# =============================================================================
# 市町村の田面積加重標高・傾斜を作る（T4: horse race の競合変数）
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   python Data_cleaning/13build_city_elevation.py
#
# T4（261023_event_2010_horse_race_T4.md）で、heat × 経度の2010年交差項が
# t = 4.1〜4.4 で頑健に効き、heat × 緯度・経度を年飽和で入れると
# delta_t の同質性が棄却されなくなった。経度が何の代理かを特定する必要があり、
# 第一候補が標高である（関東は西が山地、東が平野）。
#
# 入力: Elevation_data/G04-a-11_*_GML.zip
#   国土数値情報 G04-a「標高・傾斜度3次メッシュデータ」（平成11年、JGD2000）
#   G04a_001 メッシュコード / 002 平均標高 / 003 最高標高 / 004 最低標高
#   006 最大傾斜角度 / 008 最小傾斜角度 / 010 平均傾斜角度
#   海側などのメッシュは "unknown" が入るので数値化して落とす（関東11枚で4,431件）
#
# 気象集計と同一の重み行列 W（ピクセルレベルの田面積 valid_area）を標高に当てる。
# 11build_city_centroids.py と同じ経路で、基準年も同じ 2006 年。
#
# メッシュ農業気象データのグリッドは標準地域メッシュ（3次＝約1km）そのもので、
# ラスタの CRS も EPSG:4612（JGD2000）なので標高側と測地系が一致する。
# したがってピクセル中心の経緯度から3次メッシュコードを算術で決めて結合できる。
#
# 出力: results/city_elevation.csv
#   city_id, elev, slope, relief（田面積加重）, elev_geom, slope_geom（単純平均）,
#   elev_p10, elev_p90, n_pixels, n_matched, matched_share
# =============================================================================

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

YEAR = 2006
ZDIR = Path("Elevation_data")
OUT = Path("results/city_elevation.csv")

# -----------------------------------------------------------------------------
# 1. 3次メッシュ標高テーブル
# -----------------------------------------------------------------------------
REN = {"G04a_001": "mesh", "G04a_002": "elev_mean", "G04a_003": "elev_max",
       "G04a_004": "elev_min", "G04a_006": "slope_max", "G04a_008": "slope_min",
       "G04a_010": "slope_mean"}
NUM = ["elev_mean", "elev_max", "elev_min", "slope_max", "slope_min", "slope_mean"]

frames = []
for z in sorted(glob.glob(str(ZDIR / "G04-a-11_*_GML.zip"))):
    g = gpd.read_file("zip://" + z, encoding="cp932").drop(columns="geometry")
    frames.append(g.rename(columns=REN))
raw = pd.concat(frames, ignore_index=True)
n_raw = len(raw)
raw["mesh"] = pd.to_numeric(raw["mesh"]).astype("int64")
for c in NUM:
    raw[c] = pd.to_numeric(raw[c], errors="coerce")
mesh = raw[["mesh"] + NUM].dropna(subset=["elev_mean"])
assert mesh["mesh"].is_unique, "数値行に重複メッシュがある"
print(f"標高メッシュ: 読込 {n_raw} 行 → 数値 {len(mesh)} 行"
      f"（unknown {n_raw - len(mesh)} 件を除外）", flush=True)

# -----------------------------------------------------------------------------
# 2. 田面積の重みとピクセル座標（11build_city_centroids.py と同じ手口）
# -----------------------------------------------------------------------------
_src = Path("Data_cleaning/03b_add_tmp_min.py").read_text(encoding="utf-8")
_ns = {}
exec(compile(_src.split("# %%\nfor year in YEARS:")[0], "03b_functions", "exec"), _ns)

get_kanto_admin_gdf = _ns["get_kanto_admin_gdf"]
select_year_farm    = _ns["select_year_farm"]
rasterize_admin     = _ns["rasterize_admin"]
rasterize_farm      = _ns["rasterize_farm"]
calculate_weights   = _ns["calculate_weights"]

gdf_admin, admin_combo = get_kanto_admin_gdf(YEAR)
gdf_farm, farm_year = select_year_farm(YEAR)
print("admin_combo =", admin_combo, "/ farm_year =", farm_year, flush=True)

with xr.open_dataset("./MetData_Output/TMP_mea/data_" + str(YEAR) + ".nc") as ds_ref:
    ds_ref.rio.write_crs("EPSG:4612", inplace=True)
    grid_admin = rasterize_admin(gdf_admin, ds_ref)
    gdf_farm = gdf_farm.to_crs("EPSG:4612")
    grid_paddy, grid_other = rasterize_farm(gdf_farm, ds_ref, farm_year)
    valid_area, valid_city_id, valid_mask = calculate_weights(grid_admin, grid_paddy,
                                                              grid_other)
    lon2d, lat2d = np.meshgrid(ds_ref["x"].values, ds_ref["y"].values)
    assert lon2d.shape == grid_admin["city_id"].values.shape
    lon = lon2d.flatten()[valid_mask]
    lat = lat2d.flatten()[valid_mask]

print("有効ピクセル:", len(valid_city_id),
      "/ 市町村:", len(np.unique(valid_city_id)), flush=True)


# -----------------------------------------------------------------------------
# 3. 経緯度 → 3次メッシュコード
#    1次: p = floor(lat*1.5), u = floor(lon)-100
#    2次: q = floor((lat*1.5-p)*8), v = floor((lon-100-u)*8)
#    3次: r = floor(((lat*1.5-p)*8-q)*10), w = floor(((lon-100-u)*8-v)*10)
#    コードは p(2) u(2) q(1) v(1) r(1) w(1) の8桁。
#    ピクセル中心はメッシュ中心なので境界の取り違えは起きない。
# -----------------------------------------------------------------------------
def mesh_code(lon, lat):
    a = lat * 1.5
    p = np.floor(a)
    b = lon - 100.0
    u = np.floor(b)
    a = (a - p) * 8.0
    q = np.floor(a)
    b = (b - u) * 8.0
    v = np.floor(b)
    r = np.floor((a - q) * 10.0)
    w = np.floor((b - v) * 10.0)
    return (p * 1_000_000 + u * 10_000 + q * 1_000
            + v * 100 + r * 10 + w).astype("int64")


px = pd.DataFrame({"city_id": valid_city_id, "w": valid_area,
                   "lon": lon, "lat": lat})
px["mesh"] = mesh_code(px["lon"].to_numpy(), px["lat"].to_numpy())
px = px.merge(mesh, on="mesh", how="left")
px["matched"] = px["elev_mean"].notna()
print(f"メッシュ結合: {px['matched'].sum()} / {len(px)} "
      f"（{100 * px['matched'].mean():.2f}%）", flush=True)
if px["matched"].mean() < 0.98:
    print("警告: 結合率が低い。メッシュコードの算術かグリッドの整合を確認すること",
          flush=True)

px["relief"] = px["elev_max"] - px["elev_min"]

# -----------------------------------------------------------------------------
# 4. 田面積加重で市町村に集計
#    結合できなかったピクセルは重みから落とす（分母も合わせる）
# -----------------------------------------------------------------------------
m = px[px["matched"]].copy()
for c in ["elev_mean", "slope_mean", "relief"]:
    m["w_" + c] = m["w"] * m[c]

g = m.groupby("city_id", as_index=False).agg(
    wsum=("w", "sum"), n_matched=("w", "size"),
    s_elev=("w_elev_mean", "sum"), s_slope=("w_slope_mean", "sum"),
    s_relief=("w_relief", "sum"),
    elev_geom=("elev_mean", "mean"), slope_geom=("slope_mean", "mean"),
    elev_p10=("elev_mean", lambda s: s.quantile(0.10)),
    elev_p90=("elev_mean", lambda s: s.quantile(0.90)))
g["elev"] = g["s_elev"] / g["wsum"]
g["slope"] = g["s_slope"] / g["wsum"]
g["relief"] = g["s_relief"] / g["wsum"]
g = g.drop(columns=["s_elev", "s_slope", "s_relief"])

npx = px.groupby("city_id", as_index=False).agg(n_pixels=("w", "size"))
g = g.merge(npx, on="city_id", how="left")
g["matched_share"] = g["n_matched"] / g["n_pixels"]

cols = ["city_id", "elev", "slope", "relief", "elev_geom", "slope_geom",
        "elev_p10", "elev_p90", "n_pixels", "n_matched", "matched_share"]
g = g[cols].sort_values("city_id")

OUT.parent.mkdir(parents=True, exist_ok=True)
g.to_csv(OUT, index=False)
print(f"\n保存: {OUT}  （{len(g)} 市町村）")
print(g[["elev", "slope", "relief", "matched_share"]].describe().to_string())

# -----------------------------------------------------------------------------
# 5. 重心・処置変数との関係（T4 の問い：経度は何の代理か）
# -----------------------------------------------------------------------------
cen = Path("results/city_centroids.csv")
if cen.exists():
    c = pd.read_csv(cen)[["city_id", "lon", "lat"]]
    j = g.merge(c, on="city_id", how="inner")
    print(f"\n重心と結合: {len(j)} 市町村")
    print("相関:\n", j[["elev", "slope", "relief", "lon", "lat"]].corr().round(3)
          .to_string())
    print("\n田面積加重標高と単純平均標高の差（m）: "
          f"平均 {np.mean(j['elev'] - j['elev_geom']):.1f} / "
          f"中央 {np.median(j['elev'] - j['elev_geom']):.1f}")
    print("田の標高は市町村平均より低いはず（山地を含まないので）")

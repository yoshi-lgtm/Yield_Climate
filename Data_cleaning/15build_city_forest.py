# =============================================================================
# 市町村の森林率を作る（T4: 地理勾配の機序の特定）
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   python Data_cleaning/15build_city_forest.py
#
# T4 第2.5-2.6節で、2010年の高温被害の地理勾配は標高ではなく海岸距離
# （海洋性）だと分かった。森林は標高と別の機序（日中の遮蔽・冷気流出）を
# 持ちうるので、第3の候補として測る。
#
# 入力: merged_data.gpkg（03b と同じ土地利用細分メッシュの統合版）
#   2006年版は L03-a の 1991-2006 スキーマで、カラムは土地利用コード順に
#   L03a_002=0100 田 / L03a_003=0200 その他の農用地 / L03a_004=0500 森林 /
#   L03a_005=0600 荒地 / L03a_006=0700 建物用地 / ...
#   森林が L03a_004 であることは実証的に確認した（2006年メッシュで平均
#   580,299 m2＝58%、標高との相関 +0.653 で他カラムを大きく上回る。
#   次点は L03a_005 荒地の +0.227）。
#
# 2つの測り方を出す。
#   forest_city : 市町村の全メッシュでの森林率（森林 / 全地目）。
#                 「山がちな市町村か」を測る。elev と同じ性格
#   forest_nbr  : 田ピクセルの 3x3 近傍（約9km2）の森林率を田面積加重平均。
#                 「田の周りがどれだけ森林か」を測る。局所的な遮蔽・冷気の代理
#   同一メッシュ内の森林率は田と面積を奪い合うので機械的に負相関する。
#   3x3 近傍にするのはそれを緩和するためで、forest_same も参考に出す。
#
# 出力: results/city_forest.csv
# =============================================================================

from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from scipy.ndimage import uniform_filter
from geocube.api.core import make_geocube

YEAR = 2006
OUT = Path("results/city_forest.csv")
FOREST = "L03a_004"
# 1991-2006 スキーマの全地目（メッシュコード L03a_001 を除く）
ALL_LU = [f"L03a_{i:03d}" for i in range(2, 14)]

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
assert farm_year == 2006, "1991-2006 スキーマを前提にしている"

for c in ALL_LU:
    gdf_farm[c] = pd.to_numeric(gdf_farm[c], errors="coerce").fillna(0)
gdf_farm["_total_lu"] = gdf_farm[ALL_LU].sum(axis=1)
sh = gdf_farm[FOREST].sum() / gdf_farm["_total_lu"].sum()
print(f"森林カラム {FOREST}: 全体の {100 * sh:.1f}%", flush=True)

with xr.open_dataset("./MetData_Output/TMP_mea/data_" + str(YEAR) + ".nc") as ds_ref:
    ds_ref.rio.write_crs("EPSG:4612", inplace=True)
    grid_admin = rasterize_admin(gdf_admin, ds_ref)
    gdf_farm = gdf_farm.to_crs("EPSG:4612")
    grid_paddy, grid_other = rasterize_farm(gdf_farm, ds_ref, farm_year)
    valid_area, valid_city_id, valid_mask = calculate_weights(grid_admin, grid_paddy,
                                                              grid_other)
    gl = make_geocube(vector_data=gdf_farm, measurements=[FOREST, "_total_lu"],
                      like=ds_ref.isel(time=0), fill=0)
    forest2d = np.nan_to_num(gl[FOREST].values.astype("float64"))
    total2d = np.nan_to_num(gl["_total_lu"].values.astype("float64"))

print("有効ピクセル:", len(valid_city_id),
      "/ 市町村:", len(np.unique(valid_city_id)), flush=True)

# 3x3 近傍（uniform_filter は平均なので、比を取れば和の比と同じ）
f3 = uniform_filter(forest2d, size=3, mode="nearest")
t3 = uniform_filter(total2d, size=3, mode="nearest")
with np.errstate(invalid="ignore", divide="ignore"):
    share_nbr2d = np.where(t3 > 0, f3 / t3, np.nan)
    share_same2d = np.where(total2d > 0, forest2d / total2d, np.nan)

city2d = grid_admin["city_id"].values
flat = lambda a: a.flatten()

# (1) 市町村全体の森林率：田かどうかによらず、市町村に属する全ピクセル
allpix = pd.DataFrame({"city_id": flat(city2d), "f": flat(forest2d),
                       "t": flat(total2d)}).dropna(subset=["city_id"])
gc = allpix.groupby("city_id", as_index=False).agg(f=("f", "sum"), t=("t", "sum"),
                                                   n_city_px=("f", "size"))
gc["forest_city"] = np.where(gc["t"] > 0, gc["f"] / gc["t"], np.nan)
gc = gc[["city_id", "forest_city", "n_city_px"]]

# (2) 田面積加重：田ピクセルの 3x3 近傍／同一メッシュ
px = pd.DataFrame({"city_id": valid_city_id, "w": valid_area,
                   "nbr": flat(share_nbr2d)[valid_mask],
                   "same": flat(share_same2d)[valid_mask]})
px = px.dropna(subset=["nbr", "same"])
px["wn"] = px["w"] * px["nbr"]
px["ws"] = px["w"] * px["same"]
gp = px.groupby("city_id", as_index=False).agg(
    wsum=("w", "sum"), n_pixels=("w", "size"),
    swn=("wn", "sum"), sws=("ws", "sum"))
gp["forest_nbr"] = gp["swn"] / gp["wsum"]
gp["forest_same"] = gp["sws"] / gp["wsum"]
gp = gp[["city_id", "forest_nbr", "forest_same", "n_pixels"]]

g = gc.merge(gp, on="city_id", how="inner").sort_values("city_id")
g["city_id"] = g["city_id"].astype("int64")

OUT.parent.mkdir(parents=True, exist_ok=True)
g.to_csv(OUT, index=False)
print(f"\n保存: {OUT}  （{len(g)} 市町村）")
print(g[["forest_city", "forest_nbr", "forest_same"]].describe().to_string())

# -----------------------------------------------------------------------------
# 既存の地理変数との関係
# -----------------------------------------------------------------------------
j = g
for p, cols in [("results/city_centroids.csv", ["lon", "lat"]),
                ("results/city_elevation.csv", ["elev", "slope"]),
                ("results/city_coastdist.csv", ["coast_km"])]:
    if Path(p).exists():
        j = j.merge(pd.read_csv(p)[["city_id"] + cols], on="city_id", how="left")
cols = [c for c in ["forest_city", "forest_nbr", "forest_same", "elev", "slope",
                    "coast_km", "lon", "lat"] if c in j]
print(f"\n相関（{len(j)} 市町村）:\n", j[cols].corr().round(3).to_string())

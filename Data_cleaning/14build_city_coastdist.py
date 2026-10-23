# =============================================================================
# 市町村の田面積加重 海岸距離を作る（T4: 経度が何の代理かの特定）
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   python Data_cleaning/14build_city_coastdist.py
#
# T4 第2.5節で標高が否定された（lon と同時に入れると標高の2010年交差項は
# t = -3.64 -> -0.46、lat も入れると -0.00）。残る第一候補が海岸距離である。
# 関東の東西軸は太平洋からの距離とほぼ一致する。
#
# 入力: coastline_data/C23-06_{08,12,13,14}_GML.zip
#   国土数値情報 C23「海岸線」（平成18年、JGD2000）。茨城08・千葉12・
#   東京13・神奈川14＝関東の臨海4都県（群馬09/栃木10/埼玉11 は内陸）。
#   東京は小笠原・伊豆諸島を含むので関東の bbox でクリップする。
#
# 気象集計と同一の重み行列 W（ピクセルレベルの田面積 valid_area）を
# 「ピクセルから最寄り海岸線までの距離」に当てる。
# 11build_city_centroids.py / 13build_city_elevation.py と同じ経路・同じ
# 基準年 2006。距離は平面直角座標系 IX 系（EPSG:2451、JGD2000）で測る。
#
# 出力: results/city_coastdist.csv
#   city_id, coast_km（田面積加重）, coast_km_geom（単純平均）,
#   coast_km_min, coast_km_max, n_pixels
# =============================================================================

from pathlib import Path
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

YEAR = 2006
ZDIR = Path("coastline_data")
OUT = Path("results/city_coastdist.csv")
SRC_CRS = "EPSG:4612"   # JGD2000 地理座標（C23 も G04-a も .prj を持たない）
MET_CRS = "EPSG:2451"   # JGD2000 / 平面直角座標系 IX 系（関東）
# 関東の bbox。東京の小笠原・伊豆諸島を落とす
BBOX = (138.4, 34.7, 141.3, 37.3)

# -----------------------------------------------------------------------------
# 1. 海岸線
# -----------------------------------------------------------------------------
parts = []
for z in sorted(glob.glob(str(ZDIR / "C23-06_*_GML.zip"))):
    g = gpd.read_file("zip://" + z, encoding="cp932")[["geometry"]]
    g = g.set_crs(SRC_CRS, allow_override=True)
    n0 = len(g)
    g = g.cx[BBOX[0]:BBOX[2], BBOX[1]:BBOX[3]]
    print(f"{Path(z).name}: {n0} → {len(g)} 本（関東 bbox でクリップ）", flush=True)
    parts.append(g)
coast = pd.concat(parts, ignore_index=True)
coast = gpd.GeoDataFrame(coast, geometry="geometry", crs=SRC_CRS).to_crs(MET_CRS)
coast = coast[~coast.geometry.is_empty & coast.geometry.notna()]
print(f"海岸線: {len(coast)} 本、総延長 {coast.length.sum() / 1000:.0f} km", flush=True)

# -----------------------------------------------------------------------------
# 2. 田面積の重みとピクセル座標（13build_city_elevation.py と同じ手口）
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
    ds_ref.rio.write_crs(SRC_CRS, inplace=True)
    grid_admin = rasterize_admin(gdf_admin, ds_ref)
    gdf_farm = gdf_farm.to_crs(SRC_CRS)
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
# 3. 最寄り海岸線までの距離
# -----------------------------------------------------------------------------
px = gpd.GeoDataFrame(
    pd.DataFrame({"city_id": valid_city_id, "w": valid_area}),
    geometry=gpd.points_from_xy(lon, lat), crs=SRC_CRS).to_crs(MET_CRS)

j = gpd.sjoin_nearest(px, coast[["geometry"]], how="left", distance_col="dist_m")
# 同距離の海岸線が複数あると行が増えるので、ピクセルごとに1本に戻す
j = j[~j.index.duplicated(keep="first")]
assert len(j) == len(px), (len(j), len(px))
px["coast_km"] = j["dist_m"].to_numpy() / 1000.0
assert px["coast_km"].notna().all(), "距離が計算できなかったピクセルがある"
print(f"距離: 最小 {px.coast_km.min():.1f} / 中央 {px.coast_km.median():.1f} / "
      f"最大 {px.coast_km.max():.1f} km", flush=True)

# -----------------------------------------------------------------------------
# 4. 田面積加重で市町村に集計
# -----------------------------------------------------------------------------
d = pd.DataFrame(px.drop(columns="geometry"))
d["wd"] = d["w"] * d["coast_km"]
g = d.groupby("city_id", as_index=False).agg(
    wsum=("w", "sum"), n_pixels=("w", "size"), swd=("wd", "sum"),
    coast_km_geom=("coast_km", "mean"),
    coast_km_min=("coast_km", "min"), coast_km_max=("coast_km", "max"))
g["coast_km"] = g["swd"] / g["wsum"]
g = g[["city_id", "coast_km", "coast_km_geom", "coast_km_min", "coast_km_max",
       "n_pixels"]].sort_values("city_id")

OUT.parent.mkdir(parents=True, exist_ok=True)
g.to_csv(OUT, index=False)
print(f"\n保存: {OUT}  （{len(g)} 市町村）")
print(g[["coast_km", "coast_km_geom", "coast_km_min", "coast_km_max"]]
      .describe().to_string())

# -----------------------------------------------------------------------------
# 5. 既存の地理変数との関係（T4 の問い：経度は何の代理か）
# -----------------------------------------------------------------------------
cen = Path("results/city_centroids.csv")
ele = Path("results/city_elevation.csv")
if cen.exists():
    j2 = g.merge(pd.read_csv(cen)[["city_id", "lon", "lat"]], on="city_id")
    if ele.exists():
        j2 = j2.merge(pd.read_csv(ele)[["city_id", "elev", "slope"]], on="city_id")
    cols = [c for c in ["coast_km", "lon", "lat", "elev", "slope"] if c in j2]
    print(f"\n相関（{len(j2)} 市町村）:\n", j2[cols].corr().round(3).to_string())
    print("\n経度が高い（東）ほど海岸に近いはずなので coast_km と lon は負相関、")
    print("標高が高いほど内陸なので coast_km と elev は正相関になるはず")

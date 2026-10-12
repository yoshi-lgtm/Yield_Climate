# =============================================================================
# 市町村の田面積加重重心を作る（T12: Conley 空間HAC の入力）
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   python Data_cleaning/11build_city_centroids.py
#
# 気象集計と同一の重み行列 W（ピクセルレベルの田面積 valid_area）を
# 緯度経度に当てるだけ。03b/03c と同じラスタライズ経路を再利用する。
#
# 基準年は 2006 年。パネルの distinct city_id 209 件を単独で全部含む
# 唯一の年である（2000-04 は 174、2010 以降は 191-196）。重心は地理量
# なので年を跨いで使い回してよい。
#
# 出力: results/city_centroids.csv
#   city_id, lon, lat（田面積加重）, lon_geom, lat_geom（単純重心）,
#   n_pixels, paddy_area, shift_km（加重と単純の距離）
# =============================================================================

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

YEAR = 2006
OUT = Path("results/city_centroids.csv")

# 03b のループ手前までを exec して関数と初期化済みオブジェクトを取り出す（03c と同じ手口）
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
    valid_area, valid_city_id, valid_mask = calculate_weights(grid_admin, grid_paddy, grid_other)

    # ピクセル座標。ラスタは (y, x) なので meshgrid の flatten 順と一致する
    lon2d, lat2d = np.meshgrid(ds_ref["x"].values, ds_ref["y"].values)
    assert lon2d.shape == grid_admin["city_id"].values.shape, (lon2d.shape,
                                                               grid_admin["city_id"].values.shape)
    lon = lon2d.flatten()[valid_mask]
    lat = lat2d.flatten()[valid_mask]

print("有効ピクセル:", len(valid_city_id), "/ 市町村:", len(np.unique(valid_city_id)), flush=True)

df = pd.DataFrame({"city_id": valid_city_id, "w": valid_area, "lon": lon, "lat": lat})
df["wlon"] = df.w * df.lon
df["wlat"] = df.w * df.lat

g = df.groupby("city_id", as_index=False).agg(
    paddy_area=("w", "sum"), n_pixels=("w", "size"),
    swlon=("wlon", "sum"), swlat=("wlat", "sum"),
    lon_geom=("lon", "mean"), lat_geom=("lat", "mean"))
g["lon"] = g.swlon / g.paddy_area
g["lat"] = g.swlat / g.paddy_area
g = g.drop(columns=["swlon", "swlat"])

# 加重重心と単純重心のずれ（km）。緯度1度=111.0km、経度1度=111.0*cos(lat)
latr = np.deg2rad(g.lat)
dx = (g.lon - g.lon_geom) * 111.0 * np.cos(latr)
dy = (g.lat - g.lat_geom) * 111.0
g["shift_km"] = np.sqrt(dx ** 2 + dy ** 2)

g = g[["city_id", "lon", "lat", "lon_geom", "lat_geom", "n_pixels", "paddy_area", "shift_km"]]
g = g.sort_values("city_id").reset_index(drop=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
g.to_csv(OUT, index=False, encoding="utf-8-sig")

print("保存:", OUT, len(g), "市町村", flush=True)
print("加重重心と単純重心のずれ(km): 平均 %.2f / 中央 %.2f / 最大 %.2f"
      % (g.shift_km.mean(), g.shift_km.median(), g.shift_km.max()), flush=True)
print("田ピクセル数: 最小 %d / 中央 %.0f / 最大 %d"
      % (g.n_pixels.min(), g.n_pixels.median(), g.n_pixels.max()), flush=True)
print(g.head(5).to_string(index=False), flush=True)

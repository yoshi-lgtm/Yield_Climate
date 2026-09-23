"""
NaN伝播バグ修正版の市町村別・日次気象データ（2009-2020、5変数）。

`phenology_heat_2026-09.md` §7.2 で特定したバグの修正:
    疎行列積 W @ data.T は NaN を伝播させるため、市町村内に NaN ピクセルが
    1つでもあるとその市町村の全期間・全変数が NaN になっていた。
    281市町村中7市町村が、68〜267ピクセル中わずか1〜8ピクセル（0.4〜6.5%）の
    NaN で全損していた:
        成田(178px中1) 土浦(94中2) 鹿嶋(68中1) かすみがうら(124中8)
        行方(167中2) 茨城町(112中1) 栃木市(267中1)
    霞ヶ浦・北浦沿岸の主要水稲地帯が推定標本から丸ごと消えていた。

    さらに 05make_df_kanto_debuged.py の heat 集計は NaN >= 34 → False となるため、
    これらの市町村の heat は「欠測」ではなく「0」として kanto.csv に入っていた。

修正: 有効（非NaN）ピクセルだけで重みを張り直し、分母も時点ごとに計算する。
      重みの構築（行政区画ラスタ×田面積ラスタ）は03と完全に同一。

既存の avg_temp/ は上書きせず、avg_temp_fixed/ に書き出す。
全期間(1980-2025)の作り直しは高コストなので、推定標本の 2009-2020 のみを対象とする。

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/03c_rebuild_weather_nanfix.py
"""
# %%
import importlib.util
import os
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray  # noqa: F401

# 重み構築のロジックは 03b と同一なので、そちらの関数を再利用する。
# （03b はループ本体がモジュール末尾にあるため、ループ手前までを exec して関数だけ取り出す）
_src = Path("Data_cleaning/03b_add_tmp_min.py").read_text(encoding="utf-8")
_ns = {}
exec(compile(_src.split("# %%\nfor year in YEARS:")[0], "03b_functions", "exec"), _ns)

# 対象年。環境変数 REBUILD_YEARS で上書きできる（例: "2000-2008,2021-2023"）
def _parse_years(spec):
    out = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


YEARS = _parse_years(os.environ.get("REBUILD_YEARS", "2009-2020"))
VARIABLES = ["TMP_mea", "TMP_max", "TMP_min", "GSR", "APCP"]
OUT_DIR = Path("./avg_temp_fixed/")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# %%
def weighted_avg_nanaware(times, valid_city_id, data_body, valid_area):
    """NaN対応の田面積加重平均。NaNピクセルを重みから外して再正規化する。"""
    from scipy.sparse import csr_matrix
    unique_cities, inverse = np.unique(valid_city_id, return_inverse=True)
    n_pixels = len(valid_city_id)
    W = csr_matrix((valid_area, (inverse, np.arange(n_pixels))),
                   shape=(len(unique_cities), n_pixels))

    X = np.asarray(data_body, dtype=float).T          # (n_pixels, n_time)
    mask = np.isfinite(X)
    numerator = np.asarray(W @ np.where(mask, X, 0.0))
    denominator = np.asarray(W @ mask.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        avg = np.where(denominator > 0, numerator / denominator, np.nan).T

    # 何ピクセルが落ちたかを診断できるよう、有効面積比も返す
    total_area = np.asarray(W.sum(axis=1)).flatten()
    with np.errstate(invalid="ignore", divide="ignore"):
        cover = denominator / total_area[:, None]
    return (pd.DataFrame(avg, index=times, columns=unique_cities),
            pd.DataFrame(cover.T, index=times, columns=unique_cities))


# %%
if __name__ == "__main__":
    diag = []
    for year in YEARS:
        out_path = OUT_DIR / ("city_weighted_avg_temp_" + str(year) + ".csv")
        if out_path.exists():
            print(str(year) + ": 既存のためスキップ", flush=True)
            continue

        gdf_admin, admin_combo = _ns["get_kanto_admin_gdf"](year)
        gdf_farm, farm_year = _ns["select_year_farm"](year)
        if gdf_admin is None or gdf_farm is None:
            print(str(year) + ": 行政区画/土地利用データなし", flush=True)
            continue

        with xr.open_dataset("./MetData_Output/TMP_mea/data_" + str(year) + ".nc") as ds_ref:
            ds_ref.rio.write_crs("EPSG:4612", inplace=True)
            city_map = {k: f"{k}_{v}"
                        for k, v in gdf_admin.set_index("city_id")["N03_004"].to_dict().items()}

            if admin_combo not in _ns["_admin_raster_cache"]:
                _ns["_admin_raster_cache"][admin_combo] = _ns["rasterize_admin"](gdf_admin, ds_ref)
            grid_admin = _ns["_admin_raster_cache"][admin_combo]

            gdf_farm = gdf_farm.to_crs("EPSG:4612")
            if farm_year not in _ns["_farm_raster_cache"]:
                _ns["_farm_raster_cache"][farm_year] = _ns["rasterize_farm"](gdf_farm, ds_ref, farm_year)
            grid_paddy, grid_other = _ns["_farm_raster_cache"][farm_year]

            valid_area, valid_city_id, valid_mask = _ns["calculate_weights"](
                grid_admin, grid_paddy, grid_other)

            frames, cover_first = [], None
            for var in VARIABLES:
                with xr.open_dataset("./MetData_Output/" + var + "/data_" + str(year) + ".nc") as ds_var:
                    # APCPは他変数より格子のカバー範囲が狭いのでds_refに揃える（03と同じ）
                    ds_var = ds_var.reindex_like(ds_ref)
                    data = ds_var[var].values
                    pixels = data.reshape(data.shape[0], -1)[:, valid_mask]
                    times = ds_var["time"].values

                avg, cover = weighted_avg_nanaware(times, valid_city_id, pixels, valid_area)
                long = (avg.melt(ignore_index=False, var_name="city_id", value_name=var)
                           .reset_index().rename(columns={"index": "time"}))
                frames.append(long.set_index(["time", "city_id"]))
                if var == "TMP_mea":
                    cover_first = cover

        df = pd.concat(frames, axis=1).reset_index()
        df["city_id"] = df["city_id"].astype(int).map(city_map)
        df.to_csv(out_path, index=False, encoding="utf-8-sig")

        cov = cover_first.mean(axis=0)
        n_partial = int((cov < 0.999).sum())
        diag.append({"year": year, "市町村数": len(cov), "NaNを含む市町村": n_partial,
                     "最小有効面積比": round(float(cov.min()), 4)})
        print(f"{year}: 保存 {out_path} ({len(df)}行) / "
              f"NaNピクセルを含む市町村 {n_partial}件", flush=True)

    if diag:
        print("\n=== 診断 ===")
        print(pd.DataFrame(diag).to_string(index=False))

    files = sorted(OUT_DIR.glob("city_weighted_avg_temp_*.csv"))
    if files:
        alldf = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
        dest = Path(f"./avg_temp_all_years/city_weighted_avg_temp_fixed_"
                    f"{min(y for y in YEARS)}_{max(y for y in YEARS)}.csv")
        alldf.to_csv(dest, index=False, encoding="utf-8-sig")
        print("\n統合完了: " + str(dest) + " (" + str(len(alldf)) + "行)")
        na = alldf.groupby("city_id")["TMP_mea"].apply(lambda s: s.isna().all())
        print("全期間NaNの市町村: " + str(int(na.sum())) + "件")

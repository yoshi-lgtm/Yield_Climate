"""
03の補遺: TMP_min の市町村別・田面積加重平均を追加で構築する。

03calculate_all_def_debuged.py の target_weather は
['TMP_mea','TMP_max','GSR','APCP'] のみで TMP_min を含まないため、
avg_temp_all_years/ には日最低気温がない。
Schlenker & Roberts 型の日内気温分布（正弦補間）による気温ビン推定には
TMP_min が必須なので、03と同一の重み（行政区画ラスタ×田面積ラスタ）を
再構築して TMP_min だけを集計する。

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/03b_add_tmp_min.py
"""
# %%
from pathlib import Path
import re
import zipfile
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray  # noqa: F401  (xarray に .rio アクセサを生やすために必要)
from geocube.api.core import make_geocube
from scipy.sparse import csr_matrix

# %%
YEARS = range(2009, 2021)          # 圃場整備データとマージ可能な推定標本期間
shp_dir = Path("./ShapeFile/")
path_farm_gpkg = Path("merged_data.gpkg")
output_dir = Path("./avg_temp_min/")
output_dir.mkdir(parents=True, exist_ok=True)

PREFECTURES = ["chiba", "gunma", "ibaraki", "kanagawa", "saitama", "tochigi", "tokyo"]
ADMIN_NAME_RE = re.compile(r"N03-(\d+)_\d{2}_GML")
MEASURED_YEARS = [1976, 1987, 1991, 1997, 2006, 2009, 2014, 2016, 2021]


# %%
# --- 以下、03と同一のロジック（重みを1ピクセルも違えないためコピーして使う） ---
def parse_admin_year(name):
    m = ADMIN_NAME_RE.search(name)
    if not m:
        return None
    digits = m.group(1)
    if len(digits) >= 8:
        return int(digits[:4])
    yy = int(digits[:2])
    return 1900 + yy if yy >= 50 else 2000 + yy


def list_admin_sources(admin_dir):
    sources = []
    for item in admin_dir.iterdir():
        if item.is_dir() or item.suffix.lower() == ".zip":
            year = parse_admin_year(item.name)
            if year is not None:
                sources.append((year, item))
    return sorted(sources)


def select_admin_year(year, available_years):
    candidate_years = [y for y in available_years if y <= year]
    if candidate_years:
        return max(candidate_years)
    return None


def load_shp_auto(path):
    for enc in ["utf-8", "cp932", "shift_jis"]:
        try:
            gdf = gpd.read_file(path, encoding=enc)
            test = gdf["N03_004"].dropna().head(5).astype(str)
            if "å" not in test and "蟳" not in test:
                return gdf
        except Exception:
            continue
    raise ValueError("encoding判定失敗: " + str(path))


def admin_source_to_gdal_path(path):
    if path.is_dir():
        shp_path = next(p for p in path.rglob("*.shp") if not p.name.startswith("._"))
        return str(shp_path)
    with zipfile.ZipFile(path) as zf:
        shp_name = next(
            n for n in zf.namelist()
            if n.lower().endswith(".shp") and not Path(n).name.startswith("._")
        )
    return "/vsizip/" + path.resolve().as_posix() + "/" + shp_name


def get_crs_for_year(gdf_admin, year):
    if year < 2016:
        return gdf_admin.set_crs('EPSG:4612')
    return gdf_admin.set_crs('EPSG:6668').to_crs('EPSG:4612')


def rasterize_admin(gdf_admin, ds_ref):
    return make_geocube(
        vector_data=gdf_admin,
        measurements=["city_id"],
        like=ds_ref.isel(time=0),
        fill=np.nan,
    )


def rasterize_farm(gdf_farm, ds_ref, farm_year):
    if farm_year <= 1987:
        paddy_col, other_col = 'L03a_002', '_other_total'
        for col in ['L03a_002', 'L03a_003', 'L03a_004', 'L03a_005']:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)
        gdf_farm[other_col] = gdf_farm['L03a_003'] + gdf_farm['L03a_004'] + gdf_farm['L03a_005']
    elif farm_year <= 2006:
        paddy_col, other_col = 'L03a_002', 'L03a_003'
        for col in [paddy_col, other_col]:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)
    else:
        paddy_col, other_col = '田', '他農用地'
        for col in [paddy_col, other_col]:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)

    grid = make_geocube(
        vector_data=gdf_farm,
        measurements=[paddy_col, other_col],
        like=ds_ref.isel(time=0),
        fill=0,
    )
    return grid[paddy_col], grid[other_col]


def calculate_weights(grid_admin, grid_paddy, grid_other):
    area_flat = grid_paddy.values.flatten()
    city_id_flat = grid_admin['city_id'].values.flatten()
    valid_mask = (~np.isnan(city_id_flat)) & (area_flat > 0)
    valid_area = area_flat[valid_mask]
    valid_city_id = city_id_flat[valid_mask].astype(int)
    return valid_area, valid_city_id, valid_mask


def calculate_weighted_avg(times, valid_city_id, data_body, valid_area):
    unique_cities, inverse = np.unique(valid_city_id, return_inverse=True)
    n_pixels = len(valid_city_id)
    W = csr_matrix((valid_area, (inverse, np.arange(n_pixels))),
                   shape=(len(unique_cities), n_pixels))
    # Bug fix: 03 と同じ NaN 伝播バグへの対処。有効ピクセルだけで重みを張り直す。
    X = np.asarray(data_body, dtype=float).T
    mask = np.isfinite(X)
    numerator = np.asarray(W @ np.where(mask, X, 0.0))
    denominator = np.asarray(W @ mask.astype(float))
    with np.errstate(invalid='ignore', divide='ignore'):
        avg = np.where(denominator > 0, numerator / denominator, np.nan).T
    return pd.DataFrame(avg, index=times, columns=unique_cities)


# %%
admin_sources_by_pref = {p: list_admin_sources(shp_dir / ("admin_" + p)) for p in PREFECTURES}
gdf_master = gpd.read_file(path_farm_gpkg)
print("gpkg 読み込み完了", gdf_master.shape, flush=True)

_admin_gdf_cache = {}
_kanto_admin_cache = {}
_admin_raster_cache = {}
_farm_raster_cache = {}


def get_admin_gdf(pref, year):
    sources = admin_sources_by_pref[pref]
    admin_year = select_admin_year(year, [y for y, _ in sources])
    if admin_year is None:
        return None, None
    key = (pref, admin_year)
    if key not in _admin_gdf_cache:
        path = next(p for y, p in sources if y == admin_year)
        gdf_adm = load_shp_auto(admin_source_to_gdal_path(path))
        gdf_adm = get_crs_for_year(gdf_adm, admin_year)
        gdf_adm = gdf_adm.dropna(subset=['N03_007'])
        gdf_adm['city_id'] = gdf_adm['N03_007'].astype(int)
        gdf_adm['N03_004'] = gdf_adm['N03_004'].str.replace(r'(市|町|村)$', '', regex=True)
        _admin_gdf_cache[key] = gdf_adm
    return _admin_gdf_cache[key], admin_year


def get_kanto_admin_gdf(year):
    per_pref, combo = [], []
    for pref in PREFECTURES:
        g, ay = get_admin_gdf(pref, year)
        if g is None:
            continue
        per_pref.append(g)
        combo.append((pref, ay))
    if not per_pref:
        return None, None
    combo = tuple(combo)
    if combo not in _kanto_admin_cache:
        _kanto_admin_cache[combo] = gpd.GeoDataFrame(
            pd.concat(per_pref, ignore_index=True), crs=per_pref[0].crs)
    return _kanto_admin_cache[combo], combo


def select_year_farm(year):
    cands = [y for y in MEASURED_YEARS if y <= year]
    if not cands:
        return None, None
    ty = max(cands)
    return gdf_master[gdf_master['year'] == ty].copy(), ty


# %%
for year in YEARS:
    out_path = output_dir / ("city_weighted_tmp_min_" + str(year) + ".csv")
    if out_path.exists():
        print(str(year) + ": 既存のためスキップ", flush=True)
        continue

    gdf_admin, admin_combo = get_kanto_admin_gdf(year)
    gdf_farm, farm_year = select_year_farm(year)
    if gdf_admin is None or gdf_farm is None:
        print(str(year) + ": 行政区画/土地利用データなし", flush=True)
        continue

    with xr.open_dataset("./MetData_Output/TMP_mea/data_" + str(year) + ".nc") as ds_ref:
        ds_ref.rio.write_crs('EPSG:4612', inplace=True)
        # 03 と同一の f-string 形式（NaN は "nan" になるが既存CSVと整合させる）
        city_map = {k: f"{k}_{v}"
                    for k, v in gdf_admin.set_index('city_id')['N03_004'].to_dict().items()}

        if admin_combo not in _admin_raster_cache:
            _admin_raster_cache[admin_combo] = rasterize_admin(gdf_admin, ds_ref)
        grid_admin = _admin_raster_cache[admin_combo]

        gdf_farm = gdf_farm.to_crs('EPSG:4612')
        if farm_year not in _farm_raster_cache:
            _farm_raster_cache[farm_year] = rasterize_farm(gdf_farm, ds_ref, farm_year)
        grid_paddy, grid_other = _farm_raster_cache[farm_year]

        valid_area, valid_city_id, valid_mask = calculate_weights(grid_admin, grid_paddy, grid_other)

        with xr.open_dataset("./MetData_Output/TMP_min/data_" + str(year) + ".nc") as ds_var:
            ds_var = ds_var.reindex_like(ds_ref)
            data = ds_var['TMP_min'].values
            data_flat = data.reshape(data.shape[0], -1)
            weather_pixels = data_flat[:, valid_mask]
            times = ds_var['time'].values

    df_avg = calculate_weighted_avg(times, valid_city_id, weather_pixels, valid_area)
    df_long = (df_avg.melt(ignore_index=False, var_name='city_id', value_name='TMP_min')
                     .reset_index().rename(columns={'index': 'time'}))
    df_long['city_id'] = df_long['city_id'].astype(int).map(city_map)
    df_long.to_csv(out_path, index=False, encoding='utf-8-sig')
    print(str(year) + ": 保存 " + str(out_path) + " (" + str(len(df_long)) + "行)", flush=True)

# %%
files = sorted(output_dir.glob("city_weighted_tmp_min_*.csv"))
if files:
    all_df = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    dest = Path("./avg_temp_all_years/city_weighted_tmp_min_2009_2020.csv")
    all_df.to_csv(dest, index=False, encoding='utf-8-sig')
    print("統合完了: " + str(dest) + " (" + str(len(all_df)) + "行)", flush=True)

# %%
# Import necessary libraries
from pathlib import Path
import re
import zipfile
import geopandas as gpd
import xarray as xr
import pandas as pd
import numpy as np
from geocube.api.core import make_geocube
import rioxarray
from scipy.sparse import csr_matrix


# %%
# パス指定
years = range(1980, 2026)

shp_dir = Path("./ShapeFile/")
path_farm_gpkg = Path("merged_data.gpkg")
output_dir = Path("./avg_temp/")
output_dir_paddy = Path("./paddy/")

# Bug fix: 出力ディレクトリが存在しないとCSV保存時にエラーになるため作成
output_dir.mkdir(parents=True, exist_ok=True)
output_dir_paddy.mkdir(parents=True, exist_ok=True)

#%%
gdf_master = gpd.read_file(path_farm_gpkg)

#%%
# 関東7都県: 行政区画データ（zipのまま・展開しない）の年インデックスを作成
# 各県フォルダは admin_<pref>/N03-<日付>_<県コード>_GML.zip という命名の
# 国土数値情報N03形式。市町村合併のタイミングは県ごとに異なるため、
# 年の選定・キャッシュも県ごとに独立して行う。
PREFECTURES = ["chiba", "gunma", "ibaraki", "kanagawa", "saitama", "tochigi", "tokyo"]

ADMIN_NAME_RE = re.compile(r"N03-(\d+)_\d{2}_GML")


def parse_admin_year(name):
    """N03-<日付>_<県コード>_GML の<日付>部分(2/6/8桁)から西暦年を復元する"""
    m = ADMIN_NAME_RE.search(name)
    if not m:
        return None
    digits = m.group(1)
    if len(digits) >= 8:
        return int(digits[:4])
    yy = int(digits[:2])
    return 1900 + yy if yy >= 50 else 2000 + yy


def list_admin_sources(admin_dir):
    """zip・展開済みフォルダのどちらでも (年, パス) のリストを返す"""
    sources = []
    for item in admin_dir.iterdir():
        if item.is_dir() or item.suffix.lower() == ".zip":
            year = parse_admin_year(item.name)
            if year is not None:
                sources.append((year, item))
    return sorted(sources)


admin_sources_by_pref = {
    pref: list_admin_sources(shp_dir / f"admin_{pref}") for pref in PREFECTURES
}
for pref, sources in admin_sources_by_pref.items():
    print(pref, [y for y, _ in sources])

#%%
# 直近過去の行政区画データを選ぶ関数
def select_admin_year(year, available_years):
    candidate_years = [y for y in available_years if y <= year]
    if candidate_years:
        admin_year = max(candidate_years)
        print(f"{year}年 → 使用する市町村界: {admin_year}")
        return admin_year
    # Bug fix: return None の後にあったため到達不能だった print を前に移動
    print(f"{year}年に対応する過去データがありません")
    return None

#%%
# シェープファイルのエンコーディング自動判定関数
def load_shp_auto(path):
    for enc in ["utf-8", "cp932", "shift_jis"]:
        try:
            gdf = gpd.read_file(path, encoding=enc)
            test = gdf["N03_004"].dropna().head(5).astype(str)
            if "å" not in test and "螳" not in test:
                print(f"成功: {enc}")
                return gdf
        except:
            continue
    raise ValueError("encoding判定失敗")

#%%
# zipを展開せず、GDALのvsizip経由でzip内のshpを直接読むためのパスに変換する
# Bug fix: macOSでコピーした際のAppleDouble隠しファイル（._*.shp）が本物のshpと
#          混在しているケースがあり、name起因の走査順で誤って選ばれると
#          GDALが読めず「encoding判定失敗」に見えるエラーになるため除外する
def admin_source_to_gdal_path(path):
    if path.is_dir():
        shp_path = next(p for p in path.rglob("*.shp") if not p.name.startswith("._"))
        return str(shp_path)
    with zipfile.ZipFile(path) as zf:
        shp_name = next(
            n for n in zf.namelist()
            if n.lower().endswith(".shp") and not Path(n).name.startswith("._")
        )
    return f"/vsizip/{path.resolve().as_posix()}/{shp_name}"

#%%
# 県ごとの行政区画データ取得（(県, 採用年)単位でキャッシュ）
_admin_gdf_cache = {}  # (pref, admin_year) -> 処理済みgdf

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

#%%
# 関東7都県分をまとめた行政区画データ（県ごとの選定結果の組み合わせでキャッシュ）
_kanto_admin_cache = {}  # ((pref, admin_year), ...) -> 結合gdf

def get_kanto_admin_gdf(year):
    per_pref_gdfs = []
    combo = []
    for pref in PREFECTURES:
        gdf_pref, admin_year = get_admin_gdf(pref, year)
        if gdf_pref is None:
            print(f"{year}年: {pref} は行政区画データなしのため除外")
            continue
        per_pref_gdfs.append(gdf_pref)
        combo.append((pref, admin_year))

    if not per_pref_gdfs:
        return None, None

    combo = tuple(combo)
    if combo not in _kanto_admin_cache:
        merged = pd.concat(per_pref_gdfs, ignore_index=True)
        _kanto_admin_cache[combo] = gpd.GeoDataFrame(merged, crs=per_pref_gdfs[0].crs)

    return _kanto_admin_cache[combo], combo

#%%
# 土地利用データの直近過去を選ぶ関数
MEASURED_YEARS = [1976, 1987, 1991, 1997, 2006, 2009, 2014, 2016, 2021]

def select_year_farm(year, gdf_master):
    candidate_years = [y for y in MEASURED_YEARS if y <= year]
    if candidate_years:
        target_year = max(candidate_years)
        selected_gdf = gdf_master[gdf_master['year'] == target_year].copy()
        print(f"{year}年 → 使用する土地利用データ年度: {target_year}")
        # target_year もキャッシュキーとして返すよう変更
        return selected_gdf, target_year
    print(f"{year}年に対応する過去データがありません")
    return None, None

#%%
# CRS指定関数
# Bug fix: 元は gdf_admin をグローバル変数から暗黙参照していたが、引数として受け取るよう修正
# CRS は shapefile 自体の作成年（admin_year）で判断するのが正しい
def get_crs_for_year(gdf_admin, year):
    if year < 2016:
        return gdf_admin.set_crs('EPSG:4612')
    return gdf_admin.set_crs('EPSG:6668').to_crs('EPSG:4612')

#%%
# グリッドとベクターデータを重ねるための関数
def rasterize_admin(gdf_admin, ds_ref):
    return make_geocube(
        vector_data=gdf_admin,
        measurements=["city_id"],
        like=ds_ref.isel(time=0),
        fill=np.nan
    )

#%%
# 農地のラスタ化関数
# Performance: 元は make_geocube を田・他農用地で2回呼んでいたが、1回にまとめて高速化
# Bug fix: 各版のメタデータXML（KS-META同梱のtky.xml、gml:CompositeValueの
#          valueComponents）で列の意味を実データと突き合わせたところ、
#          L03a_XXX形式の列は "L03a_001 = 3次メッシュコード" であり、
#          実際の面積値は L03a_002 から始まっていた（従来コードは1列ズレたまま
#          L03a_001を田の面積として扱っており、メッシュコードという意味のない
#          巨大な数値を面積として集計していた）。
#          さらに年代によって面積列の分類自体も3パターンある。
#          - 1976・1987版: 田(002)/畑(003)/果樹園(004)/その他樹木畑(005) の4分類。
#            合算しないと「その他の農用地」が過小評価される。
#          - 1991・1997・2006版: 田(002)/その他の農用地(003) にすでに集約済み。
#            004以降は森林・荒地等で非農地のため合算してはいけない
#            （元コードはこの3年度を「else」節で処理し、存在しない'他農用地'列を
#              参照してしまい、有効ピクセルが0件→出力が空になっていた）。
#          - 2009版以降: 列名がGML側で日本語ラベル('田'/'他農用地')に変更され、
#            メッシュコードも別列('メッシュ')に分離されている（ズレなし）。
def rasterize_farm(gdf_farm, ds_ref, farm_year):
    if farm_year <= 1987:
        # 1976・1987: 田(002) / 畑(003)・果樹園(004)・その他樹木畑(005)を合算
        paddy_col = 'L03a_002'
        other_col = '_other_total'
        for col in ['L03a_002', 'L03a_003', 'L03a_004', 'L03a_005']:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)
        gdf_farm[other_col] = gdf_farm['L03a_003'] + gdf_farm['L03a_004'] + gdf_farm['L03a_005']
    elif farm_year <= 2006:
        # 1991・1997・2006: 田(002) / その他の農用地(003) にすでに集約済み
        paddy_col, other_col = 'L03a_002', 'L03a_003'
        for col in [paddy_col, other_col]:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)
    else:
        # 2009版以降: 列名が日本語ラベルに変更
        paddy_col, other_col = '田', '他農用地'
        for col in [paddy_col, other_col]:
            gdf_farm[col] = pd.to_numeric(gdf_farm[col], errors='coerce').fillna(0)

    grid = make_geocube(
        vector_data=gdf_farm,
        measurements=[paddy_col, other_col],
        like=ds_ref.isel(time=0),
        fill=0
    )
    return grid[paddy_col], grid[other_col]

#%%
# 重み計算の関数
def calculate_weights(grid_admin, grid_paddy, grid_other):
    area_flat = grid_paddy.values.flatten()
    nonpaddy_flat = grid_other.values.flatten()
    total_flat = area_flat + nonpaddy_flat
    city_id_flat = grid_admin['city_id'].values.flatten()

    valid_mask = (~np.isnan(city_id_flat)) & (area_flat > 0)
    total_valid_mask = (~np.isnan(city_id_flat)) & (total_flat > 0)

    valid_area = area_flat[valid_mask]
    valid_city_id = city_id_flat[valid_mask].astype(int)
    total_valid_area = total_flat[total_valid_mask]
    total_valid_city_id = city_id_flat[total_valid_mask]

    return valid_area, valid_city_id, total_valid_area, total_valid_city_id, valid_mask

#%%
# 市町村ごとの時間変化する重み付き平均を計算する関数
# Bug fix: valid_area をグローバル参照していたが引数として受け取るよう修正
# Performance: タイムステップごとに DataFrame+groupby を繰り返す O(T×N) ループを
#              scipy 疎行列による行列積 O(1) に置き換え（46年×365日×N市町村 → 大幅高速化）
def calculate_weighted_avg_temp(ds_weather, valid_city_id, data_body, valid_area):
    times = ds_weather['time'].values
    unique_cities, inverse = np.unique(valid_city_id, return_inverse=True)
    n_cities = len(unique_cities)
    n_pixels = len(valid_city_id)

    # 疎行列 W[c, p] = area[p]  (pixel p が city c に属する場合)
    W = csr_matrix((valid_area, (inverse, np.arange(n_pixels))), shape=(n_cities, n_pixels))

    # Bug fix: 疎行列積は NaN を伝播させるため、市町村内に NaN ピクセルが
    #          1つでもあるとその市町村の全期間・全変数が NaN になっていた。
    #          （実際に 281市町村中7市町村が 68〜267ピクセル中わずか1〜8ピクセルの
    #            NaN で全損しており、土浦・行方・鹿嶋・かすみがうら等の
    #            霞ヶ浦沿岸の主要水稲地帯が推定標本から消えていた。）
    #          有効ピクセルだけで重みを張り直し、分母も時点ごとに計算する。
    X = np.asarray(data_body, dtype=float).T          # (n_pixels, n_time)
    mask = np.isfinite(X)
    numerator = np.asarray(W @ np.where(mask, X, 0.0))          # (n_cities, n_time)
    denominator = np.asarray(W @ mask.astype(float))            # 時点ごとの有効面積
    with np.errstate(invalid='ignore', divide='ignore'):
        avg_temp = np.where(denominator > 0, numerator / denominator, np.nan).T

    return pd.DataFrame(avg_temp, index=times, columns=unique_cities)

#%%
# 水田率の計算の関数
def calculate_paddy_area(valid_city_id, valid_area, total_valid_city_id, total_valid_area):
    df_paddy = pd.DataFrame({'city_id': valid_city_id, 'area': valid_area})
    paddy_sum_by_city = df_paddy.groupby('city_id')['area'].sum()
    df_total = pd.DataFrame({'city_id': total_valid_city_id, 'area': total_valid_area})
    total_sum_by_city = df_total.groupby('city_id')['area'].sum()
    return paddy_sum_by_city / total_sum_by_city

# %%
target_weather = ['TMP_mea', 'TMP_max', 'GSR', 'APCP']

# Performance: 同じ admin_combo / farm_year のラスタ化は毎年同じ結果なのでキャッシュ化
# （気象データのグリッドが全年度で同一であることを前提とする）
_admin_raster_cache = {}  # admin_combo -> grid_admin
_farm_raster_cache = {}   # farm_target_year -> (grid_paddy, grid_other)

for year in years:
    gdf_admin, admin_combo = get_kanto_admin_gdf(year)
    if gdf_admin is None:
        continue

    gdf_farm, farm_target_year = select_year_farm(year, gdf_master)
    if gdf_farm is None:
        continue

    # Bug fix: ds_weather を with ブロックで管理してファイルリークを防止
    with xr.open_dataset(f"./MetData_Output/TMP_mea/data_{year}.nc") as ds_ref:
        ds_ref.rio.write_crs('EPSG:4612', inplace=True)

        city_map = {k: f"{k}_{v}" for k, v in gdf_admin.set_index('city_id')['N03_004'].to_dict().items()}

        # 行政区画ラスタ化（県の組み合わせ単位でキャッシュ）
        if admin_combo not in _admin_raster_cache:
            _admin_raster_cache[admin_combo] = rasterize_admin(gdf_admin, ds_ref)
        grid_admin = _admin_raster_cache[admin_combo]

        # 農地データ（キャッシュ）
        gdf_farm = gdf_farm.to_crs('EPSG:4612')
        if farm_target_year not in _farm_raster_cache:
            grid_paddy, grid_other = rasterize_farm(gdf_farm, ds_ref, farm_target_year)
            _farm_raster_cache[farm_target_year] = (grid_paddy, grid_other)
        grid_paddy, grid_other = _farm_raster_cache[farm_target_year]

        # 重みの計算
        valid_area, valid_city_id, total_valid_area, total_valid_city_id, valid_mask = \
            calculate_weights(grid_admin, grid_paddy, grid_other)

        # 気象変数ごとに重み付き平均を計算
        all_results = []
        for var in target_weather:
            with xr.open_dataset(f"./MetData_Output/{var}/data_{year}.nc") as ds_var:
                # Bug fix: APCPは他変数(TMP_mea/TMP_max/GSR)より格子が小さく
                # (同解像度だがカバー範囲が狭い)、そのままflattenするとds_ref基準の
                # valid_mask(65016要素)と次元が合わずIndexErrorになる。
                # 座標値はds_refと完全一致しているため、reindex_likeでds_refの
                # 格子に揃え、未カバー範囲はNaN埋めする。
                ds_var = ds_var.reindex_like(ds_ref)
                data = ds_var[var].values
                data_flat = data.reshape(data.shape[0], -1)
                weather_pixels = data_flat[:, valid_mask]  # valid_mask はすでに bool 配列

            df_avg = calculate_weighted_avg_temp(ds_ref, valid_city_id, weather_pixels, valid_area)
            df_long = df_avg.melt(ignore_index=False, var_name='city_id', value_name=var).reset_index()
            df_long['city_id'] = df_long['city_id'].astype(int)
            df_long = df_long.rename(columns={'index': 'time'})
            all_results.append(df_long.set_index(['time', 'city_id']))

    # 水田率の計算（ds_ref クローズ後も numpy 配列なので問題なし）
    paddy_ratio = calculate_paddy_area(valid_city_id, valid_area, total_valid_city_id, total_valid_area)

    # 結果をDataFrameにまとめて出力
    df_result = pd.concat(all_results, axis=1).reset_index()
    # Bug fix: df_ratio = pd.DataFrame(paddy_ratio) は次行で即座に上書きされるデッドコードを削除
    df_ratio = paddy_ratio.reset_index()
    df_ratio.columns = ['city_id', 'paddy_ratio']

    df_result['city_id'] = df_result['city_id'].map(city_map)
    df_ratio['city_id'] = df_ratio['city_id'].map(city_map)

    output_path = output_dir / f"city_weighted_avg_temp_{year}.csv"
    output_path_paddy = output_dir_paddy / f"paddy_ratio_{year}.csv"
    df_result.to_csv(output_path, index=False, encoding='utf-8-sig')
    df_ratio.to_csv(output_path_paddy, index=False, encoding='utf-8-sig')

# %%
all_df = pd.concat(
    pd.read_csv(f) for f in Path(output_dir).rglob("city_weighted_avg_temp_*.csv")
)
output_dir_all = Path("./avg_temp_all_years/")
output_dir_all.mkdir(exist_ok=True)
all_df.to_csv(output_dir_all / "city_weighted_avg_temp_all_years.csv", index=False, encoding='utf-8-sig')

#%%
# 水田率の統合
all_files = list(output_dir_paddy.rglob("paddy_ratio_*.csv"))
all_files_df = []

for f in all_files:
    match = re.search(r'\d{4}', f.name)
    if not match:
        continue
    file_year = int(match.group())
    df = pd.read_csv(f).assign(year=file_year)
    all_files_df.append(df)

if all_files_df:
    df_paddy_master = pd.concat(all_files_df, axis=0, ignore_index=True)
    cols = ['year'] + [c for c in df_paddy_master.columns if c != 'year']
    df_paddy_master = df_paddy_master[cols]
    output_dir_totalratio = Path("./total_puddyratio/")
    output_dir_totalratio.mkdir(exist_ok=True)
    df_paddy_master.to_csv(output_dir_totalratio / "paddy_ratio_master.csv", index=False, encoding='utf-8-sig')
    print(f"統合完了: {len(all_files)}ファイル")
else:
    print("ファイルが見つかりませんでした。")

# %%
# 以下はインタラクティブ確認用セル（ループ最終年の変数を参照）
print(df_result.head())
df_ratio.info()
df_result.info()

# %%
df_long.head()

# %%
df_paddy_master['city_id'] = df_paddy_master['city_id'].map(city_map)
df_paddy_master.head()

# %%
# paddy_ratioの空間的な散らばりを地図で確認
# （city_idの上2桁=都道府県コードで結合しているため、気象データがまだ揃っていない
#   都県は自動的にNaN=グレー表示になり、現在のデータカバレッジの確認にも使える）
import matplotlib.pyplot as plt

map_year = 2025

gdf_map, _ = get_kanto_admin_gdf(map_year)

df_paddy_year = pd.read_csv(output_dir_totalratio / "paddy_ratio_master.csv")
df_paddy_year = df_paddy_year[df_paddy_year['year'] == map_year].copy()
df_paddy_year['city_id'] = df_paddy_year['city_id'].str.split('_', n=1).str[0].astype(int)

gdf_plot = gdf_map.merge(df_paddy_year[['city_id', 'paddy_ratio']], on='city_id', how='left')

fig, ax = plt.subplots(figsize=(9, 9))
gdf_plot.plot(
    column='paddy_ratio', cmap='YlGnBu', legend=True,
    edgecolor='gray', linewidth=0.2, ax=ax,
    missing_kwds={'color': 'lightgrey', 'label': 'no data'}
)
ax.set_title(f'paddy_ratio spatial distribution ({map_year})')
# 東京都の伊豆・小笠原諸島が入ると余白だらけになるため、関東本土の範囲にズーム
ax.set_xlim(138.3, 141.0)
ax.set_ylim(34.8, 37.3)
ax.axis('off')
plt.savefig(f'paddy_ratio_map_{map_year}.png', dpi=150, bbox_inches='tight')
plt.show()

# %%

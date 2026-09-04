#%%
import pandas as pd
import numpy as np
from pathlib import Path

#%%
# パス指定（py_files配下の相対パス。02/03/04の出力先と揃えている）
avg_temp_path = "./avg_temp_all_years/city_weighted_avg_temp_all_years.csv"
paddy_path = "./total_puddyratio/paddy_ratio_master.csv"
yield_path = "./yield/kanto_yield.csv"
output_dir = Path("./results/")
output_path = output_dir / "kanto.csv"

output_dir.mkdir(parents=True, exist_ok=True)

# city_id（N03の全国地方公共団体コード、上2桁が都道府県コード）→ 都県名
# Bug fix: 元は city_name だけで結合しており、都県をまたいだ同名市町村
#         （例: 千葉市・さいたま市の中央区など）を誤結合する恐れがあったため、
#         city_idから都県を復元し結合キーに含める
PREF_CODE_MAP = {
    8: "茨城県", 9: "栃木県", 10: "群馬県", 11: "埼玉県",
    12: "千葉県", 13: "東京都", 14: "神奈川県",
}

#%%
# 気象データの読み込み（city_id列は "cityid_市町村名" 形式なので分割し、都県も復元する）
df = pd.read_csv(avg_temp_path, index_col=0)
df[['city_id', 'city_name']] = df['city_id'].str.split('_', n=1, expand=True)
df['city_id'] = df['city_id'].astype(int)
df['prefecture'] = (df['city_id'] // 1000).map(PREF_CODE_MAP)

df = df.reset_index()
df['time'] = pd.to_datetime(df['time'])
df['year'] = df['time'].dt.year
df['month'] = df['time'].dt.month
df.head()

#%%
# 7月・8月のTMP_maxが34度以上の日数を年・市町村ごとにカウント
df['heat'] = (df['TMP_max'] >= 34) & (df['month'].isin([7, 8]))
df_heat = df.groupby(['year', 'prefecture', 'city_id', 'city_name'])['heat'].sum().reset_index(name='heat_count')
print(sorted(df_heat['heat_count'].unique()))
df_heat

#%%
# 7-8月のデータに絞って年・市町村ごとに集計
df_78 = df[df['month'].isin([7, 8])].copy()
df_analysis = df_78.groupby(['year', 'prefecture', 'city_id', 'city_name']).agg({
    'APCP': 'sum',      # 7-8月合計降水量
    'GSR': 'mean',      # 7-8月平均日射量
    'TMP_mea': 'mean',  # 7-8月平均気温
    'heat': 'sum'       # 7-8月34度以上の日数
}).reset_index()
df_analysis

#%%
# 収量データの読み込み（04で都県・市町村名に分割済みのkanto_yield.csvを使用）
# Bug fix: value列は文字列型で「－」「X」「…」等の欠測記号が混在しており、
#         そのままgroupby().sum()すると数値ではなく文字列結合になってしまうため
#         集計前にpd.to_numeric(errors='coerce')で数値化する（変換できない値はNaN）
df_yield = pd.read_csv(yield_path)
df_yield['value'] = pd.to_numeric(df_yield['value'], errors='coerce')

df_yield_kanto = (
    df_yield.groupby(['時間軸（年次）', 'prefecture', 'city_name'])['value']
    .sum(min_count=1)
    .reset_index(name='yield')
)
df_yield_kanto = df_yield_kanto.rename(columns={'時間軸（年次）': 'year'})
df_yield_kanto['year'] = df_yield_kanto['year'].astype(int)
df_yield_kanto.head()

#%%
# heatとyieldを結合（year・prefecture・city_nameで結合）
df_final = pd.merge(df_yield_kanto, df_analysis, on=['year', 'prefecture', 'city_name'], how='inner')
cols = ['yield'] + [c for c in df_final.columns if c != 'yield']
df_final = df_final[cols]
df_final

# デバッグ用: 気象データと結合できなかった収量データを確認
# （行政区画データの市町村名と収量統計の市町村名が食い違うケースを検出するため）
_matched_keys = df_analysis[['year', 'prefecture', 'city_name']].drop_duplicates()
_unmatched = df_yield_kanto.merge(_matched_keys, on=['year', 'prefecture', 'city_name'], how='left', indicator=True)
_unmatched = _unmatched[_unmatched['_merge'] == 'left_only']
print(f"収量データ行数: {len(df_yield_kanto)} → 気象データと結合後: {len(df_final)}")
if len(_unmatched):
    print(f"気象データと結合できなかった (year, prefecture, city_name) の組み合わせ: "
          f"{_unmatched[['prefecture', 'city_name']].drop_duplicates().shape[0]}件")
    print(_unmatched[['year', 'prefecture', 'city_name']].drop_duplicates().head(20))

#%%
# 水田率データの読み込み・結合
df_paddy = pd.read_csv(paddy_path)
df_paddy[['city_id', 'city_name']] = df_paddy['city_id'].str.split('_', n=1, expand=True)
df_paddy['city_id'] = df_paddy['city_id'].astype(int)
df_paddy['prefecture'] = (df_paddy['city_id'] // 1000).map(PREF_CODE_MAP)

df_finals = pd.merge(df_final, df_paddy, on=['year', 'prefecture', 'city_id', 'city_name'], how='inner')
print(f"水田率と結合後: {len(df_finals)}行 (対象都県: {sorted(df_finals['prefecture'].unique())})")
df_finals

#%%
# 出力
df_finals.to_csv(output_path, index=False, encoding='utf-8-sig')
print(f"書き出し完了: {output_path}")

# %%

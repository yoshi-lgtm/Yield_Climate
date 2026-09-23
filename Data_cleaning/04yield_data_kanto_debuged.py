#%%
import os
import re
import urllib.request
from pathlib import Path

import pandas as pd
import numpy as np

#%%
# --- 設定 ---
# e-Stat の appId。キーはコミットしないので環境変数 ESTAT_APPID に入れて使う。
#   PowerShell: $env:ESTAT_APPID = "xxxx"
#   bash:       export ESTAT_APPID=xxxx
# 取得済みCSVが揃っていれば appid なしでも後段の整形だけ再実行できるようにする
appid = os.environ.get("ESTAT_APPID", "")

# データの置き場所。環境変数 YIELD_DATA_DIR で上書きできる
DATA_DIR = Path(os.environ.get("YIELD_DATA_DIR", "."))
YIELD_DIR = DATA_DIR / "yield"
YIELD_DIR.mkdir(parents=True, exist_ok=True)

# %%
# APIから取得したデータを読み込む
urls = {
    "2023": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?cdCat02=1005%2C1006%2C1007&appId={appid}&lang=J&statsDataId=0002112120&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0",
    "2022": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?cdCat02=1005%2C1006%2C1007&appId={appid}&lang=J&statsDataId=0002111964&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0",
    "2021": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?cdCat02=1005%2C1006%2C1007&appId={appid}&lang=J&statsDataId=0001994876&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0",
    "2007-2020": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?appId={appid}&lang=J&statsDataId=0003293480&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0",
    "2006": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?cdCat01=100&appId={appid}&lang=J&statsDataId=0003284440&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0",
    "1993-2005": f"http://api.e-stat.go.jp/rest/3.0/app/getSimpleStatsData?cdCat01=100&appId={appid}&lang=J&statsDataId=0003254280&metaGetFlg=Y&cntGetFlg=N&explanationGetFlg=Y&annotationGetFlg=Y&sectionHeaderFlg=1&replaceSpChars=0"
}
# %%
# データを保存する
output_path = "api_data"

for year, url in urls.items():
    dest = YIELD_DIR / f"{output_path}_{year}.csv"
    if dest.exists():
        print(f"{year}: 取得済みのためダウンロードをスキップ")
        continue
    if not appid:
        raise RuntimeError(
            f"{dest} が無く、環境変数 ESTAT_APPID も未設定のためダウンロードできません")
    urllib.request.urlretrieve(url, dest)

# %%
#1993-2005年のデータを読み込む
data_1993_2005 = pd.read_csv(YIELD_DIR / f"{output_path}_1993-2005.csv", encoding="utf-8", skiprows=92, index_col=0)
data_1993_2005

# %%
# 2006年のデータを読み込む
data_2006 = pd.read_csv(YIELD_DIR / f"{output_path}_2006.csv", encoding="utf-8", skiprows=91, index_col=0)
data_2006

# %%
# 2007-2020年のデータを読み込む
data_2007_2020 = pd.read_csv(YIELD_DIR / f"{output_path}_2007-2020.csv", encoding="utf-8", skiprows=91, index_col=0)
data_2007_2020

# %%
# 2021年のデータを読み込む
data_2021 = pd.read_csv(YIELD_DIR / f"{output_path}_2021.csv", encoding="utf-8", skiprows=91, index_col=0)
data_2021.rename(columns={"値": "yield"}, inplace=True)
data_2021

# %%
# 2022年のデータを読み込む
data_2022 = pd.read_csv(YIELD_DIR / f"{output_path}_2022.csv", encoding="utf-8", skiprows=91, index_col=0)
data_2022

# %%
# 2023年のデータを読み込む
data_2023 = pd.read_csv(YIELD_DIR / f"{output_path}_2023.csv", encoding="utf-8", skiprows=91, index_col=0)
data_2023

# %%
# とりあえず2021から2023のデータをそろえる
# 空の辞書を作る
data_all = {}

# 読み込みの段階で辞書に入れる
years = [2021, 2022, 2023]
for year in years:
    # 読み込んだデータをそのまま辞書に保存
    df = pd.read_csv(YIELD_DIR / f"{output_path}_{year}.csv", encoding="utf-8", skiprows=91)
    
    # 💡 列名の変更と年度の追加
    new_cols = []
    for col in df.columns:
        if "全国" in col: 
            new_cols.append("市町村")
        elif "北海道以外" in col: 
            new_cols.append("面積収量収穫量")
        else: 
            new_cols.append(col)
    
    df.columns = new_cols
    df['時間軸（年次）'] = year
    
    # 辞書に格納
    data_all[year] = df

# Bug fix: 2022年産以降のAPIは市町村名に都道府県名が付かない（2021は「茨城県_水戸市」、
#          2022以降は「水戸市」）。後段の抽出は「都道府県名+市町村名」で突き合わせるため、
#          そのままだと2022-2023が1行も通らなかった。
#          cat01_code は年をまたいで安定している（関東を含む1,719コードが完全対応）ので、
#          都道府県名付きの年から対応表を作って補完する。
_has_pref = {y: d['市町村'].astype(str).str.contains('_').mean() for y, d in data_all.items()}
_ref_year = max(_has_pref, key=_has_pref.get)
_code2name = (data_all[_ref_year].drop_duplicates('cat01_code')
                                 .set_index('cat01_code')['市町村'].astype(str))
for year, d in data_all.items():
    if _has_pref[year] < 0.5:
        filled = d['cat01_code'].map(_code2name)
        n_filled = int(filled.notna().sum())
        data_all[year] = d.assign(市町村=filled.fillna(d['市町村']))
        print(f"{year}: 市町村名に都道府県名を補完（{n_filled}/{len(d)}行、基準年{_ref_year}）")

print(data_all[2021].head())

# %%
#　2021から2023のデータを結合する
data_2021_2023 = pd.concat(data_all)
data_2021_2023

# %%
# 1993から2005と2006のデータを結合
# カラム名の統一
data_1993_2005 = data_1993_2005.rename(columns={"時間軸（年次）（長期累年）": "時間軸（年次）", "市町村（長期累年）": "市町村"})
data_2006 = data_2006.rename(columns={"市町村（2006年）": "市町村"})

df_past = [data_1993_2005, data_2006]
data_1993_2006 = pd.concat(df_past, ignore_index=True)

#%%
# 不要な列の削除
data_1993_2006 = data_1993_2006.drop(columns=["水陸稲種類", "cat02_code", "time_code"])

# %%
# 不要な列の削除
#data_2007_2020 = data_2007_2020.drop(columns=["cat01_code", "time_code"])

# %%
# 1993-2020まで結合
df =[data_1993_2006, data_2007_2020]
data_1993_2020 = pd.concat(df, ignore_index=True)

# %%
data_2021_2023 = data_2021_2023.drop(columns=["cat01_code", "cat02_code"])
df = [data_2021_2023, data_1993_2020]
df_yield = pd.concat(df, ignore_index=True)

# Bug fix: 2021年以降のAPIは市町村名が「茨城県_水戸市」とアンダースコア区切りで返るが、
#          1993-2020は「茨城県水戸市」と区切りなし。後段の city_list は1993-2020の
#          名称から作るため、正規化しないと2021-2023が1行も通らず
#          kanto_yield.csv が2020年止まりになっていた。
df_yield["市町村"] = (df_yield["市町村"].astype(str)
                          .str.replace("_", "", regex=False).str.strip())
#df_yield['時間軸（年次）'].unique()

# %%
#最終調整
df_yield = df_yield.drop(columns=["unit", "annotation", "area_code", "time_code"]) # 不要な列の削除
#収量だけ抽出する
df_yield_per_area = df_yield[df_yield['面積収量収穫量'].str.contains('収量', na=False)].copy()

# %%
# 時間軸（年次）の値から'年'を削除
# 最初の4桁の数字を抜き出す
df_yield_per_area["時間軸（年次）"] = df_yield_per_area["時間軸（年次）"].astype(str).str.extract(r'(\d+)')[0]
df_yield_per_area
#df_yield_per_area['時間軸（年次）'].unique()

# %%
# 関東地方分を抽出(茨城・栃木・群馬・埼玉・千葉・東京・神奈川)
# 関東地方の都道府県名
kanto_prefectures = ["茨城県", "栃木県", "群馬県", "埼玉県", "千葉県", "東京都", "神奈川県"]
pref_pattern = '|'.join(re.escape(pref) for pref in kanto_prefectures)

# 関東各都県の「都道府県名+市町村名」を抜き出す
city_names = data_1993_2020[data_1993_2020["市町村"].str.contains(pref_pattern, na=False)]['市町村']
city_names = city_names.str.strip().unique()
# 都道府県名だけの集計行(市町村が付かないもの)は除外
city_names = [name for name in city_names if name and name not in kanto_prefectures]
print(city_names)

#%%
#　関東地方の抽出(都道府県名+市町村名でマッチングし、県をまたいだ同名市町村の誤マッチを防ぐ)
city_list = '|'.join(re.escape(name) for name in city_names)
df_kanto_yield = df_yield_per_area[df_yield_per_area['市町村'].str.contains(city_list, na=False)].copy()
df_kanto_yield

#%%
# 都道府県名と市町村名に分割する
df_kanto_yield["prefecture"] = df_kanto_yield["市町村"].str.extract(f'({pref_pattern})')
df_kanto_yield["city_name"] = (df_kanto_yield["市町村"].str.replace(pref_pattern, "", regex=True)
                    .str.replace(r"[市町村]$", "", regex=True)
                    .str.replace("_", "")
                    .str.strip()
                    )
df_kanto_yield['時間軸（年次）'].unique()

#%%
# 出力と保存
df_kanto_yield.to_csv(YIELD_DIR / "kanto_yield.csv", index=False, encoding='utf-8-sig')



# %%

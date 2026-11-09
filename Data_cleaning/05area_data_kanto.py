#%%
# 作付面積（ha）の市町村×年パネルを作る。
#
# 04yield_data_kanto_debuged.py は 167 行目で
#     df_yield['面積収量収穫量'].str.contains('収量')
# と絞っているため、同じ API レスポンスに入っている「作付面積」が落ちていた
# （「収穫量」は部分文字列「収量」を含まないので、残るのは 10a 当たり収量の3系列だけ）。
# T17（2010年の作付面積反応）で面積が必要になったので、接合・正規化のロジックは
# 04 とそろえたまま、抽出だけ作付面積に差し替えた並行スクリプトを置く。
# 04 の出力 kanto_yield.csv は触らない（既存の結果の再現性を壊さないため）。
import os
import re
from pathlib import Path

import pandas as pd

DATA_DIR = Path(os.environ.get("YIELD_DATA_DIR", "."))
YIELD_DIR = DATA_DIR / "yield"
output_path = "api_data"

# 既存の CSV だけを使う（appid は不要）。無ければ 04 を先に走らせる。
for tag in ["1993-2005", "2006", "2007-2020", "2021", "2022", "2023"]:
    p = YIELD_DIR / f"{output_path}_{tag}.csv"
    if not p.exists():
        raise RuntimeError(f"{p} が無い。先に 04yield_data_kanto_debuged.py を実行すること")

#%%
# --- 1993-2020（長期累年・2006・2007-2020）---
data_1993_2005 = pd.read_csv(YIELD_DIR / f"{output_path}_1993-2005.csv", encoding="utf-8", skiprows=92)
data_2006 = pd.read_csv(YIELD_DIR / f"{output_path}_2006.csv", encoding="utf-8", skiprows=91)
data_2007_2020 = pd.read_csv(YIELD_DIR / f"{output_path}_2007-2020.csv", encoding="utf-8", skiprows=91)

data_1993_2005 = data_1993_2005.rename(columns={"時間軸（年次）（長期累年）": "時間軸（年次）",
                                                "市町村（長期累年）": "市町村"})
data_2006 = data_2006.rename(columns={"市町村（2006年）": "市町村"})

data_1993_2006 = pd.concat([data_1993_2005, data_2006], ignore_index=True)
data_1993_2006 = data_1993_2006.drop(columns=["水陸稲種類", "cat02_code", "time_code"])

data_1993_2020 = pd.concat([data_1993_2006, data_2007_2020], ignore_index=True)

#%%
# --- 2021-2023（年ごとに別表。列名が年で違う）---
data_all = {}
for year in [2021, 2022, 2023]:
    df = pd.read_csv(YIELD_DIR / f"{output_path}_{year}.csv", encoding="utf-8", skiprows=91)
    new_cols = []
    for col in df.columns:
        if "全国" in col:
            new_cols.append("市町村")
        elif "北海道以外" in col:
            new_cols.append("面積収量収穫量")
        else:
            new_cols.append(col)
    df.columns = new_cols
    df["時間軸（年次）"] = year
    data_all[year] = df

# 04 と同じ補完：2022年産以降は市町村名に都道府県名が付かないので、
# 年をまたいで安定な cat01_code 経由で補う。
_has_pref = {y: d["市町村"].astype(str).str.contains("_").mean() for y, d in data_all.items()}
_ref_year = max(_has_pref, key=_has_pref.get)
_code2name = (data_all[_ref_year].drop_duplicates("cat01_code")
                                 .set_index("cat01_code")["市町村"].astype(str))
for year, d in data_all.items():
    if _has_pref[year] < 0.5:
        filled = d["cat01_code"].map(_code2name)
        print(f"{year}: 市町村名に都道府県名を補完（{int(filled.notna().sum())}/{len(d)}行、基準年{_ref_year}）")
        data_all[year] = d.assign(市町村=filled.fillna(d["市町村"]))

data_2021_2023 = pd.concat(data_all).drop(columns=["cat01_code", "cat02_code"])

#%%
# --- 全年結合と正規化 ---
df_area = pd.concat([data_2021_2023, data_1993_2020], ignore_index=True)
# 2021年以降は「茨城県_水戸市」、1993-2020は「茨城県水戸市」
df_area["市町村"] = df_area["市町村"].astype(str).str.replace("_", "", regex=False).str.strip()
df_area = df_area.drop(columns=["unit", "annotation", "area_code", "time_code"], errors="ignore")

# ここが 04 との唯一の違い：作付面積だけを取る
df_area = df_area[df_area["面積収量収穫量"].str.contains("作付面積", na=False)].copy()
df_area["時間軸（年次）"] = df_area["時間軸（年次）"].astype(str).str.extract(r"(\d+)")[0]

#%%
# --- 関東の抽出（04 と同じ手順・同じ city_list の作り方）---
kanto_prefectures = ["茨城県", "栃木県", "群馬県", "埼玉県", "千葉県", "東京都", "神奈川県"]
pref_pattern = "|".join(re.escape(p) for p in kanto_prefectures)

city_names = data_1993_2020[data_1993_2020["市町村"].str.contains(pref_pattern, na=False)]["市町村"]
city_names = city_names.str.strip().unique()
city_names = [n for n in city_names if n and n not in kanto_prefectures]

city_list = "|".join(re.escape(n) for n in city_names)
df_kanto_area = df_area[df_area["市町村"].str.contains(city_list, na=False)].copy()

df_kanto_area["prefecture"] = df_kanto_area["市町村"].str.extract(f"({pref_pattern})")
df_kanto_area["city_name"] = (df_kanto_area["市町村"].str.replace(pref_pattern, "", regex=True)
                              .str.replace(r"[市町村]$", "", regex=True)
                              .str.replace("_", "")
                              .str.strip())

#%%
# --- 診断 ---
v = pd.to_numeric(df_kanto_area["value"], errors="coerce")
print(f"行数 {len(df_kanto_area)}  年 {df_kanto_area['時間軸（年次）'].min()}-{df_kanto_area['時間軸（年次）'].max()}")
print(f"数値化できない value: {int(v.isna().sum())}行")
if v.isna().any():
    print(df_kanto_area.loc[v.isna(), "value"].value_counts().to_string())
print(df_kanto_area.groupby("時間軸（年次）").size().to_string())

df_kanto_area.to_csv(YIELD_DIR / "kanto_area.csv", index=False, encoding="utf-8-sig")
print("書き出し:", YIELD_DIR / "kanto_area.csv")

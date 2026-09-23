"""
作物統計調査「水稲の耕種期日（都道府県別）」のパース。

calender_data/ に置かれた e-Stat の統計表を、
都道府県 × 年 × 作業区分（は種期/田植期/出穂期/刈取期）の
tidy な通日(DOY)テーブルに変換する。

これまで出穂日は GDD 積算による推定値だったが、この表は
**観測された出穂期（最盛期）**そのものなので、
  (1) GDD 推定の精度検証
  (2) 移植日・GDD閾値の再較正
  (3) 観測値が使える年は観測値をそのまま使う
が可能になる。

注意: 都道府県単位なので市町村内の変動は与えない。
      市町村差は引き続き GDD 側が担う。

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/07parse_crop_calendar.py
"""
# %%
from pathlib import Path
import re
import numpy as np
import pandas as pd

CAL_DIR = Path("./calender_data/")
OUT = Path("./results/crop_calendar_pref.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

KANTO = ["茨城", "栃木", "群馬", "埼玉", "千葉", "東京", "神奈川"]
PREF_CODE = {"茨城": 8, "栃木": 9, "群馬": 10, "埼玉": 11,
             "千葉": 12, "東京": 13, "神奈川": 14}

# 各ファイルの構成。列の意味は表のヘッダ（5〜8行目）から読み取ったもの。
# stages: (作業区分, その区分の最初の年の列番号) を年数ぶん横に並べる。
# 各区分の直後に「対平年差」列が1本入るので、年数+1 ずつ進む。
FILES = {
    "f002-30-013.xls": {
        "years": [2014, 2015, 2016, 2017, 2018],   # 平成26〜30年産
        "stages": ["は種期", "田植期", "出穂期", "刈取期"],
        "first_col": 2,
    },
    "f002-05-013.xls": {
        "years": [2019, 2020, 2021, 2022, 2023],   # 令和元〜5年産
        "stages": ["田植期", "出穂期", "刈取期"],
        "first_col": 2,
    },
    # 表番号087は「(3) 水稲の耕種期日 ウ 出穂期」で、出穂期のみを収録する別建ての表。
    # 013系と違って「対平年差」列を持たないが、パーサは区分ごとに
    # 年数+1 列進める作りなので、区分が1つしかないこの表では影響しない。
    "f002-25-087.xls": {
        "years": [2009, 2010, 2011, 2012, 2013],   # 平成21〜25年産
        "stages": ["出穂期"],
        "first_col": 2,
        "has_deviation_col": False,
    },
    # 表番号086は「(3) 水稲の耕種期日」で、は種期(平成6・7年産)と出穂期(平成16〜20年産)が
    # 1枚に同居する変則レイアウト。県名が1列目、出穂期は5列目から。は種期は使わない。
    "f002-20-086.xls": {
        "years": [2004, 2005, 2006, 2007, 2008],   # 平成16〜20年産
        "stages": ["出穂期"],
        "first_col": 5,
        "name_col": 1,
        "has_deviation_col": False,
    },
    # 表番号135は出穂期のみ、平成12〜16年産。2004年は086と重複するので整合性チェックに使える。
    "f002-16-135.xls": {
        "years": [2000, 2001, 2002, 2003, 2004],   # 平成12〜16年産
        "stages": ["出穂期"],
        "first_col": 2,
        "has_deviation_col": False,
    },
}


# %%
def parse_md(value, year):
    """'8. 4' / '7.29' / '10. 3' 形式の月日を、その年の通日(DOY)に変換する。"""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    # 表によっては全角ピリオド（8．8）や全角数字が混ざる
    s = str(value).replace("　", " ")
    s = s.replace("．", ".").replace("｡", ".").strip()
    s = s.translate(str.maketrans("０１２３４５６７８９", "0123456789"))
    m = re.match(r"^(\d{1,2})\s*\.\s*(\d{1,2})$", s)
    if not m:
        return np.nan
    month, day = int(m.group(1)), int(m.group(2))
    if not (1 <= month <= 12 and 1 <= day <= 31):
        return np.nan
    return pd.Timestamp(year=year, month=month, day=day).dayofyear


def parse_file(path, spec):
    raw = pd.read_excel(path, header=None)
    n_years = len(spec["years"])

    # 都道府県名が入っている行だけ拾う（全国・地域計・脚注を自然に除外できる）
    # 表によって県名の列位置が違う（013/087/135系は0列目、086系は1列目）
    names = raw[spec.get("name_col", 0)].astype(str).str.replace(r"[\s　]", "", regex=True)
    rows = raw[names.isin(KANTO)].copy()
    rows["pref"] = names[rows.index]

    records = []
    col = spec["first_col"]
    for stage in spec["stages"]:
        for k, year in enumerate(spec["years"]):
            for _, r in rows.iterrows():
                records.append({
                    "pref": r["pref"],
                    "year": year,
                    "stage": stage,
                    "doy": parse_md(r[col + k], year),
                })
        # 013系は各区分の直後に「対平年差」列が1本入る。出穂期のみの表には無い。
        col += n_years + (1 if spec.get("has_deviation_col", True) else 0)
    return pd.DataFrame(records)


# %%
def build():
    frames = [parse_file(CAL_DIR / f, spec) for f, spec in FILES.items()]
    d = pd.concat(frames, ignore_index=True)

    # 2004年(平成16年産)は 135 と 086 の両表に載る。独立した2表なので整合性チェックに使う。
    dup = d.duplicated(["pref", "year", "stage"], keep=False)
    if dup.any():
        g = (d[dup].groupby(["pref", "year", "stage"])["doy"]
                   .agg(["nunique", "min", "max"]).reset_index())
        bad = g[g["nunique"] > 1]
        n_cells = len(g)
        if len(bad):
            print("[警告] 重複セルで値が食い違う:")
            print(bad.to_string(index=False))
        else:
            print(f"重複セル {n_cells}件（2004年産）は全表で一致 — パーサの整合性OK")
        d = d.drop_duplicates(["pref", "year", "stage"])

    wide = (d.pivot_table(index=["pref", "year"], columns="stage", values="doy")
             .reset_index())
    wide.columns.name = None
    wide["pref_code"] = wide["pref"].map(PREF_CODE)
    wide["prefecture"] = wide["pref"] + ["県" if p != "東京" else "都" for p in wide["pref"]]

    # 登熟期間の実測長（出穂 -> 刈取）
    if "刈取期" in wide.columns and "出穂期" in wide.columns:
        wide["登熟日数"] = wide["刈取期"] - wide["出穂期"]
    if "田植期" in wide.columns and "出穂期" in wide.columns:
        wide["移植出穂日数"] = wide["出穂期"] - wide["田植期"]
    return wide.sort_values(["pref_code", "year"]).reset_index(drop=True)


# %%
def to_date(doy, year=2015):
    if pd.isna(doy):
        return ""
    return (pd.Timestamp(year=year, month=1, day=1)
            + pd.Timedelta(days=int(doy) - 1)).strftime("%m/%d")


if __name__ == "__main__":
    w = build()
    w.to_csv(OUT, index=False, encoding="utf-8-sig")
    print("書き出し完了: " + str(OUT) + "  " + str(w.shape))
    print("収録年: " + str(sorted(w["year"].unique())))
    print("収録県: " + str(sorted(w["pref"].unique())) + "\n")

    show = w.copy()
    for c in ["は種期", "田植期", "出穂期", "刈取期"]:
        if c in show.columns:
            show[c] = show[c].map(to_date)
    print(show[["pref", "year"] + [c for c in ["は種期", "田植期", "出穂期", "刈取期",
                                               "移植出穂日数", "登熟日数"]
                                   if c in show.columns]].to_string(index=False))

    print("\n===== 県別平均 =====")
    print(w.groupby("pref").agg(
        田植期=("田植期", "mean"), 出穂期=("出穂期", "mean"),
        出穂期SD=("出穂期", "std"), 刈取期=("刈取期", "mean"),
        移植出穂日数=("移植出穂日数", "mean"), 登熟日数=("登熟日数", "mean"),
    ).round(1).to_string())

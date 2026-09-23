"""
推定用パネルの構築（収量 × 水田率 × 圃場整備率 × 登熟期気象）。

これまで `research/model/kanto_consol.csv`（2009-2020）は対話的に作られており
構築スクリプトが無かった（analysis_review_2026-09.md §10「未整備」）。
期間を 2000-2023 に広げるにあたり、ここに一本化する。

入力:
  yield/kanto_yield.csv                     収量（1993-2023、04で作成）
  total_puddyratio/paddy_ratio_master.csv   水田率（1980-2025、03で作成）
  Consolidation_data/GA0001_2014_2010_*.xlsx 農林業センサス農業集落カード（2010年センサス）
  results/phenology_kanto_2000_2023.csv     登熟期気象（09で作成）

出力:
  research/model/kanto_panel_2000_2023.csv

注意: 圃場整備率(consol30/pipe)は2010年センサスの1時点のみで、時間不変として扱う。
      期間を2000年まで遡ると「2010年の整備率を2000年に当てる」ことになるため、
      整備率の順位が期間内で安定しているという仮定に依存する。

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/10build_panel.py
"""
# %%
import os
from pathlib import Path
import numpy as np
import pandas as pd

YIELD_PATH = Path("./yield/kanto_yield.csv")
PADDY_PATH = Path("./total_puddyratio/paddy_ratio_master.csv")
CONSOL_DIR = Path("./Consolidation_data/")
PHENO_PATH = Path(os.environ.get("PANEL_PHENO", "./results/phenology_kanto_2000_2023.csv"))
OUT = Path(os.environ.get("PANEL_OUT", "./research/model/kanto_panel_2000_2023.csv"))

YEAR_MIN = int(os.environ.get("PANEL_YEAR_MIN", 2000))
YEAR_MAX = int(os.environ.get("PANEL_YEAR_MAX", 2023))

PREF_CODE_MAP = {
    8: "茨城県", 9: "栃木県", 10: "群馬県", 11: "埼玉県",
    12: "千葉県", 13: "東京都", 14: "神奈川県",
}

# 農業集落カード「整備状況」の田の区画規模別面積。
# 「ほ区均平」は整備済、「その他」は未整備を指す。
SIZE_COLS = ["田（1ha以上）", "田（0．5～1．0ha）", "田（0．3～0．5haほ区均平）",
             "田（0．3～0．5haその他）", "田（0．2～0．3haほ区均平）", "田（0．2～0．3haその他）"]
CONSOLIDATED = SIZE_COLS[:3]      # 30a以上の大区画＝整備済とみなす


# %%
def load_consolidation():
    """農林業センサス農業集落カードから市町村別の圃場整備率を構築する。"""
    frames = []
    for path in sorted(CONSOL_DIR.glob("GA0001_2014_2010_*.xlsx")):
        frames.append(pd.read_excel(path))
    if not frames:
        raise FileNotFoundError("Consolidation_data/GA0001_2014_2010_*.xlsx が見つからない")
    d = pd.concat(frames, ignore_index=True)

    for c in ["PREF", "CITY", "KCITY", "RCOM"]:
        d[c] = d[c].astype(str).str.zfill({"PREF": 2, "CITY": 3, "KCITY": 2, "RCOM": 3}[c])

    # KCITY=='00' & RCOM=='000' が市町村計の行。CITY=='000' は都県計なので除く。
    m = d[(d.KCITY == "00") & (d.RCOM == "000") & (d.CITY != "000")].copy()

    for c in SIZE_COLS + ["用水_田_パイプライン", "用水_田_開水路"]:
        m[c] = pd.to_numeric(m[c], errors="coerce").fillna(0)

    m["den"] = m[SIZE_COLS].sum(axis=1)                     # 0.2ha以上の田面積
    m["consol30"] = m[CONSOLIDATED].sum(axis=1) / m["den"]  # 30a以上区画整備率
    pipe_den = m["用水_田_パイプライン"] + m["用水_田_開水路"]
    m["pipe"] = np.where(pipe_den > 0, m["用水_田_パイプライン"] / pipe_den, np.nan)
    m["city_id"] = (m["PREF"] + m["CITY"]).astype(int)

    m = m[m["den"] > 0]
    out = m[["city_id", "den", "consol30", "pipe"]].drop_duplicates("city_id")
    print(f"圃場整備率: {len(out)}市町村（2010年センサス）")
    return out


# %%
def load_yield():
    y = pd.read_csv(YIELD_PATH)
    y["value"] = pd.to_numeric(y["value"], errors="coerce")
    y = y.rename(columns={"時間軸（年次）": "year"})
    y["year"] = y["year"].astype(int)
    y = (y.groupby(["year", "prefecture", "city_name"])["value"]
          .sum(min_count=1).reset_index(name="yield"))
    return y[y["year"].between(YEAR_MIN, YEAR_MAX)]


def load_paddy():
    p = pd.read_csv(PADDY_PATH)
    p[["city_id", "city_name"]] = p["city_id"].str.split("_", n=1, expand=True)
    p["city_id"] = p["city_id"].astype(int)
    p["prefecture"] = (p["city_id"] // 1000).map(PREF_CODE_MAP)
    return p[p["year"].between(YEAR_MIN, YEAR_MAX)][
        ["year", "city_id", "city_name", "prefecture", "paddy_ratio"]]


def load_pheno():
    p = pd.read_csv(PHENO_PATH)
    p = p.rename(columns={"city_id": "city_key"})
    p["city_id"] = p["city_code"]
    drop = [c for c in ("city_key", "city_code", "pref_code", "prefecture", "city_name")
            if c in p.columns]
    return p.drop(columns=drop)


# %%
if __name__ == "__main__":
    consol = load_consolidation()
    y = load_yield()
    paddy = load_paddy()
    pheno = load_pheno()

    # 収量は (year, prefecture, city_name)、水田率は city_id を持つのでこれで橋渡しする
    d = paddy.merge(y, on=["year", "prefecture", "city_name"], how="inner")
    print(f"収量×水田率: {len(d)}行, {d.city_id.nunique()}市町村")

    d = d.merge(pheno, on=["city_id", "year"], how="inner")
    print(f"＋登熟期気象: {len(d)}行, {d.city_id.nunique()}市町村")

    d = d.merge(consol, on="city_id", how="inner")
    print(f"＋圃場整備率: {len(d)}行, {d.city_id.nunique()}市町村")

    d = d.dropna(subset=["yield"])
    d["log_yield"] = np.log(d["yield"])
    d["paddy_fix"] = d.groupby("city_id")["paddy_ratio"].transform("mean")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT, index=False, encoding="utf-8-sig")
    print(f"\n書き出し完了: {OUT}  {d.shape}")
    print(f"期間 {d.year.min()}-{d.year.max()}, {d.city_id.nunique()}市町村")
    print("\n年別の市町村数:")
    print(d.groupby("year")["city_id"].nunique().to_string())

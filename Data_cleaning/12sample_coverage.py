# =============================================================================
# 標本定義の根拠を再現する
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   python Data_cleaning/12sample_coverage.py
#
# 「関東6県」は設計上の前提ではなく、consol30 の分母が立つ市町村だけが
# 残った結果である。それを都県別に数値で示す。
#
# 出力:
#   research/model/output/sample_den_by_pref.csv    センサスの den の分布
#   research/model/output/sample_paddy_coverage.csv 田面積のカバー率
# =============================================================================

from pathlib import Path
import pandas as pd

PREF_NAME = {8: "茨城", 9: "栃木", 10: "群馬", 11: "埼玉", 12: "千葉",
             13: "東京", 14: "神奈川"}
# consol30 の分母 den を構成する6列（analysis_review_2026-09.md 第6.2節）
SZ = ["田（1ha以上）", "田（0．5～1．0ha）", "田（0．3～0．5haほ区均平）",
      "田（0．3～0．5haその他）", "田（0．2～0．3haほ区均平）", "田（0．2～0．3haその他）"]
OUT_DIR = Path("research/model/output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. センサスの den（0.2ha 以上の田面積）の分布
# -----------------------------------------------------------------------------
want = {s.strip() for s in SZ}
rows = []
for code, name in PREF_NAME.items():
    f = Path("Consolidation_data") / ("GA0001_2014_2010_%02d.xlsx" % code)
    d = pd.read_excel(f, dtype=str)
    cols = [c for c in d.columns if c.strip() in want]
    assert len(cols) == 6, (f, len(cols))
    m = d[(d.KCITY == "00") & (d.RCOM == "000") & (d.CITY != "000")].copy()
    for c in cols:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m["den"] = m[cols].sum(axis=1)
    rows.append(dict(pref=code, pref_name=name, census_cities=len(m),
                     den_positive=int((m.den > 0).sum()),
                     den_zero=int((m.den <= 0).sum()),
                     den_total=int(m.den.sum())))
den = pd.DataFrame(rows)
den.to_csv(OUT_DIR / "sample_den_by_pref.csv", index=False, encoding="utf-8-sig")

print("=== consol30 の分母 den（0.2ha 以上の田面積）の分布 ===")
print(den.to_string(index=False))
print("東京: 全 %d 市町村で den = 0 なので consol30 は 0/0 で定義できない"
      % int(den.loc[den.pref == 13, "census_cities"].iloc[0]))

# -----------------------------------------------------------------------------
# 2. 土地利用メッシュ由来の田面積（センサスとは独立）でカバー率を見る
# -----------------------------------------------------------------------------
cen = pd.read_csv("results/city_centroids.csv")
cen["pref"] = cen.city_id // 1000
panel_ids = set(pd.read_csv("research/model/kanto_panel_2000_2023.csv",
                            usecols=["city_id"]).city_id.unique())
cen["in_panel"] = cen.city_id.isin(panel_ids)

g = cen.groupby("pref").agg(mesh_cities=("city_id", "nunique"),
                            paddy_area=("paddy_area", "sum")).reset_index()
g2 = (cen[cen.in_panel].groupby("pref")
      .agg(panel_cities=("city_id", "nunique"), paddy_in_panel=("paddy_area", "sum"))
      .reset_index())
cov = g.merge(g2, on="pref", how="left").fillna({"panel_cities": 0, "paddy_in_panel": 0})
cov["pref_name"] = cov.pref.map(PREF_NAME)
cov["share_of_kanto_pct"] = (cov.paddy_area / cov.paddy_area.sum() * 100).round(2)
cov["covered_pct"] = (cov.paddy_in_panel / cov.paddy_area * 100).round(1)
cov["panel_cities"] = cov.panel_cities.astype(int)
cov = cov[["pref", "pref_name", "mesh_cities", "panel_cities", "paddy_area",
           "share_of_kanto_pct", "covered_pct"]]
cov.to_csv(OUT_DIR / "sample_paddy_coverage.csv", index=False, encoding="utf-8-sig")

print("\n=== 土地利用メッシュ由来の田面積とカバー率（2006年ラスタ）===")
print(cov.to_string(index=False))
tot = cen.paddy_area.sum()
inp = cen.loc[cen.in_panel, "paddy_area"].sum()
print("関東全体の田面積のうち推定パネルがカバーする割合: %.1f%%" % (inp / tot * 100))
print("保存:", OUT_DIR / "sample_den_by_pref.csv", "/", OUT_DIR / "sample_paddy_coverage.csv")

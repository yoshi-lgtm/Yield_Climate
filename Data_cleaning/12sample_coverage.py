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
# consol30 の分母 den を構成する6列（260910_analysis_review.md 第6.2節）
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

# -----------------------------------------------------------------------------
# 3. 3つの「田面積」の比較と、consol30 の分母カバー率
#
#    den_ha  : センサス「0.2ha 以上の区画の田」（2010年2月1日、ha）
#    mesh_ha : 土地利用メッシュの「田」（区画規模を問わない、ha 換算）
#    作付面積（作物統計）はその年に水稲を作付けした面積で、上の2つとは別物。
#    kanto_yield.csv には 10a 当たり収量しか入っていないので本スクリプトでは扱わない。
# -----------------------------------------------------------------------------
rows = []
for code, name in PREF_NAME.items():
    f = Path("Consolidation_data") / ("GA0001_2014_2010_%02d.xlsx" % code)
    d = pd.read_excel(f, dtype=str)
    cols = [c for c in d.columns if c.strip() in want]
    tot = d[(d.CITY == "000") & (d.KCITY == "00") & (d.RCOM == "000")]
    v = sum(pd.to_numeric(tot[c], errors="coerce").fillna(0).sum() for c in cols)
    rows.append(dict(pref=code, pref_name=name, den_ha=round(v, 1)))
meas = pd.DataFrame(rows)
mesh = cen.groupby("pref")["paddy_area"].sum().reset_index()
mesh["mesh_ha"] = (mesh.paddy_area / 10000).round(0)
meas = meas.merge(mesh[["pref", "mesh_ha"]], on="pref")
meas["ge02ha_share_pct"] = (meas.den_ha / meas.mesh_ha * 100).round(1)
meas.to_csv(OUT_DIR / "sample_area_measures.csv", index=False, encoding="utf-8-sig")

print("\n=== 3つの田面積の比較（都県計）===")
print(meas.to_string(index=False))
print("ge02ha_share_pct = 田のうち 0.2ha 以上の区画が占める割合")

# 市町村レベル: consol30 の分母が田全体のどれだけを捉えているか
pan = (pd.read_csv("research/model/kanto_panel_2000_2023.csv",
                   usecols=["city_id", "consol30", "den"])
       .drop_duplicates("city_id"))
mm = pan.merge(cen[["city_id", "paddy_area"]], on="city_id")
mm["mesh_ha"] = mm.paddy_area / 10000
mm = mm[(mm.mesh_ha > 0) & mm.den.notna() & mm.consol30.notna()].copy()
mm["cover"] = mm.den / mm.mesh_ha
mm["pref"] = mm.city_id // 1000

q = mm.cover.quantile([.1, .25, .5, .75, .9]).round(3)
print("\n=== 市町村レベル: 分母カバー率 cover = den / メッシュの田面積（n = %d）===" % len(mm))
print("10%% %.3f / 25%% %.3f / 中央 %.3f / 75%% %.3f / 90%% %.3f"
      % (q[.1], q[.25], q[.5], q[.75], q[.9]))
print("cover > 1 の件数（センサス > メッシュ。ラスタ解像度と土地利用版の差）: %d"
      % int((mm.cover > 1).sum()))
print("corr(consol30, cover): Pearson %.3f / Spearman %.3f"
      % (mm.consol30.corr(mm.cover), mm.consol30.corr(mm.cover, method="spearman")))
mm["cover_q"] = pd.qcut(mm.cover, 4, labels=["Q1", "Q2", "Q3", "Q4"])
byq = (mm.groupby("cover_q", observed=True)
       .agg(n=("consol30", "size"), consol30=("consol30", "mean"), cover=("cover", "mean"))
       .round(3).reset_index())
print(byq.to_string(index=False))
bypref = (mm.groupby("pref").agg(n=("consol30", "size"), cover=("cover", "mean"),
                                 consol30=("consol30", "mean")).round(3).reset_index())
bypref["pref_name"] = bypref.pref.map(PREF_NAME)
print(bypref[["pref_name", "n", "cover", "consol30"]].to_string(index=False))
mm[["city_id", "pref", "den", "mesh_ha", "cover", "consol30"]].to_csv(
    OUT_DIR / "sample_cover_by_city.csv", index=False, encoding="utf-8-sig")
print("保存:", OUT_DIR / "sample_area_measures.csv", "/", OUT_DIR / "sample_cover_by_city.csv")

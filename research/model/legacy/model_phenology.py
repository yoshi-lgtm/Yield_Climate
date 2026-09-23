"""
登熟期気象変数（GDD出穂日ベース）による収量回帰。

標本: 関東7都県 × 2009-2020（圃場整備率 consol30/pipe をマージできる期間）
仕様: ln(yield) を市町村FE + 年FE、city_id クラスターSE。

研究上の問い（research/docs/analysis_review_2026-09.md §9）:
  (1) 固定暦月の heat を登熟期ベースの指標に置き換えると within R^2 は上がるか
  (2) 上がった heat 指標のもとで heat × consol30 の交差項はまだゼロか

実行はリポジトリルート(Yield_Climate/)から:
    python research/model/legacy/model_phenology.py
"""
# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from twfe import feols, summary, stars  # noqa: E402

PANEL = Path("research/model/kanto_consol.csv")
PHENO = Path("results/phenology_kanto.csv")
OUT = Path("research/model/kanto_pheno.csv")


# %%
def build_panel():
    """収量・圃場整備率パネルに登熟期気象変数を結合する。"""
    d = pd.read_csv(PANEL)
    p = pd.read_csv(PHENO)

    p = p.rename(columns={"city_id": "city_key"})
    p["city_id"] = p["city_code"]

    drop = [c for c in ("prefecture", "city_name", "city_code", "city_key") if c in p.columns]
    m = d.merge(p.drop(columns=drop), on=["city_id", "year"], how="inner")

    print(f"パネル {len(d)}行 × 登熟期気象 {len(p)}行 -> 結合 {len(m)}行 "
          f"({m['city_id'].nunique()}市町村, {m['year'].min()}-{m['year'].max()})")

    m["log_yield"] = np.log(m["yield"])
    # paddy_ratio はメッシュvintage由来の見かけの within 変動を持つため市町村平均で固定
    m["paddy_fix"] = m.groupby("city_id")["paddy_ratio"].transform("mean")
    return m


# %%
def heat_horse_race(m):
    """heat 指標の定義を入れ替えて within R^2 を比較する（review §9-2 の判定基準）。"""
    specs = [
        ("現行: 7-8月 HD34（固定暦月）", ["heat_old"], ["GSR_78", "APCP_78"]),
        ("7-8月 平均気温", ["T_78"], ["GSR_78", "APCP_78"]),
        ("登熟期(年次) 平均気温 T_rip", ["r20_T_rip"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) HDD26", ["r20_HDD26"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) HDD27", ["r20_HDD27"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) HD35", ["r20_HD35"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) HD34", ["r20_HD34"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) 夜温26日数", ["r20_NIGHT26"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(年次) HDD26+CDD20", ["r20_HDD26", "r20_CDD20"], ["r20_GSR_rip", "r20_APCP_rip"]),
        ("登熟期(平年窓) 平均気温", ["cl_r20_T_rip"], ["cl_r20_GSR_rip", "cl_r20_APCP_rip"]),
        ("登熟期(平年窓) HDD26", ["cl_r20_HDD26"], ["cl_r20_GSR_rip", "cl_r20_APCP_rip"]),
        ("登熟期(平年窓) HDD27", ["cl_r20_HDD27"], ["cl_r20_GSR_rip", "cl_r20_APCP_rip"]),
        ("登熟期(平年窓) HD35", ["cl_r20_HD35"], ["cl_r20_GSR_rip", "cl_r20_APCP_rip"]),
        ("登熟40日(年次) HDD26", ["r40_HDD26"], ["r40_GSR_rip", "r40_APCP_rip"]),
        ("開花期(年次) HD35", ["fl_HD35"], ["fl_GSR_rip", "fl_APCP_rip"]),
    ]

    rows = []
    for label, heat_vars, ctrl in specs:
        cols = heat_vars + ctrl
        if any(c not in m.columns or m[c].isna().all() for c in cols):
            print(f"[skip] {label}: 変数が未構築")
            continue
        r = feols(m, "log_yield", cols, ["city_id", "year"], "city_id")
        rows.append({
            "spec": label,
            "heat": heat_vars[0],
            "coef": r["coef"][0],
            "se": r["se"][0],
            "sig": stars(r["p"][0]),
            "within_r2": r["within_r2"],
            "n": r["n"],
        })
    return pd.DataFrame(rows).sort_values("within_r2", ascending=False)


# %%
def bin_regression(m, prefix="r20_"):
    """気温ビン推定（Schlenker & Roberts 型）。基準ビンは 22-24℃。"""
    bins = [c for c in m.columns if c.startswith(prefix + "bin_")]
    if not bins:
        print("[skip] 気温ビン: TMP_min 未構築のためスキップ")
        return None, None
    ref = prefix + "bin_22_24"
    use = [b for b in bins if b != ref]
    ctrl = [prefix + "GSR_rip", prefix + "APCP_rip"]
    r = feols(m, "log_yield", use + ctrl, ["city_id", "year"], "city_id")
    tab = pd.DataFrame({
        "bin": use,
        "coef": r["coef"][:len(use)],
        "se": r["se"][:len(use)],
        "sig": [stars(x) for x in r["p"][:len(use)]],
    })
    return tab, r


# %%
def interaction_horse_race(m, heat_var, ctrl):
    """heat × 圃場整備率 / 水田率 / パイプライン率 の horse race。"""
    d = m.copy()
    d["h_consol"] = d[heat_var] * d["consol30"]
    d["h_paddy"] = d[heat_var] * d["paddy_fix"]
    d["h_pipe"] = d[heat_var] * d["pipe"]

    results = {}
    results["heat単独"] = feols(d, "log_yield", [heat_var] + ctrl,
                               ["city_id", "year"], "city_id")
    results["+ consol30"] = feols(d, "log_yield", [heat_var, "h_consol"] + ctrl,
                                 ["city_id", "year"], "city_id")
    results["+ paddy_ratio"] = feols(d, "log_yield", [heat_var, "h_paddy"] + ctrl,
                                    ["city_id", "year"], "city_id")
    results["horse race"] = feols(d, "log_yield",
                                 [heat_var, "h_consol", "h_paddy", "h_pipe"] + ctrl,
                                 ["city_id", "year"], "city_id")
    return results


# %%
if __name__ == "__main__":
    m = build_panel()
    m.to_csv(OUT, index=False, encoding="utf-8-sig")
    print("推定用パネルを書き出し: " + str(OUT) + "\n")

    print("###############  出穂日の記述統計  ###############")
    desc = m.groupby("prefecture").agg(
        heading_mean=("heading_doy", "mean"),
        heading_sd=("heading_doy", "std"),
        T_rip=("r20_T_rip", "mean"),
        HDD26=("r20_HDD26", "mean"),
        heat_old=("heat_old", "mean"),
        n=("year", "size"),
    ).round(2)
    print(desc)
    print("\n出穂日の within SD（同一市町村の年次変動）: "
          f"{m.groupby('city_id')['heading_doy'].transform(lambda s: s - s.mean()).std():.2f} 日")
    print("出穂日の between SD（市町村間の差）: "
          f"{m.groupby('city_id')['heading_doy'].mean().std():.2f} 日\n")

    print("###############  heat 指標の horse race（within R2）  ###############")
    hr = heat_horse_race(m)
    with pd.option_context("display.width", 160, "display.max_colwidth", 40):
        print(hr.to_string(index=False, float_format=lambda x: f"{x:.5f}"))

    best = hr.iloc[0]
    print("\n最良指標: " + best["spec"] + f"  (within R2 = {best['within_r2']:.4f})")

    print("\n###############  気温ビン推定  ###############")
    tab, rbin = bin_regression(m)
    if tab is not None:
        print(tab.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
        print(f"within R2 = {rbin['within_r2']:.4f}, n = {rbin['n']}")

    print("\n###############  交差項の horse race  ###############")
    for hv, ctrl in [("r20_HDD26", ["r20_GSR_rip", "r20_APCP_rip"]),
                     ("cl_r20_HDD26", ["cl_r20_GSR_rip", "cl_r20_APCP_rip"]),
                     ("heat_old", ["GSR_78", "APCP_78"])]:
        if hv not in m.columns:
            continue
        print("\n----- heat = " + hv + " -----")
        for name, r in interaction_horse_race(m, hv, ctrl).items():
            print(summary(r, title=name))
            print()

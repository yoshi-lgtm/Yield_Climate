"""
heat 指標の定義比較（公正な条件での horse race）と交差項の再検証。

model_phenology.py の最初の horse race はスペックごとに
コントロール（GSR/APCP）の窓が違っており、within R^2 の差が
heat 指標の差なのかコントロールの差なのか分離できていなかった。
ここでは

  (1) heat 変数だけを入れ替え、コントロールは共通に固定する
  (2) コントロールのみの模型からの within R^2 の増分（incremental R^2）を見る
  (3) encompassing 検定: 固定暦月 heat と登熟期 heat を同時に入れる
  (4) 出穂日の仮定（GDD閾値・移植日）に対する感度
  (5) 山間部（平年出穂日が極端に遅い市町村）を落とした頑健性

実行はリポジトリルート(Yield_Climate/)から:
    python research/model/legacy/heat_definition_test.py
"""
# %%
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from twfe import feols, summary, stars  # noqa: E402

PANEL = Path("research/model/kanto_pheno.csv")
HIGHLAND_CUTOFF = 232  # 平年出穂日 DOY 232 = 8/20。これ以降は関東の水稲として非現実的


# %%
def load(trim_highland=False):
    m = pd.read_csv(PANEL)
    m["heading_clim"] = m.groupby("city_id")["heading_doy"].transform("mean")
    if trim_highland:
        keep = m["heading_clim"] < HIGHLAND_CUTOFF
        print(f"山間部トリム: {m['city_id'].nunique()}市町村 -> "
              f"{m.loc[keep, 'city_id'].nunique()}市町村 ({len(m)} -> {int(keep.sum())}行)")
        m = m[keep].copy()
    return m


# %%
HEAT_DEFS = {
    "現行 7-8月 HD34": "heat_old",
    "7-8月 平均気温": "T_78",
    "登熟20日(年次) HDD26": "r20_HDD26",
    "登熟20日(年次) HDD27": "r20_HDD27",
    "登熟20日(年次) HD35": "r20_HD35",
    "登熟20日(年次) 平均気温": "r20_T_rip",
    "登熟20日(年次) 夜温26日数": "r20_NIGHT26",
    "登熟40日(年次) HDD26": "r40_HDD26",
    "登熟40日(年次) HD35": "r40_HD35",
    "登熟20日(平年窓) HDD26": "cl_r20_HDD26",
    "登熟40日(平年窓) HDD26": "cl_r40_HDD26",
    "開花期(年次) HD35": "fl_HD35",
}

# 共通コントロール。どの heat 指標を使うときも同じものを入れる。
COMMON_CTRL = ["GSR_78", "APCP_78"]


def horse_race(m, ctrl=COMMON_CTRL, label=""):
    """コントロールを固定し heat 変数だけ入れ替えて incremental within R^2 を比較。"""
    base = feols(m, "log_yield", ctrl, ["city_id", "year"], "city_id")
    rows = [{"spec": "（コントロールのみ）", "var": "-", "coef": np.nan, "se": np.nan,
             "sig": "", "within_r2": base["within_r2"], "delta_r2": 0.0}]
    for name, v in HEAT_DEFS.items():
        if v not in m.columns or m[v].isna().all():
            continue
        r = feols(m, "log_yield", [v] + ctrl, ["city_id", "year"], "city_id")
        rows.append({"spec": name, "var": v, "coef": r["coef"][0], "se": r["se"][0],
                     "sig": stars(r["p"][0]), "within_r2": r["within_r2"],
                     "delta_r2": r["within_r2"] - base["within_r2"]})
    out = pd.DataFrame(rows)
    print("\n##### heat 指標 horse race" + label
          + "（コントロール共通: " + ", ".join(ctrl) + "） #####")
    print(out.sort_values("delta_r2", ascending=False)
             .to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    return out


def standardized_race(m, ctrl=COMMON_CTRL):
    """heat 指標を標準化して係数の大きさを比較可能にする。"""
    d = m.copy()
    rows = []
    for name, v in HEAT_DEFS.items():
        if v not in d.columns or d[v].isna().all():
            continue
        z = "z_" + v
        d[z] = (d[v] - d[v].mean()) / d[v].std()
        r = feols(d, "log_yield", [z] + ctrl, ["city_id", "year"], "city_id")
        rows.append({"spec": name, "coef_per_SD": r["coef"][0], "se": r["se"][0],
                     "sig": stars(r["p"][0]), "t": r["t"][0]})
    out = pd.DataFrame(rows).sort_values("t")
    print("\n##### 1標準偏差あたりの ln(yield) 変化 #####")
    print(out.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    return out


def encompassing(m):
    """固定暦月 heat と登熟期 heat のどちらが情報を包含しているか。"""
    print("\n##### encompassing 検定 #####")
    for pair in [("heat_old", "r20_HDD26"), ("heat_old", "r40_HDD26"),
                 ("heat_old", "cl_r20_HDD26")]:
        a, b = pair
        r = feols(m, "log_yield", [a, b] + COMMON_CTRL, ["city_id", "year"], "city_id")
        print(summary(r, title=a + " vs " + b))
        print()


def bins(m, prefix="r20_", ref="bin_22_24"):
    cols = [c for c in m.columns if c.startswith(prefix + "bin_")]
    if not cols:
        return None
    use = [c for c in cols if c != prefix + ref]
    r = feols(m, "log_yield", use + COMMON_CTRL, ["city_id", "year"], "city_id")
    tab = pd.DataFrame({"bin": [c.replace(prefix + "bin_", "") for c in use],
                        "coef": r["coef"][:len(use)], "se": r["se"][:len(use)],
                        "sig": [stars(x) for x in r["p"][:len(use)]],
                        "mean_days": [m[c].mean() for c in use]})
    print("\n##### 気温ビン推定 (" + prefix + ", 基準 " + ref + ") #####")
    print(tab.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print(f"within R2 = {r['within_r2']:.4f}, n = {r['n']}, cities = {r['n_cities']}")
    return tab, r


def interactions(m, heat_var, ctrl=COMMON_CTRL, label=""):
    """heat × consol30 / paddy_ratio / pipe の horse race。"""
    d = m.copy()
    d["paddy_fix"] = d.groupby("city_id")["paddy_ratio"].transform("mean")
    d["h_consol"] = d[heat_var] * d["consol30"]
    d["h_paddy"] = d[heat_var] * d["paddy_fix"]
    d["h_pipe"] = d[heat_var] * d["pipe"]

    print("\n##### 交差項 horse race: heat = " + heat_var + label + " #####")
    rows = []
    specs = {
        "heat単独": [heat_var],
        "+consol30": [heat_var, "h_consol"],
        "+paddy": [heat_var, "h_paddy"],
        "+pipe": [heat_var, "h_pipe"],
        "3つ同時": [heat_var, "h_consol", "h_paddy", "h_pipe"],
    }
    for name, X in specs.items():
        r = feols(d, "log_yield", X + ctrl, ["city_id", "year"], "city_id")
        row = {"spec": name}
        for i, v in enumerate(X):
            row[v] = f"{r['coef'][i]:.5f}{stars(r['p'][i])}"
            row[v + "_se"] = f"({r['se'][i]:.5f})"
        row["within_r2"] = round(r["within_r2"], 4)
        rows.append(row)
    out = pd.DataFrame(rows)
    print(out.to_string(index=False))

    # 平均 heat で評価した consol30 0->1 の効果と95%CI
    r = feols(d, "log_yield", [heat_var, "h_consol"] + ctrl, ["city_id", "year"], "city_id")
    hbar = d[heat_var].mean()
    est, se = r["coef"][1] * hbar, r["se"][1] * hbar
    print(f"\n平均 heat (= {hbar:.2f}) で評価した consol30 0->100% の収量効果: "
          f"{100*est:+.2f}%  95%CI [{100*(est-1.96*se):+.2f}%, {100*(est+1.96*se):+.2f}%]")
    return out


# %%
def sensitivity_heading():
    """GDD閾値・移植日の仮定を振って、交差項と heat 係数がどれだけ動くか。"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("ph", "Data_cleaning/06phenology_heat.py")
    ph = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ph)

    daily = ph.load_daily(ph.YEARS)
    panel = pd.read_csv("research/model/kanto_consol.csv")
    panel["log_yield"] = np.log(panel["yield"])
    panel["paddy_fix"] = panel.groupby("city_id")["paddy_ratio"].transform("mean")

    print("\n##### 出穂日の仮定に対する感度 #####")
    rows = []
    for gdd in [800.0, 900.0, 1000.0, 1100.0]:
        for shift in [-10, 0, 10]:
            h = ph.estimate_heading(daily, gdd_threshold=gdd, shift_days=shift)
            row = {"GDD": gdd, "移植日shift": shift,
                   "平年出穂DOY": round(float(h["heading_doy"].mean()), 1)}
            for days, tag in ((ph.RIPEN_DAYS, "r20"), (ph.RIPEN_DAYS_LONG, "r40")):
                agg, _ = ph.window_aggregate(daily, h, 0, days - 1, tag + "_")
                agg["city_id"] = agg["city_id"].str.split("_", n=1).str[0].astype(int)
                d = panel.merge(agg, on=["city_id", "year"], how="inner")
                hv = tag + "_HDD26"
                d["h_consol"] = d[hv] * d["consol30"]
                # kanto_consol.csv の GSR/APCP はすでに 7-8月の平均/合計
                r = feols(d, "log_yield", [hv, "h_consol", "GSR", "APCP"],
                          ["city_id", "year"], "city_id")
                row[tag + " heat"] = f"{r['coef'][0]:.5f}{stars(r['p'][0])}"
                row[tag + " ×consol"] = f"{r['coef'][1]:.5f}{stars(r['p'][1])}"
                row[tag + " SE"] = f"({r['se'][1]:.5f})"
                row[tag + " R2"] = round(r["within_r2"], 4)
            rows.append(row)
    out = pd.DataFrame(rows)
    print(out.to_string(index=False))
    return out


# %%
if __name__ == "__main__":
    m = load()
    print(f"標本: {len(m)}行, {m['city_id'].nunique()}市町村, "
          f"{m['year'].min()}-{m['year'].max()}, {m['prefecture'].nunique()}都県")

    horse_race(m)
    standardized_race(m)
    encompassing(m)
    bins(m, "r20_")
    bins(m, "cl_r20_")

    for hv in ["heat_old", "r20_HDD26", "r40_HDD26", "cl_r20_HDD26"]:
        interactions(m, hv)

    mt = load(trim_highland=True)
    horse_race(mt, label="（山間部トリム）")
    for hv in ["heat_old", "r20_HDD26"]:
        interactions(mt, hv, label="（山間部トリム）")

    sensitivity_heading()

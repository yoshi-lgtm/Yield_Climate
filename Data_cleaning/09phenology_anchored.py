"""
実測出穂期にアンカーした出穂日と登熟期気象変数。

f002-25-087.xls（平成21〜25年産＝2009-2013の出穂期）の追加により、
**本パネル 2009-2020 の全12年で都県別の実測出穂期が利用可能**になった。

08phenology_calibrated.py は
  移植日（実測 or 平年外挿）→ 県別GDD閾値 → 出穂日
と推定していたが、推定値には都県×年で平均 -3.8 〜 +2.3 日の残差がある。
実測出穂期が全年で手に入った以上、残差を残す理由はない。

そこで都県×年ごとに
    offset = 実測出穂期 － GDD推定出穂日の都県平均
を求め、その都県×年の全市町村の出穂日を offset だけ平行移動する。

  - 都県×年の水準は実測に完全一致する
  - 市町村間の空間変動は GDD 由来のまま保たれる（offset は都県×年で共通）
  - 年次変動も実測のものになる

出力: results/phenology_kanto_anchored.csv

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/09phenology_anchored.py
"""
# %%
import importlib.util
import os
from pathlib import Path
import numpy as np
import pandas as pd


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ph = _load("ph", "Data_cleaning/06phenology_heat.py")
c8 = _load("c8", "Data_cleaning/08phenology_calibrated.py")

CAL_PATH = Path("./results/crop_calendar_pref.csv")
OUT = Path(os.environ.get("PHENO_OUT", "./results/phenology_kanto_anchored.csv"))
# 環境変数 PHENO_YEARS で期間を指定できる（例: "2000-2023"）
_y = os.environ.get("PHENO_YEARS", "2009-2020").split("-")
YEARS = range(int(_y[0]), int(_y[-1]) + 1)


# %%
def anchor_heading(heading, cal):
    """都県×年の実測出穂期に合わせて、市町村の出穂日を平行移動する。"""
    obs = cal[["pref_code", "year", "出穂期"]].dropna()
    est = (heading.dropna(subset=["heading_doy"])
                  .groupby(["pref_code", "year"], as_index=False)["heading_doy"]
                  .mean().rename(columns={"heading_doy": "gdd_pref_mean"}))
    off = obs.merge(est, on=["pref_code", "year"], how="inner")
    off["offset"] = off["出穂期"] - off["gdd_pref_mean"]

    print("=== 都県×年オフセット（実測 － GDD推定, 日）===")
    piv = off.pivot_table(index="pref_code", columns="year", values="offset").round(1)
    print(piv.to_string())
    print(f"\nオフセットの範囲: {off['offset'].min():.1f} 〜 {off['offset'].max():.1f} 日, "
          f"絶対値平均 {off['offset'].abs().mean():.2f} 日")

    h = heading.merge(off[["pref_code", "year", "offset", "出穂期"]],
                      on=["pref_code", "year"], how="left")
    h["heading_gdd"] = h["heading_doy"]
    h["heading_doy"] = (h["heading_doy"] + h["offset"]).round()
    h["heading_observed_pref"] = h["出穂期"]
    h["anchored"] = h["offset"].notna()
    return h.drop(columns=["出穂期"])


def ripening_length(cal):
    """都県別の実測登熟日数（出穂期→刈取期）。刈取期は2014年以降のみ収録。"""
    L = cal.dropna(subset=["登熟日数"]).groupby("pref_code")["登熟日数"].mean().round(0)
    print("\n=== 都県別 実測登熟日数（出穂→刈取, 2014-2023平均）===")
    print(L.to_string())
    return L.to_dict()


# %%
if __name__ == "__main__":
    daily = ph.load_daily(YEARS)
    cal_all = pd.read_csv(CAL_PATH)
    cal = cal_all[cal_all["year"].isin(YEARS)].copy()

    n_obs = cal["出穂期"].notna().sum()
    print(f"実測出穂期: {n_obs} / {len(cal)} 都県×年 "
          f"({sorted(cal.loc[cal['出穂期'].notna(), 'year'].unique())})\n")

    # 移植日とGDD閾値は 08 と同じ（田植期の実測は2014年以降のみ）
    cal_tp = cal_all[cal_all["田植期"].notna()]   # 実測田植期がある年（2014-2023）をすべて使う
    transplant = c8.transplant_table(cal_tp, YEARS)
    thresholds = c8.calibrate_gdd(daily, cal_tp)

    heading = c8.estimate_heading_calibrated(daily, transplant, thresholds)
    heading = anchor_heading(heading, cal)
    print(f"\n出穂日を確定できなかった city-year: "
          f"{int(heading['heading_doy'].isna().sum())} / {len(heading)}")

    lengths = ripening_length(cal_all)
    out = c8.build_windows(daily, heading, lengths)

    out[["city_code", "city_name"]] = out["city_id"].str.split("_", n=1, expand=True)
    out["city_code"] = out["city_code"].astype(int)
    out["prefecture"] = (out["city_code"] // 1000).map(ph.PREF_CODE_MAP)
    out.to_csv(OUT, index=False, encoding="utf-8-sig")
    print("\n書き出し完了: " + str(OUT) + "  " + str(out.shape))

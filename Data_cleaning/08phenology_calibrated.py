"""
観測された耕種期日で較正した出穂日推定と登熟期気象変数。

06phenology_heat.py は移植日を都県一律の仮定で置いていたが、
作物統計調査の実測（calender_data/, 2014-2023）と突き合わせた結果:

  - GDD閾値 1000 という較正はほぼ正しかった（実測の実効積算温度 = 996、変動係数 0.071）
  - 年次変動の再現も良好（within 相関 0.864、within SD 観測2.3日 vs 推定2.7日）
  - **移植日の仮定だけが外れていた**（群馬 -25日、東京 -27日、神奈川 -23日、埼玉 -7日）
    群馬・東京・神奈川は麦跡二毛作で6月中旬田植えであり、5月中旬という仮定が誤り

そこで
  (1) 移植日 = 観測された県×年の田植期（2014-2020）。それ以前は県別平年値で外挿
  (2) GDD閾値 = 県別に実測から較正（千葉956 〜 埼玉1152）
  (3) 登熟期の窓 = 実測の 出穂期→刈取期 日数（県別 43〜62日）
      ——固定の20日/40日ではなく、実際の登熟期間を使う

出穂日は市町村×年で推定されるため、県レベルでは実測に一致しつつ
市町村内の空間変動は保たれる。

出力: results/phenology_kanto_calibrated.csv

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/08phenology_calibrated.py
"""
# %%
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location("ph", "Data_cleaning/06phenology_heat.py")
ph = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ph)

CAL_PATH = Path("./results/crop_calendar_pref.csv")
OUT = Path("./results/phenology_kanto_calibrated.csv")
YEARS = range(2009, 2021)
CAL_YEARS = range(2014, 2021)   # 実測が本パネルと重なる年


# %%
def load_calendar():
    cal = pd.read_csv(CAL_PATH)
    return cal[cal["year"].isin(CAL_YEARS)].copy()


def transplant_table(cal, years):
    """県×年の移植日(DOY)。実測がある年は実測、無い年は県別平年値で外挿。"""
    obs = cal[["pref_code", "year", "田植期"]].dropna()
    clim = obs.groupby("pref_code")["田植期"].mean().round()
    sd = obs.groupby("pref_code")["田植期"].std().round(2)
    print("=== 観測田植期（県別）===")
    print(pd.DataFrame({"平年DOY": clim, "年次SD": sd}).to_string())

    rows = []
    for pc in clim.index:
        for y in years:
            hit = obs[(obs.pref_code == pc) & (obs.year == y)]
            rows.append({
                "pref_code": pc, "year": y,
                "transplant_doy": float(hit["田植期"].iloc[0]) if len(hit) else float(clim[pc]),
                "transplant_observed": len(hit) > 0,
            })
    return pd.DataFrame(rows)


def calibrate_gdd(daily, cal, t_base=ph.T_BASE):
    """観測田植期→観測出穂期 の実効積算温度を県別に求め、GDD閾値とする。"""
    d = daily.copy()
    d["pref_code"] = d["city_code"] // 1000
    # 県×日の平均気温（市町村単純平均）
    pref_daily = d.groupby(["pref_code", "year", "doy"], as_index=False)["TMP_mea"].mean()

    obs = cal[["pref_code", "year", "田植期", "出穂期"]].dropna()
    rows = []
    for _, r in obs.iterrows():
        s = pref_daily[(pref_daily.pref_code == r.pref_code)
                       & (pref_daily.year == r.year)
                       & (pref_daily.doy >= r.田植期)
                       & (pref_daily.doy <= r.出穂期)]
        rows.append({"pref_code": r.pref_code, "year": r.year,
                     "gdd": float(np.maximum(s["TMP_mea"] - t_base, 0).sum())})
    g = pd.DataFrame(rows)
    thr = g.groupby("pref_code")["gdd"].mean().round(0)
    print("\n=== 県別GDD閾値（実測から較正）===")
    print(pd.DataFrame({"閾値": thr,
                        "年次SD": g.groupby("pref_code")["gdd"].std().round(1)}).to_string())
    return thr.to_dict()


def ripening_length(cal):
    """県別の実測登熟日数（出穂期→刈取期）。"""
    L = cal.groupby("pref_code")["登熟日数"].mean().round(0)
    print("\n=== 県別 実測登熟日数（出穂→刈取）===")
    print(L.to_string())
    return L.to_dict()


# %%
def estimate_heading_calibrated(daily, transplant, thresholds, t_base=ph.T_BASE):
    """県×年の観測移植日と県別GDD閾値で、市町村×年の出穂日を推定する。"""
    d = daily[["city_id", "city_code", "year", "doy", "TMP_mea"]].copy()
    d["pref_code"] = d["city_code"] // 1000
    d = d.merge(transplant, on=["pref_code", "year"], how="left")
    d["gdd_threshold"] = d["pref_code"].map(thresholds)

    post = d[d["doy"] >= d["transplant_doy"]].copy()
    post["gdd_day"] = np.maximum(post["TMP_mea"] - t_base, 0.0)
    post["gdd_cum"] = post.groupby(["city_id", "year"])["gdd_day"].cumsum()

    reached = post[post["gdd_cum"] >= post["gdd_threshold"]]
    heading = (reached.groupby(["city_id", "year"], as_index=False)["doy"].min()
                      .rename(columns={"doy": "heading_doy"}))

    base = post.groupby(["city_id", "year"], as_index=False).agg(
        pref_code=("pref_code", "first"),
        transplant_doy=("transplant_doy", "first"),
        transplant_observed=("transplant_observed", "first"),
        gdd_threshold=("gdd_threshold", "first"),
    )
    return base.merge(heading, on=["city_id", "year"], how="left")


def validate(heading, cal):
    """推定出穂日を県平均に落として実測と比較する。"""
    h = heading.dropna(subset=["heading_doy"])
    e = h.groupby(["pref_code", "year"], as_index=False)["heading_doy"].mean()
    m = cal[["pref_code", "year", "出穂期"]].dropna().merge(e, on=["pref_code", "year"])
    m["err"] = m["heading_doy"] - m["出穂期"]

    print("\n=== 較正後の検証（推定出穂日 － 実測出穂期, 日）===")
    print(m.groupby("pref_code").agg(実測=("出穂期", "mean"), 推定=("heading_doy", "mean"),
                                     誤差平均=("err", "mean"), 誤差SD=("err", "std")
                                     ).round(2).to_string())
    print(f"\nMAE = {m['err'].abs().mean():.2f} 日,  RMSE = {np.sqrt((m['err']**2).mean()):.2f} 日")

    for c in ["出穂期", "heading_doy"]:
        m[c + "_d"] = m[c] - m.groupby("pref_code")[c].transform("mean")
    print(f"within 相関 = {m[['出穂期_d', 'heading_doy_d']].corr().iloc[0, 1]:.3f}"
          f"  (実測 within SD {m['出穂期_d'].std():.2f}日 / 推定 {m['heading_doy_d'].std():.2f}日)")
    return m


# %%
def build_windows(daily, heading, lengths):
    """登熟期の気象変数。窓は県別の実測登熟日数。比較用に20日/40日も作る。"""
    pieces = [heading]

    # 県別の実測登熟日数ぶんの窓。県ごとに長さが違うので県別に集計して縦結合する。
    per_pref = []
    per_pref_bins = []
    for pc, L in lengths.items():
        sub_head = heading[heading["pref_code"] == pc]
        if sub_head.empty:
            continue
        sub_daily = daily[daily["city_code"] // 1000 == pc]
        agg, w = ph.window_aggregate(sub_daily, sub_head, 0, int(L) - 1, "obs_")
        per_pref.append(agg)
        b = ph.window_bins(w, "obs_")
        if b is not None:
            per_pref_bins.append(b)
    if per_pref:
        pieces.append(pd.concat(per_pref, ignore_index=True))
    if per_pref_bins:
        pieces.append(pd.concat(per_pref_bins, ignore_index=True))

    # 比較用の固定長窓
    for days, tag in ((20, "c20_"), (40, "c40_"), (50, "c50_")):
        agg, w = ph.window_aggregate(daily, heading, 0, days - 1, tag)
        pieces.append(agg)
        if tag == "c50_":
            b = ph.window_bins(w, tag)
            if b is not None:
                pieces.append(b)

    pieces.append(ph.fixed_month_heat(daily))

    out = pieces[0]
    for p in pieces[1:]:
        out = out.merge(p, on=["city_id", "year"], how="left")
    return out


# %%
if __name__ == "__main__":
    daily = ph.load_daily(YEARS)
    cal = load_calendar()

    transplant = transplant_table(cal, YEARS)
    thresholds = calibrate_gdd(daily, cal)
    lengths = ripening_length(cal)

    heading = estimate_heading_calibrated(daily, transplant, thresholds)
    print(f"\n出穂日を推定できなかった city-year: "
          f"{int(heading['heading_doy'].isna().sum())} / {len(heading)}")
    validate(heading, cal)

    out = build_windows(daily, heading, lengths)
    out[["city_code", "city_name"]] = out["city_id"].str.split("_", n=1, expand=True)
    out["city_code"] = out["city_code"].astype(int)
    out["prefecture"] = (out["city_code"] // 1000).map(ph.PREF_CODE_MAP)
    out.to_csv(OUT, index=False, encoding="utf-8-sig")
    print("\n書き出し完了: " + str(OUT) + "  " + str(out.shape))

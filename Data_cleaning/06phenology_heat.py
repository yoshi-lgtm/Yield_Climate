"""
GDD積算による出穂日推定と、登熟期の気象変数の構築。

背景（research/docs/analysis_review_2026-09.md §7, §9-2）:
    現行の heat は「7-8月の固定暦月で TMP_max>=34 の日数」であり、
    出穂・登熟期の年次・地域変動を捉えていない（within R^2 = 0.025）。
    設計文書（research/slides/research_plan.md §4）の方針どおり、
    移植日から有効積算温度を積み上げて出穂日 t* を推定し、
    登熟クリティカル期間 [t*, t*+D] の気象変数を定義しなおす。

標本: 関東7都県、2009-2020（圃場整備率データ GA0001(2010年センサス) とマージ可能な期間）

出力: results/phenology_kanto.csv  (city_id × year の長期パネル)

実行はリポジトリルート(Yield_Climate/)から:
    python Data_cleaning/06phenology_heat.py
"""
# %%
from pathlib import Path
import numpy as np
import pandas as pd

# %%
# ---------------------------------------------------------------- 設定
YEARS = range(2009, 2021)
# NaN伝播バグ修正版（03c_rebuild_weather_nanfix.py）があればそちらを優先する。
# 修正版は5変数（TMP_min 込み）を1ファイルに持つので TMP_min の別途結合は不要。
weather_dir_fixed = Path("./avg_temp_fixed/")
weather_dir = Path("./avg_temp/")
tmp_min_path = Path("./avg_temp_all_years/city_weighted_tmp_min_2009_2020.csv")
output_dir = Path("./results/")
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "phenology_kanto.csv"

PREF_CODE_MAP = {
    8: "茨城県", 9: "栃木県", 10: "群馬県", 11: "埼玉県",
    12: "千葉県", 13: "東京都", 14: "神奈川県",
}

# 移植日（田植え）の仮定。都県コード -> 通日(DOY, 平年基準)
# 注意: これは各都県の一般的な水稲（うるち・普通期〜早期）の田植え期からの
#       外生的な仮定であり、観測値ではない。感度分析(SHIFT_DAYS)で頑健性を見る。
#       千葉は早期栽培が多く最も早い。群馬・埼玉は麦跡二毛作があり遅い。
TRANSPLANT_DOY = {
    8: 125,   # 茨城  5/5
    9: 130,   # 栃木  5/10
    10: 140,  # 群馬  5/20
    11: 135,  # 埼玉  5/15
    12: 115,  # 千葉  4/25
    13: 135,  # 東京  5/15
    14: 130,  # 神奈川 5/10
}

T_BASE = 10.0        # 有効積算温度の基準温度（設計文書は 10-12℃）
# 移植 -> 出穂 の有効積算温度閾値。
# 文献値をそのまま当てるのではなく、推定される県別平年出穂日が関東の
# コシヒカリの出穂期（千葉の早期 7月下旬 〜 群馬 8月上〜中旬）に整合するよう
# キャリブレーションした。1000 のとき7都県平均 7/31、
# 千葉 7/26 / 茨城 8/1 / 栃木 8/4 / 群馬 8/10。
# 感度分析は GDD_SENSITIVITY で行う。
GDD_HEADING = 1000.0
GDD_SENSITIVITY = [800.0, 900.0, 1000.0, 1100.0]
RIPEN_DAYS = 20      # 登熟クリティカル期間の長さ（出穂後20日）
RIPEN_DAYS_LONG = 40  # 登熟全期間（参考）
FLOWER_WINDOW = (-3, 5)  # 開花・受精期（障害不稔）の窓: t*-3 〜 t*+5


# %%
# ---------------------------------------------------------------- 読み込み
def load_daily(years):
    """市町村×日次の気象パネルを読み込む。修正版があればそちらを使う。"""
    years = list(years)
    use_fixed = all((weather_dir_fixed / ("city_weighted_avg_temp_" + str(y) + ".csv")).exists()
                    for y in years)
    src = weather_dir_fixed if use_fixed else weather_dir
    print(("NaN修正版を使用: " if use_fixed else "旧データを使用: ") + str(src))

    frames = []
    for year in years:
        frames.append(pd.read_csv(src / ("city_weighted_avg_temp_" + str(year) + ".csv")))
    df = pd.concat(frames, ignore_index=True)

    df[["city_code", "city_name"]] = df["city_id"].str.split("_", n=1, expand=True)
    df["city_code"] = df["city_code"].astype(int)

    if use_fixed:
        print("TMP_min は修正版に同梱: 欠測 "
              + str(int((df["TMP_min"].isna() & df["TMP_mea"].notna()).sum())) + "行")
    elif tmp_min_path.exists():
        dmin = pd.read_csv(tmp_min_path)
        # city_id の後半（市町村名）は shapefile 側が欠損だと "None"/"nan" と
        # 版によって表記が揺れるため、数値コードを結合キーにする
        dmin["city_code"] = dmin["city_id"].str.split("_", n=1).str[0].astype(int)
        dmin = dmin.drop(columns=["city_id"])
        df = df.merge(dmin, on=["time", "city_code"], how="left")
        n_na = int((df["TMP_min"].isna() & df["TMP_mea"].notna()).sum())
        print("TMP_min を結合: TMP_mea があるのに TMP_min が無い行 = " + str(n_na))
    else:
        print("[警告] TMP_min なし。気温ビンはスキップされる -> 先に 03b_add_tmp_min.py を実行")
        df["TMP_min"] = np.nan
    df["prefecture"] = (df["city_code"] // 1000).map(PREF_CODE_MAP)
    df["time"] = pd.to_datetime(df["time"])
    df["year"] = df["time"].dt.year
    df["doy"] = df["time"].dt.dayofyear

    # 全期間NaNの市町村（田ピクセルなし等）は落とす
    before = df["city_id"].nunique()
    ok = df.groupby("city_id")["TMP_mea"].transform(lambda s: s.notna().any())
    df = df[ok]
    print("気象データ欠測の市町村を除外: " + str(before) + " -> " + str(df["city_id"].nunique()))
    return df.sort_values(["city_id", "time"]).reset_index(drop=True)


# %%
# ---------------------------------------------------------------- 出穂日推定
def estimate_heading(df, gdd_threshold=GDD_HEADING, t_base=T_BASE, shift_days=0):
    """移植日からの有効積算温度が閾値に達した日(DOY)を出穂日として返す。

    戻り値: DataFrame[city_id, year, transplant_doy, heading_doy, heading_gdd_ok]
    """
    d = df[["city_id", "year", "doy", "TMP_mea", "city_code"]].copy()
    # city_code は市町村コード。移植日は都県単位の仮定なので上2桁に落として引く
    d["transplant_doy"] = (d["city_code"] // 1000).map(TRANSPLANT_DOY) + shift_days

    # 移植日以降のみ積算
    post = d[d["doy"] >= d["transplant_doy"]].copy()
    post["gdd_day"] = np.maximum(post["TMP_mea"] - t_base, 0.0)
    post["gdd_cum"] = post.groupby(["city_id", "year"])["gdd_day"].cumsum()

    reached = post[post["gdd_cum"] >= gdd_threshold]
    heading = (reached.groupby(["city_id", "year"], as_index=False)["doy"].min()
                      .rename(columns={"doy": "heading_doy"}))

    base = post.groupby(["city_id", "year"], as_index=False).agg(
        transplant_doy=("transplant_doy", "first"),
        gdd_season_total=("gdd_cum", "max"),
    )
    out = base.merge(heading, on=["city_id", "year"], how="left")
    out["heading_gdd_ok"] = out["heading_doy"].notna()
    return out


# %%
# ---------------------------------------------------------------- 気温ビン
# Schlenker & Roberts (2009) の日内正弦補間。
# 日最低 Tn・日最高 Tx の間を正弦カーブとみなし、閾値 c を上回る時間(日換算)を求める。
def _above_threshold_days(tn, tx, c):
    """1日のうち気温が c を超えていた時間（日単位, 0-1）と、その超過分の積分(degree-days)。"""
    tn = np.asarray(tn, dtype=float)
    tx = np.asarray(tx, dtype=float)
    M = (tx + tn) / 2.0
    W = (tx - tn) / 2.0

    frac = np.zeros_like(M)
    dd = np.zeros_like(M)

    # c <= tn : 終日超過
    full = c <= tn
    frac[full] = 1.0
    dd[full] = M[full] - c

    # tn < c < tx : 部分超過
    part = (~full) & (c < tx) & (W > 0)
    theta = np.arcsin(np.clip((c - M[part]) / W[part], -1.0, 1.0))
    frac[part] = 0.5 - theta / np.pi
    dd[part] = ((M[part] - c) * (0.5 - theta / np.pi)
                + W[part] * np.cos(theta) / np.pi)

    # c >= tx : ゼロのまま
    bad = ~np.isfinite(M) | ~np.isfinite(W)
    frac[bad] = np.nan
    dd[bad] = np.nan
    return frac, dd


def degree_days_in_bins(tn, tx, edges):
    """各ビン [edges[k], edges[k+1]) に滞在した時間（日単位）を返す。

    最上位は edges[-1] 以上の開区間。S&R と同じく
    「閾値超過時間」の差分としてビン滞在時間を構成する。
    """
    fracs = [_above_threshold_days(tn, tx, e)[0] for e in edges]
    cols = {}
    for k in range(len(edges) - 1):
        cols["bin_" + str(edges[k]) + "_" + str(edges[k + 1])] = fracs[k] - fracs[k + 1]
    cols["bin_" + str(edges[-1]) + "plus"] = fracs[-1]
    cols["bin_below_" + str(edges[0])] = 1.0 - fracs[0]
    return pd.DataFrame(cols)


BIN_EDGES = list(range(12, 40, 2))  # 12,14,...,38 -> 2℃刻み、38℃以上は開区間


# %%
# ---------------------------------------------------------------- 窓集計
def window_aggregate(df, heading, start_offset, end_offset, prefix):
    """出穂日 t* に対し [t*+start_offset, t*+end_offset] の気象変数を集計する。"""
    h = heading[["city_id", "year", "heading_doy"]].dropna(subset=["heading_doy"])
    d = df.merge(h, on=["city_id", "year"], how="inner")
    rel = d["doy"] - d["heading_doy"]
    w = d[(rel >= start_offset) & (rel <= end_offset)].copy()

    w["hdd26"] = np.maximum(w["TMP_mea"] - 26.0, 0.0)
    w["hdd27"] = np.maximum(w["TMP_mea"] - 27.0, 0.0)
    w["cdd20"] = np.maximum(20.0 - w["TMP_mea"], 0.0)
    w["hd33"] = (w["TMP_max"] >= 33.0).astype(float)
    w["hd34"] = (w["TMP_max"] >= 34.0).astype(float)
    w["hd35"] = (w["TMP_max"] >= 35.0).astype(float)
    w["night26"] = (w["TMP_min"] >= 26.0).astype(float)  # 高温登熟には夜温も効く

    agg = w.groupby(["city_id", "year"], as_index=False).agg(
        n_days=("doy", "size"),
        T_rip=("TMP_mea", "mean"),
        Tmax_rip=("TMP_max", "mean"),
        Tmin_rip=("TMP_min", "mean"),
        HDD26=("hdd26", "sum"),
        HDD27=("hdd27", "sum"),
        CDD20=("cdd20", "sum"),
        HD33=("hd33", "sum"),
        HD34=("hd34", "sum"),
        HD35=("hd35", "sum"),
        NIGHT26=("night26", "sum"),
        GSR_rip=("GSR", "mean"),
        APCP_rip=("APCP", "sum"),
    )
    agg = agg.rename(columns={c: prefix + c for c in agg.columns
                              if c not in ("city_id", "year")})
    return agg, w


def window_bins(w, prefix):
    """窓内の気温ビン（日単位の滞在時間合計）。TMP_min が無い場合は空を返す。"""
    if w["TMP_min"].isna().all():
        return None
    ww = w.dropna(subset=["TMP_min", "TMP_max"]).copy()
    bins = degree_days_in_bins(ww["TMP_min"].values, ww["TMP_max"].values, BIN_EDGES)
    bins.index = ww.index
    bins[["city_id", "year"]] = ww[["city_id", "year"]]
    out = bins.groupby(["city_id", "year"], as_index=False).sum()
    out = out.rename(columns={c: prefix + c for c in out.columns
                              if c not in ("city_id", "year")})
    return out


def fixed_month_heat(df):
    """比較対象: 現行定義（7-8月 TMP_max>=34 の日数）と 7-8月の気象。"""
    m78 = df[df["time"].dt.month.isin([7, 8])].copy()
    m78["hd34"] = (m78["TMP_max"] >= 34.0).astype(float)
    return m78.groupby(["city_id", "year"], as_index=False).agg(
        heat_old=("hd34", "sum"),
        T_78=("TMP_mea", "mean"),
        GSR_78=("GSR", "mean"),
        APCP_78=("APCP", "sum"),
    )


# %%
# ---------------------------------------------------------------- 実行
if __name__ == "__main__":
    daily = load_daily(YEARS)

    # --- 出穂日: 年次変動型 ---
    heading = estimate_heading(daily)
    print("\n=== 推定出穂日 (GDD=" + str(GDD_HEADING) + ", Tbase=" + str(T_BASE) + ") ===")
    print(heading["heading_doy"].describe().round(1))
    print("未到達(出穂日が推定できない city-year): "
          + str(int((~heading["heading_gdd_ok"]).sum())))
    print(pd.to_datetime(heading["heading_doy"].mean() - 1, unit="D",
                         origin="2015-01-01").strftime("平均出穂日 ≒ %m/%d"))

    # --- 出穂日: 平年値型（年次気象に依存しない外生的な窓） ---
    # 年次変動型は「その年の気象で窓が決まる」ため、窓自体が内生。例えば「冷夏だった年に出穂が遅れ、たまたまその遅れた窓の気温が高かった」のようになる。
    # 市町村ごとの平年出穂日を使う版を並行して構築し、両方で推定する。
    clim = (heading.groupby("city_id", as_index=False)["heading_doy"].mean()
                   .rename(columns={"heading_doy": "heading_doy_clim"}))
    clim["heading_doy_clim"] = clim["heading_doy_clim"].round()
    heading = heading.merge(clim, on="city_id", how="left")

    heading_clim = heading[["city_id", "year", "heading_doy_clim"]].rename(
        columns={"heading_doy_clim": "heading_doy"})

    # --- 窓集計 ---
    pieces = [heading.drop(columns=["heading_gdd_ok"])]

    for tag, head_df in (("", heading), ("cl_", heading_clim)):
        agg20, w20 = window_aggregate(daily, head_df, 0, RIPEN_DAYS - 1, tag + "r20_")
        agg40, w40 = window_aggregate(daily, head_df, 0, RIPEN_DAYS_LONG - 1, tag + "r40_")
        aggfl, _ = window_aggregate(daily, head_df, FLOWER_WINDOW[0], FLOWER_WINDOW[1],
                                    tag + "fl_")
        pieces += [agg20, agg40, aggfl]
        for w, sub in ((w20, "r20_"), (w40, "r40_")):
            b = window_bins(w, tag + sub)
            if b is not None:
                pieces.append(b)

    pieces.append(fixed_month_heat(daily))

    out = pieces[0]
    for p in pieces[1:]:
        out = out.merge(p, on=["city_id", "year"], how="left")

    out[["city_code", "city_name"]] = out["city_id"].str.split("_", n=1, expand=True)
    out["city_code"] = out["city_code"].astype(int)
    out["prefecture"] = (out["city_code"] // 1000).map(PREF_CODE_MAP)

    out.to_csv(output_path, index=False, encoding="utf-8-sig")
    print("\n書き出し完了: " + str(output_path) + "  " + str(out.shape))
    print(out[["heading_doy", "r20_T_rip", "r20_HDD26", "r20_HD35", "heat_old"]]
          .describe().round(2))

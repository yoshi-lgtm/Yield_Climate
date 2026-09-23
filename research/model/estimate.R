# =============================================================================
# 推定の正本: research/docs/panel_extension_2000_2023.md の全表を再現する
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/estimate.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv (Data_cleaning/10build_panel.py)
# 出力 : research/model/output/*.csv（各表）と標準出力
#
# 共通仕様
#   被説明変数 log_yield、市町村FE＋年FE、city_id クラスターロバストSE
#   コントロール GSR_78, APCP_78（7-8月日射・降水）を全モデル共通で固定
#   有意性の星は *** p<0.01, ** p<0.05, * p<0.1（本文の表と同じ）
#
# 標準誤差の自由度修正（SSC）
#   本文の数値は legacy/twfe.py で出したもので、固定効果パラメータを全数 K に数える
#   （fixest の ssc(K.fixef = "full") と同値）。fixest の既定は "nonnested"
#   （クラスターに入れ子の市町村FEを数えない）で、SE が約2%小さくなる。
#   また twfe.py は singleton を落とさないので、fixef.rm = "none" で揃える。
#   本文との一致を優先して既定を "full" にしている。切り替えは SSC_K を変える。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

PANEL   <- "Yield_Climate/research/model/kanto_panel_2000_2023.csv"
OUT_DIR <- "research/model/output"
CTRL    <- c("GSR_78", "APCP_78")
HEAT    <- "obs_HD35"      # 実測登熟期（出穂後 43–62日）の 35℃超日数
SSC_K   <- "full"          # "full" = 本文(twfe.py)と一致 / "nonnested" = fixest 既定

# 全推定で city_id クラスター
est <- function(f, dd) feols(f, dd, cluster = ~city_id, fixef.rm = "none",
                             ssc = ssc(K.fixef = SSC_K))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stars <- function(p) fifelse(p < 0.01, "***", fifelse(p < 0.05, "**", fifelse(p < 0.1, "*", "")))

# 係数表を data.table に。p 値は fixest の t(G-1) 分布
tidy <- function(m, vars = NULL) {
  ct <- coeftable(m)
  out <- data.table(var = rownames(ct), coef = ct[, 1], se = ct[, 2], t = ct[, 3], p = ct[, 4])
  if (!is.null(vars)) out <- out[var %in% vars]
  out[, sig := stars(p)]
  out[, `:=`(n = nobs(m), within_r2 = unname(r2(m, "wr2")))]
  out[]
}

fml <- function(rhs) as.formula(paste("log_yield ~", paste(c(rhs, CTRL), collapse = " + "),
                                      "| city_id + year"))

show <- function(title, dt, digits = 5) {
  cat("\n#####", title, "#####\n")
  num <- names(dt)[vapply(dt, is.double, logical(1))]
  p <- copy(dt)
  for (v in num) set(p, j = v, value = round(p[[v]], if (v == "within_r2") 4 else digits))
  print(p, row.names = FALSE)
}

save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

# -----------------------------------------------------------------------------
# データ
# -----------------------------------------------------------------------------
d <- fread(PANEL, encoding = "UTF-8")
d[, h_consol := get(HEAT) * consol30]
cat(sprintf("パネル: %d行, %d市町村, %d–%d\n", nrow(d), uniqueN(d$city_id),
            min(d$year), max(d$year)))

# -----------------------------------------------------------------------------
# §3 heat 指標の horse race（コントロール共通、ΔR² はコントロールのみ模型との差）
# -----------------------------------------------------------------------------
HEAT_DEFS <- c(
  "アンカー40日 HD35"          = "c40_HD35",
  "アンカー20日 HD35"          = "c20_HD35",
  "実測登熟期（43–62日）HD35"  = "obs_HD35",
  "アンカー40日 HDD26"         = "c40_HDD26",
  "実測登熟期 HDD26"           = "obs_HDD26",
  "現行 7-8月 HD34"            = "heat_old",
  "実測登熟期 平均気温"        = "obs_T_rip"
)

# 全指標で同じ標本を使う（heat 指標のどれかが欠ける 4 行を落とす）
d_hr <- na.omit(d, cols = c("log_yield", CTRL, HEAT_DEFS))
base <- est(fml(character(0)), d_hr)
base_r2 <- unname(r2(base, "wr2"))

hr <- rbindlist(lapply(names(HEAT_DEFS), function(nm) {
  v <- HEAT_DEFS[[nm]]
  tidy(est(fml(v), d_hr), v)[, spec := nm]
}))
hr[, delta_r2 := within_r2 - base_r2]
hr <- hr[order(-delta_r2), .(spec, var, coef, se, sig, within_r2, delta_r2, n)]
cat(sprintf("\nコントロールのみ within R2 = %.4f (n = %d)\n", base_r2, nobs(base)))
show("§3 horse race", hr)
save(hr, "s3_horse_race")

# -----------------------------------------------------------------------------
# §4 交差項 heat × consol30
# -----------------------------------------------------------------------------
inter <- function(dd, label, extra = character(0)) {
  m <- est(fml(c(HEAT, "h_consol", extra)), dd)
  ct <- tidy(m)
  data.table(sample = label, n = nobs(m),
             heat = ct[var == HEAT, coef], heat_sig = ct[var == HEAT, sig],
             delta = ct[var == "h_consol", coef], delta_se = ct[var == "h_consol", se],
             delta_t = ct[var == "h_consol", t], delta_sig = ct[var == "h_consol", sig])
}

# §4.1 leave-one-year-out
loyo <- rbindlist(lapply(sort(unique(d$year)), function(y)
  inter(d[year != y], as.character(y))))
setnames(loyo, "sample", "excluded_year")
show("§4.1 leave-one-year-out", loyo)
save(loyo, "s4_1_loyo")

# §4.2 頑健性
# バランス標本: 24年すべてに観測がある市町村
n_years  <- uniqueN(d$year)
balanced <- d[, .N, by = city_id][N == n_years, city_id]
d_bal    <- d[city_id %in% balanced]

rob <- rbindlist(list(
  inter(d, "2000–2023 全期間"),
  inter(d, "＋登熟期低温 CDD20 を統制", extra = "obs_CDD20"),
  inter(d[year >= 2009], "2009–2023"),
  inter(d[year >= 2005], "2005–2023（合併後）"),
  inter(d_bal, sprintf("バランス%d市町村×%d年", length(balanced), n_years)),
  inter(d_bal[year != 2010], "バランス − 2010年"),
  inter(d[year != 2010], "2010年除外"),
  inter(d[!year %in% c(2010, 2023)], "2010・2023年除外"),
  inter(d[prefecture != "群馬県"], "群馬除外")
))
show("§4.2 頑健性", rob)
save(rob, "s4_2_robustness")

# -----------------------------------------------------------------------------
# §5 非定常性
# -----------------------------------------------------------------------------
# §5.1 年別平均（市町村単純平均）
ym <- d[, .(obs_HD35 = mean(obs_HD35, na.rm = TRUE), heat_old = mean(heat_old),
            yield = mean(yield), heading_doy = mean(heading_doy, na.rm = TRUE)), by = year][order(year)]
show("§5.1 年別平均", ym, digits = 1)
save(ym, "s5_1_year_means")

# §5.2 2010 vs 2023 の断面比較（両年に観測がある市町村、年内 HD35 中央値で二分）
common <- intersect(d[year == 2010 & !is.na(obs_HD35), city_id],
                    d[year == 2023 & !is.na(obs_HD35), city_id])
cs <- d[year %in% c(2010, 2023) & city_id %in% common][
  , high := obs_HD35 > median(obs_HD35), by = year][
  , .(n_cities = .N, hd35_median = median(obs_HD35),
      yield_low = mean(yield[!high]), yield_high = mean(yield[high])), by = year][
  , diff := yield_high - yield_low][order(year)]
cat(sprintf("\n共通市町村: %d\n", length(common)))
show("§5.2 断面比較", cs, digits = 1)
save(cs, "s5_2_cross_section")

# §5.3 年別 heat 感応度 β_t
m_bt <- est(as.formula(paste("log_yield ~ i(year,", HEAT, ") +",
                               paste(CTRL, collapse = " + "), "| city_id + year")), d)
bt <- tidy(m_bt)[grepl("^year::", var)][, year := as.integer(sub("year::(\\d+):.*", "\\1", var))]
bt <- bt[, .(year, beta = coef, se, sig)]
cat(sprintf("\nβ_t 模型: within R2 = %.4f, n = %d\n", r2(m_bt, "wr2"), nobs(m_bt)))
show("§5.3 年別 β_t", bt)
save(bt, "s5_3_beta_by_year")

PERIODS <- list(c(2000L, 2004L), c(2005L, 2009L), c(2010L, 2014L), c(2015L, 2019L), c(2020L, 2023L))
# 注: tidy() の結果に p 列があるので、期間の変数名は pr にする
bp <- rbindlist(lapply(PERIODS, function(pr)
  tidy(est(fml(HEAT), d[year %between% pr]), HEAT)[
    , .(period = sprintf("%d–%d", pr[1], pr[2]), beta = coef, se, sig, n, within_r2)]))
show("§5.3 期間別 β", bp)
save(bp, "s5_3_beta_by_period")

# §5.4 気温ビン（実測登熟期窓、基準 22–24℃）前後比較
BIN_REF <- "obs_bin_22_24"
bins <- setdiff(grep("^obs_bin_", names(d), value = TRUE), BIN_REF)
bin_est <- function(dd, label) {
  m <- est(fml(bins), dd)
  tidy(m, bins)[, .(bin = sub("obs_bin_", "", var), coef, se, sig,
                    mean_days = vapply(var, function(v) mean(dd[[v]], na.rm = TRUE), 0),
                    within_r2, n, period = label)]
}
bn <- rbind(bin_est(d[year <= 2011], "2000–2011"), bin_est(d[year >= 2012], "2012–2023"))
show("§5.4 気温ビン", bn)
save(bn, "s5_4_bins")

cat("\n出力:", normalizePath(OUT_DIR), "\n")

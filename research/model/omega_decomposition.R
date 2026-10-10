# =============================================================================
# T11: プール推定量の estimand（omega_t 分解）
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/omega_decomposition.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
# 出力 : research/model/output/omega_*.csv と標準出力
#
# 恒等式（有限標本で厳密）
#   プール式      : y = FE + X + b * r        （r = 関心のある回帰子）
#   年別交差式    : y = FE + X + sum_t b_t * r * 1[t]
#   r~ を r の (FE, X) への射影残差とすると FWL より
#       b_hat = sum_t omega_t * b_hat_t,   omega_t = (sum_i r~_it * r_it) / (sum_it r~_it^2)
#   r~ は年別交差式の説明変数の張る空間に属するので残差との直交性が効き、
#   plim ではなく有限標本で等号が成り立つ。sum_t omega_t = 1。
#
# 注意（本スクリプトで検証する点）
#   goals_and_todo.md G0 は omega_t を「残差化処置の年内分散」
#   sum_i r~_it^2 と書いたが、厳密な重みは交差積 sum_i r~_it * r_it である。
#   両者は年ごとには一致しない（全年の合計では一致する）。非負性も自明でない。
#   どちらが恒等式を再現するかを数値で確定させる。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

PANEL <- if (file.exists("research/model/kanto_panel_2000_2023.csv")) {
  "research/model/kanto_panel_2000_2023.csv"
} else {
  "Yield_Climate/research/model/kanto_panel_2000_2023.csv"
}
OUT_DIR <- "research/model/output"
CTRL    <- c("GSR_78", "APCP_78")
HEAT    <- "obs_HD35"
SSC_K   <- "full"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
est  <- function(f, dd) feols(f, dd, cluster = ~city_id, fixef.rm = "none",
                              ssc = ssc(K.fixef = SSC_K))
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d[, h_consol := get(HEAT) * consol30]

# -----------------------------------------------------------------------------
# 分解の中核。regr = 関心のある回帰子、extra = プール式に残す他の回帰子
# -----------------------------------------------------------------------------
decompose <- function(dd, regr, extra = character(0), label = "") {
  vars <- unique(c("log_yield", regr, extra, CTRL, "city_id", "year"))
  s <- dd[complete.cases(dd[, ..vars])]

  rhs_pool <- paste(c(regr, extra, CTRL), collapse = " + ")
  m_pool <- est(as.formula(paste("log_yield ~", rhs_pool, "| city_id + year")), s)
  b_pool <- coef(m_pool)[regr]

  rhs_int <- paste(c(paste0("i(year, ", regr, ")"), extra, CTRL), collapse = " + ")
  m_int <- est(as.formula(paste("log_yield ~", rhs_int, "| city_id + year")), s)
  ci <- coef(m_int)
  ci <- ci[grepl("^year::", names(ci))]
  yr <- vapply(strsplit(names(ci), ":", fixed = TRUE), function(z) z[3], character(1))
  bt <- data.table(year = as.integer(yr), b_t = as.numeric(ci))

  # r を (city FE, year FE, extra, CTRL) に射影した残差
  proj <- unique(c(regr, extra, CTRL))
  M  <- as.data.frame(demean(s[, ..proj], s[, .(city_id, year)]))
  others <- setdiff(proj, regr)
  if (length(others)) {
    f <- as.formula(paste(regr, "~", paste(others, collapse = " + "), "- 1"))
    r_tilde <- as.numeric(residuals(lm(f, data = M)))
  } else {
    r_tilde <- as.numeric(M[[regr]])
  }

  s[, rtilde := r_tilde]
  s[, rraw   := get(regr)]
  ss_total <- sum(s$rtilde^2)   # 素の numeric にする（den はパネルの既存列なので j 内で衝突する）
  w <- s[, .(omega_cross = sum(rtilde * rraw) / ss_total,   # 厳密な重み（交差積）
             omega_var   = sum(rtilde^2)     / ss_total),   # 年内分散ベースの重み
         by = year][order(year)]
  s[, c("rtilde", "rraw") := NULL]

  w <- merge(w, bt, by = "year", all.x = TRUE)
  rec_cross <- w[, sum(omega_cross * b_t)]
  rec_var   <- w[, sum(omega_var   * b_t)]

  cat(sprintf("\n===== %s : %s =====\n", label, regr))
  cat(sprintf("n = %d, 市町村 = %d\n", nrow(s), uniqueN(s$city_id)))
  cat(sprintf("プール係数 b_hat                     = %+.8f\n", b_pool))
  cat(sprintf("sum omega_cross                      = %+.10f  (理論値 1)\n", w[, sum(omega_cross)]))
  cat(sprintf("sum omega_var                        = %+.10f  (理論値 1)\n", w[, sum(omega_var)]))
  cat(sprintf("sum omega_cross * b_t                = %+.8f  差 = %+.2e\n",
              rec_cross, rec_cross - b_pool))
  cat(sprintf("sum omega_var   * b_t                = %+.8f  差 = %+.2e\n",
              rec_var, rec_var - b_pool))
  cat(sprintf("min/max omega_cross                  = %+.5f / %+.5f  (負の年数 %d)\n",
              w[, min(omega_cross)], w[, max(omega_cross)], w[omega_cross < 0, .N]))
  cat(sprintf("min/max omega_var                    = %+.5f / %+.5f\n",
              w[, min(omega_var)], w[, max(omega_var)]))
  cat(sprintf("上位3年の omega_cross 合計           = %.4f\n",
              w[order(-omega_cross)][1:3, sum(omega_cross)]))

  w[, `:=`(regr = regr, spec = label,
           contrib_cross = omega_cross * b_t,
           share_cross   = omega_cross / sum(omega_cross))]
  print(w[, .(year, omega_cross = round(omega_cross, 5), omega_var = round(omega_var, 5),
              b_t = round(b_t, 6), contrib = round(contrib_cross, 7))], row.names = FALSE)
  list(w = w[], b_pool = b_pool, rec_cross = rec_cross, rec_var = rec_var,
       m_pool = m_pool, m_int = m_int, n = nrow(s))
}

cat(sprintf("パネル: %d行, %d市町村, %d-%d\n", nrow(d), uniqueN(d$city_id),
            min(d$year), max(d$year)))

# -----------------------------------------------------------------------------
# (1) beta による検証。交差項なしの式 = panel_extension 第3節 horse race
#     対応するプール値は -0.00327（-0.00395 は h_consol を含む式の heat 主効果）
# -----------------------------------------------------------------------------
res_b <- decompose(d, HEAT, character(0), "(1) beta 検証（交差項なし）")

# -----------------------------------------------------------------------------
# (2) delta の分解。heat はプールのまま、h_consol だけ年別にする
# -----------------------------------------------------------------------------
res_d <- decompose(d, "h_consol", HEAT, "(2) delta 分解（heat はプール）")

# -----------------------------------------------------------------------------
# (3) heat 主効果の分解（h_consol をプールで残した式）= -0.00395 に対応
# -----------------------------------------------------------------------------
res_h <- decompose(d, HEAT, "h_consol", "(3) heat 主効果（h_consol をプール）")

out <- rbindlist(list(res_b$w, res_d$w, res_h$w), use.names = TRUE)
save(out, "omega_weights")

summ <- data.table(
  spec      = c(res_b$w$spec[1], res_d$w$spec[1], res_h$w$spec[1]),
  regr      = c("obs_HD35", "h_consol", "obs_HD35"),
  n         = c(res_b$n, res_d$n, res_h$n),
  b_pool    = c(res_b$b_pool, res_d$b_pool, res_h$b_pool),
  rec_cross = c(res_b$rec_cross, res_d$rec_cross, res_h$rec_cross),
  rec_var   = c(res_b$rec_var, res_d$rec_var, res_h$rec_var))
summ[, `:=`(err_cross = rec_cross - b_pool, err_var = rec_var - b_pool)]
cat("\n##### 恒等式の検証まとめ #####\n")
print(summ, row.names = FALSE)
save(summ, "omega_identity_check")

# -----------------------------------------------------------------------------
# (4) heat と h_consol を両方とも年別にした版
#     プール式の delta は「年別 delta_t の加重平均」だけでは書けず、
#     年別 beta_t からの汚染項が乗る（複数処置の contamination）
#        delta_pool = sum_t wD_t * delta_t + sum_t wDh_t * beta_t
#        wD_t  = sum_i Dtil_it * D_it  / sum Dtil^2
#        wDh_t = sum_i Dtil_it * h_it  / sum Dtil^2
# -----------------------------------------------------------------------------
vars4 <- c("log_yield", HEAT, "h_consol", CTRL, "city_id", "year")
s4 <- d[complete.cases(d[, ..vars4])]

m4_pool <- est(as.formula(paste("log_yield ~", HEAT, "+ h_consol +",
                                paste(CTRL, collapse = " + "), "| city_id + year")), s4)
delta_pool <- coef(m4_pool)["h_consol"]

m4_int <- est(as.formula(paste("log_yield ~ i(year,", HEAT, ") + i(year, h_consol) +",
                               paste(CTRL, collapse = " + "), "| city_id + year")), s4)
c4 <- coef(m4_int)
c4 <- c4[grepl("^year::", names(c4))]
parts <- strsplit(names(c4), ":", fixed = TRUE)
k4 <- data.table(year = as.integer(vapply(parts, function(z) z[3], character(1))),
                 v    = vapply(parts, function(z) z[4], character(1)),
                 b    = as.numeric(c4))
bt4 <- dcast(k4, year ~ v, value.var = "b")
setnames(bt4, c(HEAT, "h_consol"), c("beta_t", "delta_t"))

proj4 <- c("h_consol", HEAT, CTRL)
M4 <- as.data.frame(demean(s4[, ..proj4], s4[, .(city_id, year)]))
f4 <- as.formula(paste("h_consol ~", paste(c(HEAT, CTRL), collapse = " + "), "- 1"))
s4[, Dtil := as.numeric(residuals(lm(f4, data = M4)))]
s4[, Draw := h_consol]
s4[, hraw := get(HEAT)]
ss4 <- sum(s4$Dtil^2)
w4 <- s4[, .(wD = sum(Dtil * Draw) / ss4, wDh = sum(Dtil * hraw) / ss4), by = year][order(year)]
w4 <- merge(w4, bt4, by = "year")

main_part  <- w4[, sum(wD  * delta_t)]
contam     <- w4[, sum(wDh * beta_t)]

cat("\n===== (4) heat と h_consol を両方 年別にした版 =====\n")
cat(sprintf("n = %d\n", nrow(s4)))
cat(sprintf("プール delta_hat                     = %+.8f\n", delta_pool))
cat(sprintf("sum wD  (理論値 1)                   = %+.10f\n", w4[, sum(wD)]))
cat(sprintf("sum wDh (理論値 0)                   = %+.10f\n", w4[, sum(wDh)]))
cat(sprintf("主項   sum wD  * delta_t             = %+.8f\n", main_part))
cat(sprintf("汚染項 sum wDh * beta_t              = %+.8f\n", contam))
cat(sprintf("合計                                 = %+.8f  差 = %+.2e\n",
            main_part + contam, main_part + contam - delta_pool))
cat(sprintf("汚染項が占める割合                   = %.1f%%\n", 100 * contam / delta_pool))
print(w4[, .(year, wD = round(wD, 5), wDh = round(wDh, 5),
             beta_t = round(beta_t, 6), delta_t = round(delta_t, 6),
             main = round(wD * delta_t, 7), contam = round(wDh * beta_t, 7))], row.names = FALSE)
fwrite(w4, file.path(OUT_DIR, "omega_contamination.csv"), bom = TRUE)

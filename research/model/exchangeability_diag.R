# =============================================================================
# 年 permutation の前提（交換可能性）の診断
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/exchangeability_diag.R
#
# goals_and_todo.md 第3.1節 の脅威を数値で確定させる。
#   交換可能性 = 同時分布が並べ替えに対して不変
#   2次モーメントでの必要条件は「分散が全年で等しい」かつ「対相関が全年対で等しい」
#     脅威2（精度の不均一）  -> 対角のばらつき
#     脅威4（AR 型の相関）    -> 非対角が |t-s| に依存するか
#   等相関（共通ショック）は交換可能なので無害。問題は |t-s| 依存と分散の不均一。
#
# 対象は T11 スペック(4)（beta と delta を両方 年別）の delta_t。
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

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
d <- fread(PANEL, encoding = "UTF-8")
d[, h_consol := get(HEAT) * consol30]
vars <- c("log_yield", HEAT, "h_consol", CTRL, "city_id", "year")
s <- d[complete.cases(d[, ..vars])]

m <- feols(as.formula(paste("log_yield ~ i(year,", HEAT, ") + i(year, h_consol) +",
                            paste(CTRL, collapse = " + "), "| city_id + year")),
           s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))

cat(sprintf("n = %d, 市町村 = %d, 年 = %d\n", nobs(m), uniqueN(s$city_id), uniqueN(s$year)))

# -----------------------------------------------------------------------------
# delta_t ブロックの VCV を取り出す
# -----------------------------------------------------------------------------
V  <- vcov(m)
nm <- colnames(V)
parts <- strsplit(nm, ":", fixed = TRUE)
is_d  <- vapply(parts, function(z) length(z) >= 4 && z[1] == "year" && z[4] == "h_consol", logical(1))
yr_d  <- rep(NA_integer_, length(nm))
yr_d[is_d] <- as.integer(vapply(parts[is_d], function(z) z[3], character(1)))
ord <- order(yr_d[is_d])
idx <- which(is_d)[ord]
years <- yr_d[idx]
Vd <- V[idx, idx]
dt_hat <- coef(m)[idx]
SEd <- sqrt(diag(Vd))

cat(sprintf("delta_t の本数 = %d (%d-%d)\n", length(years), min(years), max(years)))

# -----------------------------------------------------------------------------
# 脅威2：分散（精度）の不均一
# -----------------------------------------------------------------------------
cat("\n##### 脅威2：SE(delta_t) の不均一 #####\n")
tab <- data.table(year = years, delta = as.numeric(dt_hat), se = as.numeric(SEd))
tab[, tstat := delta / se]
cat(sprintf("SE の最小 %.6f (%d年) / 最大 %.6f (%d年) / 比 %.1f倍\n",
            min(SEd), years[which.min(SEd)], max(SEd), years[which.max(SEd)],
            max(SEd) / min(SEd)))
cat(sprintf("SE の変動係数 = %.3f\n", sd(SEd) / mean(SEd)))
cat("交換可能性は分散が全年で等しいことを要求するので、この比が1から離れるほど生の係数の並べ替えは無効。\n")

# -----------------------------------------------------------------------------
# 脅威4：非対角が |t-s| に依存するか
# -----------------------------------------------------------------------------
R <- Vd / outer(SEd, SEd)
K <- length(years)
pr <- CJ(i = 1:K, j = 1:K)[i < j]
pr[, `:=`(lag = abs(years[j] - years[i]), rho = R[cbind(i, j)])]

cat("\n##### 脅威4：対相関の構造 #####\n")
cat(sprintf("非対角相関 %d 個: 平均 %+.4f, sd %.4f, 範囲 [%+.4f, %+.4f]\n",
            nrow(pr), mean(pr$rho), sd(pr$rho), min(pr$rho), max(pr$rho)))
bylag <- pr[, .(n = .N, mean_rho = mean(rho), sd_rho = sd(rho)), by = lag][order(lag)]
print(bylag[, .(lag, n, mean_rho = round(mean_rho, 4), sd_rho = round(sd_rho, 4))],
      row.names = FALSE)

fit <- lm(rho ~ lag, data = pr)
cat(sprintf("\nrho を lag に回帰: 切片 %+.5f, 傾き %+.6f (対のSEは独立でないので記述統計)\n",
            coef(fit)[1], coef(fit)[2]))
cat(sprintf("lag 1 の平均 %+.4f / lag 10+ の平均 %+.4f / 差 %+.4f\n",
            pr[lag == 1, mean(rho)], pr[lag >= 10, mean(rho)],
            pr[lag == 1, mean(rho)] - pr[lag >= 10, mean(rho)]))
cat(sprintf("等相関からの逸脱 sd(rho)/|mean(rho)| = %.2f\n", sd(pr$rho) / abs(mean(pr$rho))))

# -----------------------------------------------------------------------------
# 参考：delta_t 系列そのものの ACF（素朴に見る量）
# -----------------------------------------------------------------------------
cat("\n##### 参考：delta_t 系列の ACF #####\n")
a <- acf(tab$delta, lag.max = 5, plot = FALSE)
cat(paste(sprintf("lag%d %+.3f", 1:5, as.numeric(a$acf)[2:6]), collapse = "  "), "\n")
cat("24点なので推定は粗い。真の異質性と標本誤差が混ざるため解釈は補助的。\n")

# -----------------------------------------------------------------------------
# 残差の系列相関（相関の源）
# -----------------------------------------------------------------------------
s[, e := resid(m)]
setorder(s, city_id, year)
s[, e_lag := shift(e), by = city_id]
s[, yr_lag := shift(year), by = city_id]
cons <- s[!is.na(e_lag) & year - yr_lag == 1]
rho1 <- cons[, cor(e, e_lag)]
cat("\n##### 残差の系列相関（市町村内・連続年のみ） #####\n")
cat(sprintf("n = %d, AR(1) 相関 = %+.4f\n", nrow(cons), rho1))
cat(sprintf("市町村FE 下の見かけの負バイアスは概ね -1/T = %.3f\n", -1 / uniqueN(s$year)))

# -----------------------------------------------------------------------------
# studentize すると 2010年の順位は変わるか
# -----------------------------------------------------------------------------
cat("\n##### studentize の効果 #####\n")
tab[, rank_raw := frank(-delta)]
tab[, rank_t   := frank(-tstat)]
cat(sprintf("生の delta: 2010年は %d / %d 位（最大は %d年）\n",
            tab[year == 2010, rank_raw], K, tab[which.max(delta), year]))
cat(sprintf("studentized: 2010年は %d / %d 位（最大は %d年）\n",
            tab[year == 2010, rank_t], K, tab[which.max(tstat), year]))
print(tab[order(-tstat)][1:6, .(year, delta = round(delta, 6), se = round(se, 6),
                                tstat = round(tstat, 2))], row.names = FALSE)

fwrite(tab, file.path(OUT_DIR, "exchange_delta_t.csv"), bom = TRUE)
fwrite(bylag, file.path(OUT_DIR, "exchange_corr_by_lag.csv"), bom = TRUE)
fwrite(as.data.table(R), file.path(OUT_DIR, "exchange_corr_matrix.csv"), bom = TRUE)

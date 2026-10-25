# =============================================================================
# T15: 年 × consol30 の飽和（2010年の政策ショックへの対処）
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t15_policy_saturate.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_centroids.csv
# 出力 : research/model/output/t15_*.csv と標準出力
#
# 261004_goals_and_todo.md T15、261023_estimand_window_policy.md 第3節。
#
# 2010年産は米戸別所得補償モデル事業の実施年（Yield_Climate/subsidy_ricefarmer.pdf）。
# 年FE は共通ショックを吸収するので、脅威は 2010 × consol30 の差別的効果のみ。
# 現行式には年 × consol30 の水準項がなく、2010年内で cov_i(heat*c, c) が正なので
# 脱落した水準シフトは delta_2010 に直接載る。
#
# 対処：主軸の式に sum_t 1[t]*c_i を入れる。2010年に何が起きたかを知らなくても
# 「すべての年固有ショックが consol30 に線形に効く成分」を吸収する。
# さらに 1[t]*c_i の係数列そのものが政策期の差別的水準効果の推定値になるので、
# 対処と検定が同一の回帰になる。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

PANEL <- if (file.exists("research/model/kanto_panel_2000_2023.csv")) {
  "research/model/kanto_panel_2000_2023.csv"
} else { "Yield_Climate/research/model/kanto_panel_2000_2023.csv" }
CENTR <- if (file.exists("results/city_centroids.csv")) {
  "results/city_centroids.csv"
} else { "Yield_Climate/results/city_centroids.csv" }
OUT_DIR <- "research/model/output"
CTRL    <- c("GSR_78", "APCP_78")
HEAT    <- "obs_HD35"
CONLEY  <- 50   # T12 第2節: 使える範囲は 10〜50km

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d[, h_consol := get(HEAT) * consol30]
cen <- fread(CENTR, encoding = "UTF-8")[, .(city_id, lon, lat)]
d <- merge(d, cen, by = "city_id", all.x = TRUE)

vars <- c("log_yield", HEAT, "h_consol", "consol30", CTRL, "city_id", "year", "lon", "lat")
s <- d[complete.cases(d[, ..vars])]
s[, pref := city_id %/% 1000L]
years <- sort(unique(s$year))
REF <- years[1]
cat(sprintf("n = %d, 市町村 = %d, 都県 = %d, 年 = %d (%d-%d)\n",
            nrow(s), uniqueN(s$city_id), uniqueN(s$pref), length(years),
            min(years), max(years)))

# -----------------------------------------------------------------------------
# 0. 年 × consol30 の列を作る
#    consol30 は時間不変なので sum_t 1[t]*c_i = c_i が市町村FEと共線。
#    基準年(2000)を落として 23 本にする。
# -----------------------------------------------------------------------------
zc <- character(0)
for (y in years[-1]) {
  nmz <- sprintf("zc_%d", y)
  set(s, j = nmz, value = s$consol30 * as.integer(s$year == y))
  zc <- c(zc, nmz)
}
cat(sprintf("年 × consol30: %d 本（基準年 %d を除外）\n", length(zc), REF))

# -----------------------------------------------------------------------------
# 1. 式の定義（baseline = T11/T12 スペック、saturated = T15）
# -----------------------------------------------------------------------------
f_pool_base <- as.formula(sprintf("log_yield ~ %s + h_consol + %s | city_id + year",
                                  HEAT, paste(CTRL, collapse = " + ")))
f_sat_base  <- as.formula(sprintf("log_yield ~ i(year, %s) + i(year, h_consol) + %s | city_id + year",
                                  HEAT, paste(CTRL, collapse = " + ")))
f_pool_sat  <- as.formula(sprintf("log_yield ~ %s + h_consol + %s + %s | city_id + year",
                                  HEAT, paste(zc, collapse = " + "),
                                  paste(CTRL, collapse = " + ")))
f_sat_sat   <- as.formula(sprintf("log_yield ~ i(year, %s) + i(year, h_consol) + %s + %s | city_id + year",
                                  HEAT, paste(zc, collapse = " + "),
                                  paste(CTRL, collapse = " + ")))

# -----------------------------------------------------------------------------
# 2. 共通の処理：プール推定 → omega → 年別係数 → 恒等式 → 検定
# -----------------------------------------------------------------------------
pick_idx <- function(mm, v) {
  parts <- strsplit(names(coef(mm)), ":", fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) >= 4 && z[1] == "year" && z[4] == v, logical(1))
  idx <- which(ok)
  idx[order(as.integer(vapply(parts[idx], function(z) z[3], character(1))))]
}

wald_equal <- function(b, V, lab) {
  K <- length(b)
  R <- matrix(0, K - 1, K)
  for (i in 1:(K - 1)) { R[i, i] <- 1; R[i, i + 1] <- -1 }
  Rb <- R %*% b; RV <- R %*% V %*% t(R)
  sol <- tryCatch(solve(RV), error = function(e) NULL)
  if (is.null(sol)) return(data.table(set = lab, df = K - 1, stat = NA_real_, p = NA_real_))
  st <- as.numeric(t(Rb) %*% sol %*% Rb)
  data.table(set = lab, df = K - 1, stat = st, p = pchisq(st, K - 1, lower.tail = FALSE))
}

run_spec <- function(label, f_pool, f_sat, extra_proj) {
  cat(sprintf("\n==================== %s ====================\n", label))

  mp <- feols(f_pool, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
  dpool <- unname(coef(mp)["h_consol"])
  sepool <- unname(se(mp)["h_consol"])
  cat(sprintf("プール delta = %+.6f (SE %.6f, t %+.2f)\n", dpool, sepool, dpool / sepool))

  ms <- feols(f_sat, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
  id_d <- pick_idx(ms, "h_consol"); id_b <- pick_idx(ms, HEAT)
  stopifnot(length(id_d) == length(years), length(id_b) == length(years))
  bh <- coef(ms)

  # omega: D を他の全回帰子（FE・heat・extra_proj・CTRL）に射影した残差
  proj <- c("h_consol", HEAT, extra_proj, CTRL)
  M <- as.data.frame(demean(s[, ..proj], s[, .(city_id, year)]))
  f4 <- as.formula(paste("h_consol ~",
                         paste(c(HEAT, extra_proj, CTRL), collapse = " + "), "- 1"))
  Dtil <- as.numeric(residuals(lm(f4, data = M)))
  ss4 <- sum(Dtil^2)
  wt <- data.table(year = s$year, Dtil = Dtil, Draw = s$h_consol, hraw = s[[HEAT]])[
    , .(wD = sum(Dtil * Draw) / ss4, wDh = sum(Dtil * hraw) / ss4), by = year][order(year)]
  stopifnot(identical(wt$year, years))

  main   <- sum(wt$wD  * bh[id_d])
  contam <- sum(wt$wDh * bh[id_b])
  err    <- dpool - (main + contam)
  cat(sprintf("恒等式 : 主項 %+.6f + 汚染項 %+.6f = %+.6f （プールとの差 %.2e）\n",
              main, contam, main + contam, err))
  cat(sprintf("         sum(omega) = %.8f, 汚染項のプール値比 %+.1f%%\n",
              sum(wt$wD), 100 * contam / dpool))
  cat(sprintf("         omega が負の年 = %d / %d\n", sum(wt$wD < 0), length(years)))

  out <- list()
  for (vl in c("cluster_city", sprintf("conley_%dkm", CONLEY))) {
    V <- if (vl == "cluster_city") vcov(ms) else
      suppressWarnings(vcov_conley(ms, lat = ~lat, lon = ~lon, cutoff = CONLEY,
                                   ssc = ssc(K.fixef = "full")))
    Vd <- V[id_d, id_d]
    se_d <- sqrt(pmax(diag(Vd), 0))
    i10 <- which(years == 2010L)
    d10 <- bh[id_d][i10]; s10 <- se_d[i10]
    v_main <- as.numeric(t(wt$wD)  %*% Vd %*% wt$wD)
    v_cont <- as.numeric(t(wt$wDh) %*% V[id_b, id_b] %*% wt$wDh)
    wl <- rbindlist(list(
      wald_equal(bh[id_d], Vd, "all24"),
      wald_equal(bh[id_d][years != 2009L], Vd[years != 2009L, years != 2009L], "drop2009"),
      wald_equal(bh[id_d][abs(wt$wD) > 0.01], Vd[abs(wt$wD) > 0.01, abs(wt$wD) > 0.01],
                 sprintf("absWD_gt_0.01_n%d", sum(abs(wt$wD) > 0.01)))))
    cat(sprintf("\n-- %s --\n", vl))
    cat(sprintf("(a) delta_2010 = %+.6f, SE %.6f, t %+.2f, p %.4f\n",
                d10, s10, d10 / s10, 2 * pnorm(-abs(d10 / s10))))
    cat("(c) 同質性 Wald:\n")
    print(wl[, .(set, df, stat = round(stat, 1), p = signif(p, 3))], row.names = FALSE)
    cat(sprintf("    主項 %+.6f (SE %.6f, t %+.2f) / 汚染項 %+.6f (SE %.6f, t %+.2f)\n",
                main, sqrt(max(v_main, 0)), main / sqrt(max(v_main, 1e-300)),
                contam, sqrt(max(v_cont, 0)), contam / sqrt(max(v_cont, 1e-300))))
    out[[vl]] <- data.table(spec = label, vcov = vl, delta_pool = dpool,
      delta_2010 = d10, se_2010 = s10, t_2010 = d10 / s10,
      omega_2010 = wt$wD[i10], share_2010 = wt$wD[i10] * d10 / dpool,
      n_neg_omega = sum(wt$wD < 0),
      wald_all = wl[set == "all24", stat], wald_all_p = wl[set == "all24", p],
      wald_no09 = wl[set == "drop2009", stat], wald_no09_p = wl[set == "drop2009", p],
      main = main, se_main = sqrt(max(v_main, 0)),
      contam = contam, se_contam = sqrt(max(v_cont, 0)),
      contam_pct = 100 * contam / dpool, identity_err = err)
  }
  list(sum = rbindlist(out),
       byyear = data.table(spec = label, year = years, delta = bh[id_d],
                           se = sqrt(pmax(diag(vcov(ms)[id_d, id_d]), 0)),
                           beta = bh[id_b], wD = wt$wD, wDh = wt$wDh,
                           contrib = wt$wD * bh[id_d],
                           share = 100 * wt$wD * bh[id_d] / dpool),
       m_sat = ms)
}

r_base <- run_spec("baseline", f_pool_base, f_sat_base, character(0))
r_sat  <- run_spec("saturated", f_pool_sat,  f_sat_sat,  zc)

# -----------------------------------------------------------------------------
# 3. 政策チャネルそのもの：1[t]*consol30 の係数列
#    これが T15 の「対処と検定が同一の回帰」の部分。
#    2010-2012 が揃って正で持続していれば政策チャネルを測ったことになる。
# -----------------------------------------------------------------------------
cat("\n==================== 政策チャネル：1[t] x consol30 ====================\n")
cat(sprintf("基準年 %d。係数は %d 年に対する収量(対数)の差分。\n", REF, REF))
ms <- r_sat$m_sat
cf <- coef(ms); V <- vcov(ms)
iz <- match(zc, names(cf))
ztab <- data.table(year = years[-1], coef = unname(cf[iz]),
                   se = sqrt(pmax(diag(V)[iz], 0)))
ztab[, t := coef / se]
ztab[, policy_era := year %in% 2010:2013]
print(ztab[, .(year, coef = round(coef, 5), se = round(se, 5), t = round(t, 2),
               policy_era)], row.names = FALSE)

# 2010-2012（モデル事業＋本格実施）が揃ってゼロか
pe <- which(years[-1] %in% 2010:2012)
Rz <- matrix(0, length(pe), length(cf)); for (k in seq_along(pe)) Rz[k, iz[pe[k]]] <- 1
Rb <- Rz %*% cf; RV <- Rz %*% V %*% t(Rz)
st <- as.numeric(t(Rb) %*% solve(RV) %*% Rb)
cat(sprintf("\n2010-2012 が揃ってゼロ : chi2(%d) = %.2f, p = %.4g\n",
            length(pe), st, pchisq(st, length(pe), lower.tail = FALSE)))
cat(sprintf("政策期(2010-2013)平均 %+.5f / 非政策期平均 %+.5f\n",
            ztab[policy_era == TRUE, mean(coef)], ztab[policy_era == FALSE, mean(coef)]))
cat("政策説はステップ（2010-2012 に持続、本格実施の 2011-12 の方が大きい）を予言し、\n")
cat("高温説は 2010 年だけのスパイクを予言する。\n")

# -----------------------------------------------------------------------------
# 4. Ibragimov-Muller（飽和スペックの delta_2010）
# -----------------------------------------------------------------------------
cat("\n==================== Ibragimov-Muller (G = 6, 飽和) ====================\n")
prefs <- sort(unique(s$pref))
im_one <- function(p) {
  sp <- s[pref == p]
  zp <- zc[vapply(zc, function(z) var(sp[[z]], na.rm = TRUE) > 0, logical(1))]
  fp <- as.formula(sprintf("log_yield ~ i(year, %s) + i(year, h_consol) + %s %s | city_id + year",
                           HEAT, paste(CTRL, collapse = " + "),
                           if (length(zp)) paste("+", paste(zp, collapse = " + ")) else ""))
  tryCatch({
    mm <- feols(fp, sp, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
    cfp <- coef(mm); pp <- strsplit(names(cfp), ":", fixed = TRUE)
    ok <- vapply(pp, function(z) length(z) >= 4 && z[3] == "2010" && z[4] == "h_consol",
                 logical(1))
    if (!any(ok)) NA_real_ else unname(cfp[ok][1])
  }, error = function(e) NA_real_)
}
imt <- data.table(pref = prefs, d2010 = vapply(prefs, im_one, numeric(1)))
print(imt[, .(pref, d2010 = round(d2010, 6))], row.names = FALSE)
bb <- imt$d2010[is.finite(imt$d2010)]; G <- length(bb)
tt <- mean(bb) / (sd(bb) / sqrt(G))
cat(sprintf("IM: G = %d, 平均 %+.6f, sd %.6f, t = %+.2f, p = %.3f (臨界値 %.3f)\n",
            G, mean(bb), sd(bb), tt, 2 * pt(-abs(tt), G - 1), qt(0.975, G - 1)))

# -----------------------------------------------------------------------------
# 5. 書き出しと比較表
# -----------------------------------------------------------------------------
summ <- rbindlist(list(r_base$sum, r_sat$sum))
byy  <- rbindlist(list(r_base$byyear, r_sat$byyear))
save(summ, "t15_summary"); save(byy, "t15_by_year"); save(ztab, "t15_policy_level")
save(imt, "t15_im_by_pref")

cat("\n==================== baseline 対 saturated ====================\n")
print(summ[vcov == "cluster_city",
  .(spec, delta_pool = round(delta_pool, 6), delta_2010 = round(delta_2010, 6),
    t_2010 = round(t_2010, 2), omega_2010 = round(omega_2010, 4),
    share_2010 = round(100 * share_2010, 1), n_neg_omega,
    wald = round(wald_all, 1), wald_p = signif(wald_all_p, 3),
    contam_pct = round(contam_pct, 1), id_err = signif(identity_err, 2))],
  row.names = FALSE)
cat("\n出力: research/model/output/t15_{summary,by_year,policy_level,im_by_pref}.csv\n")

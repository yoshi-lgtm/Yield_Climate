# =============================================================================
# T12: Var(delta_hat) を3通りで構築する
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t12_vcov.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_centroids.csv (Data_cleaning/11build_city_centroids.py)
# 出力 : research/model/output/t12_*.csv と標準出力
#
# 261004_goals_and_todo.md 第3.3節。帰無仮説ごとに使える手法が違う。
#   (a) delta_2010 = 0       1制約  -> 市町村クラスター / Conley / IM
#   (c) delta_t が全年同一   23制約 -> 市町村クラスター / Conley
#                                      （クラスターVCVの階数は高々 G-1）
# あわせて T11 の主項・汚染項（固定ウェイトの線形結合）の SE も出す。
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
CUTOFFS <- c(10, 25, 50, 75, 100, 200)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d[, h_consol := get(HEAT) * consol30]
cen <- fread(CENTR, encoding = "UTF-8")[, .(city_id, lon, lat)]
d <- merge(d, cen, by = "city_id", all.x = TRUE)
stopifnot(!any(is.na(d$lon)))

vars <- c("log_yield", HEAT, "h_consol", CTRL, "city_id", "year", "lon", "lat")
s <- d[complete.cases(d[, ..vars])]
s[, pref := city_id %/% 1000L]
cat(sprintf("n = %d, 市町村 = %d, 都県 = %d\n", nrow(s), uniqueN(s$city_id), uniqueN(s$pref)))
print(s[, .(cities = uniqueN(city_id), rows = .N), by = pref][order(pref)], row.names = FALSE)

# -----------------------------------------------------------------------------
# 0. 距離の分布。Conley の cutoff は対象地域の広がりに縛られる
# -----------------------------------------------------------------------------
cc0 <- unique(s[, .(city_id, lon, lat)])[order(city_id)]
K0 <- nrow(cc0)
latm <- mean(cc0$lat)
xkm <- cc0$lon * 111.0 * cos(latm * pi / 180)
ykm <- cc0$lat * 111.0
Dm <- as.matrix(dist(cbind(xkm, ykm)))
pd <- Dm[upper.tri(Dm)]
cat(sprintf("\n##### 市町村間距離（%d市町村・%d対） #####\n", K0, length(pd)))
cat(sprintf("最大 %.0f km / 中央 %.0f km / 平均 %.0f km\n", max(pd), median(pd), mean(pd)))
dshare <- data.table(cutoff = CUTOFFS)
dshare[, share_pairs := vapply(cutoff, function(z) mean(pd <= z), numeric(1))]
print(dshare[, .(cutoff, share_pairs = round(share_pairs, 3))], row.names = FALSE)
cat("対象地域の直径が約250kmなので、cutoff が100kmを超えると大半の対が含まれ\n")
cat("実質「全体を1クラスター」に近づく。Conley VCV は階数不足になる。\n")
save(dshare, "t12_distance_share")

# -----------------------------------------------------------------------------
# 1. 本体の推定（T11 スペック(4)）と T11 の重み
# -----------------------------------------------------------------------------
FML <- as.formula(paste("log_yield ~ i(year,", HEAT, ") + i(year, h_consol) +",
                        paste(CTRL, collapse = " + "), "| city_id + year"))
m <- feols(FML, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))

nm <- names(coef(m)); parts <- strsplit(nm, ":", fixed = TRUE)
pick <- function(v) {
  ok <- vapply(parts, function(z) length(z) >= 4 && z[1] == "year" && z[4] == v, logical(1))
  which(ok)[order(as.integer(vapply(parts[ok], function(z) z[3], character(1))))]
}
id_d <- pick("h_consol"); id_b <- pick(HEAT)
years <- sort(as.integer(vapply(parts[id_d], function(z) z[3], character(1))))
bhat <- coef(m); stopifnot(length(id_d) == 24, length(id_b) == 24)

proj <- c("h_consol", HEAT, CTRL)
M <- as.data.frame(demean(s[, ..proj], s[, .(city_id, year)]))
f4 <- as.formula(paste("h_consol ~", paste(c(HEAT, CTRL), collapse = " + "), "- 1"))
s[, Dtil := as.numeric(residuals(lm(f4, data = M)))]
ss4 <- sum(s$Dtil^2); s[, Draw := h_consol]; s[, hraw := get(HEAT)]
wtab <- s[, .(wD = sum(Dtil * Draw) / ss4, wDh = sum(Dtil * hraw) / ss4), by = year][order(year)]
stopifnot(identical(wtab$year, years))

# -----------------------------------------------------------------------------
# 2. 同質性 Wald（任意の年集合について）と VCV 別のレポート
# -----------------------------------------------------------------------------
wald_equal <- function(b, V, lab) {
  K <- length(b)
  R <- matrix(0, K - 1, K)
  for (i in 1:(K - 1)) { R[i, i] <- 1; R[i, i + 1] <- -1 }
  Rb <- R %*% b; RV <- R %*% V %*% t(R)
  ev <- eigen(RV, symmetric = TRUE, only.values = TRUE)$values
  sol <- tryCatch(solve(RV), error = function(e) NULL)
  if (is.null(sol)) {
    return(data.table(set = lab, df = K - 1, stat = NA_real_, p = NA_real_, min_eig = min(ev)))
  }
  st <- as.numeric(t(Rb) %*% sol %*% Rb)
  data.table(set = lab, df = K - 1, stat = st,
             p = pchisq(st, K - 1, lower.tail = FALSE), min_eig = min(ev))
}

i2010 <- which(years == 2010L)
keep_no09 <- which(years != 2009L)
keep_w <- which(abs(wtab$wD) > 0.01)
res <- list(); wrows <- list(); walds <- list()

report <- function(label, V, note) {
  Vd <- V[id_d, id_d]
  ev <- eigen(Vd, symmetric = TRUE, only.values = TRUE)$values
  se_d <- sqrt(pmax(diag(Vd), 0))
  d2010 <- bhat[id_d][i2010]; s2010 <- se_d[i2010]
  wD <- wtab$wD; wDh <- wtab$wDh
  main <- sum(wD * bhat[id_d]); contam <- sum(wDh * bhat[id_b])
  v_main <- as.numeric(t(wD) %*% V[id_d, id_d] %*% wD)
  v_cont <- as.numeric(t(wDh) %*% V[id_b, id_b] %*% wDh)

  wl <- rbindlist(list(
    wald_equal(bhat[id_d], Vd, "all24"),
    wald_equal(bhat[id_d][keep_no09], Vd[keep_no09, keep_no09], "drop2009"),
    wald_equal(bhat[id_d][keep_w], Vd[keep_w, keep_w],
               sprintf("absWD_gt_0.01_n%d", length(keep_w)))))[, vcov := label]

  cat(sprintf("\n----- %s %s -----\n", label, note))
  cat(sprintf("delta ブロックの最小固有値 %.3e （負なら非正定値）\n", min(ev)))
  cat(sprintf("(a) delta_2010 = %+.6f, SE %.6f, t %+.2f, p %.4f\n",
              d2010, s2010, d2010 / s2010, 2 * pnorm(-abs(d2010 / s2010))))
  cat(sprintf("    SE(delta_t) 最小 %.6f / 最大 %.6f / 比 %.1f倍\n",
              min(se_d), max(se_d), max(se_d) / min(se_d)))
  cat("(c) 同質性 Wald:\n")
  print(wl[, .(set, df, stat = round(stat, 1), p = signif(p, 3),
               min_eig = signif(min_eig, 2))], row.names = FALSE)
  cat(sprintf("    主項 %+.6f (SE %.6f) / 汚染項 %+.6f (SE %.6f)\n",
              main, sqrt(max(v_main, 0)), contam, sqrt(max(v_cont, 0))))

  walds[[label]] <<- wl
  wrows[[label]] <<- data.table(vcov = label, year = years, delta = bhat[id_d],
                                se = se_d, beta = bhat[id_b], wD = wD, wDh = wDh)
  res[[label]] <<- data.table(vcov = label, min_eig_d = min(ev),
    delta_2010 = d2010, se_2010 = s2010, t_2010 = d2010 / s2010,
    se_min = min(se_d), se_max = max(se_d),
    wald_all = wl[set == "all24", stat], wald_all_p = wl[set == "all24", p],
    wald_no09 = wl[set == "drop2009", stat],
    main = main, se_main = sqrt(max(v_main, 0)),
    contam = contam, se_contam = sqrt(max(v_cont, 0)))
}

report("cluster_city", vcov(m), "")
for (cc in CUTOFFS) {
  V <- suppressWarnings(vcov_conley(m, lat = ~lat, lon = ~lon, cutoff = cc,
                                    ssc = ssc(K.fixef = "full")))
  report(sprintf("conley_%dkm", cc), V,
         sprintf("(cutoff 内の対 %.1f%%)", 100 * dshare[cutoff == cc, share_pairs]))
}

# -----------------------------------------------------------------------------
# 3. Ibragimov-Muller（都県ごとに推定して6個の推定値に t 検定）
# -----------------------------------------------------------------------------
cat("\n##### Ibragimov-Muller (G = 6) #####\n")
prefs <- sort(unique(s$pref))
FML_P <- as.formula(paste("log_yield ~", HEAT, "+ h_consol +",
                          paste(CTRL, collapse = " + "), "| city_id + year"))
im_one <- function(p, what) {
  sp <- s[pref == p]
  tryCatch({
    if (what == "pooled") {
      mm <- feols(FML_P, sp, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
      unname(coef(mm)["h_consol"])
    } else {
      mm <- feols(FML, sp, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
      cf <- coef(mm); pp <- strsplit(names(cf), ":", fixed = TRUE)
      ok <- vapply(pp, function(z) length(z) >= 4 && z[3] == "2010" && z[4] == "h_consol",
                   logical(1))
      if (!any(ok)) NA_real_ else unname(cf[ok][1])
    }
  }, error = function(e) NA_real_)
}
imtab <- data.table(
  pref = prefs,
  cities = vapply(prefs, function(p) uniqueN(s[pref == p]$city_id), integer(1)),
  n = vapply(prefs, function(p) nrow(s[pref == p]), integer(1)),
  pooled = vapply(prefs, im_one, numeric(1), what = "pooled"),
  d2010 = vapply(prefs, im_one, numeric(1), what = "2010"))
print(imtab[, .(pref, cities, n, pooled = round(pooled, 6), d2010 = round(d2010, 6))],
      row.names = FALSE)

im_test <- function(b, lab) {
  b <- b[is.finite(b)]; G <- length(b)
  tt <- mean(b) / (sd(b) / sqrt(G))
  data.table(target = lab, G = G, mean = mean(b), sd = sd(b), t = tt,
             p = 2 * pt(-abs(tt), df = G - 1), crit95 = qt(0.975, G - 1))
}
imres <- rbindlist(list(im_test(imtab$pooled, "pooled_delta"),
                        im_test(imtab$d2010, "delta_2010")))
print(imres[, .(target, G, mean = round(mean, 6), sd = round(sd, 6), t = round(t, 2),
                p = round(p, 4), crit95 = round(crit95, 3))], row.names = FALSE)
cat("IM は群内で漸近正規・群間独立を仮定する。神奈川は8市町村しかなく群内推定は粗い。\n")
save(imtab, "t12_im_by_pref"); save(imres, "t12_im_test")

summ <- rbindlist(res); save(summ, "t12_vcov_summary")
save(rbindlist(wrows), "t12_delta_by_vcov"); save(rbindlist(walds), "t12_wald")
cat("\n##### まとめ #####\n")
print(summ[, .(vcov, min_eig_d = signif(min_eig_d, 2), delta_2010 = round(delta_2010, 6),
               se_2010 = round(se_2010, 6), t_2010 = round(t_2010, 2),
               wald_all = round(wald_all, 1), wald_no09 = round(wald_no09, 1),
               se_main = round(se_main, 6))], row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. 都県レベルの重み分解（IM が示した群馬集中を定量化する）
#    T11 と同じ恒等式を年ではなく都県について当てる
#    delta_pool = sum_p omega_p * delta_p,  omega_p = sum_{i in p, t} Dtil*D / sum Dtil^2
# -----------------------------------------------------------------------------
cat("\n##### 都県レベルの重み分解 #####\n")
FML_PI <- as.formula(paste("log_yield ~", HEAT, "+ i(pref, h_consol) +",
                           paste(CTRL, collapse = " + "), "| city_id + year"))
m_pi <- feols(FML_PI, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
cf <- coef(m_pi); pp <- strsplit(names(cf), ":", fixed = TRUE)
okp <- vapply(pp, function(z) length(z) >= 4 && z[1] == "pref" && z[4] == "h_consol", logical(1))
pr_lab <- as.integer(vapply(pp[okp], function(z) z[3], character(1)))
dp <- data.table(pref = pr_lab, delta_p = unname(cf[okp]))[order(pref)]

m_pool1 <- feols(as.formula(paste("log_yield ~", HEAT, "+ h_consol +",
                                  paste(CTRL, collapse = " + "), "| city_id + year")),
                 s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
d_pool <- unname(coef(m_pool1)["h_consol"])

wp <- s[, .(omega_p = sum(Dtil * Draw) / ss4), by = pref][order(pref)]
dp <- merge(wp, dp, by = "pref")
dp[, contrib := omega_p * delta_p]
dp[, share := contrib / d_pool]
print(dp[, .(pref, omega_p = round(omega_p, 4), delta_p = round(delta_p, 6),
             contrib = round(contrib, 7), share = round(share, 3))], row.names = FALSE)
cat(sprintf("sum omega_p = %.8f （理論値 1）\n", dp[, sum(omega_p)]))
cat(sprintf("sum omega_p * delta_p = %+.8f / プール delta = %+.8f / 差 %+.2e\n",
            dp[, sum(contrib)], d_pool, dp[, sum(contrib)] - d_pool))
save(dp, "t12_pref_decomp")

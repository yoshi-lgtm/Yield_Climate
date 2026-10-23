# =============================================================================
# T18: omega_t 分解を機序ベースの本番スペックで再計算する
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t18_omega_mechanism.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_centroids.csv / city_coastdist.csv / city_forest.csv
# 出力 : research/model/output/t18_*.csv と標準出力
#
# 本番スペックは「事前の機序」で選ぶ（検定統計量では選ばない）。
#   coast      : 海洋性の緩和（海風・日較差の縮小・高湿度）
#   forest_nbr : 局所の遮蔽・蒸散冷却・冷気流出（田の3x3近傍）
# 相互相関は 0.223 で、両方が有意に推定できる唯一のペアでもある。
# elev は機序が弱く（逓減率は気温の水準に効き、水準は HD35 自体に入っている）、
# lat/lon は機序がないので頑健性に回す。
#
# 恒等式の一般化：プール式に年別の競合交差項が入ると汚染項が増える。
#   y = FE + sum_t 1[t](beta_t h + delta_t D + gamma_t G + eta_t F) + X g + e
#   D = h*consol30, G = h*coast, F = h*forest
# Dtil を D の「他の全プール回帰子」への射影残差とすると
#   delta_pool = sum_t w_D,t delta_t        （主項）
#              + sum_t w_h,t beta_t         （汚染項1：heat 主効果由来）
#              + sum_t w_G,t gamma_t        （汚染項2：coast 交差項由来）
#              + sum_t w_F,t eta_t          （汚染項3：forest 交差項由来）
#   w_A,t = sum_i Dtil_it A_it / sum Dtil^2
# 年×水準項・CTRL・FE は Dtil と直交するので寄与しない。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

pick <- function(cands) cands[file.exists(cands)][1]
PANEL <- pick(c("research/model/kanto_panel_2000_2023.csv",
                "Yield_Climate/research/model/kanto_panel_2000_2023.csv"))
CENTR <- pick(c("results/city_centroids.csv", "Yield_Climate/results/city_centroids.csv"))
CDIST <- pick(c("results/city_coastdist.csv", "Yield_Climate/results/city_coastdist.csv"))
FRST  <- pick(c("results/city_forest.csv", "Yield_Climate/results/city_forest.csv"))
OUT_DIR <- "research/model/output"
CTRL    <- c("GSR_78", "APCP_78")
HEAT    <- "obs_HD35"
CONLEY  <- 50

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d <- merge(d, fread(CENTR, encoding = "UTF-8")[, .(city_id, lon, lat)], by = "city_id")
d <- merge(d, fread(CDIST, encoding = "UTF-8")[, .(city_id, coast = coast_km)], by = "city_id")
d <- merge(d, fread(FRST, encoding = "UTF-8")[, .(city_id, forest = forest_nbr)], by = "city_id")

vars <- c("log_yield", HEAT, "consol30", "coast", "forest", CTRL,
          "city_id", "year", "lon", "lat")
s <- d[complete.cases(d[, ..vars])]
s[, pref := city_id %/% 1000L]
years <- sort(unique(s$year)); REF <- years[1]
cat(sprintf("n = %d, 市町村 = %d, 年 = %d\n", nrow(s), uniqueN(s$city_id), length(years)))

# 競合は市町村断面で標準化（consol30 は既存文書と数値を揃えるため生のまま）
cs <- unique(s[, .(city_id, coast, forest)])
for (v in c("coast", "forest")) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = v, value = (s[[v]] - mu) / sg)
}
cat(sprintf("coast と forest の相関: %.3f\n",
            cor(unique(s[, .(city_id, coast, forest)])$coast,
                unique(s[, .(city_id, coast, forest)])$forest)))

s[, D := get(HEAT) * consol30]
s[, G := get(HEAT) * coast]
s[, F_ := get(HEAT) * forest]

mk_lev <- function(v, tag) {
  nms <- character(0)
  for (y in years[-1]) {
    nm <- sprintf("L%s_%d", tag, y)
    set(s, j = nm, value = s[[v]] * as.integer(s$year == y))
    nms <- c(nms, nm)
  }
  nms
}
Lc <- mk_lev("consol30", "c")
Lk <- mk_lev("coast", "k")
Lf <- mk_lev("forest", "f")

# -----------------------------------------------------------------------------
# 共通処理
#   slopes : 年別に交差させる変数（D 以外）。汚染項の源になる
#   levels : 年×水準項
# -----------------------------------------------------------------------------
pick_idx <- function(mm, v) {
  parts <- strsplit(names(coef(mm)), ":", fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) >= 4 && z[1] == "year" && z[4] == v, logical(1))
  i <- which(ok)
  i[order(as.integer(vapply(parts[i], function(z) z[3], character(1))))]
}
wald_equal <- function(b, V) {
  K <- length(b); R <- matrix(0, K - 1, K)
  for (i in 1:(K - 1)) { R[i, i] <- 1; R[i, i + 1] <- -1 }
  Rb <- R %*% b; RV <- R %*% V %*% t(R)
  sol <- tryCatch(solve(RV), error = function(e) NULL)
  if (is.null(sol)) return(c(NA_real_, NA_real_))
  st <- as.numeric(t(Rb) %*% sol %*% Rb)
  c(st, pchisq(st, K - 1, lower.tail = FALSE))
}

run <- function(label, slopes, levels) {
  cat(sprintf("\n==================== %s ====================\n", label))
  rhs_pool <- c("D", slopes, levels, CTRL)
  f_pool <- as.formula(sprintf("log_yield ~ %s | city_id + year",
                               paste(rhs_pool, collapse = " + ")))
  mp <- feols(f_pool, s, cluster = ~city_id, fixef.rm = "none",
              ssc = ssc(K.fixef = "full"))
  dpool <- unname(coef(mp)["D"]); sepool <- unname(se(mp)["D"])
  cat(sprintf("プール delta = %+.6f (SE %.6f, t %+.2f)\n", dpool, sepool, dpool / sepool))

  rhs_sat <- c("i(year, D)", sprintf("i(year, %s)", slopes), levels, CTRL)
  f_sat <- as.formula(sprintf("log_yield ~ %s | city_id + year",
                              paste(rhs_sat, collapse = " + ")))
  ms <- feols(f_sat, s, cluster = ~city_id, fixef.rm = "none",
              ssc = ssc(K.fixef = "full"))
  id_D <- pick_idx(ms, "D")
  stopifnot(length(id_D) == length(years))
  bh <- coef(ms)

  # Dtil : D を他の全プール回帰子に射影した残差
  proj <- c("D", slopes, levels, CTRL)
  M <- as.data.frame(demean(s[, ..proj], s[, .(city_id, year)]))
  f4 <- as.formula(paste("D ~", paste(c(slopes, levels, CTRL), collapse = " + "), "- 1"))
  Dtil <- as.numeric(residuals(lm(f4, data = M)))
  ss <- sum(Dtil^2)

  wt <- data.table(year = s$year, Dtil = Dtil, D = s$D)[
    , .(wD = sum(Dtil * D) / ss), by = year][order(year)]
  stopifnot(identical(wt$year, years))

  main <- sum(wt$wD * bh[id_D])
  terms <- list()
  for (v in slopes) {
    idv <- pick_idx(ms, v)
    wv <- data.table(year = s$year, Dtil = Dtil, A = s[[v]])[
      , .(w = sum(Dtil * A) / ss), by = year][order(year)]
    terms[[v]] <- list(w = wv$w, b = bh[idv], val = sum(wv$w * bh[idv]), idx = idv)
  }
  tot <- main + sum(vapply(terms, function(z) z$val, numeric(1)))
  cat(sprintf("恒等式: 主項 %+.6f", main))
  for (v in names(terms)) cat(sprintf(" + 汚染(%s) %+.6f", v, terms[[v]]$val))
  cat(sprintf(" = %+.6f （差 %.2e）\n", tot, dpool - tot))
  cat(sprintf("        sum(omega) = %.8f, omega が負の年 = %d / %d\n",
              sum(wt$wD), sum(wt$wD < 0), length(years)))
  for (v in names(terms)) {
    cat(sprintf("        汚染(%s) はプール値の %+.1f%%\n",
                v, 100 * terms[[v]]$val / dpool))
  }

  i10 <- which(years == 2010L)
  out <- list()
  for (vl in c("cluster_city", sprintf("conley_%dkm", CONLEY))) {
    V <- if (vl == "cluster_city") vcov(ms) else
      suppressWarnings(vcov_conley(ms, lat = ~lat, lon = ~lon, cutoff = CONLEY,
                                   ssc = ssc(K.fixef = "full")))
    Vd <- V[id_D, id_D]; se_d <- sqrt(pmax(diag(Vd), 0))
    w <- wald_equal(bh[id_D], Vd)
    v_main <- as.numeric(t(wt$wD) %*% Vd %*% wt$wD)
    cat(sprintf("\n-- %s --\n", vl))
    cat(sprintf("(a) delta_2010 %+.6f (SE %.6f, t %+.2f, p %.4f)\n",
                bh[id_D][i10], se_d[i10], bh[id_D][i10] / se_d[i10],
                2 * pnorm(-abs(bh[id_D][i10] / se_d[i10]))))
    cat(sprintf("(c) chi2(%d) = %.1f, p = %.3g\n", length(years) - 1, w[1], w[2]))
    cat(sprintf("    主項 %+.6f (SE %.6f, t %+.2f)\n",
                main, sqrt(max(v_main, 0)), main / sqrt(max(v_main, 1e-300))))
    for (v in names(terms)) {
      vv <- as.numeric(t(terms[[v]]$w) %*% V[terms[[v]]$idx, terms[[v]]$idx] %*%
                         terms[[v]]$w)
      cat(sprintf("    汚染(%-7s) %+.6f (SE %.6f, t %+.2f)\n", v,
                  terms[[v]]$val, sqrt(max(vv, 0)),
                  terms[[v]]$val / sqrt(max(vv, 1e-300))))
    }
    out[[vl]] <- data.table(spec = label, vcov = vl, delta_pool = dpool,
      delta_2010 = bh[id_D][i10], se_2010 = se_d[i10],
      t_2010 = bh[id_D][i10] / se_d[i10],
      omega_2010 = wt$wD[i10], share_2010 = 100 * wt$wD[i10] * bh[id_D][i10] / dpool,
      n_neg_omega = sum(wt$wD < 0), wald = w[1], wald_p = w[2],
      main = main, se_main = sqrt(max(v_main, 0)),
      contam_total = sum(vapply(terms, function(z) z$val, numeric(1))),
      contam_pct = 100 * sum(vapply(terms, function(z) z$val, numeric(1))) / dpool,
      identity_err = dpool - tot)
  }
  ct <- rbindlist(lapply(names(terms), function(v)
    data.table(spec = label, source = v, value = terms[[v]]$val,
               pct = 100 * terms[[v]]$val / dpool)))
  list(sum = rbindlist(out),
       byyear = data.table(spec = label, year = years, omega = wt$wD,
                           delta = bh[id_D],
                           contrib = wt$wD * bh[id_D],
                           share = 100 * wt$wD * bh[id_D] / dpool),
       contam = ct)
}

# -----------------------------------------------------------------------------
# 3段のはしご
# -----------------------------------------------------------------------------
r1 <- run("1. 地理なし・年×c なし（T11）", HEAT, character(0))
r2 <- run("2. ＋年×c 飽和（T15）", HEAT, Lc)
r3 <- run("3. ＋coast＋forest（機序ベース・本番候補）",
          c(HEAT, "G", "F_"), c(Lc, Lk, Lf))

sm <- rbindlist(list(r1$sum, r2$sum, r3$sum))
by <- rbindlist(list(r1$byyear, r2$byyear, r3$byyear))
ct <- rbindlist(list(r1$contam, r2$contam, r3$contam))
save(sm, "t18_summary"); save(by, "t18_by_year"); save(ct, "t18_contamination")

cat("\n==================== まとめ（市町村クラスター） ====================\n")
print(sm[vcov == "cluster_city",
  .(spec, d_pool = round(delta_pool, 6), d_2010 = round(delta_2010, 6),
    t_2010 = round(t_2010, 2), omega_2010 = round(omega_2010, 4),
    share_2010 = round(share_2010, 1), neg = n_neg_omega,
    chi2 = round(wald, 1), p = signif(wald_p, 3),
    contam_pct = round(contam_pct, 1), err = signif(identity_err, 2))],
  row.names = FALSE)
cat("\n汚染項の内訳:\n")
print(ct[, .(spec, source, value = round(value, 6), pct = round(pct, 1))],
      row.names = FALSE)
cat("\n重みの重い年（本番候補）:\n")
print(by[spec == "3. ＋coast＋forest（機序ベース・本番候補）"][order(-abs(omega))][1:8,
  .(year, omega = round(omega, 4), delta = round(delta, 5),
    contrib = round(contrib, 6), share = round(share, 1))], row.names = FALSE)
cat("\n出力: research/model/output/t18_{summary,by_year,contamination}.csv\n")

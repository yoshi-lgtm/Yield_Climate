# =============================================================================
# T22: heat の形状と用量依存 — プール係数は「どの年」でなく「どの用量域」から来るか
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t22_heat_shape.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_centroids.csv / city_coastdist.csv / city_forest.csv
# 出力 : research/model/output/t22_*.csv と標準出力
#
# -----------------------------------------------------------------------------
# 中核の着想
# -----------------------------------------------------------------------------
# T18 の恒等式は「プールモデルの Z への射影残差 Dtil を、飽和モデルの当てはめ値に
# ぶつける」というだけなので、飽和の軸は年である必要がない。
#
#   delta_pool = Dtil'Y / Dtil'Dtil
#              = Dtil'Yhat_sat / Dtil'Dtil        (飽和残差は Dtil と直交)
#              = sum_g w_D,g delta_g + sum_A sum_g w_A,g A_g
#   w_A,g = sum_{i,t in g} Dtil_it A_it / sum Dtil^2
#
# 成立条件はただ一つ、「飽和モデルの列空間がプールモデルの列を含む」こと。
# D = sum_g 1[g] D なので、g が年でも用量ビンでも年×用量セルでも成立する。
#
# したがって同じ Dtil を
#   by = "year" : T18 の再現（78.3% / 28.2% / -5.2% / -1.4%）
#   by = "hb"   : 用量ビン分解 ← T22 の本体
#   by = "cell" : 年×用量セル分解
# に分けられる。プールモデルは一切変えない（報告されている推定量のまま）。
#
# 用量ビン版では「delta_b が全ビンで等しい」= 「交差項が h について線形」なので、
# 同質性 Wald がそのまま関数形の検定になる。
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
CTRL   <- c("GSR_78", "APCP_78")
HEAT   <- "obs_HD35"          # 本番
CONLEY <- 50
# 用量ビン：ゼロ / 1-2 / 3-5 / 6-10 / 11+ （事前に決め、統計量で動かさない）
BRK <- c(-1, 0, 2, 5, 10, Inf)
LAB <- c("b0", "b1_2", "b3_5", "b6_10", "b11p")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

# -----------------------------------------------------------------------------
# データ
# -----------------------------------------------------------------------------
d0 <- fread(PANEL, encoding = "UTF-8")
d0 <- merge(d0, fread(CENTR, encoding="UTF-8")[, .(city_id, lon, lat)], by="city_id")
d0 <- merge(d0, fread(CDIST, encoding="UTF-8")[, .(city_id, coast = coast_km)], by="city_id")
d0 <- merge(d0, fread(FRST , encoding="UTF-8")[, .(city_id, forest = forest_nbr)], by="city_id")

prep <- function(H) {
  vars <- c("log_yield", H, "consol30", "coast", "forest", CTRL,
            "city_id", "year", "lon", "lat")
  s <- d0[complete.cases(d0[, ..vars])]
  s[, pref := city_id %/% 1000L]
  # 競合は市町村断面で標準化（consol30 は既存文書と揃えるため生のまま）
  cs <- unique(s[, .(city_id, coast, forest)])
  for (v in c("coast", "forest")) {
    mu <- mean(cs[[v]]); sg <- sd(cs[[v]]); set(s, j = v, value = (s[[v]] - mu) / sg)
  }
  s[, h  := get(H)]
  s[, D  := h * consol30]
  s[, G  := h * coast]
  s[, F_ := h * forest]
  s[, h2 := h^2]
  s[, D2 := h2 * consol30]; s[, G2 := h2 * coast]; s[, F2 := h2 * forest]
  s[, hb   := cut(h, breaks = BRK, labels = LAB)]
  s[, cell := factor(paste0(year, "_", as.character(hb)))]
  yrs <- sort(unique(s$year))
  for (tg in c("c","k","f")) {
    v <- c(c="consol30", k="coast", f="forest")[[tg]]
    for (y in yrs[-1]) set(s, j = sprintf("L%s_%d", tg, y),
                           value = s[[v]] * as.integer(s$year == y))
  }
  attr(s, "years") <- yrs
  attr(s, "levels") <- unlist(lapply(c("c","k","f"), function(tg) sprintf("L%s_%d", tg, yrs[-1])))
  s
}

# -----------------------------------------------------------------------------
# 汎用の恒等式分解： by の水準で飽和させる
# -----------------------------------------------------------------------------
# 係数名は fixest の "f::level:x" 形式。":" で割ると c(f, "", level, x)。
idx_of <- function(mm, fac, v) {
  parts <- strsplit(names(coef(mm)), ":", fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) >= 4 && z[1] == fac && z[4] == v, logical(1))
  i <- which(ok)
  lv <- vapply(parts[i], function(z) z[3], character(1))
  list(idx = i, lev = lv)
}
wald_equal <- function(b, V) {
  K <- length(b); if (K < 2) return(c(NA_real_, NA_real_, NA_integer_))
  R <- matrix(0, K - 1, K)
  for (i in 1:(K - 1)) { R[i, i] <- 1; R[i, i + 1] <- -1 }
  Rb <- R %*% b; RV <- R %*% V %*% t(R)
  sol <- tryCatch(solve(RV), error = function(e) NULL)
  if (is.null(sol)) return(c(NA_real_, NA_real_, K - 1L))
  st <- as.numeric(t(Rb) %*% sol %*% Rb)
  c(st, pchisq(st, K - 1, lower.tail = FALSE), K - 1L)
}

decomp <- function(s, by, label, quiet = FALSE) {
  LEV <- attr(s, "levels")
  slopes <- c("h", "G", "F_")                       # 汚染源（D 以外の年別/用量別交差）
  rhs_pool <- c("D", slopes, LEV, CTRL)

  # --- プールモデル（報告されている推定量。by に依らず不変） ---
  mp <- feols(as.formula(sprintf("log_yield ~ %s | city_id + year",
                                 paste(rhs_pool, collapse = " + "))),
              s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
  dpool <- unname(coef(mp)["D"]); sepool <- unname(se(mp)["D"])

  # --- 飽和モデル（by の水準で全スロープを飽和） ---
  rhs_sat <- c(sprintf("i(%s, D)", by),
               sprintf("i(%s, %s)", by, slopes), LEV, CTRL)
  ms <- feols(as.formula(sprintf("log_yield ~ %s | city_id + year",
                                 paste(rhs_sat, collapse = " + "))),
              s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
  bh <- coef(ms)

  # --- Dtil : D を「他の全プール回帰子」に射影した残差（by に依らず同一） ---
  proj <- c("D", slopes, LEV, CTRL)
  M <- as.data.frame(demean(s[, ..proj], s[, .(city_id, year)]))
  Dtil <- as.numeric(residuals(lm(
    as.formula(paste("D ~", paste(c(slopes, LEV, CTRL), collapse = " + "), "- 1")), data = M)))
  ss <- sum(Dtil^2)

  gv <- as.character(s[[by]])
  wsum <- function(v) {           # 水準ごとの cross-product 重み
    dt <- data.table(g = gv, Dtil = Dtil, A = s[[v]])[, .(w = sum(Dtil * A) / ss), by = g]
    setNames(dt$w, dt$g)
  }

  # --- 主項 ---
  iD <- idx_of(ms, by, "D")
  wD_all <- wsum("D")
  wD <- wD_all[iD$lev]; wD[is.na(wD)] <- 0
  bD <- bh[iD$idx]
  main <- sum(wD * bD)

  # --- 汚染項 ---
  terms <- list()
  for (v in slopes) {
    iv <- idx_of(ms, by, v)
    wv_all <- wsum(v)
    wv <- wv_all[iv$lev]; wv[is.na(wv)] <- 0
    terms[[v]] <- list(w = unname(wv), b = bh[iv$idx], idx = iv$idx, lev = iv$lev,
                       val = sum(wv * bh[iv$idx]))
  }
  tot <- main + sum(vapply(terms, function(z) z$val, numeric(1)))

  if (!quiet) {
    cat(sprintf("\n==================== %s （飽和の軸: %s） ====================\n", label, by))
    cat(sprintf("プール delta = %+.6f (SE %.6f, t %+.2f)   ※ by に依らず不変\n",
                dpool, sepool, dpool / sepool))
    cat(sprintf("恒等式: 主項 %+.6f", main))
    for (v in names(terms)) cat(sprintf(" + 汚染(%s) %+.6f", v, terms[[v]]$val))
    cat(sprintf(" = %+.6f （差 %.2e）\n", tot, dpool - tot))
    cat(sprintf("        sum(w_D) = %.8f, w_D<0 の水準 = %d / %d\n",
                sum(wD), sum(wD < 0), length(wD)))
    cat(sprintf("        主項 %+.1f%%", 100 * main / dpool))
    for (v in names(terms)) cat(sprintf(", 汚染(%s) %+.1f%%", v, 100 * terms[[v]]$val / dpool))
    cat("\n")
  }

  # --- 推論 ---
  inf <- list()
  for (vl in c("cluster_city", sprintf("conley_%dkm", CONLEY))) {
    V <- if (vl == "cluster_city") vcov(ms) else
      tryCatch(suppressWarnings(vcov_conley(ms, lat = ~lat, lon = ~lon, cutoff = CONLEY,
                                            ssc = ssc(K.fixef = "full"))),
               error = function(e) NULL)
    if (is.null(V)) next
    Vd <- V[iD$idx, iD$idx, drop = FALSE]
    w <- wald_equal(bD, Vd)
    v_main <- as.numeric(t(wD) %*% Vd %*% wD)
    row <- data.table(spec = label, by = by, vcov = vl, delta_pool = dpool,
                      main = main, main_pct = 100 * main / dpool,
                      se_main = sqrt(max(v_main, 0)),
                      t_main = main / sqrt(max(v_main, 1e-300)),
                      wald = w[1], wald_df = w[3], wald_p = w[2],
                      n_neg_w = sum(wD < 0), n_lev = length(wD),
                      identity_err = dpool - tot)
    for (v in names(terms)) {
      vv <- as.numeric(t(terms[[v]]$w) %*%
                         V[terms[[v]]$idx, terms[[v]]$idx, drop = FALSE] %*% terms[[v]]$w)
      row[[paste0("contam_", v)]] <- terms[[v]]$val
      row[[paste0("contam_", v, "_pct")]] <- 100 * terms[[v]]$val / dpool
      row[[paste0("contam_", v, "_t")]] <- terms[[v]]$val / sqrt(max(vv, 1e-300))
    }
    if (!quiet) {
      cat(sprintf("\n-- %s --\n", vl))
      cat(sprintf("    主項 %+.6f (%+6.1f%%) SE %.6f t %+.2f\n",
                  main, 100 * main / dpool, sqrt(max(v_main, 0)),
                  main / sqrt(max(v_main, 1e-300))))
      for (v in names(terms)) {
        vv <- as.numeric(t(terms[[v]]$w) %*%
                           V[terms[[v]]$idx, terms[[v]]$idx, drop = FALSE] %*% terms[[v]]$w)
        cat(sprintf("    汚染(%-3s) %+.6f (%+6.1f%%) SE %.6f t %+.2f\n", v, terms[[v]]$val,
                    100 * terms[[v]]$val / dpool, sqrt(max(vv, 0)),
                    terms[[v]]$val / sqrt(max(vv, 1e-300))))
      }
      cat(sprintf("    同質性 Wald: chi2(%d) = %.1f, p = %.4g",
                  w[3], w[1], w[2]))
      if (by == "hb") cat("   ← 「交差項が h について線形」の検定")
      cat("\n")
    }
    inf[[vl]] <- row
  }

  # --- 水準別の表 ---
  nn <- data.table(g = gv)[, .(n = .N), by = g]
  hm <- data.table(g = gv, h = s$h)[, .(mean_h = mean(h)), by = g]
  per <- merge(merge(data.table(lev = iD$lev, w = unname(wD), coef = unname(bD)),
                     nn, by.x = "lev", by.y = "g", all.x = TRUE),
               hm, by.x = "lev", by.y = "g", all.x = TRUE)
  per[, `:=`(spec = label, by = by, contrib = w * coef,
             share = 100 * w * coef / dpool)]

  list(sum = rbindlist(inf, fill = TRUE), per = per,
       dpool = dpool, Dtil = Dtil, ss = ss)
}

# =============================================================================
# PART 1  分布と年別サポート
# =============================================================================
cat("=============================================================\n")
cat(" T22  PART 1: 分布と年別サポート\n")
cat("=============================================================\n")
s <- prep(HEAT)
years <- attr(s, "years")
cat(sprintf("n = %d, 市町村 = %d, 年 = %d\n", nrow(s), uniqueN(s$city_id), length(years)))

HCAND <- c("obs_HD33","obs_HD34","obs_HD35","obs_HDD26","c40_HD35","c50_HD35","c20_HD35")
HCAND <- HCAND[HCAND %in% names(d0)]
dist <- rbindlist(lapply(HCAND, function(v) {
  x <- s[[v]]
  data.table(var = v, mean = mean(x), sd = sd(x), median = median(x),
             zero_pct = 100 * mean(x == 0), max = max(x),
             cor_prod = cor(x, s$h))
}))
print(dist[, .(var, mean = round(mean,3), sd = round(sd,3), median,
               zero_pct = round(zero_pct,1), max, cor_prod = round(cor_prod,3))],
      row.names = FALSE)
save(dist, "t22_dist")

sup <- s[, .(n = .N, nz = sum(h > 0), pct_nz = round(100*mean(h>0),1),
             n_uniq_pos = uniqueN(h[h > 0]), mean_h = round(mean(h),2),
             p90 = quantile(h, .9), max = max(h)), by = year][order(year)]
cat("\n年別サポート:\n"); print(sup, row.names = FALSE)
cat(sprintf("\n非ゼロ市町村 < 30 の年    : %s\n", paste(sup[nz < 30]$year, collapse=", ")))
cat(sprintf("非ゼロの異なる値 <= 5 の年: %s\n", paste(sup[n_uniq_pos <= 5]$year, collapse=", ")))
cat(sprintf("→ 年別の曲率が引ける年    : %d / %d\n",
            sup[n_uniq_pos >= 6, .N], nrow(sup)))
save(sup, "t22_support")

cat("\n用量ビンの分布:\n")
bn <- s[, .(n = .N, mean_h = round(mean(h),2), min_h = min(h), max_h = max(h),
            n_city = uniqueN(city_id), n_year = uniqueN(year)), by = hb][order(hb)]
print(bn, row.names = FALSE)
save(bn, "t22_bins")

# =============================================================================
# PART 2  飽和の軸を替える：年 / 用量ビン / 年×用量
# =============================================================================
cat("\n\n=============================================================\n")
cat(" T22  PART 2: 同じ Dtil を、年で分けるか用量で分けるか\n")
cat("=============================================================\n")
r_year <- decomp(s, "year", "本番スペック R3")
r_dose <- decomp(s, "hb",   "本番スペック R3")
r_cell <- decomp(s, "cell", "本番スペック R3", quiet = TRUE)

cat("\n---- 年別（T18 の再現） ----\n")
print(r_year$per[order(-abs(w))][1:8, .(year = lev, n, mean_h = round(mean_h,2),
  w = round(w,4), delta = round(coef,5), contrib = round(contrib,6),
  share = round(share,1))], row.names = FALSE)

cat("\n---- 用量ビン別（T22 の本体） ----\n")
print(r_dose$per[order(lev), .(bin = lev, n, mean_h = round(mean_h,2),
  w = round(w,4), delta_b = round(coef,6), contrib = round(contrib,6),
  share = round(share,1))], row.names = FALSE)

cat("\n---- 年×用量セル（重みの上位 12） ----\n")
print(r_cell$per[order(-abs(w))][1:12, .(cell = lev, n, mean_h = round(mean_h,2),
  w = round(w,4), delta = round(coef,6), contrib = round(contrib,6),
  share = round(share,1))], row.names = FALSE)
cat(sprintf("セル数 = %d, うち n<5 のセル = %d（係数は不安定。重みは推定誤差を持たない）\n",
            nrow(r_cell$per), r_cell$per[n < 5, .N]))

sm <- rbindlist(list(r_year$sum, r_dose$sum, r_cell$sum), fill = TRUE)
save(sm, "t22_decomp_summary")
save(rbindlist(list(r_year$per, r_dose$per, r_cell$per), fill = TRUE), "t22_decomp_by_level")

# 年×用量の重みマップ（重みのみ。推定誤差を持たないので厳密）
Dtil <- r_year$Dtil; ss <- r_year$ss
wm <- data.table(year = s$year, hb = as.character(s$hb), Dtil = Dtil, D = s$D)[
  , .(n = .N, w = sum(Dtil * D) / ss), by = .(year, hb)]
wmat <- dcast(wm, year ~ hb, value.var = "w", fill = 0)
cat("\n---- 年×用量の重みマップ w_{t,b}（合計 1） ----\n")
print(wmat[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)], row.names = FALSE)
cat(sprintf("合計 = %.8f\n", sum(wm$w)))
cat("\n用量ビンごとの重み合計（年を潰す）:\n")
print(wm[, .(w = round(sum(w), 4), n = sum(n)), by = hb][order(hb)], row.names = FALSE)
save(wm, "t22_year_dose_weights")

# =============================================================================
# PART 3  関数形の検定と delta(h) のプロファイル（閾値を替えて頑健性も見る）
# =============================================================================
cat("\n\n=============================================================\n")
cat(" T22  PART 3: 交差項の曲率と delta(h)\n")
cat("=============================================================\n")

curv <- list(); prof <- list()
for (H in c("obs_HD35", "obs_HD34", "obs_HD33")) {
  if (!H %in% names(d0)) next
  sh <- prep(H); LEV <- attr(sh, "levels")
  lin <- feols(as.formula(paste("log_yield ~ h + D + G + F_ +",
    paste(c(LEV, CTRL), collapse=" + "), "| city_id + year")),
    sh, cluster=~city_id, fixef.rm="none", ssc=ssc(K.fixef="full"))
  qua <- feols(as.formula(paste("log_yield ~ h + h2 + D + D2 + G + G2 + F_ + F2 +",
    paste(c(LEV, CTRL), collapse=" + "), "| city_id + year")),
    sh, cluster=~city_id, fixef.rm="none", ssc=ssc(K.fixef="full"))
  w4 <- wald(qua, "h2|D2|G2|F2", print = FALSE)
  b1 <- coef(qua)["D"]; b2 <- coef(qua)["D2"]; Vq <- vcov(qua)
  ii <- match(c("D","D2"), names(coef(qua)))
  hstar <- -b1 / (2 * b2)
  cat(sprintf("\n---- %s ----\n", H))
  cat(sprintf("線形プール  D  = %+.7f (SE %.7f, t %+.2f)\n",
              coef(lin)["D"], se(lin)["D"], coef(lin)["D"]/se(lin)["D"]))
  cat(sprintf("二次        D  = %+.7f (SE %.7f, t %+.2f)\n",
              b1, se(qua)["D"], b1/se(qua)["D"]))
  cat(sprintf("            D2 = %+.7f (SE %.7f, t %+.2f)  ← 交差項の曲率\n",
              b2, se(qua)["D2"], b2/se(qua)["D2"]))
  cat(sprintf("            G2 = %+.7f (SE %.7f, t %+.2f)  ← coast 側の曲率\n",
              coef(qua)["G2"], se(qua)["G2"], coef(qua)["G2"]/se(qua)["G2"]))
  cat(sprintf("二次項4本の同時ゼロ: F = %.2f, p = %.4g\n", w4$stat, w4$p))
  cat(sprintf("符号反転点 h* = %.2f 日（標本平均 %.2f 日、h>h* の標本比 %.1f%%）\n",
              hstar, mean(sh$h), 100 * mean(sh$h > hstar)))
  curv[[H]] <- data.table(var = H, lin_D = coef(lin)["D"], lin_t = coef(lin)["D"]/se(lin)["D"],
    qua_D = b1, qua_D_t = b1/se(qua)["D"], qua_D2 = b2, qua_D2_t = b2/se(qua)["D2"],
    qua_G2 = coef(qua)["G2"], qua_G2_t = coef(qua)["G2"]/se(qua)["G2"],
    F_joint = w4$stat, p_joint = w4$p, h_star = hstar,
    mean_h = mean(sh$h), pct_above = 100*mean(sh$h > hstar))
  grid <- sort(unique(c(0, 1, 3, mean(sh$h), hstar, 5, 10, 15, 25)))
  prof[[H]] <- rbindlist(lapply(grid, function(hh) {
    g <- c(1, 2*hh); v <- as.numeric(t(g) %*% Vq[ii, ii] %*% g); e <- b1 + 2*b2*hh
    data.table(var = H, h = hh, delta = e, se = sqrt(v), t = e/sqrt(max(v,1e-300)))
  }))
  if (H == HEAT) {
    cat(sprintf("%-10s %13s %11s %8s\n", "h", "delta(h)", "SE", "t"))
    for (k in seq_len(nrow(prof[[H]]))) with(prof[[H]][k],
      cat(sprintf("%-10.2f %+13.6f %11.6f %+8.2f\n", h, delta, se, t)))
  }
}
save(rbindlist(curv), "t22_curvature")
save(rbindlist(prof), "t22_delta_profile")

# 用量ビン × 修飾変数（水準シフト版。基底非依存の非線形プロファイル）
cat("\n---- 用量ビン × 修飾変数（水準シフト、参照 = h 0日） ----\n")
LEV <- attr(s, "levels")
mb <- feols(as.formula(paste(
  "log_yield ~ i(hb, ref='b0') + i(hb, consol30, ref='b0') + i(hb, coast, ref='b0')",
  "+ i(hb, forest, ref='b0') +", paste(c(LEV, CTRL), collapse=" + "), "| city_id + year")),
  s, cluster=~city_id, fixef.rm="none", ssc=ssc(K.fixef="full"))
cf <- coef(mb); sd_ <- se(mb); V <- vcov(mb)
bins <- LAB[-1]
mh <- s[hb != "b0", .(mh = mean(h)), by = hb]; setkey(mh, hb)
mv <- mh[match(bins, as.character(hb))]$mh
binprof <- rbindlist(lapply(c("", ":consol30", ":coast", ":forest"), function(tg) {
  rbindlist(lapply(seq_along(bins), function(j) {
    k <- paste0("hb::", bins[j], tg)
    if (!k %in% names(cf)) return(NULL)
    data.table(target = if (tg == "") "h_main" else sub("^:", "", tg),
               bin = bins[j], mean_h = mv[j], coef = cf[k], se = sd_[k],
               t = cf[k]/sd_[k], per_day = cf[k]/mv[j])
  }))
}))
for (tg in unique(binprof$target)) {
  cat(sprintf("\n[%s]\n", tg))
  print(binprof[target == tg, .(bin, mean_h = round(mean_h,2), coef = round(coef,6),
    se = round(se,6), t = round(t,2), per_day = round(per_day,7))], row.names = FALSE)
}
# 比例性検定：coef_b / mean_h_b が全ビンで等しい（= h について線形）
cat("\n比例性検定（係数/平均h が全ビンで等しい = h について線形）:\n")
prop <- rbindlist(lapply(c("consol30","coast","forest"), function(tg) {
  ks <- paste0("hb::", bins, ":", tg); ii <- match(ks, names(cf))
  if (anyNA(ii)) return(NULL)
  K <- length(ii); R <- matrix(0, K-1, length(cf))
  for (i in 1:(K-1)) { R[i, ii[i]] <- 1/mv[i]; R[i, ii[i+1]] <- -1/mv[i+1] }
  Rb <- R %*% cf
  st <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb)
  data.table(target = tg, chi2 = st, df = K-1L, p = pchisq(st, K-1, lower.tail=FALSE))
}))
print(prop[, .(target, chi2 = round(chi2,2), df, p = signif(p,4))], row.names = FALSE)
save(binprof, "t22_bin_profile"); save(prop, "t22_proportionality")

# =============================================================================
# PART 4  窓長のはしご（S1。T22 に取り込み）
# =============================================================================
cat("\n\n=============================================================\n")
cat(" T22  PART 4: 窓長のはしご（S1）\n")
cat("=============================================================\n")
wl <- rbindlist(lapply(c("obs_HD35","c40_HD35","c50_HD35","c20_HD35"), function(H) {
  if (!H %in% names(d0)) return(NULL)
  rr <- decomp(prep(H), "year", H, quiet = TRUE)
  rr$sum[vcov == "cluster_city", .(heat = H, delta_pool, main_pct,
    coast_pct = contam_G_pct, coast_t = contam_G_t,
    heat_pct = contam_h_pct, forest_pct = contam_F__pct,
    wald, wald_p, n_neg_w, identity_err)]
}))
print(wl[, .(heat, delta_pool = round(delta_pool,6), main = round(main_pct,1),
  coast = round(coast_pct,1), coast_t = round(coast_t,2),
  heat = round(heat_pct,1), forest = round(forest_pct,1),
  chi2 = round(wald,1), p = signif(wald_p,3), neg = n_neg_w,
  err = signif(identity_err,2))], row.names = FALSE)
cs <- unique(s[, .(city_id, obs_n_days, coast)])
cat(sprintf("\n窓長(obs_n_days) と coast の市町村断面相関 = %+.3f\n",
            cor(cs$obs_n_days, cs$coast)))
cat(sprintf("c40_HD35 と obs_HD35 の相関 = %.3f（窓長の変動は 35℃日が無い暦日に乗る）\n",
            cor(s$c40_HD35, s$obs_HD35)))
save(wl, "t22_window_ladder")

cat("\n出力: research/model/output/t22_*.csv\n")

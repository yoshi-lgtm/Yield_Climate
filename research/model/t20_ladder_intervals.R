# =============================================================================
# T20: 統制を順に追加した表と、区間推定（クラスター／Conley／IM）
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t20_ladder_intervals.R
#
# 出力 : research/model/output/t20_*.csv と標準出力
#
# 論文の軸「恒等式 → 重み集中 → 名前のある汚染 → 残ったもの」に対応して
#   A. 統制を順に追加した表（経路依存の開示）
#   B. delta_2010 の区間推定（informative null を区間で示す）
#   C. 汚染項の推論を Conley と IM の両方で
# を出す。
#
# 追加の順序は事前の機序で固定する（統計量を見て並べ替えない）。
#   R0 素のプール（T11）
#   R1 ＋年×consol30 の水準   … 年固有ショック（2010年の政策ショックを含む）
#   R2 ＋heat×coast を年飽和   … 海洋性の緩和
#   R3 ＋heat×forest を年飽和  … 局所の遮蔽・冷気流出  ← 本番
#   R4 ＋heat×elev（参考）     … 機序が弱い
#   R5 ＋heat×lat,lon（参考）  … 機序がない
#
# C の IM について：本番スペックを都県別に年飽和で回すと神奈川で41個の係数が
# 共線で落ちる（T19 §1.1）。汚染項は gamma_t の時間変動から出るので年別の
# 構造自体は捨てられない。そこで年を4期（6年ずつ）に粗視化した分割で
# 都県別推定を可能にし、まず全標本で粗視化が物語を保つことを確認してから
# IM に乗せる。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

pk <- function(x) x[file.exists(x)][1]
PANEL <- pk(c("research/model/kanto_panel_2000_2023.csv",
              "Yield_Climate/research/model/kanto_panel_2000_2023.csv"))
CEN <- pk(c("results/city_centroids.csv", "Yield_Climate/results/city_centroids.csv"))
CD  <- pk(c("results/city_coastdist.csv", "Yield_Climate/results/city_coastdist.csv"))
FR  <- pk(c("results/city_forest.csv", "Yield_Climate/results/city_forest.csv"))
EL  <- pk(c("results/city_elevation.csv", "Yield_Climate/results/city_elevation.csv"))
OUT <- "research/model/output"
CTRL <- c("GSR_78", "APCP_78")
HEAT <- "obs_HD35"
CONLEY <- 50
EVT <- 2010L
PREF_NAME <- c("8" = "茨城", "9" = "栃木", "10" = "群馬", "11" = "埼玉",
               "12" = "千葉", "14" = "神奈川")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sv <- function(d, n) fwrite(d, file.path(OUT, paste0(n, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d <- merge(d, fread(CEN, encoding = "UTF-8")[, .(city_id, lon, lat)], by = "city_id")
d <- merge(d, fread(CD, encoding = "UTF-8")[, .(city_id, coast = coast_km)], by = "city_id")
d <- merge(d, fread(FR, encoding = "UTF-8")[, .(city_id, forest = forest_nbr)], by = "city_id")
d <- merge(d, fread(EL, encoding = "UTF-8")[, .(city_id, elev)], by = "city_id")

vv <- c("log_yield", HEAT, "consol30", "coast", "forest", "elev", "lon", "lat",
        CTRL, "city_id", "year")
s <- d[complete.cases(d[, ..vv])]
s[, pref := city_id %/% 1000L]
years <- sort(unique(s$year))
cs <- unique(s[, c("city_id", c("coast", "forest", "elev", "lon", "lat")), with = FALSE])
for (v in c("coast", "forest", "elev", "lon", "lat")) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = v, value = (s[[v]] - mu) / sg)
}
s[, D := get(HEAT) * consol30]
for (v in c("coast", "forest", "elev", "lon", "lat"))
  set(s, j = paste0("h_", v), value = s[[HEAT]] * s[[v]])

# dat を受け取って水準列を作る（環境経由の副作用を避ける）
mk_lev <- function(dat, v, tag, keys) {
  nm <- character(0)
  for (k in keys[-1]) {
    n <- sprintf("L%s_%s", tag, k)
    set(dat, j = n, value = dat[[v]] * as.integer(dat$gkey == k))
    nm <- c(nm, n)
  }
  nm
}

# -----------------------------------------------------------------------------
# 分解の中核。gvar = 年または期のキー列名
# -----------------------------------------------------------------------------
decomp <- function(label, slopes, levvars, gvar = "year", dat = s, vcvs = TRUE) {
  dat <- copy(dat)
  dat[, gkey := as.character(get(gvar))]
  keys <- sort(unique(dat$gkey))
  lv <- character(0)
  for (v in levvars) lv <- c(lv, mk_lev(dat, v, v, keys))
  # 群内で定数になった水準列は落とす（都県別推定で起きる）
  keep <- vapply(lv, function(z) stats::sd(dat[[z]]) > 0, logical(1))
  lv <- lv[keep]

  rp <- c("D", slopes, lv, CTRL)
  fp <- as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rp, collapse = " + ")))
  mp <- tryCatch(feols(fp, dat, cluster = ~city_id, fixef.rm = "none",
                       ssc = ssc(K.fixef = "full"), warn = FALSE, notes = FALSE),
                 error = function(e) NULL)
  if (is.null(mp)) return(NULL)
  dp <- unname(coef(mp)["D"])

  rs <- c("i(gkey, D)", sprintf("i(gkey, %s)", slopes), lv, CTRL)
  fs <- as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rs, collapse = " + ")))
  ms <- tryCatch(feols(fs, dat, cluster = ~city_id, fixef.rm = "none",
                       ssc = ssc(K.fixef = "full"), warn = FALSE, notes = FALSE),
                 error = function(e) NULL)
  if (is.null(ms)) return(NULL)

  gx <- function(m, v) {
    p <- strsplit(names(coef(m)), ":", fixed = TRUE)
    i <- which(vapply(p, function(z) length(z) >= 4 && z[1] == "gkey" && z[4] == v,
                      logical(1)))
    i[order(vapply(p[i], function(z) z[3], character(1)))]
  }
  idD <- gx(ms, "D")
  if (length(idD) != length(keys)) return(NULL)

  pr <- c("D", slopes, lv, CTRL)
  M <- as.data.frame(demean(dat[, ..pr], dat[, .(city_id, year)]))
  f4 <- as.formula(paste("D ~", paste(c(slopes, lv, CTRL), collapse = " + "), "- 1"))
  Dt <- as.numeric(residuals(lm(f4, data = M)))
  ss <- sum(Dt^2)
  set(dat, j = "Dtil", value = Dt)   # 列にしてから group by する
  wD <- dat[, .(w = sum(Dtil * D) / ss), by = gkey][order(gkey)]$w
  bh <- coef(ms)
  main <- sum(wD * bh[idD])
  ct <- lapply(slopes, function(v) {
    iv <- gx(ms, v)
    wv <- dat[, .(w = sum(Dtil * get(v)) / ss), by = gkey][order(gkey)]$w
    # 共線で係数が落ちると長さが合わない。その場合は分解が成立しないので NA
    if (length(iv) != length(wv))
      return(list(v = v, w = wv, idx = iv, val = NA_real_))
    list(v = v, w = wv, idx = iv, val = sum(wv * bh[iv]))
  })
  names(ct) <- slopes
  tot <- main + sum(vapply(ct, function(z) z$val, numeric(1)))

  i10 <- which(keys == as.character(EVT))
  out <- list(label = label, dpool = dp, main = main, ct = ct, err = dp - tot,
              keys = keys, wD = wD, bD = bh[idD], idD = idD, ms = ms, i10 = i10,
              n_neg = sum(wD < 0))
  if (!vcvs) return(out)

  for (vl in c("cluster", sprintf("conley%d", CONLEY))) {
    V <- if (vl == "cluster") vcov(ms) else
      tryCatch(suppressWarnings(vcov_conley(ms, lat = ~lat, lon = ~lon,
                                            cutoff = CONLEY, ssc = ssc(K.fixef = "full"))),
               error = function(e) NULL)
    if (is.null(V)) next
    Vd <- V[idD, idD]
    sed <- sqrt(pmax(diag(Vd), 0))
    K <- length(keys)
    R <- matrix(0, K - 1, K); for (i in 1:(K - 1)) { R[i, i] <- 1; R[i, i + 1] <- -1 }
    Rb <- R %*% bh[idD]; RV <- R %*% Vd %*% t(R)
    sol <- tryCatch(solve(RV), error = function(e) NULL)
    st <- if (is.null(sol)) NA_real_ else as.numeric(t(Rb) %*% sol %*% Rb)
    out[[vl]] <- list(
      se_d = sed, se_main = sqrt(max(as.numeric(t(wD) %*% Vd %*% wD), 0)),
      se_ct = vapply(ct, function(z)
        sqrt(max(as.numeric(t(z$w) %*% V[z$idx, z$idx] %*% z$w), 0)), numeric(1)),
      wald = st, wald_p = if (is.na(st)) NA_real_ else
        pchisq(st, K - 1, lower.tail = FALSE),
      sd_ratio = sd(bh[idD]) / mean(sed))
  }
  out
}

ci <- function(b, se) sprintf("[%+.5f, %+.5f]", b - 1.96 * se, b + 1.96 * se)

# -----------------------------------------------------------------------------
# A. 統制を順に追加した表
# -----------------------------------------------------------------------------
RUNGS <- list(
  list("R0 素のプール（T11）", HEAT, character(0)),
  list("R1 ＋年×consol30", HEAT, "consol30"),
  list("R2 ＋coast（海洋性）", c(HEAT, "h_coast"), c("consol30", "coast")),
  list("R3 ＋forest（本番）", c(HEAT, "h_coast", "h_forest"),
       c("consol30", "coast", "forest")),
  list("R4 参考 ＋elev", c(HEAT, "h_coast", "h_forest", "h_elev"),
       c("consol30", "coast", "forest", "elev")),
  list("R5 参考 ＋lat,lon", c(HEAT, "h_coast", "h_forest", "h_lat", "h_lon"),
       c("consol30", "coast", "forest", "lat", "lon")))

cat("################ A. 統制を順に追加した表 ################\n")
res <- list(); rows <- list(); crows <- list()
for (r in RUNGS) {
  o <- decomp(r[[1]], r[[2]], r[[3]])
  if (is.null(o)) { cat(sprintf("%-26s 推定できず\n", r[[1]])); next }
  res[[r[[1]]]] <- o
  cl <- o$cluster
  cat(sprintf("\n--- %s ---\n", o$label))
  cat(sprintf("プール delta %+.6f  恒等式誤差 %.1e  負のomega %d/%d  omega_2010 %.4f\n",
              o$dpool, o$err, o$n_neg, length(o$keys), o$wD[o$i10]))
  cat(sprintf("  主項 %+.6f (SE %.6f, t %+.2f) %s\n", o$main, cl$se_main,
              o$main / cl$se_main, ci(o$main, cl$se_main)))
  for (v in names(o$ct)) {
    cat(sprintf("  汚染(%-8s) %+.6f (SE %.6f, t %+.2f) %s  プール値比 %+.1f%%\n",
                v, o$ct[[v]]$val, cl$se_ct[[v]], o$ct[[v]]$val / cl$se_ct[[v]],
                ci(o$ct[[v]]$val, cl$se_ct[[v]]), 100 * o$ct[[v]]$val / o$dpool))
  }
  cat(sprintf("  delta_2010 %+.6f (SE %.6f, t %+.2f) %s\n", o$bD[o$i10],
              cl$se_d[o$i10], o$bD[o$i10] / cl$se_d[o$i10],
              ci(o$bD[o$i10], cl$se_d[o$i10])))
  cat(sprintf("  chi2(%d) %.1f (p %.3g)   sd(delta_t)/平均SE %.2f\n",
              length(o$keys) - 1, cl$wald, cl$wald_p, cl$sd_ratio))
  cn <- o[[sprintf("conley%d", CONLEY)]]
  rows[[length(rows) + 1]] <- data.table(
    rung = o$label, delta_pool = o$dpool, main = o$main, se_main = cl$se_main,
    omega_2010 = o$wD[o$i10], share_2010 = 100 * o$wD[o$i10] * o$bD[o$i10] / o$dpool,
    n_neg = o$n_neg, delta_2010 = o$bD[o$i10], se_2010_cl = cl$se_d[o$i10],
    se_2010_cn = if (is.null(cn)) NA_real_ else cn$se_d[o$i10],
    wald = cl$wald, wald_p = cl$wald_p, sd_ratio = cl$sd_ratio,
    contam_total = sum(vapply(o$ct, function(z) z$val, numeric(1))),
    identity_err = o$err)
  for (v in names(o$ct)) crows[[length(crows) + 1]] <- data.table(
    rung = o$label, source = v, value = o$ct[[v]]$val,
    se_cluster = cl$se_ct[[v]],
    se_conley = if (is.null(cn)) NA_real_ else cn$se_ct[[v]],
    pct = 100 * o$ct[[v]]$val / o$dpool)
}
lad <- rbindlist(rows); con <- rbindlist(crows)
sv(lad, "t20_ladder"); sv(con, "t20_contamination")

# -----------------------------------------------------------------------------
# B/C. IM（都県別）— 本番スペックと、汚染項のための粗視化
# -----------------------------------------------------------------------------
s[, period := cut(year, breaks = c(1999, 2005, 2011, 2017, 2023),
                  labels = c("p1", "p2", "p3", "p4"))]
cat("\n\n################ C-0. 粗視化（4期）が物語を保つか ################\n")
o3 <- res[["R3 ＋forest（本番）"]]
o3c <- decomp("R3 粗視化（4期）", c(HEAT, "h_coast", "h_forest"),
              c("consol30", "coast", "forest"), gvar = "period")
cat(sprintf("年別   : プール %+.6f  主項 %+.6f  汚染(coast) %+.6f (%.1f%%)\n",
            o3$dpool, o3$main, o3$ct[["h_coast"]]$val,
            100 * o3$ct[["h_coast"]]$val / o3$dpool))
cat(sprintf("4期    : プール %+.6f  主項 %+.6f  汚染(coast) %+.6f (%.1f%%)\n",
            o3c$dpool, o3c$main, o3c$ct[["h_coast"]]$val,
            100 * o3c$ct[["h_coast"]]$val / o3c$dpool))
cat(sprintf("恒等式誤差 年別 %.1e / 4期 %.1e\n", o3$err, o3c$err))

cat("\n################ B,C. IM（G=6） ################\n")
prefs <- sort(unique(s$pref))
# delta_2010 の IM は倹約版（1[2010] との交差、共線ゼロ。T19 §1）を使う。
# 年飽和を都県別に回すと神奈川で41個が落ちて分解が成立しない。
im_d2010 <- function(p) {
  sp <- s[pref == p]
  sp[, I10 := as.integer(year == EVT)]
  f <- as.formula(sprintf(
    "log_yield ~ %s + D + h_coast + h_forest + I10:%s + I10:D + I10:h_coast + I10:h_forest + %s | city_id + year",
    HEAT, HEAT, paste(CTRL, collapse = " + ")))
  tryCatch({
    m <- feols(f, sp, fixef.rm = "none", ssc = ssc(K.fixef = "full"),
               warn = FALSE, notes = FALSE)
    cf <- coef(m)
    a <- cf[intersect(c("D:I10", "I10:D"), names(cf))[1]]
    if (is.na(a) || is.na(cf["D"])) NA_real_ else unname(cf["D"] + a)
  }, error = function(e) NA_real_)
}

im_one <- function(p, gvar, what) {
  if (what == "d2010") return(im_d2010(p))
  sp <- s[pref == p]
  o <- tryCatch(decomp("p", c(HEAT, "h_coast", "h_forest"),
                       c("consol30", "coast", "forest"), gvar = gvar,
                       dat = sp, vcvs = FALSE), error = function(e) NULL)
  if (is.null(o)) return(NA_real_)
  if (what == "d2010") {
    i <- which(o$keys == as.character(EVT))
    if (!length(i)) NA_real_ else unname(o$bD[i])
  } else if (what == "main") o$main
  else o$ct[["h_coast"]]$val
}
imt <- rbindlist(lapply(
  list(c("year", "d2010", "delta_2010（年別）"),
       c("period", "main", "主項（4期）"),
       c("period", "contam", "汚染 coast（4期）")),
  function(z) data.table(target = z[3], pref = prefs,
                         est = vapply(prefs, im_one, numeric(1),
                                      gvar = z[1], what = z[2]))))
imt[, name := PREF_NAME[as.character(pref)]]
print(dcast(imt, pref + name ~ target, value.var = "est")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 6) else x)], row.names = FALSE)

imres <- imt[, {
  b <- est[is.finite(est)]; G <- length(b)
  if (G < 2) .(G = G, mean = NA_real_, t = NA_real_, p = NA_real_,
               lo = NA_real_, hi = NA_real_)
  else {
    se <- sd(b) / sqrt(G); q <- qt(0.975, G - 1)
    .(G = G, mean = mean(b), t = mean(b) / se, p = 2 * pt(-abs(mean(b) / se), G - 1),
      lo = mean(b) - q * se, hi = mean(b) + q * se)
  }
}, by = target]
cat("\n##### IM 検定と区間 #####\n")
print(imres[, .(target, G, mean = round(mean, 6), t = round(t, 2),
                p = round(p, 3),
                CI95 = sprintf("[%+.5f, %+.5f]", lo, hi))], row.names = FALSE)
sv(imt, "t20_im_by_pref"); sv(imres, "t20_im_test")

# -----------------------------------------------------------------------------
# D. 区間推定のまとめ（本番スペック）
# -----------------------------------------------------------------------------
cat("\n################ D. 本番スペックの区間推定 ################\n")
cl <- o3$cluster; cn <- o3[[sprintf("conley%d", CONLEY)]]
iv <- rbindlist(list(
  data.table(target = "delta_2010", vcov = "市町村クラスター",
             est = o3$bD[o3$i10], se = cl$se_d[o3$i10]),
  data.table(target = "delta_2010", vcov = sprintf("Conley %dkm", CONLEY),
             est = o3$bD[o3$i10], se = cn$se_d[o3$i10]),
  data.table(target = "汚染 coast", vcov = "市町村クラスター",
             est = o3$ct[["h_coast"]]$val, se = cl$se_ct[["h_coast"]]),
  data.table(target = "汚染 coast", vcov = sprintf("Conley %dkm", CONLEY),
             est = o3$ct[["h_coast"]]$val, se = cn$se_ct[["h_coast"]]),
  data.table(target = "主項", vcov = "市町村クラスター",
             est = o3$main, se = cl$se_main),
  data.table(target = "主項", vcov = sprintf("Conley %dkm", CONLEY),
             est = o3$main, se = cn$se_main)))
iv[, `:=`(t = est / se, lo = est - 1.96 * se, hi = est + 1.96 * se)]
imadd <- imres[, .(target = c("delta_2010", "主項", "汚染 coast")[
  match(target, c("delta_2010（年別）", "主項（4期）", "汚染 coast（4期）"))],
  vcov = "IM (G=6)", est = mean, se = NA_real_, t = t, lo, hi)]
iv <- rbindlist(list(iv, imadd), use.names = TRUE)
print(iv[, .(target, vcov, est = round(est, 6), t = round(t, 2),
             CI95 = sprintf("[%+.5f, %+.5f]", lo, hi))], row.names = FALSE)
sv(iv, "t20_intervals")
cat("\n出力: research/model/output/t20_{ladder,contamination,im_by_pref,im_test,intervals}.csv\n")

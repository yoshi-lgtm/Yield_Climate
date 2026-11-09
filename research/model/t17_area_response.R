# =============================================================================
# T17: 2010年の作付面積反応
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t17_area_response.R
#
# 出力 : research/model/output/t17_*.csv と標準出力
#
# なぜ要るか。T15 で、2010年の戸別所得補償モデル事業が収量に効く向きが
# 理論的に決まらないことが分かった：
#   (a) 過剰作付の解消            → 面積↓ 10a収量↑
#   (b) 交付は作付面積あたり定額  → 面積↑ 10a収量↓
#   (c) 転作交付金で良田が流出    → 面積↓ 10a収量↓
# T15 の 年×consol30 飽和は水準差を吸収するだけでチャネルを特定しない。
# そこで面積そのものを左辺に置いて直接見る。
#
# 問いは二つで、別のものである。
#   問1（チャネル）      面積は consol30 に反応したか。符号が (a)(b)(c) を分ける。
#   問2（推定量への脅威）本番と同じ飽和スペックで左辺を ln(面積) にしたとき、
#                        delta^area_2010 が立つか。delta_2010 は 2010年の
#                        heat×consol30 の係数なので、2010年の面積再配分が
#                        heat×consol30 と相関して delta_2010 を作っているなら、
#                        ここに同じ山が出るはずである。出なければ脅威は閉じる。
#
# データ注記。作付面積は 04yield_data_kanto_debuged.py が 167 行目で
# 「収量」を含む系列だけに絞っていたため kanto_yield.csv に入っていなかった
# （「収穫量」は部分文字列「収量」を含まないので残るのは 10a 当たり収量だけ）。
# 05area_data_kanto.py で同じ接合ロジックのまま作付面積を抜き直した。
# 10a当たり収量と作付面積は作物統計で別系列として公表されるので、
# 面積は被説明変数 log_yield に機械的には入っていない。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

pk <- function(x) x[file.exists(x)][1]
PANEL <- pk(c("research/model/kanto_panel_2000_2023.csv",
              "Yield_Climate/research/model/kanto_panel_2000_2023.csv"))
AREA <- pk(c("yield/kanto_area.csv", "Yield_Climate/yield/kanto_area.csv"))
CEN <- pk(c("results/city_centroids.csv", "Yield_Climate/results/city_centroids.csv"))
CD  <- pk(c("results/city_coastdist.csv", "Yield_Climate/results/city_coastdist.csv"))
FR  <- pk(c("results/city_forest.csv", "Yield_Climate/results/city_forest.csv"))
EL  <- pk(c("results/city_elevation.csv", "Yield_Climate/results/city_elevation.csv"))
OUT <- "research/model/output"
CTRL <- c("GSR_78", "APCP_78")
HEAT <- "obs_HD35"
CONLEY <- 50
EVT <- 2010L
BASE <- 2009L   # 長期差分の基準年

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sv <- function(d, n) fwrite(d, file.path(OUT, paste0(n, ".csv")), bom = TRUE)

# -----------------------------------------------------------------------------
# 0. データ。T20 と同じ本番サンプルの作り方に、作付面積を足す
# -----------------------------------------------------------------------------
d <- fread(PANEL, encoding = "UTF-8")
d <- merge(d, fread(CEN, encoding = "UTF-8")[, .(city_id, lon, lat)], by = "city_id")
d <- merge(d, fread(CD, encoding = "UTF-8")[, .(city_id, coast = coast_km)], by = "city_id")
d <- merge(d, fread(FR, encoding = "UTF-8")[, .(city_id, forest = forest_nbr)], by = "city_id")
d <- merge(d, fread(EL, encoding = "UTF-8")[, .(city_id, elev)], by = "city_id")

a <- fread(AREA, encoding = "UTF-8")
setnames(a, c("時間軸（年次）", "value"), c("year", "area"))
a <- unique(a[, .(year = as.integer(year), prefecture, city_name,
                  area = suppressWarnings(as.numeric(area)))])
n_pre <- nrow(d)
d <- merge(d, a, by = c("prefecture", "city_name", "year"), all.x = TRUE)
stopifnot(nrow(d) == n_pre)   # 面積側に重複キーが無いことの確認
d[, log_area := log(area)]

vv <- c("log_yield", HEAT, "consol30", "coast", "forest", "elev", "lon", "lat",
        CTRL, "city_id", "year")
s <- d[complete.cases(d[, ..vv])]
cat(sprintf("本番サンプル %d行 %d市町村。うち作付面積が取れるのは %d行（欠損 %d行）\n",
            nrow(s), uniqueN(s$city_id), sum(!is.na(s$area)), sum(is.na(s$area))))
cat(sprintf("consol30 は市町村内で一定か: %s\n",
            all(s[, .(k = uniqueN(consol30)), by = city_id]$k == 1)))
c30 <- unique(s[, .(city_id, consol30)])$consol30
cat(sprintf("consol30: 連続値 %d通り, 中位 %.3f, 10-90%%分位 %.3f-%.3f, ==0 が %d市, ==1 が %d市\n",
            uniqueN(c30), median(c30), quantile(c30, 0.1), quantile(c30, 0.9),
            sum(c30 == 0), sum(c30 == 1)))

s[, pref := city_id %/% 1000L]
s[, `:=`(lat_deg = lat, lon_deg = lon)]
cs <- unique(s[, c("city_id", "coast", "forest", "elev", "lon", "lat"), with = FALSE])
for (v in c("coast", "forest", "elev", "lon", "lat")) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = v, value = (s[[v]] - mu) / sg)
}
s[, D := get(HEAT) * consol30]
for (v in c("coast", "forest"))
  set(s, j = paste0("h_", v), value = s[[HEAT]] * s[[v]])

# -----------------------------------------------------------------------------
# 1. 問1：面積は consol30 に反応したか（長期差分のイベント時間）
#    Delta_t = ln A_it - ln A_i,2009 を consol30 に回す。横断面なので
#    SE は HC1 と Conley 50km の両方を出す（市町村クラスターは1obs/市で無意味）。
# -----------------------------------------------------------------------------
cat("\n################ 1. 問1 チャネル: Delta ln(面積) を consol30 に回す ################\n")
# consol30 は連続変数（0.2ha 以上の区画のうち大区画・ほ区均平が占める割合、
# 中位 0.286）。係数は consol30 を 0 から 1 に動かしたときの ln(面積) の差である。
cat(sprintf("基準年 %d。係数は consol30 を 0→1 に動かしたときの ln(作付面積) の差\n", BASE))

w <- dcast(s[year %between% c(2004, 2014) & !is.na(log_area)],
           city_id ~ year, value.var = "log_area")
cv <- unique(s[, .(city_id, consol30, coast, forest, lat_deg, lon_deg)])
w <- merge(w, cv, by = "city_id")
bcol <- as.character(BASE)

rows <- list()
for (y in 2004:2014) {
  if (y == BASE) next
  yc <- as.character(y)
  ww <- w[!is.na(get(yc)) & !is.na(get(bcol))]
  ww[, dy := get(yc) - get(bcol)]
  for (sp in c("素", "＋coast,forest")) {
    f <- if (sp == "素") dy ~ consol30 else dy ~ consol30 + coast + forest
    m <- feols(f, ww, vcov = "hetero", warn = FALSE, notes = FALSE)
    vc <- tryCatch(suppressWarnings(vcov_conley(m, lat = ~lat_deg, lon = ~lon_deg,
                                                cutoff = CONLEY)),
                   error = function(e) NULL)
    b <- unname(coef(m)["consol30"]); se1 <- unname(se(m)["consol30"])
    se2 <- if (is.null(vc)) NA_real_ else sqrt(max(vc["consol30", "consol30"], 0))
    rows[[length(rows) + 1]] <- data.table(
      year = y, spec = sp, n = nobs(m), b = b, se_hc1 = se1, se_conley = se2,
      t_hc1 = b / se1, t_conley = b / se2)
  }
}
ev <- rbindlist(rows)
sv(ev, "t17_area_eventtime")
cat("\n")
print(ev[, .(year, spec, n, b = round(b, 4), se_hc1 = round(se_hc1, 4),
             t_hc1 = round(t_hc1, 2), t_conley = round(t_conley, 2))])

cat("\n--- 2010年だけ、基準の取り方を変えて（統制は coast,forest 入り）---\n")
rb <- list()
for (bs in list(c(2009, 2009), c(2007, 2009), c(2005, 2009))) {
  ww <- dcast(s[year %between% c(bs[1], 2010) & !is.na(log_area)],
              city_id ~ year, value.var = "log_area")
  bb <- as.character(bs[1]:bs[2])
  ww <- ww[complete.cases(ww[, c("2010", bb), with = FALSE])]
  ww[, dy := `2010` - rowMeans(.SD), .SDcols = bb]
  ww <- merge(ww, cv, by = "city_id")
  m <- feols(dy ~ consol30 + coast + forest, ww, vcov = "hetero",
             warn = FALSE, notes = FALSE)
  rb[[length(rb) + 1]] <- data.table(
    base = sprintf("%d-%d の平均", bs[1], bs[2]), n = nobs(m),
    b = unname(coef(m)["consol30"]), se = unname(se(m)["consol30"]),
    t = unname(coef(m)["consol30"] / se(m)["consol30"]))
}
rb <- rbindlist(rb); sv(rb, "t17_area_2010_bases")
print(rb[, .(base, n, b = round(b, 4), se = round(se, 4), t = round(t, 2))])

# -----------------------------------------------------------------------------
# 2. 問2：本番と同じ飽和スペックで、左辺を ln(面積) にする
#    delta^_t = 年 t の heat x consol30 の係数。2010年に山が出るかを見る。
# -----------------------------------------------------------------------------
cat("\n\n################ 2. 問2 脅威: 本番スペックの左辺を ln(面積) に差し替え ################\n")

keys <- sort(unique(as.character(s$year)))
mk_lev <- function(dat) {
  for (v in c("consol30", "coast", "forest")) for (k in keys[-1])
    set(dat, j = sprintf("L%s_%s", v, k),
        value = dat[[v]] * as.integer(as.character(dat$year) == k))
  nm <- grep("^L(consol30|coast|forest)_", names(dat), value = TRUE)
  nm[vapply(nm, function(z) stats::sd(dat[[z]]) > 0, logical(1))]
}
fit_sat <- function(lhs, dat, lv) {
  rs <- c("i(year, D)", sprintf("i(year, %s)", c(HEAT, "h_coast", "h_forest")),
          lv, CTRL)
  f <- as.formula(sprintf("%s ~ %s | city_id + year", lhs, paste(rs, collapse = " + ")))
  feols(f, dat, cluster = ~city_id, fixef.rm = "none",
        ssc = ssc(K.fixef = "full"), warn = FALSE, notes = FALSE)
}
idx <- function(m) {
  p <- strsplit(names(coef(m)), ":", fixed = TRUE)
  i <- which(vapply(p, function(z) length(z) >= 4 && z[1] == "year" && z[4] == "D",
                    logical(1)))
  i[order(vapply(p[i], function(z) z[3], character(1)))]
}
grab <- function(lhs, tag, dat, lv) {
  m <- fit_sat(lhs, dat, lv); i <- idx(m)
  Vc <- tryCatch(suppressWarnings(
    vcov_conley(m, lat = ~lat_deg, lon = ~lon_deg, cutoff = CONLEY,
                ssc = ssc(K.fixef = "full"))), error = function(e) NULL)
  data.table(lhs = tag, year = as.integer(keys), b = unname(coef(m)[i]),
             se_cl = sqrt(diag(vcov(m))[i]),
             se_cn = if (is.null(Vc)) NA_real_ else sqrt(pmax(diag(Vc)[i], 0)))
}

sa <- s[!is.na(log_area)]
lv_a <- mk_lev(sa)
sat <- rbind(grab("log_area", "ln(作付面積)", sa, lv_a),
             grab("log_yield", "ln(10a収量) 参照", sa, lv_a))
sat[, `:=`(t_cl = b / se_cl, t_cn = b / se_cn)]
sv(sat, "t17_area_saturated")
cat(sprintf("\n標本 %d行 %d市町村（面積欠損を落としたため本番より %d行少ない）\n",
            nrow(sa), uniqueN(sa$city_id), nrow(s) - nrow(sa)))
for (tg in unique(sat$lhs)) {
  cat(sprintf("\n--- 左辺 %s : delta^_t（年 t の heat x consol30 の係数）---\n", tg))
  x <- sat[lhs == tg]
  print(x[, .(year, b = round(b, 5), se_cl = round(se_cl, 5),
              t_cl = round(t_cl, 2), t_cn = round(t_cn, 2))])
  x10 <- x[year == EVT]
  cat(sprintf("  2010年: %+.5f (クラスターSE %.5f, t %+.2f / Conley t %+.2f) 95%%CI [%+.5f, %+.5f]\n",
              x10$b, x10$se_cl, x10$t_cl, x10$t_cn,
              x10$b - 1.96 * x10$se_cl, x10$b + 1.96 * x10$se_cl))
  cat(sprintf("  クラスターSEで |t|>1.96 の年: %s\n",
              paste(x[abs(t_cl) > 1.96]$year, collapse = ", ")))
}

# 面積欠損を落としたことで収量側の delta_2010 が動いていないかの確認
cat("\n--- 標本縮小の確認（左辺 ln(10a収量)）---\n")
lv_f <- mk_lev(s)
mf <- fit_sat("log_yield", s, lv_f); i <- idx(mf)
j <- which(as.integer(keys) == EVT)
bf <- unname(coef(mf)[i])[j]; sf <- sqrt(diag(vcov(mf))[i])[j]
cat(sprintf("  本番標本(%d行)     delta_2010 %+.6f (SE %.6f, t %+.2f)\n",
            nrow(s), bf, sf, bf / sf))
x <- sat[lhs == "ln(10a収量) 参照" & year == EVT]
cat(sprintf("  面積ありのみ(%d行) delta_2010 %+.6f (SE %.6f, t %+.2f)\n",
            nrow(sa), x$b, x$se_cl, x$t_cl))

# -----------------------------------------------------------------------------
# 3. 面積チャネルで delta_2010 を説明するのに必要な弾性値
#    delta_2010 = eps * delta^area_2010 となる eps = d ln(10a収量)/d ln(面積)。
#    面積が減れば残る田の平均質が上がる（限界地の脱落）ので eps<0 が想定される。
#    必要な |eps| が大きすぎれば、面積チャネルでは説明がつかない。
# -----------------------------------------------------------------------------
cat("

################ 3. 面積チャネルで delta_2010 を説明するのに必要な弾性値 ################
")
ya <- sat[lhs == "ln(作付面積)" & year == EVT]
yy <- sat[lhs == "ln(10a収量) 参照" & year == EVT]
lo <- ya$b - 1.96 * ya$se_cl; hi <- ya$b + 1.96 * ya$se_cl
cat(sprintf("delta_2010（収量）      %+.5f
", yy$b))
cat(sprintf("delta^area_2010（面積） %+.5f  95%%CI [%+.5f, %+.5f]
", ya$b, lo, hi))
cat(sprintf("点推定で全部を説明するのに必要な eps  %+.2f
", yy$b / ya$b))
cat(sprintf("CIの最も有利な端(%+.5f)でも必要な eps %+.2f
", lo, yy$b / lo))
sv(data.table(delta_2010 = yy$b, delta_area_2010 = ya$b,
              area_ci_lo = lo, area_ci_hi = hi,
              eps_needed_point = yy$b / ya$b, eps_needed_bound = yy$b / lo),
   "t17_elasticity_bound")

# =============================================================================
# T19: Ibragimov-Muller を本番スペックで
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t19_im_mechanism.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_coastdist.csv / city_forest.csv
# 出力 : research/model/output/t19_*.csv と標準出力
#
# IM は都県ごとに独立に推定した G=6 個の値に t 検定を当てる。仮定は
# 「都県内は任意の従属、都県間は独立、群内で漸近正規」で、Conley より弱い。
#
# 本番スペック（T18）をそのまま都県別に走らせると共線で潰れる。T12 §4.1 で
# 既に起きていた問題が、年×水準項が 3x23=69 本に増えてさらに悪化する
# （神奈川は8市町村・177行しかない）。したがって群内の式は推定可能な範囲まで
# 落とす必要がある。3通りを並べて劣化の度合いを示す。
#
#   A. プール delta（本番コントロール入り・年飽和なし）
#      log_yield ~ heat + D + G + F + CTRL | city + year        … 5 パラメータ
#   B. delta_2010（倹約版。2010年ダミーとの交差で偏差を取る）
#      ... + I10:(heat + D + G + F)                             … 9 パラメータ
#   C. delta_2010（本番と同じ年飽和。共線で落ちる係数を報告する）
#
# あわせて T12 §4 の「群馬が突出した外れ値」が地理の交絡だったかを確かめる。
# 群馬は内陸・森林が多い（海岸距離が最大・森林率が高い）ので、coast と forest を
# 統制すると群馬の値が縮むなら、あの外れ値は地理の押し込みだったことになる。
# =============================================================================

suppressMessages({
  library(fixest)
  library(data.table)
})

pick <- function(c_) c_[file.exists(c_)][1]
PANEL <- pick(c("research/model/kanto_panel_2000_2023.csv",
                "Yield_Climate/research/model/kanto_panel_2000_2023.csv"))
CDIST <- pick(c("results/city_coastdist.csv", "Yield_Climate/results/city_coastdist.csv"))
FRST  <- pick(c("results/city_forest.csv", "Yield_Climate/results/city_forest.csv"))
OUT_DIR <- "research/model/output"
CTRL <- c("GSR_78", "APCP_78")
HEAT <- "obs_HD35"
EVT  <- 2010L
PREF_NAME <- c("8" = "茨城", "9" = "栃木", "10" = "群馬", "11" = "埼玉",
               "12" = "千葉", "14" = "神奈川")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, n) fwrite(dt, file.path(OUT_DIR, paste0(n, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
d <- merge(d, fread(CDIST, encoding = "UTF-8")[, .(city_id, coast = coast_km)],
           by = "city_id")
d <- merge(d, fread(FRST, encoding = "UTF-8")[, .(city_id, forest = forest_nbr)],
           by = "city_id")
vars <- c("log_yield", HEAT, "consol30", "coast", "forest", CTRL, "city_id", "year")
s <- d[complete.cases(d[, ..vars])]
s[, pref := city_id %/% 1000L]

# 標準化は全標本で行う（変数の定義を都県間で固定するため）
cs <- unique(s[, .(city_id, pref, coast, forest, consol30)])
for (v in c("coast", "forest")) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = v, value = (s[[v]] - mu) / sg)
}
s[, `:=`(D = get(HEAT) * consol30, G = get(HEAT) * coast, F_ = get(HEAT) * forest)]
s[, I10 := as.integer(year == EVT)]
years <- sort(unique(s$year))

prefs <- sort(unique(s$pref))
cat(sprintf("n = %d, 都県 = %d\n\n", nrow(s), length(prefs)))

# -----------------------------------------------------------------------------
# 0. 都県の地理的な位置づけ（群馬の外れ値を読むために先に出す）
# -----------------------------------------------------------------------------
cs2 <- unique(s[, .(city_id, pref, coast, forest, consol30)])
geo <- cs2[, .(cities = .N, coast = mean(coast), forest = mean(forest),
               consol30 = mean(consol30)), by = pref][order(pref)]
geo[, name := PREF_NAME[as.character(pref)]]
cat("##### 都県の地理（標準化後の平均。coast が大きい＝内陸） #####\n")
print(geo[, .(pref, name, cities, coast = round(coast, 2),
              forest = round(forest, 2), consol30 = round(consol30, 3))],
      row.names = FALSE)

# -----------------------------------------------------------------------------
# 1. 群内推定
# -----------------------------------------------------------------------------
f_A <- as.formula(sprintf("log_yield ~ %s + D + G + F_ + %s | city_id + year",
                          HEAT, paste(CTRL, collapse = " + ")))
f_B <- as.formula(sprintf(
  "log_yield ~ %s + D + G + F_ + I10:%s + I10:D + I10:G + I10:F_ + %s | city_id + year",
  HEAT, HEAT, paste(CTRL, collapse = " + ")))
f_C <- as.formula(sprintf(
  "log_yield ~ i(year, %s) + i(year, D) + i(year, G) + i(year, F_) + %s | city_id + year",
  HEAT, paste(CTRL, collapse = " + ")))
# 地理なし版（T12 と同じ）を対照に置く
f_A0 <- as.formula(sprintf("log_yield ~ %s + D + %s | city_id + year",
                           HEAT, paste(CTRL, collapse = " + ")))
f_B0 <- as.formula(sprintf(
  "log_yield ~ %s + D + I10:%s + I10:D + %s | city_id + year",
  HEAT, HEAT, paste(CTRL, collapse = " + ")))

est_one <- function(p, f, what) {
  sp <- s[pref == p]
  tryCatch({
    m <- feols(f, sp, fixef.rm = "none", ssc = ssc(K.fixef = "full"),
               warn = FALSE, notes = FALSE)
    cf <- coef(m)
    nd <- length(m$collin.var)
    if (what == "pool") {
      c(unname(cf["D"]), nd)
    } else if (what == "dev") {
      # fixest は交互作用を "D:I10" の順で命名する
      b <- unname(cf["D"])
      a <- unname(cf[intersect(c("D:I10", "I10:D"), names(cf))[1]])
      c(if (is.na(a) || is.na(b)) NA_real_ else b + a, nd)
    } else {
      pp <- strsplit(names(cf), ":", fixed = TRUE)
      ok <- vapply(pp, function(z) length(z) >= 4 &&
                     z[3] == as.character(EVT) && z[4] == "D", logical(1))
      c(if (!any(ok)) NA_real_ else unname(cf[ok][1]), nd)
    }
  }, error = function(e) c(NA_real_, NA_integer_))
}

grab <- function(f, what, lab) {
  r <- vapply(prefs, est_one, numeric(2), f = f, what = what)
  data.table(spec = lab, pref = prefs, name = PREF_NAME[as.character(prefs)],
             est = r[1, ], n_collin = as.integer(r[2, ]))
}

tab <- rbindlist(list(
  grab(f_A0, "pool", "A0 プール delta（地理なし＝T12）"),
  grab(f_A,  "pool", "A  プール delta（本番コントロール）"),
  grab(f_B0, "dev",  "B0 delta_2010 倹約版（地理なし）"),
  grab(f_B,  "dev",  "B  delta_2010 倹約版（本番コントロール）"),
  grab(f_C,  "sat",  "C  delta_2010 年飽和（本番と同じ）")))

# -----------------------------------------------------------------------------
# 2. IM 検定
# -----------------------------------------------------------------------------
im <- tab[, {
  b <- est[is.finite(est)]; G <- length(b)
  if (G < 2) .(G = G, mean = NA_real_, sd = NA_real_, t = NA_real_, p = NA_real_,
               crit = NA_real_)
  else {
    tt <- mean(b) / (sd(b) / sqrt(G))
    .(G = G, mean = mean(b), sd = sd(b), t = tt,
      p = 2 * pt(-abs(tt), df = G - 1), crit = qt(0.975, G - 1))
  }
}, by = spec]

cat("\n##### 都県別の推定値 #####\n")
w <- dcast(tab, pref + name ~ spec, value.var = "est")
print(w[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 5) else x)],
      row.names = FALSE)
cat("\n共線で落ちた係数の数:\n")
print(dcast(tab, pref + name ~ spec, value.var = "n_collin"), row.names = FALSE)

cat("\n##### IM 検定 #####\n")
print(im[, .(spec, G, mean = round(mean, 6), sd = round(sd, 6),
             t = round(t, 2), p = round(p, 3), crit = round(crit, 3))],
      row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. 群馬の外れ値は地理の交絡だったか
# -----------------------------------------------------------------------------
cat("\n##### 群馬（内陸・森林多）の値が地理統制で縮むか #####\n")
gu <- tab[pref == 10L]
for (pair in list(c("A0 プール delta（地理なし＝T12）", "A  プール delta（本番コントロール）"),
                  c("B0 delta_2010 倹約版（地理なし）", "B  delta_2010 倹約版（本番コントロール）"))) {
  a <- gu[spec == pair[1], est]; b <- gu[spec == pair[2], est]
  oa <- tab[spec == pair[1] & pref != 10L, est]
  ob <- tab[spec == pair[2] & pref != 10L, est]
  cat(sprintf("%-40s 群馬 %+.5f -> %+.5f （他5県の平均 %+.5f -> %+.5f）\n",
              sub(" .*", "", pair[1]), a, b, mean(oa, na.rm = TRUE),
              mean(ob, na.rm = TRUE)))
  cat(sprintf("%-40s 群馬が他5県の平均から離れている倍率 %.1f -> %.1f\n", "",
              abs(a - mean(oa, na.rm = TRUE)) / sd(oa, na.rm = TRUE),
              abs(b - mean(ob, na.rm = TRUE)) / sd(ob, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# 4. 1県ずつ除いた IM（どれか1県が結論を作っていないか）
#    神奈川は8市町村しかなく群内漸近正規が怪しいので特に見る
# -----------------------------------------------------------------------------
cat("
##### 1県除外の IM（倹約版 B・本番コントロール） #####
")
bb <- tab[spec == "B  delta_2010 倹約版（本番コントロール）"]
loo <- rbindlist(lapply(c(NA_integer_, prefs), function(dp) {
  v <- bb[is.na(dp) | pref != dp, est]
  v <- v[is.finite(v)]
  tt <- mean(v) / (sd(v) / sqrt(length(v)))
  data.table(dropped = if (is.na(dp)) "なし" else PREF_NAME[as.character(dp)],
             G = length(v), mean = mean(v), t = tt,
             p = 2 * pt(-abs(tt), df = length(v) - 1))
}))
print(loo[, .(dropped, G, mean = round(mean, 6), t = round(t, 2),
              p = round(p, 3))], row.names = FALSE)

save(tab, "t19_im_by_pref"); save(im, "t19_im_test"); save(loo, "t19_im_loo")
cat("\n出力: research/model/output/t19_{im_by_pref,im_test}.csv\n")

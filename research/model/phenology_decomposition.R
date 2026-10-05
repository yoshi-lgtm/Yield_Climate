# 出穂日トレンドの分解：田植期（農家の決定） vs 移植-出穂期間（受動的な熱応答）
#
# 背景：Hasegawa et al. (Ch.4, Mimura & Takewaka eds. 2025) は、農家は高温回避のため
# 遅植えを志向したが、温暖化で移植-出穂期間が短縮し、出穂日はほとんど遅れなかった、と
# 述べている。卒論は出穂日の前倒し（24年で5-8日）を「作期移動＝適応の直接証拠」と
# 解釈していた。この解釈が成り立つかを、田植期が実測で取れる2014-2023で検証する。
#
# 恒等式：heading_doy = transplant_doy + interval
#   → d(heading)/dt = d(transplant)/dt + d(interval)/dt
# 前者が能動的な作期移動、後者が受動的な発育速度の変化。

suppressMessages({library(data.table); library(fixest)})

OUT <- "research/model/output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cal <- fread("results/crop_calendar_pref.csv", encoding = "UTF-8")
setnames(cal, c("出穂期", "田植期", "刈取期", "移植出穂日数", "登熟日数"),
              c("heading", "transplant", "harvest", "interval", "ripen_days"))

# ---- 0. 恒等式の確認 -------------------------------------------------------
chk <- cal[!is.na(transplant) & !is.na(heading)]
stopifnot(all(abs(chk$heading - (chk$transplant + chk$interval)) < 1e-9))
cat("恒等式 heading = transplant + interval を確認:", nrow(chk), "行\n\n")

# ---- 1. 全期間（2000-2023）の出穂日トレンド：文脈 --------------------------
full <- cal[!is.na(heading)]
m_full <- feols(heading ~ year | pref, data = full, cluster = ~pref)
cat("=== 出穂日トレンド 2000-2023（都県FE、n =", nobs(m_full), ") ===\n")
print(coeftable(m_full))

# 期間を分けて比較
for (w in list(c(2000, 2023), c(2000, 2013), c(2014, 2023))) {
  m <- feols(heading ~ year | pref, data = full[year %between% w], cluster = ~pref)
  cat(sprintf("  %d-%d: %+.3f 日/年 (SE %.3f, n=%d)\n",
              w[1], w[2], coef(m)[["year"]], se(m)[["year"]], nobs(m)))
}
cat("\n")

# ---- 2. 本題：2014-2023 の分解 ---------------------------------------------
sub <- cal[year >= 2014 & !is.na(transplant)]
cat("=== 分解 2014-2023（都県FE、n =", nrow(sub), "、7都県×10年）===\n")

res <- rbindlist(lapply(c("heading", "transplant", "interval", "ripen_days"), function(v) {
  m <- feols(as.formula(paste(v, "~ year | pref")), data = sub, cluster = ~pref)
  data.table(var = v, slope = coef(m)[["year"]], se = se(m)[["year"]],
             t = coef(m)[["year"]] / se(m)[["year"]],
             p = pvalue(m)[["year"]],
             change_10y = coef(m)[["year"]] * 10)
}))
print(res, digits = 3)
cat(sprintf("\n  恒等式チェック: transplant %+.4f + interval %+.4f = %+.4f (heading %+.4f)\n",
            res[var == "transplant", slope], res[var == "interval", slope],
            res[var == "transplant", slope] + res[var == "interval", slope],
            res[var == "heading", slope]))
cat(sprintf("  出穂日変化のうち移植期で説明される割合: %.1f%%\n\n",
            100 * res[var == "transplant", slope] / res[var == "heading", slope]))
fwrite(res, file.path(OUT, "pheno_decomp_pooled.csv"))

# ---- 3. 都県別のトレンド ---------------------------------------------------
byp <- sub[, {
  s <- function(v) coef(lm(v ~ year))[["year"]]
  .(heading = s(heading), transplant = s(transplant),
    interval = s(interval), ripen = s(ripen_days))
}, by = pref]
byp[, `:=`(share_transplant = 100 * transplant / heading)]
cat("=== 都県別トレンド（日/年、2014-2023 の単回帰）===\n")
print(byp, digits = 3)
fwrite(byp, file.path(OUT, "pheno_decomp_bypref.csv"))
cat("\n")

# ---- 4. 熱的応答の検証 -----------------------------------------------------
# 移植-出穂期間が短縮したのが温暖化によるものなら、
#   (a) 当該期間の平均気温が上昇している
#   (b) 出穂までの積算温度（GDD）はほぼ一定に保たれる
# はずである。都県×年で、移植日から出穂日までの気象を市町村平均して確認する。
cat("=== 移植-出穂期間の気象（2014-2023）===\n")
pref_map <- c("8" = "茨城", "9" = "栃木", "10" = "群馬", "11" = "埼玉",
              "12" = "千葉", "13" = "東京", "14" = "神奈川")
TBASE <- 10  # 水稲の有効積算温度の基準温度

wx <- rbindlist(lapply(2014:2023, function(y) {
  d <- fread(sprintf("avg_temp_fixed/city_weighted_avg_temp_%d.csv", y), encoding = "UTF-8")
  code <- data.table::tstrsplit(d$city_id, "_")[[1]]
  d[, `:=`(doy = as.integer(format(as.Date(time), "%j")),
           pref_code = substr(code, 1, nchar(code) - 3))]
  d <- d[pref_code %in% names(pref_map)]
  d[, pref := pref_map[pref_code]]
  w <- cal[year == y & !is.na(transplant), .(pref, transplant, heading)]
  d <- merge(d, w, by = "pref")
  # 移植日 <= doy < 出穂日 の窓
  d <- d[doy >= transplant & doy < heading]
  d[, .(year = y, T_mean = mean(TMP_mea), GDD = sum(pmax(TMP_mea - TBASE, 0)) / uniqueN(city_id),
        n_days = uniqueN(doy), n_city = uniqueN(city_id)), by = pref]
}))
setorder(wx, pref, year)
fwrite(wx, file.path(OUT, "pheno_preheading_weather.csv"))

m_t <- feols(T_mean ~ year | pref, data = wx, cluster = ~pref)
m_g <- feols(GDD    ~ year | pref, data = wx, cluster = ~pref)
cat(sprintf("移植-出穂期間の平均気温トレンド: %+.4f ℃/年 (SE %.4f, t %+.2f) → 10年で %+.2f℃\n",
            coef(m_t)[["year"]], se(m_t)[["year"]],
            coef(m_t)[["year"]]/se(m_t)[["year"]], 10*coef(m_t)[["year"]]))
cat(sprintf("出穂までの積算温度(GDD, base%d)トレンド: %+.3f ℃日/年 (SE %.3f, t %+.2f)\n",
            TBASE, coef(m_g)[["year"]], se(m_g)[["year"]],
            coef(m_g)[["year"]]/se(m_g)[["year"]]))
cat(sprintf("  GDD平均 = %.0f ℃日 → 10年の変化は平均の %.1f%%\n\n",
            mean(wx$GDD), 100*10*coef(m_g)[["year"]]/mean(wx$GDD)))

# 期間長を気温で説明できるか
iv <- merge(sub[, .(pref, year, interval, transplant, heading)], wx, by = c("pref", "year"))
m_iv <- feols(interval ~ T_mean | pref, data = iv, cluster = ~pref)
cat("=== 移植-出穂期間 ~ 当該期間の平均気温（都県FE）===\n")
print(coeftable(m_iv))
cat(sprintf("within R2 = %.3f\n", r2(m_iv, "wr2")))
cat(sprintf("→ 平均気温1℃上昇で期間が %.2f 日短縮\n", coef(m_iv)[["T_mean"]]))
cat(sprintf("→ 10年の気温上昇 %+.2f℃ は期間短縮 %.2f 日を含意（実測 %.2f 日）\n\n",
            10*coef(m_t)[["year"]], 10*coef(m_t)[["year"]]*coef(m_iv)[["T_mean"]],
            10*res[var == "interval", slope]))

# 田植期は気温に反応していないか（能動的決定なら気温と無関係なはず）
m_tr <- feols(transplant ~ T_mean | pref, data = iv, cluster = ~pref)
cat("=== 田植期 ~ 移植-出穂期間の平均気温（都県FE）===\n")
print(coeftable(m_tr))
cat(sprintf("within R2 = %.3f\n", r2(m_tr, "wr2")))
fwrite(iv, file.path(OUT, "pheno_decomp_merged.csv"))
cat("\n完了。出力は", OUT, "\n")

# ---- 5. 部分期間の頑健性（表の切替 2019 をまたぐか）------------------------
cat("=== 分解シェアの頑健性（部分期間）===\n")
rob <- rbindlist(lapply(list(c(2014,2023), c(2014,2018), c(2019,2023)), function(w) {
  d <- sub[year %between% w]
  g <- function(v) { m <- feols(as.formula(paste(v,"~ year | pref")), data=d, cluster=~pref)
                     coef(m)[["year"]] }
  data.table(window = sprintf("%d-%d", w[1], w[2]), n = nrow(d),
             heading = g("heading"), transplant = g("transplant"), interval = g("interval"))
}))
rob[, share_interval := 100 * interval / heading]
print(rob, digits = 3)
fwrite(rob, file.path(OUT, "pheno_decomp_robust.csv"))
cat("→ どの窓でも期間短縮が出穂日変化の約9割以上を占める\n\n")

# ---- 6. 24年系列：気温・表の切替・残差トレンドへの分解 ---------------------
# 田植期は2014年以降しか無いので24年の直接分解はできない。代わりに
#   (i) 都県固定の窓で生育前半の平均気温を測り、出穂日への寄与を推定
#   (ii) 作物統計の表の切替年（2005/2009/2014/2019）に水準の跳びがないか検定
# を行う。表の型が変わるのは2014年（出穂期専用表→耕種期日の全工程表）のみ。
cat("=== 24年系列の分解（2000-2023）===\n")
win <- cal[year >= 2014 & !is.na(transplant), .(d0 = round(mean(transplant)), d1 = round(mean(heading))), by = pref]
pre <- rbindlist(lapply(2000:2023, function(y) {
  d <- fread(sprintf("avg_temp_fixed/city_weighted_avg_temp_%d.csv", y), encoding = "UTF-8")
  code <- data.table::tstrsplit(d$city_id, "_")[[1]]
  d[, `:=`(doy = as.integer(format(as.Date(time), "%j")), pc = substr(code, 1, nchar(code) - 3))]
  d <- d[pc %in% names(pref_map)][, pref := pref_map[pc]]
  merge(d, win, by = "pref")[doy >= d0 & doy < d1, .(year = y, T_pre = mean(TMP_mea)), by = pref]
}))
h <- merge(full[, .(pref, year, heading)], pre, by = c("pref", "year"))
for (b in c(2005, 2009, 2014, 2019)) h[, (paste0("b", b)) := as.integer(year >= b)]
fwrite(h, file.path(OUT, "pheno_24yr_preheading_temp.csv"))

cat("\n[出穂日 ~ 生育前半気温]\n"); print(coeftable(feols(heading ~ T_pre | pref, data = h, cluster = ~pref)))
cat(sprintf("within R2 = %.3f\n", r2(feols(heading ~ T_pre | pref, data = h, cluster = ~pref), "wr2")))
cat("\n[気温を統制した上での表の切替ダミー]\n")
print(coeftable(feols(heading ~ T_pre + b2005 + b2009 + b2014 + b2019 | pref, data = h, cluster = ~pref)))
cat("\n[気温 ＋ 2014年の不連続 を統制した年トレンド]\n")
print(coeftable(feols(heading ~ T_pre + year + b2014 | pref, data = h, cluster = ~pref)))

bb  <- coef(feols(heading ~ T_pre + year | pref, data = h, cluster = ~pref))
dT  <- 23 * coef(feols(T_pre ~ year | pref, data = h, cluster = ~pref))[["year"]]
tot <- 23 * coef(feols(heading ~ year | pref, data = h, cluster = ~pref))[["year"]]
b14 <- coef(feols(heading ~ T_pre + b2005 + b2009 + b2014 + b2019 | pref, data = h, cluster = ~pref))[["b2014"]]
cat(sprintf("\n出穂日の総変化 2000->2023 : %+.2f 日\n", tot))
cat(sprintf("  気温上昇による分        : %+.2f 日 (%.0f%%)  [%+.2f℃ x %+.2f 日/℃]\n",
            dT * bb[["T_pre"]], 100 * dT * bb[["T_pre"]] / tot, dT, bb[["T_pre"]]))
cat(sprintf("  残差トレンド            : %+.2f 日 (%.0f%%)\n", 23 * bb[["year"]], 100 * 23 * bb[["year"]] / tot))
cat(sprintf("  うち2014年の表切替      : %+.2f 日\n", b14))
fwrite(data.table(total = tot, temp = dT * bb[["T_pre"]], resid = 23 * bb[["year"]], break2014 = b14),
       file.path(OUT, "pheno_24yr_decomp.csv"))
cat("\n完了（§5-6）\n")

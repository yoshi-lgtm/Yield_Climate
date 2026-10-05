# T1: heat のサプライズ分解
#
# 問い：高温感応度 beta_t が2010年以降に縮小した（panel_extension §5.3）のは、
# 被害関数が真に非定常なのか、それとも定式化が誤っているのか。
#
# 仮説：効くのは heat の水準ではなく、予想からの乖離（サプライズ）である。
#   - 2010年 = 大きな正のサプライズ → 大きな被害
#   - 2023年 = 2010年並みに暑いがトレンド上 → 被害なし
#     （Hasegawa et al. 2025 第4章：2023年の関東の MET26 はトレンド上に留まった）
#   - 冷夏年（2003・2006）= 負のサプライズ → 係数が正
#
# 分解： heat_it = E_{t-1}[heat_it] + surprise_it
#   E は市町村ごとの後ろ向き移動平均（当年を含まない）。
#
# 推定：
#   log y = a_i + t_t + bE*E + bS*S + dE*(E x consol30) + dS*(S x consol30) + gX + e
#   heat = E + S なので、これは log y = a_i + t_t + bS*heat + (bE-bS)*E + ... の
#   再パラメータ化。bS は「期待水準を固定したときの当年 heat の効果」、
#   (bE - bS) は「期待水準が高い地域ほど被害が小さいか」＝適応の符号。
#
# 検証すべき予測：
#   P1. bS < 0 で、年別 bS_t は水準の beta_t より安定している
#   P2. bE - bS > 0（期待が高い＝暑さに慣れた地域ほど heat 1日あたりの被害が小さい）
#   P3. dS > 0 が2010年以外からも識別される（LOYO で消えない）

suppressMessages({library(data.table); library(fixest)})

# リポジトリルートからでも親の sotsuron/ からでも動くように解決する
ROOT <- if (file.exists("research/model/kanto_panel_2000_2023.csv")) "." else "Yield_Climate"
stopifnot(file.exists(file.path(ROOT, "research/model/kanto_panel_2000_2023.csv")))
OUT <- file.path(ROOT, "research/model/output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setFixest_ssc(ssc(K.fixef = "full"))   # estimate.R と自由度修正を揃える

d <- fread(file.path(ROOT, "research/model/kanto_panel_2000_2023.csv"), encoding = "UTF-8")
d <- d[!is.na(obs_HD35)]
setorder(d, city_id, year)
d[, heat := obs_HD35]

# ---- 1. 期待の構築（市町村別・後ろ向き・当年を除く）------------------------
# 合併で年が欠ける市町村があるので、暦年グリッドに合わせてから窓を取る
grid <- CJ(city_id = unique(d$city_id), year = 2000:2023)
g <- merge(grid, d[, .(city_id, year, heat)], by = c("city_id", "year"), all.x = TRUE)
setorder(g, city_id, year)

trail_mean <- function(x, k) {
  # 直前 k 年の平均（当年を含まない）。有効年が k の半分未満なら NA
  n <- length(x); out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (i == 1L) next
    w <- x[max(1L, i - k):(i - 1L)]
    w <- w[!is.na(w)]
    if (length(w) >= max(2L, ceiling(k / 2))) out[i] <- mean(w)
  }
  out
}
expand_mean <- function(x) {
  n <- length(x); out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (i == 1L) next
    w <- x[seq_len(i - 1L)]; w <- w[!is.na(w)]
    if (length(w) >= 3L) out[i] <- mean(w)
  }
  out
}
g[, `:=`(E_ma3  = trail_mean(heat, 3L),
         E_ma5  = trail_mean(heat, 5L),
         E_ma10 = trail_mean(heat, 10L),
         E_exp  = expand_mean(heat)), by = city_id]

d <- merge(d, g[, .(city_id, year, E_ma3, E_ma5, E_ma10, E_exp)], by = c("city_id", "year"))
for (k in c("ma3", "ma5", "ma10", "exp")) d[[paste0("S_", k)]] <- d$heat - d[[paste0("E_", k)]]

cat("=== 期待の定義別の標本 ===\n")
for (k in c("ma3", "ma5", "ma10", "exp")) {
  e <- d[[paste0("E_", k)]]
  cat(sprintf("E_%-4s : n = %4d  開始年 %d  mean(E) = %.3f  sd(E) = %.3f  sd(S) = %.3f\n",
              k, sum(!is.na(e)), min(d$year[!is.na(e)]), mean(e, na.rm = TRUE),
              sd(e, na.rm = TRUE), sd(d[[paste0("S_", k)]], na.rm = TRUE)))
}

# ---- 2. 分散分解（市町村FE＋年FE を除去した後）------------------------------
cat("\n=== heat の within 変動のうち、期待とサプライズの寄与（MA5）===\n")
sub <- d[!is.na(E_ma5)]
dm <- function(v) resid(feols(as.formula(paste(v, "~ 1 | city_id + year")), data = sub))
vh <- var(dm("heat")); ve <- var(dm("E_ma5")); vs <- var(dm("S_ma5"))
cat(sprintf("Var(heat~) = %.3f,  Var(E~) = %.3f (%.0f%%),  Var(S~) = %.3f (%.0f%%),  2Cov = %.3f\n",
            vh, ve, 100 * ve / vh, vs, 100 * vs / vh, vh - ve - vs))
cat(sprintf("corr(E~, S~) = %+.3f\n", cor(dm("E_ma5"), dm("S_ma5"))))

# ---- 3. 主推定 -------------------------------------------------------------
CTRL <- "GSR_78 + APCP_78"
fit <- function(rhs, data = sub) feols(as.formula(paste("log_yield ~", rhs, "| city_id + year")),
                                      data = data, cluster = ~city_id)
show <- function(m, lab) {
  cat(sprintf("\n--- %s  (n = %d, within R2 = %.4f)\n", lab, nobs(m), r2(m, "wr2")))
  print(round(coeftable(m)[, 1:4], 6))
}

cat("\n\n===== 3. 主推定（期待 = MA5、標本 2005-2023）=====\n")
m0 <- fit(paste("heat +", CTRL));                  show(m0, "(0) 水準のみ")
m1 <- fit(paste("E_ma5 + S_ma5 +", CTRL));         show(m1, "(1) 期待＋サプライズ")
m2 <- fit(paste("heat + E_ma5 +", CTRL));          show(m2, "(1b) 再パラメータ化 heat + E")
m3 <- fit(paste("heat + heat:consol30 +", CTRL));  show(m3, "(2) 水準×consol30")
m4 <- fit(paste("E_ma5 + S_ma5 + E_ma5:consol30 + S_ma5:consol30 +", CTRL))
show(m4, "(3) 期待・サプライズ × consol30")

bE <- coef(m1)[["E_ma5"]]; bS <- coef(m1)[["S_ma5"]]
cat(sprintf("\nP2 検証: bE - bS = %+.5f （正なら「期待が高い地域ほど heat 1日の被害が小さい」＝適応）\n", bE - bS))
cat(sprintf("          (1b) の E の係数 = %+.5f がちょうどこれに対応\n", coef(m2)[["E_ma5"]]))

res <- rbindlist(lapply(list(list("水準のみ", m0), list("期待＋サプライズ", m1),
                             list("水準×consol30", m3), list("期待・サプライズ×consol30", m4)),
  function(z) { ct <- coeftable(z[[2]])
    data.table(spec = z[[1]], term = rownames(ct), est = ct[, 1], se = ct[, 2],
               t = ct[, 3], p = ct[, 4], n = nobs(z[[2]]), wr2 = r2(z[[2]], "wr2")) }))
fwrite(res, file.path(OUT, "surprise_main.csv"))

# ---- 4. P1：年別の係数の安定性 ---------------------------------------------
cat("\n\n===== 4. 年別係数の安定性（P1）=====\n")
yr <- sort(unique(sub$year))
bl <- feols(log_yield ~ i(year, heat) + GSR_78 + APCP_78 | city_id + year, data = sub, cluster = ~city_id)
bs <- feols(log_yield ~ i(year, S_ma5) + E_ma5 + GSR_78 + APCP_78 | city_id + year, data = sub, cluster = ~city_id)
pick <- function(m, pat) {
  ct <- coeftable(m); k <- grep(pat, rownames(ct))
  data.table(year = as.integer(sub(".*?(\\d{4}).*", "\\1", rownames(ct)[k])), est = ct[k, 1], se = ct[k, 2])
}
byyear <- merge(pick(bl, "year::.*:heat"), pick(bs, "year::.*:S_ma5"),
                by = "year", suffixes = c("_lvl", "_sur"))
print(byyear, digits = 3)
cat(sprintf("\n水準       beta_t : sd = %.5f, 範囲 [%+.5f, %+.5f], 符号が正の年 = %d/%d\n",
            sd(byyear$est_lvl), min(byyear$est_lvl), max(byyear$est_lvl),
            sum(byyear$est_lvl > 0), nrow(byyear)))
cat(sprintf("サプライズ bS_t : sd = %.5f, 範囲 [%+.5f, %+.5f], 符号が正の年 = %d/%d\n",
            sd(byyear$est_sur), min(byyear$est_sur), max(byyear$est_sur),
            sum(byyear$est_sur > 0), nrow(byyear)))
cat(sprintf("→ 年次ばらつきの比: %.2f 倍（1未満ならサプライズの方が安定）\n",
            sd(byyear$est_sur) / sd(byyear$est_lvl)))
fwrite(byyear, file.path(OUT, "surprise_byyear.csv"))

cat("\n--- 5年ブロック平均 ---\n")
sub[, blk := cut(year, c(2004, 2009, 2014, 2019, 2023),
                 labels = c("2005-09", "2010-14", "2015-19", "2020-23"))]
blk <- rbindlist(lapply(levels(sub$blk), function(b) {
  dd <- sub[blk == b]
  a  <- feols(log_yield ~ heat + GSR_78 + APCP_78 | city_id + year, data = dd, cluster = ~city_id)
  c2 <- feols(log_yield ~ S_ma5 + E_ma5 + GSR_78 + APCP_78 | city_id + year, data = dd, cluster = ~city_id)
  data.table(block = b, n = nobs(a),
             beta_lvl = coef(a)[["heat"]], se_lvl = se(a)[["heat"]],
             bS = coef(c2)[["S_ma5"]], se_S = se(c2)[["S_ma5"]])
}))
print(blk, digits = 3)
fwrite(blk, file.path(OUT, "surprise_byblock.csv"))

# ---- 5. P3：LOYO（交差項が2010年以外からも識別されるか）--------------------
cat("\n\n===== 5. leave-one-year-out（P3）=====\n")
loyo <- rbindlist(lapply(yr, function(y) {
  dd <- sub[year != y]
  a <- feols(log_yield ~ heat + heat:consol30 + GSR_78 + APCP_78 | city_id + year,
             data = dd, cluster = ~city_id)
  b <- feols(log_yield ~ E_ma5 + S_ma5 + E_ma5:consol30 + S_ma5:consol30 + GSR_78 + APCP_78 | city_id + year,
             data = dd, cluster = ~city_id)
  data.table(drop = y,
             d_lvl = coef(a)[["heat:consol30"]],
             t_lvl = coef(a)[["heat:consol30"]] / se(a)[["heat:consol30"]],
             d_sur = coef(b)[["S_ma5:consol30"]],
             t_sur = coef(b)[["S_ma5:consol30"]] / se(b)[["S_ma5:consol30"]])
}))
print(loyo, digits = 3)
cat(sprintf("\n水準       delta: 範囲 [%+.5f, %+.5f]、|t| < 1 の年 = %s\n",
            min(loyo$d_lvl), max(loyo$d_lvl), paste(loyo$drop[abs(loyo$t_lvl) < 1], collapse = ", ")))
cat(sprintf("サプライズ delta: 範囲 [%+.5f, %+.5f]、|t| < 1 の年 = %s\n",
            min(loyo$d_sur), max(loyo$d_sur), paste(loyo$drop[abs(loyo$t_sur) < 1], collapse = ", ")))
fwrite(loyo, file.path(OUT, "surprise_loyo.csv"))

# ---- 6. 期待の定義に対する頑健性 -------------------------------------------
cat("\n\n===== 6. 期待の定義に対する頑健性 =====\n")
rob <- rbindlist(lapply(c("ma3", "ma5", "ma10", "exp"), function(k) {
  E <- paste0("E_", k); S <- paste0("S_", k); SC <- paste0(S, ":consol30")
  dd <- d[!is.na(get(E))]
  f1 <- sprintf("log_yield ~ %s + %s + %s | city_id + year", E, S, CTRL)
  f2 <- sprintf("log_yield ~ %s + %s + %s:consol30 + %s:consol30 + %s | city_id + year", E, S, E, S, CTRL)
  a  <- feols(as.formula(f1), data = dd, cluster = ~city_id)
  b  <- feols(as.formula(f2), data = dd, cluster = ~city_id)
  bb <- feols(as.formula(f2), data = dd[year != 2010], cluster = ~city_id)
  data.table(E_def = k, n = nobs(a), from = min(dd$year),
             bE = coef(a)[[E]], bS = coef(a)[[S]], t_S = coef(a)[[S]] / se(a)[[S]],
             dS = coef(b)[[SC]], t_dS = coef(b)[[SC]] / se(b)[[SC]],
             dS_no2010 = coef(bb)[[SC]], t_dS_no2010 = coef(bb)[[SC]] / se(bb)[[SC]])
}))
print(rob, digits = 3)
fwrite(rob, file.path(OUT, "surprise_robust.csv"))

cat("\n完了。出力は", OUT, "\n")

# ---- 7. 仮説の直接形：年ごとの平均サプライズと、その年の被害の大きさ --------
# 年FE を入れている以上、ある年の内部では heat と S は E（緩慢に動く）の分だけしか
# 違わない。したがって年別係数の比較では両者をほとんど区別できない（§4 がそう）。
# サプライズ仮説の本来の内容は「年ごとの平均的な被害の大きさが、その年の
# 平均サプライズで説明される」である。これを直接検証する。
cat("\n\n===== 7. 年平均サプライズ vs その年の被害（仮説の直接形）=====\n")
ybar <- sub[, .(heat_bar = mean(heat), S_bar = mean(S_ma5), E_bar = mean(E_ma5)), by = year]
yb <- merge(byyear[, .(year, beta_lvl = est_lvl, bS = est_sur)], ybar, by = "year")
setorder(yb, year)
print(yb, digits = 3)

cat("\n[相関]\n")
cat(sprintf("corr(beta_t, heat_bar) = %+.3f\n", cor(yb$beta_lvl, yb$heat_bar)))
cat(sprintf("corr(beta_t, S_bar)    = %+.3f   <- 仮説は負（サプライズが大きい年ほど強い被害）\n",
            cor(yb$beta_lvl, yb$S_bar)))
print(coeftable(feols(beta_lvl ~ S_bar, data = yb)))
cat("\n[2010年と2023年の比較]\n")
print(yb[year %in% c(2003, 2006, 2010, 2018, 2020, 2023)], digits = 3)
fwrite(yb, file.path(OUT, "surprise_yearlevel.csv"))

# ---- 8. 非対称性：正のサプライズと負のサプライズ ----------------------------
# 「暑い方への外れ」だけが効き、「涼しい方への外れ」は効かない可能性。
cat("\n\n===== 8. サプライズの非対称性 =====\n")
sub[, `:=`(S_pos = pmax(S_ma5, 0), S_neg = pmin(S_ma5, 0))]
m5 <- fit(paste("E_ma5 + S_pos + S_neg +", CTRL));  show(m5, "(4) 正／負のサプライズを分離")
m6 <- fit(paste("E_ma5 + S_pos + S_neg + S_pos:consol30 + S_neg:consol30 +", CTRL))
show(m6, "(5) ＋ consol30 との交差")
cat(sprintf("\n正のサプライズ1日 = %+.5f、負のサプライズ1日 = %+.5f\n",
            coef(m5)[["S_pos"]], coef(m5)[["S_neg"]]))
cat("（負の側の係数が負なら「涼しい外れで収量が上がる」＝冷害の不在。正の側だけが\n")
cat("  高温害。両者の絶対値が近ければ、単一の線形 heat で近似できていることになる）\n")
w <- wald(m5, "S_pos|S_neg")
cat(sprintf("S_pos = S_neg の検定: "))
print(coeftable(feols(log_yield ~ E_ma5 + S_ma5 + I(S_pos) + GSR_78 + APCP_78 | city_id + year,
                      data = sub, cluster = ~city_id))["I(S_pos)", , drop = FALSE])
cat("（I(S_pos) の係数が正負の係数差。有意なら非対称）\n")

rob8 <- rbindlist(lapply(list(list("正負分離", m5), list("正負分離×consol30", m6)), function(z) {
  ct <- coeftable(z[[2]])
  data.table(spec = z[[1]], term = rownames(ct), est = ct[, 1], se = ct[, 2], t = ct[, 3], p = ct[, 4]) }))
fwrite(rob8, file.path(OUT, "surprise_asym.csv"))
cat("\n完了（§7-8）\n")

# ---- 9. 決定的な検定：非対称スペックの delta は2010年を抜いても残るか -------
# §8 で S_pos x consol30 が5%有意（t = 2.07）になった。T1 の本来の目的は
# 「delta が2010年以外からも識別されるか」なので、ここを LOYO で確かめる。
cat("\n\n===== 9. 非対称スペックの LOYO（本命の検定）=====\n")
F_ASYM <- paste("log_yield ~ E_ma5 + S_pos + S_neg + S_pos:consol30 + S_neg:consol30 +",
                CTRL, "| city_id + year")
loyo9 <- rbindlist(lapply(yr, function(y) {
  dd <- sub[year != y]
  m <- feols(as.formula(F_ASYM), data = dd, cluster = ~city_id)
  data.table(drop = y, n = nobs(m),
             S_pos = coef(m)[["S_pos"]],
             d_pos = coef(m)[["S_pos:consol30"]],
             t_pos = coef(m)[["S_pos:consol30"]] / se(m)[["S_pos:consol30"]],
             d_neg = coef(m)[["S_neg:consol30"]],
             t_neg = coef(m)[["S_neg:consol30"]] / se(m)[["S_neg:consol30"]])
}))
print(loyo9, digits = 3)
cat(sprintf("\nS_pos x consol30: 範囲 [%+.5f, %+.5f]、t の最小 = %+.2f（除外年 %d）\n",
            min(loyo9$d_pos), max(loyo9$d_pos),
            min(loyo9$t_pos), loyo9$drop[which.min(loyo9$t_pos)]))
cat(sprintf("  t > 1.64 の年 = %d/%d、t > 1.96 の年 = %d/%d\n",
            sum(loyo9$t_pos > 1.64), nrow(loyo9), sum(loyo9$t_pos > 1.96), nrow(loyo9)))
cat(sprintf("  2010年を除外したとき: %+.5f (t = %+.2f)\n",
            loyo9[drop == 2010, d_pos], loyo9[drop == 2010, t_pos]))
fwrite(loyo9, file.path(OUT, "surprise_asym_loyo.csv"))

# 2010年を抜いた本体も表示
cat("\n--- 2010年を除いた非対称スペック ---\n")
show(feols(as.formula(F_ASYM), data = sub[year != 2010], cluster = ~city_id), "(5) 2010年除外")
cat("\n--- 2010年・2023年を除いた非対称スペック ---\n")
show(feols(as.formula(F_ASYM), data = sub[!year %in% c(2010, 2023)], cluster = ~city_id),
     "(5) 2010・2023年除外")

# 期待の定義に対する頑健性（非対称スペック）
cat("\n--- 期待の定義に対する頑健性（非対称スペック、S_pos x consol30）---\n")
rob9 <- rbindlist(lapply(c("ma3", "ma5", "ma10", "exp"), function(k) {
  S <- paste0("S_", k); E <- paste0("E_", k)
  dd <- d[!is.na(get(E))]
  dd[, `:=`(sp = pmax(get(S), 0), sn = pmin(get(S), 0))]
  f <- sprintf("log_yield ~ %s + sp + sn + sp:consol30 + sn:consol30 + %s | city_id + year", E, CTRL)
  a <- feols(as.formula(f), data = dd, cluster = ~city_id)
  b <- feols(as.formula(f), data = dd[year != 2010], cluster = ~city_id)
  data.table(E_def = k, n = nobs(a),
             S_pos = coef(a)[["sp"]], S_neg = coef(a)[["sn"]],
             d_pos = coef(a)[["sp:consol30"]], t_pos = coef(a)[["sp:consol30"]] / se(a)[["sp:consol30"]],
             d_pos_no2010 = coef(b)[["sp:consol30"]],
             t_pos_no2010 = coef(b)[["sp:consol30"]] / se(b)[["sp:consol30"]])
}))
print(rob9, digits = 3)
fwrite(rob9, file.path(OUT, "surprise_asym_robust.csv"))

# 水準 heat でも正負に折ると delta は出るのか（サプライズ分解が本質か、
# 非対称性が本質かの切り分け）
cat("\n--- 切り分け：期待を使わず、heat を市町村平均で折った場合 ---\n")
sub[, hbar := mean(heat), by = city_id]
sub[, `:=`(H_pos = pmax(heat - hbar, 0), H_neg = pmin(heat - hbar, 0))]
f_h <- paste("H_pos + H_neg + H_pos:consol30 + H_neg:consol30 +", CTRL)
show(fit(f_h), "(6) 市町村平均からの偏差を正負に折る")
show(fit(f_h, sub[year != 2010]), "(6) 同上・2010年除外")
cat("\n完了（§9）\n")

# ---- 10. bE は期間について安定か -------------------------------------------
# §6 は bE が「スペック間・2010年除外に対して頑健」と述べたが、期間に対しては
# 確かめていなかった。修論の E[y] に bE を使えるかを決めるので、ここで検証する。
cat("\n\n===== 10. bE / bS の期間別安定性 =====\n")
stab <- rbindlist(lapply(c(list("全期間"), as.list(levels(sub$blk))), function(b) {
  dd <- if (identical(b, "全期間")) sub else sub[blk == b]
  m <- feols(log_yield ~ E_ma5 + S_ma5 + GSR_78 + APCP_78 | city_id + year,
             data = dd, cluster = ~city_id)
  data.table(period = b, n = nobs(m),
             bE = coef(m)[["E_ma5"]], se_E = se(m)[["E_ma5"]],
             bS = coef(m)[["S_ma5"]], se_S = se(m)[["S_ma5"]])
}))
stab[, `:=`(t_E = bE / se_E, t_S = bS / se_S)]
print(stab, digits = 3)
cat(sprintf("\n期間間の sd: bE = %.5f,  bS = %.5f  → %s の方が不安定\n",
            sd(stab[period != "全期間", bE]), sd(stab[period != "全期間", bS]),
            ifelse(sd(stab[period != "全期間", bE]) > sd(stab[period != "全期間", bS]), "bE", "bS")))
cat("\n--- 2010年除外（全期間）---\n")
m10 <- feols(log_yield ~ E_ma5 + S_ma5 + GSR_78 + APCP_78 | city_id + year,
             data = sub[year != 2010], cluster = ~city_id)
print(round(coeftable(m10)[c("E_ma5", "S_ma5"), 1:3], 5))
fwrite(stab, file.path(OUT, "surprise_bE_stability.csv"))
cat("\n完了（§10）\n")

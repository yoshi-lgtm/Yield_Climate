# =============================================================================
# T4: 2010年イベント（年別係数）と交互作用の horse race
#
# 実行はリポジトリルート(Yield_Climate/)またはその親から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t4_event_2010.R
#
# 入力 : research/model/kanto_panel_2000_2023.csv
#        results/city_centroids.csv
# 出力 : research/model/output/t4_*.csv と標準出力
#
# goals_and_todo.md G1/T4。本番スペックは T15 の飽和版（年 × consol30 入り）。
#
# T4 の本体は horse race である（頑健性ではない）。delta の識別仮定は
# 「heat 曝露を条件づけたとき、consol30 と相関する他の市町村特性が高温被害を
# 修飾しない」。単一事象設計では「複数事象で頑健」という逃げ道がないので、
# これを本体に置く。
#
# T15 で delta_2010 が t 2.07 -> 0.98 になった原因は 2010年内で heat*c と c が
# 共線なことだった。したがって次の問いは「consol30 固有の何かが効いているのか、
# 規模や地理の代理になっているだけか」である。
#
# 機序の切り分け：理論は「整備 -> 用排水インフラ -> 掛け流し・深水かんがい ->
# 水温低下 -> 高温害緩和」だが consol30 が測っているのは区画の大きさ。
# 水管理の話なら pipe の方が直接的。consol30 で効いて pipe で効かないなら
# 水温の物語は弱く、「大区画＝規模の大きい農家」という別解釈が残る。
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
EVT     <- 2010L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
save <- function(dt, name) fwrite(dt, file.path(OUT_DIR, paste0(name, ".csv")), bom = TRUE)

d <- fread(PANEL, encoding = "UTF-8")
cen <- fread(CENTR, encoding = "UTF-8")[, .(city_id, lon, lat)]
d <- merge(d, cen, by = "city_id", all.x = TRUE)

# 標高（Data_cleaning/13build_city_elevation.py）。T4 第2.4節で経度が
# t = 4.1〜4.4 で効いたので、その代理の第一候補として入れる。
ELEV <- if (file.exists("results/city_elevation.csv")) {
  "results/city_elevation.csv"
} else { "Yield_Climate/results/city_elevation.csv" }
ev <- fread(ELEV, encoding = "UTF-8")[, .(city_id, elev, slope, relief)]
d <- merge(d, ev, by = "city_id", all.x = TRUE)
d[, log_elev := log1p(elev)]

# 海岸距離（Data_cleaning/14build_city_coastdist.py）。T4 第2.5節で標高が
# 否定されたので、経度の代理候補の残り第一候補。経度との相関は -0.594。
CDIST <- if (file.exists("results/city_coastdist.csv")) {
  "results/city_coastdist.csv"
} else { "Yield_Climate/results/city_coastdist.csv" }
cd <- fread(CDIST, encoding = "UTF-8")[, .(city_id, coast = coast_km)]
d <- merge(d, cd, by = "city_id", all.x = TRUE)
d[, log_coast := log1p(coast)]

# 森林率（Data_cleaning/15build_city_forest.py）。標高・傾斜とは別の機序
# （日中の遮蔽・冷気流出）の候補。経度との相関は -0.248/-0.144 と弱い。
FRST <- if (file.exists("results/city_forest.csv")) {
  "results/city_forest.csv"
} else { "Yield_Climate/results/city_forest.csv" }
fr <- fread(FRST, encoding = "UTF-8")[, .(city_id, forest_city, forest_nbr)]
d <- merge(d, fr, by = "city_id", all.x = TRUE)

vars <- c("log_yield", HEAT, "consol30", "pipe", "paddy_fix", "den",
          "elev", "log_elev", "slope", "coast", "log_coast",
          "forest_city", "forest_nbr", CTRL, "city_id", "year", "lon", "lat")
s <- d[complete.cases(d[, ..vars])][den > 0]
s[, pref := city_id %/% 1000L]
years <- sort(unique(s$year)); REF <- years[1]
cat(sprintf("n = %d, 市町村 = %d, 都県 = %d, 年 = %d\n",
            nrow(s), uniqueN(s$city_id), uniqueN(s$pref), length(years)))

# -----------------------------------------------------------------------------
# 0. 競合変数。すべて時間不変。比較できるよう市町村断面で標準化する。
#    consol30 だけは既存文書と数値を揃えるため生のまま使う。
# -----------------------------------------------------------------------------
s[, log_den := log(den)]
COMP <- c("pipe", "paddy_fix", "log_den", "lat", "lon", "elev", "log_elev",
          "slope", "coast", "log_coast", "forest_city", "forest_nbr")
cs <- unique(s[, c("city_id", COMP), with = FALSE])
for (v in COMP) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = paste0("z_", v), value = (s[[v]] - mu) / sg)
}
Z <- paste0("z_", COMP)
cat("競合（市町村断面で標準化）:", paste(COMP, collapse = ", "), "\n")
cat("consol30 との相関:\n")
print(round(cor(unique(s[, c("consol30", Z), with = FALSE])[, -1],
                unique(s[, c("consol30", Z), with = FALSE])$consol30), 3))

# 交差項と年×水準項
s[, h_consol := get(HEAT) * consol30]
for (v in Z) set(s, j = paste0("h_", v), value = s[[HEAT]] * s[[v]])
mk_levels <- function(v, tag) {
  nms <- character(0)
  for (y in years[-1]) {
    nm <- sprintf("L%s_%d", tag, y)
    set(s, j = nm, value = s[[v]] * as.integer(s$year == y))
    nms <- c(nms, nm)
  }
  nms
}
Lc <- mk_levels("consol30", "c")
Lz <- lapply(Z, function(v) mk_levels(v, sub("^z_", "", v)))
names(Lz) <- Z

# -----------------------------------------------------------------------------
# 1. 共通の推定・集計
# -----------------------------------------------------------------------------
pick_idx <- function(mm, v) {
  parts <- strsplit(names(coef(mm)), ":", fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) >= 4 && z[1] == "year" && z[4] == v, logical(1))
  idx <- which(ok)
  idx[order(as.integer(vapply(parts[idx], function(z) z[3], character(1))))]
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

# 年飽和スペック。extra = 追加する競合（交差項を年別に入れる）
run_sat <- function(label, extra_z) {
  rhs <- c(sprintf("i(year, %s)", HEAT), "i(year, h_consol)", Lc)
  for (v in extra_z) rhs <- c(rhs, sprintf("i(year, h_%s)", v), Lz[[v]])
  rhs <- c(rhs, CTRL)
  f <- as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rhs, collapse = " + ")))
  m <- feols(f, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
  id <- pick_idx(m, "h_consol"); V <- vcov(m)
  b <- coef(m)[id]; se <- sqrt(pmax(diag(V[id, id]), 0))
  i10 <- which(years == EVT)
  w <- wald_equal(b, V[id, id])
  # 競合側の 2010 係数
  comp10 <- vapply(extra_z, function(v) {
    j <- pick_idx(m, paste0("h_", v)); unname(coef(m)[j][i10])
  }, numeric(1))
  comp10se <- vapply(extra_z, function(v) {
    j <- pick_idx(m, paste0("h_", v)); sqrt(diag(V)[j][i10])
  }, numeric(1))
  cat(sprintf("\n-- %-28s delta_2010 %+.6f (SE %.6f, t %+.2f)  chi2(23) %7.1f\n",
              label, b[i10], se[i10], b[i10] / se[i10], w[1]))
  if (length(extra_z)) {
    for (k in seq_along(extra_z)) {
      cat(sprintf("     %-14s の 2010年交差項 %+.6f (SE %.6f, t %+.2f)\n",
                  extra_z[k], comp10[k], comp10se[k], comp10[k] / comp10se[k]))
    }
  }
  list(row = data.table(spec = label, n_comp = length(extra_z),
         delta_2010 = b[i10], se_2010 = se[i10], t_2010 = b[i10] / se[i10],
         wald = w[1], wald_p = w[2],
         comp = paste(extra_z, collapse = "+"),
         comp_2010 = if (length(extra_z) == 1) comp10[1] else NA_real_,
         comp_t_2010 = if (length(extra_z) == 1) comp10[1] / comp10se[1] else NA_real_),
       byyear = data.table(spec = label, year = years, delta = b, se = se))
}

# -----------------------------------------------------------------------------
# 2. horse race
# -----------------------------------------------------------------------------
cat("\n==================== horse race（年飽和・本番スペック） ====================\n")
cat("参照は T15 の飽和版（競合なし）。競合は1つずつ入れ、最後に全部同時。\n")

out <- list(); byy <- list()
r <- run_sat("base (consol30 のみ)", character(0)); out[[1]] <- r$row; byy[[1]] <- r$byyear
for (k in seq_along(Z)) {
  r <- run_sat(sprintf("+ %s", COMP[k]), Z[k])
  out[[length(out) + 1]] <- r$row; byy[[length(byy) + 1]] <- r$byyear
}

# 地理を年飽和したまま同時に入れる。lat/lon は単独で (c) の統計量を大きく
# 下げるので、この2本の同時投入が (c) に対する最も厳しい検定になる。
# 標高は経度と -0.497 で相関する。「経度は標高の代理か」を直接問うため、
# elev と lon を同時に入れてどちらが2010年の交差項を担うかを見る。
GRP <- list(c("z_lat", "z_lon"), c("z_elev", "z_lon"),
            c("z_elev", "z_lat", "z_lon"),
            c("z_coast", "z_lon"), c("z_log_coast", "z_lon"),
            c("z_coast", "z_elev"), c("z_coast", "z_lat", "z_lon"),
            c("z_coast", "z_elev", "z_lat", "z_lon"),
            c("z_forest_nbr", "z_coast"),
            c("z_forest_nbr", "z_elev"),
            c("z_coast", "z_elev", "z_forest_nbr"),
            c("z_forest_nbr", "z_lon"), Z)
for (g in GRP) {
  r <- run_sat(sprintf("+ %s（年飽和）",
                       paste(sub("^z_", "", g), collapse = "+")), g)
  out[[length(out) + 1]] <- r$row; byy[[length(byy) + 1]] <- r$byyear
}

# 全競合同時（年飽和だと水準項が 5*23 本になり薄くなるので、競合は
# プール傾き + 2010年固有の傾き + 2010年の水準 に絞る）
cat("\n-- 全競合同時（競合はプール傾き＋2010年固有） --\n")
s[, I10 := as.integer(year == EVT)]
for (v in Z) {
  set(s, j = paste0("h10_", v), value = s[[paste0("h_", v)]] * s$I10)
  set(s, j = paste0("L10_", v), value = s[[v]] * s$I10)
}
rhs <- c(sprintf("i(year, %s)", HEAT), "i(year, h_consol)", Lc,
         paste0("h_", Z), paste0("h10_", Z), paste0("L10_", Z), CTRL)
f_all <- as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rhs, collapse = " + ")))
m_all <- feols(f_all, s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
id <- pick_idx(m_all, "h_consol"); V <- vcov(m_all)
b <- coef(m_all)[id]; se <- sqrt(pmax(diag(V[id, id]), 0)); i10 <- which(years == EVT)
w <- wald_equal(b, V[id, id])
cat(sprintf("-- %-28s delta_2010 %+.6f (SE %.6f, t %+.2f)  chi2(23) %7.1f\n",
            "+ 全競合", b[i10], se[i10], b[i10] / se[i10], w[1]))
for (v in Z) {
  cf <- coef(m_all); sv <- sqrt(diag(V))
  cat(sprintf("     %-14s プール %+.6f (t %+.2f) / 2010年追加 %+.6f (t %+.2f)\n",
              sub("^z_", "", v), cf[paste0("h_", v)],
              cf[paste0("h_", v)] / sv[paste0("h_", v)],
              cf[paste0("h10_", v)], cf[paste0("h10_", v)] / sv[paste0("h10_", v)]))
}
out[[length(out) + 1]] <- data.table(spec = "+ 全競合", n_comp = length(Z),
  delta_2010 = b[i10], se_2010 = se[i10], t_2010 = b[i10] / se[i10],
  wald = w[1], wald_p = w[2], comp = "all", comp_2010 = NA_real_, comp_t_2010 = NA_real_)
byy[[length(byy) + 1]] <- data.table(spec = "+ 全競合", year = years, delta = b, se = se)

# -----------------------------------------------------------------------------
# 2b. (c) に対する最も厳しい検定：heat の傾きを都県 × 年で自由にする
#     lat/lon を年飽和で入れると chi2 が 81.9 -> 29.9 に落ちたので、
#     「delta_t の異質性は空間 × 年の熱感応度勾配ではないか」を直接問う。
#     T12 第3.4節の「都県×年FE 自体が対処でもある」の傾き版にあたる。
# -----------------------------------------------------------------------------
cat("\n========== (c) の最も厳しい検定：heat の傾きを都県×年で自由にする ==========\n")
s[, pref_year := paste0(pref, "_", year)]
for (lab in c("heat 傾きを都県×年", "heat 傾きを都県×年 ＋ 都県×年FE")) {
  fe <- if (grepl("FE", lab)) "city_id + year + pref_year" else "city_id + year"
  f <- as.formula(sprintf(
    "log_yield ~ i(pref_year, %s) + i(year, h_consol) + %s + %s | %s",
    HEAT, paste(Lc, collapse = " + "), paste(CTRL, collapse = " + "), fe))
  mm <- tryCatch(feols(f, s, cluster = ~city_id, fixef.rm = "none",
                       ssc = ssc(K.fixef = "full")), error = function(e) NULL)
  if (is.null(mm)) { cat(sprintf("-- %s : 推定できず\n", lab)); next }
  idd <- pick_idx(mm, "h_consol"); Vv <- vcov(mm)
  bb <- coef(mm)[idd]; ss <- sqrt(pmax(diag(Vv[idd, idd]), 0))
  ww <- wald_equal(bb, Vv[idd, idd]); j10 <- which(years == EVT)
  cat(sprintf("-- %-32s delta_2010 %+.6f (t %+.2f)  chi2(%d) %7.1f  p %.3g\n",
              lab, bb[j10], bb[j10] / ss[j10], length(bb) - 1, ww[1], ww[2]))
  out[[length(out) + 1]] <- data.table(spec = lab, n_comp = NA_integer_,
    delta_2010 = bb[j10], se_2010 = ss[j10], t_2010 = bb[j10] / ss[j10],
    wald = ww[1], wald_p = ww[2], comp = "pref_year_heat",
    comp_2010 = NA_real_, comp_t_2010 = NA_real_)
}
cat(sprintf("参考: chi2(23) の5%%臨界値 %.1f / 1%%臨界値 %.1f\n",
            qchisq(0.95, 23), qchisq(0.99, 23)))

# -----------------------------------------------------------------------------
# 3. 機序の切り分け：consol30（区画規模）対 pipe（水管理）
# -----------------------------------------------------------------------------
cat("\n==================== 機序：consol30 対 pipe ====================\n")
cat("理論は水管理（掛け流し・深水かんがい）経由。consol30 が測るのは区画の大きさ。\n")
cat("consol30 で効いて pipe で効かないなら水温の物語は弱い。\n")
mp <- feols(as.formula(sprintf(
  "log_yield ~ i(year, %s) + i(year, h_consol) + i(year, h_z_pipe) + %s + %s + %s | city_id + year",
  HEAT, paste(Lc, collapse = " + "), paste(Lz[["z_pipe"]], collapse = " + "),
  paste(CTRL, collapse = " + "))),
  s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"))
ic <- pick_idx(mp, "h_consol"); ip <- pick_idx(mp, "h_z_pipe")
Vp <- vcov(mp); cfp <- coef(mp); svp <- sqrt(diag(Vp))
mech <- data.table(year = years,
  d_consol = cfp[ic], se_consol = svp[ic],
  d_pipe = cfp[ip], se_pipe = svp[ip])
mech[, `:=`(t_consol = d_consol / se_consol, t_pipe = d_pipe / se_pipe)]
print(mech[year %in% c(2007, 2010, 2012, 2013, 2015, 2018, 2023),
  .(year, d_consol = round(d_consol, 5), t_consol = round(t_consol, 2),
    d_pipe = round(d_pipe, 5), t_pipe = round(t_pipe, 2))], row.names = FALSE)
cat(sprintf("\n2010年: consol30 %+.6f (t %+.2f) / pipe %+.6f (t %+.2f)\n",
            mech[year == EVT, d_consol], mech[year == EVT, t_consol],
            mech[year == EVT, d_pipe], mech[year == EVT, t_pipe]))
cat(sprintf("|t|>1.96 の年数: consol30 %d / pipe %d（24年中）\n",
            mech[abs(t_consol) > 1.96, .N], mech[abs(t_pipe) > 1.96, .N]))

# -----------------------------------------------------------------------------
# 4. 副：long difference（base 年が任意なので副に留める）
# -----------------------------------------------------------------------------
cat("\n==================== 副：long difference ====================\n")
ld <- dcast(s[year %in% c(2007, 2008, 2009, EVT)], city_id + consol30 ~ year,
            value.var = c("log_yield", HEAT))
setnames(ld, gsub("-", "_", names(ld)))
ld <- ld[complete.cases(ld)]
ld[, dy := log_yield_2010 - rowMeans(cbind(log_yield_2007, log_yield_2008, log_yield_2009))]
ld[, h10 := obs_HD35_2010]
ld[, hc := h10 * consol30]
mld <- feols(dy ~ h10 + hc, ld, vcov = "hetero")
cat(sprintf("n = %d 市町村。base = 2007-09 平均\n", nrow(ld)))
print(mld)

# -----------------------------------------------------------------------------
# 5. 書き出し
# -----------------------------------------------------------------------------
hr <- rbindlist(out)
save(hr, "t4_horse_race"); save(rbindlist(byy), "t4_by_year"); save(mech, "t4_mechanism")
cat("\n==================== まとめ ====================\n")
print(hr[, .(spec, delta_2010 = round(delta_2010, 6), t_2010 = round(t_2010, 2),
             chi2_23 = round(wald, 1), p = signif(wald_p, 3),
             comp_2010 = round(comp_2010, 6), comp_t = round(comp_t_2010, 2))],
      row.names = FALSE)
cat("\n出力: research/model/output/t4_{horse_race,by_year,mechanism}.csv\n")

# =============================================================================
# T21: target estimand を明示し、プール推定量との距離を測る
#
# 実行はリポジトリルート(Yield_Climate/)から:
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" research/model/t21_target_estimand.R
#
# 出力 : research/model/output/t21_target_estimand.csv
#         research/model/output/t21_target_weights.csv
#         および標準出力
#
# 動機：T11・T18・T20 は「プール値の中身」を報告してきた。すなわち汚染 28.2%、
# 2010年の寄与 87.2% のように、**プール値そのものを分母にした構成比**である。
# これは「報告された数値は何で出来ているか」に答えるが、
# 「**本当に知りたかった量からどれだけ離れているか**」には答えない。
# 分母が内生的に選ばれているので、離れ方の大きさも符号も読めない。
# そこで target estimand theta を先に書き、プール値との差を測る。
#
# 恒等式（T18 §2）の副産物として、これは追加の推定を要さない：
#   delta_pool = sum_t w_D,t delta_t + sum_t w_h,t beta_t
#                                    + sum_t w_G,t gamma_t + sum_t w_F,t eta_t
# の重み w は X に条件付ければ回帰子行列だけの関数で**推定誤差を持たない**
# （T18 §2.4、outline §3.2）。target が theta = sum_t a_t delta_t の形であれば
# theta も、ギャップ delta_pool - theta も、**同じ飽和式の係数の固定ウェイト
# 線形結合**になる。したがって点推定・SE・CI はすべて同一の VCV から出る。
# delta method も bootstrap も要らない。
#
# ギャップの分解：
#   delta_pool - theta = sum_t (w_D,t - a_t) delta_t   … 重みのずれ
#                      + sum_t w_h,t beta_t + ...      … 汚染（target に依存しない）
# 汚染項は a に依存しないので、どの target でも**同じ値**が出る。
# 変わるのは重みのずれだけで、これが「target を変えると何が変わるか」の答えになる。
#
# target の選び方（post-selection を避ける）：
#   統制の選択（T18 §1）と同じ規律を target にも課す。
#   **検定統計量を見てから target を選ばない。** 下の4つは
#   (A) 推定量にもデータにも依存しない（1/24）
#   (B)(C) **曝露分布だけ**から決まる（左辺を一切見ない）
#   (D) 事前に決めてある設計年（T4・T17・T20 の 2010年）
#   のいずれかで、どれも delta_t や t 値を見て決めていない。
#   (E) は (B)(C) の中間を見るための参考で、主結果には置かない。
#
#   target のラベルは (A)-(E) とする。T## は本プロジェクトの todo 番号
#   （T1 = サプライズ分解、T4 = 2010年イベント …）に予約されているので、
#   target の番号付けには使わない。
#
# IM が使えない理由（本文に書く）：target はいずれも delta_t の線形結合なので、
# IM に乗せるには都県別の年飽和が必要になる。それは神奈川で41係数が落ちて
# 推定不能（T19 §1.1）、年を4期に粗視化すれば通るが粗視化そのものが
# 汚染項を壊す（28.2% -> 5.0%、T20 §C）。**(c) を塞ぐのと同じ G=6 の階数制約**
# である。よって本表はクラスターと Conley 50km の二本で出す。
# =============================================================================

suppressMessages({ library(fixest); library(data.table) })

pk <- function(x) x[file.exists(x)][1]
PANEL <- pk(c("research/model/kanto_panel_2000_2023.csv",
              "Yield_Climate/research/model/kanto_panel_2000_2023.csv"))
CEN <- pk(c("results/city_centroids.csv", "Yield_Climate/results/city_centroids.csv"))
CD  <- pk(c("results/city_coastdist.csv", "Yield_Climate/results/city_coastdist.csv"))
FR  <- pk(c("results/city_forest.csv",    "Yield_Climate/results/city_forest.csv"))
EL  <- pk(c("results/city_elevation.csv", "Yield_Climate/results/city_elevation.csv"))
stopifnot(!is.na(PANEL))
OUT <- file.path(dirname(PANEL), "output")
stopifnot(dir.exists(OUT))

CTRL   <- c("GSR_78", "APCP_78")   # 窓ベースの日射・降水（T20 と同一）
HEAT   <- "obs_HD35"               # 観測出穂アンカーの登熟期 >=35.0C 日数
CONLEY <- 50                       # km。lat/lon は**生の度**を渡す（outline §12）
NTOP   <- 5                        # (C) で平均する高曝露年の本数

sv <- function(x, nm) {
  p <- file.path(OUT, paste0(nm, ".csv"))
  fwrite(x, p, bom = TRUE)
  cat(sprintf("  -> %s\n", p))
}

# -----------------------------------------------------------------------------
# 0. 標本（T20 §0 と同一の作り方。ここを変えると T18・T20 と比較できない）
# -----------------------------------------------------------------------------
d <- fread(PANEL, encoding = "UTF-8")
d <- merge(d, fread(CEN, encoding = "UTF-8")[, .(city_id, lon, lat)], by = "city_id")
d <- merge(d, fread(CD,  encoding = "UTF-8")[, .(city_id, coast = coast_km)], by = "city_id")
d <- merge(d, fread(FR,  encoding = "UTF-8")[, .(city_id, forest = forest_nbr)], by = "city_id")
d <- merge(d, fread(EL,  encoding = "UTF-8")[, .(city_id, elev)], by = "city_id")

vv <- c("log_yield", HEAT, "consol30", "coast", "forest", "elev", "lon", "lat",
        CTRL, "city_id", "year")
s <- d[complete.cases(d[, ..vv])]

# Conley 用に**生の度**を退避してから標準化する（T20 のバグ修正、outline §12）
s[, `:=`(lat_deg = lat, lon_deg = lon)]
cs <- unique(s[, c("city_id", c("coast", "forest", "elev", "lon", "lat")), with = FALSE])
for (v in c("coast", "forest", "elev", "lon", "lat")) {
  mu <- mean(cs[[v]]); sg <- sd(cs[[v]])
  set(s, j = v, value = (s[[v]] - mu) / sg)
}

s[, D := get(HEAT) * consol30]
for (v in c("coast", "forest")) set(s, j = paste0("h_", v), value = s[[HEAT]] * s[[v]])
s[, gkey := as.character(year)]
keys <- sort(unique(s$gkey)); K <- length(keys)
cat(sprintf("\n標本: %d行 / %d市町村 / %d年 (%s-%s)\n",
            nrow(s), uniqueN(s$city_id), K, keys[1], keys[K]))

# -----------------------------------------------------------------------------
# 1. 本番スペック R3（T20 のはしごの R3。プール式と年飽和式を両方回す）
# -----------------------------------------------------------------------------
slopes  <- c(HEAT, "h_coast", "h_forest")       # 年飽和させる傾き
levvars <- c("consol30", "coast", "forest")     # 年飽和させる水準（T15 の政策項を含む）

lv <- character(0)
for (v in levvars) for (k in keys[-1]) {
  n <- sprintf("L%s_%s", v, k)
  set(s, j = n, value = s[[v]] * as.integer(s$gkey == k))
  lv <- c(lv, n)
}

rs <- c("i(gkey, D)", sprintf("i(gkey, %s)", slopes), lv, CTRL)
ms <- feols(as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rs, collapse = " + "))),
            s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"),
            warn = FALSE, notes = FALSE)
rp <- c("D", slopes, lv, CTRL)
mp <- feols(as.formula(sprintf("log_yield ~ %s | city_id + year", paste(rp, collapse = " + "))),
            s, cluster = ~city_id, fixef.rm = "none", ssc = ssc(K.fixef = "full"),
            warn = FALSE, notes = FALSE)
dpool  <- unname(coef(mp)["D"])
sepool <- unname(se(mp)["D"])

# 年飽和式の係数位置。名前は "gkey::2010:D" の形なので4要素に割って拾う
gx <- function(m, v) {
  p <- strsplit(names(coef(m)), ":", fixed = TRUE)
  i <- which(vapply(p, function(z) length(z) >= 4 && z[1] == "gkey" && z[4] == v, logical(1)))
  i[order(vapply(p[i], function(z) z[3], character(1)))]
}
idD <- gx(ms, "D"); idH <- gx(ms, HEAT); idG <- gx(ms, "h_coast"); idF <- gx(ms, "h_forest")
stopifnot(length(idD) == K, length(idH) == K, length(idG) == K, length(idF) == K)

# Dtil = D を「他の全プール回帰子」に射影した残差（FE は demean で落とす）
M  <- as.data.frame(demean(s[, ..rp], s[, .(city_id, year)]))
Dt <- as.numeric(residuals(lm(
        as.formula(paste("D ~", paste(c(slopes, lv, CTRL), collapse = " + "), "- 1")), data = M)))
ss <- sum(Dt^2); set(s, j = "Dtil", value = Dt)

wof <- function(v) s[, .(w = sum(Dtil * get(v)) / ss), by = gkey][order(gkey)]$w
wD <- wof("D"); wH <- wof(HEAT); wG <- wof("h_coast"); wF <- wof("h_forest")

b  <- coef(ms)
Vc <- vcov(ms)
Vk <- tryCatch(suppressWarnings(vcov_conley(ms, lat = ~lat_deg, lon = ~lon_deg,
          cutoff = CONLEY, ssc = ssc(K.fixef = "full"))), error = function(e) NULL)
if (is.null(Vk)) cat("  ! Conley VCV を構築できなかった。クラスターのみで出す\n")

# 係数空間での重みベクトル。target・主項・汚染・ギャップを同じ土俵に乗せる
vec <- function(...) { a <- rep(0, length(b)); for (p in list(...)) a[p$i] <- a[p$i] + p$w; a }
P   <- function(i, w) list(i = i, w = w)
DG  <- new.env()                 # Conley の診断。T12 §2 の基準を線形結合に当てる
DG$neg <- 0L; DG$small <- character(0)
est <- function(a, tag = NA_character_) {
  pt  <- sum(a * b)
  se1 <- sqrt(max(as.numeric(t(a) %*% Vc %*% a), 0))
  q2  <- if (is.null(Vk)) NA_real_ else as.numeric(t(a) %*% Vk %*% a)
  if (!is.na(q2) && q2 < 0) DG$neg <- DG$neg + 1L            # 二次形式が負（非正定値の露出）
  se2 <- if (is.na(q2)) NA_real_ else sqrt(max(q2, 0))
  if (!is.na(se2) && se2 < se1 && !is.na(tag)) DG$small <- c(DG$small, tag)
  c(pt = pt, se = se1, t = pt / se1, lo = pt - 1.96 * se1, hi = pt + 1.96 * se1,
    se_c = se2, t_c = pt / se2, lo_c = pt - 1.96 * se2, hi_c = pt + 1.96 * se2)
}

a_main <- vec(P(idD, wD))
a_cont <- vec(P(idH, wH), P(idG, wG), P(idF, wF))
a_pool <- vec(P(idD, wD), P(idH, wH), P(idG, wG), P(idF, wF))
r_main <- est(a_main, "主項"); r_cont <- est(a_cont, "汚染 合計"); r_pool <- est(a_pool, "プール")

cat(sprintf("\nプール delta_hat = %+.6f (SE %.6f)   恒等式の誤差 %.2e   sum w_D = %.10f\n",
            dpool, sepool, dpool - r_pool["pt"], sum(wD)))

# -----------------------------------------------------------------------------
# 2. target estimand の定義（曝露分布だけで決める。左辺は見ない）
# -----------------------------------------------------------------------------
hbar <- s[, .(h = mean(get(HEAT))), by = gkey][order(gkey)]$h   # 年平均 HD35
top  <- order(hbar, decreasing = TRUE)[1:NTOP]
med  <- order(hbar)[ceiling(K / 2)]

a_eq  <- rep(1 / K, K)                                        # (A) 等重み
a_exp <- hbar / sum(hbar)                                     # (B) 曝露加重
a_ext <- { z <- rep(0, K); z[top] <- 1 / NTOP; z }             # (C) 高曝露NTOP年の平均
a_10  <- { z <- rep(0, K); z[which(keys == "2010")] <- 1; z }  # (D) 2010年
a_med <- { z <- rep(0, K); z[med] <- 1; z }                    # (E) 中位曝露年（参考）

TAR <- list("(A) 等重み平均 1/24"  = a_eq,
            "(B) 曝露加重平均"     = a_exp,
            "(C) 高曝露5年の平均"  = a_ext,
            "(D) delta_2010"       = a_10,
            "(E) 中位曝露年 参考"  = a_med)
LAB  <- c("A", "B", "C", "D", "E")
MAIN <- names(TAR)[1:4]   # 主結果。(E) は参考

cat(sprintf("\n年平均 HD35: 上位%d年 = %s\n", NTOP,
            paste(sprintf("%s(%.2f)", keys[top], hbar[top]), collapse = "  ")))
cat(sprintf("             中位年   = %s(%.2f)   全年平均 %.3f\n",
            keys[med], hbar[med], mean(hbar)))

# -----------------------------------------------------------------------------
# 3. target estimand の推定
# -----------------------------------------------------------------------------
cat("\n################ A. target estimand ################\n")
cat(sprintf("%-22s %11s %10s %7s %24s %7s\n",
            "target", "theta_hat", "SE(cl)", "t", "95%CI(cl)", "t(Cnl)"))
RES <- list()
for (nm in names(TAR)) {
  r <- est(vec(P(idD, TAR[[nm]])), nm); RES[[nm]] <- r
  cat(sprintf("%-22s %+11.6f %10.6f %+7.2f  [%+.5f, %+.5f] %+7.2f\n",
              nm, r["pt"], r["se"], r["t"], r["lo"], r["hi"], r["t_c"]))
}
for (z in list(list("(参考) 主項 sum wD dt", r_main), list("(参考) 汚染 合計", r_cont),
               list("プール delta_hat", r_pool))) {
  r <- z[[2]]
  cat(sprintf("%-22s %+11.6f %10.6f %+7.2f  [%+.5f, %+.5f] %+7.2f\n",
              z[[1]], r["pt"], r["se"], r["t"], r["lo"], r["hi"], r["t_c"]))
}

cat("\n精度の代価 SE(theta)/SE(プール):\n")
for (nm in names(TAR))
  cat(sprintf("  %-22s %5.2f 倍\n", nm, RES[[nm]]["se"] / r_pool["se"]))

# -----------------------------------------------------------------------------
# 4. ギャップ delta_pool - theta の分解
#    汚染は a に依存しないので全 target で同一。動くのは重みのずれだけ
# -----------------------------------------------------------------------------
cat("\n################ B. ギャップ delta_pool - theta ################\n")
cat(sprintf("%-22s %11s %11s %11s %11s %7s %24s\n",
            "target", "gap", "重みのずれ", "汚染", "theta_hat", "t(gap)", "95%CI(gap)"))
GAP <- list()
for (nm in names(TAR)) {
  at <- TAR[[nm]]
  rg <- est(vec(P(idD, wD - at), P(idH, wH), P(idG, wG), P(idF, wF)),
            paste("gap", nm))                                          # ギャップ全体
  rw <- est(vec(P(idD, wD - at)), paste("重みのずれ", nm))             # 重みのずれのみ
  rt <- RES[[nm]]
  GAP[[nm]] <- list(gap = unname(rg["pt"]), gw = unname(rw["pt"]), gc = unname(r_cont["pt"]),
                    se = unname(rg["se"]), t = unname(rg["t"]), t_c = unname(rg["t_c"]),
                    lo = unname(rg["lo"]), hi = unname(rg["hi"]),
                    gw_se = unname(rw["se"]), gw_t = unname(rw["t"]))
  cat(sprintf("%-22s %+11.6f %+11.6f %+11.6f %+11.6f %+7.2f  [%+.5f, %+.5f]\n",
              nm, rg["pt"], rw["pt"], r_cont["pt"], rt["pt"], rg["t"], rg["lo"], rg["hi"]))
}

cat("\n内訳（構成比はギャップ自身を分母にする。プール値は分母に使わない）:\n")
for (nm in names(TAR)) {
  g <- GAP[[nm]]; th <- unname(RES[[nm]]["pt"])
  cat(sprintf("  %-22s gap %+.6f = 重み %+.6f (%6.1f%%) + 汚染 %+.6f (%6.1f%%)   delta_pool/theta = %s\n",
              nm, g$gap, g$gw, 100 * g$gw / g$gap, g$gc, 100 * g$gc / g$gap,
              ifelse(abs(th) < 1e-9, "-", sprintf("%+.2f", dpool / th))))
}

# -----------------------------------------------------------------------------
# 4.5 Conley の診断（T12 §2 の基準を線形結合に当てる）
#   T12 §2：「Conley SE がクラスターより小さくなるのは数値的な破綻である」。
#   この基準は**負の重みを多く含む線形結合**では単純には使えないので、
#   固有値（非正定値性）と SE の大小を**分けて**出す。
#
#   fixest の Conley は一様カーネルで、**同一市町村の異時点ペア（距離0）を含む**。
#   合成データで確認した（2ユニットを cutoff 外に置き AR(0.9) を入れると
#   Conley SE x sqrt(G/(G-1)) がクラスター SE に一致する）。したがって
#   Conley は市町村クラスターを**入れ子に含み**、差は小標本補正
#   （クラスターは G/(G-1)、G=209 なので 1.005）と**市町村間の空間共分散**だけ。
#   よって Conley < クラスター は「空間共分散が差し引きで負」を意味し、
#   VCV が正定値である限り**数値的破綻ではない**。
#   ただしその量について Conley は**保守的でない**ので、見出しの数値は
#   クラスターで報告し、Conley は参考として併記する。
# -----------------------------------------------------------------------------
cat("\n################ D. Conley の診断 ################\n")
if (is.null(Vk)) {
  cat("Conley VCV なし\n")
} else {
  ev <- eigen((Vk + t(Vk)) / 2, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("Vk の固有値: 最小 %.3e / 最大 %.3e / 負の個数 %d / %d\n",
              min(ev), max(ev), sum(ev < 0), length(ev)))
  cat(sprintf("二次形式が負になった量: %d 個（SE は 0 に丸めている）\n", DG$neg))
  if (length(DG$small)) {
    cat(sprintf("Conley SE < クラスター SE になった量: %d 個\n", length(DG$small)))
    for (x in DG$small) cat(sprintf("     - %s\n", x))
    cat("   -> VCV は正定値なので破綻ではなく、市町村間の空間共分散が差し引きで負。\n")
    cat("      ただし Conley は保守的でないので、見出しはクラスターで報告する\n")
  } else {
    cat("Conley SE < クラスター SE の量はない\n")
  }
}

# -----------------------------------------------------------------------------
# 5. 出力
# -----------------------------------------------------------------------------
out <- rbindlist(lapply(names(TAR), function(nm) {
  r <- RES[[nm]]; g <- GAP[[nm]]
  data.table(target = nm, role = ifelse(nm %in% MAIN, "main", "reference"),
             theta = unname(r["pt"]), se_cl = unname(r["se"]), t_cl = unname(r["t"]),
             ci_lo = unname(r["lo"]), ci_hi = unname(r["hi"]),
             se_conley = unname(r["se_c"]), t_conley = unname(r["t_c"]),
             se_ratio = unname(r["se"] / r_pool["se"]),
             delta_pool = dpool, gap = g$gap, gap_weight = g$gw, gap_contam = g$gc,
             gap_se = g$se, gap_t = g$t, gap_t_conley = g$t_c,
             gap_lo = g$lo, gap_hi = g$hi,
             gap_weight_se = g$gw_se, gap_weight_t = g$gw_t)
}))
ref <- data.table(
  target = c("主項 sum wD delta_t", "汚染 合計", "プール delta_hat"), role = "component",
  theta = unname(c(r_main["pt"], r_cont["pt"], r_pool["pt"])),
  se_cl = unname(c(r_main["se"], r_cont["se"], r_pool["se"])),
  t_cl  = unname(c(r_main["t"],  r_cont["t"],  r_pool["t"])),
  ci_lo = unname(c(r_main["lo"], r_cont["lo"], r_pool["lo"])),
  ci_hi = unname(c(r_main["hi"], r_cont["hi"], r_pool["hi"])),
  se_conley = unname(c(r_main["se_c"], r_cont["se_c"], r_pool["se_c"])),
  t_conley  = unname(c(r_main["t_c"],  r_cont["t_c"],  r_pool["t_c"])),
  se_ratio = unname(c(r_main["se"], r_cont["se"], r_pool["se"]) / r_pool["se"]),
  delta_pool = dpool)
out <- rbindlist(list(out, ref), use.names = TRUE, fill = TRUE)

wt <- data.table(year = keys, hbar_HD35 = hbar,
                 w_D = wD, w_h = wH, w_G = wG, w_F = wF,
                 delta_t = unname(b[idD]), beta_t = unname(b[idH]),
                 gamma_t = unname(b[idG]), eta_t = unname(b[idF]))
for (j in seq_along(TAR)) {
  set(wt, j = paste0("a_",    LAB[j]), value = TAR[[j]])
  set(wt, j = paste0("gapw_", LAB[j]), value = wD - TAR[[j]])
}
wt[, contrib_D := w_D * delta_t]

cat("\n################ C. 年別の重み（抜粋） ################\n")
print(wt[, .(year, hbar_HD35 = round(hbar_HD35, 2), w_D = round(w_D, 4),
             delta_t = round(delta_t, 5), contrib_D = round(contrib_D, 6),
             a_B = round(a_B, 4), gapw_A = round(gapw_A, 4))], row.names = FALSE)

cat("\n出力:\n")
sv(out, "t21_target_estimand")
sv(wt,  "t21_target_weights")

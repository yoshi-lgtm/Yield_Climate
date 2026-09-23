"""
二元固定効果(TWFE)推定の最小実装。

R の fixest / python の linearmodels がこの環境に入っていないため、
交互射影による within 変換 + クラスターロバスト分散 を自前で持つ。
fixest::feols(y ~ X | city + year, cluster = ~city) と同じ数値を返すことを
意図している（自由度修正も fixest の既定に合わせる）。
"""
import numpy as np
import pandas as pd
from scipy import stats


def _demean(arr, groups, tol=1e-10, maxiter=200):
    """交互射影で複数の固定効果を除去する。arr: (n, k), groups: list of (n,) int codes"""
    x = np.asarray(arr, dtype=float).copy()
    if x.ndim == 1:
        x = x[:, None]
    if len(groups) == 1:
        g = groups[0]
        means = np.zeros((g.max() + 1, x.shape[1]))
        counts = np.bincount(g, minlength=g.max() + 1)[:, None]
        np.add.at(means, g, x)
        return x - (means / counts)[g]

    for _ in range(maxiter):
        prev = x.copy()
        for g in groups:
            n_g = g.max() + 1
            means = np.zeros((n_g, x.shape[1]))
            counts = np.bincount(g, minlength=n_g)[:, None]
            np.add.at(means, g, x)
            x -= (means / counts)[g]
        if np.max(np.abs(x - prev)) < tol:
            break
    return x


def feols(df, y, X, fe, cluster, weights=None):
    """TWFE + クラスターロバストSE。

    df      : DataFrame
    y       : 被説明変数名
    X       : 説明変数名のリスト
    fe      : 固定効果にする列名のリスト（例 ["city_id", "year"]）
    cluster : クラスター変数名（例 "city_id"）

    戻り値: dict(coef, se, t, p, n, n_fe, within_r2, resid)
    """
    # cluster が fe と同じ列を指すことが多いので重複を落とす（列重複は factorize を壊す）
    cols = list(dict.fromkeys([y] + list(X) + list(fe) + [cluster]))
    d = df[cols].dropna().copy()
    n = len(d)

    yv = d[y].to_numpy(dtype=float)
    Xv = d[X].to_numpy(dtype=float)
    groups = [pd.factorize(d[f])[0] for f in fe]

    yt = _demean(yv, groups).ravel()
    Xt = _demean(Xv, groups)

    XtX = Xt.T @ Xt
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ (Xt.T @ yt)
    resid = yt - Xt @ beta

    # within R^2: 固定効果除去後の説明力
    within_r2 = 1.0 - resid.var(ddof=0) / yt.var(ddof=0)

    # クラスターロバスト分散（fixest 既定の自由度修正）
    cl = pd.factorize(d[cluster])[0]
    G = cl.max() + 1
    k = Xt.shape[1]
    n_fe_params = sum(g.max() + 1 for g in groups) - (len(groups) - 1)

    meat = np.zeros((k, k))
    u = Xt * resid[:, None]
    sums = np.zeros((G, k))
    np.add.at(sums, cl, u)
    meat = sums.T @ sums

    dof_adj = (G / (G - 1)) * ((n - 1) / (n - k - n_fe_params))
    V = XtX_inv @ meat @ XtX_inv * dof_adj

    se = np.sqrt(np.diag(V))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df=G - 1)

    return {
        "vars": list(X),
        "coef": beta,
        "se": se,
        "t": t,
        "p": p,
        "n": n,
        "n_cluster": G,
        "n_cities": int(pd.Series(d[fe[0]]).nunique()),
        "within_r2": within_r2,
        "vcov": V,
        "resid": resid,
    }


def stars(p):
    return "***" if p < 0.01 else ("**" if p < 0.05 else ("*" if p < 0.1 else ""))


def summary(res, title="", digits=5):
    lines = []
    if title:
        lines.append("=== " + title + " ===")
    lines.append(f"{'variable':<26}{'coef':>12}{'se':>12}{'t':>8}   ")
    for i, v in enumerate(res["vars"]):
        lines.append(f"{v:<26}{res['coef'][i]:>12.{digits}f}{res['se'][i]:>12.{digits}f}"
                     f"{res['t'][i]:>8.2f}   {stars(res['p'][i])}")
    lines.append(f"n = {res['n']}, cities = {res['n_cities']}, "
                 f"clusters = {res['n_cluster']}, within R2 = {res['within_r2']:.4f}")
    return "\n".join(lines)

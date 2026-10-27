#!/usr/bin/env python3
"""Power of the two-sided paired t-test, computed by the checker's own quadrature (D61).

With d_j = -s dT + e_j, e normal with sd sd_e, the statistic is T' = (Z + delta) / sqrt(V / df), Z standard normal,
V chi-square with df = n - 1 degrees of freedom, delta = sqrt(n) |s| dT / sd_e. The two-sided power at level alpha is
    P(T' > t_c) + P(T' < -t_c) = E_x[ Phi(delta - t_c x / sqrt(df)) + Phi(-delta - t_c x / sqrt(df)) ],
x = sqrt(V) chi-distributed with df degrees of freedom, t_c the upper alpha/2 point of Student's t with df degrees of
freedom. The expectation is integrated with scipy.integrate.quad over x on a window of +-14 standard deviations about
the chi mode (the chi density outside it is below 1e-40 of its peak), with the integrand built from log-densities, so
it does not overflow at large df or delta.
"""
import math

import numpy as np
from scipy import integrate, special, stats

ALPHA = 0.05
POWER = 0.80


def _logchi_pdf(x, k):
    # chi density with k degrees of freedom: x^(k-1) exp(-x^2/2) / (2^(k/2 - 1) Gamma(k/2))
    return (k - 1) * np.log(x) - 0.5 * x * x - (0.5 * k - 1) * math.log(2.0) - special.gammaln(0.5 * k)


def power_quad(n, effect, alpha=ALPHA, lower_tail=True):
    """effect = |s| dT / sd_e; returns (power, upper-tail part, lower-tail part)."""
    df = n - 1
    delta = math.sqrt(n) * effect
    tc = stats.t.ppf(1.0 - alpha / 2.0, df)
    rt = math.sqrt(df)
    mode = math.sqrt(df - 1) if df > 1 else 0.0
    lo = max(1e-300, mode - 14.0)
    hi = mode + 14.0

    def f_up(x):
        return math.exp(_logchi_pdf(x, df) + special.log_ndtr(delta - tc * x / rt))

    def f_lo(x):
        return math.exp(_logchi_pdf(x, df) + special.log_ndtr(-delta - tc * x / rt))

    pts = [p for p in (mode, (delta / tc) * rt) if lo < p < hi]
    up = integrate.quad(f_up, lo, hi, points=pts or None, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    low = integrate.quad(f_lo, lo, hi, points=pts or None, epsabs=1e-14, epsrel=1e-12, limit=400)[0] if lower_tail else 0.0
    return up + low, up, low


def power_nct(n, effect, alpha=ALPHA):
    df = n - 1
    delta = math.sqrt(n) * effect
    tc = stats.t.ppf(1.0 - alpha / 2.0, df)
    return float(stats.nct.sf(tc, df, delta) + stats.nct.cdf(-tc, df, delta))


def n_required(effect, alpha=ALPHA, target=POWER, lower_tail=True, nmax=10 ** 7):
    """Smallest integer n >= 2 with power >= target, by doubling then bisection on the monotone power."""
    pw = lambda n: power_quad(n, effect, alpha, lower_tail)[0]
    if pw(2) >= target:
        return 2
    lo, hi = 2, 4
    while pw(hi) < target:
        lo, hi = hi, hi * 2
        if hi > nmax:
            return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if pw(mid) >= target:
            hi = mid
        else:
            lo = mid
    assert pw(hi) >= target and pw(hi - 1) < target or hi == 2
    return hi


def lognormal_shape(m, q):
    """|e| log-normal with mean m and 95th percentile q (D61): sigma^2/2 - z sigma + ln(q/m) = 0, smaller root."""
    z = stats.norm.ppf(0.95)
    disc = z * z - 2.0 * math.log(q / m)
    s_small = z - math.sqrt(disc)
    s_large = z + math.sqrt(disc)
    out = {}
    for name, s in (("small", s_small), ("large", s_large)):
        mu = math.log(m) - s * s / 2.0
        out[name] = {"sigma": s, "mu": mu, "median": math.exp(mu), "mean_check": math.exp(mu + s * s / 2.0),
                     "p95_check": math.exp(mu + z * s), "E_e2": math.exp(2 * mu + 2 * s * s),
                     "sd_e": math.sqrt(math.exp(2 * mu + 2 * s * s))}
    return out


def normal_sd(m):
    return m * math.sqrt(math.pi / 2.0)


if __name__ == "__main__":
    # self-checks against scipy's nct and a textbook value
    for n, eff in ((2, 5.0), (3, 2.0), (10, 1.0), (25, 0.585), (100, 0.3), (1000, 0.09), (5000, 0.04)):
        print(n, eff, power_quad(n, eff), power_nct(n, eff))

#!/usr/bin/env python3
"""P12: does the sign of D_(N,l) (Corollary 4) depend on the Tietz scale a?

Pass line fixed before the run: for a in {0.40, 0.45, 0.50, 0.536795, 0.60, 0.65, 0.70} and (pi/8)^(2/3),
every D_(N,l) estimate is positive for every consecutive pair l, l+1 in every group 3 <= N <= 8.
Edges of the passing range are bisected if the pass line fails inside the grid. Rigorous certificate
(tf_certificate machinery, not modified, a set by assigning its module constant) only at the two ends of
the passing range found in the grid, with --certify.
Estimate: D = 2 int_0^1 h (p_l - p_(l+1)) dx, h = 1 - G(rho), G = (1+rho)^2 (P+Q)/(2A), P = phi(rho/a),
Q = phi(1/(a rho))/rho^2, rho = (1-x)/(1+x), exact rational densities (tf_certificate.densities), phi from
tf_accurate (double). A is fixed at 1.012783: it multiplies G by a constant, so D scales by a positive factor.
Run: nice -n 19 env OMP_NUM_THREADS=1 python3 a_scan.py [--certify]
"""
import os, sys, math, time
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tf_certificate as tc
import tf_accurate

_tf = tf_accurate.default()
A = 1.012783


def h_of(x, a):
    rho = (1 - x) / (1 + x)
    if rho < 1e-300:
        return 1 - 0.5 / A
    P = float(_tf.phi(rho / a))
    Q = float(_tf.phi(1 / (a * rho))) / rho ** 2
    return 1 - (1 + rho) ** 2 * (P + Q) / (2 * A)


def pairs():
    out = []
    for N in range(3, 9):
        for d in tc.densities(N):
            c = [float(v) for v in d['p'].c]
            out.append((N, d['l'], np.polynomial.polynomial.Polynomial(c)))
    return out


def Dest(a, prs):
    res = []
    for N, l, p in prs:
        f = lambda x: h_of(x, a) * p(x)
        pts = [0.5, 0.9, 0.99, 0.999, 0.9999]
        v = quad(f, 0, 1 - 1e-12, points=pts, epsabs=1e-14, epsrel=1e-12, limit=500)[0]
        res.append((N, l, 2 * v))
    return res


def row(a, prs):
    r = Dest(a, prs)
    m = min(r, key=lambda t: t[2])
    return m, all(t[2] > 0 for t in r)


def main():
    prs = pairs()
    grid = [0.40, 0.45, 0.50, 0.536795, (math.pi / 8) ** (2 / 3), 0.60, 0.65, 0.70]
    rows = {}
    for a in grid:
        m, ok = row(a, prs)
        rows[a] = (m, ok)
        print(f"a = {a:.6f}  min D = {m[2]:.6e} at (N,l) = ({m[0]},{m[1]})  all positive: {ok}", flush=True)
    passing = [a for a in grid if rows[a][1]]
    print("grid pass line:", "PASS" if len(passing) == len(grid) else "FAIL", "passing a:", passing)
    if len(passing) < len(grid):
        # bisect each edge of the passing set against the nearest failing grid point
        srt = sorted(grid)
        for lo_, hi_ in zip(srt, srt[1:]):
            if rows[lo_][1] != rows[hi_][1]:
                good, bad = (lo_, hi_) if rows[lo_][1] else (hi_, lo_)
                for _ in range(20):
                    mid = (good + bad) / 2
                    if row(mid, prs)[1]:
                        good = mid
                    else:
                        bad = mid
                print(f"edge between {lo_} and {hi_}: last passing a = {good:.6f}, first failing {bad:.6f}", flush=True)
    for a in [0.20, 0.25, 0.30, 0.35, 0.675, 0.69]:   # informational probes outside the fixed grid
        m, ok = row(a, prs)
        print(f"probe a = {a:.3f}  min D = {m[2]:.6e} at (N,l) = ({m[0]},{m[1]})  all positive: {ok}", flush=True)
    if '--certify' in sys.argv:
        ends = [float(sys.argv[sys.argv.index('--certify') + 1]), float(sys.argv[sys.argv.index('--certify') + 2])]
        sl, sh = tc.certify_sandwich(0)
        W = tc.Weak(sl, sh)
        groups = [(N, tc.densities(N)) for N in range(3, 9)]
        for a in ends:
            tc.a_TIETZ = tc.iv.mpf(repr(a))
            t0 = time.time()
            r = tc.certify(W, groups, '1/1000')
            lbs = [(e['cross'][0] + e['pole'][0], e['N'], e['l']) for e in r['res']]
            mn = min(lbs, key=lambda t: t[0])
            print(f"CERT a = {a}: min certified LB(D) = {tc.mp.nstr(mn[0], 6)} at (N,l) = ({mn[1]},{mn[2]}), all LB > 0: {all(t[0] > 0 for t in lbs)}  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == '__main__':
    main()

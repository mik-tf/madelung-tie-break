#!/usr/bin/env python3
"""Proposition 5: non-perturbative path from the nearest Tietz field to the TF field.

Pass line fixed in advance: for every N, 3 <= N <= 12, every l with
l+1 <= floor((N-1)/2) and every grid eps in (0, 1]: Z_(N,l+1)(eps) < Z_(N,l)(eps), smallest relative gap
(Z_l - Z_(l+1))/Z_l above 1e-6. Any crossing inside a group fails. Prints PASS/FAIL, on FAIL the first (N,l,eps).

Setup. Weight family w (1 + eps g_TF) = (1-eps) A T(s) + eps phi_TF(s), times s, T = (1+a s)^-2, (a, A) the nearest
Tietz field in the weighted L2 norm (recomputed here as in tietz_extremal.py). Radial equation in
sigma = ln(s)/2 (Poschl-Teller variable, as review_tiebreak_check.py part 4):
   y'' = (kap^2 - 2 Z W(sigma, eps)) y,  kap = 2l+1,  W = (4a/A) s [ (1-eps) A T + eps phi_TF ],  s = e^(2 sigma).
At eps = 0, W = sech^2(sigma - sigma0), Z_(N,l) = N(N+1)/2. Z_(N,l) = the (N-2l)-th zero-energy coupling in sector l
(radial quantum number n-l-1 = N-2l-1). Solved by two-sided shooting with the 4th order Numerov scheme,
Wronskian at sigma = 0, Illinois root finding, warm-started along eps. Step sizes h = 1/256 and 1/512 (error estimate),
domain |sigma| <= 12. Double precision (30 digits would be possible; the gap line 1e-6 is 8 orders above double errors).
Run: nice -n 19 env OMP_NUM_THREADS=1 python3 tf_path_check.py [neps=41]
"""
import os, sys, time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tf_accurate

NEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 41
XD = 12.0
NMAX = 12
t0 = time.time()
_tf = tf_accurate.default()
phi = lambda s: _tf.phi(s)

# ---- nearest Tietz field (as tietz_extremal.py c)
def ints(a):
    t = lambda s: (1 + a * s) ** -2
    pp = quad(lambda s: s * phi(s) * t(s), 0, 12, epsabs=1e-14, epsrel=1e-13, limit=400)[0] \
        + quad(lambda s: s * phi(s) * t(s), 12, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    tt = quad(lambda s: s * t(s) ** 2, 0, np.inf, epsabs=1e-14, epsrel=1e-13)[0]
    ff = quad(lambda s: s * phi(s) ** 2, 0, 12, epsabs=1e-14, epsrel=1e-13, limit=400)[0] \
        + quad(lambda s: s * phi(s) ** 2, 12, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    A = pp / tt
    return ff - 2 * A * pp + A * A * tt, A
res = minimize(lambda p: ints(p[0])[0], [0.5367], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-18))
a = float(res.x[0]); A = ints(a)[1]
print(f"nearest Tietz field: a = {a:.7f}  A = {A:.7f}   ({time.time()-t0:.0f}s)", flush=True)

class Grid:
    def __init__(self, h, X=XD):
        self.h = h
        n = int(round(X / h))
        self.sig = np.arange(-n, n + 1) * h
        self.M = n                      # index of sigma = 0
        s = np.exp(2 * self.sig)
        self.WT = 4 * a * s / (1 + a * s) ** 2
        self.WF = (4 * a / A) * s * np.array([phi(v) for v in s])

def D(Z, l, eps, G):
    """normalised Wronskian at sigma = 0 for coupling array Z"""
    Z = np.atleast_1d(np.asarray(Z, float))
    kap2 = (2 * l + 1) ** 2
    h = G.h; h2 = h * h / 12.0
    W = (1 - eps) * G.WT + eps * G.WF
    n = len(W); M = G.M; kap = 2 * l + 1
    def sweep(order, stop):
        # q_j = kap2 - 2 Z W_j ; idx = order list of node indices
        y0 = np.ones_like(Z); y1 = np.full_like(Z, np.exp(kap * h))
        ys = {order[0]: y0, order[1]: y1}
        qm = kap2 - 2 * Z * W[order[0]]; qn = kap2 - 2 * Z * W[order[1]]
        ym, yn = y0, y1
        for k in range(1, len(order) - 1):
            qp = kap2 - 2 * Z * W[order[k + 1]]
            yp = (2 * (1 + 5 * h2 * qn) * yn - (1 - h2 * qm) * ym) / (1 - h2 * qp)
            ym, yn, qm, qn = yn, yp, qn, qp
            if order[k + 1] in stop:
                ys[order[k + 1]] = yp
        return ys
    left = sweep(list(range(0, M + 3)), set(range(M - 2, M + 3)))
    right = sweep(list(range(n - 1, M - 3, -1)), set(range(M - 2, M + 3)))
    def dy(ys):
        return (-ys[M + 2] + 8 * ys[M + 1] - 8 * ys[M - 1] + ys[M - 2]) / (12 * h)
    nl = np.sqrt(sum(left[j] ** 2 for j in range(M - 2, M + 3)))
    nr = np.sqrt(sum(right[j] ** 2 for j in range(M - 2, M + 3)))
    return (dy(left) * right[M] - dy(right) * left[M]) / (nl * nr)

def roots(l, eps, G, guess, width):
    """refine one root near each guess; bracket by widening, then Illinois."""
    guess = np.asarray(guess, float)
    lo = guess - width; hi = guess + width
    flo = D(lo, l, eps, G); fhi = D(hi, l, eps, G)
    for _ in range(10):
        bad = flo * fhi > 0
        if not bad.any(): break
        lo[bad] -= width[bad] if hasattr(width, '__len__') else width
        hi[bad] += width[bad] if hasattr(width, '__len__') else width
        flo = D(lo, l, eps, G); fhi = D(hi, l, eps, G)
    else:
        raise RuntimeError(f"no bracket l={l} eps={eps}")
    side = np.zeros(len(guess), int)
    for it in range(60):
        z = hi - fhi * (hi - lo) / (fhi - flo)
        fz = D(z, l, eps, G)
        left = fz * fhi > 0       # root in [lo, z]
        hi_new = np.where(left, z, hi); fhi_new = np.where(left, fz, fhi)
        # Illinois: if same side kept twice, halve the retained value
        lo_new = np.where(left, lo, z); flo_new = np.where(left, flo, fz)
        flo_new = np.where(left & (side == -1), flo_new / 2, flo_new)
        fhi_new = np.where(~left & (side == 1), fhi_new / 2, fhi_new)
        side = np.where(left, -1, 1)
        lo, hi, flo, fhi = lo_new, hi_new, flo_new, fhi_new
        if np.max(hi - lo) < 1e-13 * 100: break
    return 0.5 * (lo + hi)

G1 = Grid(1 / 256); G2 = Grid(1 / 512)
print(f"grids built ({time.time()-t0:.0f}s)", flush=True)

LS = range(0, (NMAX - 1) // 2 + 1)            # l = 0..5
# root labelling: Z_(N,l), N = 2l+1+j, j = 0 .. NMAX-2l-1
def Nlist(l): return list(range(2 * l + 1, NMAX + 1))

Zs = {}      # eps -> {(N,l): (Z_h2, |Z_h1 - Z_h2|)}
def compute(eps, prev):
    out = {}
    for l in LS:
        Ns = Nlist(l)
        if prev is None:
            guess = np.array([N * (N + 1) / 2 for N in Ns]); width = np.full(len(Ns), 0.05)
        else:
            guess = np.array([prev[(N, l)][0] for N in Ns]); width = np.maximum(0.03, 0.03 * guess)
        r2 = roots(l, eps, G2, guess, width)
        r1 = roots(l, eps, G1, r2, np.full(len(Ns), 0.02))
        for N, z2, z1 in zip(Ns, r2, r1):
            out[(N, l)] = (z2, abs(z2 - z1))
    return out

epsgrid = [k / (NEPS - 1) for k in range(NEPS)]
prev = None
for e in epsgrid:
    # continuation: warm start from the previous eps (eps = 0 exact values)
    Zs[e] = compute(e, prev); prev = Zs[e]
    print(f"eps={e:.4f} done ({time.time()-t0:.0f}s)", flush=True)

def gaps(zd):
    g = {}
    for N in range(3, NMAX + 1):
        for l in range(0, (N - 1) // 2):       # l+1 <= floor((N-1)/2)
            zl, zl1 = zd[(N, l)][0], zd[(N, l + 1)][0]
            g[(N, l)] = (zl - zl1) / zl
    return g

def refine():
    es = sorted(Zs)
    new = []
    for e0, e1 in zip(es[:-1], es[1:]):
        if min(min(gaps(Zs[e0]).values()), min(gaps(Zs[e1]).values())) < 1e-3:
            new += [e0 + (e1 - e0) * k / 10 for k in range(1, 10)]
    for e in new:
        near = min(Zs, key=lambda x: abs(x - e))
        Zs[e] = compute(e, Zs[near])
    return len(new)
nref = refine()
print(f"refinement points added: {nref} ({time.time()-t0:.0f}s)")

# ---- checks
z0 = Zs[0.0]
dev0 = max(abs(z0[(N, l)][0] - N * (N + 1) / 2) for (N, l) in z0)
print(f"check eps=0: max |Z - N(N+1)/2| = {dev0:.2e} (double precision)")
c0 = 2 * a * Zs[1.0][(1, 0)][0] / A
print(f"check eps=1, N=1, l=0: c0 = 2 a Z / A = {c0:.8f} (expected 1.05926488)")
shoot = max(v[1] for zd in Zs.values() for v in zd.values())
shoot_rel = max(v[1] / v[0] for zd in Zs.values() for v in zd.values())
print(f"shooting error: max |Z(h=1/256) - Z(h=1/512)| = {shoot:.2e} (relative {shoot_rel:.2e}); 4th order, so error of the h=1/512 values about this/15")

# ---- pass line
es = sorted(e for e in Zs if e > 0)
fail = None; mingap = (1e9, None)
per_group = {}
for e in es:
    g = gaps(Zs[e])
    for (N, l), v in g.items():
        if fail is None and v <= 1e-6 and (N, l, e) is not None:
            fail = (N, l, e, v)
        if v < mingap[0]: mingap = (v, (N, l, e))
        if v < per_group.get(N, (1e9,))[0]: per_group[N] = (v, l, e)
print("smallest relative gap per group over eps in (0,1]  (gap, l vs l+1, eps):")
for N in sorted(per_group): print(f"  N={N}: {per_group[N][0]:.4e}  l={per_group[N][1]}  eps={per_group[N][2]:.4f}")
print(f"smallest overall: {mingap[0]:.4e} at (N,l,eps) = {mingap[1]}")
# strict crossing test, also between grid points via refined set
print("eps=1 order per group (TF ordering, lowest coupling first = binds first):")
for N in range(1, NMAX + 1):
    items = sorted(((Zs[1.0][(N, l)][0], l) for l in range(0, (N - 1) // 2 + 1)))
    print(f"  N={N}: " + "  ".join(f"l={l} (Z={z:.5f})" for z, l in items))
if fail is None:
    print("P8 PASS")
else:
    print(f"P8 FAIL: first (N, l, eps) = ({fail[0]}, {fail[1]}, {fail[2]:.4f}), relative gap {fail[3]:.3e}")
print(f"total time {time.time()-t0:.0f}s")

#!/usr/bin/env python3
"""the Tietz field as a sharp-inequality extremal.

Light: about 5 minutes at nice 19, one thread. mpmath (high precision) for (a), scipy for (b), (c).
Pass lines are fixed in advance:
 (a) E(mu) = Z0(mu) (int r v_mu^2 d3x)^(1/2) / (C/2), mu in (0,1]. PASS: minimum at
     mu = 1/2 with E = 1 to 1e-8, and E > 1 elsewhere.
 (b) E for the Thomas-Fermi field; recorded. <0.05 physics content, 0.05-0.2 weak, >0.2 none.
 (c) Tietz field nearest TF in this norm: its constant within 2 percent of 0.53625.

Setup (fixed before the run). DO family, R = 1, unit amplitude:
  v_mu(r) = 1 / (r^2 (r^mu + r^-mu)^2), so mu = 1/2 is the Tietz form 1/(r (1+r)^2).
  Z0 = inf int|grad psi|^2 / int 2 v psi^2 = smallest coupling with a zero-energy
  s-state. With x = ln r, u = r^(1/2) w: -w'' - (Z/2) sech^2(mu x) w = -(1/4) w, a
  Poschl-Teller well, giving Z0(mu) = (1 + 2 mu)/2 (checked numerically below).
  C from the sharp Sobolev constant S_4 = 8 pi / sqrt(6) through the KS map
  (d4u = (pi/4) d3x / r, |grad_4 f|^2 = 4 r |grad_3 psi|^2): C = S_4 / (2 sqrt(pi)).
  The extremal psi = 1/(1+r) is checked to attain C numerically.
TF field: V = -Z phi(r/b)/r, b = 1: v = phi(s)/s, norm^2 = 4 pi int s phi^2 ds,
  coupling c0 = smallest c with u'' = -c phi/s u, u(0)=0, flat at infinity (Z0 = c0/2).
"""
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, minimize

mp.mp.dps = 40
PI = mp.pi

# ------------------------------------------------------------------ (a)
S4 = 8 * PI / mp.sqrt(6)
C = S4 / (2 * mp.sqrt(PI))
print(f"S_4 = {mp.nstr(S4, 20)}   C = S_4/(2 sqrt(pi)) = {mp.nstr(C, 20)}")

# extremal check: psi = 1/(1+r)
I = 4 * PI * mp.quad(lambda r: r**2 * (1 / (1 + r)**2)**2, [0, 1, mp.inf])
J = 4 * PI * mp.quad(lambda r: r**2 * (1 / (1 + r))**4 / r, [0, 1, mp.inf])
print(f"extremal psi=1/(1+r): int|grad psi|^2 / (int psi^4/r)^(1/2) = {mp.nstr(I/mp.sqrt(J), 20)}  (C - that = {mp.nstr(C - I/mp.sqrt(J), 3)})")

def Z0_closed(mu):
    return (1 + 2 * mu) / 2

def Z0_numeric(mu):
    """Ground-state coupling by shooting in x = ln r, -w'' + w/4 = (Z/2) sech^2(mu x) w."""
    mu = float(mu)
    def end(Z):
        L = 60.0 / min(1.0, mu) if mu < 1 else 60.0
        f = lambda x, y: [y[1], (0.25 - 0.5 * Z / np.cosh(mu * x)**2) * y[0]]
        k = 0.5
        s = solve_ivp(f, (-L, 0), [1.0, k], method="DOP853", rtol=1e-12, atol=1e-30)
        a = s.y[:, -1]
        # even ground state: need w'(0) = 0
        return a[1] / a[0]
    return brentq(end, 0.3, 3.0, xtol=1e-13)

def norm(mu):
    mu = mp.mpf(mu)
    g = lambda x: (2 * mp.cosh(mu * x))**-4   # r^3 v^2 dr = g dx
    return mp.sqrt(4 * PI * mp.quad(g, [-mp.inf, 0, mp.inf]))

def E(mu):
    mu = mp.mpf(mu)
    return Z0_closed(mu) * norm(mu) / (C / 2)

print("check Z0(mu): numeric shooting vs (1+2mu)/2")
for mu in (0.25, 0.5, 0.75, 1.0):
    print(f"  mu={mu}: numeric {Z0_numeric(mu):.12f}  closed {float(Z0_closed(mu)):.12f}")
print("check norm(mu) vs sqrt(pi/(3 mu)) at mu=0.3:", mp.nstr(norm(0.3), 15), mp.nstr(mp.sqrt(PI / (3 * mp.mpf('0.3'))), 15))

mus = [mp.mpf(k) / 1000 for k in range(1, 1001)]
Es = [E(m) for m in mus]
imin = min(range(len(Es)), key=lambda i: Es[i])
mu_star = mp.findroot(lambda m: mp.diff(E, m), mp.mpf('0.5'))
E_star = E(mu_star)
others = min(Es[i] for i in range(len(Es)) if abs(mus[i] - mp.mpf('0.5')) > mp.mpf('0.005'))
print(f"(a) grid minimum at mu = {mp.nstr(mus[imin], 6)}; refined minimiser mu* = {mp.nstr(mu_star, 15)}")
print(f"    E(mu*) = {mp.nstr(E_star, 20)}   |E-1| = {mp.nstr(abs(E_star - 1), 3)}   smallest E off mu=1/2 (|mu-1/2|>0.005) = {mp.nstr(others, 10)}")
for m in ("0.1", "0.25", "0.4", "0.5", "0.6", "0.75", "1"):
    print(f"      E({m}) = {mp.nstr(E(m), 12)}")
passa = abs(E_star - 1) < mp.mpf('1e-8') and abs(mu_star - mp.mpf('0.5')) < mp.mpf('1e-8') and others > 1 and all(Es[i] > 1 or abs(mus[i] - mp.mpf('0.5')) < mp.mpf('1e-12') for i in range(len(Es)))
print("(a) " + ("PASS" if passa else "FAIL"))

# ------------------------------------------------------------------ (b), (c)
B_TF = 1.588071022611375
def phi_ode(smax=12.0):
    s0 = 1e-8
    y0 = [1 - B_TF * s0 + 4 / 3 * s0**1.5, -B_TF + 2 * s0**0.5]
    rhs = lambda s, y: [y[1], max(y[0], 0) ** 1.5 / np.sqrt(s)]
    return solve_ivp(rhs, (s0, smax), y0, method="DOP853", rtol=1e-13, atol=1e-16, dense_output=True)

SM = 12.0
sol = phi_ode(SM)
LAM = 0.772
def phi_somm(s):
    return (1 + (s / 12**(2 / 3))**LAM)**(-3 / LAM)
# tail: Sommerfeld form, rescaled to match the numeric value at SM (continuity)
tail_scale = sol.sol(SM)[0] / phi_somm(SM)
print(f"\nTF phi(12) numeric {sol.sol(SM)[0]:.6e}, Sommerfeld {phi_somm(SM):.6e} (tail rescale {tail_scale:.4f})")
def phi(s):
    s = np.asarray(s, float)
    return np.where(s <= SM, sol.sol(np.minimum(s, SM))[0], tail_scale * phi_somm(s))

import os, sys
TF_MODE = os.environ.get("TF_MODE", "accurate")   # "old": double-precision shot spliced to Sommerfeld at 12; default: tf_accurate
tail_fn = lambda s: tail_scale * phi_somm(s)
if TF_MODE == "accurate":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import tf_accurate
    _tf = tf_accurate.default()
    phi = lambda s: _tf.phi(s)
    tail_fn = lambda s: _tf.phi(s)
print("TF_MODE =", TF_MODE)

def tf_norm2():
    a = quad(lambda s: s * phi(s)**2, 0, SM, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    b = quad(lambda s: s * tail_fn(s)**2, SM, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    return 4 * np.pi * (a + b), 4 * np.pi * b

def c0_tf():
    S = 2e4
    def end(c):
        f = lambda s, y: [y[1], -c * phi(s) / s * y[0]]
        r = solve_ivp(f, (1e-10, S), [1e-10, 1.0], method="DOP853", rtol=1e-12, atol=1e-30,
                      first_step=1e-10, max_step=0.5)
        return r.y[1, -1]
    # smallest c with u'(inf) = 0; scan sign change from small c
    cs = np.linspace(0.5, 12, 24)
    vals = [end(c) for c in cs]
    for i in range(len(cs) - 1):
        if vals[i] * vals[i + 1] < 0:
            return brentq(end, cs[i], cs[i + 1], xtol=1e-12)
    raise RuntimeError("no threshold found")

n2, n2tail = tf_norm2()
c0 = c0_tf()
Z0tf = c0 / 2
Eb = Z0tf * np.sqrt(n2) / (float(C) / 2)
print(f"(b) TF (b=1 units): c0 = {c0:.8f}  Z0 = c0/2 = {Z0tf:.8f}  norm^2 = {n2:.8f} (tail beyond 12: {n2tail:.2e})")
print(f"    E_TF = {Eb:.6f}   E - 1 = {Eb - 1:.6f}")
cls = "physics content" if Eb - 1 < 0.05 else ("weak" if Eb - 1 < 0.2 else "no content")
print(f"(b) recorded: {cls}")

# (c) nearest Tietz field to TF: min over amplitude A and a of  int s (phi - A (1+a s)^-2)^2 ds
def dist2(a, A=None):
    t = lambda s: (1 + a * s)**-2
    pp = quad(lambda s: s * phi(s) * t(s), 0, SM, epsabs=1e-14, epsrel=1e-13, limit=400)[0] \
        + quad(lambda s: s * tail_fn(s) * t(s), SM, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    tt = quad(lambda s: s * t(s)**2, 0, np.inf, epsabs=1e-14, epsrel=1e-13)[0]
    ff = quad(lambda s: s * phi(s)**2, 0, SM, epsabs=1e-14, epsrel=1e-13, limit=400)[0] \
        + quad(lambda s: s * tail_fn(s)**2, SM, np.inf, epsabs=1e-14, epsrel=1e-13, limit=400)[0]
    if A is None:
        A = pp / tt
    return ff - 2 * A * pp + A * A * tt, A

res = minimize(lambda p: dist2(p[0])[0], [0.5], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-16))
a_best = res.x[0]
A_best = dist2(a_best)[1]
ff = n2 / (4 * np.pi)
print(f"(c) amplitude free: best a = {a_best:.6f}, A = {A_best:.6f}, relative L2 distance = {np.sqrt(dist2(a_best)[0] / ff):.4f}")
res1 = minimize(lambda p: dist2(p[0], 1.0)[0], [0.5], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-16))
print(f"    secondary, amplitude fixed A = 1 (neutral Tietz): best a = {res1.x[0]:.6f}")
dev = abs(a_best - 0.53625) / 0.53625
print(f"    deviation from Tietz a = 0.53625: {100 * dev:.2f} percent (amplitude free); {100 * abs(res1.x[0] - 0.53625) / 0.53625:.2f} percent (A = 1)")
print("(c) " + ("PASS" if dev <= 0.02 else "FAIL"))

"""Independent review checks of the tie-break proof (paper, Section 4 and Appendix A).

Written without reusing tiebreak_proof_check.py or tiebreak_dominance.py.
Run: OMP_NUM_THREADS=1 nice -n 19 python3 review_tiebreak_check.py   (about 1 minute)

Parts
 1. Lemma B with an abstract f: f'' and f''' are eliminated with the Gegenbauer ODE only
    (symbolic a, mu), so the identity is checked for every solution, not per degree.
 2. The densities rebuilt from the paper-1 R^3 eigenfunction psi_nlm (radial density
    w r^2 R_nl^2 dr, mapped by c = (1 - r)/(1 + r)), not from the S^4 block formula.
    From these, H_l(t) = int_t^1 (F_(l+1) - F_l) exactly, compared with closed form (5),
    and positivity on [0, 1) by exact real-root counting, N <= 20.
 3. Lemma C: Q_j has no real root (exact count), k <= 20, and a counterexample showing
    that the real-zero hypothesis is needed.
 4. First-order step, independent of the S^4 lift: the R^3 radial problem in Poschl-Teller
    form phi'' - (2l+1)^2 phi + 2 Z (1 + eps g(c)) sech^2(sigma) phi = 0, c = -tanh(sigma),
    solved by shooting for eps = +-1e-5; dZ/deps compared with -Z_N E_(q_2l)[g].
"""
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer

ok = True
def report(name, cond):
    global ok
    ok = ok and bool(cond)
    print(("PASS " if cond else "FAIL ") + name)

# ---------------- 1. Lemma B, abstract f ----------------
c, a, mu = sp.symbols('c a mu')
F0, F1 = sp.symbols('F0 F1')           # f, f'
s = 1 - c**2
F2 = (2*a*c*F1 - mu*F0)/s              # from s f'' - 2 a c f' + mu f = 0
def D(expr):                            # total derivative along solutions of the ODE
    return sp.diff(expr, c) + sp.diff(expr, F0)*F1 + sp.diff(expr, F1)*F2
V = s**(a+1)*(F1**2 - (a-1)/a*F0*F2)
lhs = D(V)
rhs = -(2/a)*(mu + a*(a-1))*s**a*F0*F1
report("Lemma B: V' = -(2/a)(mu + a(a-1)) s^a f f' for every solution, symbolic a, mu",
       sp.simplify(sp.powsimp(sp.expand((lhs - rhs)/s**(a-1)), force=True)) == 0)
# Lemma A, abstract
lhsA = D(-s**a*F0*F1/mu)
rhsA = s**(a-1)*F0**2 - s**a*F1**2/mu
report("Lemma A: d/dt[-(1/mu) s^a f f'] = s^(a-1) f^2 - s^a f'^2/mu, symbolic",
       sp.simplify((lhsA - rhsA)/s**(a-1)) == 0)
k_ = sp.symbols('k')
report("mu + a(a-1) = (k+1)(k+2) with a = j+2, mu = (k-j)(k+j+3)",
       sp.expand(((k_-sp.Symbol('j'))*(k_+sp.Symbol('j')+3) + (sp.Symbol('j')+2)*(sp.Symbol('j')+1))
                 - (k_+1)*(k_+2)) == 0)

# ---------------- 2. densities from the R^3 eigenfunction, closed form (5) ----------------
t, r = sp.symbols('t r', real=True)
def q_from_R3(k, l):
    """Density of c for psi_(n l m), n_r = k - 2l, R = 1, built in R^3 and mapped to c."""
    nr = k - 2*l
    Rnl = r**l*(1 + r)**(-(2*l+1))*sp.gegenbauer(nr, sp.Rational(4*l+3, 2), (1 - r)/(1 + r))
    w = 2/(r*(1 + r)**2)
    dens_r = w*r**2*Rnl**2
    rc = (1 - c)/(1 + c)
    dens_c = sp.simplify(sp.factor(dens_r.subs(r, rc)*sp.diff(rc, c)*(-1)))  # dr = -(...) dc
    dens_c = sp.Poly(sp.expand(sp.cancel(dens_c)), c)
    Z = sp.integrate(dens_c.as_expr(), (c, -1, 1))
    return sp.expand(dens_c.as_expr()/Z)

def H_direct(k, l):
    p0 = q_from_R3(k, l); p1 = q_from_R3(k, l + 1)
    u = sp.symbols('u')
    F0_ = sp.integrate(p0, (c, -t, t)); F1_ = sp.integrate(p1, (c, -t, t))
    return sp.expand(sp.integrate((F1_ - F0_).subs(t, u), (u, t, 1)))

def H_closed(k, l):
    tot = 0
    for jj in (2*l, 2*l + 1):
        lam = sp.Rational(2*jj + 3, 2); m = k - jj; aa = jj + 2
        f = sp.gegenbauer(m, lam, t)
        st = 1 - t**2
        Zj = sp.integrate((st**(jj+1)*f**2).subs(t, c), (c, -1, 1))
        muj = (k - jj)*(k + jj + 3)
        Q = sp.expand(sp.diff(f, t)**2 - sp.Rational(jj + 1, jj + 2)*f*sp.diff(f, t, 2))
        tot += sp.Rational(aa, 1)/((k + 1)*(k + 2)*muj*Zj)*st**(jj+3)*Q
    return sp.expand(tot)

npairs = 0; bad = 0
for N in range(2, 21):
    k = N - 1
    for l in range(0, k//2):
        Hd = H_direct(k, l); Hc = H_closed(k, l)
        same = sp.expand(Hd - Hc) == 0
        # positivity on [0,1): H(1)=0 and no root in [0,1), H(0) > 0
        P = sp.Poly(Hd, t)
        nroot = sp.Poly(sp.sqf_part(P.as_expr()), t).count_roots(0, 1)
        has_one = P.eval(1) == 0
        inner = nroot - (1 if has_one else 0)
        pos = (inner == 0) and P.eval(0) > 0 and has_one
        npairs += 1
        if not (same and pos):
            bad += 1; print("  problem at", N, l, same, pos)
report(f"(5) equals int_t^1 (F_(l+1) - F_l) built from psi_nlm in R^3, and H_l > 0 on [0,1), H_l(1) = 0: {npairs} pairs, N <= 20", bad == 0)

# ---------------- 3. Lemma C ----------------
badC = 0
for k in range(1, 21):
    for jj in range(0, k):
        f = sp.gegenbauer(k - jj, sp.Rational(2*jj + 3, 2), t)
        Q = sp.Poly(sp.expand(sp.diff(f, t)**2 - sp.Rational(jj + 1, jj + 2)*f*sp.diff(f, t, 2)), t)
        if Q.count_roots() != 0 or Q.eval(0) <= 0:
            badC += 1
report("Q_j has no real root and Q_j(0) > 0 for all 0 <= j < k <= 20", badC == 0)
fx = t**2 + 1
report("hypothesis needed: f = t^2 + 1 (no real zeros) gives f'^2 - f f'' = 2t^2 - 2 < 0 at t = 0",
       sp.expand(sp.diff(fx, t)**2 - fx*sp.diff(fx, t, 2)).subs(t, 0) < 0)

# ---------------- 4. first-order step, R^3 radial problem, no lift ----------------
def gfun(cc):          # weight change w -> w(1 + eps g); even part convex decreasing in |c|, odd part 0.25 c
    return -(0.6*cc**2 + 0.4*np.abs(cc)**3) - 0.25*cc   # odd part 0.25 c must drop out

def mismatch(Z, l, eps, X=14.0):
    kap = 2*l + 1
    def rhs(sig, y):
        cc = -np.tanh(sig)
        return [y[1], (kap**2 - 2*Z*(1 + eps*gfun(cc))/np.cosh(sig)**2)*y[0]]
    yl = solve_ivp(rhs, (-X, 0), [1.0, kap], rtol=1e-12, atol=1e-14).y[:, -1]
    yr = solve_ivp(rhs, (X, 0), [1.0, -kap], rtol=1e-12, atol=1e-14).y[:, -1]
    return yl[1]*yr[0] - yr[1]*yl[0]          # Wronskian at sigma = 0

def Zroot(N, l, eps):
    Z0 = N*(N + 1)/2
    return brentq(lambda Z: mismatch(Z, l, eps), Z0 - 0.3, Z0 + 0.3, xtol=1e-14)

def Eq(N, l):          # E_(q_2l)[g] with q_j ~ s^(j+1) C^(j+3/2)_(k-j)^2
    k = N - 1; jj = 2*l
    dens = lambda x: (1 - x*x)**(jj + 1)*eval_gegenbauer(k - jj, jj + 1.5, x)**2
    Zn = quad(dens, -1, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
    return quad(lambda x: gfun(x)*dens(x), -1, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0]/Zn

worst = 0.0; rows = []
eps = 1e-5
for N in (3, 4, 5, 6):
    for l in range(0, (N - 1)//2 + 1):
        if N - l < l + 1: continue
        dZ = (Zroot(N, l, eps) - Zroot(N, l, -eps))/(2*eps)
        pred = -(N*(N + 1)/2)*Eq(N, l)
        rows.append((N, l, dZ, pred)); worst = max(worst, abs(dZ - pred)/abs(pred))
for N, l, dZ, pred in rows:
    print(f"  N={N} l={l}: dZ/deps shooting {dZ:.7f}   -Z_N E_q[g] {pred:.7f}")
report(f"first order: dZ_(N,l)/deps = -Z_N E_(q_2l)[g] (R^3 shooting, no lift), worst rel. error {worst:.1e}",
       worst < 1e-5)
# order within each group: higher l has smaller dZ (binds first) for this convex weakening
order_ok = all(rows[i][2] > rows[i + 1][2] for i in range(len(rows) - 1) if rows[i][0] == rows[i + 1][0])
report("within each group dZ decreases with l (higher l = lower n needs less coupling)", order_ok)

print("ALL OK" if ok else "SOME CHECK FAILED")

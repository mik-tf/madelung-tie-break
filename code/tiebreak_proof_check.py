"""Checks every identity used in the paper (exact sympy arithmetic, no floating point).

Notation. s = 1 - x^2. For lambda > 1/2 and m >= 1, f = C^(lambda)_m (Gegenbauer), w = s^(lambda-1/2),
mu = m(m + 2 lambda). Block j of degree k (j = 0..k): q_j proportional to s^(j+1) [C^(j+3/2)_(k-j)]^2 on [-1,1],
normalised to mass 1. Phase-free block l is j = 2l. U_q(t) = E_q|c - t|.

(1) General identity (symbolic lambda, mu; uses only the Gegenbauer ODE):
    V = s^(lambda+3/2) [f'^2 - ((2 lambda - 1)/(2 lambda + 1)) f f'']  satisfies
    V'' = alpha * w * (f^2 - s f'^2 / mu),  alpha = mu (4 lambda^2 + 4 mu - 1)/(2 lambda + 1).
    With lambda = j + 3/2, m = k - j: alpha = 2 mu (k+1)(k+2)/(j+2).
    Lemma B: (s^(a+1) [f'^2 - ((a-1)/a) f f''])' = -(2/a)(mu + a(a-1)) s^a f f', a = lambda + 1/2.
(2) For every k <= 29 (N = k + 1 <= 30) and j < k, exactly:
    Lemma A: F_j - F_(j+1) = -(1/(mu nu)) s^(j+2) f f' (CDFs of c on [-1,1]);
    f satisfies the ODE; f' = 2 lambda C^(lambda+1)_(m-1); int s^(lambda+1/2) f'^2 = mu int w f^2;
    W_j := 2 V / (alpha nu) (nu = int w f^2) has W_j'' = 2 (q_j - q_(j+1)) and W_j(-1) = W_j(1) = 0,
    hence W_j = U_(q_j) - U_(q_(j+1)); Q_j = f'^2 - ((j+1)/(j+2)) f f'' has no real root and Q_j(0) > 0.
(3) For every N <= 30 and every consecutive phase-free pair l, l+1: H_l(t) := int_t^1 (F_(l+1) - F_l) ds,
    computed from the CDF definition of tiebreak_dominance.py, equals W_(2l)(t) + W_(2l+1)(t) identically.
(4) f'^2 - f f'' = sum_i (f/(x - x_i))^2 for f = prod (x - x_i), symbolic roots, degree 1..6.
Run: nice -n 19 env OMP_NUM_THREADS=1 python3 tiebreak_proof_check.py   (prints ALL OK at the end)
"""
import sys
from sympy import (symbols, Rational, Poly, QQ, gegenbauer, factor, simplify, expand, diff, prod,
                   cancel, integrate)

x, t, lam, mu, F0, F1 = symbols('x t lam mu F0 F1')
ok = True
def check(cond, msg):
    global ok
    if not cond:
        ok = False
        print("FAIL:", msg)

# ---- (1) symbolic identity
s = 1 - x**2
F2 = ((2*lam + 1)*x*F1 - mu*F0)/s          # f'' from the ODE  s f'' - (2 lam + 1) x f' + mu f = 0
Dx = lambda e: diff(e, x) + F1*diff(e, F0) + F2*diff(e, F1)
cst = (2*lam - 1)/(2*lam + 1)
R = s**2*(F1**2 - cst*F0*F2)                 # V = s^(lam - 1/2) R
dR = lambda e: Dx(e) + (lam - Rational(1, 2))*(-2*x)/s*e
R2 = dR(dR(R))                               # V'' = s^(lam - 1/2) R2
alpha = mu*(4*lam**2 + 4*mu - 1)/(2*lam + 1)
check(simplify(R2 - alpha*(F0**2 - s*F1**2/mu)) == 0, "(1) V'' identity")
j, k = symbols('j k')
m_ = k - j; lam_ = j + Rational(3, 2); mu_ = m_*(m_ + 2*lam_)
check(simplify(alpha.subs({lam: lam_, mu: mu_}) - 2*mu_*(k + 1)*(k + 2)/(j + 2)) == 0, "(1) alpha in j,k")
# Lemma B (first derivative, by hand in the md): with a = lam + 1/2,
# (s^(a+1) [f'^2 - ((a-1)/a) f f''])' = -(2/a) (mu + a(a-1)) s^a f f', and mu + a(a-1) = (k+1)(k+2).
A = lam + Rational(1, 2)
Rb = s*(F1**2 - (A - 1)/A*F0*F2)          # V = s^A Rb, so V' = s^A (Dx Rb + A (-2x)/s Rb)
dVr = Dx(Rb) + A*(-2*x)/s*Rb
check(simplify(dVr + 2/A*(mu + A*(A - 1))*F0*F1) == 0, "(1) Lemma B first derivative")
check(expand(mu_ + (lam_ + Rational(1, 2))*(lam_ - Rational(1, 2)) - (k + 1)*(k + 2)) == 0, "(1) mu + a(a-1)")
print("(1) Lemma B (V' = -(2/a)(k+1)(k+2) s^a f f') verified")
print("(1) symbolic identity V'' = alpha w (f^2 - s f'^2/mu) verified; alpha = 2 mu (k+1)(k+2)/(j+2)")

# ---- (2), (3) exact, k <= 29
P = lambda e: Poly(e, x, domain=QQ)
S = P(1 - x**2)
def geg(n, a):
    return P(gegenbauer(n, a, x)) if n > 0 else P(1)
def integ01(p, a, b):
    I = p.integrate(); return I.eval(b) - I.eval(a)
KMAX = 29
for k in range(1, KMAX + 1):
    q = []
    for jj in range(k + 1):
        a = Rational(2*jj + 3, 2)
        d = S**(jj + 1)*geg(k - jj, a)**2
        q.append(d*(1/integ01(d, -1, 1)))
    W = []
    for jj in range(k):
        a = Rational(2*jj + 3, 2); m = k - jj; muv = m*(m + 2*a)
        f = geg(m, a); f1 = f.diff(x); f2 = f1.diff(x)
        check((S*f2 - P((2*a + 1)*x)*f1 + f*muv).is_zero, f"ODE k={k} j={jj}")
        check((f1 - geg(m - 1, a + 1)*(2*a)).is_zero, f"derivative k={k} j={jj}")
        nu = integ01(S**(jj + 1)*f**2, -1, 1)
        check(integ01(S**(jj + 2)*f1**2, -1, 1) == muv*nu, f"normalisation k={k} j={jj}")
        Q = f1**2 - f*f2*Rational(jj + 1, jj + 2)
        check(Q.count_roots() == 0 and Q.eval(0) > 0, f"Q positive k={k} j={jj}")
        # Lemma A: CDF of c, F_j - F_(j+1) = -(1/(mu nu)) s^(j+2) f f'
        Fj = q[jj].integrate(); Fj = Fj - Fj.eval(-1)
        Fj1 = q[jj + 1].integrate(); Fj1 = Fj1 - Fj1.eval(-1)
        check((Fj - Fj1 + S**(jj + 2)*f*f1*(1/(muv*nu))).is_zero, f"Lemma A k={k} j={jj}")
        alv = 2*muv*(k + 1)*(k + 2)/Rational(jj + 2)
        Wj = S**(jj + 3)*Q*(2/(alv*nu))
        check((Wj.diff(x).diff(x) - (q[jj] - q[jj + 1])*2).is_zero, f"W'' k={k} j={jj}")
        check(Wj.eval(1) == 0 and Wj.eval(-1) == 0, f"W(+-1) k={k} j={jj}")
        W.append(Wj)
    # (3) CDF route, phase-free pairs, t in [0,1]
    for l in range(k//2):
        def cdf(p):                       # F(t) = P(|c| <= t) = int_0^t 2 p
            I = (p*2).integrate(); return I - I.eval(0)
        Dl = cdf(q[2*l + 2]) - cdf(q[2*l])
        G = Dl.integrate()
        H = P(G.eval(1)) - G             # int_t^1 (F_(l+1) - F_l)
        check((H - W[2*l] - W[2*l + 1]).is_zero, f"H = W+W N={k+1} l={l}")
    print(f"N={k+1}: identities checked for {k} blocks, {k//2} phase-free pairs")

# ---- (4) sum-of-squares identity
for n in range(1, 7):
    r = symbols(f'r0:{n}')
    f = prod([x - ri for ri in r])
    lhs = diff(f, x)**2 - f*diff(f, x, 2)
    rhs = sum((cancel(f/(x - ri)))**2 for ri in r)
    check(expand(lhs - rhs) == 0, f"SOS n={n}")
print("(4) f'^2 - f f'' = sum (f/(x-x_i))^2 verified for degree 1..6")
print("ALL OK" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)

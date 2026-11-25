"""tie-break dominance test, exact arithmetic.

Definitions (fixed before the run).
Group N = n + l = k + 1, k = 0..15 (N <= 16). Phase-free (circle-invariant) harmonics of
degree k on the round S^4 have SO(3) content l = 0..floor(k/2), once each. With
c = cos(theta) (pole = nucleus/infinity axis, Kelvin map c -> -c) and s^2 = 1 - c^2, the
harmonic is C^{(2l+3/2)}_{k-2l}(c) s^{2l} Y(circle-invariant, SO(3) type l), and its S^4
volume density in c is proportional to (1-c^2)^(2l+1) [C^{(2l+3/2)}_{k-2l}(c)]^2.
The density is even in c; |c| has density twice that on [0,1]. F_l(t) = P(|c| <= t).
Test, for consecutive l, l+1 in the same group: first-order stochastic dominance of
|c| (l+1 more concentrated toward the equator): F_{l+1}(t) >= F_l(t) for all t in [0,1].
Exact: Gegenbauer polynomials and CDFs are rational polynomials; nonnegativity of the
difference is decided by exact root counting (no floating point).
PASS: every pair dominates. FAIL: print first violating pair.
"""
import sys
from sympy import symbols, gegenbauer, Poly, integrate, Rational, QQ, sqf_list

c, t = symbols('c t')

def cdf(k, l):
    lam = Rational(4 * l + 3, 2)
    g = Poly(gegenbauer(k - 2 * l, lam, c), c, domain=QQ)
    dens = Poly((1 - c**2) ** (2 * l + 1), c, domain=QQ) * g * g
    cum = dens.integrate()
    tot = cum.eval(1) - cum.eval(0)
    p = (cum - cum.eval(0)) * 2 / (2 * tot)  # |c| on [0,1]: factor 2 cancels with total 2*tot
    return Poly(p.as_expr().subs(c, t), t, domain=QQ)

def nonneg_on_01(D):
    """True iff polynomial D(t) >= 0 on [0,1]; exact."""
    if D.is_zero:
        return True
    for _ in range(2):  # strip factors t and (t-1) (endpoint zeros)
        pass
    q = D
    while q.eval(0) == 0:
        q = q.quo(Poly(t, t, domain=QQ))
    while q.eval(1) == 0:
        q = q.quo(Poly(t - 1, t, domain=QQ))
    _, facs = sqf_list(q.as_expr(), t)
    for f, m in facs:
        if m % 2 == 1 and Poly(f, t, domain=QQ).count_roots(0, 1) > 0:
            return False  # sign change inside (0,1) (endpoints already stripped)
    # also sign of q at 1/2 (no sign change inside, so constant on (0,1)); D = t^a (t-1)^b q
    return D.eval(Rational(1, 2)) >= 0

def summary():
    """Extra record after a FAIL: all violating pairs, and whether mean |c| decreases in l."""
    bad, nonmono, total = [], [], 0
    for N in range(2, 17):
        k = N - 1
        prev = pm = None
        for l in range(0, k // 2 + 1):
            F = cdf(k, l)
            m = 1 - integrate(F.as_expr(), (t, 0, 1))
            if prev is not None:
                total += 1
                if not nonneg_on_01(F - prev):
                    bad.append((N, l - 1, l))
                if not m < pm:
                    nonmono.append((N, l - 1, l))
            prev, pm = F, m
    print(f"summary: {len(bad)} of {total} consecutive pairs violate dominance; "
          f"{len(nonmono)} violate the weaker order mean|c| decreasing in l")
    print("violating pairs (N,l,l+1):", bad)
    # exact sign-change count of F_(l+1) - F_l on (0,1) (odd-multiplicity roots), N <= 12: expected k - 2l - 2
    mism = []
    for N in range(3, 13):
        k = N - 1
        for l in range(0, k // 2):
            q = cdf(k, l + 1) - cdf(k, l)
            while q.eval(0) == 0:
                q = q.quo(Poly(t, t, domain=QQ))
            while q.eval(1) == 0:
                q = q.quo(Poly(t - 1, t, domain=QQ))
            _, fs = sqf_list(q.as_expr(), t)
            n = sum(Poly(g, t, domain=QQ).count_roots(0, 1) for g, m in fs if m % 2 == 1)
            if n != k - 2 * l - 2:
                mism.append((N, l, n))
    print("sign-change count k-2l-2 for N <= 12:", "OK for every pair" if not mism else f"MISMATCH {mism}")

def main():
    for N in range(2, 17):
        k = N - 1
        prev = None
        for l in range(0, k // 2 + 1):
            F = cdf(k, l)
            if prev is not None:
                if not nonneg_on_01(F - prev):
                    print(f"FAIL: N={N} (k={k}), pair l={l-1} -> l={l}")
                    summary()
                    return 1
            prev = F
        print(f"N={N}: ok ({k//2+1} harmonics)")
    print("PASS: dominance for every consecutive pair, N <= 16")
    return 0

sys.exit(main())

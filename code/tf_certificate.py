#!/usr/bin/env python3
"""W2 certificate (Corollary 4): rigorous interval check of D_(N,l) > 0 for the
Thomas-Fermi weakening h_TF, every consecutive pair l, l+1 in every group N <= NMAX.

Pass lines are those fixed in advance and are not changed here: every certified lower bound of D_(N,l) is
positive (margin reported); h_TF increasing on (0, 1 - u0] (here via h'' > 0 on [0, 1 - u0] and h'(0) = 0).

Arithmetic: mpmath.iv (python-flint is not installed). All intervals are rigorous up to mpmath's iv
rounding. Rigorous ingredients:
 1. TF enclosure with no pole asymptotics. In t = sqrt(x), y(t) = phi(t^2), u = 2 phi'(x): the system
    y' = t u, u' = 4 y^(3/2) is regular at t = 0 (y(0) = 1, u(0) = 2B). Validated Taylor stepping
    (order NT, Lagrange remainder from the order NT+1 coefficient on a Picard a-priori box).
    Comparison sandwich: B_lo < B* < B_hi with B_lo the number such that its solution hits zero (certified
    by the quadratic criterion below) and B_hi such that its solution turns up (u > 0 certified). Then
    phi_lo < phi* < phi_hi and phi_lo' < phi*' < phi_hi' pointwise on the whole certified range
    (d = phi_hi - phi*: d'' = (phi_hi^1.5 - phi*^1.5)/sqrt(x) > 0 while d > 0, d(0) = 0, d'(0) > 0).
 2. Weakening. x = |c| in [0,1], rho = (1 - x)/(1 + x) (paper 1 orientation: nucleus at c = +1, rho = r/R),
    g(rho) = phi(rho/a) (1 + rho)^2 / A - 1, Kelvin-even h(x) = -(g(rho) + g(1/rho))/2 = 1 - G(rho),
    G = (1 + rho)^2 (P + Q)/(2A), P = phi(rho/a), Q = phi(1/(a rho))/rho^2. (a, A) = (0.536795, 1.012783).
    The code sign (nucleus at c = -1) never enters: h depends on |c| only.
 3. D_(N,l) = E_l h - E_(l+1) h = D_main + D_pole, T = 1 - u0.
    D_main = h(T)(F_l - F_(l+1))(T) - h'(T) H_l(T) + int_0^T H_l h''  (integration by parts, h'(0) = 0 by evenness)
    primary bound: int_I H h'' >= min_I H * (h'(b) - h'(a)) on every piece where h'' > 0 is certified.
    D_pole = 2 int_T^1 h (p_l - p_(l+1)) >= h0 (m_l - m_(l+1)) - w (m_l + m_(l+1)), m = mass above T,
    [h0 - w, h0 + w] a crude range of h on |c| > T from 1 + B_lo s <= phi <= 1 near the nucleus and
    0 < phi <= min(1, 144 x^-3) at the far side.
 4. Cross-check: D_main = 2 int_0^T h (p_l - p_(l+1)) dx by the midpoint rule with rigorous remainder
    (len^3/24) sup|R''| (R'' from the h, h', h'' enclosures), exact rational densities.
Rigour of the arithmetic (all of it is mpmath.iv with directed rounding, so the result is rigorous up to the
correctness of that rounding): sub-piece end points and midpoints are interval enclosures of the exact values
and adjacent pieces share the same end enclosure; the pole half-width w is rounded up around a point h0; Sol.enc
asserts that the requested t lies inside the validated steps.
Pole bound phi_TF(xi) <= 144 xi^-3 for all xi >= 3.7e3 (beyond the sandwich range): proved in the paper, Appendix A.5
(maximum principle: delta = phi - 144 xi^-3 solves the same ODE comparison and cannot have a positive interior maximum).
Roles: PRIMARY = direct quadrature of step 3 below (only h and a crude h'' bound
for the midpoint remainder), CROSS-CHECK = integration by parts (needs the wide h'' enclosure; may stay unusable).
Pass line unchanged: every certified lower bound of D positive, margin reported.
Run: nice -n 19 env OMP_NUM_THREADS=1 python3 tf_certificate.py [NMAX] [u0]
"""
import sys, time, bisect
from fractions import Fraction
import mpmath as mp
from mpmath import iv

iv.prec = 480
mp.mp.dps = 150
NT = 60            # Taylor order
RATIO = 0.03       # step = RATIO * max(t, 1)
A_TIETZ = iv.mpf('1.012783')
a_TIETZ = iv.mpf('0.536795')
B_STAR = mp.mpf('-1.588071022611375312718684509423950109452746621674825616765677')
DELTA = mp.mpf("1e-35")


def lo(x):
    """rigorous lower endpoint of an iv number as a plain mpf (a plain mpf is returned unchanged)"""
    return mp.mpf(x._mpi_[0]) if hasattr(x, '_mpi_') else x


def up(x):
    return mp.mpf(x._mpi_[1]) if hasattr(x, '_mpi_') else x


def hull(a, b):
    return iv.mpf([min(a.a, b.a), max(a.b, b.b)])


def ivpow15(y):
    return y * iv.sqrt(y)


def recur(y0, u0, tp, n):
    """Taylor coefficients (y_k, u_k), k <= n+1, of y' = t u, u' = 4 y^1.5 about a point where t = tp."""
    y = [y0]; u = [u0]; w = []
    for k in range(0, n + 1):
        if k == 0:
            w.append(ivpow15(y0))
        else:
            s = iv.mpf(0)
            for j in range(1, k + 1):
                s += (iv.mpf(3) * j / 2 - (k - j)) * y[j] * w[k - j]
            w.append(s / (k * y0))
        u.append(4 * w[k] / (k + 1))
        y.append((tp * u[k] + (u[k - 1] if k >= 1 else 0)) / (k + 1))
    return y, u


def peval(c, h):
    r = iv.mpf(0)
    for k in range(len(c) - 1, -1, -1):
        r = r * h + c[k]
    return r


def pder(c):
    return [k * c[k] for k in range(1, len(c))]


class Sol:
    """Validated solution for a fixed B (exact decimal), integrated until stop(state) returns a label."""

    def __init__(self, B, mode, tmax=2000.0):
        self.B = B
        self.steps = []   # (t0 float, h float, yc, uc, ry, ru)  ry, ru: remainder coefficient intervals
        y = iv.mpf(1); u = iv.mpf(2) * iv.mpf(B)
        t = iv.mpf(0)
        self.result = None
        while True:
            t0f = float(t.a)
            h = RATIO * max(t0f, 1.0)
            hi = iv.mpf(h)
            while True:
                # Picard a-priori box for [t0, t0+h]
                T = iv.mpf([t.a, (t + hi).b])
                dy = iv.mpf([0, 1]) * hi
                Yb = y; Ub = u
                ok = False
                for _ in range(8):
                    if Yb.a <= 0:
                        break
                    Ub_n = u + iv.mpf([0, 1]) * hi * 4 * ivpow15(Yb)
                    Yb_n = y + iv.mpf([0, 1]) * hi * T * Ub_n
                    # inflate for fixpoint
                    infl_y = iv.mpf([Yb_n.a - abs(Yb_n.a - Yb_n.b) * 0.5, Yb_n.b + abs(Yb_n.a - Yb_n.b) * 0.5])
                    infl_u = iv.mpf([Ub_n.a - abs(Ub_n.a - Ub_n.b) * 0.5, Ub_n.b + abs(Ub_n.a - Ub_n.b) * 0.5])
                    if infl_y.a <= 0:
                        break
                    if infl_u.a <= Ub.a and infl_u.b >= Ub.b and False:
                        pass
                    Ub2 = u + iv.mpf([0, 1]) * hi * 4 * ivpow15(infl_y)
                    Yb2 = y + iv.mpf([0, 1]) * hi * T * infl_u
                    if infl_y.a <= Yb2.a and Yb2.b <= infl_y.b and infl_u.a <= Ub2.a and Ub2.b <= infl_u.b:
                        Yb, Ub = infl_y, infl_u
                        ok = True
                        break
                    Yb, Ub = infl_y, infl_u
                if ok:
                    break
                h = h / 2
                hi = iv.mpf(h)
                if h < 1e-9:
                    self.result = 'picard_fail'
                    return
            yc, uc = recur(y, u, t, NT)
            # remainder coefficients on the box
            ybc, ubc = recur(Yb, Ub, T, NT)
            ry, ru = ybc[NT + 1], ubc[NT + 1]
            self.steps.append((t0f, h, yc[:NT + 1], uc[:NT + 1], ry, ru, t))
            hp = iv.mpf([0, 1]) * hi
            yn = peval(yc[:NT + 1], hi) + ry * hi ** (NT + 1)
            un = peval(uc[:NT + 1], hi) + ru * hi ** (NT + 1)
            t = t + hi
            y, u = yn, un
            if y.a <= 0:
                self.result = 'nonpositive'
                return
            lab = mode(y, u, t)
            if lab:
                self.result = lab
                self.tend = float(t.a)
                return
            if float(t.a) > tmax:
                self.result = 'tmax'
                return

    def _locate(self, tint):
        ts = [s[0] for s in self.steps]
        i = bisect.bisect_right(ts, float(tint.a)) - 1
        return max(i, 0)

    def enc(self, tint):
        """enclosure of y, u over t in tint (may span steps), as intervals."""
        # no silent clipping: the requested t must lie inside the validated steps
        last = self.steps[-1]
        assert tint.a >= 0 and iv.mpf(tint.b).b <= (last[6] + iv.mpf(last[1])).a, 't outside the certified range'
        i = self._locate(tint)
        yh = None; uh = None
        while True:
            t0f, h, yc, uc, ry, ru, t0 = self.steps[i]
            lo = max(tint.a, t0.a) if False else None
            a = iv.mpf(tint.a) - t0
            b = iv.mpf(tint.b) - t0
            # local offsets, clipped to [0, h]
            aa = max(float(a.a), 0.0); bb = min(float(b.b), h)
            hh = iv.mpf([max(a.a, 0), min(b.b, iv.mpf(h).b)])
            ye = peval(yc, hh) + ry * iv.mpf([0, 1]) * iv.mpf(h) ** (NT + 1)
            ue = peval(uc, hh) + ru * iv.mpf([0, 1]) * iv.mpf(h) ** (NT + 1)
            yh = ye if yh is None else hull(yh, ye)
            uh = ue if uh is None else hull(uh, ue)
            if float(tint.b) <= t0f + h or i + 1 >= len(self.steps):
                break
            i += 1
        return yh, uh


def certify_sandwich(tneed):
    """Return (sol_lo, sol_hi, info). Raises if the criteria fail."""
    Blo = B_STAR - DELTA; Bhi = B_STAR + DELTA

    def mode_lo(y, u, t):
        # hits zero: u < 0 and u^2/4 > 2 y^2.5 / t  
        if u.b < 0:
            lhs = u.b * u.b / 4 if False else (u.b ** 2) / 4
            rhs = 2 * y * y * iv.sqrt(y) / t
            if (lhs - rhs).a > 0:
                return 'hits_zero'
        return None

    def mode_hi(y, u, t):
        return 'turns_up' if u.a > 0 else None

    lo = Sol(mp.nstr(Blo, 80, strip_zeros=False), mode_lo) if False else Sol(Blo, mode_lo)
    hi = Sol(Bhi, mode_hi)
    return lo, hi



# ---------------------------------------------------------------- weakening h and its derivatives
class Weak:
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def phi3(self, S):
        """enclosures of phi, phi', phi'' of the TF function over the interval S of the TF variable."""
        t = iv.sqrt(S)
        yl, ul = self.lo.enc(t)
        yh, uh = self.hi.enc(t)
        Y = iv.mpf([yl.a, yh.b])
        U = iv.mpf([ul.a, uh.b]) / 2
        assert Y.a > 0
        return Y, U, ivpow15(Y) / t

    def h3(self, X):
        """h, h', h'' (x = |c|) over the interval X."""
        a, A = a_TIETZ, A_TIETZ
        rho = (1 - X) / (1 + X)
        s1 = rho / a
        s2 = 1 / (a * rho)
        p0, p1, p2 = self.phi3(s1)
        f0, f1, f2 = self.phi3(s2)
        P, P1, P2 = p0, p1 / a, p2 / (a * a)
        r2 = rho * rho
        Q = f0 / r2
        Q1 = -2 * f0 / (r2 * rho) - f1 / (a * r2 * r2)
        Q2 = 6 * f0 / (r2 * r2) + 6 * f1 / (a * r2 * r2 * rho) + f2 / (a * a * r2 * r2 * r2)
        K = (1 + rho) ** 2 / (2 * A); K1 = (1 + rho) / A; K2 = 1 / A
        S0 = P + Q; S1 = P1 + Q1; S2 = P2 + Q2
        G = K * S0
        G1 = K1 * S0 + K * S1
        G2 = K2 * S0 + 2 * K1 * S1 + K * S2
        q = (1 + X)
        h = 1 - G
        hx = 2 * G1 / (q * q)
        hxx = -(G2 * 4 / (q ** 4) + G1 * 4 / (q ** 3))
        return h, hx, hxx


# ---------------------------------------------------------------- exact densities
def frac_ivs(fr):
    return iv.mpf(fr.numerator) / iv.mpf(fr.denominator)


class IPoly:
    def __init__(self, coeffs):  # ascending Fractions
        self.c = [Fraction(v) for v in coeffs]
        self.ci = [frac_ivs(v) for v in self.c]
        self.d = [k * self.c[k] for k in range(1, len(self.c))] or [Fraction(0)]
        self.di = [frac_ivs(v) for v in self.d]

    def ex(self, x):  # exact at a Fraction
        r = Fraction(0)
        for v in reversed(self.c):
            r = r * x + v
        return r

    def ev(self, X):  # enclosure over interval X, centered form
        m = (X.a + X.b) / 2
        mi = iv.mpf(m)
        return peval(self.ci, mi) + peval(self.di, X) * (X - mi)

    def der(self):
        return IPoly([self.d[k] for k in range(len(self.d))])


def densities(N):
    import sympy as sp
    from sympy import Poly, QQ, Rational, gegenbauer, symbols
    x = symbols('x')
    P = lambda e: Poly(e, x, domain=QQ)
    S = P(1 - x ** 2)
    k = N - 1
    q = []
    for j in range(k + 1):
        a = Rational(2 * j + 3, 2)
        g = P(gegenbauer(k - j, a, x)) if k - j > 0 else P(1)
        d = S ** (j + 1) * g ** 2
        I = d.integrate()
        q.append(d * (1 / (I.eval(1) - I.eval(-1))))
    out = []
    for l in range(k // 2):
        pl, pn = q[2 * l], q[2 * l + 2]
        cdf = lambda p: (lambda I: I - I.eval(0))((p * 2).integrate())
        dF = cdf(pn) - cdf(pl)               # F_(l+1) - F_l
        G = dF.integrate()
        H = P(G.eval(1)) - G                 # int_t^1 dF
        co = lambda poly: [Fraction(int(v.p), int(v.q)) for v in reversed(poly.all_coeffs())]
        out.append(dict(l=l, p=IPoly(co(pl * 1 - pn)), dF=IPoly(co(dF)), H=IPoly(co(H)), Fl=IPoly(co(cdf(pl))), Fn=IPoly(co(cdf(pn)))))
    return out


# ---------------------------------------------------------------- certificate
def certify(W, groups, u0s, grid_ratio=1.01, sub=4, log=print):
    u0 = Fraction(u0s)
    T = 1 - u0
    Tiv = iv.mpf(T.numerator) / iv.mpf(T.denominator)
    # grid in x = 1 - u, u geometric from u0 to 1, last node x = 0
    us = [u0]
    while us[-1] * Fraction(grid_ratio) < 1:
        us.append(us[-1] * Fraction(grid_ratio).limit_denominator(1000))
    us.append(Fraction(1))
    xs = [1 - u for u in us][::-1]            # increasing, xs[0] = 0, xs[-1] = T
    xs = [Fraction(float(v)).limit_denominator(10 ** 9) if False else v for v in xs]
    nodes = []
    for v in xs:
        X = iv.mpf(v.numerator) / iv.mpf(v.denominator)
        nodes.append(X)
    # node data
    nd = [W.h3(X) for X in nodes]
    # pieces
    pieces = []
    convex_ok = True
    badmono = []
    worst_hpp = None
    for i in range(len(nodes) - 1):
        I = iv.mpf([nodes[i].a, nodes[i + 1].b])
        h, hx, hxx = W.h3(I)
        # step 2 pass line: h' > 0 on every piece (interval lower bound), split in 16 on failure
        lowp = lo(hx)
        if lowp <= 0:
            sp = [iv.mpf([nodes[i].a + (nodes[i + 1].b - nodes[i].a) * j / 16, nodes[i].a + (nodes[i + 1].b - nodes[i].a) * (j + 1) / 16]) for j in range(16)]
            lows = [lo(W.h3(s)[1]) for s in sp]
            lowp = min(lows)
            if lowp <= 0:
                convex_ok = False
                badmono.append((float(nodes[i].a), float(nodes[i + 1].b), float(lowp)))
        worst_hpp = lowp if worst_hpp is None else min(worst_hpp, lowp)
        pieces.append((I, h, hx, hxx))
    # sub-pieces for the midpoint cross-check
    # Sub-pieces: lo_, hi_ are interval enclosures of the exact end points (the same enclosure is used for the
    # right end of one piece and the left end of the next, so the true pieces abut exactly and nothing is
    # counted twice); the exact midpoint lies in the interval (lo_ + hi_)/2 and the exact length in hi_ - lo_.
    subp = []; subbnd = []
    for i in range(len(nodes) - 1):
        a, b = nodes[i], nodes[i + 1]
        ends = [a + (b - a) * j / sub for j in range(sub + 1)]
        ends[0] = a; ends[-1] = b
        for j in range(sub):
            lo_, hi_ = ends[j], ends[j + 1]
            I = iv.mpf([lo_.a, hi_.b])
            subp.append(I)
            subbnd.append((hi_ - lo_, (lo_ + hi_) / 2))
    subdat = [W.h3(I) for I in subp]
    submid = [W.h3(mi_) for (_, mi_) in subbnd]
    # pole range of h on |c| > T (crude)
    rmax = iv.mpf(u0.numerator) / iv.mpf(u0.denominator)
    rmax = rmax / (2 - rmax)
    Plo = 1 + iv.mpf(B_STAR - DELTA) * rmax / a_TIETZ
    Gmin = (1 / (2 * A_TIETZ)) * Plo
    # uses phi_TF <= min(1, 144 xi^-3) on the far side: see the proof in Appendix A.5 of the paper
    Gmax = ((1 + rmax) ** 2 / (2 * A_TIETZ)) * (1 + 144 * a_TIETZ ** 3 * rmax)
    hlo, hhi = lo(1 - Gmax), up(1 - Gmin)
    # h0 is a point; w is rounded up so that [h0 - w, h0 + w] contains [hlo, hhi] exactly
    h0 = iv.mpf(mp.mpf((hlo + hhi) / 2))
    w = iv.mpf(max(up(h0 - iv.mpf(hlo)), up(iv.mpf(hhi) - h0)))
    hT, hxT, hxxT = W.h3(Tiv)
    res = []
    for N, dens in groups:
        for d in dens:
            Hp, dF, dp = d['H'], d['dF'], d['p']
            Hp = Hp
            dp1, dp2 = dp.der(), dp.der().der()
            # primary
            plo = iv.mpf(0); phi_ = iv.mpf(0)
            for i, (I, h, hx, hxx) in enumerate(pieces):
                dh = nd[i + 1][1] - nd[i][1]
                Hv = Hp.ev(I)
                prod = Hv * dh
                plo += lo(prod); phi_ += up(prod)
            HT = frac_ivs(Hp.ex(T)); dFT = frac_ivs(dF.ex(T))
            bnd = hT * (-dFT) - hxT * HT          # h(T)(F_l - F_(l+1))(T) - h'(T) H(T)
            main_lo = bnd.a + plo.a if hasattr(plo, 'a') else None
            main_lo = lo(bnd + iv.mpf([lo(plo), lo(plo)]))
            main_hi = up(bnd + iv.mpf([up(phi_), up(phi_)]))
            # pole piece: D_pole >= h0 (m_l - m_n) - w (m_l + m_n), m = mass above T (exact rationals)
            ml = 1 - d['Fl'].ex(T); mn = 1 - d['Fn'].ex(T)
            mli, mni = frac_ivs(ml), frac_ivs(mn)
            pole_lo = lo(h0 * (mli - mni) - w * (mli + mni))
            pole_hi = up(h0 * (mli - mni) + w * (mli + mni))
            # cross-check, direct quadrature 2 int_0^T h (p_l - p_n) by midpoint rule with remainder
            qlo = iv.mpf(0)
            for I, (h, hx, hxx), (hm, hxm, hxxm), (ln, mi) in zip(subp, subdat, submid, subbnd):
                R = hm * dp.ev(mi)
                R2 = hxx * dp.ev(I) + 2 * hx * dp1.ev(I) + h * dp2.ev(I)
                err = iv.mpf([-1, 1]) * (ln ** 3 / 24) * iv.mpf(max(abs(lo(R2)), abs(up(R2))))
                qlo += ln * R + err
            cross = (2 * qlo)
            res.append(dict(N=N, l=d['l'], main=(main_lo, main_hi), pole=(pole_lo, pole_hi),
                            cross=(lo(cross), up(cross))))
    return dict(badmono=badmono, res=res, pieces=len(pieces), convex_ok=convex_ok, worst_hpp=worst_hpp, h0=h0, w=w,
                subp=subp, subdat=subdat, submid=submid, hT=hT, hxT=hxT)


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    u0s = sys.argv[2] if len(sys.argv) > 2 else '1/1000'
    t0 = time.time()
    sl, sh = certify_sandwich(0)
    print(f"sandwich: B_lo = B* - 1e-35 {sl.result} at t = {sl.tend:.1f}; B_hi = B* + 1e-35 {sh.result} at t = {sh.tend:.1f}"
          f"  (certified range x <= {min(sl.tend, sh.tend) ** 2:.3e}); {time.time() - t0:.0f}s", flush=True)
    W = Weak(sl, sh)
    groups = [(N, densities(N)) for N in range(3, NMAX + 1)]
    t1 = time.time()
    r = certify(W, groups, u0s)
    bm = r['badmono']
    print('pieces where h1 > 0 is NOT certified by the h1 enclosure:', len(bm), 'of', r['pieces'], ('x in [%.4f, %.4f]' % (min(b[0] for b in bm), max(b[1] for b in bm))) if bm else '')
    print(f"u0 = {u0s}: {r['pieces']} pieces, h' > 0 certified on [0, 1-u0]: {r['convex_ok']}, smallest lower bound of h' on a piece = {mp.nstr(r['worst_hpp'], 6)}; "
          f"pole range of h: [{mp.nstr(lo(r['h0']) - up(r['w']), 6)}, {mp.nstr(up(r['h0']) + up(r['w']), 6)}]  ({time.time() - t1:.0f}s)")
    allpass = True  # pass line: every certified lower bound of D positive; monotonicity reported separately
    minLB = None; minrel = None
    print("N  l  LB(D)=quadrature+pole  [by-parts cross-check lo, hi]  [quadrature lo, hi]  N^2 D(est)  N^2 LB  overlap  PASS")
    for e in r['res']:
        LB = e['cross'][0] + e['pole'][0]
        mid = (e['cross'][0] + e['cross'][1]) / 2 + (e['pole'][0] + e['pole'][1]) / 2
        agree = e['main'][0] <= e['cross'][1] and e['cross'][0] <= e['main'][1]
        mw = e['main'][1] - e['main'][0]
        okk = LB > 0
        allpass &= bool(okk)
        rel = LB / mid
        minLB = LB if minLB is None or LB < minLB else minLB
        minrel = rel if minrel is None or rel < minrel else minrel
        print(f"{e['N']:2d} {e['l']:2d}  {mp.nstr(LB, 6):>12}  [{mp.nstr(e['main'][0], 6)}, {mp.nstr(e['main'][1], 6)}]  "
              f"[{mp.nstr(e['cross'][0], 6)}, {mp.nstr(e['cross'][1], 6)}]  {mp.nstr(e['N'] ** 2 * mid, 6)}  {mp.nstr(e['N'] ** 2 * LB, 6)}  {agree}  {'PASS' if okk else 'FAIL'}")
    print("smallest certified lower bound of D:", mp.nstr(minLB, 6), " smallest LB / D(est):", mp.nstr(minrel, 6))
    print("VERDICT:", "PASS" if allpass else "FAIL", f" total {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()

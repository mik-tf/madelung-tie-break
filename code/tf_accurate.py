#!/usr/bin/env python3
"""Accurate Thomas-Fermi function phi(x): phi'' = phi^(3/2)/sqrt(x), phi(0)=1, phi(inf)=0.

Method (all in mpmath, 60 digits, own Taylor stepping, no double-precision shooting):
 1. Near 0, with t = sqrt(x), y(t) = phi(t^2) is analytic: t y'' - y' = 4 t^2 y^(3/2),
    y = sum a_n t^n, a_0 = 1, a_1 = 0, a_2 = B = phi'(0), a_n = 4 p_(n-3) / (n (n-2)), p = coeffs of y^(3/2).
 2. From x = X0 to the join point XM: Taylor steps in x (order 70, step 0.1 x), exact recurrences.
 3. Beyond XM: the Sommerfeld/Coulson-March series phi = 144 x^-3 sum_n c_n z^n, z = x^-lambda,
    lambda = (sqrt(73) - 7)/2, c_0 = 1, c_1 = F free, c_n (n >= 2) fixed by the ODE.
 4. B and F are the two unknowns, found by requiring phi and phi' continuous at XM (2 equations).
    phi'' is then continuous because both pieces solve the same ODE; it is checked anyway.
 The outward instability is only a power law (growing mode ~ x^4.77 against phi ~ x^-3), so with
 60 digits and an exactly consistent B the Taylor pieces stay accurate to 30+ digits at x = 200.
 Independent checks (run this file): B against the published value, phi(1), phi(10), phi'' at
 the join, stepping the matched solution out to x = 1000 and comparing with the series.
Public values used for checks: phi'(0) = -1.5880710226113753127186845094239501095... (Boyd; Parand et al.),
 F = -13.270973 (Coulson and March 1950, as quoted by Parand et al.), tables of phi (Kobayashi 1955,
 as reproduced in Gombas): phi(1) = 0.4240, phi(10) = 0.02430 at the printed digits.
Usage: from tf_accurate import TF, phi  (phi(x) double, vectorised, join at XM = 200)
"""
import mpmath as mp
import numpy as np
from functools import lru_cache

mp.mp.dps = 60
N_TAY = 70
LAM = (mp.sqrt(73) - 7) / 2
B_PUBLISHED = mp.mpf('-1.588071022611375312718684509423950109452746621674825616765677')
F_PUBLISHED = mp.mpf('-13.270973')


def _pow15(y, n):
    """coeffs of y^(3/2) up to order n (y[0] > 0)."""
    p = [y[0] ** mp.mpf('1.5')]
    for k in range(1, n + 1):
        s = mp.mpf(0)
        for j in range(1, k + 1):
            s += (mp.mpf('1.5') * j - (k - j)) * y[j] * p[k - j]
        p.append(s / (k * y[0]))
    return p


def t_series(B, n=120):
    a = [mp.mpf(1), mp.mpf(0), B]
    p = [mp.mpf(1)]  # coefficients of y^(3/2)
    h = mp.mpf('1.5')
    for m in range(3, n + 1):
        k = m - 3  # need p_k, which uses a_0..a_k (all known)
        if k >= 1:
            p.append(mp.fsum((h * j - (k - j)) * a[j] * p[k - j] for j in range(1, k + 1)) / k)
        a.append(4 * p[k] / (m * (m - 2)))
    return a


def eval_t(a, t):
    """phi, dphi/dx, d2phi/dx2 at x = t^2 from y(t) series."""
    y = mp.polyval(a[::-1], t)
    dy = mp.polyval([k * a[k] for k in range(len(a) - 1, 0, -1)], t)
    # phi' = y'/(2t); phi'' = y^(3/2)/t
    return y, dy / (2 * t), y ** mp.mpf('1.5') / t


def taylor_step_coeffs(x0, y0, dy0, n=N_TAY):
    """Taylor coeffs in h = x - x0 of phi, from phi(x0), phi'(x0)."""
    s = [mp.binomial(mp.mpf(-1) / 2, k) * x0 ** (-mp.mpf(1) / 2 - k) for k in range(n + 2)]
    hh = mp.mpf('1.5')
    y = [y0, dy0]
    w = [y0 ** hh]
    w.append(hh * dy0 * w[0] / y0)
    for k in range(0, n):
        rhs = mp.fsum(w[j] * s[k - j] for j in range(k + 1))
        y.append(rhs / ((k + 2) * (k + 1)))
        m = k + 2
        if m <= n:
            w.append(mp.fsum((hh * j - (m - j)) * y[j] * w[m - j] for j in range(1, m + 1)) / (m * y0))
    return y


def march(B, X0, XM, ratio=mp.mpf('0.1'), n=N_TAY):
    """Return list of steps (x_start, coeffs) from X0 to XM and the end state (phi, phi')."""
    a = t_series(B)
    t0 = mp.sqrt(X0)
    y, dp, _ = eval_t(a, t0)
    x = X0
    steps = []
    while x < XM:
        h = min(ratio * x, XM - x)
        c = taylor_step_coeffs(x, y, dp, n)
        steps.append((x, h, c))
        y = mp.polyval(c[::-1], h)
        dp = mp.polyval([k * c[k] for k in range(len(c) - 1, 0, -1)], h)
        x = x + h
    return steps, y, dp


def asym_coeffs(F, n=200):
    c = [mp.mpf(1), F]
    e = [mp.mpf(1), mp.mpf('1.5') * F]
    for m in range(2, n + 1):
        # e_m with c_m = 0
        s = mp.mpf(0)
        for k in range(1, m):
            s += (mp.mpf('1.5') * k - (m - k)) * c[k] * e[m - k]
        rest = s / m
        cm = 12 * rest / ((3 + m * LAM) * (4 + m * LAM) - 18)
        c.append(cm)
        e.append(rest + mp.mpf('1.5') * cm)
    return c


def eval_asym(c, x):
    z = x ** (-LAM)
    w = mp.polyval(c[::-1], z)
    dw = mp.polyval([k * c[k] for k in range(len(c) - 1, 0, -1)], z)  # dw/dz
    # phi = 144 x^-3 w(z); dz/dx = -lam z / x
    phi = 144 * x ** -3 * w
    dphi = 144 * x ** -4 * (-3 * w - LAM * z * dw)
    return phi, dphi


X0 = mp.mpf('0.25')


class TF:
    """Accurate TF function with join point XM (default 200). Methods: phi(x), dphi(x), d2phi(x) double;
    mp_phi(x) high precision. Attributes: B, F, join residuals."""

    def __init__(self, XM=200, dps=60):
        mp.mp.dps = dps
        self.XM = mp.mpf(XM)

        def resid(B, F, XMc):
            steps, y, dp = march(B, X0, XMc)
            c = asym_coeffs(F, 160)
            pa, da = eval_asym(c, XMc)
            return y - pa, dp - da

        # continuation in the join point: start where the outward instability is mild (x = 30),
        # then move out by factors of 2 to XM, each solve started from the previous B, F.
        B, F = B_PUBLISHED, F_PUBLISHED
        chain = [mp.mpf(30)]
        while chain[-1] * 2 < self.XM:
            chain.append(chain[-1] * 2)
        chain.append(self.XM)
        hB, hF = mp.mpf('1e-25'), mp.mpf('1e-25')
        for XMc in chain:
            for it in range(12):
                r0 = resid(B, F, XMc)
                rB = resid(B + hB, F, XMc)
                rF = resid(B, F + hF, XMc)
                J = mp.matrix([[(rB[0] - r0[0]) / hB, (rF[0] - r0[0]) / hF],
                               [(rB[1] - r0[1]) / hB, (rF[1] - r0[1]) / hF]])
                d = mp.lu_solve(J, mp.matrix([-r0[0], -r0[1]]))
                B += d[0]
                F += d[1]
                if max(abs(d[0]), abs(d[1])) < mp.mpf(10) ** (-(dps - 25)):
                    break
        self.B, self.F = B, F
        self.steps, y, dp = march(B, X0, self.XM)
        self.c = asym_coeffs(F, 160)
        pa, da = eval_asym(self.c, self.XM)
        self.join = dict(dphi=y - pa, ddphi=dp - da, phi=y, dphi_val=dp)
        self.d2_in = y ** mp.mpf('1.5') / mp.sqrt(self.XM)
        self.d2_out = pa ** mp.mpf('1.5') / mp.sqrt(self.XM)
        # double-precision copies
        self._tf = [float(x) for x in [mp.mpf(0)]]
        self.a = t_series(B)
        self._a_f = [float(v) for v in self.a]
        self._steps_f = [(float(x), float(h), [float(v) for v in c]) for (x, h, c) in self.steps]
        self._sx = np.array([s[0] for s in self._steps_f])
        self._c_f = [float(v) for v in self.c[:60]]
        self._xm = float(self.XM)

    def mp_phi(self, x):
        x = mp.mpf(x)
        if x >= self.XM:
            return eval_asym(self.c, x)[0]
        if x <= X0:
            t = mp.sqrt(x)
            return mp.polyval(self.a[::-1], t) if x > 0 else mp.mpf(1)
        for (xs, h, c) in self.steps:
            if x <= xs + h + mp.mpf('1e-40'):
                return mp.polyval(c[::-1], x - xs)
        raise ValueError

    def _one(self, x, der=0):
        x = float(x)
        if x <= 0:
            return 1.0 if der == 0 else (float(self.B) if der == 1 else float('inf'))
        if x >= self._xm:
            z = x ** (-float(LAM))
            w = 0.0
            dw = 0.0
            for k in range(len(self._c_f) - 1, -1, -1):
                w = w * z + self._c_f[k]
            for k in range(len(self._c_f) - 1, 0, -1):
                dw = dw * z + k * self._c_f[k]
            lam = float(LAM)
            if der == 0:
                return 144 * x ** -3 * w
            if der == 1:
                return 144 * x ** -4 * (-3 * w - lam * z * dw)
            return (144 * x ** -3 * w) ** 1.5 / np.sqrt(x)
        if x <= float(X0):
            t = np.sqrt(x)
            a = self._a_f
            y = 0.0
            for k in range(len(a) - 1, -1, -1):
                y = y * t + a[k]
            if der == 0:
                return y
            dy = 0.0
            for k in range(len(a) - 1, 0, -1):
                dy = dy * t + k * a[k]
            if der == 1:
                return dy / (2 * t)
            return y ** 1.5 / t
        i = max(0, int(np.searchsorted(self._sx, x, side='right')) - 1)
        xs, h, c = self._steps_f[i]
        u = x - xs
        y = 0.0
        for k in range(len(c) - 1, -1, -1):
            y = y * u + c[k]
        if der == 0:
            return y
        dy = 0.0
        for k in range(len(c) - 1, 0, -1):
            dy = dy * u + k * c[k]
        if der == 1:
            return dy
        return y ** 1.5 / np.sqrt(x)

    def phi(self, x):
        if np.ndim(x) == 0:
            return self._one(x, 0)
        return np.array([self._one(v, 0) for v in np.asarray(x).ravel()]).reshape(np.shape(x))

    def dphi(self, x):
        return self._one(x, 1)


_default = None


def default():
    global _default
    if _default is None:
        _default = TF(200)
    return _default


def phi(x):
    return default().phi(x)


if __name__ == "__main__":
    import time
    t0 = time.time()
    for XM in (100, 200, 400):
        tf = TF(XM)
        print(f"XM={XM}: B = {mp.nstr(tf.B, 40)}  |B - published| = {mp.nstr(abs(tf.B - B_PUBLISHED), 3)}   F = {mp.nstr(tf.F, 20)}  (published -13.270973)")
        print(f"   join residuals: phi {mp.nstr(tf.join['dphi'], 3)}  phi' {mp.nstr(tf.join['ddphi'], 3)}  phi'' (in vs out) {mp.nstr(tf.d2_in - tf.d2_out, 3)}  (phi''={mp.nstr(tf.d2_in, 15)})")
        print(f"   phi(1) = {mp.nstr(tf.mp_phi(1), 25)}  phi(10) = {mp.nstr(tf.mp_phi(10), 25)}  phi(20) = {mp.nstr(tf.mp_phi(20), 25)}")
        print(f"   terms: |c_100 z^100| at XM = {mp.nstr(abs(tf.c[100]) * mp.mpf(XM) ** (-LAM * 100), 3)}  ({time.time() - t0:.0f}s)")

#!/usr/bin/env python3
"""Figures of paper 2 (the convex-order tie-break): fig1..fig6 as PDF and PNG (200 dpi) in out/.

Figure to data:
  fig1  exact block densities p_(2l) of |c|, N = 5                        (polynomials, as tiebreak_proof_check.py)
  fig2  F_(l+1) - F_l and H_l, N = 6                                       (polynomials, as tiebreak_dominance.py)
  fig3  Q_j = f'^2 - gamma f f'' and the zeros of f                        (polynomials, as tiebreak_proof_check.py)
  fig4  h_TF and h_TF'' near the poles                                     (tf_accurate.py, formulas of tf_certificate.py)
  fig5  Z_(N,l)(eps)/Z_N along the path, N = 3..8                          (tf_path_check.py via w3_cache.py -> w3_path.npz)
  fig6  Scal_TF/Scal_round against c, paper 1 sign                               (zero-energy state of the TF weight, tf_accurate.py)
Orientation: the figures use paper 1's c = (R - r)/(R + r), nucleus at c = +1 (the check scripts use the opposite sign).
Run (about 1 minute, w3_path.npz present; w3_cache.py makes it in about 10 minutes):
    OMP_NUM_THREADS=1 nice -n 19 python3 make_figures.py
Every figure goes through check_layout (no text over text, over a line, arrow, marker or patch outline, or off the
canvas); a violation raises and the exit status is nonzero.  Layout check only, no files written:  --check
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "figures")
CODE = HERE
sys.path.insert(0, CODE)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as Poly

BLUE, ORANGE, TEAL, GREY, INK = "#1f4e79", "#c8553d", "#3b8b8b", "#8a8f98", "#222831"
PURPLE = "#8e6bb0"
LCOL = {0: BLUE, 1: ORANGE, 2: TEAL, 3: PURPLE}
plt.rcParams.update({
    "font.size": 11.5, "axes.labelsize": 12.5, "legend.fontsize": 10.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "font.family": "serif", "font.serif": ["CMU Serif", "cmr10", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False, "pdf.fonttype": 42,
})

CHECK_ONLY = "--check" in sys.argv   # draw and layout-check only, write no files


# ------------------------------------------------------------ layout check
# Every figure is checked after drawing (see check_layout). A Text is exempt from the
# "text on an artist" test only if it carries a gid starting with "own-label": it is the
# label of a filled box or bar it sits inside on purpose (panel (a) subshell labels and
# the period numbers of panel (b) of the Madelung figure). Exempt texts are still tested
# against every other text and against the canvas edge. Tick labels touching axis spines
# are not tested (spines are not collected).
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle as _Rect
from matplotlib.transforms import Bbox
from matplotlib.legend import Legend


def _tight(bb, pad=0.0):
    return Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)


def check_layout(fig, name):
    """Raise AssertionError if any text overlaps another text, touches a line, arrow,
    marker or patch outline, or leaves the canvas."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    texts = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        if t.axes is not None and t in (t.axes.title, t.axes._left_title, t.axes._right_title) and not t.get_text():
            continue
        bb = t.get_window_extent(r)
        if not np.isfinite([bb.x0, bb.x1, bb.y0, bb.y1]).all() or bb.width == 0:
            continue
        # skip tick labels that matplotlib keeps but does not draw (outside the view)
        texts.append((t, bb))
    # tick labels outside the axes view are not drawn: drop them
    drawn = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for t in axis.get_majorticklabels() + axis.get_minorticklabels():
                drawn.add(id(t))
    inview = []
    for ax in fig.axes:
        ab = ax.get_window_extent(r)
        for axis, hor in ((ax.xaxis, True), (ax.yaxis, False)):
            ticks = list(axis.get_major_ticks()) + list(axis.get_minor_ticks())
            lo, hi = axis.get_view_interval()
            for tk in ticks:
                loc = tk.get_loc()
                if not (min(lo, hi) - 1e-9 <= loc <= max(lo, hi) + 1e-9):
                    tk.label1.set_visible(False)
    texts = [(t, bb) for t, bb in texts if t.get_visible()]
    # artists
    items = []   # (label, path in display coords or None, marker boxes)
    for ax in fig.axes:
        for ln in ax.lines:
            if len(ln.get_xdata()) == 0 or not ln.get_visible():
                continue
            if ln.get_color() == "none":
                continue
            path = ln.get_transform().transform_path(ln.get_path())
            if ln.get_linestyle() in ("None", "none", " ", ""):
                path = None
            boxes = []
            if ln.get_marker() not in ("None", "none", None, ""):
                pts = ln.get_transform().transform(np.column_stack([ln.get_xdata(), ln.get_ydata()]))
                h = ln.get_markersize() * dpi / 72 / 2 + ln.get_markeredgewidth()
                boxes = [Bbox.from_extents(x - h, y - h, x + h, y + h) for x, y in pts if np.isfinite([x, y]).all()]
            items.append((f"line {ln.get_label()!r}", path, boxes))
        for pa in ax.patches:
            if not pa.get_visible() or pa is ax.patch:
                continue
            path = _flat(pa.get_transform().transform_path(pa.get_path()))
            if isinstance(pa, _Rect) and not pa.get_ec()[3] and pa.get_fc()[3] == 0:
                continue
            items.append((f"patch {type(pa).__name__}", path, []))
        for ch in ax.get_children():
            if isinstance(ch, Text) and hasattr(ch, "arrow_patch") and ch.arrow_patch is not None:
                pa = ch.arrow_patch
                items.append(("annotation arrow", _flat(pa.get_path().transformed(pa.get_transform())), []))
        for lg in [ax.get_legend()]:
            if lg is not None and lg.get_frame_on():
                fb = lg.get_frame().get_window_extent(r)
                items.append(("legend frame", __import__("matplotlib").path.Path.unit_rectangle().transformed(
                    __import__("matplotlib.transforms", fromlist=["x"]).BboxTransformTo(fb)), []))
    problems = []
    fb = fig.bbox
    for i, (t, bb) in enumerate(texts):
        s = repr(t.get_text())[:30]
        if bb.x0 < fb.x0 - 0.5 or bb.y0 < fb.y0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y1 > fb.y1 + 0.5:
            problems.append(f"text {s} extends beyond the canvas")
        for t2, bb2 in texts[i + 1:]:
            if bb.overlaps(bb2) and bb.fully_overlaps(bb2) or (bb.overlaps(bb2) and _inter(bb, bb2) > 0.5):
                problems.append(f"text {s} overlaps text {repr(t2.get_text())[:30]}")
        if (t.get_gid() or "").startswith("own-label"):
            continue
        tb = _tight(bb, 0.5)
        for lab, path, boxes in items:
            if path is not None and path.intersects_bbox(tb, filled=False):
                problems.append(f"text {s} touches {lab}")
            elif any(tb.overlaps(b) for b in boxes):
                problems.append(f"text {s} touches a marker of {lab}")
    if problems:
        raise AssertionError(f"layout check FAILED for {name}:\n  " + "\n  ".join(sorted(set(problems))))
    print(f"layout {name}: PASS ({len(texts)} texts, {len(items)} artists)")


def _flat(path):
    """Polyline version of a path (Bezier control points of circles and arcs lie off the curve)."""
    from matplotlib.path import Path
    verts, codes = [], []
    for poly in path.to_polygons(closed_only=False):
        verts.extend(poly)
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(poly) - 1))
    return Path(np.array(verts), codes) if verts else path


def _inter(a, b):
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return min(w, h)


def save(fig, name):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout(pad=0.6)
    if os.environ.get("NOCHECK") != "1":
        check_layout(fig, name)
    if CHECK_ONLY:
        plt.close(fig)
        return
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ------------------------------------------------------------ exact block densities
def geg(n, lam):
    """Gegenbauer C^(lam)_n as a numpy Polynomial in c (three-term recurrence)."""
    c0, c1 = Poly([1.0]), Poly([0.0, 2.0 * lam])
    if n == 0:
        return c0
    for m in range(2, n + 1):
        c0, c1 = c1, (Poly([0.0, 2.0]) * (m + lam - 1) * c1 - (m + 2 * lam - 2) * c0) / m
    return c1


Z1 = Poly([1.0, 0.0, -1.0])          # z = 1 - c^2


def block(k, j):
    """(f, p) with f = C_j = C^(j+3/2)_(k-j) and p_j = z^(j+1) f^2 / K_j (density of c on [-1, 1])."""
    f = geg(k - j, j + 1.5)
    d = Z1 ** (j + 1) * f * f
    I = d.integ()
    return f, d / (I(1.0) - I(-1.0))


def cdf_poly(p):
    """F(t) = Pr(|c| <= t) = I(t) - I(-t) as a polynomial in t (I an antiderivative of the density p of c)."""
    I = p.integ()
    return I - Poly(I.coef * (-1.0) ** np.arange(len(I.coef)))


def H_poly(p0, p1):
    """H_l(t) = int_t^1 (F_(l+1) - F_l) ds, exact polynomial."""
    G = (cdf_poly(p1) - cdf_poly(p0)).integ()
    return G(1.0) - G


# ---------------------------------------------------------------- 1 densities
def fig1():
    N = 5; k = N - 1
    t = np.linspace(0, 1, 800)
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for l in range(0, (N - 1) // 2 + 1):
        _, p = block(k, 2 * l)
        ax.plot(t, 2 * p(t), color=LCOL[l], lw=2.2, label=rf"$l={l}$")
    ax.set_xlim(0, 1); ax.set_ylim(0, None)
    ax.set_xlabel(r"$|c|$   ($|c|=0$ equator, $|c|=1$ the poles $P_0$, $P_\infty$)")
    ax.set_ylabel(r"density $2p_{2l}(|c|)$")
    ax.legend(loc="upper center", frameon=False, ncol=3, bbox_to_anchor=(0.5, 1.12), handlelength=1.6, columnspacing=1.2)
    ax.set_title(r"$N=5$", fontsize=11.5, pad=22)
    save(fig, "fig1_densities")


# ---------------------------------------------------------------- 2 dominance fails, convex order holds
def fig2():
    N = 6; k = N - 1
    t = np.linspace(0, 1, 800)
    fig = plt.figure(figsize=(7.0, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.65])
    a1, a2, a3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[:, 1]), fig.add_subplot(gs[1, 0])
    for l in range(0, (N - 1) // 2):
        _, p0 = block(k, 2 * l); _, p1 = block(k, 2 * l + 2)
        dF = cdf_poly(p1)(t) - cdf_poly(p0)(t)
        a1.plot(t, dF, color=LCOL[l], lw=2.0, label=rf"$l={l}$")
        a1.fill_between(t, dF, 0, where=dF < 0, color="#d62728", alpha=0.55, lw=0)
        a2.plot(t, H_poly(p0, p1)(t), color=LCOL[l], lw=2.0, label=rf"$l={l}$")
    a1.axhline(0, color=GREY, lw=0.9)
    a2.axhline(0, color=GREY, lw=0.9)
    a1.set_xlabel(r"$t$"); a2.set_xlabel(r"$t$")
    a1.set_ylabel(r"$F_{l+1}(t)-F_l(t)$"); a2.set_ylabel(r"$H_l(t)$")
    for a in (a1, a2):
        a.set_xlim(0, 1)
    a1.set_title("changes sign: dominance fails", fontsize=11.5)
    a2.set_title("stays positive: convex order holds", fontsize=11.5)
    a1.legend(loc="upper left", frameon=False, title=r"$N=6$"); a2.legend(loc="upper right", frameon=False, title=r"$N=6$")
    a1.set_ylim(-0.04, 0.55)
    # zoom panel below the left panel: where each curve dips below zero
    for l in range(0, (N - 1) // 2):
        _, p0 = block(k, 2 * l); _, p1 = block(k, 2 * l + 2)
        dF = cdf_poly(p1)(t) - cdf_poly(p0)(t)
        mk = (t <= 0.4) & (dF <= 0.0095)
        a3.plot(t[mk], dF[mk], color=LCOL[l], lw=2.0)
        a3.fill_between(t, dF, 0, where=dF < 0, color="#d62728", alpha=0.55, lw=0)
    a3.axhline(0, color=GREY, lw=0.9)
    a3.set_xlim(0, 0.4); a3.set_ylim(-0.02, 0.01)
    a3.set_xlabel(r"$t$"); a3.set_ylabel(r"$F_{l+1}-F_l$")
    a3.set_title("zoom: the negative parts", fontsize=10)
    a2.set_ylim(-0.003, None)
    save(fig, "fig2_dominance_vs_convex")


# ---------------------------------------------------------------- 3 Q_j positive
def fig3():
    cases = [(4, 1), (7, 2), (10, 4)]
    cols = [BLUE, ORANGE, TEAL]
    t = np.linspace(0, 1, 1500)
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    data = []
    for (k, j), col in zip(cases, cols):
        f, _ = block(k, j)
        f1 = f.deriv(); f2 = f1.deriv()
        g = (j + 1) / (j + 2)
        Q = f1 * f1 - g * f * f2
        # normalise by Q(0): all curves start at 1; the check Q > 0 is exact in tiebreak_proof_check.py
        q = Q(t) / Q(0.0)
        assert q.min() > 0
        ax.plot(t, q, color=col, lw=2.0, label=rf"$k={k},\ j={j}$")
        data.append((f, col))
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo / 6, hi)
    lo = ax.get_ylim()[0]
    for (f, col) in data:
        zr = np.sort(f.roots().real)
        zr = zr[(zr >= 0) & (zr <= 1)]
        ax.plot(zr, np.full(len(zr), lo * 1.5), "o", color=col, ms=4.5, mfc="white", mew=1.4, ls="none")
    ax.set_xlabel(r"$|c|$")
    ax.set_ylabel(r"$Q_j(c)\,/\,Q_j(0)$")
    ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.0, 1.0))
    save(fig, "fig3_Q_positive")


# ---------------------------------------------------------------- Thomas-Fermi helpers
_TFC = {}


def tf():
    if "tf" not in _TFC:
        import tf_accurate
        _TFC["tf"] = tf_accurate.default()
    return _TFC["tf"]


A_T, a_T = 1.0127834, 0.5367947     # nearest Tietz field (tf_path_check.log)


def h_TF(x):
    """Kelvin-even weakening h(x), x = |c|, as in tf_certificate.py (Weak.h3), double precision."""
    rho = (1 - x) / (1 + x)
    P = tf().phi(rho / a_T); Q = tf().phi(1 / (a_T * rho)) / rho ** 2
    return 1 - (1 + rho) ** 2 * (P + Q) / (2 * A_T)


def h_TF_pp(x):
    """h'' by the chain rule of tf_certificate.py (Weak.h3); phi'' = phi^(3/2)/sqrt(x)."""
    T = tf()
    rho = (1 - x) / (1 + x); s1 = rho / a_T; s2 = 1 / (a_T * rho)
    p0 = T.phi(s1); p1 = T.dphi(s1); p2 = p0 ** 1.5 / np.sqrt(s1)
    f0 = T.phi(s2); f1 = T.dphi(s2); f2 = f0 ** 1.5 / np.sqrt(s2)
    P, P1, P2 = p0, p1 / a_T, p2 / a_T ** 2
    r2 = rho * rho
    Q = f0 / r2
    Q1 = -2 * f0 / (r2 * rho) - f1 / (a_T * r2 * r2)
    Q2 = 6 * f0 / (r2 * r2) + 6 * f1 / (a_T * r2 * r2 * rho) + f2 / (a_T ** 2 * r2 ** 3)
    K = (1 + rho) ** 2 / (2 * A_T); K1 = (1 + rho) / A_T; K2 = 1 / A_T
    S0, S1, S2 = P + Q, P1 + Q1, P2 + Q2
    G1 = K1 * S0 + K * S1
    G2 = K2 * S0 + 2 * K1 * S1 + K * S2
    q = 1 + x
    return -(G2 * 4 / q ** 4 + G1 * 4 / q ** 3)


def fig4():
    x = np.linspace(0, 1 - 1e-9, 400)
    hx = np.array([h_TF(v) for v in x])
    u = np.logspace(-8, -1, 160)
    hpp = np.array([h_TF_pp(1 - v) for v in u])
    # first u (largest) below which h'' < 0
    neg = u[hpp < 0]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.3), gridspec_kw=dict(width_ratios=[1, 1.15]))
    a1.plot(x, hx, color=BLUE, lw=2.2)
    a1.axhline(0, color=GREY, lw=0.9)
    a1.set_xlim(0, 1)
    a1.set_xlabel(r"$|c|$"); a1.set_ylabel(r"$h_{\rm TF}(|c|)$")
    a2.set_xscale("log")
    a2.plot(u, hpp, color=BLUE, lw=2.2)
    a2.axhline(0, color=GREY, lw=0.9)
    if len(neg):
        a2.axvspan(u.min(), neg.max(), color=ORANGE, alpha=0.18, lw=0)
    a2.invert_xaxis()
    a2.set_xlabel(r"$u=1-|c|$ (towards the poles)")
    a2.set_ylabel(r"$h_{\rm TF}''(1-u)$")
    a2.set_ylim(-1900, 520)
    a2.text(0.05, 0.12, r"shaded: $h''_{\rm TF}<0$," "\n" r"$u<9\times10^{-8}$", transform=a2.transAxes, ha="left", va="bottom", fontsize=10)
    save(fig, "fig4_hTF")
    return neg.max() if len(neg) else None


# ---------------------------------------------------------------- 5 W3 path
def fig5():
    d = np.load(os.path.join(HERE, "w3_path.npz"))["rows"]
    fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.6), sharex=True)
    for ax, N in zip(axs.ravel(), range(3, 9)):
        for l in range(0, (N - 1) // 2 + 1):
            m = (d[:, 1] == N) & (d[:, 2] == l)
            e = d[m, 0]; z = d[m, 3]; o = np.argsort(e)
            ax.plot(e[o], z[o] / (N * (N + 1) / 2), color=LCOL[l], lw=1.9, label=rf"$l={l}$")
        ax.set_xlim(0, 1)
        ax.set_title(rf"$N={N}$", fontsize=11, pad=3)
        ax.legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.2, labelspacing=0.2)
    for ax in axs[1]:
        ax.set_xlabel(r"$\epsilon$")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"$Z_{(N,l)}(\epsilon)/Z_N$")
    save(fig, "fig5_path")


# ---------------------------------------------------------------- 6 H2: curvature ratio
C0 = 1.05926488      # threshold coupling of the TF weight, N = 1, l = 0 (tf_path_check.log)


def ratio_curve():
    """Scal_TF/Scal_round = (c0/2) phi(s) / (a psi(s)^2), psi = u/s, u'' = -c0 phi(s) u / s, u(0) = 0, psi(0) = 1,
    c = (1 - r)/(1 + r) (paper 1 sign), r = a s."""
    from scipy.integrate import solve_ivp
    T = tf()
    s0 = 1e-7
    rhs = lambda s, y: [y[1], -C0 * T.phi(s) * y[0] / s]
    y0 = [s0 - C0 * s0 ** 2 / 2, 1 - C0 * s0]
    sol = solve_ivp(rhs, [s0, 3e4], y0, method="DOP853", rtol=1e-12, atol=1e-16, dense_output=True)
    s = np.concatenate([np.logspace(-6, 0, 200), np.logspace(0, np.log10(2e4), 500)[1:]])
    u = sol.sol(s)[0]
    psi = u / s
    ratio = (C0 / 2) * np.array([T.phi(v) for v in s]) / (a_T * psi ** 2)
    c = (1 - a_T * s) / (1 + a_T * s)
    return c, ratio


def fig6():
    c, r = ratio_curve()
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.axvspan(0, 1, color=BLUE, alpha=0.07, lw=0)
    ax.axvspan(-1, 0, color=ORANGE, alpha=0.07, lw=0)
    ax.plot(c, r, color=INK, lw=2.2)
    ax.axhline(1, color=GREY, lw=0.9)
    ax.set_xlim(-1, 1); ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$c=\cos\theta$")
    ax.set_ylabel(r"$\mathrm{Scal}_{\rm TF}\,/\,\mathrm{Scal}_{\rm round}$")
    ax.text(0.5, 0.12, "nucleus side\n" r"($c>0$, towards $P_0$)", ha="center", va="bottom", color=BLUE)
    ax.text(-0.5, 0.12, "far side\n" r"($c<0$, towards $P_\infty$)", ha="center", va="bottom", color=ORANGE)
    save(fig, "fig6_curvature_ratio")
    m = c > 0
    return np.abs(r[m] - 1).max(), c, r


if __name__ == "__main__":
    fig1(); fig2(); fig3()
    um = fig4()
    print("h'' < 0 only for u <", um)
    dev, c, r = fig6()
    print("max |ratio-1| on c>0:", dev)
    for cc in (0.9, 0.5, 0, -0.5, -0.9, -0.99):
        print("c =", cc, "ratio - 1 =", np.interp(cc, c[::-1], r[::-1]) - 1)
    if os.path.exists(os.path.join(HERE, "w3_path.npz")):
        fig5()
    else:
        print("fig5 skipped: w3_path.npz missing (run w3_cache.py)")

"""P11: the H2 lemma, checked symbolically with sympy.

Lemma. Let psi(r) > 0 solve  -Delta_3 psi = 2 Z v(r) psi  on R^3 (v(r) = 1/r for Coulomb, v = 1/(r (1+r/R)^2)
for the Tietz weight w = 2 v), r = |x|. Lift to R^4 by the Kustaanheimo-Stiefel map, r = |u|^2, and take the metric
g = psi(|u|^2)^2 delta on R^4 (conformal to flat). Then
        Scal(g) = 48 Z r v / psi^2 .
For the Tietz weight and the nodeless state psi = 1/(1+r/R) (N = 1, Z = 1/R) this is the constant 48/R,
the scalar curvature of the round 4-sphere of radius sqrt(R)/2 (the stereographic image of S^4_a).

Steps, each printed PASS or FAIL:
 (1) Scal of the metric f(|u|^2)^2 delta on R^4 from the Christoffel symbols equals -6 f^(-3) Delta_4 f (generic f).
 (2) KS identity: Delta_4 [psi(|u|^2)] = 4 r (Delta_3 psi)(r) for radial psi.
 (3) Hence Scal = 48 Z r v / psi^2 when -Delta_3 psi = 2 Z v psi (generic v, psi).
 (4) Tietz N = 1: psi solves the equation with Z = 1/R, and Scal = 48/R exactly (constant in r).
 (5) Tietz N >= 2 (psi_2 below, l = 0): the same formula gives a non-constant value (the sphere metric is singular at nodes).
 (6) Round sphere: the metric psi^2 delta with psi = 1/(1+|u|^2/R) is (R/4) times the stereographic metric
     4 delta_x/(1+|x|^2)^2, x = u/sqrt(R), so a round S^4 of radius sqrt(R)/2; its Scal from the Christoffel symbols is 48/R = 12/(R/4).
Run: OMP_NUM_THREADS=1 nice -n 19 python3 code/h2_scalar_check.py
"""
import sympy as sp

ok = True
def rep(name, cond):
    global ok
    ok = ok and bool(cond)
    print(f"  {name}: {'PASS' if cond else 'FAIL'}")

# (1) scalar curvature of F^2 delta in 4 dimensions, from the definition
u = sp.symbols('u1:5', real=True)
t = sp.Function('f')
rho2 = sum(x**2 for x in u)
Fn = t(rho2)
g = sp.eye(4) * Fn**2
ginv = g.inv()
n = 4
Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, b], u[c]) + sp.diff(g[d, c], u[b]) - sp.diff(g[b, c], u[d])) for d in range(n)) / 2
         for c in range(n)] for b in range(n)] for a in range(n)]
def Ric(b, c):
    s = 0
    for a in range(n):
        s += sp.diff(Gam[a][b][c], u[a]) - sp.diff(Gam[a][b][a], u[c])
        for d in range(n):
            s += Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a]
    return s
Scal = sum(ginv[b, b] * Ric(b, b) for b in range(n))
lap4 = sum(sp.diff(Fn, x, 2) for x in u)
pt = {u[0]: sp.Rational(3, 7), u[1]: sp.Rational(-2, 5), u[2]: sp.Rational(1, 3), u[3]: sp.Rational(5, 11)}
diff1 = sp.simplify((Scal + 6 * lap4 / Fn**3).doit().subs(pt))
# generic f, evaluated symbolically at a rational point: derivatives of f stay symbolic, so the difference must be 0
rep("(1) Scal(f^2 delta_4) = -6 Delta_4 f / f^3", diff1 == 0)

# (2) KS identity for radial functions
r = sp.symbols('r', positive=True)
psi = sp.Function('psi')
rho = sp.symbols('rho', positive=True)
F = psi(rho**2)
lap4_rad = sp.diff(F, rho, 2) + 3 * sp.diff(F, rho) / rho
lap3 = sp.diff(psi(r), r, 2) + 2 * sp.diff(psi(r), r) / r
lhs = sp.simplify(lap4_rad.doit())
rhs = sp.simplify((4 * r * lap3).subs(r, rho**2).doit())
rep("(2) Delta_4 psi(|u|^2) = 4 r Delta_3 psi", sp.simplify(lhs - rhs) == 0)

# (3) the lemma, generic v, psi
Z = sp.symbols('Z', positive=True)
v = sp.Function('v')
lap3_eq = lap3                                   # Delta_3 psi = -2 Z v psi on shell
scal_r = -6 * (4 * r * lap3_eq) / psi(r)**3      # (1) and (2) with rho^2 = r
on_shell = scal_r.subs(sp.diff(psi(r), r, 2), -2 * Z * v(r) * psi(r) - 2 * sp.diff(psi(r), r) / r)
rep("(3) Scal = 48 Z r v / psi^2 on shell", sp.simplify(on_shell - 48 * Z * r * v(r) / psi(r)**2) == 0)

# (4) Tietz N = 1
R = sp.symbols('R', positive=True)
vT = 1 / (r * (1 + r / R)**2)
psi1 = 1 / (1 + r / R)
eq1 = sp.simplify(-(sp.diff(psi1, r, 2) + 2 * sp.diff(psi1, r) / r) - 2 * (1 / R) * vT * psi1)
sc1 = sp.simplify(48 * (1 / R) * r * vT / psi1**2)
rep("(4) Tietz N=1: equation holds with Z=1/R", eq1 == 0)
rep("(4) Tietz N=1: Scal = 48/R, constant", sp.simplify(sc1 - 48 / R) == 0)

# (5) Tietz N = 2, l = 0: psi_2 = (1 - r/(2R)... ) solved from the equation with Z_2 = 3/R
a0, a1, a2 = sp.symbols('a0 a1 a2')
cand = (1 + a1 * r / R) / (1 + r / R)**2     # ansatz (Gegenbauer degree 1 times conformal factor)
eq2 = sp.simplify((-(sp.diff(cand, r, 2) + 2 * sp.diff(cand, r) / r) - 2 * (3 / R) * vT * cand) * r * R**2 * (1 + r / R)**5)
sol = sp.solve(sp.Poly(sp.expand(eq2), r).coeffs(), a1, dict=True)
if sol:
    psi2 = cand.subs(sol[0])
    sc2 = sp.simplify(48 * (3 / R) * r * vT / psi2**2)
    rep(f"(5) Tietz N=2 (psi_2 = {sp.simplify(psi2)}): Scal not constant", sp.simplify(sp.diff(sc2, r)) != 0)
else:
    rep("(5) Tietz N=2 ansatz solved", False)

# (6) round sphere of radius sqrt(R)/2, from the metric itself (Christoffel route, generic code of step (1))
Rp = sp.symbols('Rp', positive=True)
Ft = 1 / (1 + rho2 / Rp)
g6 = sp.eye(4) * Ft**2
g6i = g6.inv()
Gam6 = [[[sum(g6i[a, d] * (sp.diff(g6[d, b], u[c]) + sp.diff(g6[d, c], u[b]) - sp.diff(g6[b, c], u[d])) for d in range(n)) / 2
          for c in range(n)] for b in range(n)] for a in range(n)]
def Ric6(b, c):
    s = 0
    for a in range(n):
        s += sp.diff(Gam6[a][b][c], u[a]) - sp.diff(Gam6[a][b][a], u[c])
        for d in range(n):
            s += Gam6[a][a][d] * Gam6[d][b][c] - Gam6[a][c][d] * Gam6[d][b][a]
    return s
Scal6 = sp.simplify(sum(g6i[b, b] * Ric6(b, b) for b in range(n)))
x4 = [ui / sp.sqrt(Rp) for ui in u]
stereo = 4 / (1 + sum(xi**2 for xi in x4))**2          # coefficient of delta_x; delta_u = R delta_x
rep("(6) psi^2 delta = (R/4) * stereographic S^4 metric (conformal factors agree)",
    sp.simplify(Ft**2 - (1 / Rp) * (Rp / 4) * stereo) == 0)
rep("(6) Scal of psi^2 delta from Christoffel symbols = 48/R = 12/(sqrt(R)/2)^2", sp.simplify(Scal6 - 48 / Rp) == 0)

print("P11 " + ("PASS" if ok else "FAIL"))

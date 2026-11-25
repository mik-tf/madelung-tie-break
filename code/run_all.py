#!/usr/bin/env python3
"""Run the checks of the paper and print PASS/FAIL per check. Run from the repository root:
    OMP_NUM_THREADS=1 nice -n 19 python3 code/run_all.py [--certificate] [--path] [--ascan]
Default: the quick checks (about 5 minutes on one core). Opt-in long runs:
    --certificate  tf_certificate.py 17 1/1000   (Corollary 4, N <= 17; about 4 minutes)
    --path         tf_path_check.py              (Proposition 5, N <= 12; about 12 minutes)
    --ascan        a_scan.py --certify           (interval certificates at a = 0.40 and 0.65; about 2 minutes)
Needs Python 3, numpy, scipy, sympy, mpmath, matplotlib. Exit status is nonzero if anything fails.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, *args):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *args], capture_output=True, text=True, cwd=HERE)
    return p.returncode, p.stdout + p.stderr


flags = set(sys.argv[1:])
rows = []

rc, out = run("tiebreak_proof_check.py")
rows.append(("T1  exact identities (G1)-(G3), Lemmas A, B, closed forms, N <= 30 (sympy)", rc == 0 and "ALL OK" in out))

rc, out = run("review_tiebreak_check.py")
rows.append(("T2  independent rebuild from psi_nlm in R^3, 90 pairs; first-order slope by shooting", rc == 0 and "ALL OK" in out))

rc, out = run("tiebreak_dominance.py")
# the script exits 1 by design: it reports that first-order dominance FAILS for 49 of 56 pairs
rows.append(("T3  dominance fails (distribution functions cross) but mean |c| decreases in l; sign counts",
             "0 violate the weaker order" in out and "OK for every pair" in out and "FAIL: N=4" in out))

rc, out = run("tf_accurate.py")
errs = [float(x) for x in re.findall(r"\|B - published\| = ([0-9.eE+-]+)", out)]
rows.append(("T4  Thomas-Fermi function to 60 digits; slope B agrees with the published value",
             rc == 0 and len(errs) > 0 and max(errs) < 1e-20))

rc, out = run("h2_scalar_check.py")
rows.append(("T5  Lemma D, symbolic scalar curvature in six steps (sympy)", rc == 0 and "P11 PASS" in out))

rc, out = run("tietz_extremal.py")
rows.append(("T6  the Tietz field as an extremal (parts a and c)", rc == 0 and "(a) PASS" in out and "(c) PASS" in out))

rc, out = run("a_scan.py")
ok = {}
for m in re.finditer(r"a = (\d\.\d+)\s+min D = \S+ at \(N,l\) = \(\d+,\d+\)\s+all positive: (True|False)", out):
    ok[round(float(m.group(1)), 3)] = m.group(2) == "True"
inside = [a for a in ok if 0.40 <= a <= 0.65]
rows.append(("T7  sign of D_(N,l) on the grid of Tietz scales: positive for 0.40 <= a <= 0.65, fails at 0.70",
             rc == 0 and len(inside) >= 5 and all(ok[a] for a in inside) and ok.get(0.7) is False))

rc, out = run("make_figures.py", "--check")
rows.append(("T8  figure layout check (no text overlaps text, lines, arrows, patches; inside canvas)",
             rc == 0 and out.count(": PASS") >= 6))

if "--ascan" in flags:
    rc, out = run("a_scan.py", "--certify")
    rows.append(("L1  a_scan.py --certify: certificates at a = 0.40 and a = 0.65",
                 rc == 0 and len(re.findall(r"all LB > 0: True", out)) == 2))
if "--certificate" in flags:
    rc, out = run("tf_certificate.py", "17", "1/1000")
    rows.append(("L2  Corollary 4: interval certificate for every group N <= 17", rc == 0 and "VERDICT: PASS" in out))
if "--path" in flags:
    rc, out = run("tf_path_check.py")
    rows.append(("L3  Proposition 5: no level crossing along the path, N <= 12", rc == 0 and "P8 PASS" in out))

width = max(len(r[0]) for r in rows)
print()
for name, ok_ in rows:
    print(f"{name.ljust(width)}  {'PASS' if ok_ else 'FAIL'}")
skipped = [f for f in ("--certificate", "--path", "--ascan") if f not in flags]
if skipped:
    print("\nnot run (opt-in): " + " ".join(skipped))
allok = all(ok_ for _, ok_ in rows)
print("\nALL PASS" if allok else "\nSOME FAILED")
sys.exit(0 if allok else 1)

# Predictions

Each prediction is a line `T<n>:` or `L<n>:` and is tested by the script named in brackets. Run all quick checks with `python3 code/run_all.py`; the long ones with `--certificate`, `--path`, `--ascan`.

T1: the identities (G1) to (G3), Lemmas A and B and the closed forms (4) and (5) hold exactly for all N <= 30 [code/tiebreak_proof_check.py]
T2: an independent rebuild of (5) from the states psi_nlm in R^3 agrees for 90 pairs with N <= 20, Q_j has no real root, and the first-order slope of Theorem 3 agrees with shooting to relative error 1e-9 [code/review_tiebreak_check.py]
T3: first-order dominance fails for 49 of the 56 consecutive pairs with N <= 16, and none violates the weaker order (mean |c| decreasing in l); the sign-change counts are k - 2l - 2 [code/tiebreak_dominance.py]
T4: phi_TF'(0) = -1.58807102261137531271868... to all digits listed in the paper [code/tf_accurate.py]
T5: Lemma D holds, in six symbolic steps [code/h2_scalar_check.py]
T6: the Tietz field is the extremal at mu = 1/2 with E = 1 to 1e-8 [code/tietz_extremal.py]
T7: the sign of D_(N,l) is positive on the grid for 0.40 <= a <= 0.65 and fails at a = 0.70 [code/a_scan.py]
T8: in each of the six figures no text overlaps another text, a line, an arrow, a marker or a patch outline, and every text lies inside the canvas [code/make_figures.py --check]
L1: the interval certificates at a = 0.40 and a = 0.65 close [code/a_scan.py --certify]
L2: Corollary 4 holds for every group N <= 17 [code/tf_certificate.py 17 1/1000]
L3: no level crossing inside a group along the path from the Tietz to the Thomas-Fermi field, N <= 12 [code/tf_path_check.py]

Falsifiers: a pair (N, l) where the exact computation of T1 or T2 fails; a certified lower bound of D_(N,l) that is not positive for some N <= 17; a level crossing inside a group on the path; a figure with overlapping labels.

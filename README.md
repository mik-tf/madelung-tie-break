# A first-order Madelung tie-break on the round four-sphere: Kelvin parity, convex order and the Thomas-Fermi field

Author: Michaël Perreault

Version 1.1.1.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24533673.svg)](https://doi.org/10.5281/zenodo.24533673)

In the Demkov-Ostrovsky model of the Madelung rule the levels with equal N = n + l are degenerate at zero energy. The rule's second part, lower n first within a group, is not explained by the model. This paper proves that within every group the laws of |cos θ| for consecutive l are in strict increasing convex order (not first-order dominance: the distribution functions cross), and deduces a first-order tie-break for perturbations of the Tietz field that are Kelvin-even, increasing and convex. For the Thomas-Fermi field an interval computation certifies the splitting in every group with N <= 17, and a shooting check finds no level crossing along the path from the Tietz field for N <= 12. It follows [The Demkov-Ostrovsky Madelung model as the circle-invariant conformal Laplacian of the round 4-sphere](https://doi.org/10.5281/zenodo.24522163).

- **Theorems 1 and 2.** Strict increasing convex order of the laws of |cos θ| for consecutive l in every group, with a closed form from Gegenbauer polynomials and the Laguerre inequality.
- **Theorem 3.** A perturbation depending on r alone, with Kelvin-even weakening increasing and convex in |cos θ|, splits every group lower n first at first order (near the Tietz field).
- **Corollary 4.** The Thomas-Fermi field splits every group with N <= 17 at first order (interval certificate).
- **Proposition 5.** No level crossing inside a group along the path from the Tietz to the Thomas-Fermi field, N <= 12 (numerical).

Paper: [PDF](a-first-order-madelung-tie-break-on-the-round-four-sphere-kelvin-parity-convex-order-and-the-thomas-fermi-field.pdf), source in [paper.md](paper.md). Predictions: [predictions.md](predictions.md).

## Licences

- Text (`paper.md`), figures and logs: CC BY 4.0, see [LICENSE-paper.md](LICENSE-paper.md).
- Code in `code/`: MIT, see [LICENSE](LICENSE).

## Re-run the checks

Python 3 with numpy, scipy, sympy, mpmath and matplotlib. From the repository root:

```
OMP_NUM_THREADS=1 nice -n 19 python3 code/run_all.py
```

It runs the quick checks and prints PASS or FAIL for each (about 5 minutes on one core). The long runs are opt-in: `--certificate` (Corollary 4, N <= 17), `--path` (Proposition 5) and `--ascan` (the interval certificates at the two ends of the range of the Tietz scale). Each script also runs alone from the directory `code/`. Figures: `python3 code/make_figures.py` writes `figures/`; with `--check` it only runs the layout check.

| Check | Command (from `code/`) | What it checks |
|---|---|---|
| T1 | `python3 tiebreak_proof_check.py` | exact identities and closed forms, N <= 30 |
| T2 | `python3 review_tiebreak_check.py` | independent rebuild from the hydrogenic states, first-order slope |
| T3 | `python3 tiebreak_dominance.py` | the distribution functions cross (exits 1 by design) |
| T4 | `python3 tf_accurate.py` | Thomas-Fermi function to 60 digits |
| T5 | `python3 h2_scalar_check.py` | Lemma D, symbolic |
| T6 | `python3 tietz_extremal.py` | the Tietz field as an extremal |
| T7 | `python3 a_scan.py` | dependence on the Tietz scale |
| T8 | `python3 make_figures.py --check` | figure layout |
| L1 to L3 | `a_scan.py --certify`, `tf_certificate.py 17 1/1000`, `tf_path_check.py` | long runs |

## Cite as

Perreault, M. A first-order Madelung tie-break on the round four-sphere: Kelvin parity, convex order and the Thomas-Fermi field, version 1.1.1. https://github.com/mik-tf/madelung-tie-break (DOI: [10.5281/zenodo.24533673](https://doi.org/10.5281/zenodo.24533673) (all versions; v1.0 is [10.5281/zenodo.24533674](https://doi.org/10.5281/zenodo.24533674))). See [CITATION.cff](CITATION.cff).

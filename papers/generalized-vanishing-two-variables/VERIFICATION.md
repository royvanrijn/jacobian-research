# Verification matrix

This manuscript is a conceptual proof, not a computer-assisted theorem.
The matrix records where each load-bearing step is proved and what it uses.

| Claim | Paper location | Input |
|---|---|---|
| Degree at most the lowest positive operator order | Proposition 2.3 | Split-symbol polarization and Duistermaat--van der Kallen |
| Initial strict horizontal gap | Lemma 3.1 | Complete polarization, Duistermaat--van der Kallen, Hall's marriage theorem |
| Factorial valuation at a good prime | Lemma 3.2 | Legendre's formula at an unramified prime ideal |
| Shifted unequal-weight faces cannot overlap | Proposition 4.1 | Number-field specialization, Frobenius, Lemma 3.2, and unique-minimum valuation |
| Unequal-weight pure-zero faces are horizontally disjoint | Corollary 4.2 | Proposition 4.1 |
| A common unequal threshold is terminal for every fixed multiplier | Proposition 5.1 | Equal-face coordinate gap and bounded weight defect |
| The global envelopes reach a common threshold | Section 6 | Finite piecewise-linear envelopes and the initial Hall gap |
| Unrestricted GVC in two variables | Theorem 1.1 | The preceding rows |
| Long's adjacent-coefficient mechanism | Section 7 | Long's published Gaussian moment identity |
| Homogeneous GVC(3) counterexample | Theorem 8.1 | Quadric identity, exact spherical coefficients, and the Reynolds--apolar identity |
| Full winding--profile--radial failure family | Theorem 9.1 | The same phase extraction with arbitrary endpoint contact and profile |
| `GVC(n)` and homogeneous `GVC(n)` hold exactly for `n <= 2` | Theorem 10.1 | Theorems 1.1 and 8.1, de Bondt's split-symbol theorem, and unused-variable padding |

No bounded search, generated artifact, or symbolic certificate is used in
the proof.  The Ferrers support statement in Section 12 is a consequence of
horizontal separation.  The dependency-free GVC(3) checker is an independent
bounded replay, not the all-order proof.

The partial Lean audit in [`../../formal/gvc`](../../formal/gvc) checks the
literal ternary polynomials and their degree-twelve homogeneity, the cusp
identity, the algebraic beta evaluation of the endpoint moment, the complete
binomial coefficient ladder for an arbitrary rational endpoint profile, the
literal multivariate profile family and its degree/order formula, and
positivity of the exact mixed scalar.  It also proves the coefficientwise
equal-degree apolar contraction, composition of differential symbols,
the full concrete Reynolds/Laurent phase extraction, coefficientwise
algebraic change of height variable, and identification with the endpoint
kernel.  The resulting bridge is constructed in Lean, yielding the pure
identity, exact mixed scalar, and literal ternary counterexample.  Coefficient
base change and unused-variable padding then prove GVC failure over every
characteristic-zero field in every finite dimension at least three.  Thus
Theorem 8.1 and the negative half of Theorem 10.1 are fully Lean-verified.
Lean also proves abstractly that binary GVC implies unary GVC.  On the binary
side it proves the rational `p`-adic lower and exact factorial valuations
from Lemma 3.2 and the negative-final-slope plus intermediate-value envelope
step.

For Theorem 9.1, the Lean interface assumes only the multivariate phase
coefficients and the theorem's stated nonzero-moment hypothesis; from those
it derives the pure identity, the complete exact multiplier ladder, and the
literal GVC counterexample.

The audit still exposes, but does not yet prove, the full-profile phase
bridge and the Hall/shifted-ray/no-reversal/common-threshold obligations for
the binary theorem.  In particular, the transfer of the rational valuation statement to
an unramified prime ideal remains part of the written number-field proof.  A
spherical integral realization of the coefficientwise Reynolds functional is
not checked, but the formal differential identity uses the algebraic top
contraction directly and does not require that realization.
Therefore Theorems 1.1 and 9.1 and the complete biconditional in Theorem 10.1
should not yet be described as fully Lean-verified.

Build and inspect the manuscript with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  papers/generalized-vanishing-two-variables/main.tex
python3 scripts/check_latex_log.py \
  papers/generalized-vanishing-two-variables/main.log
python3 scripts/verify_gvc3_homogeneous_counterexample.py
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
make verify-gvc-lean
```

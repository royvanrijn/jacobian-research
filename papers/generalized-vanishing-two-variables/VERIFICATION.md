# Verification matrix

This manuscript is a conceptual proof, not a computer-assisted theorem.
The matrix records where each load-bearing step is proved and what it uses.

| Claim | Paper location | Input |
|---|---|---|
| Degree at most the lowest positive operator order | Proposition 2.3 | Split-symbol polarization and Duistermaat--van der Kallen |
| Initial strict horizontal gap | Lemma 3.1 | Complete polarization, Duistermaat--van der Kallen, Hall's marriage theorem |
| Shifted unequal-weight faces cannot overlap | Proposition 4.1 | Rational-point exposure and one good-prime valuation argument |
| Unequal-weight pure-zero faces are horizontally disjoint | Corollary 4.2 | Proposition 4.1 |
| A common unequal threshold is terminal for every fixed multiplier | Proposition 5.1 | Equal-face coordinate gap and bounded weight defect |
| The global envelopes reach a common threshold | Section 6 | Finite piecewise-linear envelopes and the initial Hall gap |
| Unrestricted GVC in two variables | Theorem 1.1 | The preceding rows |
| Long's adjacent-coefficient mechanism | Section 7 | Long's published Gaussian moment identity |
| Homogeneous GVC(3) counterexample | Theorem 8.1 | Quadric identity, exact spherical coefficients, and the Reynolds--apolar identity |
| Full winding--profile--radial failure family | Theorem 9.1 | The same phase extraction with arbitrary endpoint contact and profile |
| `GVC(n)` and homogeneous `GVC(n)` hold exactly for `n <= 2` | Theorem 10.1 | Theorems 1.1 and 8.1, de Bondt's split-symbol theorem, and unused-variable padding |

No bounded search, generated artifact, or symbolic certificate is used in
the proof.  The degree-eight and degree-nine Ferrers calculations mentioned
in Section 12 are discovery evidence and regressions only.  The
dependency-free GVC(3) checker is an independent bounded replay, not the
all-order proof.

Build and inspect the manuscript with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  papers/generalized-vanishing-two-variables/main.tex
python3 scripts/check_latex_log.py \
  papers/generalized-vanishing-two-variables/main.log
python3 scripts/verify_gvc3_homogeneous_counterexample.py
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
```

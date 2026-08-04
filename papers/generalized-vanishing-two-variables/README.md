# The Generalized Vanishing Conjecture: the two-variable theorem and the first failing dimension

Standalone manuscript source for the unrestricted binary GVC theorem, the
homogeneous three-variable counterexample, and the exact classification
`GVC(n)` if and only if `n <= 2`.  It also includes the complete
winding--profile--radial failure family and the specialization for every
power `Delta^k`, `k >= 6`.

Build from the repository root with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  papers/generalized-vanishing-two-variables/main.tex
python3 scripts/check_latex_log.py \
  papers/generalized-vanishing-two-variables/main.log
```

The canonical research proof is
[`../../extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md`](../../extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md).
The supporting Hall-localization, shifted-ray, and common-threshold results
are developed in
[`../../extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md`](../../extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md).
The canonical counterexample proof is
[`../../extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md`](../../extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md),
the full family is
[`../../extended-geometry/CUSP_PROFILE_SUSPENSION_THEOREM.md`](../../extended-geometry/CUSP_PROFILE_SUSPENSION_THEOREM.md),
and its dimensional consequences are recorded in
[`../../extended-geometry/GVC3_HOMOGENEOUS_SPILLOVERS.md`](../../extended-geometry/GVC3_HOMOGENEOUS_SPILLOVERS.md).
The manuscript restates their load-bearing arguments so that it can be read
without the longer research notes.

The only deep external theorem used in the proof is the Duistermaat--van der
Kallen constant-term theorem.  The manuscript
also gives the all-order spherical-coefficient proof of the GVC(3)
counterexample and explains its origin in Christopher D. Long's Gaussian
counterexample.

Replay the independent bounded audit of the three-variable formulas with:

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
```

A partial Lean audit is in
[`../../formal/gvc`](../../formal/gvc).  It checks the concrete polynomial
definitions and homogeneity, cusp identity, algebraic beta evaluation,
all-order full cusp-profile coefficient ladder, the literal multivariate
profile family and its degree/order formula, rational `p`-adic factorial
valuation lemma, exact scalar nonvanishing, coefficientwise apolar
contraction and operator composition, the full concrete Reynolds/Laurent
phase extraction and endpoint-kernel identification, characteristic-zero
base change, unused-variable padding, and the negative-final-slope plus
continuous-envelope step.  Consequently Theorem 8.1 and the negative
`n >= 3` half of Theorem 10.1 are fully Lean-verified.  The development does
not yet formalize the number-field shifted-ray transfer, the
Hall/no-reversal parts of the binary theorem, or the arbitrary-profile
phase bridge, so it is not a complete formal verification of the whole
manuscript.

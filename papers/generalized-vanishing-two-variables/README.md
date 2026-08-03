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

The proof has one substantive external input: the Duistermaat--van der
Kallen constant-term theorem.  It uses no computer algebra.  The manuscript
also gives the all-order spherical-coefficient proof of the GVC(3)
counterexample and explains its origin in Christopher D. Long's Gaussian
counterexample.  The paper is an active internal draft and has not been
externally reviewed.

Replay the independent bounded audit of the three-variable formulas with:

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
```

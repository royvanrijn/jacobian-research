# The Generalized Vanishing Conjecture: the two-variable theorem and the first failing dimension

Standalone manuscript source for the unrestricted binary GVC theorem, the
homogeneous three-variable counterexample, and the exact classification
`GVC(n)` if and only if `n <= 2`.  It also includes the complete
winding--profile--radial failure family and the specialization for every
power `Delta^k`, `k >= 6`.

<!-- status-consumer: GVC2SC f31ee48fbbecd427 -->
<!-- status-consumer: GVC2OC 75a2a340b8aa099a -->

Section 6.1 also gives the binary support-certificate equivalence over the
original field and a finite decision procedure for rational inputs. It proves
the optimal operator-independent mixed cutoff `m > C_d deg Q`, where
`d = deg P` and `C_d = 1 + floor((d+1)^2/4)`, with an equality family
in every degree. The earlier bound `m > (d + deg Lambda) deg Q` also holds;
either bound can be used. The canonical
[proof and exact checker](../../extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md)
record the algorithm and its assurance boundary.

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
base change, unused-variable padding, the arbitrary-profile quadric
restriction, all even-phase coefficient extractions, the shifted primitive
identity, the negative-final-slope plus continuous-envelope step, and the
complete finite-support common-threshold cutoff of Proposition 5.1.  It also
checks the finite Hall-deficiency and binary direction-localization core of
Lemma 3.1, including the sharp `d - e + 1` factor count and the two
coordinate-free power-divisibility normal forms.
Consequently Theorems 8.1 and 9.1 and the negative `n >= 3` half of Theorem
10.2 are fully Lean-verified.  The profile theorem uses exactly the
manuscript's declared-degree condition `S.natDegree <= e`.  The development
does not yet formalize the translated Duistermaat--van der Kallen argument
which rules out a split-factor matching, the literal field-extension and
split-factor construction and coordinate-change/exact-quotient packaging of
Lemma 3.1, the number-field shifted-ray
transfer and its equal-face ordering, or the no-reversal/global-envelope
parts of the binary theorem.  Thus it is not a complete formal verification
of the whole manuscript.
The Hall-direction descent and optimal uniform bound in Section 6.1 are written
proofs and are not currently formalized.

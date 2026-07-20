# A three-dimensional Jacobian counterexample

This repository verifies and studies the polynomial map

\[
F(x,y,z)=\left(
(1+xy)^3z+y^2(1+xy)(4+3xy),
y+3x(1+xy)^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

Its central certificate is finite and exact:

\[
\det DF=-2
\]

and

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Thus `F` is everywhere étale over `C` but is not injective. Subject only to
checking the displayed arithmetic, it is a counterexample to the classical
Jacobian conjecture in dimension three. Appending identity coordinates gives
counterexamples in every dimension at least three. Its coordinate degrees are
`(7,6,4)`.

The repository is deliberately centered on this three-dimensional result and
its direct consequences. It contains no plane-search program.

## Quick verification

The smallest verifier uses only Python's standard library and implements its
own sparse polynomial arithmetic:

```bash
python3 scripts/verify_counterexample_independent.py
```

For the full exact suite:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
make verify
```

The logical core is `make verify-core`. The remaining targets verify the
inverse geometry and the quartic member of the weighted family.

## Inverse and global geometry

Writing a target as `(a,b,c)` and setting `t=y+1/x`, the inverse problem reduces
to

\[
cT^3-2T^2+bT-2a=0.
\]

If `r=3ct^2-4t+b`, the source is reconstructed by

\[
x=2/r,\qquad y=t-r/2,\qquad
z=5r^2/4-3tr/2-cr^3/8.
\]

Consequently the generic function-field degree is three. The discriminant is
`-4Q`, where

\[
Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

The exact image, fibers, and nonproperness set are

\[
F(\mathbb C^3)=\mathbb C^3\setminus\Gamma,\qquad S_F=V(Q),
\]

with

\[
\Gamma=V(3bc-4,12a-b^2).
\]

Fibers have `3`, `1`, and `0` affine points respectively on `Q != 0`, on
`Q = 0` away from `Gamma`, and on `Gamma`. See
[CONSTRUCTION.md](notes/CONSTRUCTION.md) and
[IMAGE_AND_NONPROPERNESS.md](notes/IMAGE_AND_NONPROPERNESS.md).

## Further results retained here

- Exact finite-field fiber distributions in every characteristic:
  [FINITE_FIELD_VALUE_DISTRIBUTION.md](notes/FINITE_FIELD_VALUE_DISTRIBUTION.md).
- Commuting inverse-Jacobian vector fields:
  [COMMUTING_FLOWS.md](notes/COMMUTING_FLOWS.md).
- A degree-14 polynomial with three nondegenerate global minima and an explicit
  Palais--Smale sequence:
  [GRADIENT_INFINITY.md](notes/GRADIENT_INFINITY.md).
- Exact image, boundary geometry, and full `S_4` monodromy for a quartic-sheet
  member of the weighted family:
  [QUARTIC_WEIGHTED_GEOMETRY.md](notes/QUARTIC_WEIGHTED_GEOMETRY.md).
- Universal inverse-discriminant normalization and full `S_n` monodromy for
  every characteristic-zero weighted seed:
  [WEIGHTED_SEED_THEOREM.md](notes/WEIGHTED_SEED_THEOREM.md).
- Exact images and nonproperness sets for the full canonical family
  `H_d(W)=W^d(1-W)`:
  [CANONICAL_FAMILY_IMAGE.md](notes/CANONICAL_FAMILY_IMAGE.md).
- Boundary geometry and the exact image theorem for deformations with one
  additional simple primitive zero:
  [DEFORMED_SEED_BOUNDARY.md](notes/DEFORMED_SEED_BOUNDARY.md).
- Exact omitted-value classification for arbitrary fixed seeds and a closed
  exceptional-locus equation for two additional simple zeros:
  [OMITTED_VALUE_CLASSIFICATION.md](notes/OMITTED_VALUE_CLASSIFICATION.md).
- Exact discriminant saturation and boundary strata for repeated, including
  nonsplit, extra primitive roots:
  [REPEATED_ROOT_BOUNDARY.md](notes/REPEATED_ROOT_BOUNDARY.md).
- Explicit 95-dimensional cubic-homogeneous and 510-dimensional Drużkowski
  counterexamples:
  [CUBIC_HOMOGENEOUS_REDUCTION.md](notes/CUBIC_HOMOGENEOUS_REDUCTION.md) and
  [CUBIC_LINEAR_REDUCTION.md](notes/CUBIC_LINEAR_REDUCTION.md).
- Audited implications for related conjectures:
  [DIRECT_CONSEQUENCES.md](notes/DIRECT_CONSEQUENCES.md).

## Additional commands

```bash
make verify-normal-forms
make scan-weighted-seeds
.venv/bin/python scripts/finite_field_distribution.py
.venv/bin/python scripts/commuting_flows.py
.venv/bin/python scripts/nonproper_fiber_benchmark.py
julia --project=. scripts/nonproper_fiber_homotopy.jl
```

The normal-form generators write large ignored JSON artifacts under `results/`.
Julia is optional and is used only for the numerical nonproper-fiber benchmark.

For a claim-by-claim status and audit order, see
[FACTS.md](notes/FACTS.md), [IMPLEMENTATION_STATUS.md](notes/IMPLEMENTATION_STATUS.md),
and [FIRST_CHECKLIST.md](notes/FIRST_CHECKLIST.md). Provenance and historical
priority remain distinct from the exact certificate; see
[PROVENANCE_AUDIT.md](notes/PROVENANCE_AUDIT.md).

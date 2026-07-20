# Status of the `(72,108)` plane target

Updated 20 July 2026.  Primary source: Guccione, Guccione, Horruitiner, and
Valqui, [Increasing the degree of a possible counterexample to the Jacobian
Conjecture from 100 to 108](https://arxiv.org/abs/2204.14178).

## What is now transcribed

The unresolved initial-corner data are `A0=(8,28)` and `(m,n)=(3,2)`.  The
last rational change in Proposition 4.3 is

\[
x\mapsto x^{-1},\qquad y\mapsto x^4y,
\]

and changes the bracket to `x^2` up to a nonzero scalar.  It produces two
reduced Newton-polygon cases:

1. `P`: `(0,0),(1,0),(8,14),(8,16),(0,8)`; `Q`:
   `(0,0),(2,1),(12,21),(12,24),(0,12)`.
2. `P`: `(0,0),(1,0),(8,14),(8,16)`; `Q`:
   `(0,0),(2,1),(12,21),(12,24)`.

The broad system builder enumerates every lattice coefficient in each polygon,
forms `[P,Q]-x^2`, and saturates by every required vertex coefficient so that a
solver cannot escape by shrinking the Newton polygon.

## Current system sizes and bounded runs

| reduced case | variables | bracket equations | vertex factors | bounded modular result |
|---|---:|---:|---:|---|
| with vertical vertices | 186 | 302 | 10 | not launched |
| without vertical vertices | 72 | 92 | 8 | F4SAT timeout at 30 s |

The timeout is only a complexity measurement.  It says nothing about existence
or nonexistence over the finite field or over characteristic zero.

An experimental common-leading-form substitution replaces the square/cube
boundary coefficients by two scales and two edge-root parameters.  It reduces
the two cases to 163 and 49 variables, respectively.  The 49-variable saturated
system also exceeded 30 seconds modulo `1000003`.  This refinement is marked
experimental until every leading-power implication used in Proposition 4.3 is
independently audited; it must not be treated as an exhaustive theorem about
the reduced polygons yet.

## Important normalization warning

The original affine collision normalization is not imposed in these reduced
coordinates.  The rational transformations used to obtain Proposition 4.3 do
not preserve the selected affine points.  Any collision constraint must be
transported through the charts explicitly; copying `F(0)=F(1,0)=0` into this
system would exclude valid possibilities without justification.

## What the project has and has not reached

- The neighboring `(9,27)` Section 5 elimination is now an exact executable
  regression oracle, including its characteristic-zero terminal relation.
- The two broad `(72,108)` reduced systems are now machine-generated from the
  published polygons rather than guessed from degrees 72 and 108.
- No plane counterexample candidate has been found.
- No exclusion of either `(72,108)` reduced case has been proved.
- The next mathematical task is to audit and encode all common-leading-form,
  valuation, and dicritical constraints that turn 49--186 generic coefficients
  into the thin recursive systems used by the paper's successful exclusions.

Relevant artifacts:

- `scripts/build_72_108_reduced_system.py`
- `targets/degree_72_108.json`
- `results/reduced_72_108_systems_broad.json`
- `results/reduced_72_108_systems_edgepowers.json`
- `results/systems/reduced_72_108_*_modp.input`

# Alternate-Q80 norm-twelve chart character sweep (2026-09-04)

<!-- status-consumer: EC-K3-R17-NORM12-ALTERNATE-CHART-CHARACTER-SWEEP 16c92ccfdeeb13df -->

## Status

This is an exact negative result.  It does **not** prove generic rank at least
19 or 20.

The inherited height-four layer is now complete on all ten alternate-Q80
norm-twelve markings.  It contains 1,198 degree-two covers when counts from
the ten distinct base charts are summed.  Within every chart, exact
squareclass normalization retains the rational constant factor and finds no
equal-cover collision, no internal three-character relation, and no product
match against the eleven committed rank-28 characters or `q_103b2`.

| chart | inherited covers | result |
|---|---:|---|
| `11952` | 121 | no collision or product closure |
| `08ab4` | 131 | no collision or product closure |
| `08f72` | 86 | no collision or product closure |
| `091e4` | 155 | no collision or product closure |
| `135b7` | 127 | no collision or product closure |
| `10f72` | 118 | no collision or product closure |
| `1183a` | 120 | no collision or product closure |
| `09952` | 125 | no collision or product closure |
| `098fc` | 95 | no collision or product closure |
| `0ae21` | 120 | no collision or product closure |

The complete smooth rational-bisection layer has also been constructed and
checked on all ten charts:

```text
11952, 08ab4, 08f72, 091e4, 135b7, 10f72, 1183a, 09952, 098fc, 0ae21.
```

Each is an exact `24 I1` equation with a saturated rank-17 section basis.
Each complete atlas contains 39,147 verified lifted sections and 39,147
distinct quadratic characters.  On every chart the exact checker finds zero
equal-cover collisions, zero internal three-character relations, and zero
matches to the twelve committed characters.  Thus this layer comprises
391,470 exact cover records across ten different base coordinates, but it
does not supply the requested common quadratic cover or `V4` character
triple.

The summary certificate is
[`elkies-k3-r17-norm12-alternate-chart-character-sweep-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-alternate-chart-character-sweep-v1.json).

## Factorless exact squareclasses

Some branch-polynomial contents contain large composite integers for which
general integer factorization is irrelevant and prohibitively slow.  The
streaming checker instead factors only the geometric polynomial support and
tests rational constants by the exact condition that a numerator and
denominator are positive perfect squares.  Equal-cover and product tests are
therefore exact without assuming a factorization of the content.

Every pair product is additionally compared against the eleven certified
rank-28 quartic characters and `q_103b2`.  These are formal variable-rename
comparisons; any positive match would still require the separate
Neron--Severi degree-two compatibility gate before transferring a known twist
section.  No formal match occurs here.

## Exact saturation of the four formerly obstructed markings

The direct two-neighbor equations on `135b7`, `10f72`, `09952`, and `0ae21`
initially came with seventeen rational sections spanning a finite-index
rank-17 sublattice.  The first independent-generator choice had determinant
2 on `135b7` and determinant 4 on the other three charts.

For a degree-one old section of height `h`, the norm-twelve trace `w` satisfies

\[
\langle w,v\rangle=h,
\qquad
(v-w)^2=12-h.
\]

The published R17 lattice is rootless, so height 10 is impossible; at height
12 positive definiteness forces `v=w`.  Consequently the exact height-eight
Fincke--Pohst shell together with this single equality case is the complete
old-section search through the Cauchy bound.  On each chart these degree-one
sections have rank 17 and saturation index 2.  Reducing the generator matrix
against the complete old-section set gives an explicit determinant-2 basis
on every chart; the earlier determinant-4 values were generator artifacts.
A scan of all 39,120 committed rational bisections still finds no degree-one
curve in the missing coset.

The missing coset is instead constructed on the child equation.  For a known
section \(R=2Q\), the possible \(x(Q)\) are the roots of

\[
z^4-4x(R)z^3-2Az^2-(4Ax(R)+8B)z+A^2-4Bx(R).
\]

On each of the four charts, exact short-vector enumeration selects a
height-4 vector in the missing coset whose double is an integral combination
of the determinant-2 input basis.  The duplication quartic has a rational
linear factor.  The resulting function-field point satisfies the Weierstrass
equation and doubles exactly to the declared combination; replacing one input
section gives a coordinate determinant of 1 and a height-Gram determinant of
948.

| chart | input index | rational-half height | output index |
|---|---:|---:|---:|
| `135b7` | 2 | 4 | 1 |
| `10f72` | 2 | 4 | 1 |
| `09952` | 2 | 4 | 1 |
| `0ae21` | 2 | 4 | 1 |

These exact halves remove the marking obstruction and make all four complete
39,147-class atlases replayable from saturated equation-level bases.

## Replay

```bash
.venv/bin/python \
  elkies-k3/scripts/audit_r17_norm12_alternate_chart_character_sweep.py \
  --check

.venv/bin/python \
  elkies-k3/scripts/search_r17_norm12_quadratic_character_closure_streaming.py \
  --source-label norm12-orbit-0ae21 \
  --priority-table artifacts/generated-results/elkies-k3-r17-norm12-0ae21-alternate-bisection-priority-v1.tsv \
  --input artifacts/generated-results/elkies-k3-r17-norm12-0ae21-alternate-bisections-full-v1.json \
  --output artifacts/generated-results/elkies-k3-r17-norm12-0ae21-streaming-character-closure-v1.json \
  --check

sage -python \
  elkies-k3/scripts/saturate_r17_norm12_direct_section_basis.sage \
  --source-label norm12-orbit-135b7 \
  --check
```

## Boundary

The requested conclusions remain `UNKNOWN`.  This certificate closes the
inherited layer on all ten alternate norm-twelve charts and the complete
smooth rational-bisection layer on all ten.  It does not cover
higher-arithmetic-genus degree-two curves, singular normalizations beyond the
separately certified bounds, or uncatalogued non-coboundary twist characters.

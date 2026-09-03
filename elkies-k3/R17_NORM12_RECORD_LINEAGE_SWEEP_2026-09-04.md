# Complete norm-twelve record-lineage sweep on published R17

Status: **exact 43-chart compilation, exact rational-`PGL2` quotient, bounded
record-curve exclusion, and exact five-fibre published-R17 lineage
realization**.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 00e39f6b05c2688a -->

## Outcome

All 43 certified norm-twelve shared-zero degree-two fibrations on the
published R17 surface now have polynomial short Weierstrass equations

```text
Y^2 = X^3 + A(u) X + B(u),
deg(A,B,Delta) = (8,12,24).
```

Their normalized `j`-maps form exactly six classes under rational changes of
base.  Curves 273 and 302 have no rational preimage in any class.  One
eight-chart published-R17 class contains all five lineage curves 351, 356,
376, 377, and 385.  Every one of those 40 chart/fibre matches has trivial
quadratic twist over `QQ`.

Thus curve 356 and its four companions are proved to be fibres of one common
rootless rank-17 K3 family.  Curve 302 is proved not to be a fibre of any
chart in this complete 43-member shared-zero norm-twelve atlas, including the
alternate-Q80 classes and the hidden `norm12-orbit-103b2` published-R17 copy.
This is an atlas-bounded theorem, not an exclusion from rootless fibrations
outside these 43 certified charts.

## The six `j`-map classes

Each displayed equivalence is checked by an exact rational-function identity
after the stored `PGL2(Q)` substitution.  Different rows are separated by
exact critical-value eliminants.

| representative | size | frame type | orbit labels without prefix |
| --- | ---: | --- | --- |
| `074d9` | 8 | published R17 | `074d9,08aaa,1104c,0b4c1,05de2,10f74,10e1c,0a2f9` |
| `08234` | 5 | published R17 | `08234,1f002,04749,00678,04724` |
| `0e80b` | 18 | published R17 | `0e80b,06617,04b19,04ec4,10ba1,0681f,1f222,0bf17,183a6,090d4,103b2,1aa9d,0666f,0f983,0f09b,0698f,06867,13654` |
| `11952` | 6 | alternate Q80 | `11952,08ab4,091e4,10f72,1183a,098fc` |
| `07ca9` | 2 | published R17 | `07ca9,05592` |
| `08f72` | 4 | alternate Q80 | `08f72,135b7,09952,0ae21` |

The class sizes sum to `43` and recover the certified frame distribution
`33+10`.  For each chart and target, the atlas stores the primitive
projective degree-24 equation for `j(u)=j(E)`.  A projective modular no-root
prime proves nearly every miss; the `11952`/curve-302 case is closed by exact
`QQ` rational-root factorization because none of the predeclared small primes
alone obstructs it.

## Exact lineage parameters and twists

On the representative `norm12-orbit-074d9` chart, the five parameters are:

| curve | `u` |
| ---: | --- |
| 351 | `-1328795128332249413091778368/9355667371956741789361767413` |
| 356 | `-970521251165766412987689802944/5201452365084059419448702694373` |
| 376 | `-2885452538662711051989717312/13394519325797029034300898515` |
| 377 | `-24894098859432097011003667776/89496335058901143374489629459` |
| 385 | `-4569968803191492789881348058432/20219201289082406615406548938643` |

For every chart/fibre pair, let `(A_u,B_u)` be the specialized chart and let
`(A_E,B_E)` be the standard short target.  The checker constructs exact
`q,s in QQ` satisfying

```text
s^2 = q,
A_E = q^2 A_u,
B_E = q^3 B_u.
```

Consequently

```text
x_E = q x_u,
y_E = s^3 y_u
```

is a `QQ`-isomorphism, not merely a quadratic-twist identification.  The
artifact also records the inverse translations to every original generalized
Weierstrass model and aligns the ordinate sign with all 17 specialized
sections.

## The seventeen generic sections

Ten labelled columns interpolate literally.  Seven elementary differences
complete the basis:

```text
P2, P3, P4, P5, P8, P11, P13, P15, P16, P17,
P1-P2, P1-P6, P1-P7, P1-P9, P1-P10, P1-P12, P1-P14.
```

The resulting `17 x 17` word matrix has determinant `-1`.  Five-fibre
interpolation gives exact polynomial coordinates with `deg x <= 4` and
`deg y <= 6`; every section identity holds over `QQ(u)`.  Rational base
changes and weighted Weierstrass scalings then give explicit polynomial
coordinates for all 17 sections on each of the eight native charts.

Because the discriminant is squarefree with `24 I1` fibres, Shioda's formula
computes the height pairing directly from section intersections with zero.
The resulting Gram is positive definite of rank 17, has determinant `948`,
has no norm-two vectors, and is integrally isometric by exact PARI `qfisom`
to the pinned published-R17 frame.  The recovered section subgroup is
therefore saturated.

This explains the older first-jet misses: the literal labelled basis and the
complete declared one-shear search were genuinely obstructed, while the
successful basis uses seven simultaneous differences and lies outside that
search boundary.

## Exact displayed exceptional quotients

The current public records are hash-pinned.  Exhaustive exact
`E(F_p)/2E(F_p)` signatures, together with a modular irreducibility witness
for the 2-division cubic, prove independence of every displayed point.  Since
the generic word matrix is unimodular on `P1,...,P17`, the quotient of the
full displayed subgroup by the specialized generic subgroup is:

| curve | displayed subgroup rank | quotient | free basis modulo generic |
| ---: | ---: | --- | --- |
| 351 | 25 | `Z^8` | `P18,...,P25` |
| 356 | 29 | `Z^12` | `P18,...,P29` |
| 376 | 22 | `Z^5` | `P18,...,P22` |
| 377 | 23 | `Z^6` | `P18,...,P23` |
| 385 | 29 | `Z^12` | `P18,...,P29` |

All Smith torsion factors are trivial.  These are exact quotients inside the
subgroups generated by the displayed public points.  They do not assert that
those subgroups are the complete Mordell--Weil groups, and no rank upper bound
is inferred.

## Reproduction

From the repository root:

```bash
sage -python \
  elkies-k3/scripts/compile_r17_norm12_record_lineage_atlas.sage

sage -python \
  elkies-k3/scripts/certify_r17_norm12_wgxli_lineage_fibres.sage

sage -python \
  elkies-k3/scripts/compile_r17_norm12_record_lineage_atlas.sage --check

sage -python \
  elkies-k3/scripts/certify_r17_norm12_wgxli_lineage_fibres.sage --check
```

The first artifact contains all 43 equations, normalized `j`-maps, exact
`PGL2(Q)` quotient data, target equations, and rational-root decisions:

[`../artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json).

The second contains all 40 `QQ`-isomorphisms, all eight native 17-section
bases, the exact height and integral-isometry certificate, the public-point
independence blocks, and the five quotient decompositions:

[`../artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json).

## Claim boundary

The sweep proves a complete miss only within the 43 certified norm-twelve
shared-zero degree-two charts.  It does not exclude curve 273 or 302 from an
unlisted rootless fibration, prove exact ranks for the five lineage curves,
or produce rank 32.  Its positive result does provide five exact calibration
fibres and their exceptional complements for a targeted neighbourhood search.

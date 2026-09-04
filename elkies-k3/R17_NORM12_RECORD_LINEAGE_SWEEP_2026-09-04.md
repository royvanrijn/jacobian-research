# Complete norm-twelve record-lineage sweep on published R17

Status: **exact 43-chart compilation, exact rational-`PGL2` quotient, complete
474-curve pinned-ICARM preimage/twist sweep, and exact published-R17 and
alternate-Q80 specialization certificates**.

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 8a4c932153e2bb2d -->

## Outcome

All 43 certified norm-twelve shared-zero degree-two fibrations on the
published R17 surface now have polynomial short Weierstrass equations

```text
Y^2 = X^3 + A(u) X + B(u),
deg(A,B,Delta) = (8,12,24).
```

Their normalized `j`-maps form exactly six classes under rational changes of
base.  The complete pinned 2026-09-01 ICARM snapshot contains 474 curve
equations and 462 distinct `j`-invariants.  Exact projective preimage tests
against all six classes give 69 rational hits and 2,775 misses.  Every hit has
one rational root in one class, and all 376 induced native chart/fibre
comparisons have trivial quadratic twist over `QQ`.  All six classes occur.

In particular, every remaining wgxli component is recognized exactly:
`{363,364,378}` lies in alternate-Q80 class `08f72`, `{389,390,391}` lies in
published-R17 class `07ca9`, curve 393 lies in published-R17 class `0e80b`,
and curve 395 lies in alternate-Q80 class `11952`.  The original five curves
351, 356, 376, 377, and 385 remain in published-R17 class `074d9`.

ICARM curve 12, the independently certified Elkies--Klagsbrun
rank-at-least-29 curve, is an untwisted fibre of the native alternate-Q80
`11952` equation.  Exact specialization of its saturated generic MW17 basis
into the 29 displayed public points gives a free displayed quotient
`Z^12`.  This supplies the first high-rank native alternate-Q80 control.
There is no rank-30 or rank-31 hit: curves 273, 302, and 398 all miss the six
classes.

Thus curve 356 and its four companions are proved to be fibres of one common
rootless rank-17 K3 family, while the other wgxli numerical components are
now resolved as four further exact `j`-class lineages.  Curve 302 is proved
not to be a fibre of any chart in this complete 43-member shared-zero
norm-twelve atlas, including the alternate-Q80 classes and the hidden
`norm12-orbit-103b2` published-R17 copy.  This is an atlas-bounded theorem,
not an exclusion from rootless fibrations outside these 43 certified charts.

## The six `j`-map classes

Each displayed equivalence is checked by an exact rational-function identity
after the stored `PGL2(Q)` substitution.  Different rows are separated by
exact critical-value eliminants.

| representative | size | frame type | pinned ICARM hits | largest certified lower bound |
| --- | ---: | --- | ---: | ---: |
| `074d9` | 8 | published R17 | 5 | 29 |
| `08234` | 5 | published R17 | 54 | 28 |
| `0e80b` | 18 | published R17 | 2 | 27 |
| `11952` | 6 | alternate Q80 | 2 | 29 |
| `07ca9` | 2 | published R17 | 3 | 28 |
| `08f72` | 4 | alternate Q80 | 3 | 28 |

The class sizes sum to `43` and recover the certified frame distribution
`33+10`.  For each chart and target, the atlas stores the primitive
projective degree-24 equation for `j(u)=j(E)`.  A projective modular no-root
prime proves nearly every miss; the `11952`/curve-302 case is closed by exact
`QQ` rational-root factorization because none of the predeclared small primes
alone obstructs it.

The orbit labels in each class are unchanged from the original atlas artifact;
the representatives in table order have member counts `8,5,18,6,2,4`.

## Complete public-snapshot scan

The exact hit inventory by representative is:

```text
074d9: 351,356,376,377,385
08234: 11,65,66,67,394,396,397,412,413,414,416,417,418,419,
       420,421,422,423,427,428,429,430,431,432,433,434,435,436,
       437,438,439,440,441,442,443,445,446,447,448,449,450,451,
       452,453,454,455,456,457,458,459,462,463,464,465
0e80b: 393,404
11952: 12,395
07ca9: 389,390,391
08f72: 363,364,378
```

Of the 2,844 curve/class decisions, 2,698 misses have a stored prime at which
the projective equation has no root in `P1(F_p)`.  Exact `QQ` factorization is
used for the remaining 146 decisions: 69 have one rational root and 77 have
none.  The artifact stores a digest of every primitive equation, the full
primitive equation for every hit, all rational projective parameters, and the
twist comparison on every native member chart.

The live public endpoint grew after the hash-pinned 474-curve response.  The
certificate therefore retains precisely ids 1 through 474, checks that they
predate 2026-09-01, and requires the pinned SHA-256 of their
`(id,curve_key,ainvs,created_at)` projection.  That complete equation projection
is embedded in the output and is now the default offline replay source.  The
separate 69-fibre artifact retains each recognized curve's pinned point prefix
and is the default offline source for the five-fibre section and quotient
checker.  `--live-source` and `--live-pinned-source` are explicit source-drift
audits; later records are deliberately outside this theorem.

## Remaining wgxli components

The requested eight new wgxli recognitions are:

| curve | rank at least | representative | frame | representative parameter `u` |
| ---: | ---: | --- | --- | --- |
| 363 | 27 | `08f72` | alternate Q80 | `-1158455904970405917875280423497790960/7008997530131828358330671985973631951` |
| 364 | 28 | `08f72` | alternate Q80 | `-4889478984105233975746550959971939120/31980538646115718907687342513498498947` |
| 378 | 24 | `08f72` | alternate Q80 | `-566882880896877014074302259297950480/3961009332576056965142591776688723713` |
| 389 | 24 | `07ca9` | published R17 | `-2894205664468994496955959240/13458275699231869785327455027` |
| 390 | 26 | `07ca9` | published R17 | `-15016892529059225953186513980120/107318212874126569523668597319401` |
| 391 | 28 | `07ca9` | published R17 | `-45893899045288706850898992590520/317459242949710578550437826245181` |
| 393 | 26 | `0e80b` | published R17 | `468216440440228603200/743386878825686635399` |
| 395 | 28 | `11952` | alternate Q80 | `-205030389514850913908989849692/84586289990466564607930545463` |

Every member-chart twist in these four classes is trivial.  Thus the earlier
height-residual components were detecting four genuine lineages, but not one
common family: they land in four different rational-`PGL2` `j`-classes.

## Native alternate-Q80 rank-29 control

Curve 12 occurs on the direct `norm12-orbit-11952` equation at

```text
u = -517322777134156443281121957597/93938104767880054864772987570.
```

The twist is trivial on all six member charts.  On the representative, the
seventeen stored saturated alternate-Q80 sections specialize to exact integer
combinations of the 29 independently certified public points.  The resulting
`29 x 17` coordinate matrix has rank 17 and Smith factors all equal to one.
Augmenting it by the public points

```text
P2, P11, P4, P3, P6, P8, P17, P10, P28, P24, P19, P15
```

has determinant `+/-1`.  Hence, inside the displayed subgroup,

```text
<P1,...,P29> / MW17_specialized = Z^12.
```

All seventeen recovered relations are replayed by exact elliptic-curve group
law.  The numerical height solve was used only to discover the integer matrix,
not to prove it.

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
  elkies-k3/scripts/certify_r17_norm12_icarm_database_sweep.sage

sage -python \
  elkies-k3/scripts/certify_r17_norm12_curve12_alternate_q80_quotient.sage

sage -python \
  elkies-k3/scripts/compile_r17_norm12_record_lineage_atlas.sage --check

sage -python \
  elkies-k3/scripts/certify_r17_norm12_wgxli_lineage_fibres.sage --check

sage -python \
  elkies-k3/scripts/certify_r17_norm12_icarm_database_sweep.sage --check

sage -python \
  elkies-k3/scripts/certify_r17_norm12_curve12_alternate_q80_quotient.sage --check
```

The combined status checker for the last two commands is:

```bash
sage -python \
  elkies-k3/scripts/verify_r17_norm12_icarm_database_and_curve12.sage
```

All four commands above now replay from committed exact projections by
default.  The repaired offline checks reproduce the 2,844 curve/class
decisions, all 69 hits, and the five exact section/quotient certificates
without consulting the mutable endpoint.

The first artifact contains all 43 equations, normalized `j`-maps, exact
`PGL2(Q)` quotient data, target equations, and rational-root decisions:

[`../artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json).

The second contains all 40 `QQ`-isomorphisms, all eight native 17-section
bases, the exact height and integral-isometry certificate, the public-point
independence blocks, and the five quotient decompositions:

[`../artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json).

The complete 474-curve decision ledger, including the 69 hits, 2,775 misses,
all 376 native twist tests, and the pinned equation projection, is:

[`../artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json).

The exact curve-12 native alternate-Q80 specialization and displayed
`Z^12` quotient are:

[`../artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json).

## Claim boundary

The sweep proves a complete decision only against the 43 certified
norm-twelve shared-zero degree-two charts and only for the 474 equations in
the pinned snapshot.  It does not exclude a missed curve from an unlisted
rootless fibration, cover later ICARM records, prove an exact rank, or produce
rank 32.  Exact displayed generic/exceptional quotients are currently proved
for the original five `074d9` curves and curve 12; the other new fibre
recognitions do not yet carry such quotient certificates.

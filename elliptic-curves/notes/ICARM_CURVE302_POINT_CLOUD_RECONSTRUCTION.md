# ICARM curve 302: exact point-cloud reconstruction probes

Status: **exact bounded negative reconstruction result, with curve 273 and a
known Fermigier--Mestre negative control**.  No construction provenance,
Mordell--Weil upper bound, or exclusion of a moving-section family is claimed.

## Bottom line

Curve 302 should not currently be modelled as a visible `17+14` specialization.
Exact good-reduction Kummer codes through prime 1000 have full column rank 31
over both `F_2` and `F_3`.  Modulo the span of the submitted first seventeen
points, the remaining fourteen points contribute all fourteen dimensions in
both codes.  Pairwise local-class similarity does not separate the first
seventeen, the remaining fourteen, and cross-block pairs.

The elementary coordinate fossils tested here are negative.  There are no
repeated pair sums, oriented pair differences or pair products among the
public `x`-coordinates; no three-term arithmetic or geometric progressions;
no pair with the same rational squareclass of `x`; and no pair whose absolute
`x`-difference is a rational square.  Rational interpolation in the public
submission index, fitted on seven frozen points and tested only on held-out
points, gets no exact held-out `X` or `Y` value at total degree at most six.

The largest raw congruence clusters are not evidence for a hidden global
family.  They occur at bad primes.  Of the 31 public `x`-numerators, 25 are
divisible by 5, 23 by 11, and 15 by 23; all three primes are bad for curve
302.  This agrees with strong local conditioning, but the exact component
calculation in
[`RECORD_CURVES_273_302_FIRST17_SUBGROUPS.md`](RECORD_CURVES_273_302_FIRST17_SUBGROUPS.md)
already shows that the first seventeen points surject onto the complete
bad-component product.  The fourteen-point quotient remains invisible there.

These misses rule out only the small architectures stated below.  The same
coordinate and fixed-`X` probes also miss the public point basis of curve 245,
whose Fermigier--Mestre parent is known exactly.  They therefore do not argue
against a construction with moving section coordinates or a nontrivial
change of Mordell--Weil basis.

## Canonical mod-2 and mod-3 finite Kummer codes

For every usable good prime `p <= 1000`, the replay enumerates `E(F_p)` and
forms

```text
E(F_p) / ell E(F_p),  ell = 2 or 3.
```

The temporary basis of this quotient is discarded.  The artifact retains the
reduced row-echelon form of the submitted-point image in point order, which is
invariant under a change of local quotient basis.  It also retains the
canonical relation-space dimension and the exact local row spaces.

| curve | code | informative primes | full image rank | first-17 rank | remaining rank modulo first 17 | last rank-increase prime |
|---|---:|---:|---:|---:|---:|---:|
| 302 | mod 2 | 91 | 31 | 17 | 14 | 311 |
| 302 | mod 3 | 64 | 31 | 17 | 14 | 523 |
| 273 | mod 2 | 102 | 30 | 17 | 13 | 257 |
| 273 | mod 3 | 74 | 30 | 17 | 13 | 353 |
| 245 control | mod 2 | 95 | 20 | -- | -- | 149 |
| 245 control | mod 3 | 69 | 20 | -- | -- | 211 |

Thus the combined relation space is zero in all six rows.  For curves 302 and
273 this is a second prime-descent independence witness; its purpose here is
structural.  Full rank means that there is no nonzero linear relation in
these finite Kummer products from which to read a construction.

For curve 302 the mean fraction of informative primes separating a pair is:

| pair class | mod 2 | mod 3 |
|---|---:|---:|
| within first 17 | 0.58282 | 0.66153 |
| within remaining 14 | 0.58761 | 0.67119 |
| across the boundary | 0.58076 | 0.67155 |

These are descriptive finite-code statistics, not probabilistic tests.  Their
useful conclusion is narrow: the submitted index-17 boundary is not visible
as a local Kummer-similarity boundary.

## A basis-invariant Mestre 2-cover fossil

A soluble genus-one quartic is a 2-cover of its Jacobian.  After choosing one
rational quartic point, the covariant map sends all rational quartic points
into one affine coset of `2E(Q)`.  This gives a sharper diagnostic than raw
coordinate similarity and survives a change of Mordell--Weil basis.

The reconstructed curve-245 Fermigier--Mestre control verifies the diagnostic
exactly.  In the twenty-point public basis, all twelve visible Mestre points
and Fermigier's extra quartic point have the common parity vector

```text
(1,0,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0).
```

The first twelve coordinate rows also sum to zero exactly, as required by the
known visible-point relation.

By contrast, the full-rank curve-302 mod-2 code separates every pair of
submitted points.  Therefore no two **raw submitted points** can be images of
rational points on one soluble quartic 2-cover.  This does not exclude a
Mestre-like origin: the relevant visible sections could be nontrivial
combinations of the submitted basis, or the construction could use a
different map.  It does show that a direct search for a twelve-point visible
subset of the printed list is impossible and should be replaced by a
same-mod-2-coset search among controlled Mordell--Weil combinations.

## Denominators and rational squareclasses

Eighteen curve-302 points have integral public `x`, including points 20 and
23.  The repeated nontrivial denominator-root class is

```text
sqrt(den(x)) = 3  at points 16, 21, 25.
```

It crosses the first-17 boundary.  The other denominator roots are

```text
5, 2, 22, 31, 7, 29, 43, 60257, 229, 2987,
```

with `2987=29*103`; each occurs once.  Hence the denominator data show a
simple-to-complicated ordering tendency, not a clean seventeen-section block.

The replay computes exact local squareclasses of the public `x`-coordinates
at every prime through 97.  Local repetitions are retained in the artifact,
but no two global `x` values have the same rational squareclass.  Among
small-prime numerator-divisibility clusters of size at least three, the only
good-reduction curve-302 cluster is five points at 17.  The much larger
clusters are precisely the bad-prime rows at 5, 11 and 23 described above.

## Genuine held-out interpolation

The interpolation protocol is fixed before evaluation:

1. use the one-based public submission index as parameter;
2. take seven alternating points as training data;
3. fit every rational function `P(i)/Q(i)` with
   `deg(P)+deg(Q) <= 6` separately to the canonical integral-short-model
   coordinates `X` and `Y`;
4. evaluate exact rational equality only on points not used for fitting.

The protocol is applied to all 31 points, the first 17, and the remaining 14.
Every curve-302 model has zero held-out hits for `X`, zero for `Y`, and hence
zero joint hits.  The same is true for curve 273 and the curve-245 control.
This excludes low-degree dependence on the submitted ordering; it says
nothing about a latent parameter after reordering or a nonlinear section
label.

## A first nontrivial fixed-X deformation obstruction

Unrestricted first-order deformation is vacuous.  On a short model

```text
Y^2 = X^3 + A X + B,
```

every point with nonzero `Y` lifts formally under arbitrary infinitesimal
changes of `A` and `B` by the implicit-function theorem.  The replay instead
tests the first exact algebraic closure condition

```text
A(t) = A + a t + c t^2,
B(t) = B + b t + d t^2,
X_i(t) = X_i,
Y_i(t) = Y_i + h_i t.
```

Coefficient comparison gives

```text
h_i = (a X_i + b)/(2 Y_i),
h_i^2 = c X_i + d.
```

For a fixed rational projective direction `[a:b]`, the points
`(X_i,h_i^2)` must therefore be collinear.  Each triple gives a quadratic
condition on `a/b`; the direction at infinity is checked separately.  This
is an exact held-out test: a triple proposes a direction, and every further
point on the line is an independent hit.

All `C(31,3)=4495` curve-302 triples were tested.  None admits a rational
projective direction, so no pencil of this form preserves even three supplied
points.  The analogous maxima are also two for curve 273 and the curve-245
control.  This rules out a constant-`X`, linear-section Mestre-like fossil in
the canonical short chart, but the negative control shows why it must not be
generalized to moving sections.

## Interpretation and next reconstruction layer

The present exact evidence says:

- both prime Kummer layers see 31 undivided directions rather than a visible
  structural core plus a tail;
- the denominator boundary is an ordering/simplicity boundary and crosses
  point 17;
- the conspicuous numerator clusters are predominantly bad-reduction data;
- no degree-six submitted-index interpolation or fixed-`X` quadratic pencil
  survives a single held-out point.

The next useful attack should allow both a Mordell--Weil basis change and
moving section coordinates.  It should first be calibrated on the **actual
transported twelve-dimensional generic subgroup** of curve 245, rather than
on its unrelated submitted basis.  Two concrete exact routes are now
well-defined:

1. differentiate the reconstructed Fermigier--Mestre family and its twelve
   transported sections at the curve-245 parameter, then determine which
   low-degree moving-`X` jet invariants survive a change of Weierstrass chart;
2. search bounded low-height combinations of the 31 curve-302 points for a
   latent twelve-point/six-pair quartic configuration, using the transported
   curve-245 generic subgroup as a positive control and curve 273 as a second
   target.

Neither route should be dimensioned at 17 or scored against the R17 lattice.

## Reproduction

From the repository root:

```bash
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_icarm_curve302_point_cloud.py --check
```

The compact replay artifact is
[`icarm_curve302_point_cloud_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve302_point_cloud_v1.json).
Its SHA-256 is

```text
26b7b646cb3d282bf0b5811a9e923a7143303ff56dafe59c69f735ca5b55fa3d
```

<!-- status-consumer: EC-ICARM-CURVE302-POINT-CLOUD 1e1eb37dd6d4350f -->

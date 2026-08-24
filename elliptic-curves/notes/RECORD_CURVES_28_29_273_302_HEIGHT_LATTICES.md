# Record curves 28, 29, 273 and 302: comparative height-lattice audit

Status: **high-precision numerical computation and calibrated bounded
structural evidence**.  This note does not identify a K3 specialization or
prove an isometry with the generic rootless Mordell--Weil lattice.

## Bottom line

The four public point lattices have been recomputed with PARI/GP at 100 decimal
digits and LLL-reduced.  A common, reproducible short-vector procedure finds a
primitive rank-17 candidate space in every curve.  Those spaces are real
structural concentrations: an out-of-sample check, not used in selecting them,
shows substantially more integral points inside the candidates than outside.
The enrichment is strongest for curves 273 and 302.

The stronger claim initially suggested by this computation does **not**
survive calibration.  A forced rank-17 space in ICARM curve 245, whose exact
parent is independently known to be the Fermigier--Mestre generic-rank-12
family, has almost the same normalized first-1,311 profile and accepts an
exact unimodular determinant-948 R17 basis-entry fit almost as well.  Mapping
all 1,311 exact R17 minimal lines through those fitted bases produces a very
broad height distribution, not an approximate specialized minimal shell.
Thus neither the nearest-vector profile nor the optimized R17 Gram fit
distinguishes R17 provenance.  They must not be cited as strong common-H3
evidence.

The normalized nearest-rank quantiles of the first 1,311 unoriented vectors
are:

| curve | minimum | 10% | 25% | median | 75% | 90% | vector 1,311 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rank 28 | 3.379 | 3.689 | 3.892 | 4.140 | 4.387 | 4.575 | 4.762 |
| rank 29 | 3.115 | 3.691 | 3.896 | 4.172 | 4.442 | 4.618 | 4.763 |
| curve 273 | 3.203 | 3.678 | 3.879 | 4.150 | 4.429 | 4.669 | 4.821 |
| curve 302 | 3.160 | 3.689 | 3.890 | 4.142 | 4.409 | 4.641 | 4.809 |

Here every height is divided by the determinant-forced scale `lambda` defined
by

```text
det(candidate core) = 948 * lambda^17.
```

In particular, the medians are all within `0.032` of one another and within
`0.172` of the generic minimal norm `4`; all six non-minimum quantile spreads
are below `0.10`.  This is a reproducible compatibility statistic, but the
negative control below shows that it is not a provenance discriminator.

It is not proof.  The first 1,311 specialized vectors are not separated from
vector 1,312 by a visible height gap, and no integral change of basis taking
their Gram exactly to the pinned `R17` Gram has been found.  Bounded
specialization error can smear the norm-four and higher generic shells into
one another.

There is selection bias because the search dimension was set to 17.  The
curve-245 negative control now measures one concrete consequence of that
bias.  The useful surviving observation is the out-of-sample integrality
enrichment, which detects a structured core but does not name its lattice.

## Public-point height lattices

The supplied rational point lists were used with their printed signs.  PARI's
`ellheightmatrix` computed the Neron--Tate pairing, and `qflllgram` returned a
unimodular column transformation `U` with reduced Gram `U^t H U`.

| curve | displayed independent points | shortest height | determinant of displayed subgroup | LLL diagonal range |
|---|---:|---:|---:|---:|
| rank 28 | 28 | 47.2659171617 | `3.857298234011609...e34` | 47.266--77.113 |
| rank 29 | 29 | 46.3601302163 | `1.433744182671713...e36` | 46.360--75.010 |
| curve 273 | 30 | 43.4922918319 | `1.072068660411565...e37` | 43.492--69.947 |
| curve 302 | 31 | 51.3285105293 | `5.520367374821894...e39` | 52.778--85.997 |

The complete 100-digit matrices, LLL transformations, reduced Grams and
determinants are in
[`record_height_lattices_28_29_273_302_v1.json`](../../artifacts/generated-results/elliptic-curves/record_height_lattices_28_29_273_302_v1.json).
Its SHA-256 is

```text
aa89dea9eb7cf547633a67b522bd9bef67868b0fe686942853953f0258a1472f
```

These are lattices of the displayed finite-index candidates.  For curves 273
and 302 only rank lower bounds are unconditional, so the phrase
"Mordell--Weil lattice" here must not be read as a certified basis of the full
free group.

## How the rank-17 candidates were selected

For each LLL-reduced ambient lattice, PARI enumerated all unoriented vectors up
to a declared height bound.  The search counted exact coordinate identities
of the form `a+b=c` and `a-b=c`; vectors with many visible additive relations
were sampled first.  RANSAC proposed a rational 17-space.  The winning space
was then treated exactly:

1. its integer kernel was computed;
2. the double integer kernel
   `S=matkerint(matkerint(B)^t)` saturated it in `Z^r`;
3. every reported short-vector membership was checked by the resulting exact
   integer equations;
4. its canonical Gram `S^t H S` was independently LLL-reduced.

| curve | height bound | all ambient lines | exact lines in rank-17 candidate | core determinant | `lambda` | log-RMS of LLL lengths from `4 lambda` |
|---|---:|---:|---:|---:|---:|---:|
| rank 28 | 60 | 2,423 | 787 | `3.546378114169868...e22` | 14.1694 | 0.1188 |
| rank 29 | 60 | 1,507 | 485 | `8.170831310196924...e22` | 14.8824 | 0.1324 |
| curve 273 | 65 | 5,936 | 971 | `6.794075535535679...e22` | 14.7218 | 0.0992 |
| curve 302 | 70 | 3,458 | 864 | `3.622820375387032...e23` | 16.2450 | 0.1091 |

The candidate Grams, saturated bases in public-point coordinates, core LLL
transforms, exact high-precision determinants and shell quantiles are in
[`record_rank17_core_candidates_v1.json`](../../artifacts/generated-results/elliptic-curves/record_rank17_core_candidates_v1.json),
SHA-256

```text
0e1f32f73bb80033aeb1e2bee55402e685f0bb9a9e4e84c53e3a8026791c2d55.
```

The rank-29 result is a numerical candidate for the specialized generic
subgroup: a primitive rank-17 space containing 485 of the 1,507 ambient
vectors below height 60, with determinant scale `lambda=14.882435...`.  Of
those 485 vectors, 346 (71.3%) give integral points, compared with 551 of
1,022 (53.9%) outside the space.  The public announcement did not give the 17
transported section coordinates, so this is not an identification of the
actual transported basis.

## Dimension elbow and exceptional directions

A common bounded dimension scan gave the following best counts.  These are
search diagnostics; only the selected rank-17 memberships in the preceding
table were subsequently replayed through exact kernel equations.

| curve | dim 15 | dim 16 | dim 17 | dim 18 | dim 19 |
|---|---:|---:|---:|---:|---:|
| rank 28 | 381 | 546 | 787 | 0* | 0* |
| rank 29 | 236 | 350 | 485 | 514 | 545 |
| curve 273 | 435 | 654 | 971 | 1,071 | 1,200 |
| curve 302 | 372 | 595 | 864 | 991 | 1,055 |

`*` For rank 28 the top additive sampling pool itself had rank 17; zero is not
a statement that no rank-18 ambient subspace exists.

All four profiles have a marked gain at dimension 17 followed by smaller
marginal gains.  Saturating the candidate gives primitive decompositions of
dimensions

```text
28 = 17 + 11
29 = 17 + 12
30 = 17 + 13
31 = 17 + 14.
```

The determinant ratios `det(ambient)/det(core)`, equivalently the determinant
of the orthogonally projected exceptional quotient, are approximately

```text
1.0877e12, 1.7547e13, 1.5779e14, 1.5238e16.
```

Their per-exceptional-dimension geometric scales are `12.42, 12.70, 12.36,
14.32`, again placing the first three curves together and curve 302 at a
somewhat larger specialization scale.

## Negative-control calibration and the surviving fingerprint

ICARM curve 245 is an appropriate adverse control because its public rank-20
subgroup is large enough to force a rank-17 subspace, while its parent is
reconstructed exactly in this repository as the Fermigier--Mestre
generic-rank-12 family.  Running the declared rank-17 search at height 28
finds a saturated candidate containing 473 of 1,928 short lines.  Its
determinant-forced scale is `6.741934...`, and its normalized first-1,311
quantiles are

```text
2.689, 3.589, 3.931, 4.356, 4.593, 4.685, 4.752.
```

These overlap the four record profiles.  A second calibration optimizes over
unimodular bases whose columns are drawn from the exact 1,311-line minimal
shell of R17.  Every fitted integral Gram is therefore exactly isometric to
R17; only the scalar and basis are optimized.  This fits the 289 entries of a
chosen Gram.  It does **not** assert that all 1,311 minimal lines map to short
vectors of the numerical target.

| curve | fitted scale | normalized entrywise RMS | maximum entry error |
|---|---:|---:|---:|
| rank 28 | 12.586 | 0.462 | 1.430 |
| rank 29 | 13.724 | 0.460 | 1.346 |
| curve 273 | 13.479 | 0.457 | 1.353 |
| curve 302 | 14.817 | 0.428 | 1.245 |
| curve 245 negative control | 5.594 | 0.474 | 1.176 |

The control RMS is only `1.049` times the four-record mean.  At the declared
`1.20` separation threshold, this is a failed negative control.

The full-shell replay exposes an additional failure.  After normalizing each
set of 1,311 transported heights to mean 4, the fitted-shell medians are
`2.425, 3.151, 2.359, 2.534, 1.882` for rank 28, rank 29, curve 273, curve 302
and the control, respectively.  Their RMS deviations from 4 are `4.556,
2.996, 4.569, 4.135, 4.850`; only `12.2%`--`23.6%` of the mapped lines lie
between normalized heights 3 and 5.  These are dispersed clouds, not recovered
specialized R17 shells.  The entrywise fit therefore does not recover
transported generic sections or even the declared minimal-shell geometry.

The control can now be checked against exact ground truth rather than merely
its known family label.  Transporting the twelve visible Mestre quartic
points and Fermigier's extra point at canonical `T=5801/10` to the public
curve, then projecting with 150-digit canonical pairings, gives integral
coordinates in the twenty published points.  Exact PARI group-law replay
verifies all thirteen identities.  The only relation is

```text
Q1 + Q2 + ... + Q12 = O,
```

and `Q13` is independent, so the specialized known generic subgroup has rank
12.  Its sum with the forced rank-17 candidate has rank 20, hence their
intersection has rank `12+17-20=9`, the smallest dimension possible.  The
heuristic candidate therefore misses the true generic subgroup as strongly
as the ambient dimension permits.  At height 28 the true subgroup contains
144 short lines, of which 112 (77.8%) are integral; outside it 1,077 of 1,784
(60.4%) are integral, giving odds ratio `2.30`.  This shows that integrality
enrichment can validate structure while neither recovering the true generic
dimension nor identifying the lattice.

Fitting the numerical cores directly, without naming R17, fails even more
clearly.  The relative entrywise RMS values into curve 302 are `0.303` from
curve 273, `0.302` from rank 29, `0.295` from rank 28, and `0.300` from the
curve-245 control.  The control/record-source mean ratio is `1.00009`.
Likewise, the control fits curve 273 with RMS `0.318`, slightly better than
rank 29 at `0.323`.  Thus bounded GL(17,Z)-plus-scale alignment cannot decide
whether 273 and 302 have the same unknown generic lattice.

A third search directly minimizes the coefficient of variation of all 1,311
mapped R17 minimal-vector heights, rather than the entries of a selected
basis Gram.  Sixteen deterministic restarts over unimodular bases drawn from
each bounded target short-vector cloud give:

| curve | best shell-height CV | normalized RMS from norm 4 | fraction in [3,5] |
|---|---:|---:|---:|
| rank 28 | 0.3447 | 1.3788 | 51.6% |
| rank 29, known R17 positive | 0.3420 | 1.3679 | 50.7% |
| curve 273 | 0.3378 | 1.3513 | 52.7% |
| curve 302 | 0.3514 | 1.4054 | 51.0% |
| curve 245 negative control | 0.3335 | 1.3341 | 54.6% |

The negative control is numerically *better* than rank 29 and all four record
candidates under this objective.  Consequently the shell-aware search also
fails calibration.  Curve 273's value cannot be promoted to R17 evidence, and
curve 302's slightly worse value is not evidence against R17.  The bounded
search, bases, restart distributions and shell quantiles are in
[`record_rank17_shell_embedding_search_v1.json`](../../artifacts/generated-results/elliptic-curves/record_rank17_shell_embedding_search_v1.json).
Its SHA-256 is

```text
34d9d287d37ee07a7caf862f0db53785a3f197164b3862fb8bb24fe1152f48b3.
```

One out-of-sample fingerprint does separate the candidates.  Point
integrality was not used by the additive-RANSAC selection:

| curve | integral inside core | integral outside | rate ratio | odds ratio |
|---|---:|---:|---:|---:|
| rank 28 | 441/787 = 56.0% | 666/1,636 = 40.7% | 1.38 | 1.86 |
| rank 29 | 346/485 = 71.3% | 551/1,022 = 53.9% | 1.32 | 2.13 |
| curve 273 | 446/971 = 45.9% | 905/4,965 = 18.2% | 2.52 | 3.81 |
| curve 302 | 396/864 = 45.8% | 570/2,594 = 22.0% | 2.09 | 3.00 |
| curve 245 control | 310/473 = 65.5% | 879/1,455 = 60.4% | 1.08 | 1.25 |

This says that the selected 17-spaces in the records, especially 273 and 302,
are not merely arbitrary coordinate slices of the short-vector clouds.  It
supports a stable section-core interpretation.  Integrality is not an
isometry invariant and can occur in many elliptic-surface constructions, so
it supplies no R17 or H3 identification.

The complete fitted unimodular bases, exact isometric R17 Grams, all mapped
shell dispersions, exact thirteen-point control reconstruction, control Gram,
and integrality counts are in
[`record_rank17_fingerprint_calibration_v1.json`](../../artifacts/generated-results/elliptic-curves/record_rank17_fingerprint_calibration_v1.json),
SHA-256

```text
6294fe759b7c03eb908afea12546c896d7137566d884c38caa13541de42ee763.
```

## Coefficients, discriminants, conductors and j

| curve | digits in `(a4,a6)` | digits in displayed `Delta` | `log |Delta|` | `log N` |
|---|---:|---:|---:|---:|
| rank 28 | (56,83) | 166 | 381.202 | 325.904 |
| rank 29 | (62,92) | 185 | 425.442 | 343.720 |
| curve 273 | (63,94) | 188 | 432.125 | 339.348 |
| curve 302 | (67,99) | 197 | 453.047 | 375.222 |

All four displayed discriminants contain `2,3,5,7,13`.  Additional overlaps
are consistent with the previously recorded CRT-shaped lineage:

- rank 28 and rank 29 also share `11,17`;
- rank 29 and curve 273 also share `31,41`;
- rank 29 and curve 302 also share `11,41`;
- curves 273 and 302 also share `41`.

This is not a family invariant.  Pairwise `j`-invariants are different; every
pair has `gcd(c6_i,c6_j)=1`, and every pair except rank 28/rank 29 has
`gcd(c4_i,c4_j)=1` (that exceptional gcd is only `37`).  Thus there is no
simple equality or twist pattern hidden by the large coefficients.

The published rank-28 Tate calculation reports exponent two at `3` and
exponent one at its other eleven bad primes; that long conductor replay was
not completed locally in this audit.  The public rank-29 account reports a
squarefree conductor on its seventeen bad primes.  Curve 273's exact
global-minimal conductor is the separately pinned 148-digit integer, while
curve 302's is the pinned 163-digit integer with conductor exponent two only
at `5`.  These local differences prevent a simple "same minimal model up to
scaling" explanation.

## Mestre, Nagao and family recognition

The exact normalized six-root Mestre census now includes curve 302 as well as
curve 273.  For each curve it tested all 2,329 nonsingular normalized root
tuples of diameter at most 300 plus the larger Fermigier control tuple.  Curve
273 had 113 exact-factorization survivors and curve 302 had 146; neither had a
rational-square parameter or an exact `j` match.  The full artifact is
[`icarm_construction_fingerprints_v2.json`](../../artifacts/generated-results/elliptic-curves/icarm_construction_fingerprints_v2.json),
SHA-256

```text
ced63ed67c61bb23484039237259127ffd0864426ae41429cd005e6989bfdc4a.
```

This excludes only the declared bounded fixed-root Mestre census.  There is
no exact Nagao, Kihara, generalized-Mestre or other-family recognition for
curves 273 or 302.  Their certified trivial torsion also excludes direct
specialization of the implemented Elkies--Klagsbrun model having forced point
`(0,0)` of order two, but does not exclude another fibration, an isogenous
image, or the rootless rank-17 descendant.

## Interpretation and exact next gate

The height-lattice evidence motivates using the H3 machinery against both
current record curves, but it no longer supports an R17 provenance claim by
itself.  Curves 273 and 302 contain unusually integral, additive rank-17
candidate spaces.  Their determinant scaling, LLL-length spread and
normalized 1,311-vector profiles are compatible with R17, while the calibrated
Gram-entry test is demonstrably non-specific and does not recover the exact
minimal shell.

The conclusion must remain:

```text
structured rank-17-core evidence; R17 provenance unresolved.
```

An exact decision still requires the generic rootless equation and `j`-map,
solving `j_R17(t)=j_273` and `j_R17(t)=j_302`, transporting all seventeen
sections, and proving the resulting curves Q-isomorphic to the public minimal
models.  If either equation has no rational solution, the candidate core in
this note becomes a concrete lattice target for the reverse-neighbour/RR/
quartic/Jacobian reconstruction path rather than evidence for H3 itself.

## Public sources

- [Dujella's rank-28 page](https://web.math.pmf.unizg.hr/~duje/tors/rk28.html)
  supplies the 2006 curve and 28 points.
- [Dujella's rank-29 page](https://web.math.pmf.unizg.hr/~duje/tors/rk29.html)
  supplies the 2024 curve and 29 points.
- [Dujella's rank-30 page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html)
  supplies curve 273 and its 30 points.
- [ICARM curve 302](https://elliptic-rank.icarm.cloud/curve/302) supplies the
  rank-at-least-31 model and 31 points.
- The [public rank-29 construction account](https://mathoverflow.net/questions/477849/background-for-the-elkies-klagsbrun-curve-of-rank-29/478509)
  states that the search used the same K3 as rank 28, with generic
  `Z^17` and twelve additional independent points.
- Klagsbrun--Sherman--Weigandt,
  [*The Elkies Curve Has Rank 28 Subject Only to GRH*](https://arxiv.org/abs/1606.07178),
  is the source of the published rank-28 local/conductor calculation.

## Reproduction

See [`../REPRODUCE.md`](../REPRODUCE.md).  The principal commands require
PARI/GP; the Mestre census additionally uses the pinned Python dependencies.

# Visibility, selection and a point-supplied construction

Completed bounded follow-up to the
[structural reassessment](RANK_JUMP_REASSESSMENT_2026-09-05.md).
The [portable exact replay](../../artifacts/generated-results/elliptic-curves/rank_jump_diagnostics_replay_v1.json)
passes. No parameter census is enlarged, no record search or descent is run,
and no existing gate or rank-status entry is changed. Mathematical status
remains in [MATH_STATUS.json](../../MATH_STATUS.json).

## 1. Explain the MW18 misses before increasing exposure

The [literal representative audit](../../artifacts/generated-results/elliptic-curves/mw18_retained_visibility_v1.json)
applies `search_observability.point_visibility` to both signs of the eleven
published MW17-complement points in every one of the 600 retained charts:
13,200 point/chart checks. One direction has become generic on each cover;
these eleven labels are not eleven independent directions beyond MW18.
Every literal representative is outside the height-100,000 box. The nearest
requires height 1,440,705; many require vastly more.

That alone would be misleading: the trial recovered other representatives
of some of these directions. The
[translated audit](../../artifacts/generated-results/elliptic-curves/mw18_translated_visibility_v1.json)
uses the specialized canonical height Gram to propose **one generic translate
per sign, representative and retained centre**, another 13,200 checks.
For centre Q and oracle R, the proposed translate is R+m, with m near
Q/2 minus the projection of R onto the generic subgroup. Floating CVP proposes
the word; exact group law, the complete chart transport and the final square
identity verify every retained witness. This is an oracle-only diagnostic,
not a new centre-selection rule or a proof of an optimal translate.

Exact expressions of the recovered points in the independent public basis
also distinguish a missing representative from a missing rational direction.
There are **no omissions inside claimed completed coverage**. All translated
points lying in it occur in the output, and the retained discoveries pass
the same independent observability check.

The following are **additional rational-span dimensions with explicit
translated witnesses inside the indicated boxes**, beyond the original
deepest-policy recovery. They are not results of larger searches.

| Presentation | Original gain beyond MW18 | H=200,000 | H=1,000,000 | H=10,000,000 |
|---|---:|---:|---:|---:|
| 531 / `12f61` | 6 | 0 | 2 | 3 |
| 534 / `13d7a` | 2 | 0 | 1 | 2 |
| 534 / `1a371` | 3 | 0 | 1 | 2 |
| 545 / `08c1e` | 4 | 1 | 2 | 4 |
| historical / `15a68` | 7 | 0 | 0 | 2 |
| Presentation-weighted total | 22 | 1 | 6 | 13 |

The two curve-534 presentations are correlated. A tenfold height bound can
cost roughly one hundred times the primitive numerator/denominator work;
this table is not an instruction to enumerate those boxes.

There are mixed causes, not a single detector bug:

- **Modest extensions suffice for some representatives.** The closest still
  missing translated representative on 545 has exact height 160,634; other
  new directions have witnesses below one million.
- **Others retain a large coordinate cost.** For the hardest missing basis
  representative among the forty deepest charts, the smallest proposed
  height is 297,016,773,067,303,433 on `534/13d7a`, and
  100,391,649,163,954,452,831 on `534/1a371`. These are minima over the
  declared proposals, not lower bounds over every generic translate or
  horizontal coordinate system. Further changes of representative or
  coordinates remain worth diagnosing before a box increase on this scale.
- **No missing direction is already inside the claimed box in this audit.**
  The rational-span extension at the original height is zero on every arm.

Diverse-deep differs in detail: its witnessed additions at one million are
`2,0,1,3,1`, relative to original gains `6,2,3,3,7`. In particular the
historical presentation has a missing witness at height 237,446 in that arm.
The certificates retain every label, generic word, signed point, chart index
and primitive coordinate. No adaptive wave has been run and the original
35/50, minimum-five gate remains failed.

## 2. Mask generic directions on the actual ordinary fibres

The experiment uses exactly the **31 previously searched distinct curves**
of the [height-population pilot](FIBRE_HEIGHT_POPULATION_2026-09-05.md):
15 MW16 and 16 MW18 fibres. It neither samples more parameters nor claims
that this previously selected subset is an unbiased population.

Before opening the oracle, section index zero is removed, and the principal
block of the original generic Gram is retained. For each family, a fixed
hash samples 256 parity classes of the remaining subgroup. Every sampled
coset minimum is checked exactly; the twelve deepest sampled classes are
frozen. This is not the complete deepest stratum. All curves use those same
twelve family-specific centres, `metric:16`, H=100,000 and 20 seconds per
chart. Four supervised workers have 900 seconds and 2 GiB per curve.
Separate files hold the blind input and oracle. The worker reads only the
blind points; centres remain fixed and there is no adaptive wave.

| Family | Curves | Completed charts | Withheld directions recovered | Distinct returned points, summed by curve |
|---|---:|---:|---:|---:|
| MW16 | 15 | 180 | 15 | 256 |
| MW18 | 16 | 192 | 16 | 32 |

All **31/31 withheld directions are recovered**. No literal signed withheld
point appears: its nearest retained-chart height has 91–447 bits for MW16
and 430–913 bits for MW18. The output instead supplies generic translates
or combinations. This is the same representative issue seen retrospectively
on the anchors.

The [search transcripts](../../artifacts/generated-results/elliptic-curves/ordinary_masked_controls_v1.json.gz)
contain the frozen protocol, exact inputs and all 372 completed charts.
The [relation certificate](../../artifacts/generated-results/elliptic-curves/ordinary_masked_relations_v1.json)
expresses **all 288 returned points** exactly in the original independent
subgroups. Finite reductions at moduli 3,5,7,11 propose bounded words;
exact group law accepts them. A pointed partner also supplies one word from
its centre and an already verified point. Every curve has a returned point
with nonzero coefficient on the withheld basis direction. No unresolved
relation remains. The portable replay separately verifies the original
31 independence certificates.

This endpoint is **WITHHELD_KNOWN_DIRECTIONS**, not new rank. Total worker
time is 927.84 seconds for MW16 and 1,003.02 for MW18, excluding preparation
and replay. The result rules out complete insensitivity to generic escapes
on these particular large models. A subgroup of rank r-1 with a withheld
generic point has different geometry from rank r with an exceptional point;
31/31 does not promise comparable exceptional-point sensitivity or decide
whether the earlier prospective null was caused by rare incidence.

## 3. Test the selector without enlarging the population

[Elkies, Lecture III](https://arxiv.org/html/0709.2908v1) describes thousands
of primes, millions of specializations and precomputed residue counts. That
account does not validate our much shorter screens. This experiment freezes
the existing **91 MW16 and 93 MW18** eligible classes and tests natural
prefixes of 25,64,128,256 primes starting at 5. The selection primes end at
1627; two disjoint 64-prime validation blocks run from 1637 through 2663.
A fixed hash also splits each family into development and curve holdout
subsets (45/46 and 46/47). These are within-family holdouts, not independent
fibration holdouts.

Two policies are tested without fitting:

\[
 S_N=\sum_{p\ \mathrm{good}}\frac{2-a_p}{p+1-a_p}\log p,
 \qquad S_M=\sum_{p\ \mathrm{good}}\log\frac{p+1-a_p}{p}.
\]

Each prime contribution is rounded to 10^-12 before summation. Bad primes
of the retained global minimal model contribute zero and their counts are
retained. The [trace artifact](../../artifacts/generated-results/elliptic-curves/bounded_prime_selector_traces_v1.json.gz)
stores 68,667 distinct `(prime, reduced equation)` counts for 70,656
candidate/prime lookups. Counts are shared for equal residue equations;
there is no unreported full residue or parameter census. Sage independently
recomputes every distinct finite-field count in the portable replay.

For the Nagao score, the
[comparison](../../artifacts/generated-results/elliptic-curves/bounded_prime_selector_comparison_v1.json)
gives the following Spearman association with the combined validation blocks.
The last column is the mean validation percentile of the **four** highest
scoring held-out curves; 0.5 is the centre of the holdout distribution.

| Family | Prefix | Full-population rho | Holdout rho | Selected holdout validation percentile |
|---|---:|---:|---:|---:|
| MW16 | 25 | 0.000 | -0.155 | 0.411 |
| MW16 | 64 | 0.083 | 0.016 | 0.311 |
| MW16 | 128 | 0.039 | 0.040 | 0.350 |
| MW16 | 256 | -0.022 | 0.045 | 0.506 |
| MW18 | 25 | -0.002 | 0.057 | 0.402 |
| MW18 | 64 | 0.005 | 0.167 | 0.647 |
| MW18 | 128 | -0.019 | 0.177 | 0.712 |
| MW18 | 256 | -0.078 | 0.163 | 0.560 |

The MW18 holdout improves through 128 and then weakens; the full-population
association does not support the same improvement. MW16 has little rank
association. The Mestre-log score likewise has no uniform monotone benefit;
all results, both individual validation blocks, selected IDs and coefficient
sizes are retained. These small deterministic comparisons do not establish
a rank-enrichment effect or select a winner for production.

The natural first-25-prime arm is **not** the old sparse 25-prime capped
Pareto selector. Its exact previously retained selections are also evaluated
as a baseline: their full-population validation percentiles average 0.615
(MW16) and 0.346 (MW18) under the Nagao validation score. The baseline is
not retuned and its original selections are preserved.

Point recovery is compared only where it was actually measured. Within each
fixed 15/16-curve measured subset, every prefix and both scores select four
curves: each has **4/4 masked direction recoveries in 48 completed charts**,
and **zero original unmasked quotient gains**. The
[joined replay/report](../../artifacts/generated-results/elliptic-curves/rank_jump_diagnostics_replay_v1.json)
records these exact endpoints, including the legacy baseline. Both outcomes
have zero variance, at opposite extremes. They cannot distinguish selector
quality. Unsearched population rows are never assigned zero recovery.

## 4. A construction that supplies the additional section

One fixed construction pairs retained `531/08234-12f61` with the lexically
first other cover at that anchor, `08234-0a9bf`. No pair enumeration or
Selmer calculation is involved. This is an instance of the standard V4
construction, not a claim of a novel rank-19 method.

For its two exact coprime, squarefree quadratics q1,q2, use the actual base

\[
 C:\quad u_1^2=q_1(t),\quad u_2^2=q_2(t).
\]

It is connected, of degree four over the old parameter and genus one. The
common anchor supplies a rational point. Parametrizing the first conic
produces a pointed quartic model of **C itself**, rather than just its
degree-two product quotient. Its opposite-anchor point P has reductions
of orders dividing 32,54,60 at good primes 23,41,47. Their gcd is two,
and exact computation gives 2P != O. Thus C has rank at least one; no
rank upper bound, generator saturation or point search is required.

Both covers supply exact lifted sections T1,T2. Their distinct quadratic
characters make the anti-invariant parts orthogonal. If G is the old
rank-17 Gram and q_i the two trace words, the inherited block is 4G,
the cross columns are 2Gq_i, both section heights are 16, and
`<T1,T2>=q1^T G q2`. The Schur complement is diag(6,6), with determinant
`4^17 * 948 * 36 > 0` for the full Gram. Consequently the displayed
subgroup has generic rank **at least 19**, supplying **one independent
direction beyond the retained MW18** after the further base change.
The same rank-28 anchor then has nine displayed directions remaining,
instead of ten; its rank has not increased.

The [construction certificate](../../artifacts/generated-results/elliptic-curves/point_supplied_mw19_diagnostic_v1.json)
also maps the predetermined base points 2P and 3P to exact fibres and both
supplied points, recording the arithmetic cost:

| Base point | Native t-height, bits | log2 j-height | Largest raw fibre coefficient, bits | Largest supplied point coordinate, bits |
|---|---:|---:|---:|---:|
| 2P | 370 | 6,806.2 | 5,990 | 3,013 |
| 3P | 882 | 19,201.8 | 12,179 | 6,101 |

The j-height is intrinsic. The other sizes refer to the exact retained raw
models and maps, without a claim of minimal presentation or a primitive base
generator. The sampled points are verified on their fibres; independence
on those two particular fibres is **not tested**. The achieved endpoint is
the independent generic section and positive-rank base, not a specialized
rank certificate. This supplies the missing point-producing implication,
but its height cost does not establish an advantage for a record search.

Small-point-conditioned twist construction and staged descent have explicit
precedents in [Watkins et al., §§3–6](https://pmb.centre-mersenne.org/item/10.5802/pmb.9.pdf).
Here the next useful construction comparison should constrain j-height and
point-coordinate cost alongside the new section; a large Selmer space would
not substitute for these witnesses.

## Reproduction

The portable check reads only committed inputs and witnesses:

```sh
sage -python elliptic-curves/cas/replay_rank_jump_diagnostics.sage
python3 -m unittest discover -s elliptic-curves/tests -p test_rank_jump_diagnostics.py
python3 -m unittest discover -s elliptic-curves/tests -p test_search_observability.py
python3 -m unittest discover -s elliptic-curves/tests -p test_rank_jump_audit_cli.py
```

All 25 focused tests pass. Replay checks 31 original independent subgroups,
372 search charts, 288 masked-point relations, 622 recovered MW18/public
relations, 162 non-endpoint translated minima and all 68,667 distinct local
counts. Search coverage still trusts the pinned enumeration worker; replay
does not independently enumerate the large boxes or prove global CVP minima.

Discovery commands, with per-chart/per-prime checkpoints in ignored local
directories, are:

```sh
python3 elliptic-curves/cas/audit_mw18_retained_visibility.py
sage -python elliptic-curves/cas/audit_mw18_translated_visibility.sage
sage -python elliptic-curves/cas/run_ordinary_masked_controls.sage freeze
sage -python elliptic-curves/cas/run_ordinary_masked_controls.sage run
sage -python elliptic-curves/cas/run_ordinary_masked_controls.sage audit
sage -python elliptic-curves/cas/audit_ordinary_masked_relations.sage
python3 elliptic-curves/cas/compare_bounded_prime_selectors.py freeze
python3 elliptic-curves/cas/compare_bounded_prime_selectors.py compute
python3 elliptic-curves/cas/compare_bounded_prime_selectors.py analyze
sage -python elliptic-curves/cas/construct_point_supplied_mw19_diagnostic.sage
```

`freeze` requires a new protocol path and refuses to replace an existing
protocol. The masked worker never reads the separate oracle files. Results
and source hashes are retained; the original 600-chart and height-population
certificates are not rewritten.

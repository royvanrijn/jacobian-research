# Elkies 2026 rank-17 paper: certified integration and search handoff

Date: 2026-08-31

Source: Noam D. Elkies, *An elliptic K3 surface X/Q(t) with Mordell-Weil
rank 17, I: Formulas for X and base changes of ranks 18 and 19*,
[arXiv:2608.25406v1](https://arxiv.org/abs/2608.25406), submitted 2026-08-26.

## Outcome

The paper is now integrated into the exact H3 theorem package and the active
search framework. The repository has:

1. the compact published `A(t),B(t)` model;
2. all seventeen published sections, stored as the displayed x-coordinates,
   `P1`'s y-coordinate, and the sixteen quadratic chord descriptions;
3. an exact Mobius and Weierstrass identification with q12/orbit5867;
4. exact positive controls at all four disclosed high-rank parameters, with
   quotient gains `8,9,10,11` beyond the generic rank 17;
5. a complete compact-`t`, height-10000, three-ensemble Nagao calibration;
6. a fail-closed residual 2-Selmer gate before expensive point search;
7. an explicit eighteenth section and rational parameter on the first conic;
8. both new sections on the paired cover, explicit `E0 -> t` maps, four
   independent `E0` generators, and a modular Mordell--Weil sieve.

The independent q8/q12 construction remains the canonical marked route proof.
The paper supplies the preferred arithmetic coordinate and the base changes.

## Published rank-17 data

The compact equation is

```text
y^2 = x^3 + A(t)*x + B(t),  deg(A,B,Delta)=(8,12,24).
```

The exact verifier reconstructs every ordinate from the published chords and
checks all seventeen Weierstrass identities. It also proves:

```text
height Gram determinant                    948
height-2 vector pairs                        0
unoriented height-4 pairs                 1311
published-basis -> pinned-R17 determinant    1
```

The coordinate matcher proves exact identities

```text
t = (a*u+b)/(c*u+d),
A_q12(u) = s^4 (c*u+d)^8  A(t),
B_q12(u) = s^6 (c*u+d)^12 B(t)
```

over `QQ`, with trivial twist. The compact published `t` chart is therefore
the default specialization and Nagao chart. The raw q12 coordinate remains a
construction regression; its images of the four disclosed fibres have
221--234-bit numerators.

## Calibration anchors

Use these before assessing a search change:

```text
certified rank >= 25: t = -2/377      quotient gain >= 8
certified rank >= 26: t = -308/251   quotient gain >= 9
certified rank >= 27: t = 2456/135   quotient gain >= 10
certified rank >= 28: t = -9529/5471 quotient gain >= 11
```

All four exact fibres equal the corresponding public minimal models in
Dujella's rank-record table. The imported exact public point sets have lengths
25, 26, 27 and 28. Finite-reduction pivoting keeps the generic seventeen
first, selects public complements of dimensions `8,9,10,11`, and certifies
each combined list in one matrix. These are unconditional rank lower bounds,
not exact-rank results. The rank-28 fibre is a positive control, not permission
for an indefinite four-point search. ICARM curves 302/351/356 are no longer
needed to infer the family.

## Compact-t scoring calibration

The accepted scorer partitions the 102 primes from 19 through 599
round-robin into three disjoint ensembles of 34. Each local Nagao contribution
is centered and population-standardized over the good fibres of
`P^1(F_p)`; singular fibres receive the prime mean and remain counted. The
primary ranking key is the weakest normalized block.

The complete height-10000 scan contains 121,589,944 primitive projective
parameters. The four positive-control ranks are:

```text
t=-2/377      rank  54,624   fraction 0.0449%
t=-308/251    rank 593,936   fraction 0.4885%
t=2456/135    rank 422,873   fraction 0.3478%
t=-9529/5471  rank  55,387   fraction 0.0456%
```

All pass the declared top-one-percent gate. The ranking remains heuristic; it
only decides which candidates deserve descent.

## Residual 2-Selmer gate

For a candidate fibre the required object is

```text
Sel_2(E_t) / <P1,...,P17>.
```

Rank 32 requires residual dimension at least 15. A completed unconditional
2-descent with residual dimension below 15 is therefore an exact rejection.
Only a result at least 15 authorizes two-cover solving or expensive point
search on the same minimal model. The BNF-free Kummer signature, norm-one
cubic elements, incomplete relation ledgers, `K(S,2)` envelopes and candidate
local classes do not pass this gate.

The first exact-backend attempt used PARI `ellrank` through Sage on the public
rank-28 fibre, with all 28 certified points supplied. It reached the strict
300-second limit at 230,338,560 bytes peak observed RSS without returning a
Selmer dimension. An independent eclib invariant-quartic descent, with
`selmer_only=True` and both point-search bounds zero, also reached 300 seconds
without a result at 232,099,840 bytes peak observed RSS. Both pinned results are
therefore
`INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN`: it is not an upper bound and it
does not authorize a point search. The rank-28 control already supplies eleven
quotient directions and would need four more for rank 32.

The hidden first cost in PARI was factorization of the monic 2-division cubic
discriminant. That discriminant is now completely factored as

```text
2^23 3^6 5^6 7^4 11^2 13^4 17^5 19^3
* 48463 * 20650099
* 315574902691581877528345013999136728634663121
* 376018840263193489397987439236873583997122096511452343225772113000611087671413.
```

Every factor is independently proved prime and the product is checked against
the closed cubic-discriminant formula. A factor-supplied, certified number
field then computes the generic-seventeen Kummer rows at every bad odd prime,
at 2, and at the real places. All thirteen local blocks complete; their
concatenated coordinate matrix has rank 15. This closes the factorization and
known-image layer only: it does not enumerate ambient `K(S,2)`, compute the
`S`-class quotient, classify all locally soluble ambient classes, or give a
Selmer upper bound.

Supplying that proved factor table to PARI removes the repeated factorization.
The first 8 GB-stack run entered the descent immediately and reached a strict
600-second wall limit at 5,698,514,944 bytes peak observed RSS, without a
Selmer dimension. Its status therefore remains
`INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN`. A second run with the same
8 GB PARI stack reached its strict 1,800-second limit at 6,040,723,456 bytes
peak observed RSS and likewise returned no dimension. This demonstrates that
factorization is no longer the active delay, but does not complete the
class-group/local-solubility layer or authorize search.

The class-group layer is now independently isolated. A stage-aware PARI
worker supplies the same proved factor support to `nfinit`, completes
unconditional maximal-order `nfcertify`, and then runs `bnfinit(...,0)` before
the one-sided `bnfcertify(...,1)` class-quotient certificate. With
`c1=0.01,c2=4,nrpid=20`, the pinned 120-second run stops inside `bnfinit` at
265,261,056 bytes peak observed RSS; PARI's diagnostic stream shows a
243-ideal factor base followed by a random-relation plateau at 153 requested
relations. Since `bnfcertify` is never reached, this is implementation
evidence only.

An exact `polredabs` variant replaces the large translated cubic by a
depressed cubic and records the original-generator map `theta=-3*x+1`. The
polynomial discriminant drops by `3^6` and the defining-order index by 27;
factor-supplied `nfinit` and `nfcertify` still prove the same maximal field.
In the matched 120-second envelope this cleaner model reaches the same
243-ideal setup and 153-request random-relation plateau, timing out at
264,839,168 bytes peak observed RSS. It is therefore an exact coordinate
improvement but not a class-group, Selmer, or rank bound.

The BNF-free collector now also has an exact rank-28 preset and reuses the
proved factor hints in collection, canonical-row augmentation, and audit. A
factor-base-1000 paired-special-ideal pilot retains 10,288 sampled algebraic
integers; its augmented ledger adds and verifies 172 rational principal rows.
No noncanonical relation closes. The bounded model has 327 ideal columns, 26
`S` columns, exact relation rank 172, and displayed quotient dimension 141,
but its factor base is below the Bach/ERH generation bound 1,202,640. The
pinned classification is therefore `UNCERTIFIED_FACTOR_BASE`. This makes the
short-vector relation route reproducible and gives it a scale-specific stop
result; it does not compute an `S`-class quotient, ambient `K(S,2)`, local
solubility, Selmer dimension, or rank bound.

The local layer now has an exact rank-28 positive control. A fresh fixture
recomputes all 53 bad-place coordinates for both the generic seventeen and
the certified public complement of eleven. The generic local-signature rank
is 15; after adjoining all eleven globally independent directions it remains
15, and the incremental rank is zero in every local block. This is direct
calibration evidence that known Kummer signatures do not measure the
Mordell--Weil quotient and cannot replace the ambient Selmer calculation.

An independent coverage audit proves that the generic images span the full
local Kummer image at four of the eleven odd bad primes (`3`, `19`,
`20650099`, and
`315574902691581877528345013999136728634663121`) and at infinity. Seven odd
places and the two-adic place remain unresolved. A bounded exact norm-one
generator then produces 49 candidate cubic classes and their two-cover
quadrics. The resumable local supervisor tests the first twelve covers at
seven odd primes in 84 owned cover/place workers: 60 return smooth-reduction
local points, 19 singular lift trees hit the state cap, and five workers time
out. No local obstruction is found. This pilot validates the one-sided local
witness path while leaving 24 cases inconclusive; it proves neither
everywhere-local solubility nor Selmer membership and authorizes no search.

The public complement now supplies the matching genuine positive control for
this layer. For each of its eleven points `Q`, an exact builder forms
`alpha=X(Q)-theta`, verifies `Norm(alpha)=Z(Q)^2`, and materializes the two-cover
with rational point `[1:0:0:1]`. The generic cover builder and local checker
replay all eleven witnesses. Combined with the finite-reduction certificate
that makes these point classes independent modulo the generic seventeen, this
proves a residual 2-Selmer lower bound of 11 at the rank-28 fibre. The result
shows that the explicit-cover path sees every known exceptional direction even
though the 53-coordinate signature span sees none of their quotient gain. It
is still a lower bound: the four additional directions needed for rank 32 and
the complete ambient Selmer upper bound remain open.

The external complete path is
`build_elkies_2026_rank28_relative_descent_magma.py`. Before emitting code it
replays the exact generic-17 and generic-plus-public-11 finite-reduction
certificates. The generated job calls unconditional
`TwoSelmerGroup(E : Bound := -1)` at basis level, computes
both `dim Sel_2(E)-17` and `dim Sel_2(E)-28`, and exits if the former is below
15 (equivalently if the latter is below four). Only the passing branch may
call relative `TwoDescent`, with `RemoveGens` equal to the generic seventeen
plus the public complement eleven; it therefore materializes only genuinely
unexplained covers. There is no point-search primitive in the job. Its parser
rejects partial logs and emits the same model-bound gate schema. Magma is
unavailable on the current host, so no external transcript or Selmer dimension
is claimed.

## First rank-18 cover

The published conic is

```text
u^2 = 4225*t^2 + 38636*t + 289444.
```

Branch reductions isolate a single nonzero mod-2 trace coset. A shortest trace
has height 10. Exact chord recovery factors the residual discriminant as

```text
(127170526080*h(t))^2 * (4225*t^2+38636*t+289444).
```

The stored section has the form

```text
x = x0(t) + x1(t)*u,
y = y0(t) + y1(t)*u,
```

and passes both coefficient identities modulo the conic equation and the
Galois trace identity. The conic parameterization is

```text
t = (289444-r^2)/(130*r-38636),
u = 65*t+r.
```

The paper's anti-invariance lemma then proves generic rank at least 18 over
`QQ(r)`. `search_elkies_2026_rank18_conic_nagao.py` pulls the existing local
R17 tables back along this rational map and sieves primitive rational `r`.

## Paired rank-19 cover

The second cover is

```text
u2^2 = 54756*t^2 - 3269604*t + 22473889.
```

Its second height-10 trace and polynomial section are also recovered exactly.
After parameterizing the first conic, the fibre product becomes

```text
v^2 = 54756*r^4 + 425048520*r^3 + 221786712628*r^2
      - 348786049427920*r + 74698868489239696.
```

A rational plane-cubic transformation gives the paper's curve

```text
E0: y^2 = x^3 + 1029367969*x^2 - 42900734074705920*x.
```

The exact data file records rational functions `r(x,y)` and `v(x,y)`, followed
by

```text
t  = (289444-r^2)/(130*r-38636),
u1 = 65*t+r,
u2 = v/(130*r-38636).
```

All identities reduce to zero modulo the `E0` equation. Four displayed points
are unconditionally independent, proving `rank(E0(Q)) >= 4`; the paper states
exact rank 4. The two anti-invariant characters give generic rank at least 19.

`search_elkies_2026_E0_mw_nagao.py` reduces the four-generator lattice once
per prime, enumerates its finite image in at most `#E0(F_p)` steps, caches the
`E0 -> t` Nagao symbol, and scans a bounded coefficient box. Finalists include
the exact rational `E0` point and exact `t` value.

## Active search order

1. Keep the four exact fibres in every score calibration report.
2. Search compact `t=a/b` with three or more disjoint prime ensembles and rank
   by weakest-block performance.
3. Compute the actual residual 2-Selmer quotient for a survivor. Reject it
   exactly when the residual dimension is below 15.
4. Solve two-covers or run expensive point search only after the same minimal
   curve passes that gate.
5. Promote rank 32 only after fifteen certified quotient directions and one
   finite-reduction independence matrix of rank 32.

The first-conic and `E0(Q)` lattice sieves remain supporting base-change
routes. They do not bypass the residual gate.

Nagao scores, bounded coefficient boxes, and absence of points are experiments,
not rank bounds.

## Parked paths

The following are provenance or regression paths, not current priorities:

- q8/q12 coefficient discovery after the exact endpoint;
- q12/orbit4484 equation lifting;
- fixed-corridor reverse lifts, q323, changed-zero reranking, and compiler
  optimization without a direct specialization use;
- ICARM fingerprint fitting as a way to infer the now-published family;
- ungated raw `ratpoints`, slope-box, and two-cover point searches; the direct
  q12 ratpoints, affine-chart, eclib-search, and slope-slice commands now refuse
  to start without a same-parameter, same-minimal-model passing gate;
- the old low-complexity x-ansatz parameter/point search, now hard parked in
  favor of compact-`t` calibrated scoring followed by descent;
- further searching of the rank-28 calibration fibre merely to rediscover its
  known 28 independent points.

## Reproduction

```bash
SAGE=/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python

$SAGE elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
$SAGE elkies-k3/scripts/match_h92_q12o5867_to_elkies_2026_qq.sage
python3 elliptic-curves/scripts/verify_elkies_2026_high_rank_calibrations.py
python3 elkies-k3/scripts/calibrate_elkies_2026_positive_controls_nagao.py
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --timeout 300 --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend eclib --timeout 300 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_eclib_v1.json \
  --overwrite
python3 elliptic-curves/scripts/specialize_q12o5867_candidate.py \
  --a -9529 --b 5471 --overwrite
python3 elliptic-curves/cas/build_elkies_2026_rank28_bad_place_ledger.py \
  --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend pari-factored --timeout 600 \
  --pari-stack-bytes 8000000000 --rss-limit-bytes 12000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_pari_factored_8g_v1.json \
  --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend pari-factored --timeout 1800 \
  --pari-stack-bytes 8000000000 --rss-limit-bytes 12000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_pari_factored_8g_30min_v1.json \
  --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_s_class_pari.py \
  --timeout 120 --c1 0.01 --c2 4 --nrpid 20 \
  --pari-stack-bytes 2000000000 --rss-limit-bytes 3000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_s_class_pari_v1.json \
  --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_s_class_pari.py \
  --field-model polredabs --timeout 120 --c1 0.01 --c2 4 --nrpid 20 \
  --pari-stack-bytes 2000000000 --rss-limit-bytes 3000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_s_class_pari_polredabs_v1.json \
  --overwrite
mkdir -p artifacts/local/elliptic-curves/elkies-rank28-bnf-free
$SAGE elliptic-curves/cas/run_fermigier_rank20_minkowski_specialq.py \
  --elkies-rank28 --factor-base-bound 1000 \
  --special-q-min 1009 --special-q-max 5000 --max-special-q 20 \
  --special-ideal-mode cycle-pairs --pair-cycle-length 20 \
  --trial-prime-bound 1000 --lattice-combination-bound 2 \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_v1.json
$SAGE elliptic-curves/cas/augment_bnf_free_canonical_principal_relations.py \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_v1.json \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_canonical_v1.json
$SAGE elliptic-curves/cas/audit_bnf_free_s_class_quotient.py \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_canonical_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_bnf_free_s_class_pilot_v1.json
$SAGE elliptic-curves/cas/build_elkies_2026_rank28_local_coverage.py \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --overwrite
$SAGE elliptic-curves/cas/audit_bnf_free_local_kummer_coverage.py \
  --signature-map artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_coverage_v1.json
$SAGE elliptic-curves/cas/generate_q12o5867_norm_one_cover_candidates.py \
  --signature-map artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --coefficient-bound 2 \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_v1.json
$SAGE elliptic-curves/cas/build_bnf_free_two_covers.py \
  --candidates artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_v1.json \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_covers_v1.json
python3 elliptic-curves/cas/run_bnf_free_two_cover_local_supervisor.py \
  --covers artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_covers_v1.json \
  --primes 3,5,7,11,13,17,19 --max-covers 12 \
  --max-enumeration-prime 19 --max-lift-precision 6 --max-lift-states 5000 \
  --timeout-per-place 2 \
  --cache-dir artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_local_pilot_blocks_v1 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_norm_one_local_pilot12_v1.json
python3 elliptic-curves/cas/build_elkies_2026_rank28_public_selmer_controls.py \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_selmer_candidates_v1.json \
  --overwrite
$SAGE elliptic-curves/cas/build_bnf_free_two_covers.py \
  --candidates artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_selmer_candidates_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_two_cover_controls_v1.json
$SAGE elliptic-curves/cas/audit_bnf_free_two_cover_reduction.py \
  --covers artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_two_cover_controls_v1.json \
  --primes 2,3,5,7,11,13,17,19,48463,20650099,315574902691581877528345013999136728634663121,376018840263193489397987439236873583997122096511452343225772113000611087671413 \
  --max-enumeration-prime 2 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_global_cover_witness_audit_v1.json
python3 elliptic-curves/cas/build_elkies_2026_rank28_relative_descent_magma.py \
  --output artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.m
# Run the generated program with an unconditional Magma installation, then:
python3 elliptic-curves/cas/parse_elkies_2026_rank28_relative_descent.py \
  --program artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.m \
  --log artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.log \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_magma_v1.json
$SAGE elkies-k3/scripts/verify_elkies_2026_rank18_first_cover.sage
$SAGE elkies-k3/scripts/verify_elkies_2026_rank19_paired_cover.sage

python3 elkies-k3/scripts/search_elkies_2026_rank18_conic_nagao.py \
  --numerator-bound 1000 --denominator-bound 1000 \
  --output artifacts/local/elkies-k3/elkies-r18-conic-nagao-h1000.json

python3 elkies-k3/scripts/search_elkies_2026_E0_mw_nagao.py \
  --coefficient-bound 6 --primes 11-97 \
  --output artifacts/local/elkies-k3/elkies-E0-mw-nagao-b6.json
```

The model, coordinate, control, and cover verifiers are exact replays. The
Nagao command is a complete bounded heuristic ranking. The two local
residual-descent commands are exact backend attempts, but their pinned timeouts
are incomplete and authorize no search. The Magma builder/parser is an exact,
fail-closed external path, not a completed computation on this host. The last
two commands are bounded supporting heuristic searches.

<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-ELKIES-2026-HIGH-RANK-CALIBRATIONS 345b9fb977057133 -->
<!-- status-consumer: EC-K3-ELKIES-2026-NAGAO-POSITIVE-CONTROL f99c98cdb6b8cd7d -->
<!-- status-consumer: EC-K3-ELKIES-2026-R28-S-CLASS-PILOT a791713dc40f7caf -->
<!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE 56298144d268ab70 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R18-COVER 6b4ee5bbc1afc01e -->
<!-- status-consumer: EC-K3-ELKIES-2026-R19-PAIRED f1e135d2ba803e80 -->

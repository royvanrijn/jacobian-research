# Prospective R17 CRT rank-jump experiment

Date: 2026-09-04  
Status: complete frozen first experiment; detector-limited; replacement half-lattice protocol frozen and running; rank-32 promotion forbidden

<!-- status-consumer: EC-K3-R17-074D9-PROSPECTIVE-CRT-LOCAL-STABILITY 0edaaa6f05041634 -->
<!-- status-consumer: EC-K3-R17-074D9-PROSPECTIVE-CRT-ESCAPE-EXPERIMENT 021a952efb9ea0f4 -->

## Result

The original `p^3` premise fails, but a pre-outcome local-only refinement does
produce two empirically stable finite-sample cylinders.  The complete frozen
experiment then finds no Mordell--Weil escape in any cohort under the uniform
direct rational-point bound:

| cohort | intended condition | completed | certified MW17 escapes | certified extra directions |
|---|---|---:|---:|---:|
| A | full 356 fingerprint | 256 | 0 | 0 |
| B | full 385 fingerprint | 256 | 0 | 0 |
| C | matched ordinary | 512 | 0 | 0 |
| D | `p=2` only | 512 | 0 | 0 |
| E | four odd places only | 512 | 0 | 0 |
| F | random equal-codimension residues | 512 | 0 | 0 |

There were no timeouts, backend failures, or structurally invalid rows in the
frozen 2,560-row experiment.  A post-experiment sensitivity audit then applied
the unchanged search call to the known +12 fibres 356 and 385 and rediscovered
no escape on either.  Thus the experiment is detector-limited and supplies no evidence that
full multi-place alignment, `p=2` alone, or the odd conditions alone enriches
the probability of an escape at this search bound.  It does **not** show that
the fibres have rank 17, that the true escape probabilities are equal, or that
the local fingerprints lack predictive value at a deeper bound.

For the primary pooled comparison, full A+B is `0/512` and ordinary C is
`0/512`: risk difference `0`, risk ratio and odds ratio undefined, and
two-sided Fisher exact `p=1`.  Each arm's two-sided 95% Clopper--Pearson upper
limit is `0.007179`; the corresponding conservative Cartesian bounds for the
risk difference are `[-0.007179, 0.007179]`.  Each individual 256-row anchor
arm has zero events and upper limit `0.014306`.  Because the response has zero
variance, no outcome predictor is fit.

## Phase 1: target-blind premise audit

The audit specializes the exact `norm12-orbit-074d9` R17 equation and verifies
all seventeen generic sections before doing local arithmetic.  Its only
structural rejection conditions are exact singularity and failure of an exact
section identity or local-reduction computation.  Each selected rational
prime is handled independently in a `p`-maximal 2-division order.  The
comparison payload includes reduction type, Tamagawa and section-component
data, the ambient local squareclass dimension, the actual known-MW17 local
image presentation, and its source kernel.  Odd-prime raw uniformizer
coordinates are retained diagnostically, while cross-fibre matching uses the
source-kernel invariant because a PARI uniformizer may change by a nonsquare.

On 64 salted samples from each original `p^3` cylinder, every specialization
was valid, but the intended local conditions survived as follows:

| anchor | `p=2` | `p=13` | other three intended places | matched-place totals |
|---|---:|---:|---:|---|
| 356 | 0/64 | 1/64 | 64/64 each | 63 rows at 3/5; 1 row at 4/5 |
| 385 | 0/64 | 32/64 | 64/64 each | 32 rows at 3/5; 32 rows at 4/5 |

In particular, `t mod p^3` does not force the tested local Kummer fingerprint.
No Mordell--Weil search result was read during refinement.  For each
anchor/place the script selected the first exponent `k=3,...,20` passing 16
salted discovery draws and then 64 disjoint confirmation draws.  All selected
places passed `16/16` and `64/64`:

| anchor | prime-power exponents in prime order | frozen exact class |
|---|---|---|
| 356 | `2^15,13^5,37^3,53^3,71^3` | `t = 14503094594046060605284928 + 32837835553385972844560384*n` |
| 385 | `2^18,13^4,37^3,53^3,67^3` | `t = -4622816536543600156190144 + 16981283803546391486922752*n` |

The frozen cylinder-definition hash is
`500dc6931c5aeaf3d6d9982bb994286d7aee36e7c87b9e414e8b7e0ef8aef15c`.
This is an exact commitment to empirically selected residue classes, not a
proof of constancy on either complete `p`-adic cylinder.

## Phase 2: matched frozen cohorts

The manifest uses the salt `r17-prospective-crt-frozen-cohorts-v1` and freezes
256 rows per anchor lane before any point-search outcome is opened.  A and B
contain the two full classes.  C avoids every individual target residue.  D
preserves only the refined 2-adic residue and avoids every target odd residue.
E preserves all four target odd residues and excludes the refined 2-adic
residue.  F uses salted random residues at exactly the same prime powers as
the corresponding full anchor.

All parameters are integers in the common 111-bit absolute-height shell
`[2^110,2^111-1]`, with denominator one and paired salted target heights.
Exact parameter and exact `j` duplicates are excluded.  Across every cohort,
the numerator bit ranges are 1904--1945 for `A(t)`, at most 2851--2915 for
`B(t)`, and at most 5696--5818 for the discriminant core.  No Nagao score,
point, public rank, cover split, or later outcome enters selection.

The complete candidate-list commitment is
`5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb`.
The original manifest remains unopened; outcomes live in separate ledgers.

## Phase 3: pre-search arithmetic panel

All 2,560 frozen rows pass exact specialization, all seventeen section
identities, and the predeclared arithmetic computation at
`2,13,37,53,67,71`.  The compressed artifact stores the actual local image
rows/source kernels and component data, not only dimensions, together with
cumulative intersection and leave-one-place-out matrices and the fixed
three-block Nagao comparison feature.

The intended fingerprint survival is:

| cohort | exact matched-place distribution |
|---|---|
| A | `5:256` |
| B | `5:256` |
| C | `0:493, 1:19` |
| D | `1:495, 2:17` |
| E | `4:512` |
| F | `0:512` |

So the refined full and odd-only congruences do reproduce the intended local
data on every frozen row, while D has 17 accidental odd-place matches and C
has 19 accidental single-place matches.  These are retained without
rebalancing.

The full known-MW17 stacked-localization ranks range from 9--11 in A, 7--9 in
B, 4--11 in C, 3--11 in D, 7--11 in E, and 8--10 in F.  These are dimensions
of images of the known specialized MW17 subgroup; they are **not** residual
Selmer dimensions or Mordell--Weil upper bounds.

The all-pairs Hilbert/Tate panel was uniformly deferred as infeasible for this
moderate 2,560-fibre pass and was not used to tune a cylinder.  The monotone
sieve authorized all 2,560 bounded searches but produced zero complete
2-Selmer groups and zero finite proved residual upper bounds.  Every absent
upper bound remains `null`, never an inferred value.

## Phase 4: identical bounded search

The manifest's original v1 backend was Sage/eclib at height 12.  A single
predeclared canary spent its full 300 seconds inside
`mwrank_EllipticCurve` initialization, never reached the search call, and
returned no point.  This is recorded as an operational feasibility failure,
not a bounded miss.

Before any point-search call completed or returned a point, a uniform v2
amendment was frozen without changing candidates, cohorts, matching, features,
or comparisons.  Its protocol hash is
`63d6b9e83f52bc7208b9057298e05941dfcedc85d53f5681186c953498947d4b`.
For the deterministic exact integral 2-minimal model

```text
y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
```

PARI `hyperellratpoints` searches the exact completed-square cubic

```text
Y^2 = 4*x^3 + (a1^2+4*a2)*x^2 + (2*a1*a3+4*a4)*x + (a3^2+4*a6)
```

at rational `x` numerator/denominator height 10,000, with 30 wall seconds
including setup, 8 GB memory, zero retries, no adaptive stopping, and the same
backend for every row.  Any returned point must transport exactly to the
original specialization and make the finite mod-2 reduction certificate for
MW17 plus all counted extras attain full column rank.  No returned candidates
occurred; all 2,560 rows are genuine bounded-protocol misses.

After the ledger was frozen, the identical completed-square call was applied
to historical positive controls 356 and 385.  Both have twelve exactly
certified known directions beyond MW17, but the height-10,000 call returned
zero points on both.  This post-experiment diagnostic changes no candidate,
protocol, or contrast.  It shows that the detector is not sensitive to even
the motivating positive controls at the declared bound; accordingly the
2,560 zeroes cannot be read as strong evidence against the local-arithmetic
hypothesis.

## Phase 5: interpretation

The first arrow in the proposed chain is supported only in the finite sampled
sense after substantial exponent refinement:

```text
refined parameter congruences
        -> intended local MW17/Kummer fingerprints on all frozen rows.
```

The second arrow remains open because no complete Selmer group or finite
proved residual upper bound was obtained.  The third arrow is unresolved
because the bounded search has zero events in every cohort and fails its
post-experiment positive-control sensitivity check.  Consequently the
`random -> p=2 -> odd -> full` ablation has no observable ordering at this
bound, and the first experiment does not authorize a large rank-32 campaign.

Any continuation must be a separately specified, separately hashed deeper
search or second replication.  The present candidates and outcomes must not
be extended or rebalanced in place.  A sensible next gate would first validate
that a stronger uniform search has nonzero sensitivity on predeclared positive
controls, then freeze its application to this existing ledger or a new salted
replication.

The fixture-blind half-lattice replay passes this detector-sensitivity gate.
The fixed generic-deepest 43 recover exact quotient gains 12 and 3 on the +12
controls 356 and 385; the generic/specialized unions recover 12 and 4.  A new
protocol therefore applies that detector to this same committed cohort without
altering the first experiment or its zero-event interpretation.  This is not a
rank-32 promotion gate: sensitivity to `jump>=1` does not show that the detector
score tracks jump size or extreme-tail incidence.

## Phase 6: frozen half-lattice detector replacement

The replacement protocol is separate from the original direct-search ledger
and was frozen before any half-lattice outcome on the 2,560 prospective fibres
was opened.  It retains candidate-list hash
`5df03637d4db0baa95cb9e5f697fe35e5e897838676b6370c0e08bdae5aa9aeb`
and defines:

- Stage A: the fixed 43 exact generic R17 norm-12 classes, with a specialized
  shortest representative on each fibre;
- one identical minimize/reduce/search call per class at height 100,000 and a
  15-second whole-cover timeout;
- Stage B only after an exactly certified Stage-A quotient direction, followed
  by the full specialized ranking and the incremental generic/specialized
  top-43 union;
- exact transport to the original specialization and a full finite-reduction
  independence certificate for every counted point;
- pooled A+B versus C Stage-A detector yield per scheduled fibre as the primary
  contrast, with Stage B reported only as conditional incremental recovery.

The protocol-definition hash is
`9584174de7625031e5f95ce73d0117a9caf8341d91063061ea672f2e4e36e521`.
The compact protocol, deterministic chunk runner, and frozen analyzer are
generated by the commands below.  The run is checkpointed under
`artifacts/local/`; no cohort comparison is made before all scheduled rows are
merged.

### Rank-32 promotion boundary

The post-freeze promotion audit is deliberately restrictive and applies
regardless of the eventual v3 result.  Protocol v3's primary event is binary:
one exactly certified Stage-A direction is enough.  It may therefore compare
bounded detector-visible escape yield, but no v3 success, effect size, or
`p`-value may promote a fibre for a rank-32 campaign.

The calibration does not fill that gap.  Development scores are `8,9,10,9`
on displayed jumps `8,9,10,11`, including an order reversal.  The three sealed
controls all have displayed jump 12 but score `10,12,3`, so the independent
panel has no between-stratum magnitude information.  Its largest displayed
jump is 12, while rank 32 requires 15 independent directions beyond `MW17`.

A successor promotion rule must freeze the exact certified quotient gain as an
ordinal score, validate directional score--jump association on an independent
multi-stratum panel, and predeclare and pass an upper-tail enrichment endpoint.
The separate completed residual 2-Selmer gate on the same minimal curve remains
mandatory before expensive follow-up.  Finding fifteen exactly certified
directions is different: that directly proves rank at least 32 and needs no
heuristic promotion.

The same guard pins the repository-wide chart-order policy.  The pointed
quartics are birational search charts, not nontrivial 2-coverings.  Frozen
fields named `depth`, `deepest`, `old-deep-43`, or quotient Hamming `weight`
are search-order metadata for their exact recorded bases; they are not
arithmetic filtrations or Selmer data.  Any lattice or finite-index
enlargement, generator/basis change, height-form change, or quotient-basis
change invalidates the order.  A successor must recompute chart identities,
representatives, scores, and order, bind them to the new state fingerprint,
and revalidate any efficiency claim on blinded controls.  No miss can imply
point absence, local insolubility of a 2-covering, Selmer structure,
Mordell--Weil saturation, or a rank upper bound.

It also repairs the prospective denominator without changing the frozen v3
search.  The original scheduled-row analyzer remains at its protocol-pinned
hash for reproduction, but it is not an inferential endpoint.  The restrictive
v4 analyzer uses only complete Stage-A rows and authorizes a contrast only
when each distinct censor-status proportion is exactly balanced between its
arms.  If that gate fails, all effect estimates and the Fisher p-value are
null; exact events observed on censored rows remain descriptive lower bounds.

<!-- status-consumer: EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE 9a1f080523c9ecae -->

## Artifacts and commitments

| artifact | SHA-256 |
|---|---|
| `elkies-k3-r17-prospective-crt-local-stability-v1.json` | `bb573fbe66bede1625afc57e777ec37a41e58a560162a4207359968f745f173f` |
| `elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json` | `7e8c43a6f67eac96dd9dede333f94e0cce139fa685b421f83ad7e4d69c1a75d4` |
| `elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz` | `709ac16370fd5ced34857068102d57d634f375ad2a4ab98c73c221f0baf48d6b` |
| `elkies-k3-r17-prospective-crt-search-protocol-v2.json` | `1863d2b7b0ec31127cbe8b93beb203ff802c874b5d01dabc12bc58cac6b870d0` |
| `elkies-k3-r17-prospective-crt-point-search-ledger-v2.json` | `f99bdf697e9d9781237011eaaabbaf35cbb55ab9fe58b8ed2809c9210da7277b` |
| `elkies-k3-r17-prospective-crt-search-sensitivity-v1.json` | `9787d6010c8384b7ce7f13915345b03cff30c87bdc7fea64b3c32861036a7a01` |
| `elkies-k3-r17-prospective-crt-statistical-analysis-v1.json` | `e947ecc5ea0210f329afdd4d856f1e1c20d1e964a54c2aeb2e10f2818903fd8e` |
| `elkies-k3-r17-prospective-crt-half-lattice-protocol-v3.json` | `a402b1a286dd72ad579c753315a55309a92d03886f60ac8fee84e434119da626` |
| `elkies-k3-r17-prospective-crt-half-lattice-promotion-gate-v1.json` | `68350e6d099027e2817ee2969825fb57d6f215f713b61758fd6bdd78a4213812` |

The large arithmetic artifact stores every exact local matrix.  The point
ledger stores every frozen outcome and retains separate timeout/backend
fields even though their counts are zero.

## Reproduction

```bash
# Phase 1: exact p^3 audit and local-only refinement
sage -python elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage

# Phase 2: freeze the target-blind matched cohorts
python3 elkies-k3/scripts/build_r17_prospective_crt_cohorts.py

# Phase 3: eight deterministic checkpoints and canonical merge
for i in $(seq 0 7); do
  sage -python elkies-k3/scripts/run_r17_prospective_crt_arithmetic_features.sage \
    --chunk-index "$i" --chunk-count 8
done
sage -python elkies-k3/scripts/run_r17_prospective_crt_arithmetic_features.sage \
  --merge --chunk-count 8

# Freeze the pre-outcome feasibility amendment
python3 elkies-k3/scripts/build_r17_prospective_crt_search_protocol.py

# Phase 4: 32 deterministic checkpoints and canonical merge
for i in $(seq 0 31); do
  sage -python elkies-k3/scripts/run_r17_prospective_crt_direct_point_search.sage \
    --chunk-index "$i" --chunk-count 32
done
sage -python elkies-k3/scripts/run_r17_prospective_crt_direct_point_search.sage \
  --merge --chunk-count 32

# Phase 5 and regression checks
sage -python elkies-k3/scripts/audit_r17_prospective_crt_search_sensitivity.sage
python3 elkies-k3/scripts/analyze_r17_prospective_crt_experiment.py
python3 -m unittest elliptic-curves/tests/test_r17_prospective_crt_experiment.py

# Freeze and run the separately gated replacement detector
sage -python elkies-k3/scripts/build_r17_prospective_crt_half_lattice_protocol.sage --check
for i in $(seq 0 31); do
  sage -python elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage \
    --chunk-index "$i" --chunk-count 32
done
sage -python elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage \
  --merge --chunk-count 32
python3 elkies-k3/scripts/analyze_r17_prospective_crt_half_lattice_experiment.py
python3 -m unittest elliptic-curves/tests/test_r17_prospective_crt_half_lattice_protocol.py

# Enforce the separate rank-32 interpretation boundary
python3 elkies-k3/scripts/build_r17_prospective_crt_half_lattice_promotion_gate.py --check
python3 -m unittest elliptic-curves/tests/test_r17_prospective_crt_half_lattice_promotion_gate.py
```

For immutable replay checks of the deterministic commitments:

```bash
sage -python elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage --check
python3 elkies-k3/scripts/build_r17_prospective_crt_cohorts.py --check
python3 elkies-k3/scripts/build_r17_prospective_crt_search_protocol.py --check
sage -python elkies-k3/scripts/audit_r17_prospective_crt_search_sensitivity.sage --check
python3 elkies-k3/scripts/analyze_r17_prospective_crt_experiment.py --check
```

The Phase-1 artifact deliberately retains the upstream fingerprint file's
whole-file hash at freeze time,
`54548e6b7110d0b53ae3bd86a97bbd06fd1159836d19c3fc3ad4e23b77320fbc`.
That upstream artifact has since received additive diagnostics.  Phase-1
`--check` accepts this provenance drift only after recomputing and confirming
that every mathematical payload field is identical; it does not rewrite or
refreeze the original artifact.

## Claim boundary

- No complete Selmer group was computed.
- No incomplete descent datum is called an upper bound.
- No Selmer dimension is identified with Mordell--Weil rank.
- No bounded miss is interpreted as rank 17.
- No local finite-sample stability is promoted to a cylinder theorem.
- No outcome was used to refine a residue condition or rebalance a cohort.
- The 69 public R17 fibres are not used as a population denominator.
- Zero events at the declared bound are evidence only about that bounded
  protocol; its failure on both known +12 controls makes the result
  detector-limited, not a proof against the underlying arithmetic hypothesis.
- The replacement protocol's binary `jump>=1` endpoint is a detector-yield
  response only.  Without independent magnitude and upper-tail validation it
  cannot promote a rank-32 candidate, even if its cohort contrast is positive.

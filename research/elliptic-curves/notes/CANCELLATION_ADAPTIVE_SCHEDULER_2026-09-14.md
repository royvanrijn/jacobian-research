# Adaptive residue-box search with complete certification costs

This is a new, user-authorized next-direction experiment following the failed
[41-curve fitted-selector validation](FINITE_CANCELLATION_VALIDATION_2026-09-14.md).
Its objective is independently certified directions per complete CPU, including
independent rank verification. A residue-membership score cannot pass that gate.
The independently held-out successor recovers **22/24 directions under both
arms**, reducing complete CPU from **325.34 to 204.30 seconds**. Its declared
gate passes. The subsequent M27/M29/M30 control recovers all three directions
and reduces CPU from **106.95 to 42.82 seconds**. These are measured search
and independent certification results, not height correlations. The first
144-call comparison still has a failed promotion gate. Rank32 remains **UNKNOWN**.

The subsequent four-fibre cold pilot produces certified new directions on
all four previously unsearched fibres. Adaptive full-cloud lower bounds are
**17,17,18,17**, from explicitly certified rank16 starting subgroups. Its
complete CPU saving is only **4.25%**; V3 returns one additional cloud
direction overall. This is a positive first-direction policy result with
an unresolved deep-amplification boundary.

The [protocol](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/protocol.json),
[inputs](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/inputs.json),
[fit](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/fit.json)
and [training replay](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/training-replay.json)
are retained separately from both earlier experiments. The failed V3
promotion gate remains failed. No old result is retagged.

## What cancellation can make cheap

Use the jointly primitive pointed forms `N,D`, with `D∘M = rho*q`, from the
[finite-cancellation law](FINITE_CANCELLATION_PREDICTOR_2026-09-14.md).
For an old primitive address `w=(m,n)`, write `H0=||w||∞` and `theta=w/H0`.
Suppose `T=[[p,r],[0,1]] U`, with `U` integral unimodular, is an admitted
prime-neighbour transformation. Define

\[
 \chi_p(w)=1_{m\equiv rn\pmod p},\qquad
 z=\frac{\operatorname{adj}(T)w}{p^{\chi_p(w)}}.
\]

The numerator has content exactly `p` in the distinguished ball and `1`
outside. Invertibility away from `p` and the adjugate identity prove this;
primitivity excludes any further factor at `p`. Thus the new address is
primitive and its height is exactly

\[
 \boxed{\frac{H_{\rm new}}{H_0}
 =\frac{\|\operatorname{adj}(T)\theta\|_\infty}{p^{\chi_p(w)}}.}
\]

If `p^e` is the common coefficient content of `N∘T,D∘T`, the previously
proved neighbour law also gives

\[
 \frac{g_{\rm new}}{g_{\rm old}}=p^{4(1-\chi_p(w))-e}.
\]

This identifies the local event that can compress a coordinate: membership
in the multiple-root ball, together with a favourable real orientation.
The radius affordable in a new box of height `H` is
`H*p^chi / ||adj(T)*theta||∞`. Large *absolute* `g` is not the reward:
`H0^4=Hx*g/S` still has `g` in the numerator. None of these identities says
the ball contains an independent rational point.

The Hessian/derivative corollary in the earlier validation supplies the local
interpretation: outside the stated exceptional primes, positive primitive
cancellation requires a multiple root of `q` after its content correction.
The new implementation uses those roots to construct integral neighbours;
it does not factor discriminants, optimize `C/L`, or compute a gcd tree on
every candidate model at every small prime.

The policy uses information available from `q`; it does **not** establish
additional predictive information in a gcd feature beyond `q` and its
derivatives. Its improvement concerns how that local information schedules
actual search boxes.

## Target-blind features and fitted inputs

The baseline chart and at most two integral neighbours are constructed with
the existing factor-free mapper. Their actual boxes are deduplicated by the
retained exact signed-permutation criterion. The new target-blind input is
the old quartic's distribution across the distinguished residue balls.

For each prime actually used by a neighbour, a q-only partition has depth
four and a 512-node refinement cap. Already-open siblings receive terminal
leaves. On resolved square leaves the reference weight is exactly
`projective_mass * p^(v(q)/2)`. Unresolved leaves remain explicit; normalizing
the resolved mass is a scheduling assumption and supplies no completeness
or global solubility assertion.

A single clipped least-squares mixture of uniform projective mass and that
differential reference is trained on **1,866 j groups** of the finite-cancellation
corpus. Each j group has weight one, regardless of its number of signs,
anchors or target representatives. All 53 previous CPU curves, Curve302 and
the 41 new CPU controls are excluded. The fitted differential weight is
`0.9441799733822807`, with at least 2% uniform mass required by the policy.
Five j-fold fits have weights between 0.9398 and 0.9487. Withheld membership
Brier error is 0.16257, versus 0.42813 for uniform mass. This is a membership
calibration, not a claim of useful search performance.

The radius prior comes from **545 completed development V2 factor-free
calls**, with **39 certified positive calls**. For each positive call, use
the smaller exact coordinate height of the two signs of its first certified
witness. It is not necessarily the smallest independent representative in
the entire returned cloud. Unsuccessful completed calls remain censored at
125000; their unknown target heights are not manufactured. Anchor hazards
use the fixed buckets `0..3`, `4..15`, and `16..47`, with fixed smoothing.
This conditional radius prior does not model an unknown tail beyond the
development box or prove a probability of new directions on unseen fibres.

The independent training replay checks **157,180 residue leaves** on
**3,816 anchors**, reconstructs exact rational differential weights, refits
the mixture, checks withheld calibration, and reconstructs all 39 witness
heights against retained independent V2 endpoints. Fitting takes 8.27 CPU
seconds; this replay takes 9.81. Both are charged in the campaign-cost
sensitivity, separately from per-arm execution.

## Adaptive allocation

[`cancellation_scheduler.py`](../cas/cancellation_scheduler.py) builds a
finite real quadrature on both square-boundary charts. The weights approximate
`|m dn-n dm| / sqrt(q(m,n))`; a quadrature miss falls back to the full boundary
instead of declaring real insolubility. This numerical grid is not a bound.

Neighbour roots at the **same prime are mutually exclusive categories**.
Only distinct primes are combined by a product-reference assumption. Each
real/local state determines exact height ratios for every model. The empirical
radius CDF then assigns a surrogate detection mass to a box. The candidate
levels are `8000`, `32000` and `125000`.

The scheduler chooses the highest estimated *additional* detection mass per
CPU among existing boxes and opening the next retained anchor. It uses the
union already searched on that anchor: a completed miss raises the covered
radius statewise, while a timeout removes **no** exposure. Completed call
costs update a model with fixed overhead and squared-height scaling, shrunk
toward four development observations. Candidate preparation, scoring,
discarded models, repeated inner boxes and unsuccessful calls are all actually
charged. The cost model itself is not a timing certificate.

The worker stops on its first independently certified direction. It does not
continue with an enlarged basis and a stale landscape. Basis changes require
a new compatible bank and a separately frozen continuation.

## New withheld CPU experiment

The roster reuses the V2 eligibility audit without repeating its census.
All previous CPU groups and Curve302 are excluded. It selects every remaining
eligible deep group within the fixed six-per-family cap, then six shallow
groups per family using a new deterministic j-hash order. The result is
**36 shallow and five deep curves**. Missing deeper banks are not rebuilt.

All 41 starting subgroups and larger known endpoints pass portable and
independent Sage finite-group checks. Endpoint points and ranks live in a
separate oracle opened only by preflight and post-execution analysis. The
worker sees the equation, known subgroup proof, retained centre words and
identifying metadata. This is code-level target separation and new CPU
holdout, not an externally pristine corpus.

The three arms are:

| Arm | Allocation |
|---|---|
| `factor_free` | Original retained V3 order, one factor-free box at 125000 |
| `adaptive_uniform` | Adaptive union of the same neighbour boxes, uniform residue mass |
| `adaptive_local` | Identical adaptive policy with the fitted local residue mass |

The uniform arm still uses multiple-root information to construct neighbours.
It ablates the fitted residue **probabilities**, not all use of local
arithmetic. A comparison with V3 can also reflect multiple box sizes and
anchor reordering; the two adaptive arms isolate the additional weighting.

Each arm has 40 CPU seconds for search, at most 48 retained anchors and 144
point calls, and a five-second wall cap per call. An outer 120-second wall
limit and a 100-second process CPU safeguard bound certification work too.
All six arm permutations cycle through the roster, one arm at a time, with
BLAS/OpenMP thread counts fixed to one. The
[supervisor](../cas/run_cancellation_scheduler.py) uses an exclusive lock,
retains each call and atomic event checkpoints, resumes only unstarted arms,
and refuses to silently rerun an interrupted or unreceipted arm.

The outer CPU meter includes interpreter startup, known-subgroup verification,
all preparation, policy scoring, backend calls, admission, exact transcript
replay and **independent Sage finite-group rank certification**. Successful
results are published only after this independent check. Unsuccessful arms
also pay for replay and their starting-subgroup check. Common retained
landscapes define the input interface: this is not a measurement of cold
whole-fibre construction. Earlier corpus creation and development are not
invented zero-cost inputs.

Promotion requires the local arm to have at least as many successes, at least
10% more successes per complete CPU, and a paired family/stratum bootstrap
interval wholly above parity against **both** comparators. There are 10,000
fixed-seed draws; each comparison uses a central 97.5% interval for the two
predeclared promotion comparisons. Zero-success denominators leave the gate
UNKNOWN. No outcome-based extra cases, timing repeats or policy changes are
allowed. The report also charges all fitting, training replay and shared
preflight CPU to the local candidate alone as a conservative campaign-cost
sensitivity.

In the first protocol, fresh-fibre execution was gated on both comparisons passing without
infrastructure ambiguity. Its next input must be a source-bound generic-only
fibre packet at a previously unsearched address, with cold construction charged
and a new paired finite budget frozen before search. The present known-control
run cannot be relabelled as fresh discovery or rank32 progress. That gate
failed; the independent successor below supplies the separate later gate.

## Replay and implementation boundaries

The first comparison completed **123 arms, 3,416 point calls and 106 new
control-direction proofs**, with no map failure, backend failure or timeout.
Every arm paid for its independent rank verification inside the outer meter.

| First comparison | Directions | Complete CPU seconds | Point calls |
|---|---:|---:|---:|
| Factor-free V3 | 36/41 | 410.21836 | 463 |
| Adaptive uniform | 35/41 | 258.55665 | 1,597 |
| Adaptive local | 35/41 | 223.64384 | 1,356 |

Local weighting has **1.15611 times** the recoveries per CPU of uniform
weighting, with paired central-97.5% interval **[1.03828,1.29015]**. This
passes the probability-weighting ablation gate. Relative to V3 the rate ratio
is 1.78330, interval [1.47318,2.09145], but the gate fails because one V3
success is missed. Neither adaptive arm gains on a V3 miss. Literal withheld
representatives occur in 30,14,20 cases respectively; independently certified
alternative representatives count as directions, not literal recoveries.
These measurements are actual search and certification costs, not height
correlations. They still do not establish a universal speedup or a new record.

The [complete report](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/summary.json)
retains stratum breakdowns and all paired rows. Charging all 65.34 seconds of
fitting, training replay and shared preflight to the local candidate alone
reduces its rate ratio to V3 to 1.38009, and to uniform to 0.89471. Thus the
local increment pays at the measured execution interface, but not under this
deliberately conservative one-batch allocation of every shared offline cost.
Amortization cannot be omitted from a later production claim.

The [byte-checked replay bundle](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/replay-manifest.json)
retains 6,336 files, including raw primary and separately labelled development
receipts. The initial source lock omitted some transitive repository helpers.
A [supplementary audit](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/dependency-audit.json)
binds all 44 additional helpers to their unchanged start-commit bytes. It is
explicitly a post-launch audit, not a retagged original protocol. The successor
locks those helpers before execution.

## Independent successor at the same CPU allowance

The sole successor change is to allow every declared anchor/model/height job:
`48*3*3=432`, instead of stopping at 144 calls while CPU remains. The original
144-call experiment and its misses are preserved. The best development arm,
`adaptive_local`, becomes a single preselected candidate, tested against V3
on **24 different curves**, four from each of six families. It retains the
same 40-second search CPU and 125000 height ceilings. It does not retry or
extend any of the previous 41 cases.

The [new protocol](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4/protocol.json)
excludes all earlier CPU groups and the 24 new groups from fitting. Its
identical mixture fit on 1,842 j groups has weight 0.9435116837720751.
Independent input/endpoint checks and 154,130 training-leaf checks pass.
The new decision rule tests the complete policy against V3, without claiming
that every part of an improvement is due to the probability weights. A pass
permits a separately frozen transfer test on the complete retained
Curve302 M27/M29/M30 banks. A failed gate permits no automatic expansion.

Commands for the sealed successor use
[`cancellation_scheduler_round4.py`](../cas/cancellation_scheduler_round4.py),
with `run`, `report`, `replay` and `pack` as distinct operations. No point
search is invoked by reporting, replay or packaging.

The [completed successor report](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4/summary.json)
records:

| Independent successor | Directions | Complete CPU seconds | Point calls |
|---|---:|---:|---:|
| Factor-free V3 | 22/24 | 325.341976 | 412 |
| Adaptive local | 22/24 | 204.298177 | 1,317 |

The rate ratio is **1.59249**, with the predeclared central-97.5% paired
family bootstrap interval **[1.15894,2.36396]**. There are no discordant
recoveries. The candidate saves **37.20% CPU** and passes every successor
gate. Charging all 43.60 seconds of fitting, training replay and shared
preflight to the candidate still gives a rate ratio of **1.31242**.
This is evidence for this finite policy on independent CPU holdouts, not
a guarantee on all subgroups or families.

All 48 arms complete and all 44 enlarged subgroups independently certify.
There are 1,728 complete point calls and one timeout on a last call with
0.0293 seconds left; both arms fail on that curve. The timeout contributes
no completed prefix and its CPU is charged. There are no preparation
unknowns. The [execution audit](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4/execution-audit.json)
retains that boundary. The byte-checked successor archive contains 3,200 files.

## High-rank transfer

The passing successor enables a separately
[frozen transfer protocol](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_transfer_v1/protocol.json).
It loads the **complete 256-centre banks** for Curve302 M27, M29 and M30,
checks their original input hashes and independently verifies all three
starting subgroups. No winning centre, target coordinate or truncated bank
chooses the inputs. The curve was excluded from fitting. Each arm has 150
search CPU seconds and up to 2,304 declared jobs; heights are unchanged.

| Known control | Certified output | V3 complete CPU | Adaptive complete CPU |
|---|---:|---:|---:|
| M27 | 28 | 38.208403 | 5.449909 |
| M29 | 30 | 22.788957 | 23.568068 |
| M30 | 31 | 45.949575 | 13.797233 |
| Total | All three recovered | 106.946935 | 42.815210 |

The [transfer report](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_transfer_v1/summary.json)
passes its declared aggregate saving and per-case slowdown gates. CPU falls
**59.97%** overall; M29 is **3.42% slower**, which remains visible. All six
rank proofs are independent and included in these costs. The archive contains
1,571 files. Three subgroups on one known curve do not provide three
independent-population samples or a new rank31 discovery.

The transferred policy is a measured candidate for subsequent bounded search.
This result does not change the baseline of every existing runner: cold
generic-only transfer must be measured separately, and the original V3,
uniform-weight ablation, earlier failed selectors and all misses remain
available as regressions.

## Fresh fibres and cold costs

The [fresh protocol](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_fresh_v1/protocol.json)
is frozen only after the transfer gate passes. It selects four nonbaseline
admitted X948 parents by a fixed hash of the family name, each at zero-based
address 64 of the old parameter plan, `t=-11/18`. These addresses have no
retained exposure or checkpoint. The selection reads no specialized rank,
point or j value. It uses existing certified generic parents and does not
restart their stopped foundry.

Each isolated arm specializes the 16 generic sections, certifies their
independence, and constructs banks 0,1,2 using the existing exact generic
CVP recipe. Both exact CVP implementations agree on all 192 proposed classes;
48 selected anchors enter search. This complete cold preparation is repeated
and charged for each arm, with a 90-second wall limit. Both arms reconstruct
byte-identical seed/bank inputs. The warm search retains V4's 40-CPU-second,
432-call and H125000 limits. Outer safeguards are 300 wall seconds and 240
process CPU seconds, including independent proof work.

Search stops at the first independently certified direction. All points in
the returned clouds are then admitted against the original subgroup, with
no further point calls. Their larger subgroup receives a standalone portable
certificate and independent Sage finite-group verification. Both phases
and all rejected/unused work are inside the outer CPU receipt. Generic
parent construction and the earlier corpus are retained development inputs;
their historical costs are not presented as zero or charged anew per fibre.

| Fresh X948 family suffix | Initial subgroup | V3 full-cloud lower bound | Adaptive full-cloud lower bound | V3 complete CPU | Adaptive complete CPU |
|---|---:|---:|---:|---:|---:|
| 04380 | 16 | 17 | 17 | 15.720259 | 11.149112 |
| 02021 | 16 | 18 | 17 | 11.603136 | 14.751798 |
| 01191 | 16 | 18 | 18 | 15.057154 | 11.291808 |
| 01615 | 16 | 17 | 17 | 13.705256 | 16.508786 |
| Total | Four fibres | Six added directions | Five added directions | 56.085805 | 53.701504 |

All eight arms recover a first direction and all **54 point calls complete**,
with no map or cold-preparation unknown. Cold setup alone costs 28.36938 CPU
seconds for V3 and 28.22274 for the candidate. The
[fresh report](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_fresh_v1/summary.json)
retains each equation identifier, rank certificate location and paired cost.
Only after execution are j values compared with the retained corpus:
all four are distinct and absent from its 1,961 groups. Freshness is relative
to that corpus and the recorded parent/address history, not a global novelty
or rank-record claim. The byte-checked bundle retains 462 files, including
generic sections, cold derivations, input seals and all eight cloud proofs.

The candidate saves **4.25%** complete CPU for four first gains, but yields
**five** full-cloud directions to V3's **six**. Those directions are measured
relative to each fibre's own certified rank16 subgroup, not added together
on one elliptic curve. This tiny pilot supplies no population speed claim.
It also explains why cheaper first recovery alone cannot establish a better
route to rank32: early small-box success can sacrifice a richer cloud.

## Current policy and next gate

The finite cancellation data support a useful scheduling policy at the
retained-bank interface; they do not support a universal negative result
about predictability. The policy has produced actual independently certified
directions on fresh fibres. It has not produced rank32, a fresh rank31 curve,
or evidence of a better deep-amplifier rate.

The [exported policy](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_programme_v1/policy.json)
retains the fitted parameters, finite profile, input contract and selected
protocol hash. Use the frozen V4 fit for this candidate and retain factor-free as the paired
reference. Further work should price independently certified **cloud gains
and later gains**, rebuild compatible banks after every basis enlargement,
and validate on another disjoint cohort before scaling. The four fresh
receipts are now development data, not reusable holdouts. A continuation
must declare which current subgroup it extends, its complete CPU allowance
and a certificate endpoint. An unresolved local state or finite-column test
cannot exclude any fibre.

This completed programme does not schedule another campaign. Its current
limits, retained failures, source locks and next unscheduled gate are recorded
in the mathematical ledger and work ledger. No old experiment is overwritten.

The separate [two-direction successor](CANCELLATION_CLOUD_SCHEDULER_2026-09-14.md)
has now failed its frozen gate: 39 certified directions in460.734 CPU seconds
versus V3's40 in429.056 on24 new rank20 controls. Its full-cloud fit, partial
gains and complete costs are retained. A new-basis integration check on one
fresh rank18 output passes in5.176 CPU seconds with no point searches; later
basis-aware efficacy still needs a new disjoint control cohort.

## Replay commands

The fixed backend is unchanged. PARI's
[official rational-point search documentation](https://pari.math.u-bordeaux.fr/dochtml/ref/Hyperelliptic_curves.html#hyperellratpoints)
defines its finite naive-height boxes. Returned squares, coordinate maps,
selected actions and independence witnesses are replayed; completed negative
box coverage retains the stated trust in the pinned PARI executable.

From the repository root:

```sh
sage -python -m pytest -q research/elliptic-curves/tests/test_cancellation_scheduler.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 research/elliptic-curves/cas/run_cancellation_scheduler.py
sage -python research/elliptic-curves/cas/verify_cancellation_scheduler.py
python3 research/elliptic-curves/cas/report_cancellation_scheduler.py
```

The supervisor's command resumes unstarted arms of this exact protocol; it
does not authorize a new run or increase its allowance. Replaying does not
repeat point searches. The four focused regressions check the exact prime
height/gcd law, same-prime categorical mass, timeout/overlap handling and
failure on a false prime-neighbour input. Isolated negative and positive
development smoke runs are retained outside the primary benchmark.

The final [programme audit](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_programme_v1/audit.json)
checks all four stages' complete receipts, training/CPU separation, retained
gates, source hashes, point files and replay bundles. Its default command is
an integrity/accounting audit of arithmetic already independently verified
inside the timed arms; `--arithmetic` additionally reconstructs the rank
proofs, policy traces and cold banks without making point-search calls:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  sage -python research/elliptic-curves/cas/verify_cancellation_scheduler_programme.py
# Separate, optional arithmetic replay of all completed stages:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  sage -python research/elliptic-curves/cas/verify_cancellation_scheduler_programme.py --arithmetic
```

The primary machine was WSL Linux on an AMD Ryzen9 9950X3D, Sage10.9 and
PARI2.15.4; the [environment receipt](../../artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3/environment.json)
retains versions and binary hashes. All stage workers use one BLAS/OpenMP
thread. Reproducibility means replayable inputs, scheduling rules, recorded
timing decisions and exact proofs, not bit-identical CPU times on another
machine. The cost adaptation depends on measured times; a new execution can
legitimately follow a different trace while the retained trace replays exactly.

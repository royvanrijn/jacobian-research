# Orbit8044: a parametrized M18 seed factory

**Implemented and validated alongside the running determinant1092 funnel.**
All seven reduced rational parameters of height at most2 produced independently
certified M18 subgroups in57.998 supervised seconds. One is the already known
small conic fibre; exact rational-isomorphism deduplication admitted the other
six to a durable V3 queue. All six bounded V3 runs completed at18 after114 charts each, with independent
replay; see the [frozen M18 comparison](M18_LANDSCAPE_COMPARISON_2026-09-08.md). A
larger4,096-address batch is frozen and prepared, with no producer launched.

The [validation package](../../artifacts/generated-results/elliptic-curves/orbit8044_seed_factory_v1/manifest.json)
contains the chart identity, equations, all seven ordered18-point packets,
standalone finite-group proofs, exact dependence-filter records, deduplication
decisions, queue entries and restart checks. These are subgroup lower bounds,
not exact ranks or a claim that the six curves are new to the literature.

## Parameter map and height order

The factory uses a reduced rational coordinate on the existing orbit8044 conic:

\[
\boxed{S(u)=
\frac{5193-125501628u-193042839045u^2}
     {35630-130270680u-1324582766150u^2}.}
\]

This factory coordinate `u` differs from the old slope coordinate in the
original conic certificate. It comes from PARI `qfparam`, reducing the third
form of the primitive conic at the previously constructed rational conic
point. This is an exact reparametrization of the same bisection. It does not
change the generic rank of the original K3 family.

Writing the numerator and denominator as `X(u),Y(u)`, the retained quadratic
`Z(u)` satisfies the coefficient identity

\[
Z^2=4940136120823281Y^2-21931265870123820XY
       +51437496474840100X^2.
\]

The producer verifies this identity and that `X,Y` have no common root. Their
ratio is a degree-two map to the reduced parent parameter `s`. The previously
verified parameter matrix transports `s` to the original parent coordinate.
All coefficient arrays and input hashes are retained in the
[chart certificate](../../artifacts/generated-results/elliptic-curves/orbit8044_seed_factory_v1/parameter-chart.json).

Enumerate `u=a/b` with `b>0`, `gcd(a,b)=1`, in increasing
`(max(abs(a),b),b,a)` order, including `0/1`. Evaluate the quadratic forms
homogeneously, then reduce `s=X(a,b)/Y(a,b)` exactly. The larger prepared
batch is the first4,096 addresses with height at most64; it is a finite prefix,
not a claim to have processed every rational parameter. Poles, singular fibres
and exceptional residual maps receive explicit statuses.

## Seed and queue gates

Each smooth ordinary fibre follows the requested sequence:

1. Specialize and certify the17 inherited sections on the exact reduced model.
2. Solve the frozen residual quadratic and line to construct both orbit8044
   points. Check both original equations and their exact reduced-model transports.
3. Check signed generic singles and pairs by exact rational addition. If one
   conic branch has an inherited word, the exact trace identity supplies the
   other word. No finite-mod2 failure or numerical score proves dependence.
4. Run the finite-group extra-direction certificate. Stop admission at the
   first proved extra point and export exactly the inherited17 plus that point.
5. In a separate Sage process, independently verify the rational points,
   specialization, conic incidence and18 independent finite-quotient columns.
6. Bucket by exact rational `j`, then test rational isomorphism. Same-`j` twists
   remain distinct. Every independently certified class without a duplicate
   receives an immutable V3 queue entry, without a score cutoff or seed quota.

The dependence filter is checked against the earlier exceptional split fibre
`s=-528/3635`, whose two conic points have exact inherited words, and against
a conic M18 where neither branch may be rejected. The residual quadratic
constructor agrees with four independently retained two-branch controls;
the nonsplitting quartic-seed fibre is correctly rejected by that constructor.
Unresolved extra directions retain `UNKNOWN` semantics.

No quartic point search is used in seed production. V3 receives exactly one
certified extra direction and keeps its original numerical policy. Deduplication
uses the factory's earlier admissions and a frozen snapshot of the other
certified seed equations. Later discoveries in other campaigns are not
silently added to that snapshot.

## Pilot result and cost

| `u` | Reduced parent `s` | Certified subgroup | Admission |
|---|---|---:|---|
| `-1` | `1004777772/6898189895` | 18 | V3 queue |
| `0` | `5193/35630` | 18 | Known rational-isomorphism duplicate |
| `1` | `1609736129/11039275010` | 18 | V3 queue |
| `-2` | `257306782577/1766023495870` | 18 | V3 queue |
| `2` | `257474118081/1766197190110` | 18 | V3 queue |
| `-1/2` | `2073030269/14240022390` | 18 | V3 queue |
| `1/2` | `7159030427/49068265370` | 18 | V3 queue |

The six admitted fibres have `j` numerator heights974–1,183 bits. Their seed
proofs are affordable; whether V3 finds further directions at these heights
is a separate experiment. The first exact M18 V3-map preflight passes.
No amplification outcome is asserted by this seed-validation package.

## Checkpoints and finite budgets

The [controller](../cas/run_orbit8044_seed_factory.py) has separate producer and
V3 locks, ledgers and budgets. A seed attempt and its independent verification
each have120 seconds and1GiB RSS. The producer has a four-hour aggregate
supervised allocation. V3 uses one worker with3GiB RSS, up to14,400 seconds per
fibre including replay, and a twelve-hour aggregate allocation per frozen run.
Queue admission does not promise every queued fibre finishes within that budget.
Unprocessed entries remain queued.

Completed stage receipts and immutable outputs are checked on resume. An
interrupted stage without a completed receipt consumes its reservation and
remains censored. A crash between queue publication and admission recording
is idempotent. The real pilot resumes with all14 completed proof stages and
all six queue files unchanged. The V3 consumer drains the queue snapshot
present at invocation; another invocation can consume later arrivals within
the same remaining aggregate budget.

The [pure arithmetic/queue contracts](../cas/orbit8044_seed_factory.py) and
[Sage worker](../cas/orbit8044_seed_factory_worker.py) are new files. The live
original90 campaign, completed conic follow-up and all frozen V3 sources are
preserved.

From the repository root, inspect the running pilot consumer:

```sh
python3 research/elliptic-curves/cas/run_orbit8044_seed_factory.py status \
  --directory research/artifacts/local/elliptic-curves/orbit8044-seed-factory-pilot-v1
```

The larger batch is already frozen. These separate commands start its producer
and consume its resulting queue; they are not implicitly run by `freeze`:

```sh
python3 research/elliptic-curves/cas/run_orbit8044_seed_factory.py launch \
  --directory research/artifacts/local/elliptic-curves/orbit8044-seed-factory-production-v1
python3 research/elliptic-curves/cas/run_orbit8044_seed_factory.py launch-v3 \
  --directory research/artifacts/local/elliptic-curves/orbit8044-seed-factory-production-v1
```

The same commands resume within existing budgets. A fresh protocol uses
`freeze --parent-run ... --max-height H --maximum-parameters N`, with explicit
`--exclude-run` snapshots as needed. Source/input/runtime changes fail closed.

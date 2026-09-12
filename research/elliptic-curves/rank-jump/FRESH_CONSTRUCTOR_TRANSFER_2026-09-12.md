# Fresh dependency-to-point transfer: frozen inputs, commissioning pending

The next milestone is one fresh MW16-05 fibre on which equation and generic
sections produce a fresh dependency, an independently new strict class, a
blind rational lift and a certified rank gain. Freeze the algorithm across
parameters; rebuild the arithmetic objects. The fixed-word finiteness theorem,
column 7 and the completed carrier bank are not new search inputs.

**Only input selection and a historical cost audit have run.** No new number
field, principal relation, strict class, cover-point search or V3 run is claimed.
The [nine-case stage ledger](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/initial-stage-ledger.json)
records every stage as `NOT_STARTED`, including the separate commissioning
control. These are not preparation timeouts or negative class results.

## Frozen fresh inputs

For reduced `t=m/n`, `n>0`, the declared band is

```
8 <= max(abs(m),n) <= 32.
```

The [selection ledger](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/selection.json)
orders the entire finite band by SHA256 of the fixed seed, a newline and the
canonical rational string, with `(m,n)` breaking a hash tie. The seed is
`MW16-05 fresh dependency transfer v1; 2026-09-12`.
No rank, score, conductor, norm-smoothness estimate or distance from `3/17`
enters the order. Both signs are included.

The eight selected parameters, in frozen order, are

\[
 27/4,\quad -29/4,\quad18/17,\quad-4/19,\quad
 32/3,\quad-17/14,\quad-7/17,\quad4/11.
\]

The first eight ordered candidates all pass the declared exclusions. Exact
Q-isomorphism checks use both `c4` and `c6`, retaining constant twists; equality
of `j` alone is insufficient. The
[alias snapshot](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/equation-aliases.json)
contains the 491 equations in the frozen curve database and its source hash,
with historical points, ranks, scores and parameter labels omitted. The control
equation and earlier selected equations are also excluded. This proves freshness
relative to that declared snapshot, not absence from every historical exposure.
There will be no replacements after arithmetic begins.

Each input contains the exact specialized short equation and all sixteen
generic points, checked by rational substitution. An
[independent Sage replay](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/independent-input-replay.json)
reconstructs the order, uses Sage's Q-isomorphism test for exclusions, checks
the 491-equation projection and all 128 fresh point equations, and replays
3,884 timing receipts. Generic independence on these fibres remains pending;
an equation identity is not an independence certificate.

## The complete reference cost is not yet known

The [cost audit](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/historical-cost-audit.json)
totals **3,202.922 recorded wall seconds**, about 53.38 minutes, in disjoint
retained stages. It includes the early relation searches and their available
replays, nine initial target waves, failed BNF attempts, the fixed target wave,
the adaptive completion and the small-representative stage. Enclosing worker
times and their component times are not added twice.

This is an incomplete subtotal. Cold factorization and prime proofs, some
order/factor-base/local preparation, target setup and audits, standalone
strict extraction and verification, compaction and the full lifting route
lack a complete common meter. Implementation time is another separate item.
The audit does not convert an absent timing into zero or equate wall time with
CPU time. The later 18.401-second cover recovery is not a constructor budget.

The proposed resource rule is in the
[design record](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/design.json):
first meter a successful fresh commissioning route at `3/17`, including
preparation, dependencies, failed attempts, compaction, reduction, lifting
and independent verification. If its charged cost is `C`, freeze
`min(7200,ceil(2*C))` seconds per route per fresh fibre. V3 receives the same
allowance. Proposed commissioning cap: two hours; one worker and 2 GiB RSS.
The worst-case aggregate ceiling, including both routes on all nine fibres,
is 36 worker-hours. This is a ceiling, not a feasibility or runtime forecast.

**The resource question is pending.** No numeric fresh-fibre allowance has
been substituted for the missing complete measurement. A failed or censored
commissioning run cannot supply `C`. The executable constructor policy and
all source hashes must also be frozen before commissioning, then remain
unchanged across the fresh panel.

## Reuse and implementation gates

The successful reference is
[the adaptive small-representative constructor](TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md),
with its complete principal dependencies. Reuse its norm-form reduction,
targeted smoothness testing, exact ideal relations and strict extraction.
Protect only independently certified generic anchors; include small
representatives whose ideal columns have already become pivots. Keep every
outside coordinate until exact cancellation. Evaluate factored classes away
from their full support.

The original entry points are not ready-made multi-fibre workers. They read
fixed reference paths, assume generic strict dimension six in the extractor,
and use residue-degree-one dyadic generators in some verification routines.
Those assumptions require a parameterized adapter with complete local data
and independent tests. Copying the reference's arrays or changing only its
parameter string would not execute the proposed experiment. This adapter
remains to be implemented and commissioned; no hidden ready-runner claim is made.

Every fibre rebuilds its field, maximal order, local maps, generic classes,
ideal labels and principal dependencies. No expected strict dimension, final
rank or class-group quotient size is a stopping target. A full class-group
upper bound is not a prerequisite.

Retain at most the first two independently new strict classes. Order them
by completed extraction checkpoint and deterministic coefficient-kernel order,
admitting each only after independent strictness and incremental independence
against the generic span and earlier retained candidates. Freeze this order
before arithmetic; do not select classes by observed lifting difficulty.

The lift worker receives only the frozen covers, rational norm data and model
maps. Reuse the verified reduction route from
[blind recovery](BLIND_CONSTRUCTED_CLASS_RECOVERY_2026-09-12.md), retaining every
coordinate transformation. No exceptional points or sibling V3 output may
enter it. A missed search leaves rational solubility `UNKNOWN`, not Sha.

## Separate endpoints

| Stage | Required evidence |
|---|---|
| Arithmetic preparation | Exact field/model transport, certified maximal order and the complete local data used by the constructor |
| Class construction | Explicit factored class, all strict conditions and independent nontriviality modulo the generic span |
| Rational lifting | Exact cover point, cubic-field identity, full transport and original elliptic equation |
| Rank gain | Fresh certificate for the generic sixteen plus the recovered point |
| Ordinary ideal-class gain | Separate half-ideal/Artin detection; the possible unit kernel is not ignored |

Preparation failure, no class within budget, unresolved candidate verification
and unresolved rational solubility must have different ledger states. Stop and
retain completed mathematics at every bound. No enlargement or replacement
is implied by a miss.

V3 starts independently from the same generic subgroup. Seal both outputs
before a comparison worker reads either together. **Transfer success** is a
fresh independently certified class-to-point chain even if V3 finds the same
direction faster. **Complementarity success** requires an additional certified
direction in the combined subgroup beyond V3's sealed matched endpoint. If
the combined independence test is inconclusive, complementarity stays `UNKNOWN`;
matching separate rank totals proves neither dependence nor independence.

All preparation, dependency work, reduction, failures and verification count
towards the class route. Keep remote and local accounting distinguishable.
Eight fibres do not establish a success rate, an infinite-family theorem, or
uniform strict/ordinary ideal-class transfer.

Replay the completed preflight only:

```sh
sage -python research/elliptic-curves/rank-jump/verify_fresh_constructor_transfer.py
```

This does not launch commissioning or any fresh arithmetic experiment.

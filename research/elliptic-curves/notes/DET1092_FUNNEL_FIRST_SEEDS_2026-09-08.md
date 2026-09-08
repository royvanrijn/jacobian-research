# First certified seeds and a rank21 cascade from the determinant1092 funnel

**Verified examples.** Three fibres from the audited ten-million-parameter
population have independently proved M18 subgroups. The fibre at
`s=1926/2699` has subsequently reached a separately verified rank21 subgroup
in unchanged V3. This certifies a completed epoch; its cascade remains active.

| Reduced parameter `s` | Selection | Certified lower bound in this package | `j` numerator/denominator bits |
|---|---|---:|---:|
| `1117/2193` | Original90, bounded first-point constructor | 18 | 486/469 |
| `1926/2699` | Exact conic-splitting follow-up | 21 | 460/449 |
| `2953/1671` | Exact conic-splitting follow-up | 18 | 509/492 |

The [certificate package](../../artifacts/generated-results/elliptic-curves/det1092_funnel_first_seeds_v1/manifest.json)
retains equations, ordered points, exact generic-section prefixes and complete
finite-group proofs. All three equations are unmatched under rational
isomorphism in the pinned630 public and201 local equations. This comparison
does not establish literature-wide novelty, a conductor improvement, or an
exact rank. The original90-fibre search continues under its frozen protocol.

## Arithmetic selection and conic confirmation

The [funnel](DET1092_SEARCH_FUNNEL_2026-09-08.md) completed all ten million
equation records and their independent replay. Its fixed trace-only retention
gave8,999 arithmetic candidates and90 seed inputs. None of those retained
candidates splits the frozen orbit8044 conic.

A complete read-only audit of the same ten million records found exactly
three rational conic-splitting fibres. A separate protocol selected all three,
without replacing any of the original90. Exact specialization of the conic
maps supplied candidate points; this follow-up allows zero quartic seed-search
charts. Two fibres passed standalone M17-plus-one-point independence proofs.
The first admitted point on each fibre is its complete exported seed.

The remaining split address, `s=-528/3635`, has two conic points with exact
words in the inherited ordered basis:

\[
P_0=-M_2+M_{15},\qquad P_1=M_{10}-M_{14}-M_{15}.
\]

These rational group identities prove that these two witnesses add no rank.
They do not bound the rank of the elliptic curve. The original finite-mod2
failure alone was insufficient to prove dependence. A bounded height-based
word proposal subsequently passed exact rational addition. Its initial run
retained `UNKNOWN` because default decimal printing lost proposal precision;
the retained second invocation sets PARI output precision to120 digits.
Only the exact group identities support the conclusion. A separate one-layer
halving attempt produced one rational half but no certified extra direction;
its unchanged earlier result remains available in local evidence.

## Independent proofs and amplification

The [seed checker](../cas/verify_det1092_funnel_seed.sage) recomputes each
equation and all18 rational points, checks the exact inherited17-point prefix
and original parameter transport, and checks conic incidence where applicable.
Complete finite elliptic groups and their doubled subgroups yield18 independent
quotient columns. An odd-order good reduction excludes rational2-torsion;
infinite descent proves independence. Positive controls pass on both a quartic
seed and the separate small conic seed. Altered off-curve and dependent points
are rejected.

The [epoch checker](../cas/verify_det1092_funnel_epoch.sage) independently
performs the same elementary finite-group proof on the21 retained rational
points of `s=1926/2699`, preserving both the generic17 and initial18 prefixes.
It imports no point constructor or repository rank/Kummer backend. The
[standalone proof](../../artifacts/generated-results/elliptic-curves/det1092_funnel_first_seeds_v1/rank21/standalone-proof.json)
binds the completed epoch, its point cloud, and the frozen input protocol.
This establishes the rank lower bound without waiting for the rest of V3.
Full search-policy replay is a separate terminal check.

The first main-funnel seed, `s=1117/2193`, completed114 V3 charts with no further
certified gain. Its full independent landscape, map and point-cloud replay
passed, including the mod3 and mod5 audits. This is a bounded no-gain result.
The split follow-up is sequential: the second independently proved M18 waits
for the current cascade. Each cascade retains unchanged V3 settings, a
14,400-second combined search/replay cap, and3GiB RSS. Each seed proof and
standalone seed verification has a120-second cap. There are exactly three
fixed follow-up fibres, with one attempt per stage.

## Reproduction

The [follow-up controller](../cas/run_det1092_split_followup.py) freezes the
three extracted equation rows, confirms all seeds, independently verifies
successful seeds, and then runs their cascades. Completed stages and raw
charts stay under `artifacts/local/elliptic-curves/`; the compact package
preserves the mathematical witnesses and their hashes.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_funnel_seed.sage \
  --seed research/artifacts/local/elliptic-curves/det1092-funnel-conic-split-seeds-v1/seeds/funnel-002537010 \
  --output research/artifacts/local/elliptic-curves/det1092-funnel-conic-split-seeds-v1/seeds/funnel-002537010/standalone-replay.json
```

This package records verified lower bounds at a fixed checkpoint. Later
progress belongs to the live controller and subsequently verified certificates;
the rank32 and conductor-above-rank22 objectives remain open.

# Broad-rank certified rank-22 ledger

This note records the separately frozen rank-22 cohort from the completed
broad-rank-v1 campaign. It is distinct from the earlier rank-greater-than-22
publication cutoff: its 46 source packets were replayed and checked against the
pinned 445-curve inventory before being rendered into the curve list.

The [canonical snapshot](../../artifacts/generated-results/elliptic-curves/broad_rank22_ledger_snapshot_v1.json)
and [fresh two-path replay](../../artifacts/generated-results/elliptic-curves/broad_rank22_ledger_replay_v1.json)
verify all point memberships, independence, rational 2-torsion exclusions,
specialized parent equations, generic point prefixes, and Q-nonisomorphism
against the baseline and within the cohort. The inventory therefore grows to
491 curves. These are certified subgroup lower bounds, not exact ranks,
conductor claims, or a worldwide novelty claim.

The rank-22 cohort has no conductor calculation in this ledger; its conductor
and minimal-model fields remain `UNKNOWN` until separately certified.

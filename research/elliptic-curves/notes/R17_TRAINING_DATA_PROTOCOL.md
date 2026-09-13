# R17 learned-score protocol — historical boundary

<!-- status-consumer: EC-K3-R17-TRAINING-EXACT-ARITHMETIC-GROUP-GATE 427bf822e774c81e -->

The complete 2 September 2026 protocol, frozen cohorts, commands and outcome
tables are preserved byte-for-byte in [the archive](../../archive/elliptic-curves/notes/R17_TRAINING_DATA_PROTOCOL.md.txt)
(`sha256: 07c89fcd33d712bdc3a97861db3ce3eed87423ae4711d73726b6c1e28c42a622`).
It does not authorize new training, feature selection, point search, or
retuning. Mathematical status remains in [MATH_STATUS.json](../../MATH_STATUS.json).

## Canonical arithmetic-grouping gate

The outcome-free audit groups the immutable 100,000-row development universe,
4,922 selected rows, 5,000-row prospective holdout and four quarantined
controls by their exact reduced rational `j`-invariant. This conservatively
groups every rational isomorphism and twist class, including `j=0` and
`j=1728`. It finds 100,000 singleton development groups, no selected/holdout
parameter or twist-class overlap, and no control overlap with development or
holdout.

The [compact input](../../artifacts/generated-results/elliptic-curves/r17_training_arithmetic_group_inputs_v1.json.gz),
[exact audit](../../artifacts/generated-results/elliptic-curves/r17_training_arithmetic_group_audit_v1.json),
and [`audit_r17_training_arithmetic_groups.py`](../scripts/audit_r17_training_arithmetic_groups.py)
authorize only the pinned v1 learned-score artifact. Any changed population,
split, holdout, controls, family equation or fitted score requires a new
outcome-free grouping audit. This closes leakage only; it is not a rank,
Selmer, or prediction theorem.

## Retained decision

The frozen bisection-gain score generalized to its 5,000-row bisection-label
holdout but missed the published rank-28 control at a one-percent review
budget. It learned bisection visibility, not extreme total rank. The
quarantine is opened, so retuning on those four controls cannot be described
as a fresh replay. More labels from this atlas are retired for the total-rank
objective; use the [programme map](../README.md) and
[algorithm memory](../../knowledge/ALGORITHMS.md) for separately scoped work.

The detailed frozen inputs, censored-label rules, score comparisons and
literature context remain in the archive. The current script interfaces are
listed in [`scripts/README.md`](../scripts/README.md) and the broader
historical evidence is summarized in
[SEED_AND_AMPLIFICATION_HISTORY_2026-09-13.md](SEED_AND_AMPLIFICATION_HISTORY_2026-09-13.md).

# Fermigier quotient fingerprints and frozen-score replay

## Historical calibration

This completed retrospective replay ranked all 60,815,684 primitive
`T=a/b` parameters in its declared box with four frozen local score orderings.
The known E22 and rank-20 controls appeared after position 2.7 million in every
ordering, so each had zero recall through the historical budget 100,000.

This says the recorded local scores do not retrieve generic-lattice escape in
the Fermigier family. It does not classify the other censored fibres, establish
an accuracy metric, or authorize a new candidate search. The controls were
known during score design and are not prospective holdouts.

The compact [quotient fingerprints](../../artifacts/generated-results/elliptic-curves/fermigier_rank_jump_fingerprints_v1.json)
and [score replay](../../artifacts/generated-results/elliptic-curves/fermigier_rank_jump_replay_v1.json)
remain available. The full original protocol, exact control positions and
slow replay commands are preserved in the
[archive](../../archive/elliptic-curves/notes/FERMIGIER_RANK_JUMP_REPLAY.md.txt).

New score work requires a separately scoped prospective protocol with exact
endpoints; this historical complete-box replay is not a work queue.

# Five new catalogue-unmatched curves with rank at least26

These five curves have globally minimal integral equations and26 exactly
independent rational points each. They are mutually nonisomorphic over Q,
absent from the pinned593-equation ICARM snapshot and absent from549 previously
measured equations. This is a bounded novelty comparison, not universal novelty.
None is being claimed to have exact rank26 or to set a rank record.

The [model certificate](../../artifacts/generated-results/elliptic-curves/compact192_rank26_models_v1.json)
contains every equation, point, exact coordinate transport, finite independence
certificate and local minimality witness. The executable
[Sage export](../../artifacts/generated-results/elliptic-curves/new_compact192_rank26_curves.sage)
provides all five curves and their130 independent-point witnesses.

| Candidate | R17 family parameter | Certified lower bound | gcd(c4,c6) |
| --- | --- | ---: | ---: |
| 07ca9-001 | 760/211 | 26 | 3136 |
| 07ca9-006 | -575/1207 | 26 | 64 |
| 07ca9-011 | 3307/1128 | 26 | 1 |
| 103b2-029 | 726/761 | 26 | 121 |
| 11952-010 | -1826/2583 | 26 | 3 |

## Exact verification

The complete initial192-curve point batch and its post-terminal point proof
and catalogue replay pass. These five are all its catalogue-unmatched results
with lower bound at least26. For each, an integral normalized Weierstrass model
is constructed with the same c4,c6 invariants as the discovery model. All point
transports and their inverses are checked exactly; the original finite
independence certificate is recomputed on the transported-back points.

The small invariant gcds permit complete cheap prime factorization. At every
prime dividing that gcd, the exact local invariant valuations or complete
normalized-coefficient exclusions rule out an integral model with invariants
divided by the fourth and sixth powers of the prime. Other primes are excluded
by their invariant valuations. This proves global minimality without factoring
the full discriminant. No conductor is claimed.

`export_compact192_rank26_models.py --check` replays these model, point and
comparison proofs and verifies the exact Sage export. Both the checker and
execution of that export pass locally. The [standalone proof bundle](../../artifacts/generated-results/elliptic-curves/compact192_five_rank26_evidence_v1.json)
requires no base archive. Its [isolated replay report](../../artifacts/generated-results/elliptic-curves/compact192_five_rank26_portable_replay_v1.json),
when present, records both the exact model-proof check and Sage-export execution.

## Subsequent work

The [complete compact192 experiment](COMPACT192_UNSEARCHED_TRIAL_2026-09-06.md)
governs the cohort's remaining evidence and inventory promotion. The separately
bounded [own-subgroup follow-up](COMPACT192_SPECIALIZED_FOLLOWUP_2026-09-06.md)
uses the additional points to seek a lower bound of28. Its outcome must not be
inferred from these26-point certificates. Whole-curve upper bounds, full
saturation, exact ranks and rank28/32 discovery remain separate questions.

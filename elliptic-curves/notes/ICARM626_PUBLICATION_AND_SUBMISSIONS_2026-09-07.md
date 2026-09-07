# ICARM #626 and further submission priorities

Our curve `new-20260905-36`, family `a1-fibration-05` at `3/17`, is now
[ICARM #626](https://elliptic-rank.icarm.cloud/curve/626), submitted by
Roy van Rijn at the catalogue timestamp `2026-09-07 07:58:17`.
Its public equation, 22-point set and exact conductor match the
[existing proof](NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md) exactly.
The unconditional statement remains rank at least 22. The separate
[conditional upper-bound proof](SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md)
gives exact rank 22 under its stated GRH assumption.

## Current publication status

The database downloaded at `2026-09-07T08:08:39.902705+00:00` contains 626
curves; raw SHA256 is
`6a7ebb045f371832f54ba7509c7c802e484de14ef0fcbf0b2b10d0d019bc96d6`.
The [retained projection and public entry](../../artifacts/generated-results/elliptic-curves/icarm626_publication_snapshot_v1.json)
include ICARM's attribution and NSF DMS 2425401 acknowledgement.
Raw downloads and their metadata are retained under
`artifacts/local/elliptic-curves/icarm626-publication-20260907/`.

The [replayable publication supplement](../../artifacts/generated-results/elliptic-curves/icarm626_publication_audit_v1.json)
compares every equation in the 201-curve V22 inventory over Q:

| Local discovery ID | ICARM ID | Certified lower bound in inventory |
|---|---:|---:|
| `new-20260905-12` | 600 | 23 |
| `new-20260906-188` | 619 | 28, from separately reproduced public points |
| `new-20260905-36` | 626 | 22 |

The other **198 equations are unmatched in this snapshot**: seven at least
27, eighteen at least 26, thirty-eight at least 25, fifty at least 24,
forty-six at least 23 and thirty-nine at least 22, each counted once.
This is not a universal novelty or discovery-priority claim. The supplement
updates publication metadata; V22 and all earlier discovery inputs and
certificates retain their original snapshot semantics.

For #626, #376, #575 and the newer #625 have smaller recorded
rank-at-least-22 conductors. Nine relevant entries lack conductors:
#537, #543, #545, #581, #615, #616, #617, #619 and #620. Thus its placement
is now fourth among recorded values in this snapshot, not a conductor record.
The third-place comparisons in the earlier notes remain historical.

## Submission advice

ICARM's [submission form](https://elliptic-rank.icarm.cloud/) and
[API documentation](https://elliptic-rank.icarm.cloud/api), checked on
2026-09-07, do not require a record. They accept equations with independently
certified rational points. Complete bad-prime lists are optional and permit
recording the conductor. An improved point basis can update an existing
entry; supplying missing bad primes can backfill its conductor without a
rank gain. New submissions and improvements receive distinct attribution.

There is no obligation to submit all 198 curves. Our recommended first batch
is the seven unmatched curves with our strongest local bound, rank at least 27:

| Local ID | Family | Parameter |
|---|---|---|
| `new-20260906-40` | `074d9` | `2818/1535` |
| `new-20260906-71` | `103b2` | `3726/881` |
| `new-20260906-41` | `11952` | `-2448/11` |
| `new-20260906-72` | `11952` | `2012/211` |
| `new-20260906-48` | `11952` | `2828/2015` |
| `new-20260906-186` | `11952` | `4286/1881` |
| `new-20260906-90` | `a1-fibration-01` | `-1867/270` |

Equations, point sets and certificate references are in the
[V22 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v22.json).
The [seven ready-to-paste submission sheets](../../artifacts/generated-results/elliptic-curves/icarm_rank27_submissions_v1/README.md)
provide the proved minimal equations, all 27 points, commentary and JSON
payloads, also in one [ZIP](../../artifacts/generated-results/elliptic-curves/icarm_rank27_submissions_v1.zip).
Their separate 08:14 UTC catalogue check still finds no matches. All seven
existing finite-reduction proofs and exact point transports replay; the optional
bad-prime fields are left blank because the complete lists remain unknown.
This package makes no external submission.

The [187-curve conductor audit](INVENTORY187_CONDUCTOR_BOUNDS_2026-09-06.md)
establishes no record improvement; bounds above a benchmark leave actual
conductors unresolved. These are optional high-rank examples, not established
new records. Additional curves can be selected for size or construction
interest; a bulk upload is not necessary. Recheck the live catalogue before
submission, preserve construction and AI-assistance provenance, and never
submit a conductor upper bound as an exact conductor. No further submission
was made by this update.

Replay without network access or a new point search:

```sh
python3 elliptic-curves/cas/audit_icarm626_publication.py --check
```

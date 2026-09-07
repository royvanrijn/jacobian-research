# Current ICARM database and recovery controls

The [current database manifest](../data/icarm_current.json) now contains **630 curves**,
including627–630. Exact checks pass630 discriminants and9,484 point memberships;
194 of201 local equations remain unmatched. The latest conductor comparison is
recorded in the [201-curve audit](../../artifacts/generated-results/elliptic-curves/inventory201_conductor_bounds_v1.json).

The recovery controls below retain their frozen
[626-entry manifest](../data/icarm_snapshots/626_manifest.json). Its raw SHA256 is
`6a7ebb045f371832f54ba7509c7c802e484de14ef0fcbf0b2b10d0d019bc96d6`.
The fetch at **2026-09-07 08:18:18 UTC** still contains **626 curves**.
Our seven prepared rank-27 submissions are not yet visible in this snapshot.
The database credits ICARM and NSF Grant DMS 2425401.

## Database and publication update

Compared with the earlier full 620-entry database, IDs **621–626** are new.
All other entries are unchanged, including their points and metadata; no entry
was removed. The [refresh audit](../../artifacts/generated-results/elliptic-curves/icarm_catalogue_626_6a7ebb045f37_delta.json)
checks all **626 discriminants and 9,376 point memberships** exactly.
This does not prove every published point list independent or independently
verify the catalogue's reported conductors.

The frozen626 publication view joins an exact Q-isomorphism comparison to the
frozen V22 research inventory. Local IDs 12, 188 and 36 match ICARM 600, 619
and 626 respectively. **198 of 201 inventory equations remain unmatched**,
including seven with local certified lower bound 27. The publication overlay
updates the stale V22 status for #626 without changing any discovery equation,
point, rank certificate, rank bound or public-point provenance.

Use `load_catalogue()` and `load_inventory()` from
[`refresh_icarm_local_database.py`](../cas/refresh_icarm_local_database.py)
for this current view. Frozen experiments continue to use their declared
snapshots. The [publication history and submission packets](ICARM626_PUBLICATION_AND_SUBMISSIONS_2026-09-07.md)
remain indexed separately.

## Stronger size benchmarks

The new entries change the recorded conductor minima for four rank thresholds:

| Rank threshold | Previous ID | Current ID | Previous recorded minimum / current |
|---|---:|---:|---:|
| at least 22 | 376 | 625 | about 822.2 |
| at least 23 | 539 | 624 | about 839,012 |
| at least 24 | 548 | 623 | about 3.469 |
| at least 25 | 542 | 622 | about 20,659 |

The rank-26 and rank-27 thresholds are unchanged. These are comparisons of
**reported conductors**; the new conductor computations have not been replayed
here, and nine relevant entries still have no conductor. The exact current
integers and missing-entry IDs are retained in the refresh audit.

For our #626, a hypothetical 23rd independent point would give a conductor
about **20.08 times smaller** than #624's recorded rank-at-least-23 minimum.
The historical factor of about 16.8 million used the older minimum. This is
only an updated comparison: no additional point exists in our certificates,
and the separate exact-rank-22 conclusion under the stated GRH assumption
remains unchanged. No new search on #626 is proposed by this refresh.

## Locally certified calibration inputs

The [bounded public-point replay](../../artifacts/generated-results/elliptic-curves/icarm_recent_controls_v1/replay.json)
proves all displayed points independent for **eleven curves**:
#615–620 each have rank at least 28; #621–625 have respective bounds
20, 25, 24, 23 and 22. These are public-data reproductions, not discoveries.
Each proof transports the points exactly to a short equation, computes finite
quotients at good primes no larger than 997 and excludes rational 2-torsion.
The finite signature replay shares arithmetic code with the construction;
it is not claimed as an independent second implementation.

Six controls now separate worker inputs from the withheld witnesses:

- **#619: genuine historical missed direction.** The worker receives precisely
  our old 27-point basis for local ID188. The validator retains the published
  point that extends it to an exactly certified 28-point subgroup. This control
  was already used in earlier development and is not a new holdout.
- **#615, #616, #617, #618 and #620: synthetic deletions.** The worker receives
  the first 27 points of the public basis, and the validator retains the last.
  All five 28-point groups and their 27-point seeds are certified. Artificial
  deletion does not model the distribution of directions missed by our search.

The [seed file](../../artifacts/generated-results/elliptic-curves/icarm_recent_controls_v1/seeds.json)
contains only curves, 27-point bases and seed proofs. The separate
[oracle file](../../artifacts/generated-results/elliptic-curves/icarm_recent_controls_v1/oracle.json)
holds withheld witnesses and full 28-point proofs. Workers must not read that
file, the public proof files or the public database to construct their centres.
These are retrospective controls, not target-free prospective inputs.

The concrete next experiment is an equal-budget comparison of our current
point-recovery method on these six fixed seeds. Freeze its algorithm, chart
rule, height, time limits and scoring before running; record completion and
runtime separately from success. Success requires an exactly independent 28th
direction, not recovery of one particular representative. Report the genuine
miss separately from the five deletions. No search was run by this update and
no improvement in recovery rate has yet been demonstrated.

The public commentary on #615–620 describes rescoring through `2^20` and an
iterated 2-covering search over 3,000 cosets of the growing subgroup. It motivates
testing adaptive subgroup coverage, but is not a complete reproducible algorithm
or evidence that increasing either budget improves our method. Construction
provenance for #621–625 is absent from their commentary and remains UNKNOWN;
recognizing those small models is a separate possible follow-up.

## Replay

```sh
python3 elliptic-curves/cas/refresh_icarm_local_database.py --check
python3 elliptic-curves/cas/replay_icarm626_controls.py
```

Raw database and metadata inputs remain under the paths recorded in the
manifest and refresh audit. The complete current snapshot is also retained
as compressed raw JSON under `artifacts/generated-results/elliptic-curves/`.
No network, factorization, descent or rational-point search runs during replay.

The initial uncached proof builder hit its 60-second wrapper limit after
checkpointing #615–617. Its source is preserved beside the raw intake.
The completed builder uses the existing in-memory finite quotient cache,
reproduces those same three proofs, and completes all eleven proofs and six
controls with the same prime bound. The timeout implies no mathematical miss.

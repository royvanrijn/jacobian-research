# Two new rank-at-least-27 curves from extended prime selection

Two independently certified curves improve the strongest new result of this
programme from rank at least 26 to **rank at least 27**. Both have proved global
minimal integral models and no matching rational isomorphism class in the pinned
586-equation ICARM snapshot or 299 earlier measured equations. This is a precise
catalogue comparison, not a proof of universal novelty or exact rank.

| Inventory ID | Compact family | Parameter | Certified lower bound |
|---|---|---|---:|
| `new-20260906-40` | `074d9` | `2818/1535` | 27 |
| `new-20260906-41` | `11952` | `-2448/11` | 27 |
| `new-20260906-42` | `074d9` | `808/2259` | 26 |

The first two equations are

```
y² = x³ + x²
     - 478084759064998390933143695734260558414948830620*x
     + 119142731212438304178867546335237459944511190723779503318854859356288944

y² + x*y = x³
     - 1534141051320185336576223734099764498891473676600*x
     + 724098778739678960943055805127196924900556657734889009686013364096928832
```

The [first Sage file](../../artifacts/generated-results/elliptic-curves/new_paired_rank27_curve.sage)
and [second Sage file](../../artifacts/generated-results/elliptic-curves/new_paired_rank27_curve_11952.sage)
contain all 27 rational points on their respective integral models. The
[rank-26 file](../../artifacts/generated-results/elliptic-curves/new_paired_rank26_curve.sage)
contains the third equation and 26 independent points.
The [minimal-model proof collection](../../artifacts/generated-results/elliptic-curves/paired_high_rank_minimal_proofs_v2.json)
recomputes exact finite quotient signatures and point transports without Sage.
For the first curve, `gcd(c4,c6)=16` excludes nonminimality at every odd prime,
and the discriminant valuation 8 excludes it at 2. The second curve and the
new rank-26 curve have `gcd(c4,c6)=1`. These are global minimality proofs,
not assertions from a partial factorization.

## What changed in selection

The earlier 562-prime score retained 128 addresses per compact family. It was
an avoidable mistake to treat this small-prime truncation as an adequate stand-in
for a selector using thousands of primes. Earlier controls already showed
losses, but direct character-sum table construction made extensions appear
costly. A new fixed six-curve benchmark uses PARI's finite-field elliptic point
counting: 5,978 further primes per curve completed in 1.297 seconds overall.
All 48 predetermined direct character-sum checks agree.

The new experiment freezes all **768 already retained addresses**, then computes
4,591,104 additional prime traces. It completes in 126.967 seconds, with exact
program/prime rosters and score replay. Selection adds primes 4099–32749 to the
original score. Primes 32771–65521 are a disjoint validation band and never
enter the selector. The original indices 0–3 are excluded from the new point
roster; two extended-score finalists per family are compared with original
indices 4 and 5. One shared address gives **23 distinct curve attempts**.
No catalogue equation, rank label or public exceptional point enters selection
or point searching. Catalogue comparison waits for the terminal fixed batch.

All **1,037 generic-centre boxes** complete at height 100,000, using every exact
maximum generic parity class: 43 or 49 per family. All search/admission/archive
histories and all complete point-cloud certificates replay. The observed
certified lower bounds average 24.75 in the extended-score arm and 18.25 in the
original-score arm, with 12 arm memberships each. The extended arm also has
higher mean validation-prime score in all six families.

These are finite discovery and validation outcomes. They do not prove a general
rank classifier, actual-rank enrichment, or completeness outside the previously
truncated 768-address pool. Numerical heights choose centres only; exact points
and finite reduction certificates establish lower bounds.

## New curves, known controls and retained failures

The [complete batch certificate](../../artifacts/generated-results/elliptic-curves/fresh_r17_paired_results_v1.json)
labels four known catalogue rediscoveries: 376, 390, 396 and 498. Another result
matches an earlier curve already in this repository. They do not count as new
curves. Eight new rank-at-least-22 curves remain, including the three strongest
listed above. The independently replayed
[47-curve inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v4.json)
preserves all previous IDs; an
[equation CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v4.csv)
is available.

The first four-model certificate builder incorrectly required all four high-rank
results to be new. It rejected catalogue curve 390. Its source and failed log
remain intact. Version 2 verifies the actual comparison lists and explicitly
retains that known curve alongside the three unmatched curves.

A separate performance gap appeared in standalone certification: the default
finite-reduction cache publishes individual point facts to disk. One diagnostic
was observed waiting in `jbd2_log_wait_commit`. The new explicit
[memory certificate helper](../cas/memory_rank_certificate.py) uses the same
finite arithmetic and returns the same proof schema, without publishing
individual ephemeral facts. A focused test checks exact agreement with the
original backend and rejects repeated and off-curve points. The durable
certificates retain all points, primes, signatures and torsion witnesses.
No existing frozen worker or checker was changed.

## Bounded follow-ups and remaining target

The first rank-27 curve completed a separately frozen 301-centre adaptive
follow-up using its ten discovered directions beyond the generic seventeen.
All maps and archived histories replay. Its complete **1,832-point cloud** still
certifies 27 modulo 2, 3 and 5. The
[coverage certificate](../../artifacts/generated-results/elliptic-curves/paired_rank27_adaptive_coverage_v1.json)
records the finite search scope. This does not prove exact rank or saturation.

The second rank-27 curve also completed all 301 fixed adaptive boxes. All
archived histories replay; its complete 1,791-point cloud still certifies 27 modulo
2, 3 and 5. Its [coverage certificate](../../artifacts/generated-results/elliptic-curves/paired_second27_adaptive_coverage_v1.json)
records the same bounded scope. Neither follow-up found a certified 28th direction.

A cheap local conductor audit on this curve leaves a 110-digit unfactored
cofactor and gives a 120-digit conductor upper bound. It does not beat the
listed rank-27 minimum as an upper bound; the exact conductor remains unknown.
No full discriminant factorization or descent campaign was started.

New rank-at-least-28/32 curves and unconditional exact ranks remain open.
The next selection question is how far to extend the previously truncated
candidate pool, with equal exposure and held-out validation retained.

## Reproduction

```sh
python3 elliptic-curves/cas/extend_retained_r17_prime_scores.py replay
python3 elliptic-curves/cas/certify_paired_high_rank_minimal_v2.py --check
python3 elliptic-curves/cas/certify_fresh_r17_paired_results.py --check
python3 elliptic-curves/cas/export_new_high_rank_curve_index_v4.py \
  --check artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v4.json
```

`MATH_STATUS.json` remains the status authority. Earlier rank-26 and low-conductor
proofs remain valid and retained. Full search histories, failures and raw point
outputs accompany the portable evidence; no claim depends on treating a bounded
miss as a theorem of nonexistence.

The [portable evidence supplement](../../artifacts/generated-results/elliptic-curves/paired_rank27_discovery_evidence_v1.json)
names five pinned base archives. Its
[35-stage isolated verifier](../cas/verify_paired_rank27_discovery_bundle.py)
checks exact centres and point-cloud provenance for all 1,639 charts, the
768 trace rosters, independent-point and minimal-model proofs, the 47-curve
inventory and targeted tests. It does not repeat the separately passed
admission/archive histories or run another point search. The
[completed-run summary](../../artifacts/generated-results/elliptic-curves/paired_rank27_portable_replay_v1.json)
records the actual isolated outcomes.

The initial isolated run passed 34 stages. Its inventory stage hit the
180-second limit after 26 curves, so that run is retained as incomplete for
the inventory. The separately versioned
[memory inventory replayer](../cas/replay_inventory_v4_memory.py) checks the
same 47 certificates, source bindings, stable IDs, distinctness and catalogue
exclusions, and additionally checks the CSV against the proved JSON. It
completed locally in 6.655 seconds. The
[portable completion record](../../artifacts/generated-results/elliptic-curves/paired_rank27_portable_completion_v2.json)
records its isolated result; the
[versioned evidence delta](../../artifacts/generated-results/elliptic-curves/paired_rank27_inventory_replay_evidence_v2.json)
retains both the original timeout and the replacement, without rewriting the
frozen exporter or rerunning the 34 passed stages.

The [latest-eight cross-family extension](LATEST_EIGHT_CROSS_FAMILY_INCIDENCE_2026-09-06.md)
also closes the remaining incidence coverage: all 564 pairs for the 47 curves
and twelve recorded presentations are now checked. The eight additions have
only their original parameters in their original families, so this route
supplies no additional recorded generic subgroup.

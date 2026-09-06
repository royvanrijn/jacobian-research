# A third rank-at-least-27 curve and the next fixed24 experiment

The next fixed candidate batch gives a third catalogue-unmatched curve with
**27 exactly independent rational points**, together with two further curves
with 26. Their global minimal equations and point transports are proved in the
[minimal-model certificate](../../artifacts/generated-results/elliptic-curves/next24_high_rank_minimal_proofs_v4.json).
No exact rank, universal novelty, or rank-at-least-28/32 result is claimed.

| Inventory ID | Family | Parameter | Certified lower bound |
|---|---|---|---:|
| `new-20260906-48` | `11952` | `2828/2015` | 27 |
| `new-20260906-49` | `103b2` | `2773/962` | 26 |
| `new-20260906-50` | `11952` | `-1173/127` | 26 |

The new rank-27 equation is

```
y² + x*y = x³
  - 32856125399473211842388268920550481481461685620655*x
  + 50166999670818394109703619550634496881463755230595821387380798956334398777
```

Its [Sage file with all 27 points](../../artifacts/generated-results/elliptic-curves/new_next24_rank27_11952_041_curve.sage)
and the two rank-26 files
([103b2](../../artifacts/generated-results/elliptic-curves/new_next24_rank26_103b2_021_curve.sage),
[11952](../../artifacts/generated-results/elliptic-curves/new_next24_rank26_11952_030_curve.sage))
load the exact equations and rational point coordinates.
All three are unmatched in the pinned 586-equation catalogue and 322 earlier
measured equations. These comparisons occur after the fixed batch and its
independent history and cloud replays terminate.

## Fixed selection and actual yield

The [earlier extended-prime experiment](EXTENDED_PRIME_RANK27_DISCOVERIES_2026-09-06.md)
remains unchanged. This continuation freezes the next four untried addresses
per family under exactly the same selection score, excluding original indices
0–3 and the previous 23 addresses. It adds no parameters to the existing
768-address pool. Primes 32771–65521 remain outside selection; catalogue labels
and public exceptional points enter neither selection nor the workers.

All 24 attempts complete all **1,080** declared generic-centre boxes at height
100,000, with 43 or 49 exact maximum generic parity classes per family. The
same limits apply to every candidate. Every recorded admission/archive history
and complete retained point cloud replays. The
[batch certificate](../../artifacts/generated-results/elliptic-curves/next24_r17_results_v1.json)
records lower bounds, exact generic transports, and post-terminal comparisons.

The batch contains three rediscoveries of earlier repository curves. They are
retained and labelled, without new inventory IDs. Fifteen additional distinct
catalogue-unmatched curves have lower bounds at least 22. The independently
replayed [62-curve inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v5.json)
and [CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v5.csv)
preserve all 47 earlier IDs. The lower-bound buckets are 3 at 27, 4 at 26,
11 at 25, 13 at 24, 14 at 23, and 17 at 22. All 62 have distinct j-invariants.
The earlier cross-family incidence result covers its original 47 curves;
it has not yet been extended to these fifteen additions.

## Minimality and replay gaps corrected

For the new rank-27 curve, `gcd(c4,c6)=19`, and the valuations at 19 are
`(1,2,3)`. For the `103b2` rank-26 curve the gcd is 1. Those necessary
nonminimality conditions fail immediately.

For the `11952` rank-26 curve the gcd is 1600. At 2 the valuations are
`(6,9,13)`: they pass necessary conditions for a smaller model, which does
**not** prove such a model exists. Builder version 1 stops at that gate.
Version 2 attempts scaling by 2 and fails to find an integral normalized
model. Both frozen sources and the captured version-2 failure remain retained.

Versions 3 and 4 resolve the gap by exhausting normalized equations with the
scaled invariants. Version 4 records a denominator divisible by 2 in at least
one coefficient of every candidate. This is a local obstruction, not an
assumption from an incomplete factorization. Over the p-integral rationals,
Weierstrass translations can normalize a1 and a3 to 0 or 1 and a2 to -1, 0 or 1:
use residues modulo 2 or 3 when that modulus is p, and its inverse otherwise.
For each of the twelve triples, c4 and c6 uniquely determine a4 and a6.
The exhaustive obstruction therefore excludes a smaller model at 2. At 5
the invariant valuations `(2,2,4)` already exclude nonminimality.

The full-cloud auditor now has a separately frozen version 3 using the explicit
memory reduction cache. A predetermined existing rank-27 cloud reproduces
all mathematical fields of version 2, including its chosen independent subset
and finite signatures; only schema/source bindings change. Its build and
standalone check take about 1.35 and 0.28 seconds in that regression.
Earlier frozen checkers remain unchanged.

## What the known rank-29 control says about height

A separate [retrospective visibility proof](../../artifacts/generated-results/elliptic-curves/native11952_published_visibility_v2.json)
transports the 29 published points exactly to the native control model and
rechecks their independence. It locates both signs in all 49 already frozen
PARI coordinate charts, checking every quartic-square identity and the
completed height-100,000 transcript.

None of these published basis representatives lies at height at most 100,000
or 1,000,000; the smallest height across all points, signs and charts is
1,472,464. The blind control already certifies 27 using other points, and
its union with the published subgroup certifies 29. Thus direct visibility
of a particular published basis is an inadequate proxy for visibility of
its Mordell–Weil cosets. This audit supplies no gate for a tenfold height
campaign. It does not exclude smaller translated representatives, unrecorded
points, or gains from a different chart policy. No oracle point enters a
prospective search. The remaining useful diagnostic is reduction modulo the
blindly recovered subgroup, with that separation maintained.

The first visibility builder succeeds, but its checker compares in-memory
tuples to JSON lists and rejects. Version 2 normalizes that comparison and
passes both build and replay; the first certificate and failed replay log
remain retained. No arithmetic failure is relabelled as a success.

## Bounded follow-up on the new rank-27 curve

The separately frozen 301-centre adaptive follow-up uses only this curve's
own 27 discovered points. All 301 boxes complete at height 100,000; every
archived admission history replays. The initial and adaptive outputs contain
**915 distinct points up to sign**, whose finite quotients still certify
27 modulo 2, 3 and 5. The
[coverage certificate](../../artifacts/generated-results/elliptic-curves/next24_rank27_adaptive_coverage_v1.json)
records that bounded outcome. No 28th direction was established, and this is
not an upper bound, saturation proof, or absence theorem.

## Reproduction

```sh
python3 elliptic-curves/cas/certify_next24_r17_results.py --check
python3 elliptic-curves/cas/certify_next24_high_rank_minimal_v4.py --check
python3 elliptic-curves/cas/export_new_high_rank_curve_index_v5.py \
  --check artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v5.json
python3 elliptic-curves/cas/audit_native11952_published_visibility_v2.py --check
```

`MATH_STATUS.json` remains authoritative. The new lower bounds are exact;
search completeness applies only to the declared finite boxes.

The [portable evidence supplement](../../artifacts/generated-results/elliptic-curves/next24_discovery_evidence_v1.json)
names eight pinned base archives and preserves all new search histories,
failures and replay inputs. Its
[isolated verifier](../cas/verify_next24_discovery_bundle.py) performs 32
bounded checks, including exact geometry and point provenance for all 1,381
new charts, the 62-curve inventory and CSV, all finite point proofs, local
minimality regressions, and the separate published-point visibility audit.
It repeats no point search. Admission/archive histories were replayed locally
and are retained; that distinct check is not repeated by the isolated verifier.
The [portable replay outcome](../../artifacts/generated-results/elliptic-curves/next24_portable_replay_v1.json)
records the actual stage results.

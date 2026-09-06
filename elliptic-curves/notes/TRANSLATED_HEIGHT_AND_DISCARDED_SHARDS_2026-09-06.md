# Translated visibility, a measured height increase, and eight more curves

The saved-shard experiment adds **eight catalogue-unmatched curves** with
certified lower bounds 23–26. The independently replayed inventory now contains
**70 curves**, including the three earlier rank-at-least-27 examples. Exact
ranks, universal novelty, and new rank-at-least-28/32 curves remain open.

The strongest addition, `new-20260906-63`, is compact family `074d9` at
`2588/1603`. Its proved global minimal equation is

```
y² + x*y + y = x³
  - 24458430718461499382493264906272310854066801448*x
  + 1404034490168971387859499604490268925902658396966107000566121492008506
```

The [Sage file](../../artifacts/generated-results/elliptic-curves/new_discarded_rank26_curve.sage)
contains all 26 independent rational points. The
[minimal-model certificate](../../artifacts/generated-results/elliptic-curves/discarded_rank26_minimal_proof_v1.json)
replays their independence, exact transport and rational-isomorphism comparisons.
Here `gcd(c4,c6)=121`, and the valuations at 11 are `(2,2,4)`, excluding
nonminimality. There is no match among the pinned 586 catalogue equations or
346 earlier measured equations. Catalogue absence is not universal novelty.

The [70-curve JSON](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v6.json)
and [CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v6.csv)
preserve all 62 earlier IDs. The exact lower-bound buckets are 3 at 27, 5 at 26,
14 at 25, 14 at 24, 17 at 23, and 17 at 22. All 70 have distinct j-invariants.
Cross-family incidence remains checked for the earlier 47 curves; the latest
23 additions have not yet received that audit.

## The missed control directions become visible after translation

The [previous direct visibility audit](NEXT24_RANK27_DISCOVERY_2026-09-06.md)
found no supplied published basis representative below height one million
in the native rank-29 control's 49 fixed charts. This was a statement about
those representatives, not their quotient classes.

The exact union certificate identifies published points 2 and 27, with zero-based
indices, as independent of the blindly recovered rank-27 subgroup. A separate
retrospective audit uses a 384-bit numerical height Gram to propose one translate
per missing direction, sign and existing chart, aiming near half the pointed
centre. A unimodular LLL change and floating CVP choose integer group words;
neither numerical optimality nor a complete coset search is claimed.

All **196 translations** replay with independent rational group arithmetic and
exact pointed-quartic square identities. Their best representatives require
heights **113,933** and **918,522**, rather than the original approximately
5.54 trillion and 3.14 trillion. The first is in chart 20 at coordinate
`113933/373`; the second is in chart 12 at `156623/918522`.
The [exact translation replay](../../artifacts/generated-results/elliptic-curves/native11952_translated_visibility_replay_v1.json)
retains every integer word and rational point. These are oracle-only outputs;
they enter neither a prospective worker nor its centre selection.

## A matched control distinguishes useful height from censored work

Both matched arms start afresh from generic rank 17, using the same frozen
49 maps and the PARI point engine. They do not receive the public points or
translated words. All chart outputs and archived admissions replay.

| Height | Completed boxes | Timed-out boxes | Certified bound from that arm |
|---:|---:|---:|---:|
| 100,000 | 49 | 0 | 27 |
| 125,000 | 49 | 0 | 28 |
| 1,000,000 | 0 | 49 | 17 generic points only |

The first and second workers take approximately 37.1 and 46.2 seconds. The
million-height arm reaches its ten-second per-chart cap on every chart and
finishes its declared attempts in 509.4 seconds. Its retained generic bound
17 is **not a rank loss** or a point-search exclusion. No larger-height
completeness is claimed. The 125,000 arm was frozen separately after the
first ten million-height charts timed out: its bound is the next 25,000
increment beyond the exact 113,933 witness, and it changes no prior protocol.

The [height-control certificate](../../artifacts/generated-results/elliptic-curves/native29_height_control_v2.json)
recomputes the finite point proofs and compares exact maps and coverage.
Version 2 checks supervisor logs relative to the extracted workspace; version 1
followed their recorded absolute paths and remains retained as a local checker.
Recovering 28 on this known rank-29 curve is a detector calibration, not a
new rank record or proof that 125,000 suffices on other curves.

A separately frozen prospective follow-up then applies 125,000 to the original
141 generic maps of all three newly found rank-27 curves, starting from their
own certified 27-point subgroups. All boxes complete and all histories replay.
Their old-plus-new clouds contain 1,858, 1,815 and 940 points up to sign and
still certify 27 modulo 2, 3 and 5. The
[coverage certificate](../../artifacts/generated-results/elliptic-curves/new27_height125_followup_coverage_v1.json)
records these finite misses without asserting exact ranks or saturation.

## Reuse the candidates discarded by the earlier merge

The old H4096 scan covered 20,400,078 signed primitive parameters per family,
122,400,468 addresses altogether. It saved 128 finalists per sign, then kept
only 128 overall per family under the short-prime score. The first extended
prime experiments therefore revisited only half of the saved addresses.

This new protocol takes exactly the other **768 already saved addresses**,
verifies their original scanner logs and hashes, and computes the same
4,591,104 additional prime traces. It completes in 126.4 seconds; all trace
rosters and score extraction replay. No new broad parameter enumeration occurs.
Selection uses primes through 32749, while 32771–65521 remain outside selection.
The sign-shard population is still truncated; this does not recover every
candidate excluded by the original short-prime scan.

The best two per family are frozen before point measurements or catalogue
comparison. All twelve start from exactly 17 generic sections and receive the
same 125,000 height and ten-second per-chart cap, justified by the separate
completed control. All **540 boxes** complete. Every archived admission
history and complete point cloud replays. The
[batch certificate](../../artifacts/generated-results/elliptic-curves/discarded12_r17_results_v1.json)
finds no catalogue or earlier-equation match among the twelve; eight have
certified lower bounds at least 22 and enter the inventory. The others retain
lower bounds 17, 17, 19 and 21. Those lower bounds are not exact ranks.

One discarded candidate's extended score exceeds every previously retained
candidate in its family, demonstrating a real selection omission. A score is
not a rank theorem. Because this batch also increases point height, its yield
is not a clean causal estimate of changing the selector alone. No validation
score or public exceptional point influences candidate ordering.

## Reproduction

```sh
python3 elliptic-curves/cas/certify_native29_height_control_v2.py --check
python3 elliptic-curves/cas/extend_r17_discarded_shard_scores.py replay
python3 elliptic-curves/cas/certify_discarded12_r17_results.py --check
python3 elliptic-curves/cas/certify_discarded_rank26_minimal.py --check
```

The [portable supplement](../../artifacts/generated-results/elliptic-curves/height_and_discarded_evidence_v1.json)
retains the fixed experiments, censored outputs and replay inputs, with pinned
base archives. Its [isolated verifier](../cas/verify_height_and_discarded_bundle.py)
checks exact geometry and raw-point provenance for 828 charts, the 196 oracle
translations, standalone rank proofs and the 70-curve inventory/CSV. It repeats
no point search; separately passed local admission histories remain retained.
The [replay outcome](../../artifacts/generated-results/elliptic-curves/height_and_discarded_portable_replay_v1.json)
records actual stage results. `MATH_STATUS.json` remains authoritative.

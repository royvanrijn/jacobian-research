# Further native M27 short passes: historical boundary

The complete native M27 short-pass roster, including maps, parent banks,
controller records, receipts, and cursors, is preserved byte-for-byte in
[the archive](../../archive/elliptic-curves/notes/HIGH_RANK_SHORT_PASS_ROSTER_2026-09-09.md.txt)
(`sha256: 8131b97d0094d70034c131b5d2f7e686dac6ce504fe8f3bb45bf1693e71c3297`).
It does not authorize a new V3 pass or continuation.

All completed point-search branches retained their certified M27 lower bounds:

- Curve40's first attempt stopped during map construction before any point
  call. Its bounded-map replacement completed 100 calls with no gain and six
  map-censored boxes; a separate six-chart preconditioned policy completed
  without gain. The raw minimized boxes and remaining suffix stay `UNKNOWN`.
- Curve71 completed a 100-call productive-parent pass and a separate 100-call
  complementary-parent pass, both without timeout or gain. Their unvisited
  selected centres remain `UNKNOWN`.
- Curve48 completed 300 productive-parent calls and 600 complementary-parent
  calls, all without timeout or gain. Those are separate bounded policies;
  their untouched suffixes remain `UNKNOWN`.

The [Curve40 bounded-map result](../../artifacts/generated-results/elliptic-curves/curve40_bounded_maps_v3_v1/result.json),
[Curve71 preconditioned result](../../artifacts/generated-results/elliptic-curves/curve71_preconditioned_v3_v1/result.json),
and [Curve48 complement result](../../artifacts/generated-results/elliptic-curves/curve48_complement_cached_v3_v5/result.json)
retain the outcomes. A raw map timeout is not a point-search miss; a completed
alternative coordinate box adds only that declared exposure. The full archive
preserves the Curve52 calibration in which 2-and-3-only minimization moved a
known gain outside the height bound, so it was not adopted as a replacement.

Reuse the [bounded-map rule](../../knowledge/ALGORITHMS.md#method-ec-lean-maps-bound-map-construction-separately-from-point-search),
[box-deduplication rule](../../knowledge/ALGORITHMS.md#method-ec-box-dedup-check-actual-bounded-box-equivalence-before-another-point-call),
and [parent-span rule](../../knowledge/ALGORITHMS.md#method-ec-parent-span-diversify-outside-the-earlier-productive-parent-span)
only in separately scoped work selected by the current [programme map](../README.md).

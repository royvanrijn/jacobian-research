# Curve92 productive-anchor search: historical boundary

The complete Curve92 V3 experiments, parent construction, maps, receipts,
controller histories, and terminal cursors are preserved byte-for-byte in
[the archive](../../archive/elliptic-curves/notes/CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md.txt)
(`sha256: df031d6b669c1f8da1e9875a81cd3a58d2012516c107c3d5282f5714375c1fcb`).
They do not authorize another V3 pass or continuation.

The initial three-anchor policy completed all 312 planned calls at certified
M26 without timeout or gain. A pairwise-parent policy exhausted its selected
146-centre/two-map domain after 292 calls with no gain; a 16-centre
preconditioned-map audit found only already completed height-125,000 boxes.
An outside-span complementary-parent branch then completed 400 calls without
timeout or gain. Its remaining 508 selected centres are untested and remain
`UNKNOWN`; none of these bounded misses proves exact rank, rules out a 27th
point, or establishes a conductor record.

The [initial result](../../artifacts/generated-results/elliptic-curves/curve92_productive_v3_v1/result.json),
[exhausted pairwise result](../../artifacts/generated-results/elliptic-curves/curve92_pairwise_cached_v3_v2/result.json),
and [third complement continuation](../../artifacts/generated-results/elliptic-curves/curve92_complement_cached_v3_v3/result.json)
retain the finite evidence. The old conductor and root-number observations are
scheduling diagnostics only. Reuse the
[box-deduplication](../../knowledge/ALGORITHMS.md#method-ec-box-dedup-check-actual-bounded-box-equivalence-before-another-point-call)
and [parent-span](../../knowledge/ALGORITHMS.md#method-ec-parent-span-diversify-outside-the-earlier-productive-parent-span)
methods only in separately scoped work selected by the [programme map](../README.md).

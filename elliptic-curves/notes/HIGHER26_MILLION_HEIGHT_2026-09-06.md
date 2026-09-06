# Million-height follow-up on the new higher-parameter rank26 curve

The new curve at `11952`, `7460/32309`, completes all 49 original generic
charts at height 1,000,000 with a sixty-second cap per chart. Its independent
admission replay still certifies 26. The complete union of the initial 49,
both adaptive 301-chart waves, and these 49 larger boxes contains 1,233 points
up to sign. Exact checks modulo 2, 3 and 5 each give lower bound 26.

The [coverage and point-union certificate](../../artifacts/generated-results/elliptic-curves/higher26_million_higher_11952_coverage_v1.json)
binds the exact chart geometry, raw point witnesses, all 700 chart inputs and
the full-cloud proofs. All 49 larger boxes completed; no censored interval is
reported as covered. The worker and admission replay together took 1,466.654
seconds, with one point worker. Full-cloud construction/checks took 11.349
seconds; the separate geometry and provenance check took 0.481 seconds.

The protocol was frozen before execution. This particular curve had already
gained a 26th direction in its first adaptive wave; its second new-direction
wave completed without gain. The known29 retrospective control recovered its
29th direction in a fixed million-height chart within 39 seconds, supplying a
finite cost and visibility gate. The earlier million-height pilot on a
different new27 curve was null and remains an explicit negative control.
None of these observations guaranteed a further point on this curve.

The seed is the curve's own 26-point proof. The original 49 generic maps,
their coordinates and order are unchanged. There is no public point, new
centre geometry, refill or automatic retry. A provisional bound28 would have
stopped the worker pending replay; that endpoint was not reached.

The curve remains `new-20260906-99` in the 100-curve inventory. This bounded
follow-up establishes no new27/28/32 result, exact rank, saturation or absence
of additional rational points. The [minimal equation and point proof](HIGHER_PARAMETER_RANK26_2026-09-06.md)
remain canonical for its certified bound26.

Replay the recorded admission history with
`higher26_million_height_pilot.py replay --index 0`; replay exact geometry and
full-cloud provenance with `audit_higher26_million_height.py --kind million
--check`. The coverage certificate links both standalone quotient checkers.

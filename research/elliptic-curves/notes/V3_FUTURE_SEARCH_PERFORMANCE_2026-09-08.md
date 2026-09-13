# V3 implementation methods and historical benchmarks

The dated V3 runs are historical calibration evidence. They completed or
stopped under their recorded bounds and do not authorize a restart, larger box,
or new campaign. Current work is selected only by the
[programme map](../README.md).

## Reusable implementation rules

- Reuse exact arithmetic contexts, finite-reduction signatures and immutable
  receipts through the [shared runtime](SHARED_RESEARCH_RUNTIME.md).
- Keep lean map construction separate from point search. The recorded cypari2
  map-worker improvement is a fixed component comparison, not an end-to-end
  search speed theorem.
- A cached continuation reuses only a bound, independently verified landscape;
  it must replay inherited map identities and new point witnesses. A changed
  basis requires a new landscape and rational-CVP verification.
- Deduplicate actual equal bounded boxes, including infinity, before another
  point call. A coordinate presentation alone does not establish new exposure.
- Exact finite-quotient admission establishes lower bounds; a failed finite
  column remains `UNKNOWN`.

[ALGORITHMS.md](../../knowledge/ALGORITHMS.md) contains the scoped methods and
current source links. The [next-direction benchmark](NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md)
and [seed history](SEED_AND_AMPLIFICATION_HISTORY_2026-09-13.md) give the
current policy boundary: cheaper known-control recovery and bounded no-gain
results do not justify scaling a stalled route.

The full historical preparation, timing tables, frozen run protocols and
replay commands are preserved unchanged in the
[archive](../../archive/elliptic-curves/notes/V3_FUTURE_SEARCH_PERFORMANCE_2026-09-08.md.txt).

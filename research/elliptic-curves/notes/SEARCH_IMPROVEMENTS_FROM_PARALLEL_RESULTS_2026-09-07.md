# Retained norm preflight — canonical boundary

The complete 7 September 2026 chronological review is preserved byte-for-byte
in [the archive](../../archive/elliptic-curves/notes/SEARCH_IMPROVEMENTS_FROM_PARALLEL_RESULTS_2026-09-07.md.txt)
(`sha256: f6d14f2c292f0509b70d3ba29840e305aca7cd9a6458b69d205bf221c9b9e7f3`).
It does not authorize a score change, parameter scan, point search, class-group
campaign, or automatic follow-up. Mathematical status remains in
[MATH_STATUS.json](../../MATH_STATUS.json).

## Exact preflight

For a finite dictionary in a separable monic rational cubic, recompute each
norm exactly. Remove a coefficient only when its norm has an isolated
nonsquare remainder outside the forbidden support; peel again after removals
and retain generators with shared support because their product can cancel.
For elliptic Kummer use, the forbidden support must include every bad elliptic
place. A nonempty residual dictionary remains `UNKNOWN`.

The frozen retrospective protocol covers all 132 original-panel and 296
fixed-box generators. It proves 428 forced-zero coefficients, zero remaining
dictionary capacity, and no curve or rank exclusion. The
[result](../../artifacts/generated-results/elliptic-curves/retained_norm_preflight_v1.json),
[independent Sage replay](../../artifacts/generated-results/elliptic-curves/retained_norm_preflight_sage_v1.json),
and [`audit_retained_norm_preflight.py`](../cas/audit_retained_norm_preflight.py)
retain all exact witnesses. The preflight excludes these constructed
dictionaries only; the known-soluble reference confirms that genuine additional
directions can exist despite a zero result for a particular dictionary.

The separate [public28 integration](../../artifacts/generated-results/elliptic-curves/search_result_integration_v1.json)
joins an exact equation and finite-rank replay to the existing `m=17`, `k=0`,
`a=4` boundary for curve188, giving a necessary additional strict rational
dimension of at least seven. It is outcome-derived, not a prospective feature,
class-group bound, or exact-rank result.

Use [ALGORITHMS.md](../../knowledge/ALGORITHMS.md) and
[`norm_ramification.py`](../cas/research_runtime/norm_ramification.py) before
any separately scoped class construction. The detailed review, its dated
construction restrictions and all historical links remain in the archive.

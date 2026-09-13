# Historical elliptic-curve machinery audit — 5 September 2026

This is a navigation record for the September 5 audit, whose roster counts,
prospective outcomes, operational commands and programme gaps are historical.
The exact original is retained in
[`archive/elliptic-curves/notes/ELLIPTIC_BREAKTHROUGH_AUDIT_2026-09-05.md.txt`](../../archive/elliptic-curves/notes/ELLIPTIC_BREAKTHROUGH_AUDIT_2026-09-05.md.txt)
(SHA-256 `779c00babae7397b6e8f58e34dc29d907f905f0e6f91f5aee84c6ebf0be2e40d`).

Mathematical status is now solely in [MATH_STATUS.json](../../MATH_STATUS.json)
and its canonical sources. The active programme map is
[README.md](../README.md); the current cleanup boundary is
[CLEANUP_REPORT_2026-09-13.md](../CLEANUP_REPORT_2026-09-13.md).

## Preserved implementation lesson

Bind a resumable pointed-chart checkpoint to its protocol and complete immutable
`MWState`, including observations, rather than only its curve, basis, centre
and budget. A new observation can leave chart coordinates unchanged while
changing the state a replay must verify. Separate state keys must therefore
receive separate checkpoints; an unmatched state is a cache miss, never a
reusable result. The regression covers two such states without repeating a
search. The shared runtime documents its broader cache and replay contracts in
[SHARED_RESEARCH_RUNTIME.md](SHARED_RESEARCH_RUNTIME.md).

The historical audit also records old bounded selection and cache measurements.
They are evidence for their original protocols only and do not authorize a
restart, extend the current candidate queue, or establish a rank bound.

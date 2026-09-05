# Shared research runtime

Use `ArithmeticContext` / `TwoTorsionContext` for exact arithmetic identity,
`MWState` for subgroup transitions, and the common supervisor for workers.
Discovery fills immutable caches; replay requires retained witnesses and never
silently regenerates missing arithmetic. Heuristic policies do not prove bounds.

`VoronoiIterator` orders cosets from shallow to deep. Maximum-depth selection
uses the separate [`deep_centres`](deep_centres.py) finite audit and
[MW18 calibration](../../notes/MW18_DEEP_CENTRE_CALIBRATION_2026-09-05.md).

[Contracts, connected callers, controls and measured policy sweep](../../notes/SHARED_RESEARCH_RUNTIME.md).

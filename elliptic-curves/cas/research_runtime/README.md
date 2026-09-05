# Shared research runtime

Use `ArithmeticContext` / `TwoTorsionContext` for exact arithmetic identity,
`MWState` for subgroup transitions, and the common supervisor for workers.
Discovery fills immutable caches; replay requires retained witnesses and never
silently regenerates missing arithmetic. Heuristic policies do not prove bounds.

[Contracts, connected callers, controls and measured policy sweep](../../notes/SHARED_RESEARCH_RUNTIME.md).

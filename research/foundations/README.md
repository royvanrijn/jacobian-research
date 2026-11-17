# Foundational proof library

Earlier source proofs, exact checkers, formal projects and failed routes supporting the current [research topics](../README.md) live here. Published manuscripts have moved to [results](../../results/README.md); drafts live in [unpublished](../../unpublished/README.md).

| Collection | Scope and status |
|---|---|
| [Verified core](verified/README.md) | Scoped Keller, cancellation and arithmetic results; global extensions remain open. |
| [Geometry sources](extended-geometry/README.md) | GVC, Gaussian moments and Image-space proofs, with partial frontiers retained. |
| [Plane Jacobian](plane-jc/README.md) | Earlier plane reductions and closed branches; the full conjecture is open. |
| [Formal projects](formal/README.md) | Exact Lean scopes and remaining bridges. |
| [Arithmetic](arithmetic/) | Construction specifications and certificates. |
| [Checkers](scripts/) | Explicit replay programs, including bounded controls. |
| [Hessian catalogue](../index/hessian.md) | Canonical HC4 and Schur source notes in this directory. |

[MATH_STATUS.json](../MATH_STATUS.json) is the sole current claim authority. Historical notes and outputs do not change a claim's state or certify external review.

For a selected retained proof target, inspect [Makefile](Makefile) and run `make -C research/foundations TARGET`. Default repository maintenance is described in the [replay guide](../REPRODUCE.md).

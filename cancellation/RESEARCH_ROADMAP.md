# Cancellation continuation pointers

This file no longer maintains an independent roadmap.  The sole continuation
queue is generated in [`STATUS.md`](../STATUS.md) from
[`MATH_STATUS.json`](../MATH_STATUS.json).

The primary cancellation continuation is `OP-CR`; cancellation arithmetic is
parked as `OP-ARITH`, and controlled-boundary suspension classification is
parked as `OP-SUSP`.  The former roadmap is preserved in
[`archive/legacy-notes`](../archive/legacy-notes/CANCELLATION_RESEARCH_ROADMAP_2026-07-22.md).

The geometric deck group of the cancellation inverse polynomial is treated
separately in [`INVERSE_MONODROMY.md`](INVERSE_MONODROMY.md).  Its generic
`A_N/S_N` theorem and branch-cycle search are independent of the parameter-
polynomial Galois problem in `OP-ARITH`.

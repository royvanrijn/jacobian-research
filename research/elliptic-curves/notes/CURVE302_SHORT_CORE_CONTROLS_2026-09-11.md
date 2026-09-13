# Curve302 historical short-core control proposal

This 128-panel short-vector null model was **not run**.  Its synthetic tests
exercise input validation, exact lattice admission, random-stream isolation and
censoring behavior; they do not supply a Curve302 control result.  It has no
claim in [`MATH_STATUS.json`](../../MATH_STATUS.json), no generated certificate,
and is not a current computation queue.

The completed [short-vector core experiment](CURVE302_SHORT_VECTOR_CORE_EXPERIMENTS_2026-09-11.md)
is the canonical source for the finite 1,288,441-direction enumeration and its
partial match to the observed common-core chain.  The unexecuted design proposed
three distinct checks: literal first-29 lattice generation, deficit-one basin
summaries, and a conditional short-vector reference ensemble.  None may be
reported as measured, calibrated, or predictive.

The full proposal, its frozen bands, primitive-extension rule, censorship policy
and historical command surface are preserved byte-for-byte in the
[archive](../../archive/elliptic-curves/notes/CURVE302_SHORT_CORE_CONTROLS_2026-09-11.md.txt).
Its runner and synthetic regressions remain retained at
[`run_curve302_short_core_controls.py`](../cas/run_curve302_short_core_controls.py)
and [`test_curve302_short_core_controls.py`](../tests/test_curve302_short_core_controls.py).
They require sealed closure and short-vector result directories that are not
available as a complete replay bundle in this checkout.

Do not launch or extend this null model merely because its code exists.  Any
future use needs a separately scoped statistical question, immutable inputs,
an explicit treatment of censoring and an evidence publication plan.  New
elliptic-curve work follows the [current programme map](../README.md).

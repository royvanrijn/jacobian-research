# Curve302 historical chart-exposure control proposal

This retrospective chart-exposure design was **not run**.  It has no
[`MATH_STATUS.json`](../../MATH_STATUS.json) claim, no normalized chart ledger,
no generated output, and no current work authorization.  Its fail-closed
requirements are retained: a chart needs an explicit quotient word for positive
exposure, and a negative comparison needs complete chart coverage; point-only
transcripts and multi-gain ordering remain `UNKNOWN`.

The existing [historical replay adapter](CURVE302_CHART_REPLAY_ADAPTER_2026-09-11.md)
and [canonical integral-core proof](CURVE302_EXACT_SHARED_CORE_2026-09-11.md)
are the relevant completed evidence.  They do not provide the missing
chart-exposure ledger or establish this proposal's multiplicity or scheduling
comparisons.

The full candidate definition, batch semantics, counterfactual boundary,
schema and command surface are preserved byte-for-byte in the
[archive](../../archive/elliptic-curves/notes/CURVE302_CHART_EXPOSURE_CONTROLS_2026-09-11.md.txt).
The retained [`run_curve302_chart_exposure.py`](../cas/run_curve302_chart_exposure.py)
does not make the experiment live.  A future use requires a separately scoped
question, immutable complete transcripts or an exact normalized ledger, explicit
censoring treatment, and a portable evidence plan.  Do not infer zero exposure
from the absent inputs or launch the proposal by default.

# Curve302 historical closure follow-up

This completed 14-run, 180-acquisition retrospective analysis has no claim in
[`MATH_STATUS.json`](../../MATH_STATUS.json) and is not a current experiment.
The canonical finite result is the [integral shared-core proof](CURVE302_EXACT_SHARED_CORE_2026-09-11.md).
The full protocol, commands, execution history and frozen-source hashes are
preserved byte-for-byte in the [archive](../../archive/elliptic-curves/notes/CURVE302_CLOSURE_FOLLOWUP_2026-09-11.md.txt).

## Retained findings

- In the supplied fourteen-axis quotient atlas, seven singleton seeds lowered
  the full-closure threshold.  The best two reduced its recorded value from
  30.94463075 to 25.63565575 (17.1564%).  This is a weighted-hypergraph result
  for that fixed atlas, not intrinsic rational-direction closure.
- Every run filled its declared local quotient by dimensions 4--7.  The
  rank-29 exception ended with strict/local dimensions (8,4); 125 of the 128
  recorded strict increments came from individually mixed vectors.  These
  observations motivated the later integral-core audit; they do not give a
  propagation theorem or a rank bound.
- Two bounded literal-vector vocabularies did not make the tested adaptive
  residual or unlock scores beat static quotient norm.  Their coverage was only
  49/180 and 70/180 respectively, and absent targets were charged in the loss.
  The result is a calibration failure for those scores and vocabularies, not a
  prospective predictor comparison or an assertion about all possible points.

## Evidence and replay boundary

The runner and its synthetic regressions remain at
[`run_curve302_closure_followup.py`](../cas/run_curve302_closure_followup.py)
and [`test_curve302_closure_followup.py`](../tests/test_curve302_closure_followup.py).
They require the earlier sealed `closure-laws.json`, `quotient-relations.json`,
`trajectories.json`, and `REPORT.json` bundle.  That bundle and the completed
follow-up output folder are ignored local artifacts and are absent from this
checkout; the generated landscape input alone is insufficient.  Therefore a
new run or check is not a reproducible action here.  If the full immutable
bundle is recovered for a separately scoped integrity purpose, the runner can
deterministically recompute the three retrospective summaries, but it still
does not repeat the original elliptic-curve coordinate proofs.

Use the [current programme map](../README.md) for new work.  Do not reopen this
analysis merely to retune a score, enlarge its vocabulary, or re-run its
historical controller surface.

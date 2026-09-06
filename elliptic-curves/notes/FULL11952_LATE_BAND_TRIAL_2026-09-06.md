# A second-prime-band selection trial on4096 fresh equations

**Active bounded experiment. No new rank result is certified here yet.**

The completed [full11952 cohort](FULL11952_RETENTION_TRIAL_2026-09-06.md)
has scalar traces in a32771..65521 band that did not enter its selection.
Combining this band with the original score moves the known29 control from
position29 to1 within those64 selected curves. Five of the top six combined
entries have measured lower bounds at least19, versus one among the first
six original entries. This is retrospective evidence on an already selected
cohort, not a rank predictor, incidence estimate or calibrated sensitivity.

The new trial freezes4096 previously unmeasured equations from the existing
retained population. It excludes all853 prior cohort equations and593 pinned
catalogue equations by exact rational isomorphism, then excludes duplicates
within the new roster. The first4096 eligible entries occur by old extended
score position4194, after98 exclusions. The available prefix was fixed at
32768 rows; no new parameter scan, prefix enlargement or refill is allowed.
This4096 roster is an explicit further truncation, so the trial does not
test the entire million-row population with the second band.

For every frozen equation:

1. Compute fresh scalar traces4099..65521 and require agreement with the
   earlier cached4099..32749 score and good-prime count.
2. Independently check the character sums at4099 and32771.
3. Select64 by the full quantized score through65521, then good-prime count,
   denominator and signed numerator. Public ranks and points do not enter.
4. After selection freezes, compute wholly disjoint65537..131071 validation
   and three direct character sums per finalist. These values do not alter
   selection or geometry.

The first eight predetermined rows pass the runtime gate, projecting
514.921 seconds serial against a2400-second gate. They are reused once.
The scoring run has two workers,20 seconds per curve, a1800-second outer
cap and checkpoints every16 rows. Complete score replay and an additional
equation-roster/selection replay precede any point protocol.

The separately frozen point controller allows exactly64 generic17-only
attempts, at most3136 initial charts, height125000 and ten seconds per chart.
All64 exact map files must precede every point search. Two point workers and
two independent proof workers run within fixed caps. A provisional28 stops
that curve pending proof; there is no adaptive wave. Every terminal history,
rational map and retained point cloud is checked, including modulo3,5.
Novelty, model export and inventory promotion require their separate proofs.

Frozen sources are `../cas/full11952_late_band_selection.py`,
`../cas/verify_full11952_late_band_selection.py`,
`../cas/full11952_late64_r17_pari_batch.py` and
`../cas/finish_full11952_late64_points.py`.
Protocols, ledgers and raw calls are under
`artifacts/local/elliptic-curves/full11952-late-band-selection-v1`,
`full11952-late64-controller-v1` and `full11952-late64-r17-pari-v1`.

One premature point-controller startup is preserved separately: it read the
protocol before its setup command finished and exited before any stage ran.
The prepared execution begins only after setup completed. No point or trace
attempt was repeated by this startup correction.

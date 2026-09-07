# Calibrated point exposure at record-scale arithmetic height

The recovered determinant1092 MW17 parent is the main search family for this
campaign. The [302 calibration](CURVE302_RECOVERED_SUBGROUP_CALIBRATION_2026-09-07.md)
recovers seven of fourteen known exceptional directions before the prospective
intake starts. Selection is complete; point exposure is active and has no
terminal outcome yet. No new near-record curve is claimed.

The [immutable intake package](../../artifacts/generated-results/elliptic-curves/det1092_record_scale_intake_v1/manifest.json)
preserves the frozen protocol, full projective trace tables, population block
hashes,2048 retained score summaries, final48 choices, and completed stage ledger.
SHA-generated primitive fractions give1048576 distinct parameters with height
65536..1048576, requiring1133153 draws. Every equation has its reduced rational
j-height computed exactly. The admitted640..767-bit numerator bands contain
905827 equations, compared with302's667-bit numerator and the old sibling
selection's256..319-bit band. Arithmetic-height preference is explicitly
requested;302's parameter and exceptional points are not selection inputs.

The166 first-stage primes5..997 use complete projective trace tables checked
independently against PARI. At ambiguous short reductions, exact integer
coefficients undergo repeated p^4/p^6 scaling before good/bad classification.
Synthetic5/13 scaling regressions pass. No restored scaling occurs in this
particular population. The strongest1024 first scores in each64-bit band
receive corrected scalar scoring through32749. Every retained first score
matches its scalar prefix; direct character-sum checks also pass. Validation
primes65537..131071 are neither computed nor read.

Each band contributes sixteen strongest scores, four consecutive rows starting
one-third down its ranking, and four fixed SHA choices from the bottom third.
All48 choices have distinct exact j-invariants and numerator heights644..734bits.
The selector uses no point outcome, rank label, catalogue or point-based refill.
Distinct j-invariants distinguish these48 geometric isomorphism classes;
novelty against external catalogues is a separate post-search question.

Selection costs400.966280993 recorded seconds:5.441532995 for tables,
15.475669559 for the population and380.049078439 for score extension. The
replay checker binds raw score hashes and reconstructs equations, heights and
selection from retained scores; it does not repeat the entire population or
all original character sums.

A subsequent [complete intake audit](../../artifacts/generated-results/elliptic-curves/det1092_record_scale_population_audit_v1.json)
reconstructs every SHA draw, including rejected fractions and deduplication,
then rechecks all1048576 rational j-heights, every admitted first-stage table
score, and both1024-row retention heaps. Direct polynomial evaluation checks
all2048 retained models independently of homogeneous Horner evaluation. It
finds no actual removable local scaling and no selection mismatch. This costs
11.070738775 supervised seconds, separately from the original selection and
point exposure. It reuses the previously independently checked trace tables;
it does not repeat every finite-field point count or add a parameter search.

## Frozen prospective point policy

Every fibre first supplies its seventeen specialized generic sections and
passes exact finite-group independence. Each initial wave uses49 charts from
those sections only. A certified gain permits another49-chart wave using the
enlarged independent subgroup and excluding centre parities in the preceding
input subgroup. The calibrated2048-mask,384-bit metric, largest-norm policy,
factor-free maps and125000-height/ten-second point boxes remain unchanged.

There are at most five waves per curve: one initial and four follow-ups.
A curve stops after a no-gain wave, certified lower bound32, or the fifth wave.
Every selected curve is attempted; there is no replacement or additional
population. The maximum is11760 boxes, not a claim of completed exposure.
One worker uses2GiB RSS,180-second geometry/individual-proof limits and
1200-second worker/history limits. A failed or censored stage halts the driver
with its evidence intact. Exact histories, geometry and independent cloud
proofs modulo2,3,5 precede every accepted gain.

Report discovered directions and completed exposure by score stratum, including
proof cost and any failures. The32/8/8 stratum sizes and two height bands make
this a bounded portfolio comparison, not a universal score law. Adaptive
exposure is identical as a policy but may consume more boxes on gaining curves.

Live checkpoints are in
`artifacts/local/elliptic-curves/det1092-record-scale-points-v1/ledger.json`.
The driver is `elliptic-curves/cas/det1092_record_scale_points.py`.
`report_det1092_record_scale_points.py --snapshot` checks completed proofs and
writes a local progress report with completed boxes, certified gains and stage
costs by height band and score stratum. Its final mode refuses an unfinished
campaign; it never treats unsearched choices as null results.
The completed selection replays with:

```sh
python3 elliptic-curves/cas/report_det1092_record_scale_selection.py --check
```

Keep the frozen campaign sources and gate unchanged while it runs. Update this
note with independently completed outcomes when the fixed48-curve experiment
terminates; intermediate bounds do not stand in for a completed comparison.

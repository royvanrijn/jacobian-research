# Curve 302: one bisection-orbit arm lifts the recovered subgroup to 25

This is a deliberately small, retrospectively calibrated recovery test.  It
does not change the sealed generic17-to24 calibration policy, and it makes no
prospective claim about how to rank rational-bisection orbits.  Mathematical
status remains in [`MATH_STATUS.json`](../../MATH_STATUS.json).

## Exact reconstruction

The cheapest residual strict presentation in the prior diagnostic's declared
32-chart vetted set uses the M17 half-lattice centre

\[
(1,0,-1,2,0,1,2,-1,0,-1,1,0,1,-3,-1,0,1).
\]

Under the explicit M17-to-degree-two LLL transport it is the canonical
minimum-norm representative of rational-bisection orbit `117420`.  Its
generic norm is `10`, and the corresponding degree-two divisor label is
`2O+4F+phi(w)`.  The M24 centre is that word with seven trailing zeroes, so
its M24 parity mask is `95909`.

The exact factor-free reduced quartic has coefficients, in ascending order,

\[
\begin{split}
(&498834529132684723375303442265823300,\\
&1200839569316129875617237773991311840,\\
&498473661155487017705111260021173720,\\
&-250662428412975696824929816386794080,\\
&21255117339473537836188155458220644).
\end{split}
\]

Its bounded box contains the reduced coordinate `2945/581` (with the exact
square witness retained in the transcript).  The full exact map and all
returned points are in the [generated report](../../artifacts/generated-results/elliptic-curves/curve302_m24_bisection_orbit117420_v1.json).

## Why the frozen 17→24 policy omitted it

The old M17 arm first generated `2048` SHA masks and then selected only the
`49` *largest* specialized rounded-height norms.  Parity `95909` was not in
that sampled set, hence could not be selected.  The later recovery waves also
required a nonzero coordinate beyond the previous subgroup.  At the M24 wave
the floor is M22; this padded M17 word has `95909 >> 22 = 0`, so it is
ineligible before deep-centre ranking.

Thus the omission is a direct consequence of both parts of the frozen policy:
initial deep sampling missed the label, and later novelty filtering excludes
every M17-supported label.  It is not arithmetic absence.

## Frozen amendment and result

The smallest possible addition is one new chart: take the predeclared complete
degree-two orbit `117420`, use its exported canonical norm-ten M17 word, and
append zeroes in the seven recovered M24 coordinates.  The one-chart arm
starts from the independently certified M24 cloud and retains the historical
height `125000` and ten-second limit.

The selector was chosen after retrospective inspection, so it is calibrated,
not prospective.  Once frozen, however, the geometry, worker, and exact replay
were artifact-read guarded: their arithmetic inputs were only the M24 cloud,
the M17 parent, the degree-two certificate, and its orbit table.  They did not
read the residual diagnostic, a public missing-point list, or any post-search
rank result.

The one chart completed in under a second.  Its complete returned point cloud
has exact finite-quotient rank `25` modulo `2`, `3`, and `5`; the independently
replayed certificates therefore prove rank at least `25`.  The state-machine
admission remained at24 because it attempts incremental mod-2 admission before
the complete cloud is audited, but the retained complete-cloud mod-2 proof
does certify the new twenty-fifth direction.

## Replay

From `research/`:

```sh
python3 elliptic-curves/cas/audit_recorded_point_mod2_rank_v3.py --check artifacts/generated-results/elliptic-curves/curve302_m24_bisection_orbit117420_mod2_v1.json
python3 elliptic-curves/cas/audit_retained_cloud_modl.py --check artifacts/generated-results/elliptic-curves/curve302_m24_bisection_orbit117420_modl_v1.json
python3 -m unittest elliptic-curves.tests.test_curve302_m24_bisection_orbit117420
```

The producer is
`elliptic-curves/cas/run_curve302_m24_bisection_orbit117420.py`; it refuses to
overwrite the frozen local run or generated certificates.

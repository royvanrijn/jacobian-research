# MW16 sensitivity recovery

Date: 2026-09-05.

This is a complement-blind detector calibration on 398/400/401/542/548,
followed by a separately gated replay of the frozen 104 prospective fibres.
The historical [54/55 ladder](ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md)
and [weight-one pointed sieve](POINTED_QUARTIC_SIEVE.md) remain preserved
baselines. Coordinate presentations within a curve are nested trials, not
independent observations.

The frozen policy recovers **55/55** historical control directions, improving
on the historical 54/55 result. Its initial wave gives `5,5,11,10,8`;
adapting the two controls with initial gains below eight gives
`14,12,11,10,8`. Both adaptive centre constructions separately reach the same
curve-398 and curve-400 ranks. This is calibration recovery, not a new record
curve or an upper bound on any control's rank. Here 55/55 compares certified
quotient dimensions; it does not identify the recovered basis with a public
complement, whose points are not loaded.

The subsequent 104-fibre rerun completes all 856 height-100,000 boxes with
zero finite points and zero certified gains. No prospective gain was found.

## Calibration policy

The sole control fixture supplies the short equation, sixteen specialized
generic sections, and the exact generic Gram. No public exceptional point,
public complement, or per-curve target rank enters search, centre selection,
or coordinate selection. The aggregate threshold of 54 is the declared gate,
not a source of points or a per-curve stopping rule.

The same complete maximum-depth generic half-lattice stratum is used on every
control: respectively 12, 4, 8, 10, and 8 masks. The sensitivity menu varies:

- projective height: 10,000, 100,000, and 200,000;
- centres: generic minimum-norm representatives and representatives shortened
  using the specialized height form;
- horizontal coordinates: the original Gauss metric, weights 1/16 and 16,
  and exact PARI horizontal reduction of the denominator-cleared quartic;
- rational slope parametrization after reduction: `z=s/2`, `z=2s`, and
  `z=s+1/2`;
- adaptive centres: representatives shortened in the recovered point lattice
  and literal zero/one quotient residues.

Every finite enumeration uses the same C++/GMP modular sieve. The horizontal
reduction comparison calls `hyperellred` to choose coordinates; no variant
calls `hyperellratpoints` or `hyperellminimalmodel` for point enumeration.
The selected metric policy calls neither quartic minimization nor reduction.

For integral `Q=(a/d^2,b/d^3)` on `y^2=x^3+A*x+B`, the selected lattice
metric is

```text
N^2 + 16*(|a| + d^2*(floor(sqrt(|A|))+1))*D^2.
```

The original denominator congruence is unchanged. Every horizontal matrix and
ordinate scale is explicit. The checker verifies the five coefficients by
substitution into the original rational pointed quartic, transports both
square-root signs, and checks infinity. Rational parametrizations change the
finite height box; they do not change the rational curve.

The common height-100,000 generic-centre menu is scored by the sum of the five
**certified quotient ranks**. Ties use height and the setting key; runtime is
not an input to selection. Weight 16 wins that menu with initial gains
`5,5,11,10,8`, totalling 39. In particular, the missing curve-401 direction is
already recovered in this initial wave.

| initial configuration | 398 | 400 | 401 | 542 | 548 | total |
|---|---:|---:|---:|---:|---:|---:|
| preserved weight-one, specialized centres, H=10,000 | 0 | 0 | 2 | 10 | 8 | 20 |
| weight-one, specialized centres, H=100,000 | 2 | 2 | 10 | 10 | 8 | 32 |
| horizontal reduction, specialized centres, H=100,000 | 4 | 3 | 10 | 10 | 8 | 35 |
| three rational parametrizations, specialized centres, H=100,000, union | 5 | 4 | 10 | 10 | 8 | 37 |
| horizontal reduction, generic centres, H=100,000 | 5 | 4 | 9 | 10 | 8 | 36 |
| weight 16, generic centres, H=100,000 | 5 | 5 | 11 | 10 | 8 | 39 |
| horizontal reduction, generic centres, H=200,000 | 9 | 6 | 11 | 10 | 8 | 44 |
| horizontal reduction, both centre choices, H=200,000, union | 12 | 7 | 11 | 10 | 8 | 48 |

With the same adaptive trigger, the horizontal-reduction/generic-centre
height-100,000 policy reaches only `14,12,9,10,8` (53 directions). The selected
weight-16 policy reaches 55.

The higher-height comparisons are retained. Selection of the final policy is
judged after adaptive recovery, rather than treating initial rank or execution
time as the final objective.

## Adaptive construction

The frozen uniform trigger is `0 < initial certified quotient gain < 8`.
A gain of zero supplies no quotient direction to lift; gains of eight or more
end this bounded calibration policy. This is a scheduling rule, not an upper
rank bound. It activates curves 398 and 400 from their blind initial results,
without consulting their atlas ranks.

For each active curve, retain at most five already discovered independent directions,
ordered by the bit heights of their rational `x` and `y` coordinates, then by
the point. The original sixteen generic sections remain the basis prefix.
No public point or target rank enters this choice.

Enumerate every nonzero word in this five-bit quotient slice above every
frozen deepest generic mask. Compare the canonical-height CVP representative
with the literal binary residue at height 10,000. The selected metric run
therefore has 372 and 124 charts per centre construction on 398 and 400:
992 bounded searches altogether. The other three controls retain their
completed initial certificates. This is a complete enumeration of the
**declared five-bit slice**, not of all quotient classes when the initial gain
exceeds five. It avoids the historical curve-401 8,184-chart next wave without
claiming to have searched that wave.

Each setting is scored separately, and its union with the initial discoveries
is classified by exact group law. Finite reductions and a no-rational-2-torsion
witness certify point independence. Every reported dependent point has an
exact integral relation in the retained basis. These are recovered-subgroup
certificates, not saturation certificates for all of `E(Q)` or rank upper
bounds.

## Prospective result

Only after the 55-direction control certificate passed did the frozen 104
fibres receive the selected generic-centre, weight-16 policy at height
100,000. All **856/856** deepest-mask boxes complete with **zero finite
points, zero certified quotient gains, and zero timeouts**. Thus no
prospective adaptive wave is triggered and no candidate reaches residual
Selmer or unrestricted continuation.

The declared boxes cover 17,120,085,600,000 integer numerator/denominator pairs
and perform 2,803,616 exact square tests. Integral quartic coefficients have
1,539--1,790 bits. Modular workers use 4,637.80 seconds in total, with a maximum
of 7.31 seconds per chart under the 20-second cap; campaigns use at most eight
workers. These are recorded timings, not a runtime guarantee.

The [prospective summary](../../artifacts/generated-results/elliptic-curves/mw16_sensitivity_prospective_summary_v1.json)
and [self-contained evidence](../../artifacts/generated-results/elliptic-curves/mw16_sensitivity_prospective_v1.json.gz)
retain the exact frozen inputs and calibration gate. Every generic MW16 basis
receives a fresh finite-reduction independence check in the verifier. This
recovers the requested control sensitivity but finds **no prospective gain**;
the null result proves no rank upper bound or absence of further points.

## Evidence and replay

The [control summary](../../artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_summary_v1.json)
and [control evidence bundle](../../artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_v1.json.gz)
retain every completed comparison and the frozen policy. Bundles retain executed source snapshots, seeds, chart identities,
returned points, exact relations, finite-reduction witnesses, limits, and
software identifiers. Separate exact replay uses the bundles and does not require
ignored local checkpoints. The checker recomputes the polynomial identities
and finite signatures but shares arithmetic helpers with the search; this
is not a claim of a fully independent software implementation.
`check_mw16_sensitivity_policy.py` additionally binds every recorded height,
coordinate label, deepest mask and quotient word to the declared complete
policy, then invokes that exact checker.

```bash
python3 -m unittest discover -s elliptic-curves/tests -p test_mw16_sensitivity_backend.py
python3 -m unittest discover -s elliptic-curves/tests -p test_half_lattice_pointed_sieve.py

bash elliptic-curves/cas/replay_mw16_sensitivity.sh
```

Exploratory adaptive calls on the larger-initial-gain controls are retained
under ignored local checkpoints. Their expensive whole-point-set relation
classification was stopped after the 55-direction target had already been
reached. No completion or extra rank claim is made for those interrupted
exploratory ledgers. The published policy is replayed separately with exactly
the two eligible controls.

The shell replay declares and checkpoints the whole calibration before
constructing the prospective gate. The new menu calls have a 20-second limit (15 seconds for the height-only
weight-one baseline);
the height-200,000 comparison uses 60 seconds. Exact relation chunks retain
the 180-second limit. At most eight modular workers run within a campaign.
Completed chart checkpoints are content-addressed; interrupted charts cannot
be promoted to completed boxes. If a campaign records an incomplete setting,
use a fresh campaign output to retry it; completed chart cache entries are
reused and incomplete entries are searched again.

The gate requires at least 54 certified directions on all five controls and a
passing independent replay of the exact source-linked initial/adaptive
pipeline. The prospective checker also requires all 104 frozen candidate IDs
and the same initial centre/coordinate/height policy. A bounded prospective
miss supplies no rank upper bound, structural rejection, or Selmer result.

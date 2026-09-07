# Retained outer and native fibres: exposure and a preparation fix

Four retained curves receive 196 completed point-search boxes, with no newly
certified direction. The two outer curves remain at lower bound26 and the
two native carrier fibres at19. The useful machinery change is a separate
factor-free chart policy: it prepares all49 charts on a 3,875-bit displayed
model in **2.579 seconds**, after the historical pipeline had reached its
300-second preparation cap without one chart. This is preparation feasibility,
not improved rank yield or a new record. No parameter sweep follows.

## Outer26 exposure gap closed

The [source-exposure audit](../../artifacts/generated-results/elliptic-curves/recent_outer26_exposure_audit_v1.json)
verified that the original92 charts on IDs189 and192 used only their generic
17/16-point prefixes. This separate experiment uses each certified26-point
subgroup. Exactly2,048 deterministic parity masks, nonzero above the generic
prefix, are selected using the existing `full11952-specialized-followup-v1`
domain. Canonical heights use384-bit precision and rounding at10^6; unimodular
LLL and numerical CVP propose representatives. Exact parity and rounded norms
are checked. The49 largest computed norms are retained, with both map files
frozen before points. CVP optimality and complete coset coverage are not claimed.

Each chart receives height125000 and ten seconds, one worker, no rank stop.

| Curve | Family and parameter | Completed boxes | Point-worker seconds | Cloud / certified lower bound |
|---|---|---:|---:|---:|
| 189 | 11952, 2618/26913 | 49 | 47.334 | 30 / 26 |
| 192 | MW16-01, −4444/217 | 49 | 44.871 | 58 / 26 |

The [report](../../artifacts/generated-results/elliptic-curves/recent_outer26_followup_v1.json)
binds histories, maps, point provenance and finite proofs modulo2,3,5.
Its supervised stages total116.222 seconds, excluding later Sage and standalone
checks. [Independent Sage](../../artifacts/generated-results/elliptic-curves/recent_outer26_sage_replay_v1.json)
enumerates complete finite groups and quotients by doubles on all88 points.
The [three-file standalone bundle](../../artifacts/generated-results/elliptic-curves/recent_outer26_rank_standalone_v1.zip)
passes in a fresh directory. Neither26 is an exact-rank claim.

## Twelve retained native19 fibres, then two fixed trials

The [native carrier construction](NATIVE_RANK3_CARRIER_IMAGES_2026-09-06.md)
already supplied twelve08234 fibres with17 inherited and two independent native
directions. Their original400-bit model budget gate remains a valid record of
that earlier experiment. It was not a mathematical exclusion.

The [retained score audit](../../artifacts/generated-results/elliptic-curves/retained_native19_scores_v1.json)
evaluates all12 equations at3,510 primes from5 through32749:42,120 traces and
576 independent character-sum checks. Prime-local short-model scaling is tested;
none is needed on this population. The score is the frozen rounded sum of
`(2-a_p) log(p)/(p+1-a_p)` over good primes. Validation primes65537..131071 are
not computed or used. [Exact intake](../../artifacts/generated-results/elliptic-curves/retained_native19_intake_v1.json)
replays all twelve19-point proofs and finds twelve distinct j-invariants,
unmatched against the pinned620 catalogue equations and201 inventory equations.
This is relative novelty, not proof of unpublished curves.

The new bounded trial admits exactly the strongest score, `native19-08`,
word[-2,0], and the smallest model, `native19-01`, word[1,0]. They have3,875
and543 displayed coefficient bits. The smaller fibre's parameter is
−7119612289/2394065174; the full larger parameter is retained in the report.
Their scores are126.642309679919 and124.016921795943. This pair tests feasibility
outside the old height budget; it is not a height-matched score comparison.

Preparation failures are preserved:

- V1 rejects a1.430511474609375e−6 numerical discrepancy at rounded norm
  2183933217 using its absolute1e−6 guard. It spends0.994 supervised seconds
  and attempts no point box.
- The separate V3 geometry helper accepts finite discrepancies within
  `max(1e−6,64 ulps of each compared value)`. It still verifies exact parity
  and rounded norms; four focused tests pass. This is a numerical consistency
  guard, not an error-bound or closest-vector theorem. V2 of the trial uses
  this guard and completes all2,048 representatives, but its combined
  `hyperellminimalmodel` / `hyperellred` call does not produce the first chart
  before the300-second cap. Supervision records301.043 seconds and zero boxes.

## Factor-free mapping, separately calibrated

[`factor_free_pari_mapping.sage`](../cas/factor_free_pari_mapping.sage) constructs
an integral denominator/Gauss quartic, removes only verified square content,
then calls `hyperellred` alone. It verifies the complete rational horizontal
map and a positive rational square multiplier against the original pointed
quartic. Global minimalization is not needed. Old helpers and protocols remain
unchanged. The new matrices define different finite coordinate boxes.

Two fixed geometry probes pass: the pending native chart takes0.064 seconds,
and a historical control chart0.037 seconds. The latter probe's local ID
`known28-chart5` denotes **zero-based index5**, hence the sixth chart; it is not
the historical first-gain chart. A separate point control explicitly uses
zero-based index4, the historical fifth chart, with the old27-point seed only.
Its single125000-height box completes in0.435 seconds and supplies a cloud
of29 points with [certified lower bound28](../../artifacts/generated-results/elliptic-curves/factor_free_known28_control_v1.json).
This retrospectively selected control does not estimate prospective sensitivity.

The control's original wrapper incorrectly demanded independence of its whole
point cloud and failed after search completion. A separate audit selects
independent columns from the retained cloud. Its original CLI comparison also
rejected descriptive timing floats in an exact-only digest; the
[corrected replay entry point](../cas/replay_factor_free_control.py) compares the
complete JSON. Both sources remain preserved, no control point search is repeated,
and the later independent Sage and fresh-directory proofs confirm rank28.

V3 of the native trial freezes all98 new maps before any points, with the same
two curves,19-point seeds,2,048-mask samples,49-centre rule,125000 height and
ten-second per-chart cap. No public points enter the native workers.

| Fibre | Map seconds | Point-worker seconds | Completed boxes | Cloud / certified lower bound |
|---|---:|---:|---:|---:|
| native19-08 | 2.579 | 46.781 | 49 | 19 / 19 |
| native19-01 | 2.016 | 42.736 | 49 | 19 / 19 |

The [completed report](../../artifacts/generated-results/elliptic-curves/retained_native19_trial_v3_results_v1.json)
replays histories, all4,096 sampled exact parity/norm transports,98 rational maps
and point provenance, and both clouds modulo2,3,5. Supervised trial stages total
113.583 seconds, or415.620 including the two failed preparations. This excludes
the separately recorded score/intake, method controls, report generation and
later proofs; it is not an all-session CPU total. Independent Sage checks add
0.891 seconds, and a fresh-directory standalone replay0.944 seconds.

The [four-file standalone bundle](../../artifacts/generated-results/elliptic-curves/retained_native19_rank_standalone_v1.zip)
contains both native certificates, the known28 control certificate and a Sage
verifier, requiring no repository imports. It independently verifies all67
points and lower bounds19,19,28 by full finite-group enumeration. A bounded
null does not exclude further directions or prove that larger native models
are unproductive. The remaining ten retained native fibres are not automatically
scheduled. No high-rank inventory version changes.

## Reproduction

These checks consume retained evidence and launch no elliptic point search:

```sh
python3 elliptic-curves/cas/report_recent_outer26_followup.py --check
python3 elliptic-curves/cas/score_retained_native19.py check
python3 elliptic-curves/cas/audit_retained_native19_intake.py --check
python3 elliptic-curves/cas/replay_factor_free_control.py
python3 elliptic-curves/cas/report_retained_native19_trial_v3.py --check
python3 elliptic-curves/cas/package_retained_native19_proofs.py --check
```

Local protocols, failed attempts, raw outputs and supervision remain under
`artifacts/local/elliptic-curves/recent-outer26-followup-v1`,
`retained-native19-scores-v1`, `retained-native19-point-trial-v1/v2/v3`,
`factor-free-map-probe-v1` and `factor-free-known28-control-v1`.

# The saved own27 charts recover the known28 direction

The previously frozen 49 own27 charts for the11952 fibre110314/102227
**recover rank28 on chart5**, using only the earlier local27-point seed.
All49 charts complete at height125000 with ten seconds per chart and no
rank stop. History replay, whole-cloud proofs modulo2,3,5 and independent
Sage group enumeration pass. This is a known-curve detector calibration:
the curve is ICARM619, already published with28 independent points. There
is no new curve,29th direction, exact-rank or first-discovery claim.

The [point outcome](../../artifacts/generated-results/elliptic-curves/inventory188_fixed49_point_control_v1.json)
records49 retained points and rank28. Point execution takes46.4273 seconds;
the six search/history/cloud stages total55.2589 seconds. The additional
independent Sage rank check takes0.8447 seconds. A
[9,980-byte standalone bundle](../../artifacts/generated-results/elliptic-curves/inventory188_point_search_standalone_v1.zip)
contains the49-point certificate and its Sage checker. Its
[fresh-directory replay](../../artifacts/generated-results/elliptic-curves/inventory188_point_search_standalone_v1.json)
passes with no repository imports or public point file. That small bundle
proves the rank bound; complete search-history provenance remains separate.

## Why the preceding representative tests missed it

The earlier [fixed chart-coverage control](INVENTORY188_CHART_COVERAGE_2026-09-07.md)
located two public extra representatives and one-basis translates. All
tested coordinates were outside125000. Those statements remain exact,
but they do not imply detector failure on their quotient direction.

A separate fixed nearest-translate diagnostic now uses both public
representatives against the old27 subgroup. It computes384-bit height
matrices, rounds at10^6, reduces only the old27 principal block by unimodular
LLL, and projects the negative public point into that block. Floating CVP
proposes a translate. The roster is the original point, that translate,
and its54 neighbours obtained by adding or subtracting one reduced basis
vector. No search-coordinate feedback changes this roster.

For each witness there are56 distinct representatives, both signs and98
saved charts: **21,952 exact coordinate checks**, all outside125000.
The two nearest translates are the same rational point; hence the two
public witnesses represent the same class modulo the old subgroup. Its
best coordinate is−119747429439/44329639748, in own27 chart44, with height
119,747,429,439. The [V2 diagnostic](../../artifacts/generated-results/elliptic-curves/inventory188_nearest_translate_visibility_v2.json)
and [independent rational replay](../../artifacts/generated-results/elliptic-curves/inventory188_nearest_translate_replay_v2.json)
check all112 exact group sums, metric/parity transports and coordinates.
No CVP-optimality or exhaustive representative claim is made.

The first diagnostic worker incorrectly reused a mutable Sage vector for
its neighbours. The independent roster check rejected its two-representative
output. Original sources, output and failure log remain preserved under
`inventory188-nearest-translates-v1`; its `invalid-roster-journal.json`
overrides the worker's premature completion label. V2 copies from a fresh
Python list and checks the55-or56 roster size before coordinates. It uses
the original declared policy and passes independent replay. V1 supplies
no completed-neighbour-exposure claim.

## Exact recovered representative and subgroup identity

The successful point is

```
x = 295093036892844166031650709690689/3468
y = 468599499225192255806797732011201343935008184960/4913.
```

On chart5 it has coordinate **−94237/33087**, height94,237, and is present
in the recorded PARI output and exact square witnesses. The
[independent explanation audit](../../artifacts/generated-results/elliptic-curves/inventory188_recovery_explanation_v1.json)
checks the raw map against the earlier frozen map and finds no
visible-but-unrecorded discrepancy in its98 signed chart observations.

Let P be public point26 and B_0,...,B_26 the old basis in its stored order.
The recovered point is exactly P+sum(w_i B_i), with

```
w = [-1,0,1,-1,0,1,-1,2,-1,-1,-1,-1,-1,0,1,1,0,0,-1,0,0,0,1,0,0,0,0].
```

A bounded height-matrix proposal was accepted only after the Sage group
identity passed. A separate ordinary-Python group-law replay proves the
same identity. The first relation wrapper failed while parsing a decimal
string as a Sage rational, before producing any certificate; its source
and log remain preserved. V2 parses through exact `Fraction` and retains
the original denominator-at-most64 proposal limit. The successful relation
has denominator1 and coefficient1 on P, so the recovered point and P
generate **the same subgroup extension** of the old27 group. Public points
are read only after point execution, to explain the result.

In the same rounded height matrix the recovered representative has norm
56,546,824, while the nearest translate has norm42,364,637. The latter
has a much larger best search-coordinate height. Thus minimizing canonical
height does not by itself minimize the particular coordinate enumerated by
this chart atlas. Both norms are scheduling measurements, not rigorously
enclosed canonical heights or minimum theorems.

## Consequence for further searches

This closes one missing exceptional-direction calibration of the existing
own-subgroup policy. It strengthens the earlier known29 recovery evidence;
it does not establish a universal detection rate. The same policy already
completed49-box trials on all seven currently unmatched27-point curves,
with no gain, so those boxes should not be repeated unchanged.

The useful next candidates are retained high-rank curves whose exceptional
points have not yet entered their chart construction, especially the recent
outer26-point discoveries. That is a different exposure from their generic
initial searches. Their completed histories must be checked before allocating
new boxes. No parameter rescan, score change or automatic next campaign is
part of this control.

The [retained exposure audit](../../artifacts/generated-results/elliptic-curves/recent_outer26_exposure_audit_v1.json) now replays both26-point certificates and all92 source charts for two such candidates:

| Inventory ID | Family / parameter | Parameter height | Generic centre dimension | Known directions beyond generic |
|---|---|---:|---:|---:|
| 189 | 11952, 2618/26913 | 26913 | 17 | 9 |
| 192 | MW16-01, −4444/217 | 4444 | 16 | 10 |

Their source runs use no centre coefficient beyond the generic prefix, even
after admitting new points. Both are unmatched in the pinned620 catalogue.
This verifies a possible change of exposure on retained curves beyond4096;
it does not certify absence of every historical follow-up or promise a gain.
No new point campaign is launched by that audit.

## Replay

```sh
python3 elliptic-curves/cas/verify_inventory188_nearest_translates_v2.py --check
python3 elliptic-curves/cas/inventory188_fixed49_point_control.py check
python3 elliptic-curves/cas/audit_inventory188_recovery_explanation.py --check
```

Full history replay is the point-control script's `replay` mode. The49-point
cloud has separate `audit_recorded_point_mod2_rank_v3.py` and
`audit_retained_cloud_modl.py` checkers. The Sage standalone checker runs
in a fresh extracted copy to preserve its generated output. Local protocols,
logs, raw PARI output and checkpoints are under
`artifacts/local/elliptic-curves/inventory188-fixed49-point-control-v1/`.

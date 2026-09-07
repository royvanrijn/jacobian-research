# Focus on302 and its determinant1092 siblings

The fixed196 point boxes have completed. Final independent cloud and combined
302 certificates are pending; the worker outputs alone are not promoted here.

## Why this parent needs a different test

The user explicitly prioritized302 and comparable good fibres from its
[recovered MW17 parent](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md).
The [completed score-stratum comparison](MIXED_REDUCED_PARENT_POINT_EXPOSURE_2026-09-07.md)
produced four improved sibling lower bounds18/19/19/18, with six new directions
in total. Both strong and moderate score pairs gained three directions;
the moderate pair used less computation. No sibling yet approaches302's31.

There is a substantial unsampled height region. The previous six sibling
selections have reduced j-numerator heights256..319bits. The equation of302
has reduced j-height667/645bits. In the small integral parameter chart its
coordinate is−164518/143797, far outside the previous height8 panel. This
coordinate is a retrospective diagnostic, not an input to new fibre selection.
The parent was chosen using302; its sibling addresses and point searches
remain free of302's exceptional points.

The [invariant quartic bound](BOUNDED_PRIME_EXPOSURE_AND_HEIGHT_FLOOR_2026-09-07.md)
requires at least107bits in any integral completed-square binary quartic for302.
The two chosen siblings have corresponding bounds45 and42bits. Thus the
height gap cannot be erased merely by changing integral quartic coordinates.
It does not itself bound rational point height or rank.

## Frozen experiment

| Row | Supplied independent subgroup | Endpoint |
|---|---:|---|
|302 generic-only control|17|Recover directions absent from the specialized generic subgroup|
|302 direct target|31|Certify any direction beyond the full supplied public31|
|Sibling s=4/3|19|Improve the strongest completed strong-stratum sibling|
|Sibling s=7/2|19|Improve the strongest completed moderate-stratum sibling|

The control exports only the generic equation and17 section formulas from the
parent artifact, specializes at0, and uses the exact short-model map
X=36x+15,Y=108(2y+x+1). The other public directions do not enter its centre
geometry or point admission. Manifest verification hashes other input files;
this is an input-restricted mathematical control, not an OS-isolated blind
process. The earlier full original27-only blind28 regression remains separately
certified and unchanged.

The direct target uses the original31 public points. These are not ordered
with the recovered generic basis as a prefix, so its `generic_points` field
is empty solely to indicate that lack of alignment. The parent still has
generic rank17. Its centre policy samples all nonzero parities of the full
supplied31-dimensional subgroup.

Each sibling starts with its17 specialized sections plus two points found and
independently certified in the completed comparison. Only centre parities
outside the earlier17-dimensional span are admissible. The exact
[centre-class theorem](POINTED_QUARTIC_CENTRE_CLASSES_2026-09-07.md) proves these
degree-two presentations inequivalent to every earlier-subgroup centre under
rational curve automorphisms and rational base changes of degree one. It does
not guarantee new rational points or directions.

All four inputs pass copied-input Sage checks of membership and complete
finite-group independence before the protocol freezes. Each row samples2048
SHA parities, uses384-bit heights rounded at10^6, and takes49 largest computed
norm centres with exact norm/parity transport. All196 factor-free maps precede
every point worker. Point boxes remain height125000 and ten seconds; geometry
has180seconds per row, point processing and history replay1200seconds each,
one worker and2GiB RSS. There is no rank stop, refill, adaptive wave or new
parameter population. Numerical CVP choices carry no optimality claim.

The seed-preparation V1 attempt stopped before any maps or point boxes because
an equality compared serialized coordinate strings with Fractions. V2 converts
both sides to exact rationals. The failed source and0.281173210seconds of
supervision are retained and charged; no mathematical rank failure is hidden.

## Independent endpoints

Each complete cloud receives mod2 and mod3/mod5 construction/replay, exact
geometry replay, and a copied-input standalone finite-group rank proof.
The post-search union separately includes all public31 points, all17 generic
images, and every point returned by both302 runs. This prevents a control
point escaping the public31 span from being overlooked merely because the
control alone detects fewer than31 directions. Agreement of finite ranks
still supplies no upper bound or universal rational-span assertion.

Only after the runs and their proofs finish, a retrospective audit places all
62 signed original public representatives into all49 charts of each302 run.
It checks exact quartic lifting and completed-box visibility. The original
representatives are not all their translates or quotient representatives;
their misses are not a proof that their directions were invisible.

## Reproduction

Run from `research/`. Immutable checkpoints are in
`artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2/`.
The producer is `run_curve302_focused_point_exposure_v2.py`; the independent
combined-cloud supplement is `certify_curve302_focused_union.py`.
Neither command overwrites existing evidence. The final reporter fails closed
until all required proofs finish:

```
python3 elliptic-curves/cas/report_curve302_focused_point_exposure.py --check
```

The [302 lower-bound note](ICARM_CURVE302_RANK31.md) records a public conditional
rank31 statement, whose upper-bound computation is not reproduced locally.
It is not an unconditional exclusion of32. A finite direct attempt is useful,
but302 alone should not consume the whole record-search programme. The next
production decision should depend on improved point exposure at302-sized
arithmetic heights before a larger target-free sibling population is launched.

# Generic17-only recovery on curve302

The [completed calibration](../../artifacts/generated-results/elliptic-curves/curve302_recovery_calibration_v1.json)
recovers **seven of the fourteen known exceptional directions** of302, starting
only from the seventeen generic sections specialized at original parameter0.
The final independently certified subgroup has rank24. This is recovery on a
known curve, not a new curve or rank record.

| Exposure | Input rank | Certified output rank | Completed boxes |
|---|---:|---:|---:|
|Original generic-only control|17|19|49|
|Recovered-subgroup wave1|19|22|49|
|Recovered-subgroup wave2|22|24|49|
|Recovered-subgroup wave3|24|24|49|

The original control is part of the preserved
[focused experiment](CURVE302_FOCUSED_POINT_EXPOSURE_2026-09-07.md).
The follow-up protocol fixes at most four waves before running. It rebuilds
factor-free degree-two charts from the independently certified recovered
subgroup, admitting centre parities outside the preceding input subgroup.
Each wave uses2048 deterministic SHA parity samples,384-bit heights rounded
at10^6, and the49 largest computed norms. Every box has height125000 and a
ten-second search limit. Numerical lattice choices have no optimality claim.
The no-gain third wave triggers the frozen stop; no fourth wave is run.

Geometry, worker and history processes install an artifact-read guard before
research imports. Their allowed inputs are the original17 seed, its recovered19
cloud, and this campaign's own files. Recorded accesses contain no public-span
or public-union files. The public comparison is constructed only after all
point waves and their independent proofs terminate. These input restrictions
do not make an OS sandbox or erase the researcher's knowledge of the control.

All three histories and eighteen cloud certificate stages pass. Exact
independence agrees modulo2,3,5 and in copied-input standalone Sage proofs.
Afterward, numerical heights propose rational expressions for all24 recovered
basis points in the public31 basis; exact rational elliptic-curve group sums
verify every identity. A separate checker uses only those exact sums and
rational matrix ranks. The expressions span dimension24, with the original17
prefix spanning17, hence a seven-dimensional quotient in the known
fourteen-dimensional exceptional space. This counts directions rather than
requiring the original published representatives to occur in a search box.

A [self-contained identity and rank package](../../artifacts/generated-results/elliptic-curves/curve302_recovered_public_span_v1/manifest.json)
retains both point clouds, exact words and standalone Sage checkers. Its basis
lists agree exactly with the identity inputs; the packaged identities and both
rank proofs have also been replayed independently. These extra packaging checks
are recorded separately from the frozen calibration cost below.

The181-point union of public31 and all retained calibration outputs still
certifies31 modulo2,3,5 and independently in Sage. That union test alone is
not an upper bound or proof that every cloud point belongs to the public span;
the exact span identities above concern the recovered24 basis.

Follow-ups and their post-search proofs consume217.990320965 recorded stage
seconds, in addition to88.153354279 for the original generic17 arm. All196
calibration boxes complete. Accepting half the known quotient as the requested
“significant fraction” is an operational decision made with the recovery result
available, not a prospectively registered statistical threshold or a general
sensitivity theorem.

The successful policy has been deployed in the
[record-scale determinant1092 campaign](DET1092_RECORD_SCALE_CAMPAIGN_2026-09-07.md).
Its source roster, parent inputs and calibration gate are frozen. Public points
and the known302 address do not choose prospective parameters or centres.

Replay from the research directory:

```sh
python3 elliptic-curves/cas/report_curve302_recovered_subgroup.py --check
```

Raw replay inputs and all stage ledgers remain under
`artifacts/local/elliptic-curves/curve302-recovered-subgroup-followup-v1/`.
Exact span producer/checker sources are
`certify_curve302_recovery_in_public_span.sage` and
`verify_curve302_recovered_public_span.sage`; the cloud union uses
`build_recorded_point_cloud_union.py` and the independent rank checkers.

# Quotient-aware rank-escape detector v2

Date: 2026-09-04  
Status: **Outcome D — fixture-separated MW17-only controls frozen but unrun; complete descents blocked; prospective sample unopened**

<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 f07ee569c95bf3a1 -->

## Result

Detector v2 now has a fail-closed exact measurement layer, a separately
blinded five-cohort sample, and fixture-separated MW17-only record-control
programs.  It is **not calibrated**: neither curve 356 nor curve 385 has a
completed MW17-relative global 2-descent, so the frozen prospective rows have
not been submitted to the pipeline.

The exact current control table is:

| curve | actual `im(MW17/2)` | exceptional image mod MW17 | actual displayed MW29 image | `dim Sel_2` | `s_res` | `dim Sel_2/im(MW29/2)` |
|---:|---:|---:|---:|---:|---:|---:|
| 356 | 17 | 12 | 29 | `UNKNOWN` | `UNKNOWN` | `UNKNOWN` |
| 385 | 17 | 12 | 29 | `UNKNOWN` | `UNKNOWN` | `UNKNOWN` |

Both exact root numbers are `-1`, so the 2-Selmer dimensions are odd under
the proved absence of rational 2-torsion.  Together with the certified
29-dimensional rational-point image this only restricts each dimension to
`29,31,33,...`; it does not choose 29.  If a future complete descent returns
29, the existing rank lower bound 29 immediately gives exact rank 29.  No
saturation hypothesis is needed for that lower-bound/upper-bound equality.

The specialized MW17 group is exactly the subgroup on `P1,...,P17`, and is
primitive inside the displayed subgroup on `P1,...,P29`.  The combined
good-reduction `E(F_p)/2E(F_p)` signature matrix has rank 29 on the displayed
points; hence MW29, and its MW17 subspace, are 2-saturated in `E(Q)`.  The
underlying generic function-field MW17 basis is saturated as well.  What
remains unproved is all-prime saturation of either specialized subgroup
inside the as-yet-unknown full `E(Q)`.

The compact certificate is
[`../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-controls-v1.json`](../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-controls-v1.json).

## Prospective blinding boundary

The earlier hard record-fibre route embeds and quotients all 29 public point
classes before its local calculation.  It is useful for post-discovery exact-rank
closure, but it removes precisely the twelve directions that a new candidate
would not possess.  It therefore cannot calibrate an operational candidate
gate.

The replacement replay is pinned in
[`../artifacts/generated-results/elliptic-curves/r17_mw17_only_selmer_control_inputs_v1.json`](../artifacts/generated-results/elliptic-curves/r17_mw17_only_selmer_control_inputs_v1.json).
Each generated Magma executable contains only the global minimal curve and the
seventeen specialized generic points.  A mechanical source audit requires
exactly 17 point declarations and forbids every `P18,...,P29` coordinate row,
held-out label, MW29 token, half-ideal shortcut, external file read, cover
search, or point search.  It computes the complete unconditional 2-Selmer group
first, verifies the MW17 Kummer rank, emits a basis of
`Sel_2(E)/im(MW17/2)`, and terminates at `blind_freeze`.

The separately pinned run ledger is
[`../artifacts/generated-results/elliptic-curves/r17_mw17_only_selmer_control_run_v1.json`](../artifacts/generated-results/elliptic-curves/r17_mw17_only_selmer_control_run_v1.json).
It currently records zero completed replays.  The Selmer machinery is not an
operational candidate gate until both source-hash-matched transcripts complete,
each certifies MW17 image dimension 17, and each reports residual dimension at
least 12 before the committed public control truth is consulted.  The control
truth is bound by coordinate hashes outside the executables.  Even a passing
two-record replay would calibrate the Selmer gate only; it would not by itself
authorize the prospective sample.

## Exact record models and places

The replay compares the public-fibre, exact-lineage, and local-Kummer
certificates coefficient for coefficient.  The global minimal models are

```text
356: [0,1,0,
 -24391876744717707263532695900840552395172973498186560300,
 46943906433780620456844832699051340439698711588743845207309557656274241785479710000]

385: [1,-1,1,
 -331827496674406562041164370816053963496434513510649284995530467,
 2322853282053688692296179831887155386042997843373920023045057766047634710076537198987220653859]
```

Every finite bad place has exact ambient local Kummer dimension, reduction
kind, Kodaira symbol, Tamagawa number, conductor exponent, minimal
discriminant valuation, exceptional-subgroup image dimension, and component
image-order multiset in the certificate.  The complete bad-prime lists are:

```text
356: 2,3,5,13,23,29,37,41,139,751,28960331,
     1204882855601765528877267647500895974865482613,
     197980272243427555346397293722916980361535459279712115031762027678304939

385: 2,5,7,11,13,29,37,41,43,47,73,89,109,127,
     4678955899327531799218956351083,
     3466418496046307687463380088830029870757187849126892374834076735615029496431422958053227129921575121
```

Infinity is explicit.  Curve 356 has one real component and local Kummer
dimension zero.  Curve 385 has two components and dimension one; the MW17
image already has dimension one there, while all twelve exceptional points
add zero modulo MW17.  Exact rational root-isolation intervals and the 29
component bits are stored.

These are exact local curve/known-point data, not the global Selmer condition
matrix.  Until the global `S`-squareclass domain exists, the local condition
row spaces, local codimensions, auxiliary descent places, and
`s_res^(-v)` remain null.

## Measurement implementation

[`../elliptic-curves/cas/quotient_rank_escape_detector_v2.py`](../elliptic-curves/cas/quotient_rank_escape_detector_v2.py)
implements the post-descent measurement over `F_2`.  Given a completed global
norm-squareclass basis, every all-place condition row, and the MW17 Kummer
rows, it computes:

- the canonical left-pivot RREF of the stacked global local-condition matrix;
- actual MW17 image dimension and `s_res`;
- a deterministic presentation of `Sel_2/im(MW17/2)`;
- summed and independent local codimensions;
- `s_res^(-v)` and single-place suppression for every place;
- the held-out exceptional image and remaining quotient after MW29 is loaded;
- a predeclared rank-profile distance for cross-curve comparison.

The quotient presentation is coordinate-deterministic only.  The code sets
the pairing claim to null and does not infer a pairing on a complement.
Complete row spaces for different cubic fields do not share an ambient vector
space, so a direct Grassmann distance is undefined.  The declared comparison
distance uses only the basis-independent sorted per-place triple
`(local codimension, leave-one-out suppression, local Kummer dimension)`.

The checkpointed Simon worker now records infinity even when its condition is
vacuous, loops over every elliptic bad prime in addition to descent support,
stores all local Kummer dimensions, and exports the Selmer basis in global
norm coordinates.  A class-number-one smoke curve completes through this
interface and includes its otherwise redundant bad prime 37.

## Frozen stratified sample

The original six lanes induce five comparison cohorts by pooling the two full
cylinders while preserving their 356/385 anchor strata:

```text
full cylinders, matched ordinary, 2-only, odd-only, random equal-codimension.
```

Within every cohort/anchor stratum, selection is solely lexicographic order
of the already frozen 24-hex `sample_id`.  Stage 1 takes hash rank one in each
anchor stratum: two fibres per cohort, ten total.  Stage 2 takes ranks one
through three: six per cohort, thirty total.  Stage 2 contains Stage 1.

The descent input has no cohort or anchor labels:
[`../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-sample-v1.json`](../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-sample-v1.json).
The separately frozen key is
[`../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-unblinding-key-v1.json`](../artifacts/generated-results/elkies-k3-r17-quotient-rank-escape-detector-v2-unblinding-key-v1.json).
Neither stage is authorized until both record controls complete and blindly
recover all twelve exceptional quotient directions.

## Exact obstruction and narrowed task

The completed-square cubic maximal orders and every bad-prime factorization
are certified and cached.  The obstruction is the global class/unit step.
Known Kummer classes already force full auxiliary cubic class-group 2-rank
lower bounds 21 for curve 356 and 15 for curve 385; after the recorded unit
allowance, the known exceptional blocks force residual class-image lower
bounds 11 and 10.  A generic full-BNF computation therefore spends most of
its effort reconstructing classes that the calibrated quotient will remove.

Preserved threaded PARI runs include strict 1,200-second timeouts for both
fields.  A broader curve-385 run reached six factor-base stages and ended
after 2,693.58 seconds with a requested relation deficit of 5,224 relations
for 5,219 ideals.  No run produced a certified BNF checkpoint.  Independent
five-minute eclib Selmer-only canaries on both curves also returned no Selmer
dimension.  These failures are operational measurements only.

The prospective technical task is therefore:

1. run the sealed curve-356 and curve-385 executables with no access to the
   held-out coordinates or half-ideals;
2. compute and certify the complete unconditional 2-Selmer group, including
   every finite and real local condition;
3. quotient only by the certified specialized MW17 image and freeze the
   residual dimension and basis;
4. require residual dimension at least twelve on both records before consulting
   the committed public control truth;
5. retain quotient-by-MW29 collectors only for post-discovery closure work.

This is Outcome D.  There is no cohort comparison, cylinder conclusion,
rank-enrichment claim, or candidate promotion.

## Replay

```bash
python3 elkies-k3/scripts/build_r17_quotient_rank_escape_detector_v2_sample.py --check
python3 elliptic-curves/cas/build_r17_mw17_only_selmer_replay.py --check
python3 elliptic-curves/cas/run_r17_mw17_only_selmer_replay.py --check
sage -python elkies-k3/scripts/certify_r17_quotient_rank_escape_detector_v2_controls.sage --check
python3 -m unittest \
  elliptic-curves/tests/test_quotient_rank_escape_detector_v2.py \
  elliptic-curves/tests/test_elkies_relative_2selmer_checkpointed.py \
  elliptic-curves/tests/test_r17_mw17_only_selmer_replay.py
```

To materialize the two ignored local executables and, on a licensed host, run
them under the declared one-day/16-GB per-case envelope:

```bash
python3 elliptic-curves/cas/build_r17_mw17_only_selmer_replay.py --overwrite
python3 elliptic-curves/cas/run_r17_mw17_only_selmer_replay.py \
  --execute --overwrite
```

The older fixture-sequenced Magma inputs for the record pair are retained in
[`../artifacts/generated-results/elliptic-curves/elkies_2026_record_pair_relative_2selmer_inputs_v1.json`](../artifacts/generated-results/elliptic-curves/elkies_2026_record_pair_relative_2selmer_inputs_v1.json).
They remain reproducibility evidence, not the prospective calibration input.

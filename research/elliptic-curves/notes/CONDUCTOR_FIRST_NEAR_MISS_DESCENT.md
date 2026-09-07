# Conductor-first descent on the four low-conductor near misses

## Status

This is an active exact-descent programme, not a new rank theorem.  The four
fixed targets now have source-pinned models and exact known subgroups whose
images in `E(Q)/2E(Q)` have the full certified dimensions `20,20,19,19`.
Complete 2-Selmer computations are still open.

| target | certified rank lower bound | `log(N)` | directions needed for rank 21 | known mod-2 image |
| --- | ---: | ---: | ---: | ---: |
| ICARM 245 | 20 | `150.668907...` | 1 | 20 |
| Fermigier--Mestre `u=28917/20` | 20 | `159.934825...` | 1 | 20 |
| split-infinity family 2, `u=483` | 19 | `157.759935...` | 2 | 19 |
| split-infinity family 3, `u=660` | 19 | `164.053646...` | 2 | 19 |

The pinned input is
[`conductor_first_near_miss_descent_targets_v1.json`](../../artifacts/generated-results/elliptic-curves/conductor_first_near_miss_descent_targets_v1.json).
It records every point, the finite-quotient matrices, the trivial rational
2-torsion witnesses, the exact 2-division cubics, source hashes, conductor
data, and the required residual quotient dimensions.

## Removing a false mod-2 bottleneck on the Mestre fibres

The rank-19 bases originally selected by the mod-3 certificates have
one-dimensional images in the bounded product of mod-2 finite quotients.  A
relative 2-descent against those lists would therefore spend most of its
output rediscovering subgroup-index information.

A deterministic PARI `ellsaturation(E,P,3)` pass returns much smaller rational
points on each short model.  The target builder does not treat this call as a
global saturation proof.  Instead it independently checks every returned
point and proves that each 19-point replacement list has binary rank 19 in
exact products of `E(F_p)/2E(F_p)`.  Thus each replacement list is an exact
rank-19 subgroup with full 19-dimensional Kummer image, regardless of what
PARI proves about saturation.

The retained certificate primes are:

```text
family2_u483:
13,17,31,47,59,67,71,73,79,83,97,101,127,151,179

family3_u660:
17,23,31,37,41,47,53,59,67,71,73,79,83,97,101,109,131,151
```

The Fermigier target similarly uses the already pinned small-prime-saturation
candidate basis, not the older imported list.  Its exact mod-2 certificate has
rank 20.  ICARM 245 already had a full rank-20 mod-2 certificate.

## Exact Selmer jobs

Run the dimension-only 2-Selmer pass first:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --target icarm-245 \
  --mode selmer-dimension \
  --output artifacts/local/elliptic-curves/icarm245-2selmer-dimension.m
```

The other target identifiers are:

```text
fermigier-u28917-20
family2-u483
family3-u660
```

This calls `TwoSelmerGroup(E)` without constructing quartics.  Because the
pinned known subgroup has full mod-2 image and rational 2-torsion is trivial,
the exact residual dimension is

```text
dim Sel_2(E) - certified known Kummer dimension.
```

A zero value closes the fibre before any cover construction.  For a positive
value, generate the relative-cover pass with `--mode relative-covers`.  That
program calls

```magma
TwoDescent(E :
    RemoveGens := SequenceToSet(known),
    RemoveTorsion := true
);
```

and uses `CFNMSEL|` protocol records.  It prints every nonzero residual
2-cover and, below a declared cap, the complete Cassels--Tate pairing.  No GRH
class-group bounds are enabled.

An independent unconditional 3-Selmer upper-bound job is also available:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --target icarm-245 \
  --mode three-selmer-dimension \
  --output artifacts/local/elliptic-curves/icarm245-3selmer-dimension.m
```

It calls `ThreeSelmerGroup(E)` but does not construct ternary cubics.  If its
rank upper bound equals the certified lower bound, the fibre is independently
closed even when the 2-Selmer computation is obstructed.  A positive gap does
not identify a rational direction and does not replace the relative 2-cover
job.

A bounded local PARI diagnostic is available separately:

```sh
.venv/bin/python elliptic-curves/cas/run_conductor_first_pari_diagnostic.py \
  --target icarm-245 \
  --output artifacts/local/elliptic-curves/icarm245-pari-diagnostic.json
```

It verifies that every pinned point lies on the exact `descent_model` before
starting GP, supervises wall time and resident memory, and records raw
`CFNMPARI` protocol output.  Its upper endpoint is not an unconditional
Selmer/rank bound: PARI's cubic-field `Buchall` path uses provisional BNF data
and does not call `bnfcertify`.  Any returned rational point is exact-checked
and admitted as a new unconditional lower-bound direction only after a fresh
full mod-2 finite-reduction certificate.

## BNF-free S-class collector priority

The exact planning artifact
[`conductor_first_s_class_envelopes_v1.json`](../../artifacts/generated-results/elliptic-curves/conductor_first_s_class_envelopes_v1.json)
constructs each integral monic 2-division cubic from the global minimal model,
computes its field discriminant, and materializes every prime ideal over the
Bach bound together with every prime ideal over `2*Delta(E)`.  The resulting
collector order is:

| target | signature | Bach/ERH bound | prime ideals including `S` |
| --- | ---: | ---: | ---: |
| ICARM 245 | `(1,1)` | `245977` | `39904` |
| family 2, `u=483` | `(3,0)` | `252689` | `40754` |
| Fermigier, `u=28917/20` | `(3,0)` | `262523` | `42251` |
| family 3, `u=660` | `(3,0)` | `271542` | `43512` |

Thus ICARM 245, rather than Fermigier, is the first target for the generalized
Minkowski/special-q relation collector.  These counts and fields are exact;
Bach's prime-generation conclusion is conditional on ERH/GRH.  No principal
relation completeness, S-class quotient, Selmer dimension, or rank bound is
claimed by this comparison.

Interpretation is exact:

- `CLOSED_EXACT_RANK_20` or `CLOSED_EXACT_RANK_19` means the residual
  2-Selmer quotient is zero and the fixed fibre is unconditionally closed.
- The corresponding `_CT` outcome means residual Selmer classes exist but
  the pairing radical is zero; the fixed fibre is again closed.
- `HIGHER_DESCENT_REQUIRED` identifies the exact radical classes on which to
  run 4-descent.  It is not evidence for another rational point.
- `RESIDUAL_COVERS_UNPAIRED` means only that the configured pairing cap was
  reached.

If a radical class remains, use `FourDescent` only on that class, remove the
same known elliptic-curve generators again, minimize and reduce the resulting
intersection-of-quadrics models, and pair them against the full residual
2-cover set before any point search.  A mapped point must still pass the
repository's exact finite-quotient escape test before promotion.

## Exact-engine calibration and current resource boundary

The following calibration was run on 2026-08-31.  It records resource failures,
not Selmer evidence:

- PARI/GP 2.15.4: strict 30-second `ell2cover` probes on ICARM 245 and
  both Mestre fibres did not finish.  A two-GB, 300-second family-2 probe also
  timed out.  Both `bnfinit(f,1)` and `bnfinit(f,0)` on the reduced family-2
  cubic field timed out at the same two-GB/300-second boundary.  A subsequent
  unconditional `bnfinit(f,1)` run with a 3.5-GB PARI stack, five-GB address
  limit, and 600-second wall limit also timed out without returning a BNF
  object; it remained CPU-bound and used roughly one GB resident memory.  The
  faster `bnfinit(f,0)` route likewise returned no provisional BNF within 900
  seconds, so the cheaper `bnfcertify(b,1)` quotient certification could not
  begin.
- eclib/mwrank 20231212: all four targets reach the same machine-integer
  conversion failure in Selmer-only 2-descent (`lower bound on c too large`)
  and return no rank bound.  In repository-local replays with 22 auxiliary
  primes, ICARM 245 fails when converting a quartic coefficient lower bound
  of about `-1.72e28`, while family 3 at `u=660` fails at about `-2.24e31`;
  the earlier Fermigier and family-2 runs have the identical failure class.
  A source build at tag `v20231212`, commit
  `c4a9d5be304c1cc80cb5020ebeeecd3f237ca4c6`, was then patched to use NTL
  integers for the quartic-enumeration variable `c` and its bounds; the pinned
  patch has SHA-256
  `eacaf6a9667d5f0ea70652913cdb7dce681cb66f8f349fdbdeee78a686cdac3d`.
  The patched executable agrees with the packaged executable on the rank-5
  smoke curve `[0,0,1,-79,342]` and passes the previous ICARM 245 conversion
  failure at 256-bit real precision.  It then remains in the first Type-3
  quartic enumeration until the strict 30-second wall stop (peak RSS 19,520
  KiB), returning no rank or Selmer bound.  An independent 100-decimal PARI
  evaluation gives `5,728,687,860,386,985,887,994,542,644` actual `c` loop
  values in the first `a=1,b=0` slice before the auxiliary local flags.  Thus
  widening the integer type removes this representation failure but does not
  make the existing exhaustive enumeration viable.  Current eclib master at
  commit `8e0f64171a663f11e859afa15b29e589328e1738` was also attempted, but its
  FLINT interface does not compile against Ubuntu's FLINT 3.0.1 headers; no
  arithmetic conclusion is drawn from that build failure.
- Official Magma calculator V2.29-9: `TwoSelmerGroup(E)` reaches the service
  memory limit on all four targets before assigning a Selmer group.  The
  dimension-only call has the same failure as `TwoDescent`, locating the
  bottleneck before quartic construction.  `ThreeSelmerGroup(E)` passes the
  initial torsion stage on all four targets but reaches the service's
  60-second limit without returning a group.
- PARI `ellisomat(E,0,1)` finds one rational isogeny class representative for
  each target, so none of the four has a nontrivial rational isogeny route to a
  cheaper 2- or 3-isogeny descent.
- A temporary build of PARI 2.19.0-development at upstream commit `6af5b91`
  did not return a provisional BNF for the family-2 `u=483` field within 900
  seconds.  For the smaller `u=481` field, `bnfinit(f,0)` completed after
  `209496` ms with provisional class group order `32` and elementary divisors
  `[2,2,2,2,2]`.  This remains GRH-conditional: `bnfcertify(b,1)` warned that
  the Zimmert bound was `886251290812616611594240` and did not finish before
  the same 900-second total cutoff.  No unconditional class-group or Selmer
  claim is extracted from it.
- One-hour fail-closed `ellrank(...,0,known)` diagnostics with that PARI 2.19
  build returned no interval or point on Fermigier, family 2 at `u=483`, or
  family 3 at `u=660`.  The durable family-2 and family-3 supervisors recorded
  peak observed RSS of `150757376` and `395194368` bytes respectively before
  their wall stops.  ICARM 245 first overflowed a 12-GB PARI stack; its 16-GB
  retry was stopped at the explicit 18-GB RSS boundary after `2774.79`
  seconds, with peak observed RSS `18001264640` bytes and still no interval.
  These are resource terminals, not Selmer/rank bounds or bounded point-search
  evidence.
- The generalized BNF-free Minkowski/special-q collector was calibrated on
  ICARM 245 using the complete 39,904-column Bach-plus-`S` factor base.  After
  appending and independently verifying all 21,713 canonical principal rows
  `(p)`, the ledger auditor certifies under ERH/GRH that
  `dim Cl(O_K[S^-1])/2 <= 18181`, hence the unfiltered envelope has
  `dim K(S,2) <= 18202`.  This is not a Selmer bound: norm and local Kummer
  conditions remain unapplied.  At lattice-combination bound five, 46,450
  single-ideal candidates outside the factor base produced 56 exact closed
  cycles and raw relation rank three, but zero rank after quotienting by the
  canonical rows and `S`.  The matched in-factor-base run produced no closed
  relation, and the 46,470-candidate paired-special cycle run likewise
  produced none.  A matched degree-two-special-ideal run produced 65 exact
  cycles and raw rank three, again with zero post-quotient gain.  These
  bounded zero-gain results close the collector's current single degree-one,
  paired degree-one, in-base, and degree-two geometries; they do not establish
  relation-lattice completeness.
- The matched family-2 `u=483` calibration sampled 59,600 single-ideal
  candidates outside its complete 40,754-column Bach-plus-`S` factor base.
  It produced 153 exact closed cycles and raw relation rank three, again with
  zero rank beyond the 22,257 canonical `(p)` rows and `S`.  The independent
  ledger audit consequently certifies under ERH/GRH only
  `dim Cl(O_K[S^-1])/2 <= 18483` and `dim K(S,2) <= 18511`.  This second
  zero-gain field confirms that the initial special-q lattice geometry, not
  the ordering of the targets, is the part that must change.
- The quadratic `a+b*theta+c*theta^2` collector now accepts an arbitrary
  monic cubic and declared Selmer primes, and its optional sparse-hypergraph
  eliminator retains exact generator witnesses while cancelling any bounded
  number of large-prime ideal columns.  On ICARM 245, 20,000 hybrid candidates
  trial-divided through 10,000 and a matched 2,000-candidate full-factor-base
  trial produced no partial relation.  Exact factorization of 200 samples
  found respectively `30,80,63,21,5,1` norms with one through six odd prime
  factors outside the Bach base.  Retaining the 30 one-residual and 80
  two-residual cases as exact two- and three-vertex hyperedges produced no
  dependency.  A subsequent adaptive walk reused degree-one residual ideals
  as new special ideals, excluded duplicate principal generators, and ran two
  capped generations.  Its 1,345 candidates supplied 796 accepted exact
  hyperedges (`110,450,236` by depths zero, one, and two), with zero
  dependency and hence exact incidence rank 796.  Thus naive residual reuse
  still gives a forest-like graph in this bounded sample.
- A follow-up paired-ideal experiment deliberately arranged the ten smallest
  reused residual ideals in a cycle and sampled 324,000 elements at
  lattice-combination bound ten.  It returned zero closed cycles and zero
  post-quotient gain, but the old hybrid cofactor gate retained only one
  partial edge.  Therefore this run does **not** show that a multiple-residual
  cycle geometry fails; it only shows that pairing the ideals without exact
  residual factorization does not populate the intended graph.  This corrects
  the earlier, too-strong prescription that merely forcing paired residual
  ideals would be the next productive geometry.
- The full-ideal Minkowski collector now reuses the same exact
  sparse-hypergraph eliminator as the quadratic collector.  Its opt-in exact
  mode fully factors the remaining norm, retains a bounded number of odd
  residual prime-ideal columns, and records the exact incidence rank and
  nullity.  A pre-optimization full-Bach pilot spent its 600-second wall limit
  constructing the 39,904-column field/factor-base state and emitted no
  collection record.  A smaller-factor-base setup was then interrupted on the
  user's stop request while PARI was reconstructing the maximal order.  The
  collector now supplies the already declared discriminant primes to PARI
  before that construction, but no post-change exact-factor pilot has been
  run.  Consequently there is no new principal relation, S-class bound,
  Selmer bound, or rank conclusion from this unfinished path.

No output from these bounded probes changes a mathematical status.  Complete
jobs require a full Magma installation with substantially more memory/time,
or an independent class-group computation that can be imported and certified.
The durable output is therefore a three-stage, source-pinned Magma workflow,
not another increase of a raw `ell2cover` or `ratpoints` bound.

## Replay

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_targets.py --check

.venv/bin/python -m unittest \
  elliptic-curves/tests/test_conductor_first_near_miss_descent.py
```

The first command reruns the quick PARI basis replacement on the two Mestre
fibres and replays all four exact Kummer-image certificates.  It does not run
a complete Selmer computation.

## Family-search order

The surrounding Mestre/Fermigier search should now use this strict gate order:

1. discriminant radical and squareful-value sieve;
2. exact global minimalization and local Tate data;
3. full-dimensional known-subgroup Kummer certificate;
4. residual Selmer dimension relative to that exact subgroup;
5. explicit residual-cover recovery;
6. point search only on surviving covers.

Candidate comparison is lexicographic and Pareto-based: maximize the certified
rank lower bound first and minimize exact conductor second.  Local scores,
Nagao scores, point counts, and discriminant radicals may order work within a
gate, but they are not blended into the mathematical frontier.

[`conductor_first_pipeline.py`](../cas/conductor_first_pipeline.py) enforces
this order at the record boundary.  In particular, it rejects a point-recovery
record without completed exact Tate data, a complete residual Selmer result,
and an explicitly constructed locally surviving cover.  Stage a family batch
with:

```sh
.venv/bin/python elliptic-curves/cas/stage_conductor_first_family_candidates.py \
  --input artifacts/local/elliptic-curves/conductor-first-candidates.json \
  --output artifacts/local/elliptic-curves/conductor-first-queues.json
```

The output separates the Tate, residual-Selmer, residual-cover, and
point-recovery queues and includes both the rank-first ordering and the exact
rank/conductor Pareto frontier.

## Closed anchor-neighborhood pilot

The first applied family ledger is
[`conductor_first_family_anchor_pilot_v1.json`](../../artifacts/generated-results/elliptic-curves/conductor_first_family_anchor_pilot_v1.json).
Its closed population consists of offsets `-4..4` in each of three lanes:
Fermigier `u=(28917+k)/20`, family 2 `u=483+k`, and family 3 `u=660+k`.
All 27 fibres receive the cheap discriminant-factorization sieve.  Per lane,
the anchor and the two best other sieve rows receive exact global
minimalization, local Tate data, and a full-dimensional known mod-2 Kummer
certificate.  The ledger builder runs no standalone point search.  Its
`u=481` record additionally pins two exact points returned incidentally by a
separate provisional descent diagnostic, while withholding that diagnostic's
upper bound.

The nine exact residual-Selmer inputs are:

| fibre | certified known rank | Kummer dimension | `log(N)` | root number |
| --- | ---: | ---: | ---: | ---: |
| family 2, `u=481` | 14 | 14 | `128.494548...` | `+1` |
| family 2, `u=483` | 19 | 19 | `157.759935...` | `-1` |
| family 2, `u=484` | 12 | 12 | `188.515120...` | `+1` |
| family 3, `u=660` | 19 | 19 | `164.053646...` | `-1` |
| family 3, `u=662` | 12 | 12 | `210.776762...` | `+1` |
| family 3, `u=663` | 12 | 12 | `204.814202...` | `+1` |
| Fermigier `u=28913/20` | 12 | 12 | `174.193694...` | `+1` |
| Fermigier `u=28915/20` | 12 | 12 | `169.827590...` | `+1` |
| Fermigier `u=28917/20` | 20 | 20 | `159.934825...` | `+1` |

The new `u=481` conductor is below ICARM 245, but its current certified rank
is only 14.  It is a low-conductor Selmer-triage survivor, not a rank-21 near
miss.  Both its unconditional 2- and 3-Selmer dimension jobs reached the
public calculator's 60-second limit without returning a group.  PARI 2.19's
provisional `ell2cover` basis has dimension 14, and helped `ellrank` reported
`[14,14]` with Sha indicator zero, but those upper-bound data depend on the
uncertified BNF.  The two new rational points recovered during that run are
independently exact and raise the unconditional lower bound from 12 to 14;
their full 14-dimensional mod-2 certificate uses primes
`23,31,37,53,83,89,101,103,109`.

Replay the ledger and generate its next job with:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_family_pilot.py --check

.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --manifest artifacts/generated-results/elliptic-curves/conductor_first_family_anchor_pilot_v1.json \
  --target family2-u481 \
  --mode selmer-dimension \
  --output artifacts/local/elliptic-curves/family2-u481-2selmer-dimension.m
```

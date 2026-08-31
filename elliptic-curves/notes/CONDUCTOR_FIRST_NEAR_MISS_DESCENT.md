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

## Relative 2-Selmer job

Generate one unconditional Magma job with:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --target icarm-245 \
  --output artifacts/local/elliptic-curves/icarm245-relative-2selmer.m
```

The other target identifiers are:

```text
fermigier-u28917-20
family2-u483
family3-u660
```

The generated program calls

```magma
TwoDescent(E :
    RemoveGens := SequenceToSet(known),
    RemoveTorsion := true
);
```

and uses `CFNMSEL|` protocol records.  It prints every nonzero residual
2-cover and, below a declared cap, the complete Cassels--Tate pairing.  No GRH
class-group bounds are enabled.

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

## Local PARI calibration

PARI/GP 2.15.4 is available in the current repository environment; Magma,
Sage, and mwrank are not.  Strict 30-second probes of `ell2cover` on ICARM 245
and both Mestre fibres did not finish.  A single two-GB, 300-second follow-up
on `family2_u483` also remained inside `ell2cover` until the fixed timeout.
These timeouts supply no Selmer dimension, cover, or rank evidence.  The
earlier Fermigier work identifies the same cubic class-group bottleneck.  This
is why the durable output is a relative, source-pinned job rather than another
widening of a raw `ell2cover` or `ratpoints` bound.

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
3. known-subgroup Kummer certificate and residual Selmer dimension;
4. explicit residual-cover recovery;
5. point search only on surviving covers.

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

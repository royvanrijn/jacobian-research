# X1092: nineteen rootless MW17 frame types

The complete pruned Niemeier census has **19 integral-isometry classes** of
rootless rank-17 complements of determinant 1092. The recovered curve302 frame
is class **6**. Thus there are **18 additional J2/frame types** on this lattice.
This is not a J1/surface-automorphism classification and does not supply
rational marked nef-U realizations, elliptic equations, or specialized ranks.

The scientific change is that the parent-construction queue is now finite and
explicit. X1092 has geometric Picard rank 19, so Shioda--Tate caps its generic
MW rank at 17. We now know all frame types attaining that ceiling. This gives
eighteen alternatives to curve302's frame on the same surface, each with an
exact lattice witness. It does not establish a new surface or any better
specialization tail. No new rational curve of rank 28 or higher was found by
this computation.

The contrast with X948's two rootless J2 types is a structural observation
worth retaining. It is not evidence that the number of frame types predicts
specialization rank, or an explanation of curve302's rank-at-least-31 fibre.

## Exact coverage and algorithm

The certified rank-seven Nishiyama auxiliary contains the fixed D5. The
complete D5-anchor catalogue covers sixteen anchor orbits in thirteen rooted
Niemeier types; other rooted types contain no D5, and the Leech lattice has
no roots. Every primitive auxiliary embedding therefore appears in this cover.
The stabilizer's residual Weyl group puts the sixth vector in a dominant
chamber. The seventh vector must be strictly dominant on every remaining root
component precisely when its final complement is rootless.

For current base B with prescribed pairing row b, the projection has norm
`b (B G B^t)^(-1) b^t`. Residual ADE components are orthogonal both to B and
to each other. Their labels l contribute `l C^(-1) l^t` to a **shared** norm
budget. An exact branch is discarded as soon as this sum, plus the smallest
remaining component costs, exceeds the budget. Strict positive labels have
minimum at the all-one vector because ADE inverse-Cartan entries are positive.
The implementation checks the resulting block-diagonal Schur-complement
identity exactly in every chamber reaching the enumeration stage.

Remaining affine shells use rational LDL enumeration, exact integrality and
Smith-form primitivity checks. Orthogonal complements use saturated integer
kernels. There is no score cutoff, numerical epsilon, or bounded-search
inference. Counts below are representatives in the declared Weyl cover,
**not** all embeddings or full ambient-automorphism embedding-orbit counts.

| Anchor | Primitive sixth vectors | Rootless embeddings |
|---|---:|---:|
| D24:1 | 4 | 0 |
| D16_E8:1 | 20 | 0 |
| D16_E8:2 | 56 | 0 |
| 3E8:1 | 75 | 0 |
| 2D12:1 | 59 | 0 |
| A17_E7:1 | 350 | 0 |
| D10_2E7:1 | 411 | 0 |
| D10_2E7:2 | 162 | 0 |
| A15_D9:1 | 230 | 0 |
| 3D8:1 | 437 | 0 |
| A11_D7_E6:1 | 874 | 0 |
| A11_D7_E6:2 | 854 | 0 |
| 4E6:1 | 969 | 0 |
| 2A9_D6:1 | 1794 | 0 |
| 4D6:1 | 1878 | 0 |
| 2A7_2D5:1 | 1925 | 208 |

The total is **10,098 primitive sixth-vector representatives**. The first
fifteen anchors yield no rootless complement; the last supplies all 208.
Their order is an implementation choice. Partial zero counts carried no
evidence for the eventual classification.

Every class has minimum 4. Minimum-vector counts include both signs.

| Class | Minimum vectors | Automorphism order | Embeddings in cover |
|---|---:|---:|---:|
| 1 | 2436 | 2 | 16 |
| 2 | 2446 | 4 | 8 |
| 3 | 2438 | 2 | 16 |
| 4 | 2438 | 2 | 16 |
| 5 | 2442 | 2 | 16 |
| 6 (curve302) | 2436 | 2 | 16 |
| 7 | 2436 | 2 | 16 |
| 8 | 2442 | 2 | 16 |
| 9 | 2434 | 2 | 16 |
| 10 | 2440 | 4 | 8 |
| 11 | 2440 | 4 | 8 |
| 12 | 2438 | 8 | 4 |
| 13 | 2440 | 4 | 8 |
| 14 | 2452 | 8 | 4 |
| 15 | 2438 | 4 | 8 |
| 16 | 2446 | 4 | 8 |
| 17 | 2448 | 4 | 8 |
| 18 | 2444 | 4 | 8 |
| 19 | 2438 | 4 | 8 |

The certificate retains all nineteen exact representative Grams, every
embedding-to-class integral isometry, and the explicit known-frame isometry.
Class comparison uses exact PARI integral isometry, not minimum-vector counts
or other fingerprints. In particular, the known frame is tested only after
all complements have been retained and deduplicated.

## Validation and performance

Two completed small anchors agree vector-for-vector with the old algorithm,
including every subsequent seventh-vector calculation. Four small positive
shell controls check strict/non-strict and root-free residual cases. Those
controls exposed an old LDL indexing error, repaired in the successor source.
A full separate run reproduces **all sixteen anchor packets byte-for-byte**;
the separate witness-assembly pass also reproduces the final certificate.
This is a clean replay of the same enumerator, not an independently written
proof of enumeration completeness.

Enumeration used about **104 worker-seconds**, split over three detached
workers. Anchors 11 and 12 took 8.3 and 8.8 seconds; anchor 16 took 15.2 seconds.
Isometry classification and replay are additional costs. The original process
was left untouched during replacement and validation, then stopped after
explicit user authorization. Its log and termination receipt are preserved in
`artifacts/local/elkies-k3/det1092-rootless-j2-census-v1b*`. Its hours of
computation do not contribute to acceptance of this certificate.

The slow implementation first formed Cartesian products of independent
component-label lists, granting each component the whole norm budget. For
anchor 12 this produced 2,580,928 sixth-stage label combinations before the
shared norm test. Expensive rational matrix calculations then rejected most
combinations. It also rebuilt root data and LDL factorizations and logged
only completed anchors. The replacement applies the exact shared budget
before that work, caches repeated data, and checkpoints within anchors.

Three correctness defects were repaired in the successor source:

1. Every rootless complement was required to match the known curve302 frame.
   A new type would therefore abort the run. All types are now retained and
   deduplicated before the expected known-class check.
2. The minimum was initially hard-coded as four; an attempted replacement
   also reversed PARI's count/minimum return fields. Exact minimum and signed
   minimum-vector count are now computed and checked on a small control.
3. The LDL off-diagonal update used an out-of-scope index in its denominator.
   A positive residual-space test exposed this failure; the denominator now
   uses the current diagonal entry.

No end-to-end speedup ratio is claimed: the old implementation never finished,
and the 104 worker-seconds measure enumeration, excluding startup, final class
assembly, validation and the clean replay. Future performance reports should
include those costs and failed attempts explicitly.

## Prospective realization order

Each additional class has a representative sharing its first six auxiliary
vectors with a known-class embedding in the same Niemeier model. Their frames
therefore share an exact primitive rank-16 core. The priority packet orders
these candidates by common-core determinant and auxiliary-coordinate size,
with class index breaking ties. It selects **class 1**, with core determinant
4100. Thirteen types tie on the first two criteria; this is a tractability
heuristic, not a theorem about equation size or specialized rank.

The next gate is an explicit transport to the known rational NS marking,
followed by a primitive nef U and effective rational zero. Only an exact
rational elliptic equation and generic MW17 certificate may release a
specialization panel. Carrier and constructed-strict-class machinery follows
those arithmetic gates. Broad A1/MW16 expansion remains stopped.

## Recommended continuation

**Build one usable new MW17 parent, starting with class 1.** Do not enlarge
the census or collect increasingly expensive equivalent embeddings. The
current priority is a deterministic starting order, not a demonstrated best
parent: classes 1, 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 17 and 19 tie on the
two cost criteria. A concrete obstruction or a materially cheaper realization
can move another tied class ahead.

### 1. Turn the common core into a source-marked construction problem

Use known embedding **25** and new embedding **33** in the anchor-16 packet,
sharing sixth index **155**. These are zero-based packet indices. Their
primitive common rank-16 core has determinant **4100**, and its basis and Gram
are in the priority certificate. Bind both ambient bases and the class-6
isometry to the recovered generic MW17 basis before any transport.

Compute the core inclusions, orthogonal rank-one factors and finite-index
gluing on both sides. Then attempt a compatible transport between the
stabilized lattices `U + W_known(-1)` and `U + W_1(-1)`, or use that data to
give an explicit request to the existing
[marked-U planner](../../elkies-k3/scripts/plan_marked_u_realizations.sage).
The shared Niemeier core is input to this step, not already a marked NS
isometry. Existence of an abstract stable isometry does not prove that a
chosen core identification extends as requested.

**Deliverable:** an exact integral transport placing the target primitive U
and its frame in the known rational NS marking. Verify the full Gram identity,
unimodularity and an exact isometry of the target complement to class 1.
A failed bounded construction leaves that realization UNKNOWN and retains the
failed route; it does not exclude class 1 or justify an unbounded deeper run.

### 2. Prove the physical and arithmetic realization gates

Choose a positive effective primitive fibre class, establish nefness, and
produce an effective rational zero with fibre intersection one. Preserve any
Weyl reflections and zero changes in the exact source marking. A lattice
vector satisfying square zero alone does not pass this stage.

**Deliverable:** a rational marked nef-U certificate tied to the known K3 and
the selected class. Keep its J1 orbit UNKNOWN unless separately established.

### 3. Compile and certify one equation

Construct the Riemann--Roch pencil, derive its rational elliptic equation and
exact map to the source model, then compactify while retaining the transport.
Recover seventeen rational sections and prove their independence over Q(t)
with an exact height Gram. Check the fibre configuration/rootlessness and
identify its frame with the selected class. The Picard-19 ceiling then closes
the generic rank at exactly 17.

**Deliverable:** one self-contained parent packet containing equation, maps,
sections, height pairing and a replay. If another presentation is cheaper,
retain it after exact deduplication; the class index alone is not a coordinate
quality measure.

### 4. Run carrier and strict-class construction before parameter search

On the new marked equation, compute the accessible carrier characters and
their certified dependencies modulo the marked generic subgroup. Do not copy
curve302's strict dimension or an old chart's carrier quotient into the new
fibration. Record complete versus bounded carrier coverage and all cold
construction costs. Use independently supplied prospective selectors from
Agent 1 when their inputs and exact output semantics are pinned.

The [current strict-seed dispatcher](../rank-jump/prospective_constructed_strict_seed_gate.py)
is an interface prototype: its tests check state transitions, but it does
**not** replay certificate contents or verify hashes of named certificates.
It is not yet a production admission checker. Before connection, add
hash-bound certificates, producer/verifier receipts, matching curve/parameter/
model/subgroup identities, and provenance checks. A principal dependency
should enter a construction work queue without being labelled a strict class.

Keep these four endpoints separate:

| Gate | Required evidence | Permitted next action |
|---|---|---|
| Class constructed | Explicit nonzero strict representative and exact local/global class checks | Build and analyze its cover |
| Cover soluble | A rational witness or an unconditional rational-solubility proof | Extract an explicit point through the certified map |
| New rational point | Exact point and map checks, with novelty relative to the declared input point set | Test independence modulo the declared subgroup |
| Independent quotient direction | Exact independence certificate tied to that subgroup and model | Submit the enlarged subgroup to adaptive V3 |

On the **specified curve302 fibre and marked generic subgroup**, generic
strict dimension is zero. Its first certified nonzero constructed strict
class is therefore already outside that generic strict subgroup, and should
be preserved as a meaningful construction event. It still proves neither a
rational point nor a Mordell--Weil rank gain. A soluble point can be new as a
coordinate while remaining in the existing rational span. Local solubility,
zero CT pairing, split counts and principal-relation counts cannot replace
the corresponding later gates.

### 5. Commission a bounded prospective panel and amplify certified gains

Only after the parent certificate and carrier analysis are ready, freeze a
small initial exposure on rational parameters. For example, use twelve
addresses per newly realized parent, with a matched known-parent baseline,
and predeclare heights, expensive-call limits and candidate provenance.
Include unguided controls so the selector is judged against comparable
exposure. Twelve is a commissioning choice, not a significance threshold.

Report the four gate counts separately and compare certified quotient gains
and high-rank output per total CPU cost, including parent construction,
compaction, cold carrier traces, failures, full-cloud reconciliation and rank
replays. Preserve both the one-time construction cost and recurring marginal
cost. Retain bounded UNKNOWN outcomes. Promote productive cascades and
promising parents; permit fresh representation classes and breadth instead of
repeatedly extending a stalled fixed fibre. No daily budget is proposed.

The immediate completion milestone is **one new rational MW17 equation with
exact sections and carrier-ready marking**. A rank-28-or-higher specialization
is the next search milestone; rank 32 remains the scientific target. Broad
A1/MW16 parameter expansion is not supported by its current weak panels.

## Evidence and replay

- [Complete census certificate](../../artifacts/generated-results/elliptic-curves/det1092_pruned_rootless_j2_census_v1.json)
- [Portable full anchor packets](../../artifacts/generated-results/elliptic-curves/det1092_pruned_anchor_packets_v1.zip)
- [Prospective realization priorities and shared-core witnesses](../../artifacts/generated-results/elliptic-curves/det1092_frame_realization_priority_v1.json)
- [Exact pruned enumerator](../../elkies-k3/scripts/classify_det1092_rootless_j2_pruned.sage)
- [Witness and class checker](../../elkies-k3/scripts/assemble_det1092_pruned_census.sage)
- [Shared-core priority builder](../../elkies-k3/scripts/rank_det1092_frame_realizations.sage)

From the repository root, use a fresh work directory for enumeration:

```sh
sage -python research/elkies-k3/scripts/classify_det1092_rootless_j2_pruned.sage --regression
sage -python research/elkies-k3/scripts/classify_det1092_rootless_j2_pruned.sage --work /tmp/x1092-replay --anchors 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
sage -python research/elkies-k3/scripts/assemble_det1092_pruned_census.sage --work /tmp/x1092-replay
```

The assembler requires all sixteen completed packets and their source/input
hashes. Missing recovery of the known frame produces a contradiction artifact
and a nonzero exit, never a nonattainment theorem. Partial progress is saved
every 25 sixth vectors and after every rootless hit; resuming rechecks hashes.

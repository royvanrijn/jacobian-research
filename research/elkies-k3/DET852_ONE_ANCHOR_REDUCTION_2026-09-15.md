# Determinant 852: only one Niemeier anchor can support a rootless frame

The subsequent [global root obstruction](DET852_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md)
now excludes every MW17 fibration on this NS. Its full rational marking
and non-CM period remain valid; any rootless UNKNOWN below belongs to
this earlier gate alone.

For the [arithmetically admitted determinant852 NS](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md),
any rootless rank17 frame must arise from the D5 anchor in the Niemeier
lattice with root system2A7+2D5. The other15 complete D5 anchor cases are
excluded by exact energy bounds. **Existence in the remaining case is
UNKNOWN.** This is neither a rootless witness nor a global exclusion.

## Auxiliary and primitive glue

Use the positive even rank7 Gram

```text
 2 -1 -1 -1  1 -1  1
-1  2  0  0 -1  0  0
-1  0  2  0  0  1  0
-1  0  0  2  0  0 -1
 1 -1  0  0  2 -1  1
-1  0  1  0 -1 14  0
 1  0  0 -1  1  0 18
```

Its determinant is852 and its first five vectors span primitive D5.
The remaining two projected vectors have Gram

```text
51/4  3/4
 3/4 67/4
```

Their trace is59/2. The checker verifies the opposite discriminant form
against the actual admitted frame, constructs the index852 primitive
graph glue and checks its even unimodular positive rank24 Gram. Every
frame of this genus can therefore occur as the complement of this auxiliary
in a Niemeier lattice. Leech is excluded by the D5 roots.

## Complete anchor reduction

The retained [D5 orbit theorem](ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md)
covers16 anchors in13 rooted Niemeier classes. The checker enumerates each
anchor's orthogonal roots and applies the exact two-vector Weyl-support
bounds from the [388 proof](DET388_GLOBAL_ROOT_OBSTRUCTION_2026-09-14.md)
with the stronger large-A bounds and sharp A7 bound from the
[622 proof](DET622_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md).
All15 cases other than2A7+2D5 require total projected norm at least30,
strictly greater than59/2. They cannot have rootless complement.

The remaining anchor has orthogonal root system2A7+D5. Its sharp component
bound is2*(39/4)+6=51/2, below the available59/2 by4. The equality argument
that excluded622 therefore does not exclude852. The full projected Gram,
actual Niemeier glue, primitive auxiliary embedding and simultaneous root
avoidance still need to be checked. Matching this norm budget is not a
positive witness.

## Bounded discovery and replay

The retained rank7 source auxiliary was enumerated using even5-neighbours,
exact isometry deduplication and genus checks. The first20000-neighbour
checkpoint contained147 classes and mass7169/2560, below target21541/7680.
A checkpointed extension capped at40000 total neighbours reached152 classes
and exact target mass after34483 neighbour evaluations. The last partly
expanded parent was revisited and exact deduplication retained prior work.
Both execution limits were120 seconds and256 classes. Class151 supplies
the displayed auxiliary. The initial and resumed checkpoints, plans and
logs remain in the packet. This auxiliary genus completeness was used for
discovery, not as a premise of the sixteen-anchor reduction.

```sh
sage -python research/elkies-k3/scripts/certify_det852_one_anchor_reduction.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-root-gate-v1/certificate.json)
contains the explicit auxiliary, primitive glue, projected Gram and all16
root decompositions and bounds. Independent implementation, formal
verification, external review and novelty are unclaimed. The next task is
the exact remaining2A7+2D5 primitive-embedding problem, not an enlarged
undirected frame search. An MW17 equation remains unconstructed.

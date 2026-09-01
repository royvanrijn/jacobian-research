# Source-first optimization for the Picard-19 lattice foundry

## Outcome

The foundry objective is now ordered as

```text
same NS/T class with a catalogued MW15--17 frame
-> MW0--2 source fibration
-> rational source marking and low equation complexity
-> cheap certified elliptic-neighbour corridor
-> low-degree multisection richness on the target.
```

Rootlessness is no longer a source-search gate.  A rootless MW17 fibration is
one possible endpoint, not the object from which equation work must start.
The search starts at any exact MW15, MW16, or MW17 frame, fixes its underlying
Neron--Severi class, and searches other primitive fibrations of that same K3.

The first direct prescribed-root production slice has now replaced random
neighbour scouting as the discovery workflow.  It fixes each of the 48 exact
rank-seven auxiliaries, works in `N(3E8)`, and enumerates all-`A` complement
root systems of rank 15--17 with two or three supports.  The exact finite run
finds 97 deterministic reduced-Gram representatives of MW2 sources in 23 NS
classes: 64 have root type `A2+A6+A7`, and 33 have root type `A1+2A7`.
Every row is attached to all catalogued MW15--17 targets in its NS class.  In
particular NS0005 now has preferred-band MW2 sources attached to its 40 target
frames, and the noncyclic NS0007 class skipped by the old cyclic-glue scout
also has direct MW2 hits.

These are exact lattice results inside the declared `3E8`/all-`A`/two-or-three
support slice.  Equal deterministic reduced Grams are merged, but distinct
reduced Grams may still represent the same integral-isometry or `J2` class.
No rational marking, equation, nef corridor, or arithmetic descent is claimed.
A miss in this slice is not a theorem that the NS class lacks an MW0--2
fibration.

The subsequent full rank-16/17 pass changes the practical conclusion.  Its
thirteen rooted-Niemeier ambients produce 2,134 MW1 reduced-Gram rows and hit
all 48 foundry NS classes; no MW0 row occurs in that declared embedding cover.
Thus the present bottleneck is no longer finding a low-rank source.  It is
selecting a source whose generator has a small equation footprint and then
descending its marking.

As a positive control, the complete NS0001 rank-15/two-support run across all
sixteen D5-anchor orbits finds 19 `E7+E8/MW2` representatives in `N(3E8)`.
Their binary height forms recover the three known Kumar classes; seven
representatives have the pinned H3 form `[[21/2,3],[3,46]]` up to integral
isometry.  This is a root-and-height regression, not a full-frame identity
claim.

## Source score and proof boundary

The exact ranking artifact orders candidates lexicographically by:

1. the preferred band `MW<=2`, then MW rank;
2. number of reducible-fibre supports;
3. compatibility with a semistable all-`A` configuration;
4. expected fibre-stratum dimension;
5. minimum nonzero-section pole order;
6. known rational marking before unknown marking, then Galois orbit size;
7. expected number of additional coefficient conditions;
8. certified neighbour cost, with an unknown route ranked last;
9. the five audited low-degree multisection coordinates as a final tie-break.

Root rank, MW rank, support count, all-`A` compatibility, and the displayed
minimum pole order are exact lattice computations.  The deformation count

```text
expected fibre-stratum dimension = 18 - root rank = 1 + MW rank
```

and the resulting estimate of `MW rank` additional section conditions needed
to isolate a Picard-19 locus are heuristics until an equation ansatz is
constructed.  A one-dimensional complex lattice-polarized moduli space does
not imply a rational parameter over `QQ`.  Rational source marking, Galois
orbit size, and rational parametrization remain explicitly unknown unless an
arithmetic certificate supplies them.

The final multisection tie-break maximizes, in order, rational bisections,
genus-one bisection candidates, rational trisection candidates, genus-one
trisection candidates, and sampled low-genus quadrisection candidates.
Degree two is complete on every audited frame.  Degree three is complete on
the six-frame census below and sampled elsewhere; degree four remains sampled.
This last coordinate cannot outrank source feasibility or a certified
corridor.

Before the prescribed-root inventory, the route-aware ranking had these
leading rows:

| NS | source root type | MW | supports | all-`A` | high-rank endpoints | route |
|---|---:|---:|---:|---:|---:|---|
| NS0024 | `A3+A4+A6` | 4 | 3 | yes | 5 | certified 13-edge degree-two route to MW17 |
| NS0005 | `A1+2A3+A6` | 4 | 4 | yes | 40 | unknown |
| NS0022 | `A1+A2+A3+A6` | 5 | 4 | yes | 13 | unknown |
| NS0005 | `A1+2A3+A5` | 5 | 4 | yes | 40 | unknown |
| NS0033 | `2A2+A3+D5` | 5 | 4 | no | 40 | unknown |

The ranking is reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_lattice_foundry_sources.sage --check
```

from
[`../artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json).

The ranking now consumes 75 individual rootful-source certificates, the 97
MW2 rows of the first direct slice, and all 2,134 MW1 rows of the full-support
census, for 2,306 stored candidates in 48 NS classes.  A prescribed-root row
retains its declared finite-slice scope; the scorer does not promote it to a
complete fibration classification or invent a rational marking or route.

Restricted to that pinned 97-row MW2 inventory and classes with an audited
rootless MW17 target, the surface leaders are

| rank | NS | source root type | source MW | audited target |
|---:|---|---|---:|---|
| 1 | NS0028 | `A2+A6+A7` | 2 | NS0028-F005 |
| 2 | NS0011 | `A1+2A7` | 2 | NS0011-F002 |
| 3 | NS0022 | `A1+2A7` | 2 | NS0022-F011 |
| 4 | NS0005 | `A1+2A7` | 2 | NS0005-F008 |
| 5 | NS0001 | `A2+A6+A7` | 2 | NS0001-F001 (R17) |

All five have three all-`A` supports.  Minimum section pole, rational marking,
and neighbour cost remain open on these prescribed-root rows, so their present
order is decided by the final audited multisection coordinates, beginning with
the exact degree-two count.  It is a ranking of the declared stored inventory,
not a claim that unsearched slices cannot contain a better source.

### Minimum-section pole audit and the current equation candidate

Support count alone is a poor proxy for the source equation.  The exact
rank-one pole audit solves an affine CVP in each primitive rank-16 root lattice
for every multiple capable of improving the current norm.  Of 2,134 MW1 rows,
1,342 have primitive root lattice and pass exact rational norm recomputation
plus an independent 256-bit MPFR closest-vector audit.  The remaining 792 are
left open because their root lattice is nonprimitive and torsion/glue can
change the section lattice.

The sharp warning is `NS0022`: its one-support `A16/MW1` source wins the
support-first lexicographic score, but its generator has height `1018/17` and
minimum pole order 30.  It is not an equation-friendly source merely because
it has one reducible support.  The useful Pareto alternatives among classes
with a rootless MW17 endpoint are:

| supports | example source | exact minimum pole |
|---:|---|---:|
| 1 | NS0022 `A16/MW1` | 30 |
| 2 | NS0001 `A11+A5/MW1` | 5 |
| 3 | NS0011 `A2+A6+A8/MW1` | 2 |
| 4 | NS0007 `A1+A3+2A6/MW1` | 0 |

The current equation-balanced test case is therefore the primitive NS0011
row `NS0011-S005` from the group-a shard, with Gram SHA-256
`0e57d9632ce02df3ff4778ed25b3e66a3d7985f614cc198a884851ca6d9ca86c`.
It has semistable-compatible profile `I3+I7+I9+5I1`, MW height `352/63`,
minimum pole two, and the exact component corrections

```text
A2/I3: 0,   A6/I7: 6/7,   A8/I9: 14/9.
```

Its same-NS rootless endpoint `NS0011-F002` has 40,507 exact rational
bisection candidates and 19,023,996 exact rational trisection candidates.
This choice is a Pareto/equation judgment; it does not silently replace the
declared support-first score by a new theorem.

The first equation gate is informative but negative.  An exhaustive split
`GF(5)` Hermite scan finds seven normalized squarefree fibre models with
profile `I9+I7+I3+5I1`.  Component-adapted pole-two systems then cover all
12 nonzero-`y` and seven smooth zero-`y` infinity charts.  Exact tensor-product
evaluation checks 1,484,375 affine tuples and finds no section.  The nodal
zero-`y` point is excluded because it lies on a nonidentity `I3` component.
This is an exact local obstruction for the displayed normalized split
characteristic-five chart, not a characteristic-zero nonexistence theorem.

A bounded characteristic-seven discovery run checks 500,000 of `7^8` fibre
coefficients, finds two squarefree fibre models, and exhausts all seven of
their pole-two infinity charts (5,764,801 tuples) without a section.  Because
the fibre scan is not exhaustive, this is only a bounded negative pilot.
Other primes, nonsplit reductions, and characteristic-zero descent remain
open.

The next equation target is now NS0007 rather than a longer NS0011 prime
search.  The primitive group-c row `NS0007-S025` has root type
`A1+A3+2A6`, MW height `11/4`, and exact minimum pole zero.  Its generator
has corrections `1/2,3/4,0,0`, hence only depth-one node conditions at the
`I2` and `I4` fibres and identity-component specialization at both `I7`
fibres.  The expected semistable profile is
`I2+I4+2I7+4I1`, and the same NS class already contains the rootless MW17
frames `NS0007-F003` and `NS0007-F014`.

This candidate also has a characteristic-five obstruction, but at the earlier
fibre gate: the complete normalized split scan covers all three cross-ratios
and `3*5^8 = 1,171,875` normalized `A` polynomials.  It finds 966
Hermite-compatible sign branches and no squarefree fibre model.  A bounded
characteristic-seven scan of 100,000 coefficient/cross-ratio cases finds 26
compatible branches and again no squarefree model.  The latter is only a
pilot.  The pole-zero geometry justified testing NS0007 first, but the empty
complete characteristic-five fibre chart moves it behind the later NS0034,
NS0043, and NS0030 tests below.

The first arithmetic-support extension is also negative.  Since the `I2` and
`I4` multiplicities are unique, Frobenius fixes their supports.  After putting
them at zero and infinity, respectively, the two equal `I7` supports may form
one irreducible quadratic orbit.  Scaling the base leaves the trace-one
quadratics `t^2-t+n` and one trace-zero squareclass; these give four support
orbits over `GF(7)` and six over `GF(11)`.  The Hermite jets are evaluated in
the quadratic residue field and descended coefficientwise to the prime field.
No tangent-cone split condition is imposed at the `I4` or `I7` fibres.  The
normalization `A(0)=-3` fixes one local twist at the `I2`: it is nonsplit at
7 and split at 11, so this is not an exhaustive two-twist scan at either
prime.

A coprime-stride sample of 100,000 normalized degree-eight `A` polynomials in
each support orbit checks 400,000 rows at 7 and 600,000 rows at 11.  At 7,
336,032 signed branches pass the local square gate and 133 pass the descended
Hermite compatibility equations.  At 11 the corresponding counts are
540,420 and 24.  In both runs no compatible branch has the exact prescribed
orders `(2,7,4)` at zero, the quadratic support, and infinity; consequently
neither run produces a squarefree residual quartic or a finite-field fibre
model.  These are exact results for the displayed deterministic samples, not
exhaustive prime obstructions and not evidence that the characteristic-zero
NS0007 source does not exist.  No `GF(13)` or characteristic-zero elimination
was run.  The pinned sample artifacts have SHA-256
`a6f9d7933678a1a4f34e322f06efb5c82d9bddededc435d947d5e2a3f7437996`
and
`bc683b8e87e4d589d308a002070a898a53f4a02e86ab22d69ee301403ac216b7`.

### Later pole-zero candidates: NS0034, NS0043, and NS0030

The next exact pole-zero source is `NS0034-S008` from the group-a shard.  It
has root type `A2+A3+A4+A7`, MW height `19/8`, primitive root lattice, and
component corrections

```text
I4: 3/4,   I8: 7/8,   I3: 0,   I5: 0.
```

The same NS class contains rootless MW17 frames `NS0034-F006`,
`NS0034-F014`, and `NS0034-F018`.  The complete normalized split `GF(5)`
fibre census finds five exact `I4+I8+I3+I5+4I1` models, but the exhaustive
component-adapted `5^8` section chart on each model finds no marked section.
This is an exact obstruction only for that displayed characteristic-five
chart.

At `GF(7)`, the fixed slice `lambda=2,A8=1` is now complete.  It checks all
`7^7=823543` remaining `A` polynomials; 1,210,104 signed branches pass the
local square gate, 505 satisfy Hermite compatibility, and exactly one has the
prescribed fibre orders and squarefree residual quartic.  Its exhaustive
`7^8=5764801` section chart contains eight polynomial sections, but none has
both required depth-one I4/I8 contacts and smooth identity specialization at
I3/I5.  The fibre and section artifacts have SHA-256
`f671cbf5e08f329d118eee74956bad3a3d8c6be456a442f7ba129a661cbe8883`
and
`6f5cb66f9d2f42ed3729b168066412237858c514d657a144d039b2e65e1c1131`.
The A-eliminated 20-variable nodal-Hermite ideal is positive-dimensional
before exact-order saturation; its boundary components are not a source
family certificate.

The source-first ledger then reveals a cleaner target that was hidden by the
earlier rootless-MW17 emphasis.  `NS0043-S005` has the same three-support
semistable profile `I9+I7+I3+5I1` as the NS0011 equation gate, but its MW
height is exactly four and its minimum pole is zero.  Therefore every
component correction is zero.  Its determinant is 756, and the same NS class
has four catalogued MW15 frames `NS0043-F001` through `NS0043-F004`.

This profile lets the complete seven-model NS0011 `GF(5)` fibre census be
reused without a new fibre search.  Exhausting all `5^5` possible polynomial
X-coordinates on every model finds no polynomial section in the square local
twist.  The nonsquare quadratic twist contains 54 polynomial sections, but no
section meets the smooth identity component at all three reducible fibres.
The failures split among seven component patterns; 26 sections hit
nonidentity components at I9, I7, and I3 simultaneously.  Thus the union of
the two local twist classes is an exact obstruction for the displayed split
`GF(5)` normalization.  The corresponding square/nonsquare artifacts have
SHA-256
`8d69b9a1b3137a6c95dfb76525bcd509bd593c89ee2850d18e3ecaf07c89d5db`
and
`cf56b3fe520be6369cfeacffd33b1dc16b939fc518a5edc8c4bffc18b81fd7f2`.
The two stored `GF(7)` fibre models give the same two-twist negative section
result, but their parent fibre scan is bounded.

The next promoted semistable source is `NS0030-S001`.  It sacrifices support
count for a stronger endpoint balance: root type `2A1+A2+2A6`, MW height
`17/14`, pole zero, determinant 714, one same-NS MW15 frame, and one same-NS
MW16 frame `NS0030-F002`.  Exact projection to the five root components gives
corrections

```text
I3: 0,   one I2: 0,   other I2: 1/2,
one I7: 6/7,   other I7: 10/7.
```

The ordered five-support ansatz normalizes the I3 and identity I2 at zero and
one, the depth-one I2/I7 at `lambda,mu`, and the depth-two I7 at infinity.
A coprime-stride `GF(5)` pilot checks 100,000 normalized `A` polynomials in
each of the six ordered support pairs.  Among 600,000 rows it finds 346 signed
Hermite-compatible branches and no branch with the exact prescribed orders
and squarefree residual cubic.  This is bounded routing evidence only; its
artifact SHA-256 is
`8cd522e1f7144b18083b9c92a99f155aa595861f0fcabcf114a69b80559f3bbb`.
The remaining coefficient rows, the second local twist, characteristic-zero
descent, and a neighbour corridor remain open.

## Direct prescribed-root enumeration

The source search is reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage \
  --all-ns --ambient-label 3E8 \
  --source-support-min 2 --source-support-max 3 --all-a-only \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json
```

The exact ledger is
[`elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json).
The determinant-948 all-ambient positive control is
[`elkies-k3-lattice-foundry-prescribed-root-sources-ns0001-all-ambients-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-ns0001-all-ambients-v1.json).

The enumeration uses every certified D5 anchor in the selected ambients,
enumerates the sixth auxiliary generator modulo the residual Weyl group, and
prescribes zero Dynkin labels of rank 15--17 before solving the remaining
seventh-generator ellipsoid.  The full Niemeier lattice, including glue
cosets, is used throughout.  Acceptance recomputes the complete norm-two root
system of the saturated complement rather than trusting the prescribed face.

### Full MW0--1 prescribed-root census

The next pass removed the `N(3E8)`, all-`A`, and two/three-support
restrictions for source root ranks 16 and 17.  Four disjoint shards cover all
48 foundry NS classes, all thirteen rooted Niemeier classes admitting a D5
anchor, and all sixteen stored D5 anchor orbits.  The compact audited summary
is
[`elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-all-ambients-summary-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-all-ambients-summary-v1.json),
with SHA-256
`c5e610ac5baf12e01f86d506a6b42b6593a48f8949311eee095dfc27b55f9ad6`.
The four source ledgers are the corresponding `group-a` through `group-d`
JSON files documented in [`../REPRODUCE.md`](../REPRODUCE.md).

The exact accounting inside this declared embedding cover is:

- 2,134 shard-local deterministic reduced-Gram source representatives;
- all 2,134 have root rank 16, hence geometric MW rank 1 at Picard rank 19;
- no root-rank-17/MW0 source occurs in the cover;
- every one of the 48 NS classes has at least one MW1 source already carrying
  its catalogued same-NS MW15--17 targets;
- no `(NS, reduced-Gram digest)` identity repeats across the four ambient
  shards.

This is much stronger source availability than the bounded random scout
showed, but `2,134` is not a count of pairwise integral-isometry or `J2`
classes.  General rank-17 integral-isometry testing exhausted multi-GiB PARI
stacks, so the discovery ledger merges only equal deterministic LLL-reduced
Grams and states that limitation explicitly.

For equation design, 245 of the MW1 representatives have at most two fibre
supports and only `A`-type components, across 33 NS classes:

| source root type | representatives | NS classes |
|---|---:|---:|
| `A14+A2` | 93 | 18 |
| `A11+A5` | 46 | 9 |
| `2A8` | 38 | 8 |
| `A1+A15` | 20 | 9 |
| `A16` | 16 | 7 |
| `A10+A6` | 15 | 7 |
| `A7+A9` | 13 | 7 |
| `A13+A3` | 4 | 1 |

These support counts are an equation-ansatz convenience, not a proof that
one row is geometrically easier than another.  In particular, neither this
census nor its same-NS target attachment constructs a rational K3 marking, a
source Weierstrass equation, or a physical elliptic-neighbour corridor.  The
MW0 miss is exact only for the declared sequential D5-anchor embedding cover,
sixth-vector norm bound, and prescribed-root window; it is not a global
non-existence theorem for all K3 fibrations.

## Historical random high-rank-frame search

The source hunter now accepts any exact foundry frame, rather than requiring a
rootless start.  The first direct trials used one catalogued MW15 or MW16 frame
in each of eight Neron--Severi classes, twelve generations, beam width twelve,
60 sampled admissible Kneser neighbours per parent, and 7,981 reduced keys per
run.  The target root rank was fifteen, equivalently source MW at most two.

| starting frame | starting MW | best exact source | source MW |
|---|---:|---|---:|
| NS0002-F003 | 15 | `2A1+2A2+2A3` | 5 |
| NS0005-F001 | 15 | `A1+2A3+A6` | 4 |
| NS0011-F003 | 16 | `2A1+A2+A3+D5` | 5 |
| NS0022-F003 | 15 | `3A2+2A3` | 5 |
| NS0024-F003 | 15 | `2A1+2A3+A4` | 5 |
| NS0028-F001 | 16 | `A1+A2+A3+A5` | 6 |
| NS0032-F001 | 16 | `2A1+A3+A6` | 6 |
| NS0033-F001 | 15 | `2A1+2A2+D5` | 6 |

Every retained row is an exact primitive root-lattice/MW computation.  The
Kneser walk is discovery provenance, not an elliptic-neighbour corridor, so
its edge count is not used as equation cost.  The negative result is complete
only for the declared deterministic beams and samples.

For example, replay the strongest new high-rank-start row with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage \
  --ns-id NS0005 --target-frame-id NS0005-F001 \
  --generations 12 --beam 12 --samples-per-parent 60 \
  --primes 3,7,11,13,17,23 --seed 20262906 \
  --target-root-rank 15 --allow-below-target \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-v1.json \
  --root-adapted-frame-output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-root-adapted.txt
```

The other seven JSON artifacts use the same
`elkies-k3-lattice-foundry-nsNNNN-mw2-source-from-high-mw-scout-v1.json`
naming pattern and record their seed, admissible prime list, generation
accounting, and visited-key count.

The prescribed-root implementation above supersedes lengthening these random
beams for source discovery.  The old artifacts remain exact bounded
provenance.  Only after a direct candidate passes rational marking and
source-equation gates should a physical neighbour corridor be optimized.

## Low-degree multisection spectrum

Proposition F5 in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
reduces rootless degree-`d`, genus-`g` section-nonnegative classes to coset
minima in `M/dM`, with threshold

```text
2*d^2 - 2*g + 2.
```

The degree-two calculation is complete through norm ten.  It exactly
reproduces the published R17 count of 39,120 geometrically rational bisection
translation orbits and finds several foundry endpoints with more; the largest
in this nine-frame batch is NS0032-F011 with 41,421, about 5.9 percent above
R17.  NS0028-F005 has 41,376 and NS0033-F026 has 40,912.  This confirms that
R17 is not extremal even for the exact minimal rational-bisection coordinate.

Genus-one bisection counts are exact lattice-candidate counts, but global
nefness, irreducibility, and arithmetic descent are not yet certified.

The primary complete degree-three run covers the current top five from the
pinned MW2 source ranking above.  For every frame it visits all
`3^17 = 129,140,163` translation cosets.  Inversion reduces the CVP work to
`64,570,082` representatives per frame, so the certificate accounts for
645,700,815 cosets and 322,850,410 inversion representatives.

| rank | frame | rational trisection cosets | genus-one trisection cosets | maximum coset minimum |
|---:|---|---:|---:|---:|
| 1 | NS0028-F005 | 19,645,256 | 34,294,400 | 26 |
| 2 | NS0011-F002 | 19,023,996 | 33,978,764 | 26 |
| 3 | NS0022-F011 | 18,774,826 | 33,788,528 | 26 |
| 4 | NS0005-F008 | 18,446,258 | 33,705,930 | 26 |
| 5 | NS0001-F001 (R17) | 18,024,616 | 33,484,468 | 26 |

Here is the complete minimum-norm spectrum; every column sums to `3^17`.

| norm | NS0028 | NS0011 | NS0022 | NS0005 | R17 |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 1 | 1 |
| 4 | 2,394 | 2,480 | 2,536 | 2,560 | 2,622 |
| 6 | 48,832 | 50,570 | 51,322 | 52,168 | 53,344 |
| 8 | 421,150 | 435,848 | 445,570 | 449,780 | 460,080 |
| 10 | 2,261,764 | 2,342,806 | 2,376,486 | 2,416,496 | 2,472,628 |
| 12 | 8,703,488 | 9,004,264 | 9,206,870 | 9,288,622 | 9,495,786 |
| 14 | 22,973,754 | 23,593,438 | 23,819,764 | 24,157,244 | 24,568,586 |
| 16 | 36,883,630 | 37,271,412 | 37,465,808 | 37,559,396 | 37,786,028 |
| 18 | 34,282,756 | 33,967,910 | 33,779,846 | 33,701,210 | 33,481,080 |
| 20 | 19,644,508 | 19,023,834 | 18,774,542 | 18,446,058 | 18,024,296 |
| 22 | 3,905,494 | 3,436,584 | 3,208,452 | 3,061,708 | 2,792,004 |
| 24 | 11,644 | 10,854 | 8,682 | 4,720 | 3,388 |
| 26 | 748 | 162 | 284 | 200 | 320 |

The exact degree-three ordering agrees with the degree-two order on these five
surfaces, and NS0028 is the clear leader in both displayed trisection
coordinates.  The 256-coset pilot nevertheless misestimated the magnitudes
substantially: it observed 49 rational hits for NS0028 where the exact expected
count per 256 is 38.944, and 83 genus-one hits for NS0022 where the exact
expected count is 66.980.  Thus a small sample was useful for triage but is not
reliable for close quantitative comparisons.

A separate completed six-frame certificate retains the five leaders from the
older pre-prescribed-root route-aware ledger plus R17.  Its three additional
surfaces have exact `(rational, genus-one)` counts
NS0033-F026 `(19,287,006, 34,122,336)`, NS0002-F007
`(18,771,452, 33,869,098)`, and NS0024-F005
`(18,133,774, 33,469,500)`.  NS0028 exceeds even NS0033, the richest surface in
that older batch.

Every returned CVP candidate has its norm recomputed with the integral Gram
matrix.  A deterministic subset of about 15,760 residues per frame is
independently repeated using 256-bit MPFR GSO arithmetic; the largest recorded
distance-to-integral-norm discrepancy is below `5.0e-14`.  This is a complete
computational lattice census, not a formal verification of the floating CVP
branch decisions.

The artifact is reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage \
  --sample-count 256 --height-slack 4 \
  --frame-id NS0001-F001 --frame-id NS0002-F007 \
  --frame-id NS0005-F008 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0024-F005 \
  --frame-id NS0028-F005 --frame-id NS0032-F011 \
  --frame-id NS0033-F026 --check
```

and stored at
[`../artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json).

Run the complete degree-three census, or check its completed checkpoint without
rerunning the CVP enumeration, with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0028-F005 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0005-F008 \
  --frame-id NS0001-F001 \
  --workers 8 --chunk-size 1000000 --float-type dd \
  --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0028-F005 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0005-F008 \
  --frame-id NS0001-F001 \
  --workers 8 --chunk-size 1000000 --float-type dd \
  --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json \
  --check
```

The primary exact certificate is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json).
Its stable whole-file SHA-256 is
`8be0e881f5c170366dada6319aed9a09fed689eacc032fcaf5ee70878d735fd0`.
The older-batch certificate is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-top5-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-top5-v1.json),
and each has its complete chunk checkpoint in the adjacent `.partial` file.
The old degree-four entries remain exact-CVP results only for the declared 256
sampled cosets and must not be promoted to a complete census.

Multisection richness is a secondary discovery coordinate, not a specialization
rank theorem.  The R17 positive controls already show that 39,120 bisections
can leave an extreme specialization largely invisible.  The geometric
motivation for retaining this coordinate is the relation between multisections
and rank jumps studied by Garbagnati--Salgado, while the use of alternative
elliptic fibrations to obtain rank jumps is consistent with Salgado's earlier
two-fibration method:

- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).
- C. Salgado,
  [*On the rank of the fibers of rational elliptic surfaces*](https://arxiv.org/abs/1307.3994).

The same-K3 fibration search is grounded in the Kneser--Nishiyama framework;
the bounded foundry catalogue is not a replacement for a complete fibration
classification:

- K. Nishiyama,
  [*The Jacobian fibrations on some K3 surfaces and their Mordell--Weil groups*](https://doi.org/10.4099/math1924.22.293).
- V. Braun, Y. Kimura, and T. Watari,
  [*On the Classification of Elliptic Fibrations modulo Isomorphism on K3 Surfaces with large Picard Number*](https://arxiv.org/abs/1312.4421),
  especially Section 4.1 for the surjection from primitive auxiliary
  embeddings to frame-isometry classes.
- M.-J. Bertin and O. Lecacheux,
  [*Elliptic Fibrations of a certain K3 surface of the Apéry--Fermi pencil*](https://doi.org/10.5802/pmb.44),
  for the warning that a non-root auxiliary must be embedded in the full
  Niemeier lattice rather than only its root lattice.
- I. Shimada,
  [*On elliptic K3 surfaces*](https://arxiv.org/abs/math/0505140).

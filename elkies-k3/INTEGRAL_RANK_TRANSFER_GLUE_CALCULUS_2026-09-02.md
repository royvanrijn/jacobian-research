# Target-directed fibration hopping: an equation-free glue calculus (2026-09-02)

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CALCULUS 7eeeeaa80d9b2bf3 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-THETA-CONVOLUTION 5ebbd3d242fdb3db -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CORE-GENERATION d0d78c49b44f55ac -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-REVERSE-THETA eee16ce986ec0a1f -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-WEIL-COMPRESSION 34d2abea91a265f4 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MODULAR-DIMENSION-SIEVE 9622c6eb4d8522bd -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-GENERATION 9a7a1e01cb22f62e -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-CONTROLS 3cbde45fb2cb0f17 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-DIRECTED-Q80 80de8b6727cd3409 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-BIRTH-DEATH a755a3956c4c97cb -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-GRAPH-REACHABILITY e02f950eba79b32a -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-ROOT-SYSTEM-SIGNATURE d32b35b66a35627c -->
<!-- status-consumer: EC-K3-NS0024-INVERSE-ADE-MUTATION 5c56f07d14129837 -->
<!-- status-consumer: EC-K3-INTEGRAL-CHARACTER-GLUE 0b76d65366279037 -->
<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL 648ec884ce7152bb -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-INTEGRAL-GLUE 52de13c8443f2b7d -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-BRIDGE-PREDICTOR-BENCHMARK 3127e24cc505f646 -->

## Outcome

Terminology and novelty boundary: **integral rank transfer** remains the
project label.  A same-surface move is a change of primitive `U`-embedding,
or a fibration hop, and its rank change is the classical Shioda--Tate rank
balance.  Nikulin graph gluing, Kneser neighbours, Kneser--Nishiyama frame
classification, Weyl reduction, and explicit fibration hopping are
established.  The contribution of this note is the target-directed inverse
use of those tools and the determinant-specific exact computations below.
See the
[`literature and novelty map`](LITERATURE_AND_NOVELTY_MAP_2026-09-03.md)
for claim-level provenance.  “No explicit antecedent located” is used for
plausibly new algorithms; no priority claim is made.

The proposed experiment has a sharp first answer.

There are three recurring integral operations in the successful examples:

1. replace the primitive hyperbolic plane `U` inside a fixed `NS` and take
   the new saturated orthogonal frame;
2. glue an auxiliary lattice to a frame by the graph of a discriminant-form
   anti-isometry, then take the complementary lattice;
3. split into rational finite-group isotypic components and recover the integral
   lattice by `|G|`-primary saturation glue.

The first operation accounts for every edge in the H3, Q80, NS0024, and
Golay-720 corridors.  It does **not** change the finite quadratic form:
because `U` is unimodular, the frame discriminant form is fixed along the
whole corridor.  The finite form is an admissibility and genus invariant, not
the object which records which roots disappeared.

There is now an exact local replacement statement behind these pivots.  For two
frames `W,W'` in the same `NS`, their common bridge core

```text
K=W intersect W'=(U+U')^perp(-1)
```

contains precisely the roots that survive the hop.  Each frame is a finite
graph-glue extension of `K+C_i`, where `C_i=K^perp` inside that frame.  On all
42 selected H3, Q80, NS0024, and Golay-720 edges, `K` has rank 15, both `C_i`
have rank two, and the two glue groups are cyclic of the same order and
project onto the full bridge discriminant.  Thus all 42 selected hops are
exact rank-two cyclic bridge replacements.  This 42-edge corpus result is a
new computation; the common-core and graph-glue formalism is a tailored
application of Nikulin's established theory.  The certificate is
[`elkies-k3-integral-rank-transfer-bridge-reglue-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-reglue-v1.json).

The maximality condition also inverts core generation at finite-form level.
For any proposed binary bridge `C` and target frame form `q_W=-q_NS`, the
rank-15 core genus is forced by

```text
q_K = q_W orthogonal_sum (-q_C),
det K = det W * det C.
```

Within a declared finite bridge universe, the decorated signature consisting
of `q_K` and every discriminant-coset theta coefficient through norm two then
determines the complete zero-support completion spectrum.  The resulting
algorithm enumerates the forced core genera first, rejects any class already
containing a root, and only then enumerates compatible graph anti-isometries.
The exact replay checks the splitting on all 84 old/new bridge presentations.
Its cheapest gate rejects all 277 primitive cores in the complete held-out E6
source shell, while the four known positive terminal cores have minimum four
and admit five rootless bridge classes among fourteen.  This supplies a
bounded core-first calculus, but still no universal bridge-determinant bound,
general fast genus enumerator, or speedup theorem.  See Theorems H0b--H0c in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) and
the generated
[`core-generation certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-core-generation-v1.json).

The theta condition itself can now be inverted as well.  A bridge and graph
define a finite reverse mask

```text
F_(C,H)={(a,2-nu):(a,b) in H and theta_C(b,nu)>0},
```

and a core completes rootlessly exactly when its theta coefficients vanish on
that mask.  The terminal masks contain 18--87 oriented cells and reduce under
`a <-> -a` symmetry to a nonredundant antichain of 14 masks with 10--44 cells.
Lazy exact coset
queries reproduce all 28 graph decisions; the complete core tables, computed
only afterward as a truth check, contain 4,418--22,579 occupied cells.  For a
rootless rank-15 core every coefficient through norm two is at most 30 by the
orthoplex bound, so the bridge masks cut out a finite explicitly enumerable
set of allowed low-order signatures.  The remaining problem is no longer
unrestricted genus enumeration but inverse theta realization: construct only
lattices whose vector-valued theta series has an allowed signature.  The
exact mask certificate is
[`elkies-k3-integral-rank-transfer-reverse-theta-masks-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-reverse-theta-masks-v1.json).

The modular inverse-theta step also has an exact quotient before any modular
forms are enumerated.  Under the forced splitting
`A_K=A_W orthogonal-sum (-A_C)`, every mask class has zero `A_W` coordinate.
Averaging over `O(q_W)` commutes with the Weil representation and leaves that
zero component unchanged.  Hence masked modular feasibility lives in

```text
C[A_W]^{O(q_W)} tensor C[A_C]
```

and then in its theta-symmetric part, not in all of `C[A_K]`.  For Golay-720,
H3, NS0024, and Q80 this reduces the coefficient dimensions respectively
from `16,560, 44,556, 181,450, 21,804` to
`864, 5,760, 24,960, 2,880`.  Exact good-prime `S,T` closure shows that the
cyclic Weil submodule generated by the zero class already fills the complete
orthogonal-orbit quotient in every control, so that particular refinement
cannot compress further.  This is a factor-of-`7.27`--`19.17` exact reduction
of the next modular prefilter, not yet an enumeration or lattice realization.
See Theorem H0e and the
[`zero-orbit compression certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-weil-compression-v1.json).

Running the invariant Riemann--Roch dimension calculation changes the next
priority.  The four compressed weight-`15/2` modular spaces have dimensions
`476, 3,121, 13,488, 1,563`, with cusp dimensions
`472, 3,120, 13,485, 1,562`.  A mask imposes at most its 10--44 coefficient
conditions, so every terminal mask leaves a cusp kernel of dimension at least
461.  Linear modularity is therefore much too underdetermined to be the
desired classifier.  This does not settle the affine constant-term condition
for rejected masks.  It says that the next useful sieve must enter the
arithmetic theta cone: bounded integral nonnegative coefficients, local
density restrictions, and ultimately lattice realization.  The same exact
certificate records this Corollary H0f audit.

The first constructive inverse-theta step now bypasses a complete core-class
enumeration in one control.  For the Golay-720 target and its class-2 binary
bridge, the finite form generates a unique rank-15 genus and that genus
supplies its own canonical representative, which is rootful.  A fixed-seed
good-prime Kneser beam uses root count and then occupied reverse-mask cells as
its exact score.  It reaches a zero-mask core after seven neighbour steps and
34,571 unique raw neighbours.  The resulting minimum-four core is not
isometric to the historical Golay core, but its order-23 completion is
rootless and isometric to the declared rank-17 Golay target.  Thus the method
has generated a genuinely new compatible core class from the forced genus,
not merely recognized an input harvested from a successful corridor.

The short certificate replays all seven neighbours, the empty class-2 mask,
and the final glue; the full mode reruns the bounded beam.  See Theorem H0g
and the
[`mask-aware core-generation certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-masked-core-neighbors-v1.json).
This is a proved construction for one control, not a genus-completeness,
success-probability, uniform-running-time, or speedup theorem.

Prospective application to the other terminal genera gives a mixed but
substantive answer.  Exact good-prime paths from their canonical 280-root
genus representatives construct nonhistorical minimum-four cores for H3 and
NS0024.  H3 completes through bridge class 2 to the declared determinant-948
rootless target.  NS0024 completes through bridge class 4 to a new
determinant-950 rootless rank-17 isometry class with the target discriminant
form.  The corresponding Q80 path reaches a rootless core but retains one
sign-paired obstruction, hence two occupied cells, on its only viable bridge
mask.

These experiments also corrected the discovery rule.  A bridge containing a
root contributes the compulsory occupied cell `(0,0)` and must be removed
before cores are ranked; bridge class 1 is excluded in all three controls.
On NS0024, preserving distinct occupied support signatures finds the new core
after 7,477 neighbours.  The same corrected eight-generation rule examines
42,300 neighbours without a hit on H3 or Q80.  Therefore masked-support
diversity is useful but not a complete surrogate for `Sigma_2` realization or
core isometry class.  See Corollary H0h and the
[`prospective control certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-masked-core-controls-v1.json).

The Q80 near miss is now resolved by passing from masked cells to their
physical dual-vector witnesses.  If the good-prime neighbour is defined by an
isotropic line `ell=<y>`, an old dual vector `x` survives exactly when
`<x,y>=0 mod p`.  Lines nonorthogonal to every forbidden witness therefore
remove the complete old defect before a neighbour is built.  New forbidden
vectors can still appear, and this replacement is the dominant behavior:
among 1,397 rootless one-step Q80 neighbours which kill all four initial
witnesses, none has zero mask.

An isometry-diverse directed beam transports the two-cell defect through
physical witness populations `4 -> 6 -> 4 -> 4 -> 0`.  It reaches zero in
four directed steps after constructing 30,228 distinct neighbours.  Appending
these steps to the canonical eight-step prefix gives a new rootless Q80 core
whose determinant-948 rank-17 completion is exactly the alternate Q80 frame:
it has 1,313 norm-four pairs and automorphism-group order four, is explicitly
integrally isometric to alternate Q80, and is not isometric to published R17.
The previous `declared_target_frame` label referred to published R17, so its
failed isometry test did not indicate a third class.  Exact local symbols at
`2`, `3`, and `79` independently place the completion in the common target
genus.  Thus the directed core flow connects the Q80 near-miss region around
the published target to the other mass-complete rootless `J2` class.  The
exact survival law and path are Theorem H0i.  The result identifies the
missing state variable: masked support must be decorated by the incidence of
its individual witnesses against isotropic neighbour lines.

The replacement side is also explicit.  For the standard lift with
`y^2=0 mod 2*p^2`, put

```text
K_y dual={x in K dual:<x,y>=0 mod p}.
```

Then the child dual is the disjoint union of
`K_y dual+j*y/p`, `0<=j<p`.  The `j=0` layer is exactly the old-survivor
theorem; nonzero layers contain every possible birth.  Under the canonical
prime-to-`p` identification of discriminant groups, the class of
`x+j*y/p` is the old class of `x`, so every masked child coefficient is a
finite affine-CVP sum.  The exact Q80 replay predicts all four child witness
sets and their complete norm-at-most-two theta profiles before construction,
then independently recovers the same profiles and vectors from the children
afterward.  This is Theorem H0i.1 and its
[`birth--death certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-q80-defect-birth-death-v1.json).
The counted abstract `Sigma_2` table alone still omits the required physical
line incidence, and no runtime speedup theorem is claimed.

Complete small-genus reachability graphs now calibrate the global dynamics.
In three mass-closed even ternary genera, three good primes are analyzed both
separately and in every pair/triple union.  Exact SCC labelling finds a closed
two-state directed `3`-trap at determinant 112 and two singleton `3`-traps at
determinant 316, despite unrestricted `3`-connectivity.  Every tested
two/three-prime union reaches zero from every state.  Exhaustive set selection
shows that one of primes `5`, `7`, `11`, or `13` already suffices in each
control, with the exact minimizers depending on the genus.  All three genera
have one proper spinor genus, so the fixed-prime traps are finer than spinor
separation.  The resulting candidate global invariant is the marked
spinor/level component for eventual access, refined at small primes by the
physical-witness incidence profile.  This is Theorem H0i.2 and the
[`small-genus reachability certificate`](../artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v2.json),
not yet a marked rank-15 theorem.

Those same physical witnesses now transfer more than zero versus nonzero
support.  Retaining every completion root `r=k+c` and the pairwise metric
`<r,r'>=<k,k'>+<c,c'>` recovers the complete signed root graph, ADE
decomposition, rank, and root discriminants.  Retaining their coordinates in
the completed frame additionally recovers the primitive closure and exact
Mordell--Weil torsion quotient.  The metric without the marked embedding does
not determine those last two invariants.  The exact NS0024 controls recover
`D5+E8`, `3A1+A2`, and rootless stages, including twelve roots spread over
five nonzero graph-glue labels; Q80 controls recognize `4A1` and `A1`.  See
Theorem H0k and the
[`metric physical-witness certificate`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-root-system-signature-v1.json).

The forward metric signature now has an exact inverse at a fixed good-prime
line.  Given `K,C,H` and `y`, old physical roots impose modular survival or
death conditions `<k,y>=0` or nonzero modulo `p`; every possible birth is a
vector in one of the finite affine shells

```text
M+r_a+k_a+j*y/p,       k^2=2-c^2,
```

joined to a bridge vector `c` through `H`.  Requiring the resulting complete
metric to equal the desired ADE root system, and every unselected shell to be
empty, is necessary and sufficient.  For an abstract ADE label this is a
finite disjunction over marked witness templates; for a selected template it
is a direct conjunction of modular incidence and affine-CVP constraints.

The first NS0024 `p=17` edge supplies the initial exact control.  Of the 140
`D5+E8` root lines, six prescribed modular forms vanish and 134 are nonzero;
the six survivors span `3A1+A2`, while exhaustive nonzero-layer and order-191
graph-glue enumeration produces no birth or extra root.  Independent child
construction afterward gives the same physical root set.  See Theorem H0l
and the
[`target-root-system certificate`](../artifacts/generated-results/elkies-k3-ns0024-inverse-ade-mutation-v1.json).

Literal Nikulin discriminant-form glue occurs in the Niemeier complement
constructions.  A closely related finite-index saturation quotient occurs in
E6 `2+2`.  Its original height blocks are rational; the census multiplies the
pairing by 12 to obtain an even integral bookkeeping lattice before computing
`A_L`, `q_L`, and the corresponding isotropic subgroup.  Before that explicit
rescaling the quotient must not be called a Nikulin glue subgroup.  It is the
cleanest local character model:

```text
L+ = <P,Q>,              det(L+) = 4/3,
L- = <T1,T2>,            det(L-) = 52,
2R1=P+T1,  2R2=Q+T2.
```

The pure character sum has determinant `208/3`; adjoining the two half-sums
has index four and gives determinant `13/3`.  Thus the exact saturation
quotient is `(Z/2)^2`.  Exhausting all six graph isomorphisms between the two
`F_2^2` eigensubgroups shows that they form one integral isometry class after
the factor-12 scaling.  The E6 `2+1` case has exactly two graph-glue types:
the actual index-one sum and one alternative index-two class, the latter with
smaller minimum.  See
[`elkies-k3-integral-character-glue-calculus-v1.json`](../artifacts/generated-results/elkies-k3-integral-character-glue-calculus-v1.json).

The requested equation-free tuples and all 42 selected same-`NS` corridor
edges are stored in
[`elkies-k3-integral-rank-transfer-glue-census-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-glue-census-v1.json).
The local bridge and character statements, and the mass-complete decorated
neighbour algorithm, are now theorems recorded in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
They do not assert uniform small-move connectivity or an equation-level lift.

## The tuple convention

Each experiment is represented in the order

```text
(L, G, {L_chi}, A_L, q_L, glue subgroups).
```

Here `L` is one integral lattice, and every `L_chi`, `A_L`, `q_L`, and glue
entry is typed relative to that same `L`.  For same-surface corridors it is
`NS`, so trivial `G` has `L_triv=NS`; the varying frames `U_i^perp(-1)` are
recorded separately.  For E6 it is the pure MW-character sum with its height
pairing multiplied by 12.  For the bisection experiments it is the certified
orthogonal character sublattice, because the full Mordell--Weil lattice after
base change is not known.

The `glue_subgroups` entries carry a context label.  In particular, a
Niemeier graph subgroup records a construction of a related complement; it is
not asserted to be an overlattice subgroup of `A_NS`.

Finite quadratic forms are normalized in `Q/2Z` directly from the pinned Gram
matrices.  A displayed list such as `(2,2,6)` means the nontrivial Smith
invariants of `A_L`.

## Exact lattice lemmas

### 1. A primitive `U` has no glue defect

Let `U` be a copy of the unimodular hyperbolic plane in an integral lattice
`L`, and put `W=U^perp`.  Then

```text
L = U direct_sum W.
```

Indeed, unimodularity makes the orthogonal projection `L -> U dual = U`
integral, its kernel is `W`, and subtracting the projection gives the claimed
direct sum.  In particular the embedding is automatically primitive.
Consequently

```text
A_W is isomorphic to A_L,       q_W is isomorphic to q_L.
```

For the positive elliptic frame convention `NS=U+W(-1)`, the sign is reversed.
Thus changing `U` can change the integral isometry class and root system of
`W`, but cannot change its finite quadratic form.

This is the lattice core of the Kneser--Nishiyama classification; see
[Nishiyama](https://www.jstage.jst.go.jp/article/math1924/22/2/22_2_293/_article)
and Sections 12.9--12.21 of
[Schuett--Shioda](https://arxiv.org/abs/0907.0298).

### 2. An overlattice cannot annihilate an existing root

If `L0` is contained in an overlattice `L1`, every norm-two vector of `L0` is
still a norm-two vector of `L1`.  Therefore adjoining one isotropic glue
subgroup can add roots but cannot remove an old root.

Root annihilation in the present examples must mean one of:

- a root leaves `U^perp` when `U` changes; or
- old cosets are removed and new cosets are added in a genuine **reglue**.

This rules out a naive calculus whose atomic command is only “enlarge by an
isotropic subgroup.”  For an **even** lattice, Nikulin's quadratic-overlattice
formula remains the bookkeeping law:

```text
H isotropic in (A_L,q_L),
[L_H:L]=|H|,
(A_LH,q_LH)=(H^perp/H, q_L restricted),
|det(L_H)|=|det(L)|/|H|^2.
```

See Proposition 1.4.1 of
[Nikulin](https://www.mathnet.ru/eng/im1677).

### 3. Finite-character glue is supported on primes dividing `|G|`

For an involution `sigma` of an integral lattice `L`, put

```text
L+ = L intersect ker(sigma-1),
L- = L intersect ker(sigma+1).
```

The two lattices are orthogonal.  Moreover, for every `x` in `L`,

```text
2x = (x+sigma*x) + (x-sigma*x),
```

so the quotient formed from the **full** lattices `L+` and `L-` is killed by
two.  More generally, decomposition into rational isotypic components has a
finite saturation quotient supported only at primes dividing `|G|`, because
the rational central idempotents have denominators dividing `|G|`.  This does
not constrain the index of a sum of merely chosen or visible sublattices.
The E6 `2+2` half-sums are an exact finite-index realization in the rational
MW height lattice; E6 `2+1` has index one.

This rational character decomposition and its finite-index integral defect
also appear in explicit base-change constructions; see
[Kloosterman](https://arxiv.org/abs/math/0502017) and
[Hulek--Schuett](https://arxiv.org/abs/0912.0608).

## The normalized successful examples

| experiment | `L` and finite form | `G` and character lattices | exact glue operation | status boundary |
| --- | --- | --- | --- | --- |
| historical H3 / J2 / R17 | `L=NS_H3`, rank 19, determinant 948; `A_NS=Z/948`, normalized `q_NS=diag(7/4,4/3,2/79)` | trivial, hence `L_triv=NS`; the separate positive frame family has rank 17 and the opposite finite form | `U` split has glue zero; the J2 classifier instead glues the rank-seven auxiliary and either rootless frame by an order-948 graph subgroup in `N(2A7+2D5)` | exactly two rootless J2 classes; no J1 claim |
| Q80 low-q | the same `L=NS_H3` and finite form | trivial, hence `L_triv=NS` | ten primitive-`U` changes, all determinant-one transports; no edgewise finite-form change | alternate rootless J2 equation remains open over characteristic zero |
| NS0024 | `L=NS0024`, rank 19, determinant 950; `A_NS=Z/950`, normalized `q_NS=diag(1/2,2/25,4/19)` | trivial, hence `L_triv=NS` | thirteen primitive-`U` changes; source also arises by cyclic order-950 graph glue in `N(A15+D9)` | exact lattice route, equation and rational descent open |
| Golay-720 | `L=NS`, rank 19, determinant 720; `A_NS=(2,6,60)` with exact form in the census | trivial, hence `L_triv=NS` | order-720 maximal graph glue in `N(24A1)`; an exact three-generator graph among 96 anti-isometries gives the rootful `N(4A5+D4)` companion; six later changes of `U` | binary-octad construction and route exact; rational source equation open |
| E6 `2+1` | `L=12(MW_+ direct_sum MW_-)`, rank 3, determinant 1536; exact `A_L,q_L` in the census | `C2`; scaled blocks have ranks `2+1`; unscaled determinants `4/3`, `2/3` | character glue zero, index one | ambient `NS` has determinant 24 and roots `2E6+A2`; same-NS rootless MW17 is impossible by the determinant bound |
| E6 `2+2` | `L=12(MW_+ direct_sum MW_-)`, rank 4, determinant 1,437,696; exact `A_L,q_L` in the census | `C2`; scaled blocks have ranks `2+2`; unscaled determinants `4/3`, `52` | two half-sums give `(Z/2)^2`, index four; this is literal even glue only in the scaled model | ambient `NS` has determinant 78 and roots `2E6+A1`; ordered incidence has `2+2`, rational quotient `1+1`; a complete Niemeier residual-rank obstruction proves its entire frame genus rootful |
| R17 genus-one bisection cover | `L=R17(2)+<16>`, the orthogonal character sublattice, rank 18 | `C2`; invariant rank 17 and one certified anti-invariant height-16 line | `2R=tau+T` gives a certified `Z/2` half-sum extension of index two; only further saturation is unknown | rank at least 18, not an exact total-rank or full-MW claim |
| paired rational R17 norm-ten bisections | `L=R17(4)+<24>+<24>`, the orthogonal character sublattice, rank 19 | `V4`; two distinct nontrivial characters give the height-24 lines | the two independent half-sums give `(Z/2)^2`, index four; further saturation and the product character are unknown | all 765,167,640 pairs are geometrically genus-one `V4` covers; rational points are separate |
| norm-12 / norm-8 simultaneous split | on the norm-12 cover, `L=R17(2)+<16>`, rank 18, determinant 1,988,100,096; exact `A_L,q_L` in the census | `C2`; invariant `R17(2)` and primitive anti-invariant `<16>` | `2R=tau+T` gives the full `Z/2` graph saturation in the displayed rational span; the saturated visible lattice has determinant 497,025,024, minimum 8, and no roots | at `t=1/25`, `0x103b2` gives an independent direction and rank at least 18; total cover anti-rank and specialization rank upper bound remain open |

The H3 entry is deliberately not written as three arrows.  `J2` means an
integral frame-isometry class.  The historical H3 corridor, the Niemeier J2
classification, and identification of the published R17 endpoint are three
different proof layers.

Here “Golay” means the binary-octad determinant-720 construction and its
rootful source corridor.  The ternary-Golay `N(12A2)` foundry backend is a
large census (151 frame classes after local deduplication), not one of these
named transition chains.  It should be expanded as a family of tuples rather
than misrepresented by one aggregate tuple.

## What happened at the same-`NS` fibration hops

Across the four selected corridors the census contains 42 exact edges:

```text
old fibre degree 2:  42 / 42
q=4:                 27
q=6:                 12
q=8:                  2
q=24:                 1
q=4 or 6:            39 / 42
nonunimodular NS transports: 0
finite-form changes:          0
```

The common operation is therefore:

```text
(NS,U,W,R,MW)  ->  (NS,U',W',R',MW')
```

with a determinant-one integral marking transport and

```text
rank(R')-rank(R) = -(rank(MW')-rank(MW)).
```

For all 42 normalized marked edges, the finer operation is

```text
W  = overlattice(K+C_old, H_old)
  -> overlattice(K+C_new, H_new) = W',
```

where `rank(K)=15`, both bridges have rank two, and `H_old,H_new` are cyclic
maximal bridge graphs of equal order.  The observed orders are

```text
15, 23, 31, 47, 63, 95, 119, 127, 143, 159,
191, 215, 303, 359, 799, 991, 1231, 1535, 2447, 3231.
```

This proves that the geometric route labels `q=4,6` are not the finite glue
orders.  They measure the chosen isotropic divisor, while the bridge core
determines the actual reglue module.

The H3, Q80, NS0024, and Golay chains in the artifact retain every root type,
rank change, old-fibre degree, and `q`.  Most selected edges shed one root;
some shed several.  Their overwhelming concentration at old degree two and
`q=4,6` is strong routing evidence, not a bounded-connectivity or monotonicity
theorem.  The route parameter `q=ab` is also **not** a Kneser neighbor prime;
these labels must not be conflated.

## What happened at the literal glue replacements

### Maximal graph glue and complementary frames

Let `N` be even and unimodular, let `K` be primitive, and require
`W=K^perp_N` with complementary rank (equivalently `K=W^perp_N`).  Then their
finite forms are anti-isometric.  The ambient lattice is obtained from `K+W`
by the maximal isotropic graph

```text
H_phi = {(x,phi(x)): x in A_K}.
```

This is the exact H3-J2 and Golay construction.  Mutating the primitive
embedding, complement class, or Niemeier ambient can change the root system
of the saturated complement.  Changing only the anti-isometry for fixed
abstract `K,W` changes the ambient gluing, not the intrinsic roots of `W`.
Rootlessness is checked in the actual complement, not inferred from the
finite form.

The H3 computation is particularly decisive: the published R17 and alternate
Q80 frames have the same rank, determinant, minimum, and discriminant form,
but are not integrally isometric.  The complete Niemeier enumeration finds
exactly these two rootless J2 classes.  Conversely, 65 old rootless
determinant-948 seeds fail the 2-adic and 79-adic form and lie in a different
genus.  Both facts show why the calculus needs the integral lattice **and**
its finite form.

### Character half-saturation

E6 `2+2` gives a completely identified nontrivial character saturation
quotient.  The census's factor-12 model makes the height blocks even integral
and verifies the finite-form statement there.  The unscaled rational height
quotient should still be described as saturation, not silently as Nikulin
glue.  The R17 bisection identity

```text
R + R' = tau,       T=R-R',       2R=tau+T
```

has the same half-sum shape.  It certifies an order-two subgroup and an
index-two extension of the orthogonal character sum; what remains unknown is
whether that extension exhausts saturation and the full anti-invariant
lattice.  For a pair of rational norm-ten bisections, the separate character
coordinates make the two half-sums independent, giving a certified
`(Z/2)^2` subgroup and index four.  Product-character rank and further
saturation remain open.

The rational Golay `3I6` model is the warning example.  A half-section and
rational 3-torsion give index `2*3=6`, changing determinant `720` to `20`.
The equation is simple and the geometric ranks look right, but it belongs to
the wrong integral `NS` class.

## The inverted calculus

The data support the following constructive program.

### State

Use a marked integral state

```text
S=(NS,U,W,R,A_W,q_W,H_root,minimum,physical roots,root metric,saturation indices).
```

The root glue `H_root` records the primitive closure of the fibre-root lattice
inside the frame.  For a finite symmetry also store the integral character
lattices and their `|G|`-primary graph glue.

### Move types

1. **U-pivot.** Choose a primitive nef isotropic `F'` and a mate defining a
   primitive `U'`; take the exact saturated complement `W'`.
2. **Complement reglue.** Change a primitive auxiliary embedding or graph
   anti-isometry in a Niemeier lattice and take the new saturated complement.
3. **Character saturation.** Compute the full integral isotypic lattices (or
   specify and verify a denominator-clearing integral model), then enumerate
   the finite `p`-primary extensions for primes dividing `|G|`.  Invoke
   isotropic discriminant subgroups only after integrality and evenness are
   established.
4. **Kneser reglue.** Across a common index-`p` sublattice, remove the old
   isotropic coset and add a new one.  This is the atomic finite operation
   capable of removing roots while staying in a genus.

For the usual clean `p`-neighbor construction, take `p` odd and prime to the
determinant; `p=2` and primes dividing the discriminant require separate local
care.  Neighbor connectivity also uses the hypotheses and allowed-prime set
in the relevant spinor-genus theorem, not an unconditional fixed-`p`
statement for every genus; see Theorem 3.18 of
[Voight](https://arxiv.org/abs/2308.11566).  The statistics in
[Chenevier](https://arxiv.org/abs/2104.06846) give complementary computational
context.  More pointedly, [Chenevier--Lannes](https://arxiv.org/abs/1409.7616)
prove that a rooted Niemeier lattice is a `p`-neighbor of the Leech lattice
exactly when the prime `p` is at least its root-system Coxeter number.  This is
a rigorous precedent for a root-annihilating neighbor calculus.

### Acceptance gate

For a proposed reglue `H_old -> H_new`, compute rather than predict:

```text
H_new isotropic,
A_new = H_new^perp/H_new,
determinant and signature,
primitive closure and saturation index,
minimum of every new glue coset,
complete norm-two root set.
```

An edge is root-annihilating only when old roots lie in removed cosets or leave
the new frame, and no new coset has norm two.  Minimum at least four is the
terminal rootless gate for the positive rank-17 frame.

For repeated reglues over a fixed core, this gate can be inverted exactly.
Precompute the low-norm discriminant-coset theta coefficients of `K` and each
rank-two bridge `C`.  The norm-two count in a glue class `(a,b)` is their
convolution

```text
rho_KC(a,b)=sum_nu theta_K(a,nu)*theta_C(b,2-nu).
```

An isotropic graph subgroup produces a rootless overlattice exactly when it
is contained in the zero support of `rho_KC`.  Thus the enumerator may list
zero-support graph subgroups first and construct only their children.  This
is an exact acceptance calculus, not a statistical predictor; its practical
value depends on reusing the core table across enough bridge and graph
choices.

The implementation certificate covers the complete fourteen-class terminal
binary-bridge universe for H3, Q80, NS0024, and Golay-720.  From the four
cyclic bridge determinants it independently enumerates the fourteen reduced
positive even binary bridges, computes their tables and four core theta
tables, derives all 28 oriented graph multipliers from finite-form isotropy,
and tests them without constructing a rank-17 child.  It reproduces every
stored signed root count and accepts exactly the five rootless bridge classes:
[`elkies-k3-integral-rank-transfer-theta-convolution-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-theta-convolution-v1.json).

### Completed bridge computation

The bridge-core computation is complete for H3, Q80, NS0024, and Golay-720.
For each of their 42 selected edges the verifier computes the common rank-15
lattice, both rank-two bridges, Smith-coordinate glue generators, and the
complete root-transfer partition.  It verifies

```text
Phi(W_old) intersect Phi(W_new) = Phi(K)
```

and applies the finite gate: the new frame is rootless exactly when neither
the core nor any selected new glue coset contains a norm-two vector.

The reusable per-edge algorithm is:

1. embed consecutive frames through their exact markings in the common `NS`;
2. compute `K=W_old intersect W_new` directly in the common `NS` marking;
3. compute both bridge lattices and graph subgroups;
4. label every norm-two vector by surviving, removed, or newly introduced
   glue coset;
5. deduplicate successful labels under the relevant automorphism groups when
   enumerating new target-directed transitions.

No separate rational isometry or Kneser chain is required for the local
theorem: the common `NS` marking supplies the canonical carrier `K`.  A later
good-prime Kneser decomposition remains useful for connectivity, but is not
the definition of the observed hop.

### First predictive test: the split bound is not selective

The normalization makes one exact pre-classification screen available.  For a
proposed new `U'`, both the roots of the common core `K` and the roots of the
rank-two bridge `C_new` survive in the child, because

```text
K + C_new subset W_new.
```

Thus their orthogonal sum gives a mandatory lower bound on both child root
rank and signed root count.  Rejecting a candidate when this lower bound
exceeds a predeclared target budget has no false negatives for that budget and
does not enumerate roots of the glued rank-17 child.

The retrospective benchmark
[`elkies-k3-integral-rank-transfer-bridge-predictor-benchmark-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-predictor-benchmark-v1.json)
replays five exact H3 first-hit histories.  Four use complete MW-quotient
shells; the terminal q6 case uses its pinned capped stream window.  The result
is negative:

```text
raw child classifications through the five first hits:  2,892
pass the core-only bound:                              2,892
pass the K + C_new bound:                             2,714
classifications rejected:                               178  (6.15%)
projected classification speedup:                    1.066x
```

The per-edge retained counts are `52/52`, `111/114`, `496/498`, `808/981`,
and `1247/1247`; in particular the rootless terminal stream receives no
benefit.  This is not substantially better than raw q-neighbor enumeration,
so the evidence does **not** support a new construction algorithm.

### Exploratory terminal reglue test: bridge minimum is enriched

The same benchmark also asks the narrower question suggested by the completed
normalization.  For each observed rootless terminal core, it fixes that core
and its observed prime bridge determinant, then exhausts every
Minkowski-reduced positive even binary bridge class and every compatible
oriented cyclic graph label in the target discriminant genus.  Full root
enumeration labels the outcome but is not part of the tested predictor.

```text
admissible binary bridge classes:                    14
rootless classes without screening:                   5  (35.7%)
classes of maximum bridge minimum:                     5
rootless classes retained:                             4  (80.0%)
rootless recall:                                       4/5
precision enrichment:                               2.24x
projected full-classification speedup:               2.80x
```

The rule selects the unique rootless class in H3, Q80, and Golay-720.  For
NS0024 the maximum minimum is shared by two classes, one rootless and one
rootful, and a second rootless class of smaller minimum is missed.  Thus the
rank-two bridge data do predict successful reglues substantially better than
unfiltered enumeration **inside these four fixed successful cores**.

This is exploratory enrichment, not yet a new construction algorithm.  The
cores and determinants were selected from the observed successful edges, so
the test does not measure whether bridge data can choose a promising core
from a raw q-neighbor shell.  The broader H3 test also shows why that gap
matters: most spoiling roots occur in nonzero graph-glue cosets.  The next
test must therefore declare a core-generation rule and a cheap coset-minimum
score before examining an untouched shell, then compare both success rate and
runtime with direct child-root enumeration.  The retained Q80 score tables
contain only already rank-growing children and do not supply an unbiased
negative corpus for that prospective test.

### Held-out prospective test: bridge minimum is not selective

<!-- status-consumer: EC-K3-E6-DET78-PROSPECTIVE-BRIDGE-NEGATIVE d23a0abd146c2ed9 -->

The predeclared protocol was run without successful-corridor inputs on the
complete determinant-78 E6 old-degree-two shell.  Candidate generation uses
only the source `2E6+A1/MW4` frame and its root Weyl group.  For each child,
the common rank-15 core, new rank-two bridge, and every nonzero graph-glue
coset minimum are computed before root enumeration or truth-catalogue lookup.

```text
dominant source-Weyl classes:                         280
primitive candidates:                                277
distinct mass-closed J2 classes reached:              31
bridge glue-coset minimum distribution:           {2:277}
child root-rank distribution: {12:18,13:168,14:71,15:20}
```

Thus the least nonzero glue-coset minimum is constant and supplies no ranking
signal on this genuinely held-out shell.  The recorded scoring time is 51.01
seconds, versus 0.48 seconds for direct norm-two enumeration and root-rank
classification: the score is about `106x` slower on this workload.  The
separate 145.06-second J2 isometry lookup is truth-set evaluation and is
excluded from that comparison; all timings are non-certificate metadata.  The
E6 genus is globally rootful, so this is a negative control, not a positive
recall test.  It rules out the scalar minimum score here while leaving the
exact decorated profile `rho_K` as the necessary next object.
See
[`elkies-k3-e6-rank4-det78-prospective-bridge-predictor-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-det78-prospective-bridge-predictor-v1.json).

## Recent evidence and failure map

The last three commits explain why the character/glue layer is now visible.

- The deep-cover and relative 2-Selmer work found no exceptional directions
  beyond the tested degree-two atlas in its declared bounded universes; those
  are negative experiments, not nonexistence theorems.
- The E6+A1 dissection fixed a new determinant-36 Picard-19 lattice and made
  the determinant obstruction operational: determinant at most 28 cannot
  support a rootless rank-17 frame.
- The E6 ordered-incidence closeout separated geometric `2+2` from arithmetic
  descent `1+1`, while the rank-28 genus-one bisection work produced exact
  height-16 anti-invariant directions but not full saturation.
- The subsequent high-throughput genus-one search found one simultaneous
  split at `t=1/25`: the norm-twelve class `0x103b2` contributes an exact
  quotient direction beyond specialized generic MW17.  Its cover-level
  character tuple is now exact: `R17(2)+<16>` acquires one order-two graph
  class and becomes a rootless visible rank-18 lattice.  Isolating the large
  eclib prime checks makes the full specialization saturation fit in memory;
  the displayed rank-18 subgroup is primitive, though no rank upper bound is
  known.

Other preserved failures constrain the theory:

- determinant, parity, minimum, and rootlessness do not determine the local
  genus;
- a pure overlattice cannot delete roots;
- an unchecked half-section or torsion class can change the NS determinant;
- a rational unordered quotient need not descend individual eigensections;
- the 39,120 rational R17 bisections have distinct squareclasses, so two
  anti-invariant gains do not occur on one quadratic cover by collision;
- bounded neighbor shells, timeouts, and sampled `p`-neighbor walks are not
  global obstructions.

## Reproduction

Replay the bridge-core theorem on all normalized marked edges:

```bash
sage -python \
  elkies-k3/scripts/certify_integral_rank_transfer_bridge_reglue.sage --check
```

Replay the retrospective bridge-split predictor benchmark:

```bash
sage -python \
  elkies-k3/scripts/benchmark_integral_rank_transfer_bridge_predictor.sage --check
```

Replay the blind determinant-78 prospective negative control:

```bash
sage -python \
  elkies-k3/scripts/benchmark_e6_det78_prospective_bridge_predictor.sage --check
```

Replay the exact inverse theta-convolution enumerator:

```bash
sage -python \
  elkies-k3/scripts/certify_integral_rank_transfer_theta_convolution.sage --check
```

Replay the first exact target-root-system control:

```bash
sage -python \
  elkies-k3/scripts/certify_ns0024_inverse_ade_mutation.sage --check
```

Byte-check the norm-12 `0x103b2` character tuple while reusing its pinned full
specialization-saturation record:

```bash
sage -python \
  elkies-k3/scripts/certify_r17_norm12_103b2_mw_glue.sage \
  --skip-specialization-saturation --check
```

Replay the exhaustive involution-graph classification:

```bash
sage -python \
  elkies-k3/scripts/certify_integral_character_glue_calculus.sage --check
```

Replay the genus-wide determinant-78 E6 obstruction:

```bash
sage -python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage \
  --rootless-obstruction --check
```

Then generate and byte-check the deterministic equation-free census:

```bash
sage -python \
  elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage
sage -python \
  elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage --check
```

The builder recomputes Smith invariants and normalized finite quadratic forms
from the pinned integral Gram matrices.  It also checks the route rank budget,
the NS0024 and Golay edge sequences against their exact route certificates,
and all input hashes.  It intentionally contains no Weierstrass coefficients.

## Literature boundary

The underlying overlattice, embedding, elliptic-surface, and neighbor tools are
standard: [Nikulin](https://www.mathnet.ru/eng/im1677),
[Nishiyama](https://www.jstage.jst.go.jp/article/math1924/22/2/22_2_293/_article),
[Shioda](https://rikkyo.repo.nii.ac.jp/records/10027), and
[Schuett--Shioda](https://arxiv.org/abs/0907.0298).  Rootless coinvariant
lattices for finite symplectic K3 actions provide a separate precedent for
symmetry **selecting a rootless coinvariant sublattice**; see
[Nikulin](https://arxiv.org/abs/1109.2879) and
[Garbagnati--Sarti](https://arxiv.org/abs/math/0603742).

The equation-free `O(NS)`/J2 calculus is complete in the following precise
sense: a mass-closing neighbour enumeration terminates with the whole frame
genus, and the decorated glue-coset root profile decides every edge exactly.
The bare tuple `(A,q,H)` is not enough--the Leech and rooted Niemeier lattices
already disprove that--so coset minima or root counts are essential data.
What remains open is a uniform small-prime or short-path bound, J1
classification modulo surface automorphisms, and a general lift of a lattice
path to explicit elliptic equations over a prescribed ground field.

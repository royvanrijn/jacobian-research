# Integral rank transfer: an equation-free glue calculus (2026-09-02)

## Outcome

The proposed experiment has a sharp first answer.

There are three recurring integral operations in the successful examples:

1. replace the primitive hyperbolic plane `U` inside a fixed `NS` and take
   the new saturated orthogonal frame;
2. glue an auxiliary lattice to a frame by the graph of a discriminant-form
   anti-isometry, then take the complementary lattice;
3. split into rational finite-group isotypic components and recover the integral
   lattice by `|G|`-primary saturation glue.

The first operation accounts for every edge in the H3, Q80, NS0024, and
Golay-720 corridors.  It is **not** a mutation of the finite quadratic form:
because `U` is unimodular, the frame discriminant form is fixed along the
whole corridor.  The finite form is an admissibility and genus invariant, not
the object which records which roots disappeared.

There is now an exact local replacement theorem behind these pivots.  For two
frames `W,W'` in the same `NS`, their common bridge core

```text
K=W intersect W'=(U+U')^perp(-1)
```

contains precisely the roots that survive the click.  Each frame is a finite
graph-glue extension of `K+C_i`, where `C_i=K^perp` inside that frame.  On all
13 NS0024 edges and all six Golay-720 edges, `K` has rank 15, both `C_i` have
rank two, and the two glue groups are cyclic of the same order and project
onto the full bridge discriminant.  Thus all 19 checked clicks are exact
rank-two cyclic bridge replacements.  The certificate is
[`elkies-k3-integral-rank-transfer-bridge-reglue-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-reglue-v1.json).

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
quotient is `(Z/2)^2`.

The requested equation-free tuples and all 42 selected same-`NS` corridor
edges are stored in
[`elkies-k3-integral-rank-transfer-glue-census-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-glue-census-v1.json).
This note is a theory-extraction experiment.  It does not change
`MATH_STATUS.json` and does not assert the conjectural neighbor calculus below.

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
| E6 `2+2` | `L=12(MW_+ direct_sum MW_-)`, rank 4, determinant 1,437,696; exact `A_L,q_L` in the census | `C2`; scaled blocks have ranks `2+2`; unscaled determinants `4/3`, `52` | two half-sums give `(Z/2)^2`, index four; this is literal even glue only in the scaled model | ambient `NS` has determinant 78 and roots `2E6+A1`; ordered incidence has `2+2`, rational quotient `1+1` |
| R17 genus-one bisection cover | `L=R17(2)+<16>`, the orthogonal character sublattice, rank 18 | `C2`; invariant rank 17 and one certified anti-invariant height-16 line | `2R=tau+T` gives a certified `Z/2` half-sum extension of index two; only further saturation is unknown | rank at least 18, not an exact total-rank or full-MW claim |
| paired rational R17 norm-ten bisections | `L=R17(4)+<24>+<24>`, the orthogonal character sublattice, rank 19 | `V4`; two distinct nontrivial characters give the height-24 lines | the two independent half-sums give `(Z/2)^2`, index four; further saturation and the product character are unknown | all 765,167,640 pairs are geometrically genus-one `V4` covers; rational points are separate |
| norm-12 / norm-8 simultaneous split | specialized section lattice at `t=1/25`, certified rank at least 18; Gram and finite form not determined | no integral deck-action lattice extracted | norm-12 class `0x103b2` escapes specialized generic MW17; no glue mutation yet identified | exact specialization quotient direction, not a generic rank-transfer theorem |

The H3 entry is deliberately not written as three arrows.  `J2` means an
integral frame-isometry class.  The historical H3 corridor, the Niemeier J2
classification, and identification of the published R17 endpoint are three
different proof layers.

Here “Golay” means the binary-octad determinant-720 construction and its
rootful source corridor.  The ternary-Golay `N(12A2)` foundry backend is a
large census (151 frame classes after local deduplication), not one of these
named transition chains.  It should be expanded as a family of tuples rather
than misrepresented by one aggregate tuple.

## What happened at the same-`NS` clicks

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

For the 19 normalized marked edges, the finer operation is

```text
W  = overlattice(K+C_old, H_old)
  -> overlattice(K+C_new, H_new) = W',
```

where `rank(K)=15`, both bridges have rank two, and `H_old,H_new` are cyclic
maximal bridge graphs of equal order.  The observed orders are

```text
23, 31, 47, 63, 119, 143, 191, 215, 303, 359, 1231.
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

## What happened at the literal glue clicks

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
S=(NS,U,W,R,A_W,q_W,H_root,minimum,roots,saturation indices).
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

### Completed bridge computation and remaining normalization

The bridge-core computation is complete for NS0024 and Golay-720.  For each of
their 19 edges the verifier computes the common rank-15 lattice, both rank-two
bridges, Smith-coordinate glue generators, and the complete root-transfer
partition.  It verifies

```text
Phi(W_old) intersect Phi(W_new) = Phi(K)
```

and applies the finite gate: the new frame is rootless exactly when neither
the core nor any selected new glue coset contains a norm-two vector.

The remaining exact normalization is to put the 13 H3 and ten Q80 edges into
the same marked-transport interface.  For each remaining pair:

1. choose abstract integral representatives of consecutive frames in their
   common genus and record a rational isometry embedding both into one common
   rational quadratic space;
2. compute `K=W_old intersect W_new` directly in the common `NS` marking;
3. compute both bridge lattices and graph subgroups;
4. label every norm-two vector by surviving, removed, or newly introduced
   glue coset;
5. deduplicate successful labels under the relevant automorphism groups.

No separate rational isometry or Kneser chain is required for the local
theorem: the common `NS` marking supplies the canonical carrier `K`.  A later
good-prime Kneser decomposition remains useful for connectivity, but is not
the definition of the observed click.

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
  quotient direction beyond specialized generic MW17.  Its integral Gram and
  glue mutation have not yet been extracted, so the census records the
  success with explicit unknown tuple fields rather than guessing them.

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

Generate the deterministic equation-free census:

```bash
sage -python \
  elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage
```

Verify the pinned artifact byte for byte:

```bash
sage -python \
  elkies-k3/scripts/build_integral_rank_transfer_glue_census.sage --check
```

Replay the bridge-core theorem on all normalized marked edges:

```bash
sage -python \
  elkies-k3/scripts/certify_integral_rank_transfer_bridge_reglue.sage --check
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

What is new and still open here is not the finite-form formalism.  It is a
local, computable criterion which predicts from a reglue label that roots will
be annihilated, together with a theorem that the relevant marked elliptic K3
frames are connected by a controlled set of such moves.

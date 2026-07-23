# Relation-complex deformation theory

> **Status and scope.** This note is a reusable language and research
> protocol, not a theorem asserting that the four theories below are
> equivalent.  Its collision and degree-thirty Ritt examples are established
> repository calculations.  The derived packaging, cycle-neutrality
> criterion, and proposed transfers to GMC and Kuranishi theory are a
> programme.  Any application must still prove that its chosen combinatorial
> complex and local algebras represent the original moduli problem.

Several parts of the repository have the same two-layer shape:

1. a combinatorial incidence object records which presentations can agree,
   collide, or be related; and
2. local algebra records the scheme structure of those interactions.

The recurring examples are:

* zero-sum weight relations in the two-real Gaussian moment problem;
* allocation partitions and primitive mergers in omitted-component
  intersections;
* decomposition words, Ritt moves, and Coxeter cells in Hessian--Ritt
  intersections; and
* tangent directions, obstruction directions, and binary-form pencils in
  Kuranishi calculations.

The purpose of the framework is operational.  It should say which local
algebras to compute, which comparisons are meaningful, when cycles are
irrelevant to length, and where higher coherence can create or detect
nonreduced structure.

## 1. Two objects, not one hypergraph

A bare hypergraph is insufficient.  There are two related combinatorial
objects.

### 1.1 The constraint hypergraph

Let `V` be a finite set of presentations, sheets, weights, or deformation
coordinates.  A hyperedge records one primitive simultaneous relation among
vertices.  Examples include a zero-sum weight circuit, a collision block at
one root, or a collection of decomposition cuts imposed at once.

If every constraint is represented by a map of affine schemes

\[
 Z_e\longrightarrow M
\]

to a common ambient moduli space, the algebra of a face `sigma` is the
derived intersection

\[
 A_\sigma
 =
 \mathop{\bigotimes\nolimits^{\mathbf L}}_
        {e\in\sigma,\mathcal O(M)}
 \mathcal O(Z_e).
\tag{1.1}
\]

This formula is appropriate for simultaneous constraints in one ambient
space.  Its degree-zero homology is the ordinary scheme-theoretic
intersection, while its negative homology records excess Tor.

### 1.2 The coherence complex

Some relations are correspondences between different presentation spaces,
not subschemes of one fixed ambient space.  Let `C` be a graph or cell
complex with

* a presentation space `X_v` at each vertex;
* a derived correspondence `X_e` at each oriented edge; and
* a comparison between composite correspondences on the two sides of each
  two-cell.

Formally, this is a correspondence-valued weak two-functor

\[
 \mathcal X:C\longrightarrow\operatorname{Corr}(\mathrm{dAff}_k).
\tag{1.2}
\]

A path means an iterated derived fiber product of its edge
correspondences.  A two-cell does not merely identify its boundary paths
combinatorially: it supplies, or asks for, a comparison of their full path
schemes.

The constraint hypergraph and coherence complex may be derived from the
same presentations, but they answer different questions.  Conflating them
is exactly what makes a graph cycle appear to predict nilpotent length when
it need not do so.

## 2. Local decorations

The useful object is a pair

\[
 (C,\mathcal A),
\tag{2.1}
\]

where \(C\) is the relevant incidence complex and \(\mathcal A\) is a
constructible diagram of augmented derived commutative algebras, or
equivalently a diagram of derived affine correspondences.  The decoration
should retain:

* the base algebra and its specialization strata;
* the full local defining ideal or a quasi-isomorphic Koszul cdga;
* face and specialization maps;
* augmentations to a reduced support or chosen geometric point;
* filtrations such as the nilradical, conductor, or parameter-adic
  filtration; and
* comparison maps attached to two-cells.

The phrase **relation hypergraph deformation algebra** may be used for a
single simultaneous-intersection algebra (1.1).  The broader phrase
**relation-complex deformation theory** refers to the full decorated
diagram (2.1).  The latter is the object needed for Ritt coherence.

## 3. Path and cell defects

Let \(\gamma_0\) and \(\gamma_1\) be the two boundary paths of a two-cell.
Write \(A_{\gamma_i}\) for their derived path algebras.  Often both paths
have a
chosen common reduced image or normalization `B`, giving augmentations

\[
 A_{\gamma_0}\longrightarrow B
 \longleftarrow A_{\gamma_1}.
\tag{3.1}
\]

The individual path defects are

\[
 \Delta_{\gamma_i}
 =
 \operatorname{fib}(A_{\gamma_i}\longrightarrow B).
\tag{3.2}
\]

If the geometry supplies a comparison
\(\phi_c:A_{\gamma_0}\to A_{\gamma_1}\), define the oriented **cell
defect**

\[
 \Delta_c=\operatorname{fib}(\phi_c).
\tag{3.3}
\]

Here `fib` is taken on the underlying derived modules or complexes; the
defect need not itself be a unital algebra.  In an underived surjective
comparison it is the ordinary kernel module.  Without a canonical
comparison, retain the span (3.1) itself.  Choosing an arbitrary map between
the two path algebras would destroy invariance.

For every path or cell defect, record the smallest useful fingerprint:

\[
 \left(
 \operatorname{Supp}\Delta,\,
 \operatorname{Ann}\Delta,\,
 H_*(\Delta),\,
 L_{A/k},\,
 P_k^A(T),\,
 \operatorname{Hilb}(A),\,
 \operatorname{Soc}(A)
 \right).
\tag{3.4}
\]

The support and annihilator locate the failure.  The cotangent complex
separates tangent excess from relation excess.  Tor and the residue-field
Poincare series detect derived complexity.  Hilbert, filtration, and socle
data distinguish local algebras of the same length.

## 4. Three meanings of a cycle

The following must be kept separate:

\[
 H_1(C;k),\qquad
 \operatorname{Tor}^{R}_*(A_1,\ldots,A_s),\qquad
 m_n\ \text{or Massey products}.
\tag{4.1}
\]

The first is combinatorial path homology.  The second measures
nontransversality of algebraic constraints.  The third measures higher
multiplication or coherence in a dg or \(A_\infty\) model.  None implies
either of the others without additional hypotheses.

If a global derived algebra can be expressed as derived sections of a local
diagram, one expects a descent or hypercohomology spectral sequence of the
form

\[
 E_2^{p,q}
 =
 H^p\!\left(C;\mathcal H^q(\mathcal A)\right)
 \Longrightarrow
 H^{p+q}\!\left(R\Gamma(C,\mathcal A)\right).
\tag{4.2}
\]

Equation (4.2) is a target formalism, not yet a theorem for every example in
this repository.  It explains what a useful theorem would need to separate:
cellular cohomology, local Tor, and differentials or higher operations that
couple them.

## 5. The four dictionaries

### 5.1 Gaussian moment support

For a rotational weight support `S`, the canonical combinatorial object is
the nonnegative zero-sum semigroup

\[
 \Lambda_S
 =
 \left\{
 \alpha\in\mathbb N^S:
 \sum_{w\in S}\alpha_w w=0
 \right\}.
\tag{5.1}
\]

Primitive support-minimal elements give the circuit-incidence hypergraph.
The semigroup, rather than only the unoriented circuit matroid, is needed:
positive multiplicities and the Gaussian factorial functional affect the
moment equations.

For a moment cutoff `N`, the first derived algebra to compute is the Koszul
cdga

\[
 K_S^{(N)}
 =
 K(M_1,\ldots,M_N;
   k[\text{coefficient coordinates}]),
\tag{5.2}
\]

after the same localization, centering, and saturation used in the
underived chart.  Filter (5.2) by elements or factorizations of
`Lambda_S`.  This retains information discarded by the binary output
“the reduced Groebner basis is `[1]`.”

For the three-level support `{-1,0,1}`, constant-term extraction packages
the local interaction as the Bessel--factorial block

\[
 \mathcal L\!\left[
 e^{tC(U)}I_0(2t\sqrt{D(U)})
 \right].
\tag{5.3}
\]

The first comparison experiment is between star supports, whose canonical
circuit graph is \(K_{1,q}\), and the first cyclic support \(K_{2,2}\).

An Orlik--Solomon algebra may describe a squarefree linearized shadow, but
it generally loses the multiplicities in (5.1) and the functional in
(5.3).  Affine semigroup, toric face, or arithmetic-oriented-matroid
language is therefore safer.

### 5.2 Collision and allocation geometry

Vertices are normalization sheets.  A collision root gives a hyperedge
partitioning the sheets by their allocations.  The local algebra at an
active root is a transfer block `Z_h`; the primitive case is

\[
 Z_1=k[\epsilon]/(\epsilon^2).
\tag{5.4}
\]

At distinct collision roots, formal Hensel factorization gives the proved
product decomposition

\[
 A_{\mathrm{tr}}
 \simeq
 \widehat{\bigotimes}_{\rho\ \mathrm{active}}Z_{\rho}.
\tag{5.5}
\]

The branch-merger graph instead records paths between allocation words.
Its cycles are not the tensor factors in (5.5).  This is why a triangle of
sheets over three active roots has one graph cycle but transverse algebra

\[
 k[\epsilon_1,\epsilon_2,\epsilon_3]/
   (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2)
\tag{5.6}
\]

of length eight.

When Hensel clusters merge, multiplication can change while length remains
constant, as in

\[
 k[X,Y]/(X^3,XY,Y^2),
\tag{5.7}
\]

or even length can jump through a spectator collision.  These are failures
of local factorization, not consequences of an abstract graph cycle.

### 5.3 Ritt intersections

Vertices are normalized complete decomposition words.  Edges are
elementary power or Chebyshev Ritt correspondences.  Commuting squares and
braid hexagons are genuine coherence cells.

In the degree-thirty braid, the two half-boundary path ideals on a fixed
endpoint chart satisfy

\[
 I_L\subsetneq K=I_R.
\tag{5.8}
\]

Thus the quotient map

\[
 S/I_L\longrightarrow S/K
\tag{5.9}
\]

has defect `K/I_L`.  It is supported on the monomial divisor `z=0`, and
adding the missing composite cut kills it.  The three labelled sectors have
different nilpotence indices, annihilator filtrations, and Artin slice
algebras.  Therefore the unlabelled Coxeter hexagon and the common
normalization do not determine the scheme-level two-cell.

This is the primary model for (3.1)--(3.4).  A higher-degree theory should
test whether its braid defects are universal local blocks under smooth base
change, tensor products with separated spectator blocks, or controlled
extensions when spectators collide.

### 5.4 Kuranishi geometry

Let

\[
 \kappa:T^1\longrightarrow T^2
\tag{5.10}
\]

be a Kuranishi map.  Vertices may represent tangent coordinates and
obstruction coordinates, while a hyperedge records the support of one
component of the quadratic obstruction.  Its local algebra is the completed
Kuranishi algebra

\[
 A_\kappa
 =
 k[[{T^1}^{\!*}]]/
   (\text{components of }\kappa).
\tag{5.11}
\]

The support hypergraph alone is not enough.  Ranks, discriminants, and
common isotropic loci of the binary-form pencils determine the quadratic
branches.  Cubic and higher continuation obstructions are naturally the
higher brackets of a controlling \(L_\infty\)-algebra.  This is the example
in which Massey products or higher operations are structurally native,
rather than merely suggested by a combinatorial cycle.

## 6. Reusable calculation protocol

Every application should answer the following questions in order.

### Step A: specify the moduli problem

1. What is the ambient ring or presentation space?
2. Which data are intrinsic, and which depend on a chart or marking?
3. Is the problem a simultaneous intersection or a composition of
   correspondences?
4. What saturation, completion, or localization is required?

### Step B: build the correct incidence object

1. List vertices and primitive relations.
2. Separate constraint hyperedges from path edges.
3. Add commuting, braid, or other coherence cells.
4. Label cells by degree, multiplicity, allocation, or weight data that can
   alter the local algebra.

### Step C: compute local algebras

For every relevant stratum or cell:

1. compute the ordinary ideal and its reduction;
2. construct a Koszul or other derived presentation;
3. compute tangent and obstruction dimensions from the cotangent complex;
4. after eliminating smooth variables, find a finite homogeneous initial
   envelope and record its top degree \(N\) and length \(L\);
5. apply the
   [Kuranishi nilpotence cutoff theorem](DEFECT_SYMBOL_APOLARITY.md#2-kuranishi-nilpotence-cutoff-theorem)
   immediately when a comparison defect is already known to have order
   \(d>N\);
6. compute the full length, Hilbert vector, socle, and multiplication table
   only when the exact Artin algebra, rather than vanishing of one defect,
   is needed;
7. compute annihilator and support of every comparison defect; and
8. retain parameter-adic, conductor, or nilradical filtrations.

Equal lengths or equal Tor ranks are not sufficient grounds for identifying
two blocks.  Conversely, exact identification is unnecessary when the
finite envelope already lies below the certified defect order.

### Step D: test factorization before invoking cycles

1. Do completed local factors split by Hensel uniqueness?
2. Are the relevant maps flat or Tor-independent?
3. Is a spectator relation contained in a distinct factor?
4. Does specialization merge factors?

Only after these tests should graph or cell homology be compared with the
derived algebra.

### Step E: compare paths and fill cells

For every two-cell:

1. compute both path schemes;
2. compare their sets, reductions, normalizations, and full scheme
   structures separately;
3. construct the canonical comparison span or map;
4. compute the defect fingerprint (3.4); and
5. check what changes when the entire boundary constraint is imposed.

### Step F: state the conclusion at the correct strength

Distinguish:

* combinatorial connectivity;
* equality of point sets;
* equality of reduced schemes;
* equality after normalization;
* equality of ordinary schemes;
* derived equivalence; and
* coherent equivalence across all higher cells.

## 7. Calculation record

New computations can use the following compact record.  Fields that have
not been proved should be marked `unknown`, not inferred from nearby
examples.

```text
name:
ambient_ring:
base_stratum:
completion_or_saturation:

incidence_object:
  vertices:
  edges_or_hyperedges:
  higher_cells:
  labels:

local_algebra:
  presentation:
  reduced_algebra:
  normalization:
  hilbert_vector:
  socle_dimension:
  nilpotence_index:
  cotangent_homology:
  poincare_series:

comparison:
  left_path:
  right_path:
  canonical_map_or_span:
  defect_module:
  defect_support:
  defect_annihilator:
  boundary_filling_effect:

factorization:
  hensel_factors:
  tor_independence:
  spectator_behavior:

status:
  proved:
  conjectural:
  executable_certificate:
```

## 8. Three working principles

The following are the useful theorem targets.

### Hensel factorization principle

Separated relation blocks should tensor in the derived category.  In the
collision setting, the underived completed factorization is proved.  A
derived version should identify hypotheses guaranteeing that no additional
cross-factor Tor appears.

### Cycle-neutrality principle

If local blocks split into independent completed factors and all gluing
maps are derived-flat and homotopy-Cartesian, combinatorial path cycles
should not alter the degree-zero local length.  They may still contribute
to global cellular cohomology.  This principle must not be applied when
specialization merges local factors.

### Braid-defect locality principle

A Ritt two-cell defect should be supported on the locus where the
correspondence comparison loses transversality, flatness, or effective
descent.  Higher-degree defects should be functorial under smooth base
change and separated spectators, with new blocks appearing when spectators
collide.

These are conjectural principles outside the cases already computed.

## 9. Highest-value next tests

1. **Degree-forty-two braid universality.** Compare the degree-thirty local
   defects with the cut-`6` degree-forty-two spectator sector.  Determine
   whether the latter is a base change, a tensor product with a smooth
   factor, or a genuinely new extension.
2. **GMC forest versus first cycle.** At matched radial degree and moment
   cutoff, compare localized Koszul homology and multigraded Betti data for
   a star support \(K_{1,q}\) and the first circuit cycle \(K_{2,2}\).
3. **Kuranishi incidence versus higher obstruction.** On the established
   normal branches, compare the quadratic support hypergraph and pencil
   discriminant with the rank jump of the next continuation equation.
   This tests whether the higher obstruction is visible combinatorially or
   only in the \(L_\infty\) decoration.

## 10. Guardrails

This framework does not justify any of the following shortcuts:

* a graph cycle automatically creates Tor or a nilpotent;
* absence of graph cycles implies a reduced intersection;
* equal local length implies equal local algebra;
* equal normalization implies equal path scheme;
* nonzero Tor implies a nontrivial Massey product;
* an incidence hypergraph determines a Kuranishi map without its
  coefficients; or
* bounded GMC moment calculations establish an all-order statement.

The value of the framework is precisely that it records the extra local
algebra and coherence data needed to decide when such implications do or do
not hold.

## 11. Existing anchors

The primary repository inputs are:

* [support graphs for the two-real Gaussian moment
  problem](TWO_REAL_GMC_SUPPORT_GRAPH_EXPLORATION.md);
* [primitive merger hypergraphs and omitted-component intersection
  algebras](OMITTED_COMPONENT_INTERSECTION_ALGEBRA.md);
* [the Ritt move two-complex and degree-thirty
  braid](RITT_MOVE_2_COMPLEX.md);
* [general Hessian--Ritt intersections](GENERAL_HESSIAN_RITT_INTERSECTIONS.md);
* [filtered rank-two Kuranishi
  calculations](RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md); and
* [defect-symbol and apolarity
  package](DEFECT_SYMBOL_APOLARITY.md).

These documents remain authoritative for their individual theorems and
scope restrictions.  The present note supplies only the common language and
calculation protocol.

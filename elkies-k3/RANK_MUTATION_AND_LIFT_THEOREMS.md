# Rank mutation and lift theorems

This note extracts general mathematics from the Elkies--K3 calculations. It
separates statements that follow from standard K3 and elliptic-surface
theorems, conditional correctness theorems for the equation compiler, and
genuinely open navigation conjectures.

Status boundary: the proofs below are a theorem-development package, not yet
new entries in `MATH_STATUS.json`. They do not promote the active orbit42
artifact or prove that the selected route is optimal.

## 1. Setup

Let `X` be a smooth projective K3 surface over an algebraically closed field of
characteristic zero. A Jacobian elliptic fibration `pi` has fibre class `F` and
zero section `O`. The classes

```text
F, O+F
```

span a copy `U_pi` of the hyperbolic plane inside `NS(X)`. Write

```text
W_pi = orthogonal complement of U_pi in NS(X),
R_pi = root lattice from non-identity reducible-fibre components,
r_pi = rank MW(pi),
t_pi = size of MW(pi)_tors,
Reg_pi = determinant of the free MW height lattice.
```

Changing the elliptic fibration means changing the embedded copy of `U`; it
does not change `NS(X)`.

## 2. The exact rank-mutation law

### Theorem A: conservation of the divisor budget

For two Jacobian elliptic fibrations `pi_1` and `pi_2` on the same K3 surface,

```text
r_2 - r_1 = rank(R_1) - rank(R_2).
```

Equivalently,

```text
rank(R_i) + r_i = rho(X) - 2.
```

#### Proof

For each fibration, the trivial lattice has rank `2+rank(R_i)`. The
Shioda--Tate formula gives

```text
rho(X) = 2 + rank(R_i) + r_i.
```

Subtract the two formulas. Nothing about the equation, neighbour degree, or
chosen route is needed. QED.

### Corollary A1: rank cannot appear from nowhere

Along any chain of fibrations on one fixed K3,

```text
r_n = r_0 + rank(R_0) - rank(R_n).
```

Thus removing one independent fibre root creates exactly one MW rank; removing
two creates two. Conversely, a reverse neighbour stores MW rank in reducible
fibres. For a rank-19 K3, a rootless Jacobian fibration automatically has MW
rank 17.

This proves the rank changes in the H3 and Q80 lattice corridors once the
marked fibrations, Picard rank, and root systems are certified. It does not
construct explicit coordinates for the new sections.

## 3. The determinant and saturation laws

### Theorem B: determinant mutation

For a Jacobian elliptic fibration,

```text
abs(disc NS(X)) = abs(disc R_pi) * Reg_pi / t_pi^2.
```

Consequently, for two fibrations on the same surface,

```text
Reg_2 / Reg_1
  = abs(disc R_1) / abs(disc R_2) * (t_2 / t_1)^2.
```

#### Proof

The trivial lattice is `U + R_pi`. Shioda's orthogonal projection identifies
the free Mordell--Weil group with the height lattice in the rational
orthogonal complement. The primitive-closure defect of the trivial lattice is
exactly MW torsion. Taking lattice discriminants gives the first formula; the
second follows by cancelling the fixed discriminant of `NS(X)`. QED.

This is stronger than rank conservation: it predicts the determinant of the
new MW lattice before its generators are explicitly lifted.

### Lemma B1: every saturation error is a square

If a full-rank lattice `L_0` has index `n` in its saturation `L`, then

```text
abs(det L_0) = n^2 * abs(det L).
```

#### Proof

Choose bases related by an integer matrix `A` of determinant `n`. Their Gram
matrices satisfy `G_0=A^T G A`, so `det(G_0)=det(A)^2 det(G)`. QED.

Hence a regulator mismatch by `81=9^2` is an index-9 warning, exactly as in the
rank-3 Elkies--K3 calculation. A non-square mismatch cannot be repaired by
finite-index saturation alone; one of the roots, torsion, heights, Picard rank,
or NS discriminant is wrong.

## 4. When an isotropic lattice vector is an elliptic fibration

### Theorem C: fibration from a primitive nef isotropic class

Let `D` be a nonzero primitive class in `NS(X)` with

```text
D^2 = 0 and D nef.
```

Then `|D|` is a base-point-free genus-one pencil and `h0(X,O(D))=2`. If there
is an effective divisor `O` with

```text
O^2=-2 and O.D=1,
```

then the pencil is Jacobian: some irreducible component of `O` has degree one
over the base and is a section. If `O` is irreducible, it is that section.

#### Proof

Riemann--Roch on a K3 gives `chi(O(D))=2+D^2/2=2`. The standard K3
base-point-freeness theorem for primitive nef isotropic classes says that
`|D|` is a genus-one pencil; primitivity rules out a multiple pencil. Since
the total horizontal degree of `O` is one, exactly one component maps with
degree one to the base and is therefore a section. QED.

### Proposition C1: Weyl reduction is a proof step, not a heuristic

The standard K3 Weyl-chamber theorem moves either sign of an isotropic class
in the positive cone to the nef cone. Algorithmically, reflection across an
effective `(-2)` wall having negative pairing strictly lowers intersection
with a fixed ample class after retaining the effective sign. Integral ample
degree makes this descent terminate. The reflection record is the exact
chamber/fixed-component correction.

The important qualification is global: nonnegative intersection with a
supplied finite list proves only `nef_on_declared_walls`. Global nefness needs
either all effective `(-2)` walls or an independent effective-cone theorem.
This is precisely the boundary already recorded by the exact-neighbour engine.

### Proposition C2: finite horizontal-wall test at fixed old-fibre degree

Write a marked Neron--Severi lattice as `U + L(-1)`, with positive-definite
Gram matrix `M` on `L`, and let

```text
D=(a,b,w),  D^2=0,  b>0.
```

For an old-horizontal `(-2)` class `C=(k,m,x)`, set `y=b*x-m*w`.  Then

```text
y.M.y = 2*b*m*(D.C) + 2*b^2.
```

Consequently, after vertical walls have been checked, every horizontal wall
with `D.C<0` occurs for some `1<=m<=b` among the finite vectors

```text
y == -m*w (mod b*L),   y.M.y < 2*b^2,
```

subject to `x=(y+m*w)/b` being integral and
`k=(x.M.x-2)/(2*m)` being integral.  Thus nefness at any fixed horizontal
degree `b` has an exact finite lattice test; it is not restricted to the
section-only closest-vector gate.

#### Proof

The root equation gives `x.M.x=2*k*m+2`, while isotropy gives
`w.M.w=2*a*b`. Expanding `(b*x-m*w).M.(b*x-m*w)` yields the displayed
identity. If an effective irreducible curve has negative intersection with
the effective class `D`, it is a fixed component, so `D-C` is effective and
`0<=F.(D-C)=b-m`; hence `0<=m<=b`. The case `m=0` is vertical. For `m>0`,
negative intersection forces the stated strict norm bound, and the congruence
and divisibility conditions reconstruct exactly the possible root classes.
Finiteness follows from positive-definiteness of `M`. QED.

#### Corollary C2.1: equation cost is scored after physical Weyl reduction

If the walls used in Proposition C1 are components of the old reducible
fibres, each reflection preserves the old-fibre degree `b` and the horizontal
class modulo the trivial lattice, while preserving isotropy and primitivity.
It need not preserve the first `U` coordinate `a`, the presentation value
`q=a*b`, the vertical layering, or a resolved-RR cost estimate.  Therefore a
root-dominant class in an abstract adapted basis is not a compiler-cost object
until it has been reduced against the physical affine cycles and has passed
the finite horizontal-wall test of Proposition C2.

This is an exact consequence of the reflection formula and the fixed-component
argument in C1--C2, not a claim that the cost always decreases.  In the H3
component-9-zero `2A5` marking, the stored q104 representative has negative
physical degrees.  Sixty-one recorded reflections produce a q10 degree-two
representative with the same horizontal quotient, `P.O=5`, three vertical
layers, and expected RR ambient 15.  Its complete physical, all-section, and
finite-horizontal-wall audit is
[`../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-certificate.json`](../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-certificate.json).

For reproducible enumeration one may use the positive-definite augmented
form

```text
Q_m(x,z)=(b*x-m*w*z).M.(b*x-m*w*z)+z^2
```

through norm `2*b^2-1` and retain `z=+/-1`. This is the gate used by
`probe_h92_pinned_r17_targeted_shell_cvp.sage` for the degree-three and
degree-four reverse searches. It certifies each retained candidate; the
target-directed ray/scale sampling that proposes candidates remains bounded.

### Proposition C3: certified neighbour loops can change the zero cheaply

Let `pi_0` and `pi_1` be marked Jacobian fibrations on the same K3 surface.
Suppose exact neighbour transports give

```text
pi_0 --D--> pi_1 --F_0--> pi_0',
```

where the second fibre class is the original fibre ray `F_0`, but the marked
section of `pi_0'` differs from that of `pi_0`.  If a curve already explicit
on `pi_0` has degree one over `D`, it may be used as the zero of `pi_1` by the
unimodular basis

```text
D, S+D, orthogonal complement of <D,S+D>.
```

If both neighbour classes pass Proposition C2 and the component walls, this
is an exact zero-changing loop, not merely an ADE recurrence.  A following
fibre can therefore have much smaller horizontal degree, pole order, or RR
dimension even though the loop temporarily revisits the same fibre ray.

#### Proof

The displayed classes have Gram `U` because `D^2=0`, `S^2=-2`, and `S.D=1`.
Primitivity of `D` and the section pairing make this `U` primitive, so its
orthogonal complement and any determinant-one completion give the full NS
transport in both directions.  The return identifies the same primitive ray
`F_0`; a different second row changes only its Jacobian marking.  Proposition
C2 plus the vertical/component audit proves that the two fibre classes define
the asserted pencils. QED.

The equation-cost consequence is route-specific: a loop is useful only when
the explicit section and the composed cost are checked.  The q6/orbit1307
H3 loop is one certified instance; it does not imply that every lattice
neighbour recurrence is compiler-cheap.

## 5. Integral marking transport

### Proposition D: lossless neighbour transport

Let `G` be a Gram matrix for `NS(X)` and let `A` be an integral matrix such
that

```text
det(A)=+1 or -1, and A^T G A = G.
```

Then `A` is an automorphism of the full integral NS lattice. It preserves
primitivity, intersections, discriminant form, and all class identities. If it
maps one marked `U` to another, it gives a lossless change of elliptic
fibration.

If `A` is merely rational, or integral with determinant other than `+/-1`, it
only describes a sublattice relation. Treating it as a full transport can
manufacture false MW generators or hide glue.

#### Proof

Unimodularity makes `A^{-1}` integral. The Gram identity makes it an isometry;
the remaining assertions follow functorially. QED.

This explains why an ADE/MW label is not enough: it omits the embedded `U` and
the integral transport that identifies the actual fibration.

## 6. Specialization mutation

### Theorem E: specialization balance law

Consider a smooth characteristic-zero family of K3 surfaces with compatible
Jacobian fibrations, and a generic-to-special NS specialization map. Whenever
the Picard ranks and fibre root systems are known,

```text
Delta(MW rank) = Delta(rho) - Delta(root rank).
```

#### Proof

Apply Shioda--Tate to the generic and special fibres and subtract. QED.

Consequences:

- if `rho` stays fixed, root growth forces equal MW-rank loss;
- if `rho` jumps by one, one extra algebraic class is available, but it may
  become a root, a section, or part of their glue;
- generic MW coordinates, pole orders, and component labels need not remain
  valid after specialization.

This is why the Q80 CM24 child is a typed specialization node rather than the
generic rootless endpoint.

## 7. Correctness of an equation lift

### Theorem F: conditional lattice-to-equation correctness

Let `D` satisfy Theorem C and let `O.D=1`. Suppose an exact compiler provides:

1. a complete resolved-chart cover for `O_X(D)`;
2. an exact rank calculation and two displayed independent sections
   `f_0,f_1` spanning `H0(X,O(D))`;
3. exact function-field elimination showing that `t=f_1/f_0` has generic
   fibre `C/k(t)` birational to the generic fibre of `X -> P1`;
4. a transported rational point giving the origin;
5. exact birational maps to a Weierstrass model, plus local minimality and
   fibre checks.

Then the output Weierstrass surface is the same marked Jacobian elliptic
fibration determined by `D`. Its root system and MW rank may be read from the
new fibres and Theorem A once `rho(X)` is known.

#### Proof

The complete cover and rank calculation identify the displayed two-plane with
the full `H0(X,O(D))`, so its ratio defines exactly the pencil `|D|`. The
function-field isomorphism identifies `C` with its generic fibre. The marked
rational point makes `C` an elliptic curve and identifies it with its Jacobian,
not merely with a torsor having the same invariants. Exact birational changes
produce its Weierstrass model. Relatively minimal smooth K3 models that are
birational are isomorphic, and the transported origin fixes the marking. QED.

### Why each hypothesis matters

- Matrix nullity without a complete chart cover is only a local upper bound.
- A binary quartic and its Jacobian can encode a 2-cover; point transport must
  record that degree.
- Matching `c4,c6,Delta` without the scalar-square and twist check can select a
  quadratic twist.
- Finite-place minimization alone does not classify the fibre at infinity.
- A child with the right ADE/MW label but no transported origin or inverse NS
  map is not the same marked node.

The q8 missing-`Dx` and double-2-cover failures are concrete examples of why
this theorem needs exact denominators and point-map degrees.

## 8. What a bounded neighbour search really proves

### Theorem G: completeness inside a declared lattice box

Suppose

```text
NS(X) = U + M(-1)
```

with `M` positive definite. Write a candidate isotropic vector as

```text
x = a*e + b*f + v,  with v in M(-1).
```

For fixed positive integers `a,b`, isotropy is the finite norm equation

```text
norm_M(v) = 2ab.
```

Therefore exact enumeration of all vectors of that norm, followed by
primitivity, chamber, section and marking tests, is complete for those `a,b`.
A finite box `a<=A`, `b<=B` is likewise decidable and complete.

#### Proof

Positive-definite lattices have finitely many vectors of any fixed norm. All
remaining gates are exact integer/rational tests once their wall and marking
inputs are declared. QED.

The negative boundary is essential: a completed box is not a theorem that no
larger neighbour, different chamber, or cheaper multi-step route exists.

## 9. A composite certificate theorem

The preceding results give a reusable theorem engine. A node is certified when
it contains:

```text
full NS lattice
+ primitive marked U
+ global nef/effectivity certificate
+ complete fibre-root classification
+ saturated MW/torsion/glue data.
```

An edge is certified at lattice level by a primitive isotropic target, a
replayable Weyl reduction, and unimodular forward/inverse transports. It is
certified at equation level only after the hypotheses of Theorem F pass.

Under those hypotheses, rank mutation, regulator mutation, endpoint identity,
and the equation lift are consequences rather than repeated discoveries.

## 10. What remains open

The following are useful research conjectures, not consequences of the current
examples:

1. **Low-degree connectivity:** all desired marked `U` embeddings in a fixed NS
   orbit are connected by neighbours of uniformly bounded old-fibre degree.
2. **Monotone root shedding:** a rootless fibration, when it exists, can always
   be reached by a path whose root rank never increases.
3. **Controlled equation cost:** such a path can be chosen with uniformly
   bounded pole order, resolved-RR dimension, and coefficient growth.
4. **Uniform ADE compiler:** the saturated local module is determined by a
   finite combinatorial package of resolved component and marking data.
   ADE type alone is demonstrably insufficient.
5. **Specialization transfer:** a generic high-rank K3 route yields useful
   rational specializations with independently retained section rank.

The first three would turn the present route finder into a general navigation
theorem. The fourth would turn the resolved-RR work into a reusable compiler.
The fifth is the bridge from a high generic K3 rank to new rational elliptic
curves.

## 11. Formalization order

The easiest results to formalize first are purely integral:

1. Proposition D (unimodular marked transport);
2. Lemma B1 (index-square saturation);
3. Theorem G (finite norm-shell enumeration);
4. Theorem A as an algebraic corollary once Shioda--Tate data are supplied.

Theorem C needs K3 linear-system geometry. Theorem F needs a function-field and
birational-model layer. This split lets a proof assistant verify the route
ledger now without pretending to formalize the whole equation compiler at
once.

## References

- M. Schuett and T. Shioda,
  [*Elliptic surfaces*](https://arxiv.org/abs/0907.0298), especially Sections 6
  and 11 for Shioda--Tate, heights, torsion and discriminants.
- T. Shioda,
  [*On the Mordell--Weil lattices*](https://rikkyo.repo.nii.ac.jp/records/10027),
  for the height-lattice and discriminant machinery.
- A. Kumar,
  [*Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715),
  especially Section 3.2 for primitive isotropic classes, Weyl reduction,
  genus-one pencils and section tests.
- N. Elkies and A. Kumar,
  [*K3 surfaces and equations for Hilbert modular surfaces*](https://arxiv.org/abs/1209.3527),
  for explicit K3 moduli navigation by elliptic fibrations.

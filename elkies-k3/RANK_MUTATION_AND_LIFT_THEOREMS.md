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

### Proposition B2: determinant obstruction to a rootless fibration

<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->

Suppose a Picard-rank `rho` K3 surface has a rootless Jacobian elliptic
fibration.  Put `n=rho-2` and let `D=abs(disc NS(X))`.  If `B_n` is any
proved upper bound for the Hermite constant `gamma_n`, then

```text
D >= (4/B_n)^n.
```

In particular Blichfeldt's bound

```text
B_n=(2/pi)*Gamma(2+n/2)^(2/n)
```

gives `(4/B_17)^17=28.8658...`.  A Picard-rank-19 K3 with `D<=28` therefore
cannot carry a rootless MW17 fibration.

#### Proof

The rootless MW frame is even and positive definite, hence its minimum is at
least four.  Rootlessness also makes every fibre correction in Shioda's
height formula zero.  A torsion section would then have height
`4+2(P.O)>0`, a contradiction, so MW torsion is automatically trivial.
Theorem B identifies the frame determinant with `D`.  Its Hermite invariant
is therefore at least `4/D^(1/n)`, while the definition of the Hermite
constant bounds this above by `B_n`.  Rearranging gives the claim. QED.

The determinant-24 `2E6+A2/MW3` family in
[`E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md`](E6_II_RANK3_QUADRATIC_BASE_CHANGE_2026-09-02.md)
is an exact application: the requested same-NS rootless search terminates
negatively before any neighbour enumeration.

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

### Proposition C2.2: the old-zero coefficient-swap obstruction

<!-- status-consumer: EC-K3-E6A1-RHO19-GENUINE-Q2-MW3 cd4314040bb028f7 -->

In split coordinates

```text
NS(X)=U+M(-1),     F=e,     O=f-e,
```

write an isotropic candidate as

```text
D=a*e+b*f+w,     w.M.w=2*a*b,     a,b>0.
```

Then

```text
D.O=a-b.
```

If `a<b` and `D` has the effective sign, the old zero is a fixed component.
Removing it with its exact multiplicity exchanges the two hyperbolic
coefficients:

```text
D-(b-a)O = b*e+a*f+w.
```

In particular the apparent degree-`b` presentation reduces to old-fibre
degree `a`.  A zero-neutral search at old degree `q` must therefore start at
`a>=q`; its smallest norm shell is `w.M.w=2q^2`, not `2q`.  For `q=2` and
`q=3` the first possible shells have norms eight and eighteen.

#### Proof

The identities follow immediately from `e^2=f^2=0`, `e.f=1`, and
`O=f-e`.  If `D.O<0`, the irreducible effective `(-2)`-curve `O` is fixed.
The reflection/fixed-component update is

```text
D <- D+(D.O)O = D-(b-a)O,
```

which is the displayed coefficient swap and preserves `D^2`.  Its pairing
with `F=e` is `a`, proving the degree reduction. QED.

The determinant-36 `E6+A1` K3 supplies an exact regression: all fourteen
nominal `e+2f+w`, `w^2=4`, Weyl orbits reduce to degree one, and the nominal
norm-six cubic layer has the same obstruction.  Its complete genuine
norm-eight quadratic census is
[`E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md`](E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md).

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

A related equation-level warning is essential.  Specializing a section to the
singular point of a Weierstrass `I_n` fibre does not by itself determine which
resolved component it meets: distinct tangent branches can have the same raw
node fingerprint.  Consequently, a component profile inferred only from
singular-node incidence can corrupt a Shioda height.  One exact audit is to
multiply by the exponent of the component groups and recover the canonical
height from compact pole-degree growth; a resolved local chart is still needed
when the oriented component label itself matters.  On q4/orbit164 this
fourfold audit corrects one coarse `I4` label and changes the affected height
from `3` to `13/4`.

This is why the Q80 CM24 child is a typed specialization node rather than the
generic rootless endpoint.

### Theorem E2: non-thin jumps from a second elliptic fibration

Let `K` be a number field and let `pi:X->P1_K` be a non-isotrivial elliptic
K3 fibration without non-reduced fibres. If `X` has a different elliptic
fibration over `K`, then Pasten--Salgado prove that the following are
equivalent:

```text
X(K) is Zariski dense;
pi has infinitely many rank-jump fibres;
{t in P1(K) : rank X_t(K) > rank MW(X,pi)} is not thin.
```

For the published R17 fibration, the exact `24I1` certificate gives
non-isotriviality and reduced fibres, the H3 `E7+E8/MW2` model is a different
elliptic fibration over `QQ` on the same K3, and the positive R17 section rank
gives Zariski density. Since its generic Mordell--Weil rank is exactly 17,
the rank-at-least-18 specialization locus is not thin. The complete
hypothesis audit is
[`PASTEN_SALGADO_NONTHIN_RANK_JUMPS_2026-08-31.md`](PASTEN_SALGADO_NONTHIN_RANK_JUMPS_2026-08-31.md).

<!-- status-consumer: EC-K3-R17-NONTHIN-RANK-JUMPS c9ed2e62cc456bdb -->

## 7. Correctness of an equation lift

### Proposition F0: section-first Tate charts

Let `k` be a field of characteristic different from two and let `K=k(t)`.
Suppose

```text
E: y^2=x^3+A*x+B
```

has a `K`-point `P=(xp,yp)` with `yp != 0`.  The unit Weierstrass change

```text
x=X+xp,  y=Y+m*X+yp,  m=(3*xp^2+A)/(2*yp)
```

puts `P` at `(0,0)` and gives

```text
Y^2+a1*X*Y+a3*Y=X^3+a2*X^2,
a1=2*m,  a2=3*xp-m^2,  a3=2*yp.
```

Thus a one-marked-section equation has no residual section condition.

More generally let `R=k[t]` and choose `a1,h,r,s,kappa in R` satisfying

```text
gcd(r,h)=gcd(r,s)=1.
```

Choose `alpha,beta in R` with `alpha*s+beta*r^2=1` and put

```text
T  = h*(r^3-h*s^2-a1*r*s),
a3 = alpha*T+kappa*r^2,
a2 = -beta*T+kappa*s.
```

Then

```text
Y^2+a1*X*Y+a3*Y=X^3+a2*X^2
```

contains the two marked points

```text
P=(0,0),  Q=(h*r,h^2*s),
```

and their affine coincidence ideal has gcd `h`.  If `h` is coprime to the
discriminant and the chart has no omitted pole or infinity intersections,
then `P.Q=deg(h)` on this affine chart.  Prescribed semistable fibres may
therefore be imposed on the compiled discriminant: exact order `n` of
`Delta` together with a `c4` unit is the local `I_n` gate.

#### Proof

Direct substitution after the first change of variables cancels the constant
and linear `X` terms and gives the displayed coefficients.  It is a unit
Weierstrass transformation, so it preserves `c4,c6,Delta`.

For the two-point construction, substitution of `Q` and division by `h^2`
reduces its equation to

```text
a3*s-a2*r^2 = h*(r^3-h*s^2-a1*r*s)=T.
```

The definitions of `a2,a3` make the left side

```text
T*(alpha*s+beta*r^2)=T.
```

The first point is immediate.  Finally
`gcd(h*r,h^2*s)=h` follows from the two coprimality hypotheses.  At a smooth
zero of `h`, `r` is a unit, so the local coincidence ideal is generated by
`h`; multiplicities give `deg(h)`.  The semistable fibre statement is the
corresponding case of Tate's algorithm. QED.

This proposition is an ansatz compiler, not an exact-rank theorem.  It builds
in marked points and their affine intersection relation, but independence,
resolved component labels, global minimality, torsion/divisibility, Picard
rank, and NS saturation remain separate gates.  Polynomial elliptic-K3
searches use `deg(a_i)<=2i`; translating an already known short model can
produce rational-function `a_i`, so new searches should begin in the chart.
The implementation and Golay/NS0031 controls are in
[`SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md`](SECTION_FIRST_NORMAL_FORM_COMPILER_2026-09-02.md).

A route-specific obstruction illustrates why the chosen section chart is part
of the hypotheses.  In the rationalized `D6` polynomial marked-section chart,
the two nontrivial equal-leading-coefficient correspondences are both
birational to `Y^2=X^3+X^2+X`.  This elliptic curve has rank zero over `QQ`
and only the degenerate boundary torsion, so that chart contains no
nontrivial rational pair.  The result is exact for the declared chart and says
nothing about a larger rational-function or section-first `D6` chart.  See
[`LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md`](LOWER_ROOT_TWO_TWIST_SEARCH_2026-09-02.md).

<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->

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

### Proposition F0b: quadratic-elimination parent gate

Let `R` be an integral domain of characteristic different from two, and let

```text
q(X)=a*X^2+b*X+c in R[X].
```

An equation lift may pass the quadratic-discriminant gate only after checking
both

```text
parent(q)=R[X],
q.discriminant()=b^2-4*a*c in R.
```

In particular, a coefficient that simplifies into `R` must be explicitly
coerced back to `R` before forming `q`.  Leaving it in `Frac(R)` can promote
the ambient polynomial to a different coefficient tower, where a
zero-argument `discriminant()` method need not dispatch to the intended
variable or parent.  Factorization, square stripping, and binary-quartic
invariants computed before this gate do not certify the neighbour.

#### Proof

For a quadratic over `R`, the discriminant is identically `b^2-4ac`; equality
after coercion is an exact algebraic check.  Every later squarefree-quartic and
Jacobian operation is functorial in that element, so a different element
computes a different genus-one curve. QED.

The orbit-96 `A7+D7` lift is the motivating counterexample.  Its tangent slope
is polynomial in the old base but was retained in a fraction-field parent.
The resulting spurious degree-three residual reproduced the old `2E6+A3`
frame.  The parent gate gives the genuine degree-four residual and the
`I8+I3*+7I1` Jacobian; see
[`E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md`](E6A1_RHO19_ORBIT96_WEIERSTRASS_GALOIS_2026-09-02.md).

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT96-A7D7-GALOIS ba008502f0e5533f -->

### Proposition F0c: unordered incidence does not certify section descent

<!-- status-consumer: EC-K3-RES-QBC-E6-RANK4-LINEAR-CHORD 3bcfe3534656b26f -->

Let `C -> C0` be a generically quadratic ordered-incidence cover over a field
`K`, with involution `sigma`, and suppose an elliptic surface descends to
`K(C0)` while displayed sections on `C` occur in exchanged pairs

```text
sigma(P)=Q,       sigma(T1)=T2.
```

Then rationality of `C0` does not imply that the individual sections descend.
If the four displayed sections span the full geometric Mordell--Weil space and
the two exchanged pairs span independent two-planes, the Mordell--Weil rank
over `K(C0)` is exactly two, with rational span

```text
P+Q,       T1+T2.
```

Rank four is defined only over `K(C)`.  In particular, a genus computation or
rational point search performed after quotienting by the section-label
involution must be followed by the squareclass and fixed-subspace descent
gate before it is advertised as a rational one-parameter rank source.

#### Proof

Galois descent identifies the rational Mordell--Weil group tensored with
`QQ` with the `sigma`-fixed subspace of the geometric group.  The fixed
subspace of each exchanged two-dimensional permutation representation is the
line generated by the sum.  Independence of the two planes gives the
displayed two-dimensional fixed subspace. QED.

The E6 linear-chord incidence is the exact counterexample motivating the
gate.  Its unordered quotient is `P1_QQ`, but recovering `v,w` introduces

```text
r^2=k^4+6*k^2+13.
```

This ordered cover is the genus-one curve `52a2`, not a conic.  It has rank
zero and rational torsion `Z/2`; its only rational points lie at infinity.
The complete parameterization, ordered-cover map, quadratic witness, and
rank-descent calculation are in
[`E6_RANK4_LINEAR_CHORD_INCIDENCE_2026-09-02.md`](E6_RANK4_LINEAR_CHORD_INCIDENCE_2026-09-02.md).

### Proposition F1: direct bisection compilation from a height-ten trace

Let

```text
E: y^2=x^3+A(t)x+B(t)
```

be an integral rootless elliptic K3 over a characteristic-zero field, in the
standard degree bounds `deg(A)<=8`, `deg(B)<=12`.  Let a height-ten section be
written in coprime form

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=3,       gcd(Nx,h)=1.
```

There is a unique polynomial `M` of degree below six satisfying

```text
M*Nx+Ny == 0 mod h^2.
```

Put

```text
U = M^2-Nx,
R = M*Nx+Ny,
N = M^4-6*M^2*Nx-8*M*Ny-3*Nx^2-4*A*h^4.
```

Then `h^2` divides `U`, `h^6` divides `N`, and the residual intersections of
the line through `-tau` with slope `M/h` satisfy

```text
x^2-(U/h^2)*x+(R^2-B*h^6)/(h^4*Nx)=0.
```

Its discriminant is

```text
h^2*q(t),       q=N/h^6,       deg(q)<=2.
```

If the class of `tau` modulo twice the Mordell--Weil lattice is one of the
section-nonnegative norm-ten bisection classes, then the residual curve is the
corresponding irreducible rational bisection.  In particular `q` is a
non-square squarefree quadratic after removing a rational square, and its
class is the exact quadratic branch squareclass.

#### Proof

Invertibility of `Nx` modulo `h^2` gives existence and uniqueness of `M`.
The section identity is

```text
Ny^2=Nx^3+A*Nx*h^4+B*h^6.
```

The congruence for `M` and this identity first give `h^2 | U`.  In the
localization at `h`, write `M=-Ny/Nx+h^2*k`; then

```text
U/h^2 == -2*(Ny/Nx)*k mod h^2,
R/h^2 == k*Nx mod h^2.
```

The exact identity

```text
Nx*N = Nx*U^2-4*R^2+4*B*h^6
```

therefore shows `h^6 | N`; coprimality removes `Nx`.  Substituting the line
through `-tau` into the cubic and removing its known root gives the displayed
quadratic.  Its discriminant is `N/h^4=h^2*q`.  The K3 degree bounds give
`deg(N)<=20`, hence `deg(q)<=2`.  Finally, the lattice argument in
[`BISECTION_COLLISION_SEARCH.md`](BISECTION_COLLISION_SEARCH.md) proves that a
surviving class is an irreducible smooth rational bisection, excluding a split
or constant residual cover and identifying its trace class with `tau mod 2M`.
QED.

This proposition replaces a nonlinear generic Riemann--Roch solve at the
rootless endpoint by one exact elliptic group-law computation and one linear
polynomial inversion modulo `h^2`.  It does not predict collisions between the
resulting quadratic squareclasses.

If a pole of `tau` lies over infinity, apply the proposition after the base
chart change `s=1/t`, with `x_s=s^4*x` and `y_s=s^6*y`, and transport the
result back.  The cover coordinate transforms by
`q_t=t^2*q_s(1/t)`.  This is again multiplication by the required square under
`u_t=t*u_s`, so it preserves the quadratic squareclass and all displayed
coefficient identities.  In the complete published-R17 batch this chart is
needed only for orbit `0x0c54f`.

### Proposition F1.1: the height-eight genus-one bisection pencil

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-BISECTION-PILOT 80fa6e59107cc9e6 -->

In the setup of Proposition F1, let `tau` instead have height eight and a
finite-pole presentation

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=2,       gcd(Nx,h)=1.
```

Let `M0` be the unique polynomial of degree below four satisfying

```text
M0*Nx+Ny == 0 mod h^2.
```

For `lambda` in the ground field put

```text
M_lambda=M0+lambda*h^2
```

and define `N_lambda` by the formula for `N` in Proposition F1.  Then

```text
h^6 divides N_lambda,
q_lambda=N_lambda/h^6,
deg(q_lambda)<=4.
```

The line of slope `M_lambda/h` through `-tau` cuts out the residual quadratic
from Proposition F1.  If `q_lambda` is squarefree of degree four, its
normalization is a genus-one bisection.  Over
`s^2=q_lambda(t)` its two lifts are

```text
x=(M_lambda^2-Nx)/(2*h^2) + (h/2)*s,
y=y0+(M_lambda/2)*s,
```

where `y0` is obtained by substituting the constant part of `x` in the line.
Their sum is the pullback of `tau`.

More precisely, let `t0` be rational with `h(t0)!=0`, and let
`Q=(x_Q,y_Q)` be a rational point of the smooth fibre with
`x_Q!=x(tau(t0))`.  There is a unique member of the pencil through `Q`, namely

```text
lambda_Q =
  (h(t0)*(y_Q+y(tau(t0)))/(x_Q-x(tau(t0)))-M0(t0))/h(t0)^2.
```

At the corresponding rational point of the cover,

```text
s_Q=(2*x_Q-(M_lambda^2-Nx)/h^2 evaluated at t0)/h(t0),
```

one has `s_Q^2=q_lambda(t0)` and the lifted point is literally `Q`.
Consequently its Kummer barcode is exactly `x_Q-theta`, not merely a class
selected by a numerical or local proxy.

If in addition `q_lambda` is coprime to the surface discriminant and to `h`,
and the displayed lift has the standard integral degree bounds, the
anti-invariant section is independent of the invariant Mordell--Weil group
and has height sixteen.

#### Proof

The congruence and divisibility proof is unchanged from Proposition F1,
because every `M_lambda` has the same residue modulo `h^2`.  A height-eight
section on a rootless elliptic K3 has `tau.O=2`; in the finite chart this gives
`deg(h)=2`, `deg(Nx)<=8`, and `deg(Ny)<=12`.  Thus
`deg(M_lambda)<=4`, every term of `N_lambda` has degree at most sixteen, and
division by `h^6` leaves degree at most four.  The residual quadratic and
lift formulas follow by direct substitution.  Riemann--Hurwitz gives genus
one when the branch divisor consists of four simple points.

Solving the line-incidence equation at `t0` gives the displayed unique
`lambda_Q`; substituting it gives both the cover witness and literal equality
with `Q`.  Finally the degree-two pullback has `chi=4`, so an integral section
disjoint from zero has self-intersection `-4`.  Coprimality with `h` makes the
two conjugate lifts meet transversely at exactly the four branch points and
nowhere else.  Their difference therefore has height

```text
2*(4-(-4))=16.
```

It is anti-invariant under the deck involution and hence orthogonal to, and
independent of, the invariant Mordell--Weil subgroup. QED.

The first published-R17 positive-control application uses the
equation-cheapest trace `tau=-P2-P5`.  All eleven exceptional rank-28 targets
select rational pencil parameters; their branch polynomials are irreducible
squarefree quartics coprime to the 24 singular fibres.  This `11/11` incidence
is the expected base-point-free pencil behavior.  The equation-level content
is the exact quartic, smooth-branch, Kummer, lift, and height certificate, not
a claim that fitting a pencil member through a known point discovers that
point.  See
[`R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md`](R17_RANK28_GENUS_ONE_BISECTIONS_2026-09-02.md).

### Proposition F1.2: the height-twelve regular quartic

<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->

In the setup of Proposition F1, let `tau` have height twelve and a finite-pole
presentation

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=4,       gcd(Nx,h)=1.
```

Let `M0` be the unique polynomial of degree below eight satisfying
`M0*Nx+Ny=0 mod h^2`, and form `N` as in Proposition F1.  Then

```text
h^6 divides N,       q=N/h^6,       deg(q)<=4.
```

Thus, when `q` is squarefree of degree four, the unique regular residual chord
is a genus-one bisection.  Unlike Proposition F1.1, there is no scalar pencil
parameter within the regular degree bound.

#### Proof

The congruence, divisibility, residual quadratic, and lift identities are the
same as in Proposition F1.  The elliptic-K3 weights give
`deg(Nx)<=12`, `deg(Ny)<=18`, and `deg(M0)<=7`.  Hence every term in `N` has
degree at most 28, while `deg(h^6)=24`; therefore `deg(q)<=4`.  Adding a
nonzero multiple of `h^2` raises the slope numerator to degree at least eight
and leaves the regular degree bound, proving uniqueness.  Riemann--Hurwitz
gives genus one for a squarefree quartic. QED.

On the published R17 frame the complete 43-element norm-twelve deep set from
Proposition F5 has `deg(h)=4`; exact construction gives 43 irreducible
squarefree quartics, all coprime to the surface discriminant and their trace
denominators.  The complete calculation and its bounded simultaneous-splitting
search are recorded in
[`R17_GENUS_ONE_BISECTION_SPLITTING_SEARCH_2026-09-02.md`](R17_GENUS_ONE_BISECTION_SPLITTING_SEARCH_2026-09-02.md).

### Theorem F2: complete injectivity on the published rootless R17 survivor set

For the published rootless R17 elliptic K3, let `C` be the 39,120
section-translation classes of section-nonnegative degree-two `(-2)`-curves
enumerated in [`BISECTION_COLLISION_SEARCH.md`](BISECTION_COLLISION_SEARCH.md).
The map

```text
C -> QQ(t)^*/QQ(t)^{*2},       [B] -> [q_B]
```

which sends a bisection to its quadratic branch extension is injective.
Every `q_B` is a squarefree quadratic coprime to the surface discriminant.
Consequently every class gives an explicit smooth quadratic base change of
generic Mordell--Weil rank at least 18, while no two distinct classes give a
common quadratic base change.  In particular this complete bisection set
cannot yield a rank-two anti-invariant collision on one quadratic cover.  This
does not exclude the distinct-extension composita in Theorem F3 below.

#### Proof

The complete norm-ten shell contains 806,238 unoriented representatives and
maps onto exactly the 39,120 surviving classes.  The exact priority replay
tests every representative and retains one published-basis trace per class.
Proposition F1 constructs its quadratic relation; coefficientwise identities
verify both lifted points on the Weierstrass equation.  The single trace with
a pole at infinity is handled by the reciprocal chart above.  The complete
coverage gate checks the pinned integral vector and its class modulo `2R17`
for every record.

For all 39,120 records the computed `q_B` has degree two, is squarefree, and
is coprime to the degree-24 surface discriminant.  Independent exact
normalization of the displayed quadratic discriminants gives 39,120 distinct
keys in `QQ(t)^*/QQ(t)^{*2}`.  This proves injectivity.

On each double cover the two lifted sections meet transversely at the two
simple branch points.  The pullback fibration remains rootless, has `chi=4`,
and for one lift `P` and the deck involution `sigma` one has
`P^2=-4` and `P.sigma(P)=2`.  Hence the anti-invariant section has height

```text
<P-sigma(P),P-sigma(P)> = 2*(2-(-4)) = 12.
```

It is non-torsion and adds one direction to the invariant rank-17 lattice.
The absence of equal squareclasses excludes a two-bisection common-cover
height matrix.  QED.

<!-- status-consumer: EC-K3-BISECT-EQUATION-BATCH a0570a5a4ea8e02b -->

### Corollary F2.1: translated trace shells cannot enlarge specialization visibility

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES 1c39220ee5fedc77 -->

Let `B` be one of the complete degree-two bisection classes of Theorem F2 and
let `S` be a section of the generic `R17` subgroup.  Translating `B` fibrewise
by `S` preserves its quadratic branch extension.  If `t_0` is a rational good
fibre at which the extension splits and one branch gives `P` in `E_{t_0}(QQ)`,
then the translated branch gives

```text
P + S(t_0).
```

Consequently `P` and its translate define the same class modulo the specialized
generic subgroup.  Inversion changes the class to its negative, which is the
same class in every quotient modulo two.  Therefore any higher-height trace
shell consisting only of translations or inversions of the complete 39,120
classes has exactly the same split-extension set and the same finite-quotient
visibility span at every rational good fibre.

#### Proof

Fibrewise translation by a section is an automorphism of the smooth locus over
the base and does not change the degree-two map from the normalization of `B`
to the parameter line.  It therefore does not change the corresponding
quadratic function field.  On a split fibre it sends each rational branch
point to its elliptic sum with `S(t_0)`.  Passing to the quotient by the
specialized generic subgroup removes this summand.  Finally `-P` and `P` have
the same image after tensoring the quotient with `F_2`.  Theorem F2 says that
the 39,120 stored classes already exhaust the relevant translation classes,
so no translated trace shell can add a new one. QED.

This corollary is a mechanism boundary, not a point-search obstruction.  A
higher-degree multisection, a different covering construction, or a direct
specialization point can still occupy a quotient direction invisible to the
bisection atlas.

### Theorem F3: distinct bisection extensions give genus-one rank-19 bases

For any two distinct classes `B_1,B_2` in the complete published-R17 survivor
set, let

```text
C_12: u^2=q_1(t),  v^2=q_2(t).
```

Then `C_12` is a geometrically connected genus-one `V4` cover of `P1`.  Over
its function field the two pulled anti-invariant sections have exact height
matrix

```text
[24  0]
[ 0 24],
```

and the pulled elliptic surface has generic Mordell--Weil rank at least 19.
All 39,120 individual conics are rational over `QQ`, and the complete set
therefore gives exactly `binomial(39120,2)=765167640` such paired bases.

#### Proof

Exact factorization shows every primitive `q_i` is an irreducible quadratic
and no two are proportional.  Hence distinct `q_i,q_j` are independent in
`QQ(t)^*/QQ(t)^{*2}`, their geometric branch sets are disjoint, and the
compositum has Galois group `V4`.  It has four branch points with inertia order
two, so Riemann--Hurwitz gives

```text
2g(C_12)-2 = 4*(-2)+4*(4-2)=0.
```

The height-12 direction from each double cover doubles to height 24 after the
other degree-two pullback.  The two directions occupy distinct nontrivial
`V4` characters; Galois invariance makes their cross-height zero.  They are
also orthogonal to the invariant rank-17 space, proving the rank bound.

Exact Hasse--Minkowski computation supplies a rational point on every
individual conic.  It is not claimed that every paired genus-one curve has a
rational point.  Exactly 5,566 pairs have an immediate common point at zero or
infinity.  A complete bounded point ledger for those curves has two certified
rank-at-least-nine bases; an empty bounded search remains only lower bound
zero, not an exact-rank statement.  QED.

<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->

### Theorem F4: multiquadratic character decomposition and base genus

<!-- status-consumer: EC-K3-BISECT-MULTIQUADRATIC-CHARACTERS dc58103d8d2494cf -->

Let `K` be a field of characteristic different from two, let `E/K` be an
elliptic curve, and let `q_1,...,q_k` be independent elements of
`K^*/K^{*2}`.  Put

```text
L=K(sqrt(q_1),...,sqrt(q_k)),
q_S=product(q_i, i in S),
```

with `q_empty=1`.  Then there is an exact rational character decomposition

```text
E(L) tensor QQ
  = direct_sum over S subset {1,...,k} of E^{q_S}(K) tensor QQ,
```

and in particular

```text
rank E(L) = sum over S rank E^{q_S}(K).
```

Distinct summands are orthogonal for every Galois-invariant canonical height
pairing.  Thus ranks on product twists are genuinely new character
contributions: they cannot be absorbed into the original curve or the
singleton twists.

Now take `K=QQ(t)` and suppose each `q_i` is a squarefree quadratic whose
reduced geometric branch divisor on `P1` is disjoint from every other one.
Let `C_k` be the smooth projective curve with function field `L`.  Then `C_k`
is a geometrically connected `2^k`-sheeted cover of `P1` and

```text
g(C_k) = 1 + 2^(k-1)*(k-2).
```

If `E(QQ(t))` has rank `r`, and each singleton twist has one known non-torsion
direction, then

```text
rank E(QQ(C_k)) >= r+k.
```

If the known direction on each individual double cover has height `h`, its
pullback to `C_k` has height `2^(k-1)*h`.  The `k` known directions therefore
have diagonal height block

```text
2^(k-1)*h * I_k.
```

For the published rootless R17 surface, `r=17` and `h=12`, so the base has
genus `1+2^(k-1)*(k-2)`, rank at least `17+k`, and new height block
`12*2^(k-1)*I_k`.  More precisely, every nonempty product twist contributes
its full rank to the corresponding additional character.  For two covers,

```text
rank E(QQ(t)(sqrt(q_i),sqrt(q_j)))
  = 17 + rank E^{q_i}(QQ(t)) + rank E^{q_j}(QQ(t))
       + rank E^{q_i*q_j}(QQ(t)).
```

Consequently either of the following would improve the current constructions:

- `rank E^{q_i}(QQ(t))>=2` gives a rational `P1` base of generic rank at
  least 19;
- `rank E^{q_i*q_j}(QQ(t))>=1` gives the associated genus-one paired base
  generic rank at least 20.

#### Proof

Let `G=Gal(L/K)`.  Since `G` is an elementary abelian two-group, the rational
group algebra has the orthogonal idempotents

```text
e_chi = 2^(-k) * sum(g in G) chi(g)*g.
```

They split `E(L) tensor QQ` into its `2^k` character eigenspaces.  Over `L`,
the standard isomorphism from `E^{q_S}` to `E` identifies
`E^{q_S}(K) tensor QQ` with the eigenspace for the character attached to
`sqrt(q_S)`.  This proves the direct sum and rank formula.  If points `P,Q`
belong to different characters, choose `g` on which the characters differ.
Galois invariance and bilinearity give

```text
<P,Q> = <gP,gQ> = -<P,Q>,
```

so their cross-height is zero.

Disjoint nonempty branch divisors make the squareclasses geometrically
independent: every nonempty product has a branch point of odd valuation.
Thus the cover is geometrically connected.  It has `2k` branch points.  Over
each one there are `2^(k-1)` points of ramification index two, so
Riemann--Hurwitz gives

```text
2g(C_k)-2 = 2^k*(-2) + 2k*2^(k-1) = 2^k*(k-2).
```

This is the displayed genus.  The rank lower bound follows by retaining the
trivial character and the `k` singleton characters.  Canonical heights on an
elliptic surface multiply by the degree of a finite base change.  Pulling an
individual double-cover direction through the remaining degree `2^(k-1)`
base change therefore multiplies its height by `2^(k-1)`; character
orthogonality makes the resulting block diagonal. QED.

The theorem is exact, but it does not determine any twist rank.  The
Frobenius-character census in
[`QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md`](QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md)
is only a candidate-ranking mechanism until an additional rational section
and its independence are certified.

### Proposition F5: rootless low-degree multisections are coset minima

Let a rootless elliptic K3 have

```text
NS(X) = U + M(-1),
```

where `F=e` is the fibre, `e.f=1`, and `M` is positive definite and even.
Every divisor class with fibre degree `d>0`, arithmetic genus `g`, and
`M(-1)` coordinate `w` has the form

```text
D = ((norm_M(w)+2g-2)/(2d))*e + d*f + w.
```

It is integral precisely when the displayed first coefficient is integral.
Translation by the section indexed by `x in M` replaces `w` by `w+d*x`.
Moreover, for the section

```text
S_x = ((norm_M(x)-2)/2)*e + f + x,
```

one has

```text
2d*(D.S_x) = norm_M(w-d*x) - (2d^2-2g+2).
```

Consequently `D` is nonnegative on every section if and only if the exact
minimum of its coset in `M/dM` is at least

```text
2d^2-2g+2.
```

For `(d,g)=(2,0)` the threshold is ten.  Every such effective class is an
irreducible smooth rational bisection: rootlessness removes vertical root
components, while a decomposition into two sections would have intersection
one by adjunction and hence negative intersection with either component.
For `g>=1` or `d>=3`, the same calculation certifies the lattice class and
all-section nonnegativity, but not global nefness, irreducibility, arithmetic
descent, or a Mordell--Weil rank gain.

#### Proof

The formula for `D` is exactly the adjunction equation `D^2=2g-2`.  Substituting
the displayed section class gives the completed-square identity.  Translation
therefore preserves the residue class of `w` modulo `dM` and the minimum over
all sections is the corresponding positive-definite coset minimum.  When
`d=2,g=0`, Riemann--Roch and `D.F>0` make `D`, rather than `-D`, effective.
If it split horizontally, adjunction forces two rational degree-one
components meeting once, so `D` would meet either component in `-1`, contrary
to section nonnegativity.  An irreducible arithmetic-genus-zero curve on a
smooth K3 is smooth and rational. QED.

The complete degree-two, complete selected-frame degree-three, and bounded
sampled degree-four applications are recorded in
[`LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md`](LATTICE_FOUNDRY_SOURCE_FIRST_OBJECTIVE_2026-09-01.md).
The degree-three certificate exhausts all `3^17` cosets on each selected
rootless frame; it does not strengthen the geometric boundary in the
proposition.  Coset abundance is only a discovery coordinate; the published
R17 experience shows that it is not by itself a predictor of exceptional
specialization rank.

### Proposition F6: the intrinsic multisection-coset metric and degree overlap

In the setup of Proposition F5, define

```text
mu_d(c) = min { norm_M(w) : w mod dM = c },       c in M/dM.
```

Let `C` and `D` be degree-`d` divisor classes of arithmetic genera `g` and
`h`, with horizontal coordinates in cosets `c` and `c'`.  Then their minimum
intersection under independent translations by sections is

```text
min C.D = mu_d(c-c')/2 + g + h - 2.
```

Thus `mu_d(c-c')` is an intrinsic translation-quotient metric.  Any threshold
graph or hypergraph defined from it is preserved by `Aut(M)` and does not
depend on chosen shortest representatives.  In particular, representative
angle distributions are useful equation gauges but are not quotient
invariants unless a representative-selection rule is pinned.

Regard `M/dM` as the `d`-torsion subgroup `(1/d)M/M` of the real lattice
torus.  If `d` divides `e`, the natural inclusion is

```text
c mod dM  |->  (e/d)c mod eM,
```

and its coset minima satisfy the exact scaling law

```text
mu_e((e/d)c) = (e/d)^2 * mu_d(c).
```

For arbitrary positive `d,e`, the literal intersection of their torsion
subgroups is the `gcd(d,e)`-torsion subgroup.  Hence coprime degree structures
meet only at zero; a stronger comparison between them requires an explicitly
defined common-modulus or CRT compatibility relation rather than a
representative-dependent overlap count.

#### Proof

Choose representatives `w+d*x` and `v+d*y`.  Substitution of the adjunction
coefficients from Proposition F5 gives

```text
C.D = norm_M((w+d*x)-(v+d*y))/2 + g + h - 2.
```

As `x-y` ranges over `M`, minimizing gives the first formula.  Lattice
automorphisms preserve norms, differences and congruence classes, proving the
invariance statement.

For `d|e`, every representative of `(e/d)c mod eM` has the form

```text
(e/d)w + e*x = (e/d)(w+d*x),
```

so taking norms and minima proves the scaling law in both directions.  The
torsion-intersection statement follows coordinatewise from Bezout, or from
the elementary identity between the `d`- and `e`-torsion subgroups of a free
real torus. QED.

The first complete R17 application is
[`R17_MULTISECTION_DIVERSITY_2026-09-02.md`](R17_MULTISECTION_DIVERSITY_2026-09-02.md).
It finds, among other things, that the 39,120 rational bisection vertices form
one connected zero-intersection graph, while their natural degree-four images
are genus-one quadrisection vertices of minimum norm 40.  These are exact
lattice statements; they do not promote the sampled degree-three or
degree-four graph data to geometric curves.

## 7A. Integral bridge reglue

### Theorem H-1: relative-`U` bridge lifting

Let

```text
L = U direct_sum W(-1),        J = [0 1; 1 0],
```

where `W` is positive definite, and fix the ordered basis `(u_1,u_2)` of
`U` with Gram matrix `J`.  Let `(u'_1,u'_2)` be an ordered basis of another
copy `U'` of `U` in `L`, and put

```text
A_ij = <u_i,u'_j>.
```

There are unique vectors `w_1,w_2 in W` such that, column by column,

```text
u'_j = (J*A)_(bullet j) + w_j.                       (H-1.1)
```

Here the first term is written in the fixed `U` basis and the second is
viewed in `W(-1)`.  Their positive Gram matrix is

```text
Gram_W(w_1,w_2) = G_A := A^t*J*A-J.                 (H-1.2)
```

Conversely, for any `A in M_2(ZZ)` and any ordered pair `(w_1,w_2)` in `W`
with Gram matrix `G_A`, formula (H-1.1) gives a sublattice with Gram `J`.
It is automatically primitive and splits off integrally.  Thus, for fixed
`A`, (H-1.1) is a bijection between ordered representations of `G_A` in `W`
and ordered embedded copies `U'` having cross-pairing matrix `A`.

Put `B=<w_1,w_2>`.  If `G_A` is positive definite, equivalently
`rank(U+U')=4`, then

```text
K = W intersect U'^perp(-1) = B^perp in W,
C = K^perp in W = saturation_W(B).                  (H-1.3)
```

Consequently, if `m=[C:B]`, then

```text
det(G_A) = m^2*det(C).                               (H-1.4)
```

The reverse projection has cross matrix `A^t` and raw positive Gram
`A*J*A^t-J`; its saturation is the bridge on the `W'=U'^perp(-1)` side.
When `G_A` has rank below two, the same projection identities hold, but this
is a degenerate relative position and there is no rank-two bridge.

#### Proof

Write the `U`-coordinate column of `u'_j` as `c_j`.  Pairing with the fixed
basis gives `J*c_j=A_(bullet j)`, hence `c_j=J*A_(bullet j)`.  The `U` part
therefore has Gram

```text
(J*A)^t*J*(J*A)=A^t*J*A.
```

Since the full Gram is `J` and `W(-1)` has the negative of the form on `W`,
subtraction gives (H-1.2).  The converse is the same calculation backwards.
The resulting `U'` is unimodular.  Orthogonal projection to `U'` is integral,
so `L=U' direct_sum U'^perp`; in particular the embedding is primitive.

For `x in W`, equation (H-1.1) gives

```text
<x,u'_j>_L = -(x,w_j)_W.
```

This proves the first equality in (H-1.3).  Positive definiteness of `G_A`
makes `B` rank two.  Taking the orthogonal complement of `K` over `QQ` gives
`B tensor QQ`, and intersecting back with `W` is exactly the saturation of
`B`.  The determinant-index identity is Lemma B1. QED.

For elliptic markings

```text
u_1=F,       u_2=O+F,
u'_1=F',     u'_2=O'+F',
```

the cross matrix is exactly

```text
    [ d       d+s       ]
A = [ d+t     d+s+t+z   ],                          (H-1.5)
```

where

```text
d=F.F',       s=F.O',       t=O.F',       z=O.O'.
```

This interpretation requires `O=u_2-u_1` and `O'=u'_2-u'_1` to be the
actual effective zero curves.  An arbitrary unimodular splitting mate still
satisfies (H-1.1)--(H-1.4), but its derived `(-2)` class can be a chamber
pseudo-zero; then `s,t,z` are lattice intersections rather than physical
compiler degrees.

The opt-in exact replay in
[`certify_integral_rank_transfer_bridge_reglue.sage`](scripts/certify_integral_rank_transfer_bridge_reglue.sage)
checks (H-1.1)--(H-1.4) in both orientations on all 42 stored H3, Q80,
NS0024, and Golay-720 clicks.  It recovers every stored old-fibre degree from
`A_11` and recovers each stored bridge as the saturation of the two projected
vectors.  All 84 oriented projected pairs happen to be saturated (`m=1`), so
these controls verify the square-index law but do not exhibit a nontrivial
index.  Some Q80 transported mates give negative values among `s,t,z`, which
is an exact warning that the 42-edge lattice ledger does not always carry an
equation-effective zero.  The generated record is
[`elkies-k3-relative-u-bridge-lifting-regression-v1.json`](../artifacts/generated-results/elkies-k3-relative-u-bridge-lifting-regression-v1.json).

#### Corollary H-1a: bounded relative-marking completeness

Fix finite sets of allowed values for

```text
F.F',       F.O',       O.F',       O.O'.
```

Then there are finitely many matrices (H-1.5).  For each matrix, there are
finitely many ordered representations of `G_A` in `W`, because each `w_j`
lies on one fixed norm shell of the positive-definite lattice.  Formula
(H-1.1) therefore gives a terminating complete enumeration of all ordered
`U'` in the declared intersection box.  For every output one can exactly:

1. construct the primitive `U'` and its full integral marking;
2. compute `W'=U'^perp(-1)`, its roots, and its saturated bridge data;
3. use Shioda--Tate to determine MW rank once `rho` and the fibre roots are
   certified;
4. test integral isometry against a declared foundry target;
5. Weyl-reduce `F'`, run the finite horizontal-wall test of Proposition C2,
   and audit whether `O'` is an effective zero;
6. retain the exact marking for the conditional equation compiler of
   Theorem F.

The same conclusion holds for a declared physical cost function only when a
cost bound implies a computable finite set of the four intersection tuples.
This coercivity hypothesis is essential: a bound on resolved-RR dimension,
coefficient height, or another downstream score does not by itself bound
`s,t,z`.  The result is completeness in the declared relative-marking box,
not global finiteness of literal `U` embeddings, an equation construction, or
a universal compiler-cost bound.

The relative projection is elementary lattice linear algebra.  Its geometric
use relies on the standard correspondence between Jacobian elliptic
fibrations and primitive `U` embeddings, while primitivity, overlattices, and
saturation are governed by Nikulin's discriminant-form formalism cited below.

### Theorem H: bridge-core root transfer

Let `L` be a nondegenerate even integral lattice and let `U_0,U_1` be copies
of the unimodular hyperbolic plane such that the positive-sign frame lattices

```text
W_i=U_i^perp(-1)
```

are positive definite.  Regard both frames as sublattices of the common
ambient `L` and put

```text
K=W_0 intersect W_1.
```

Then `K` is primitive in both frames.  For

```text
C_i=K^perp inside W_i
```

the lattice `C_i` is primitive in `W_i`, `K+C_i` has finite index in `W_i`,
and `W_i` is the even overlattice encoded by an isotropic subgroup

```text
H_i subset A_K + A_C_i.
```

Both projections of `H_i` to the two discriminant groups are injective, so
`H_i` is the graph of an anti-isometry between its two images.  Moreover the
root sets in the common ambient lattice satisfy the exact identity

```text
Phi(W_0) intersect Phi(W_1) = Phi(K).
```

Consequently every root outside `K` that disappears or appears under the
change of `U` is confined to the finite bridge/glue presentations
`K+C_0 subset W_0` and `K+C_1 subset W_1`.  In particular `W_i` is rootless
if and only if `K` and every glue coset selected by `H_i` contain no
norm-two vector.

If `U_0+U_1` has rank four, then `K` has rank `rank(L)-4`; for a Picard-rank
19 elliptic K3 this is rank 15 and each `C_i` has rank two.  If in addition

```text
|H_i|=|A_C_i|,
```

then the projection `H_i -> A_C_i` is an isomorphism: the frame is obtained
by maximal graph glue over the entire bridge discriminant group.

#### Proof

Unimodularity of `U_i` gives the integral splitting

```text
L=U_i direct_sum U_i^perp.
```

Thus each `W_i` is primitive in `L`.  If `n*x` lies in `K` for nonzero
integer `n` and `x` lies in `W_i`, then `n*x` lies in the other primitive
frame, hence `x` lies in that frame and therefore in `K`.  This proves that
`K` is primitive in each `W_i`.  The same divisibility argument for an
orthogonal complement proves that `C_i` is primitive.

Positive definiteness makes `K` nondegenerate, so `K+C_i` has full rank in
`W_i` and finite index.  The standard even-overlattice correspondence gives
the isotropic subgroup `H_i` in `A_K+A_C_i`.  If an element of `H_i` projects
to zero in `A_C_i`, its representative in `W_i` lies in `K tensor QQ`;
primitivity of `K` in `W_i` forces it to lie in `K`.  Hence the other
projection is also zero.  The same argument with `K,C_i` interchanged proves
injectivity of both projections, and isotropy says that the resulting graph
map reverses the two quadratic forms.

A vector is a root of both frames exactly when it has norm two and belongs to
their intersection, which is `K`.  Finally `W_i` is the disjoint union of the
finitely many cosets of `K+C_i` indexed by `H_i`; testing norm two in those
cosets is therefore necessary and sufficient for rootlessness.  The rank
statement is linear algebra, and the final assertion follows from an
injective map between finite groups of equal order. QED.

The exact checker
[`certify_integral_rank_transfer_bridge_reglue.sage`](scripts/certify_integral_rank_transfer_bridge_reglue.sage)
applies this theorem to all 42 selected H3, Q80, NS0024, and Golay-720 edges.
Every edge has rank-15 core and rank-two bridges.  In all 42
cases both bridge discriminant groups and both glue groups are cyclic, the
old and new glue orders agree, and each glue projects isomorphically onto the
full bridge discriminant group.  Thus these clicks are exactly **rank-two
cyclic bridge replacements**, with glue orders

```text
15, 23, 31, 47, 63, 95, 119, 127, 143, 159,
191, 215, 303, 359, 799, 991, 1231, 1535, 2447, 3231.
```

The generated certificate is
[`elkies-k3-integral-rank-transfer-bridge-reglue-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-bridge-reglue-v1.json).
It records explicit Smith-coordinate glue generators and the complete
removed, surviving, and introduced norm-two vectors on every edge.  This is
an exact theorem and replay for the selected marked edges, not a completeness
claim for all primitive `U` embeddings or a proof that every rootless frame is
reachable by such edges.

### Corollary H0: theta convolution gives an inverse glue enumerator

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-THETA-CONVOLUTION 5ebbd3d242fdb3db -->

Let `K` and `C` be positive-definite even lattices.  For a discriminant class
`a in A_K` and a nonnegative rational number `nu`, define the finite coset
theta coefficient

```text
theta_K(a,nu) = #{x in K dual : x modulo K=a and x^2=nu},
```

and define `theta_C` similarly.  If `H` is an isotropic subgroup of
`A_K direct_sum A_C` and `W_H` is the corresponding even overlattice of
`K direct_sum C`, then

```text
#Phi(W_H)
 = sum over (a,b) in H
     sum over nu in QQ intersect [0,2]
       theta_K(a,nu) * theta_C(b,2-nu).
```

Only finitely many `nu` contribute.  Put

```text
rho_KC(a,b)
 = sum_nu theta_K(a,nu) * theta_C(b,2-nu).
```

Then `W_H` is rootless if and only if every element of `H` lies in the zero
support of `rho_KC`.  For graph glue, allowed root-annihilating mutations can
therefore be generated directly:

1. compute the low-norm coset theta tables of `K` once;
2. compute the rank-two table of each proposed bridge `C`;
3. convolve the two tables;
4. enumerate isotropic graph subgroups contained in the zero support;
5. construct and classify only the surviving overlattices.

#### Proof

The overlattice is the disjoint union

```text
W_H = union over (a,b) in H of (K+a) direct_sum (C+b).
```

Orthogonality makes the norm of `(x,y)` equal to `x^2+y^2`.  Counting the
vectors of total norm two in each disjoint coset gives the convolution
formula.  Every summand is a nonnegative integer, so the total vanishes
exactly when the convolution vanishes on every selected glue class. QED.

This is the exact inverted operation sought by the foundry: for a fixed core
and bridge universe it enumerates root-annihilating glue before constructing
rank-17 children.  It does not solve core generation, bound the number of
bridge classes, or guarantee a speedup.  Negative experiment H0a below shows
that computing only the least coset norm can be slower than direct root
enumeration and loses all ranking information on one held-out shell; cached
theta convolution must be benchmarked separately.

The exact implementation
[`certify_integral_rank_transfer_theta_convolution.sage`](scripts/certify_integral_rank_transfer_theta_convolution.sage)
computes and hash-locks the complete norm-at-most-two theta tables for the
four terminal H3, Q80, NS0024, and Golay-720 cores.  Starting only from each
cyclic bridge determinant, it independently enumerates every reduced positive
even binary bridge and recovers the complete fourteen-class universe.  It
then derives all 28 oriented graph multipliers from finite-form isotropy alone
and evaluates their convolutions before reading child outcomes or constructing
any rank-17 child.  The convolution count agrees with the stored independent
child-root count in every case and selects exactly the five rootless bridge
classes.  This proves the inverse enumerator on that complete fixed-core
universe; it does not yet prove a speedup or a core-generation rule.

### Theorem H0b: maximal graph glue generates the core genus

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CORE-GENERATION d0d78c49b44f55ac -->

Retain the notation of Theorem H, and suppose that `W` is a maximal graph
overlattice of `K direct_sum C`: the projection of its isotropic glue group
`H` to `A_C` is an isomorphism.  Put `B=pr_K(H)`.  Then `B` is a
nondegenerate finite quadratic submodule of `A_K`, and

```text
(B,q_B) is isomorphic to (A_C,-q_C),
A_K = B orthogonal_sum B^perp,
(A_W,q_W) is isomorphic to (B^perp,q_K restricted to B^perp).
```

Equivalently,

```text
q_K is isomorphic to q_W orthogonal_sum (-q_C),
|det K| = |det W| * |det C|.                         (H0b.1)
```

Conversely, any finite-form splitting in (H0b.1) constructs a maximal
isotropic graph in `A_K direct_sum A_C`; the corresponding even overlattice
has discriminant form `q_W`.  Thus, for a prescribed Picard-rank-19
Neron--Severi genus and a proposed positive even binary bridge `C`, the
rank-15 **core genus** is not guessed from a fibration or equation.  It is
forced by

```text
q_W = -q_NS,
q_K = q_W orthogonal_sum (-q_C),
det K = |det NS| * det C.                            (H0b.2)
```

#### Proof

The graph map identifies `A_C` anti-isometrically with `B`.  Since `q_C` is
nondegenerate, so is `q_B`, and therefore `A_K=B orthogonal_sum B^perp`.
For `z in B^perp`, the class `(z,0)` belongs to `H^perp` and maps into
`H^perp/H`.  This map is injective.  Moreover

```text
|H^perp/H| = |A_K|*|A_C|/|H|^2 = |A_K|/|A_C| = |B^perp|,
```

so it is an isomorphism of finite quadratic forms.  This proves the first
three identities, the orthogonal splitting, and the determinant identity.
Conversely, identify the `-q_C` summand of `q_K` with `-q_C`.  Its diagonal
with `q_C` is isotropic and projects isomorphically to `A_C`; quotienting its
orthogonal complement by that diagonal leaves exactly `q_W`. QED.

This is the missing algebraic inversion one level above theta convolution.
It changes core discovery from an unstructured search over changes of `U`
to enumeration of lattice classes in explicitly generated rank-15 genera.
It does **not** say that a genus has a rootless class or provide a general
fast way to enumerate and mass-close that genus.
The finite-form step is the graph-specialized consequence of Nikulin's
[even-overlattice correspondence](https://www.mathnet.ru/eng/im1677).
Completeness of a positive-definite genus enumeration still requires a
separate argument, such as neighbour closure together with the exact mass
check used by
[Conway--Sloane](https://doi.org/10.1016/0022-314X(82)90084-1).

### Corollary H0c: theta-decorated core forms classify bounded completions

For a positive even rank-15 core define its norm-two decorated discriminant
form

```text
Sigma_2(K) = (A_K, q_K, {theta_K(a,nu): a in A_K, 0 <= nu <= 2}).
```

Fix a finite list of positive even binary bridges.  Then `Sigma_2(K)`
determines, up to relabelling by decorated finite-form isomorphisms, the full
list of maximal graph completions and which of them are rootless.  More
explicitly, `K` admits a rootless completion by `C` exactly when there are a
nondegenerate summand `B subset A_K` and an anti-isometry

```text
alpha:(A_C,q_C) -> (B,-q_K restricted to B)
```

such that

```text
rho_KC(alpha(c),c)=0 for every c in A_C.             (H0c.1)
```

If a target frame form is prescribed, add `q_K|B^perp isomorphic to q_W`.
The zero class in (H0c.1) gives the cheapest necessary gate:

```text
theta_K(0,2)=0 and theta_C(0,2)=0,
```

so any core or bridge containing a root is rejected before enumerating a
single graph or isotropic divisor.  The rest follows immediately from H0b
and the nonnegative theta convolution in H0.

The resulting bounded core-first procedure is exact:

1. enumerate the declared binary bridges `C`;
2. generate the possible core genera from (H0b.2);
3. enumerate and mass-close the rank-15 classes `K` in those genera;
4. reject every `K` with `theta_K(0,2)>0`;
5. deduplicate survivors by `Sigma_2(K)` and enumerate only the
   anti-isometries satisfying (H0c.1);
6. construct rank-17 frames only for zero-support survivors.

The checker
[`certify_integral_rank_transfer_core_generation.sage`](scripts/certify_integral_rank_transfer_core_generation.sage)
verifies (H0b.1) on both sides of all 42 recorded H3, Q80, NS0024, and
Golay-720 clicks, hence on 84 maximal graph presentations.  It records four
distinct positive terminal core signatures, all of minimum four.  At their
four determinants it independently enumerates respectively `48, 8, 8, 8`
even rank-15 genera; the form generated by `q_W orthogonal_sum (-q_C)` picks
exactly one genus in every case.  Their fourteen binary bridge classes contain
exactly five zero-support classes.
As a held-out negative control, the minimum-two gate rejects all 277 primitive
rank-15 cores in the complete determinant-78 E6 source shell before graph
enumeration.  The generated record is
[`elkies-k3-integral-rank-transfer-core-generation-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-core-generation-v1.json).

The first positive core catalogue is:

| terminal core | `det W` | `det C` | forced `det K` | `min K` | binary classes | zero-support classes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Golay-720 | 720 | 23 | 16,560 | 4 | 2 | 1 |
| H3 | 948 | 47 | 44,556 | 4 | 3 | 1 |
| NS0024 | 950 | 191 | 181,450 | 4 | 7 | 2 |
| Q80 | 948 | 23 | 21,804 | 4 | 2 | 1 |

The artifact records the normalized finite quadratic form and the digest of
the complete norm-at-most-two coset theta table for every row, so these are
reusable core signatures rather than corridor names used as proxies.

This baseline procedure is exact relative to a declared finite bridge
universe and a complete core-genus enumerator.  The next theorem removes full
genus enumeration from the *screening* layer.

### Theorem H0d: reverse theta masks and the finite signature sieve

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-REVERSE-THETA eee16ce986ec0a1f -->

For a fixed binary bridge `C` and isotropic graph `H`, define its reverse mask
on core theta cells by

```text
F_(C,H) = {(a,2-nu) : (a,b) in H and theta_C(b,nu)>0}.
```

Then

```text
W_H is rootless
  iff theta_K(a,mu)=0 for every (a,mu) in F_(C,H).   (H0d.1)
```

Thus `C` and `H` prescribe all required vanishings before `K` is constructed.
The full table `Sigma_2(K)` need not be computed: a rootless decision queries
only the finite cells in `F_(C,H)` and may stop at its first occupied cell.

#### Proof

For `(a,b) in H`, Corollary H0 writes its root contribution as

```text
sum_nu theta_K(a,2-nu) * theta_C(b,nu).
```

Every factor is a nonnegative integer.  This sum vanishes exactly when each
core coefficient paired with a positive bridge coefficient vanishes.  Doing
this for every graph element is precisely (H0d.1). QED.

There is also a finite signature statement in the direction needed for core
generation.  Suppose `K` is rootless of rank `r`.  For a fixed discriminant
class `a` and fixed `0 < nu <= 2`, distinct vectors in

```text
{x in K dual : x modulo K=a and x^2=nu}
```

differ by a nonzero vector of `K`, hence have squared distance at least four.
Their pairwise inner products are therefore at most `nu-2 <= 0`.  The
orthoplex bound for vectors with pairwise nonpositive inner products gives

```text
0 <= theta_K(a,nu) <= 2r.                            (H0d.2)
```

For rank 15 every nonconstant coefficient of `Sigma_2(K)` is consequently in
`{0,...,30}`.  Since `A_K` is finite and norms in one even discriminant coset
are fixed modulo two, the set of norm-at-most-two signatures satisfying at
least one bridge mask is finite and explicitly enumerable.  More precisely,
for a declared finite bridge/graph universe `B`, every completable core lies
in the finite union of coordinate faces

```text
T(q_K,B) = union over (C,H) in B of
  {t : I(q_K)->{0,...,30} : t restricted to F_(C,H)=0},
```

after imposing `t(0,0)=1`, `t(a,nu)=t(-a,nu)`, and the quadratic-form support
congruences.  This is the requested implication

```text
zero-support requirements -> allowed Sigma_2(K) -> compatible cores.
```

The union above is deliberately an outer approximation: not every bounded
integer signature is realized by a lattice.  A much sharper signature-first
enumerator can use that the vector of coset theta series of a positive even
rank-15 lattice is a vector-valued modular form of weight `15/2` for the Weil
representation of `q_K`.  The reverse-mask vanishings are linear conditions
on its Fourier coefficients.  Computing that modular-form space, applying a
Sturm bound, integrality and nonnegativity, and then realizing only surviving
signatures is the next exact algorithmic layer.  Vector-valued modularity of
lattice theta series is standard; see, for example, the
[Weil-representation basis theorem](https://arxiv.org/abs/2407.01205).

The checker
[`certify_integral_rank_transfer_reverse_theta_masks.sage`](scripts/certify_integral_rank_transfer_reverse_theta_masks.sage)
derives the masks for all 28 oriented graphs in the terminal bridge census.
Each oriented mask has only 18--87 core cells.  After the compulsory symmetry
`theta_K(a,nu)=theta_K(-a,nu)`, they compile to an antichain of exactly 14
unoriented masks of 10--44 cells, with no containment redundancy.  This is
compared with 4,418--22,579 occupied cells in the independently computed full
core tables.  A lazy exact CVP oracle queries only 13, 25, 94, and 13 distinct
cells across all graphs for
the Golay-720, H3, NS0024, and Q80 cores respectively.  It reproduces every
accept/reject decision; a later independent full-table phase reproduces all
28 signed root counts.  Ten orientations, representing five unoriented binary
bridge classes, pass.  No rank-17 child is constructed.  On the recorded
workstation the lazy phases took 0.010--0.089 seconds per core, versus
0.74--2.68 seconds for the independent full tables; these timings are
informative, not theorem fields.  The generated record is
[`elkies-k3-integral-rank-transfer-reverse-theta-masks-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-reverse-theta-masks-v1.json).

What remains open is now narrower: construct or enumerate exactly the lattice
classes realizing the allowed vector-valued theta signatures without walking
through every rejected class of the genus.  No constrained-realization
algorithm, constrained mass formula, unbounded bridge-determinant cutoff,
uniform complexity bound, or equation-level lift is proved here.

### Theorem H0e: zero-orbit Weil compression

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-WEIL-COMPRESSION 34d2abea91a265f4 -->

The modular-form sieve proposed after Theorem H0d does **not** require the
full coefficient space `C[A_K]`.  Write the forced splitting from Theorem H0b
as

```text
A_K = A_W orthogonal-sum B,       (B,q_B) isomorphic to (A_C,-q_C).
```

After choosing this splitting, every reverse-mask class is `(0,b)`: the graph
uses precisely the distinguished copy `B` and has zero `A_W` component.  Let

```text
E_W = C[A_W]^{O(q_W)}
```

be the space of functions constant on orthogonal-group orbits.  If `Theta_K`
is the vector-valued theta series of a candidate core, then there is a
vector-valued modular form

```text
Theta_K^av in E_W tensor C[B]
```

with exactly the same Fourier coefficients as `Theta_K` in every component
`(0,b)`.  Consequently all reverse-mask constraints can be imposed in
`E_W tensor C[B]`, followed by the compulsory theta symmetry
`(a,b) ~ (-a,-b)`, without losing or creating any solution at a masked
coordinate.

#### Proof

The Weil representation of an orthogonal sum is the tensor product of the
two Weil representations.  Every `g in O(q_W)` commutes with the `S` and `T`
operators of the first factor.  Hence the Reynolds operator

```text
P_W = (1/|O(q_W)|) sum_{g in O(q_W)} g
```

commutes with the Weil representation, and
`Theta_K^av=(P_W tensor 1)Theta_K` is again a modular form of the same weight
and type.  Since zero is fixed by every orthogonal transformation,

```text
(Theta_K^av)_(0,b)
  = (1/|O(q_W)|) sum_g (Theta_K)_(g^{-1}0,b)
  = (Theta_K)_(0,b).
```

These are exactly the components read by the masks. QED.

There is a canonical possible refinement.  Let

```text
R_W = span {rho_W(g)e_0 : g in Mp_2(Z)}.
```

It is contained in `E_W`, because `e_0` is orthogonally invariant and the two
actions commute.  Orthogonal projection to `R_W tensor C[B]` also preserves
the `(0,b)` coefficients.  Thus `dim R_W`, rather than the number of
`O(q_W)`-orbits, is the sharp representation-theoretic coefficient count for
this zero-slice problem.

The exact checker
[`certify_integral_rank_transfer_weil_compression.sage`](scripts/certify_integral_rank_transfer_weil_compression.sage)
computes both quantities for the four terminal target forms.  It constructs
the complete `O(q_W)` orbit quotient and closes `e_0` under its orbit-level
Weil `S,T` matrices.  The closure rank is certified after specialization at a
prime congruent to one modulo the complete phase level.  A full-rank minor
there cannot come from the zero cyclotomic integer, so this gives an exact
characteristic-zero rank certificate.  The common nonzero scalar in `S` is
irrelevant to invariant-subspace generation.

| corridor | `|A_W|` | `|O(q_W)|` | `dim E_W=dim R_W` | `|A_K|` | after `O(q_W)` | after theta symmetry | exact reduction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Golay-720 | 720 | 96 | 72 | 16,560 | 1,656 | 864 | `115/6` |
| H3 | 948 | 8 | 240 | 44,556 | 11,280 | 5,760 | `3713/480` |
| NS0024 | 950 | 4 | 260 | 181,450 | 49,660 | 24,960 | `18145/2496` |
| Q80 | 948 | 8 | 240 | 21,804 | 5,520 | 2,880 | `1817/240` |

Thus the exact modular feasibility problem is smaller by factors between
about `7.27` and `19.17` on the observed terminal modules.  The equality
`R_W=E_W` is also a useful negative result: for these four controls there is
no further compression obtained merely by replacing the orbit quotient with
the cyclic submodule of `e_0`.

The generated certificate is
[`elkies-k3-integral-rank-transfer-weil-compression-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-weil-compression-v1.json).
This theorem removes an avoidable full-discriminant-space cost from the next
step.  It does not yet compute the compatible modular forms, impose a Sturm
bound, realize a signature by a lattice, or generate a new core.

### Corollary H0f: the linear modular mask sieve is nonselective

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MODULAR-DIMENSION-SIEVE 9622c6eb4d8522bd -->

Let `M` and `S` be the weight-`15/2` modular-form and cusp-form spaces for the
theta-symmetric zero-orbit representation of Theorem H0e.  The invariant
Riemann--Roch trace formula gives:

| corridor | coefficient dimension | `dim M` | `dim S` | largest mask | smallest certified `dim ker(S -> mask)` |
| --- | ---: | ---: | ---: | ---: | ---: |
| Golay-720 | 864 | 476 | 472 | 11 | 461 |
| H3 | 5,760 | 3,121 | 3,120 | 19 | 3,101 |
| NS0024 | 24,960 | 13,488 | 13,485 | 44 | 13,441 |
| Q80 | 2,880 | 1,563 | 1,562 | 11 | 1,551 |

In particular, every one of the fourteen reverse masks leaves a nonzero cusp
subspace, and the smallest certified kernel dimension is `461`.

#### Proof

For the dual core form

```text
-q_K = (-q_W) orthogonal-sum q_C,
```

the checker evaluates the finite-image Riemann--Roch trace formula after
averaging by

```text
O(q_W) times {+1,-1 on A_C}.
```

The `S` and `ST` traces factor into finite twisted Gauss sums on the two
discriminant modules.  The `T` exponent sum and isotropic-cusp count are
computed directly on the product of `O(q_W)`-orbit representatives and
`b ~ -b` bridge representatives.  All phases are rational.  Complex ball
evaluation at 192 bits encloses exactly one integer in every Riemann--Roch
dimension interval, certifying the displayed dimensions.  Subtracting the
isotropic-orbit count gives `dim S`.

For a mask `F` of `m` cells, coefficient extraction defines a linear map

```text
ev_F : S -> C^m.
```

Its rank is at most `m`, so rank--nullity gives

```text
dim ker(ev_F) >= dim S - m.
```

Substitution of the fourteen exact mask sizes gives the table. QED.

This closes the first proposed modular step negatively.  Modularity is still
a valid necessary condition, but the zero masks are far too small to make its
linear space rigid.  The calculation does **not** say that a rejected mask
admits a modular form with the affine normalization
`theta_K(0,0)=1`: the constant functional might interact with the mask image.
Still less does it produce integral nonnegative coefficients or a lattice.
The next potentially selective object is therefore the arithmetic theta cone

```text
{vector-valued modular forms with theta_K(0,0)=1,
 integral nonnegative coefficients, coefficient bound <=30,
 and lattice-realizable local densities},
```

not the ambient complex modular-form space.

### Theorem H0g: forced-genus mask-aware neighbour generation

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-GENERATION 9a7a1e01cb22f62e -->

Fix a target frame form `q_W`, a positive even binary bridge `C`, and a
maximal graph type.  Let `G_K` be the rank-15 genus forced by Theorem H0b:

```text
q_K = q_W orthogonal-sum (-q_C),
det K = det W * det C.
```

Choose any representative `K_0` of `G_K`.  For primes `p` not dividing
`det K`, Kneser `p`-neighbours remain in `G_K`.  Hence the following is a
sound core-generation procedure:

1. obtain `K_0` from the generated genus, without using a historical core;
2. walk through good-prime neighbours;
3. discard rootful neighbours by the zero-class theta coefficient;
4. score each rootless neighbour by the number of occupied cells in the
   reverse masks of Theorem H0d;
5. stop when one mask has zero occupied support, and only then construct its
   rank-17 completion.

Every returned lattice is a core in the required genus, and the completion
selected at step 5 is rootless.  The first assertion is the defining local
property of a good-prime neighbour; the second is exactly Theorem H0d.  This
is a constructive use of the reverse theta constraints: they act during core
generation rather than only after a genus catalogue has been built.

For the Golay-720 control this procedure is an exact bounded construction.
The generated rank-15 genus has determinant `16,560`; Sage's canonical genus
representative has 96 signed roots.  With random seed `314159`, good primes
`7,11,13,17,19`, 300 samples per parent, and a beam of 12 elite plus 8
diversity slots, the search reaches after seven neighbours a core with

```text
rank 15, determinant 16560, minimum 4, roots 0,
class-1/class-2 mask violations (3,0).
```

The full run examines 34,571 unique raw neighbours before the hit.  The
seven prime/vector pairs are retained as a short exact certificate.  Replaying
them reconstructs the same Gram matrix, verifies every neighbour condition,
and proves that this core is not isometric to the historical Golay core.  Its
class-2 order-23 graph glue produces a rank-17 lattice of determinant 720,
minimum four and no roots; an independent integral isometry test identifies
it with the declared Golay target frame.

The checker
[`generate_integral_rank_transfer_masked_core_neighbors.sage`](scripts/generate_integral_rank_transfer_masked_core_neighbors.sage)
provides both the short replay and the fixed-seed full search.  Its generated
record is
[`elkies-k3-integral-rank-transfer-masked-core-neighbors-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-masked-core-neighbors-v1.json).
This proves core generation for one bounded control and, importantly, a new
core isometry class leading to the known target.  It does not prove that the
beam closes the genus, that it will find every compatible class, or that its
expected running time improves on complete genus enumeration.  Extending the
test prospectively to the H3, NS0024, and Q80 forced genera is the next
computational question.

### Corollary H0h: prospective masked-core controls

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-CONTROLS 3cbde45fb2cb0f17 -->

The forced-genus neighbour construction extends beyond the Golay control,
but not yet uniformly.  Before scoring a core, discard every bridge with
`theta_C(0,2)>0`: its reverse mask contains `(0,0)`, which is occupied for
every core because `theta_K(0,0)=1`.  This compulsory gate removes bridge
class 1 in each of H3, NS0024, and Q80.

Starting at the canonical representative of each forced genus, exact
good-prime paths give:

| corridor | seed roots | path length | `det K` | viable mask profile | outcome |
| --- | ---: | ---: | ---: | --- | --- |
| H3 | 280 | 8 | 44,556 | `(0,2)` | new core; class-2 completion is the declared rootless target |
| NS0024 | 280 | 3 | 181,450 | `(8,4,0,10,6,10)` | new core; class-4 completion is a new rootless target-genus class |
| Q80 | 280 | 8 | 21,804 | `(2)` | rootless core, but a two-cell near miss |

Both successful cores have minimum four, no roots, automorphism group of
order two, and are not isometric to their historical corridor cores.  The H3
child has determinant 948 and is integrally isometric to the declared target
frame.  The NS0024 child has determinant 950 and the same discriminant form
as the declared target, but is not integrally isometric to it; thus the
calculus constructs a new rootless rank-17 class rather than merely finding a
second presentation of the known one.

The NS0024 path was found by capping nonzero mask counts at three while
retaining exact zero tests, then reserving beam slots for distinct occupied
support signatures.  It required 7,477 unique raw neighbours and hit in
generation three.  On the same fixed eight-generation support-diversity
trial, H3 and Q80 each examined 42,300 unique raw neighbours and missed; the
best Q80 core retained the sign-paired two-cell obstruction.  A separate H3
root-descent trial did find the certified eight-step path, but its score still
included the subsequently identified impossible rootful bridge, so that run
is construction provenance rather than evidence for the corrected uniform
ranking rule.

The exact short checker
[`certify_integral_rank_transfer_masked_core_controls.sage`](scripts/certify_integral_rank_transfer_masked_core_controls.sage)
replays all three paths, computes the uncapped viable-bridge masks, constructs
the two successful children, and checks their integral classification.  The
experimental driver is
[`search_integral_rank_transfer_masked_core_controls.sage`](scripts/search_integral_rank_transfer_masked_core_controls.sage),
and the pinned exact record is
[`elkies-k3-integral-rank-transfer-masked-core-controls-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-masked-core-controls-v1.json).

This proves two additional core constructions and one exact near miss.  It
does not prove a Q80 nonexistence result, beam completeness, or that truncated
masked support determines which support mutations are reachable.  The
failure of count-only and support-only archives on H3/Q80 shows that the next
classifier must retain more of `Sigma_2` or actual core isometry information.

### Theorem H0i: masked-witness survival and the Q80 route to the alternate frame

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-DIRECTED-Q80 80de8b6727cd3409 -->

Let `K` be integral, let `p` be a good prime, and let `ell` be an isotropic
line in `K/pK` represented by `y`.  Write its Kneser neighbour as

```text
M = {z in K : <z,y> = 0 mod p},
K_ell = M + Z*(y/p).
```

For every `x in K dual`, there is an exact survival criterion

```text
x in K_ell dual  iff  <x,y> = 0 mod p.              (H0i.1)
```

Indeed, `x` already pairs integrally with the sublattice `M` of `K`.  The
only additional generator is `y/p`, and its pairing with `x` is integral
exactly when `p` divides `<x,y>`.  Consequently, if `V_F(K)` is the finite set
of all physical dual vectors realizing occupied forbidden theta cells, then
an isotropic line satisfying

```text
<x,y> != 0 mod p for every x in V_F(K)              (H0i.2)
```

provably removes every current masked witness before the neighbour is
constructed.  This is only a removal theorem: `K_ell dual` may contain new
vectors in the same forbidden cells.
The isotropic-line parametrization and displayed construction of a
`p`-neighbour are standard; see Chenevier's
[*Statistics for Kneser p-neighbors*](https://arxiv.org/abs/2104.06846),
equations (1.2)--(1.3).  The dual-witness survival criterion above is the
elementary finite-incidence consequence needed here.

The Q80 control exhibits both sides sharply.  Its H0h near-miss core has two
occupied signed mask cells realized by four physical dual vectors.  Among
10,000 sampled isotropic lines at
`p=7,11,13,17,19,29,31,37,41,43`, exactly 8,919 satisfy (H0i.2).  Of their
neighbours, 1,397 are rootless, but every one regenerates a forbidden vector;
the new defect distribution is

```text
2:28, 4:266, 6:584, 8:413, 10:99, 12:7.
```

Thus one-step witness killing is not sufficient.  A multistep beam retaining
the lowest-defect cores in distinct integral isometry classes succeeds.  It
constructs 30,228 distinct directed neighbours through four generations;
4,643 are rootless and 116 have the minimal nonzero two-cell defect.  The
four selected transitions remove every parent witness.  The physical witness
counts change

```text
4 -> 6 -> 4 -> 4 -> 0,
```

so the first three steps are genuine defect replacement and the fourth is
annihilation.  Combined with the eight-step canonical-seed prefix of H0h,
this gives a twelve-step exact construction of a new rootless rank-15 Q80
core.  Its class-2 graph completion has rank 17, determinant 948, minimum
four, no roots, and the target local genus.  Direct comparison with both
mass-complete rootless controls identifies it exactly: it has 1,313
norm-four pairs, automorphism-group order four, is not integrally isometric
to published R17, and is integrally isometric to the alternate Q80 frame.
An explicit determinant-minus-one integral isometry is stored in the
certificate.

This also resolves an ambiguity in the earlier wording.  The
`declared_target_frame` constructed by the Q80 corridor preparation is the
published R17 control, not the alternate Q80 control.  Its failed isometry
test was therefore already pointing to the other class.  Independently of
the compressed discriminant-form key, Sage's exact genus machinery gives
the same signature `(17,0)` and local symbols for the completion and both
controls:

```text
2-adic:    1^-16:[4^1]_1
3-adic:    1^-16 3^-1
79-adic:   1^16 79^-1
```

Together with the mass-complete determinant-948 `J2` classification, this
promotes the construction from an unnamed completion to an explicit core
flow from the Q80 near-miss region around the published target to the other
rootless class, the alternate Q80 frame.

The short checker
[`certify_integral_rank_transfer_q80_defect_completion.sage`](scripts/certify_integral_rank_transfer_q80_defect_completion.sage)
verifies (H0i.1) on every stored transition, proves removal of every old
witness, recomputes every replacement shell, and constructs the final child.
The one-step and multistep discovery drivers are
[`search_integral_rank_transfer_q80_defect_neighbors.sage`](scripts/search_integral_rank_transfer_q80_defect_neighbors.sage)
and
[`search_integral_rank_transfer_q80_defect_beam.sage`](scripts/search_integral_rank_transfer_q80_defect_beam.sage).
The exact certificate is
[`elkies-k3-integral-rank-transfer-q80-defect-completion-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-q80-defect-completion-v1.json).

This is the first finite transition law in the desired constructive calculus.
It also shows why defect cardinality and occupied support are insufficient:
neighbour reachability depends on the pairings of the individual physical
witnesses with isotropic lines.  No monotone scalar descent, universal path,
neighbour-graph completeness, or speedup theorem follows.

### Theorem H0i.1: complete defect birth--death law

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-BIRTH-DEATH a755a3956c4c97cb -->

Retain the notation of Theorem H0i, and now assume that `K` is positive
definite and even, that `p` is odd and does not divide `det(K)`, and that the
chosen lift `y` of the nonzero isotropic line satisfies

```text
y^2 = 0 mod 2*p^2.
```

Put

```text
K_y dual = {x in K dual : <x,y> = 0 mod p}.
```

Then there is a disjoint coset decomposition

```text
K_ell dual
  = disjoint_union over 0 <= j < p of (K_y dual + j*y/p).       (H0i.3)
```

In particular every old survivor is in the `j=0` layer, and every vector in
a nonzero layer is new relative to `K dual`.  Its norm is

```text
(x+j*y/p)^2
  = x^2 + 2*j*<x,y>/p + j^2*y^2/p^2.                           (H0i.4)
```

#### Proof

Theorem H0i identifies

```text
K_y dual = K dual intersect K_ell dual.
```

The pairing map `K dual -> Z/p`, `x |-> <x,y>`, is onto.  Indeed, because
`p` does not divide `det(K)`, the pairing on `K/pK` is nondegenerate, and
`y mod p` is nonzero.  Hence `K_y dual` has index `p` in `K dual`.

The vector `y/p` belongs to `K_ell dual`: it pairs integrally with `M` by
the definition of `M`, and its pairing with the remaining generator `y/p`
is integral because `y^2` is divisible by `p^2`.  Moreover `p*(y/p)=y`
belongs to `K_y dual`, while `y/p` does not belong to `K dual`; otherwise
`y mod p` would lie in the radical of `K/pK`.  Thus `y/p` has exact order
`p` modulo `K_y dual`.

Finally `K` and its `p`-neighbour have the same determinant.  Their duals
therefore have the same covolume, so the index of `K_y dual` in
`K_ell dual` is also `p`.  The `p` displayed cosets generated by `y/p` are
distinct and exhaust the child dual.  Expanding the square gives (H0i.4).
QED.

There is a canonical prime-to-`p` identification of discriminant groups

```text
iota: A_(K_ell) -> A_K,
iota(v mod K_ell) = p^(-1)*(p*v mod K),                         (H0i.5)
```

where multiplication by `p` is invertible on `A_K`.  Formula (H0i.3) makes
this concrete: if `v=x+j*y/p`, then `iota(v)=[x]`.  Thus the layer shift does
not change the discriminant label under `iota`.

This gives an exact finite transition operator, but its input must retain the
physical lattice and the line, not only the abstract counted signature
`Sigma_2(K)`.  For every class `a in A_K` and `0 <= mu <= 2`,

```text
theta_(K_ell)(iota^(-1)(a),mu)
  = sum over 0 <= j < p of
      #{x in K_y dual : [x]=a and (x+j*y/p)^2=mu}.              (H0i.6)
```

Every set on the right is finite by positive definiteness and is an affine
CVP query.  An implementation tailored to a reverse mask can work in `M`
instead.  Choose a representative `r` of `a` and one `k_0 in K` satisfying

```text
<r+k_0,y> = 0 mod p.
```

The vectors in the requested child cell are the disjoint union

```text
disjoint_union over 0 <= j < p of
  {v in M+r+k_0+j*y/p : v^2=mu}.                               (H0i.7)
```

Consequently a line gives a zero-defect child exactly when both conditions
hold:

1. it is nonorthogonal modulo `p` to every current physical forbidden
   witness, equivalently every forbidden `j=0` query is empty;
2. every forbidden query in the layers `j=1,...,p-1` is empty.

The checker
[`certify_integral_rank_transfer_q80_defect_birth_death.sage`](scripts/certify_integral_rank_transfer_q80_defect_birth_death.sage)
performs these affine queries before constructing each of the four stored Q80
children.  It also evaluates the full sum (H0i.6): the four child profiles
contain respectively `10,219`, `10,201`, `10,287`, and `10,121` dual vectors
through norm two in `5,377`, `5,435`, `5,397`, and `5,485` occupied theta
cells.  Independent child-dual enumeration then gives exactly the same full
profiles and the same physical forbidden vectors.  The predicted
witness-count regression is

```text
4 -> 6 -> 4 -> 4 -> 0,
```

and the first three replacements occur entirely in nonzero affine layers.
The exact record is
[`elkies-k3-integral-rank-transfer-q80-defect-birth-death-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-q80-defect-birth-death-v1.json).

This proves the transition law and the strong zero-defect criterion.  It does
not make `Sigma_2(K)` without physical representatives a complete state, and
it does not by itself prove that separate layer queries are faster than
forming a child Gram matrix and querying its mask directly.  Runtime is an
implementation question, not a consequence of (H0i.3).

### Corollary H0j: the NS0024 completed-core path realizes ranks `4,12,12,17`

<!-- status-consumer: EC-K3-NS0024-COMPLETED-CORE-RANK-TRANSFER 16b64051fb648d66 -->

Let `S` be the determinant-950 Picard-rank-19 Neron--Severi lattice
`NS0024`.  Complete the canonical rank-15 core and its three exact
good-prime Kneser neighbours with the order-191 binary bridge class 4.  In
the successive core bases the integral graph multipliers may be chosen as

```text
59, 50, 76, 83
```

(with sign-paired alternatives `132,141,115,108`).  The exact core and
completed-frame profiles are

| stage | incoming core prime | core root rank/count | completed roots | completed root rank | MW rank |
| --- | ---: | ---: | --- | ---: | ---: |
| 0 | -- | `13/280` | `D5+E8` | 13 | 4 |
| 1 | 17 | `5/12` | `3A1+A2` | 5 | 12 |
| 2 | 13 | `0/0` | `3A1+A2` | 5 | 12 |
| 3 | 7 | `0/0` | rootless | 0 | 17 |

Every completed frame has determinant 950 and discriminant form
anti-isometric to `q_S`.  Theorem H2 therefore realizes it as the positive
frame of a primitive hyperbolic-plane embedding in `S`, at `O(S)/J2` level.
After Weyl reduction, Theorem C gives the corresponding Jacobian fibration.
Theorem A then gives

```text
rank MW = 19 - 2 - rank R = 17 - rank R,
```

so the Mordell--Weil ranks are exactly

```text
4 -> 12 -> 12 -> 17,
```

with successive mutations `+8,0,+5`.

This path separates the two mechanisms in the glue calculus.  At stage 2
the core is already rootless, but the completed frame still has twelve
roots.  Since roots of the core embed as the zero glue coset, all twelve
must lie in nonzero graph-glue cosets.  The rank plateau is therefore exact:
annihilating the remaining core roots did not increase MW rank because roots
were present in the completion glue.  The final step keeps the core rootless
while removing those glue-coset roots, and hence creates the last five MW
directions.  The physical-witness survival criterion of Theorem H0i is the
local transition law controlling such removals; Corollaries H0 and H0c turn
the surviving decorated cosets into the completed root system, and Theorem A
turns its rank into MW rank.

The checker
[`certify_ns0024_new_rootless_source_route.sage`](scripts/certify_ns0024_new_rootless_source_route.sage)
replays all four completions, verifies the core and completed root profiles,
and checks every determinant and discriminant form.  Its exact record is
[`elkies-k3-ns0024-new-rootless-source-route-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-new-rootless-source-route-v1.json).

The primes `17,13,7` are core Kneser-neighbour primes, not elliptic-neighbour
degrees.  This corollary proves the ranks and existence of the four
fibrations over a complex K3 with Neron--Severi lattice `S`; it does not give
a marked elliptic-neighbour corridor, equations, rational maps, a field of
definition, or a universal monotone rank law.

### Theorem H0k: metric physical witnesses transfer the full root system

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-ROOT-SYSTEM-SIGNATURE d32b35b66a35627c -->

Let `K,C` be positive-definite even integral lattices, let

```text
H subset A_K direct_sum A_C
```

be isotropic, and let `W_H` be the corresponding overlattice of `K+C`.  Define
the finite set of physical completion witnesses

```text
Omega_2(K,C,H)
 = {(k,c) in (K dual)+(C dual) :
      (k mod K,c mod C) in H and k^2+c^2=2}.
```

For two witnesses `r=(k,c)` and `r'=(k',c')`, put

```text
g(r,r')=<k,k'>+<c,c'>.                              (H0k.1)
```

Then `Omega_2(K,C,H)=Phi(W_H)`, and its complete matrix `g` determines the
signed norm-two root system as an abstract metric set.  In particular it
determines:

1. the complete labelled root graph;
2. the orthogonal ADE decomposition and root rank;
3. the root lattice `R=<Phi(W_H)>`, its component and total discriminants,
   and its discriminant form;
4. the finite list of *possible* even root-overlattice, hence torsion-glue,
   subgroups, namely the isotropic subgroups of `A_R` which pass the relevant
   no-new-root test.

If the signature is **marked**, meaning that it also retains the coordinates
of every physical witness in one integral basis of `W_H`, then it additionally
determines

```text
Rbar=(R tensor QQ) intersect W_H,
H_root=Rbar/R,
```

including the invariant factors and exact index.  For a Jacobian elliptic
fibration with `NS=U direct_sum W_H(-1)`, this quotient is exactly
`MW_tors`.  Thus a marked metric physical-witness signature determines the
primitive closure and actual torsion contribution, while the unmarked metric
signature determines only the possible torsion glue.  Rootlessness is the
case `Omega_2=empty`.

#### Proof

Nikulin's overlattice correspondence writes `W_H` as the disjoint union of
the cosets `(K+k)+(C+c)` selected by `H`.  Orthogonality gives

```text
(k+c)^2=k^2+c^2,
```

so the displayed witness set is literally the complete norm-two set of
`W_H`; polarization gives (H0k.1).  In an integral lattice the reflection

```text
x |-> x-<x,r>r
```

in a norm-two vector preserves the lattice and permutes its norm-two vectors.
Hence `Phi(W_H)` is a reduced simply-laced crystallographic root system.
The connected components of its nonorthogonality graph are its irreducible
components, so the finite ADE classification recovers their types.  The rank
is the rank of `g`; choosing any integral basis of the span of the roots gives
the root Gram matrix, its determinant and discriminant form.  Nikulin's same
correspondence classifies its even finite-index extensions by isotropic
subgroups of `A_R`, with the norm-two witness test deciding which extensions
add roots.

For the marked assertion, place the root-coordinate rows in a matrix in the
fixed `W_H` basis.  Smith saturation computes their primitive row span
`Rbar`, and Smith normal form computes `Rbar/R`.  Finally Shioda's canonical
identification

```text
MW = NS/(U+R) = W_H/R
```

shows that the torsion subgroup is the torsion of `W_H/R`, namely `Rbar/R`.
QED.

The word **marked** cannot be dropped from the primitive-closure conclusion.
The abstract root lattices `A1^24` and the Niemeier lattice `N(24A1)` have
the same complete metric root system `24A1`.  In the first ambient lattice
the root lattice is already primitive; in `N(24A1)` its primitive closure is
the whole Niemeier lattice and the quotient is the binary Golay code
`(Z/2)^12`.  Pairwise root products therefore cannot recover ambient
saturation or exact torsion without the witness embedding.  This is also why
the foundry state must retain physical coordinates rather than only a
canonical ADE label.

The exact checker
[`certify_integral_rank_transfer_root_system_signature.sage`](scripts/certify_integral_rank_transfer_root_system_signature.sage)
constructs the marked signature on all four NS0024 completed-core stages.  It
recovers

```text
D5+E8 -> 3A1+A2 -> 3A1+A2 -> rootless,
```

with root ranks `13,5,5,0`, root discriminants `4,24,24,1`, primitive roots,
and trivial torsion.  At the third stage the core is rootless while the twelve
completion roots occupy five nonzero order-191 graph-glue labels, so the
pairwise metric reconstructs root gluing which zero support alone cannot see.
The same classifier independently recognizes the existing Q80 route controls
`4A1` and `A1`, of root discriminants `16` and `2`.  The generated record,
including every physical `k+c` root line and every pairwise inner product, is
[`elkies-k3-integral-rank-transfer-root-system-signature-v1.json`](../artifacts/generated-results/elkies-k3-integral-rank-transfer-root-system-signature-v1.json).

This theorem supplies an exact acceptance and targeting signature.  It does
not prove that a prescribed ADE signature is realized in a forced core genus,
that a neighbour path reaches it, or that its fibration has an equation or a
specified field of definition.

### Negative experiment H0: the orthogonal split is not a useful predictor

For any proposed replacement, `K+C_new` is a sublattice of `W_new`.
Consequently

```text
rank Phi(K) + rank Phi(C_new) <= rank Phi(W_new)
```

and the analogous inequality holds for signed root counts.  These are exact,
no-false-negative screens for a specified child root budget.  They are not,
however, selective on the preserved H3 first-hit data.  The exact benchmark
[`benchmark_integral_rank_transfer_bridge_predictor.sage`](scripts/benchmark_integral_rank_transfer_bridge_predictor.sage)
tests 2,892 historical candidates through five first hits.  The core-only
bound rejects zero; adding the rank-two bridge rejects 178, leaving 2,714 and
giving only a `2892/2714 = 1.066...` projected reduction in full child
classifications.  The terminal rootless q6 window retains all 1,247
candidates.

This is a bounded retrospective counterexample to the proposed **screening
heuristic**, not to Theorem H.  It shows that the missing information is the
norm-two profile of the nonzero graph-glue cosets.  No new construction
algorithm is claimed until a predeclared coset score succeeds on an untouched
shell and is demonstrably cheaper than direct child-root enumeration.

There is nevertheless exploratory fixed-core evidence in the rank-two data.
For the four observed rootless terminal cores, the same checker exhausts all
14 admissible binary bridge classes of the observed prime determinants and
their compatible oriented graph labels.  Five classes are rootless.  The
rank-two-only rule "retain maximum bridge minimum" keeps five classes, four
rootless: precision rises from `5/14` to `4/5`, a `56/25 = 2.24` enrichment,
with `4/5` rootless recall and a projected `14/5 = 2.8` reduction in full
classifications.  This selected-core census does not repair the prospective
gap: both the core and determinant were learned from successful edges.

### Negative experiment H0a: bridge minimum is constant on an untouched shell

<!-- status-consumer: EC-K3-E6-DET78-PROSPECTIVE-BRIDGE-NEGATIVE d23a0abd146c2ed9 -->

The predeclared score from H0 was next tested prospectively on the complete
zero-neutral old-degree-two shell of the determinant-78 E6
`2E6+A1/MW4` frame.  The core-generation rule starts only from that source
frame; the successful H3, Q80, NS0024, and Golay corridors are excluded.
For each primitive child, the score

```text
min { minimum of a nonzero graph-glue coset of K+C_new }
```

is computed before enumerating child roots.  The mass-complete 1,549-class
E6 catalogue is consulted only afterward to label the outcome.

The complete shell has 280 dominant source-Weyl classes, of which 277 are
primitive.  Every one of those 277 candidates has score exactly two, while
their child root ranks range from 12 through 15 and they occupy 31 `J2`
classes.  Therefore this scalar score has zero ranking power on the declared
untouched shell.  On the recorded workstation the exact coset scoring took
51.01 seconds, versus 0.48 seconds for direct norm-two enumeration and
root-rank classification, so it was about `106x` slower.  The separate J2
isometry lookup belongs to truth-set evaluation and is excluded from that
comparison; these wall-clock timings are not theorem fields.

This exact negative control does not refute Theorem H or the decorated
root-profile calculus.  It identifies the missing datum more sharply:
prospective inversion must retain the distribution of norm-two vectors among
the candidate glue classes (the function `rho_K`, or a proved bound on it),
not merely the least lattice norm.  Because the determinant-78 genus is
globally rootful, the experiment cannot measure rootless recall; a positive
prospective gate still needs an untouched mass-complete genus containing a
rootless class.  The replay is
[`benchmark_e6_det78_prospective_bridge_predictor.sage`](scripts/benchmark_e6_det78_prospective_bridge_predictor.sage).

### Corollary H1: a good-prime neighbour is a finite hyperbolic-line swap

Let `W_0,W_1` be distinct even `p`-neighbours in one rational quadratic
space, where `p` is odd and does not divide `det(W_i)`, and put
`K=W_0 intersect W_1`.  Then

```text
W_i/K = H_i is cyclic of order p,
(A_K)_p is isomorphic to (Z/p)^2,
H_0 and H_1 are distinct isotropic lines.
```

The `p`-primary discriminant form of `K` is split, and the neighbour move is
exactly the line swap `H_0 -> H_1`; all prime-to-`p` glue is unchanged.  If

```text
rho_K(h)=number of norm-two vectors in K dual with residue h,
```

then surviving roots lie in `h=0`, removed roots lie on
`H_0 minus {0}`, and introduced roots lie on `H_1 minus {0}`.  Thus a
good-prime root-annihilating click is decided by the finite decorated form

```text
((A_K)_p, q_K, H_0, H_1, rho_K restricted to H_0 union H_1).
```

More generally, for any two isotropic glue subgroups in `A_K`,

```text
#Phi(L_H) = sum over h in H of rho_K(h),

#Phi(L_H1)-#Phi(L_H0)
  = sum over h in H1 minus H0 of rho_K(h)
    - sum over h in H0 minus H1 of rho_K(h).
```

#### Proof

The defining intersection of distinct `p`-neighbours has index `p` in each,
so the two quotient groups are distinct order-`p` isotropic subgroups of
`A_K`.  Since `p` is prime to the neighbour determinants, the determinant
formula gives `|(A_K)_p|=p^2`.  A nondegenerate quadratic plane containing an
isotropic line is split.  Nikulin's overlattice correspondence reconstructs
the two neighbours from the two lines.  Decomposing the norm-two vectors in
`K dual` by their residue classes proves the root assertions. QED.

The hypotheses matter.  For `p>2` the move swaps all `p-1` nonzero cosets on
one line, not one coset.  At `p=2` or `p` dividing the determinant, the local
module and isotropy rules require a separate calculation.

### Theorem H2: rank-19 frame-genus realization

Let `S` be an even lattice of signature `(1,r-1)` with

```text
r >= length(A_S)+2.
```

Then the `O(S)`-orbits of embedded hyperbolic planes `U` are in bijection
with the integral isometry classes of positive-definite even lattices `P` of
rank `r-2` satisfying

```text
q_P is isomorphic to -q_S.
```

The map sends `U` to `P=U^perp(-1)`.  For the Neron--Severi lattice of a
complex Picard-rank-19 K3, the length hypothesis is automatic because
`A_NS` is anti-isometric to the discriminant group of the rank-three
transcendental lattice.  Therefore a rootless rank-17 lattice in the required
finite-form genus is equivalent, at the lattice/fibration level, to a
rootless Jacobian fibration and hence to Mordell--Weil rank 17.

#### Proof

An embedded `U` splits off integrally because it is unimodular, giving the
stated complement and finite form.  Conversely `U+P(-1)` has the same
signature and discriminant form as `S`.  Nikulin's uniqueness theorem for
indefinite even lattices under the displayed length bound makes it isometric
to `S`, producing the embedding.  An isometry between two complements
extends by the identity on `U`, proving injectivity on `O(S)`-orbits.

For a complex K3, discriminant duality with the transcendental lattice gives
the length bound.  Weyl reduction moves a primitive isotropic generator to a
nef class; its companion in `U` supplies a degree-one `(-2)` class, so Theorem
C produces a Jacobian fibration.  If `P` is rootless, Shioda--Tate gives MW
rank 17. QED.

This is an `O(NS)`/lattice-fibration classification.  It is not a
classification modulo `Aut(X)` and does not remove the need for chamber or
equation certificates.  It does, however, make the foundry's target exact:
enumerating rootless classes in the one prescribed finite-form genus is
equivalent to enumerating rank-17 frame classes at J2 level.

### Corollary H2a: rank-three Hodge rigidity and a finite rootless J1 bound

<!-- status-consumer: EC-K3-H3-ROOTLESS-J1-UNIFORM-BOUND b71330a75ad2c9ad -->

Let `X` be a complex projective K3 surface of Picard rank 19.  Then

```text
Isom(T_X)^Hodge = {+identity,-identity}.
```

Consequently the Braun--Kimura--Watari uniform multiplicity bound for every
fixed frame class `[P]` in `J2(X)` is

```text
# inverse_image_J1([P])
    <= | image(Isom(T_X)^Hodge) backslash Isom(q_T) |.
```

For the pinned determinant-948 H3/R17 surface, the cyclic discriminant form
has

```text
Isom(q_T) = {1,157,317,473,475,631,791,947} modulo 948,
image(Isom(T_X)^Hodge) = {1,947}.
```

Thus each rootless `J2` frame has `J1` multiplicity at most four.  The
mass-complete rootless `J2` classification has exactly two frame classes, so
the number `n_rootless_J1` of rootless Jacobian fibrations modulo surface
automorphisms satisfies the unconditional finite interval

```text
2 <= n_rootless_J1 <= 8.
```

#### Proof

The real transcendental space has signature `(2,1)`.  A Hodge isometry `g`
acts on `H^(2,0)` by an eigenvalue `zeta`, on `H^(0,2)` by its conjugate, and
on the real orthogonal line by `epsilon=+1` or `-1`.  All three eigenvalues
have absolute value one.  Since `g` is integral, they are roots of unity.  If
`zeta` is nonreal, the rational `epsilon`-eigenspace is a nonzero rational
`(1,1)` subspace of `T_X`, contradicting
`T_X intersect NS(X)_QQ=0`.  If `zeta=+1` or `-1` and `epsilon` has the
opposite sign, the same contradiction applies.  Hence `g` is scalar and is
`+identity` or `-identity`; both occur.

Proposition C' of Braun--Kimura--Watari bounds the `J1` multiplicity of every
`J2` class by the displayed discriminant-form coset count.  Direct exact
enumeration of units `u modulo 948` satisfying

```text
q(u*g)=q(g) modulo 2*ZZ
```

gives the eight listed units for either pinned generator of the form.  The
image of `-identity` is `947`, so the quotient has four cosets.  The complete
rootless `J2` classification supplies two distinct realized frame classes;
the map from `J1` to `J2` is surjective on them.  Summing the two upper bounds
gives eight, while distinct `J2` classes give the lower bound two. QED.

The finite calculation and both input hashes are locked by
[`certify_rootless_j1_uniform_bound.py`](scripts/certify_rootless_j1_uniform_bound.py)
and
[`elkies-k3-rootless-j1-uniform-bound-v1.json`](../artifacts/generated-results/elkies-k3-rootless-j1-uniform-bound-v1.json).
The external multiplicity theorem is
[Braun--Kimura--Watari, Proposition C'](https://arxiv.org/abs/1312.4421).

This is a strict advance from an unbounded `J1` frontier to a finite list, but
it is not the exact `J1` classification.  In particular, the recorded images
of the two frame automorphism groups in the discriminant form do not by
themselves justify replacing the uniform quotient by a smaller
frame-dependent quotient: that requires ample-cone stabilizer control or an
equivalent surface-automorphism computation.

### Theorem H3: conditional one-click root annihilation

Let `G` be a genus of positive-definite integral lattices of rank greater
than two, assume `G` is one spinor genus, and suppose it contains a rootless
class `P_star`.  For every `P` in `G` and every sufficiently large prime
`p` not dividing `det(P)`, there is a `p`-neighbour of `P` integrally
isometric to `P_star`.

Consequently, after the finite hypotheses “one spinor genus” and “a certified
rootless class” have been checked, the discriminant-form calculus has an
honest root-annihilating one-click existence theorem.  In the common
intersection `K=P intersect P_star`, Corollary H1 identifies the click with a
split finite quadratic-plane line swap: all old roots outside `K` occur on
the removed line, while the new line has no norm-two vector.

#### Proof

Chenevier's equidistribution theorem for Kneser neighbours says that, in one
spinor genus, the normalized number of `p`-neighbours of `P` isometric to a
fixed class tends as `p` grows to the positive mass of that class.  Taking
the fixed class to be `P_star` makes the count positive for every sufficiently
large good prime.  Corollary H1 gives the line-swap and root description.
QED.

This theorem is non-effective in the threshold prime and its one-spinor-genus
hypothesis has not been certified for every frame genus in this repository.
A concrete result still requires an explicit neighbour witness or a separate
spinor-genus/mass computation.  Combined with Theorem H2, however, it proves
that any such rootless neighbour is realized by another `U` on the same
Picard-rank-19 K3 at the lattice/fibration level.

### Theorem H4: mass-complete decorated glue calculus

Let `G` be a genus of positive-definite even lattices of rank at least three,
let `W_0` lie in `G`, and fix a finite set of good odd primes.  Run an exact
Kneser-neighbour breadth-first search from `W_0`, deduplicating vertices by
integral isometry and computing `|O(W)|` exactly.  If the visited vertices
satisfy

```text
sum over visited [W] of 1/|O(W)| = mass(G),
```

then the list is the complete genus.  In particular, `G` contains a rootless
lattice if and only if the list contains one, and every such lattice is
reached from `W_0` by a finite sequence of the discriminant-line mutations in
Corollary H1.  Decorating every edge by

```text
(A_K,q_K,H_old,H_new,rho_K on H_old union H_new)
```

makes the change in its complete root set an exact finite calculation.
Combined with Theorem H2, this is a terminating constructive calculus for
rootless MW-rank-17 frames at `O(NS)`/J2 level whenever the mass closes.

#### Proof

Every enumerated neighbour remains in `G`.  The Minkowski--Siegel mass is the
sum of `1/|O(W)|` over all integral isometry classes in `G`, with every term
positive.  Equality therefore leaves no missing class.  The search tree gives
a neighbour path from `W_0` to each listed vertex, and Corollary H1 identifies
each edge and its root transfer with the displayed finite decorated form.
Theorem H2 converts the frame classes to `O(NS)`-orbits of hyperbolic-plane
embeddings. QED.

The mass equality is a certificate, not a stopping heuristic.  Standard
Kneser theory supplies finite prime sets which connect the relevant class set,
but a particular implementation must still prove connectivity or close the
mass.  Nor can the decoration be reduced to the bare finite quadratic form:
the Leech lattice and the 23 rooted Niemeier lattices all have trivial
discriminant form, while only the Leech lattice is rootless.  Coset minima or
the equivalent function `rho_K` are indispensable state.

### Theorem H5: the determinant-78 E6 frame genus is globally rootful

Let `P` be the positive rank-17 frame genus of the saturated determinant-78
E6 `2+2` Neron--Severi lattice.  Every lattice in this genus contains a root.
Consequently this Picard-rank-19 lattice has no rootless MW-rank-17 fibration
at `O(NS)`/J2 level.

#### Proof

The rank-seven Nishiyama auxiliary `K` has determinant 78 and contains the
primitive root lattice

```text
S=A3 direct_sum A2 direct_sum A1.
```

In the standard basis its final generator `v` has square four and pairs only
as `(-1,0,0)` with the `A3` factor.  Hence the projection of `v` to `S` has
square `3/4`, leaving the exact orthogonal norm budget

```text
v_perp^2 = 4-3/4 = 13/4.
```

The discriminant form of `K` is anti-isometric to that of `P`.  If a rootless
`W` lay in the genus of `P`, their full graph glue would give a positive even
unimodular rank-24 lattice `N` containing `K` and `W=K^perp` primitively.
Because `K` contains roots, `N` is one of the 23 rooted Niemeier lattices.

The exact residual-Weyl enumeration described below covers every primitive
embedding of `S` in those 23 lattices.  Across its 1,591 root-anchor
representatives, the residual root system

```text
R(N) intersect S^perp
```

always has rank at least 14.  Apply its Weyl group, which fixes `S`, to move
`v` into the dominant chamber.  If `W` were rootless, `v` could be orthogonal
to no residual root, so all Dynkin labels on the residual simple roots would
be positive integers.  For every simply-laced Cartan matrix `C`, `C^-1` is
entrywise nonnegative and every diagonal entry is at least `1/2`.  Therefore
the squared norm of the residual projection is at least

```text
(1,...,1) C^-1 (1,...,1)^t >= rank(C)/2 >= 7,
```

contradicting the available norm `13/4`. QED.

The focused mode of
[`classify_e6_rank4_det78_niemeier_frames.sage`](scripts/classify_e6_rank4_det78_niemeier_frames.sage)
checks the opposite discriminant forms, enumerates 5,325 fixed-root A3
subsystems and the complete 1,591 primitive `A3+A2+A1` residual-Weyl anchor
cover, and obtains residual-rank distribution

```text
14:24, 15:721, 16:619, 17:198, 18:29.
```

The replay artifact is
[`elkies-k3-e6-rank4-det78-rootless-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-det78-rootless-obstruction-v1.json).
This closes the earlier degree-at-most-four search globally at J2 level; it
does not classify all frame isometry classes or construct equations.

## 7B. Integral character glue

### Theorem I: involution eigensublattices are joined by a 2-primary graph

Let `L` be an even integral lattice with an isometric involution `sigma`, and
let

```text
L+ = L intersect ker(sigma-1),
L- = L intersect ker(sigma+1).
```

Then `L+` and `L-` are primitive and orthogonal, and

```text
H=L/(L+ direct_sum L-)
```

is killed by two.  Under the even-overlattice correspondence, `H` embeds in
both `A_L+ [2]` and `A_L- [2]` and is the graph of an anti-isometry between
2-elementary subgroups of the two eigendiscriminant forms.  Conversely, any
such isotropic graph defines an even involution-stable overlattice of
`L+ direct_sum L-`.

#### Proof

If a nonzero multiple of `x` lies in either eigensublattice, torsion-freeness
of `L` shows that `x` has the same eigenvalue, proving primitivity.  For
`x` in `L+` and `y` in `L-`, invariance of the pairing gives

```text
(x,y)=(sigma*x,sigma*y)=-(x,y),
```

so the eigensublattices are orthogonal.  Every `z` in `L` satisfies

```text
2z=(z+sigma*z)+(z-sigma*z),
```

which proves the exponent-two assertion.  Nikulin's correspondence realizes
`L` by an isotropic subgroup of `A_L+ + A_L-`.  Primitivity makes both
projections injective, exactly as in Theorem H, so this subgroup is an
anti-isometry graph.  Conversely an isotropic graph gives an even
overlattice.  On 2-torsion the actions `+1` and `-1` coincide, hence the graph
is stable under the induced involution and the involution extends. QED.

The checker
[`certify_integral_character_glue_calculus.sage`](scripts/certify_integral_character_glue_calculus.sage)
exhausts the graph possibilities in the two E6 examples after multiplying
the rational MW height pairing by 12 to obtain an even integral model.

- For E6 `2+1`, the actual index is one.  There are three possible nonzero
  order-two graphs, all in one integral isometry class; that alternative
  lowers the scaled minimum from eight to six and is a concrete control that
  adding glue can create shorter vectors but cannot remove old ones.
- For E6 `2+2`, full index-four saturation forces a graph isomorphism between
  two copies of `F_2^2`.  All six graph isomorphisms give one integral
  isometry class, of determinant 89,856, minimum 16, six signed minimal
  vectors, and automorphism-group order eight.  Thus the observed half-sums
  exhaust the full graph choices up to integral isometry.

The exact finite enumeration is stored in
[`elkies-k3-integral-character-glue-calculus-v1.json`](../artifacts/generated-results/elkies-k3-integral-character-glue-calculus-v1.json).
It classifies the integral character glue; the input equation certificates
remain responsible for geometric existence and descent of the sections.

<!-- status-consumer: EC-K3-R17-RANK28-INTEGRAL-CHARACTER-GLUE 617f1838d8581fcd -->

### Corollary I1: the eleven fitted rank-28 lifts admit one character-glue type

On each of the eleven exact genus-one double covers fitted through the public
rank-28 specialization, let `R_i` be the lifted section, let `sigma_i` be the
deck involution, and put

```text
tau=R_i+sigma_i(R_i),       T_i=R_i-sigma_i(R_i).
```

All eleven covers have the same trace `tau=-P2-P5`.  In the pulled-back
Mordell--Weil height lattice,

```text
tau^2=T_i^2=16,       tau.T_i=0,       2R_i=tau+T_i.
```

Hence every direction is obtained from the pure rank-two character lattice

```text
<16>_+ direct_sum <16>_-
```

by adjoining the diagonal half-sum.  The glue subgroup has order two and is
represented by `(8,8)` in `Z/16+Z/16` in the Smith coordinates induced by
`(tau,T_i)`; the saturated carrier is isometric to `<8>+<8>`.  The eleven
quartic branch polynomials define eleven distinct squareclasses, so this is
one repeated integral pattern on eleven distinct characters, not one cover
with an eleven-dimensional anti-invariant space.  The pattern is relative to
the chosen trace pencil and is not an intrinsic invariant of an isolated
fibre point.

#### Proof

The trace has norm eight in `R17`, and height pairings double under the
degree-two base change, giving `tau^2=16`.  The equation certificate gives
`T_i^2=16`; invariant and anti-invariant characters are orthogonal.  The
half-sum identity follows from the definitions.  Theorem I then identifies
the order-two graph glue.  Direct basis change gives Gram matrices

```text
[16 8]          [8 0]
[ 8 8]    and   [0 8]
```

in bases `(tau,R_i)` and `(R_i,sigma_i(R_i))`.  Finally, two squarefree
irreducible quartics represent the same class in `QQ(t)^*/QQ(t)^{*2}` only
if they are proportional; primitive normalization makes the eleven stored
polynomials pairwise distinct. QED.

The exact derived replay is
[`certify_rank28_integral_character_glue.py`](scripts/certify_rank28_integral_character_glue.py),
with certificate
[`elkies-k3-r17-rank28-integral-character-glue-v1.json`](../artifacts/generated-results/elkies-k3-r17-rank28-integral-character-glue-v1.json).
It inherits geometric existence and height 16 from the genus-one-bisection
certificate.  It does not find a new specialization, make the eleven covers
split together away from the fitted fibre, or prove rank 32.

### Corollary I2: the norm-twelve `0x103b2` cover has exact visible glue

For the genus-one double cover attached to the R17 trace class `0x103b2`, let
`R` be the displayed lift, `sigma` the deck involution, and

```text
tau=R+sigma(R),       T=R-sigma(R).
```

Then the invariant and anti-invariant lattices in their displayed rational
span are

```text
L+=R17(2),       L-=<16>,       tau^2=24,       2R=tau+T.
```

The full integral saturation in this span is obtained by adjoining the one
order-two graph class `(tau/2,T/2)`.  It has rank 18, determinant
`497025024`, minimum eight, and no roots.  Its pure character lattice has
determinant `1988100096`; its complete Smith group and finite quadratic form
are stored in the certificate.  At `t=1/25` the cover splits and the resulting
point is independent of the specialized generic MW17 subgroup, proving rank
at least 18 on that fibre.

#### Proof

Heights double under the quadratic pullback, so the norm-twelve trace has
square 24.  The exact bisection intersection computation gives `T^2=16`, and
the two characters are orthogonal.  Thus `R=(tau+T)/2` is an isotropic
order-two discriminant class.  Theorem I and primitivity of both integral
eigensublattices show that this is the full saturation inside the declared
rational span.  The determinant follows by division by the square of the
index; exact shortest-vector enumeration gives minimum eight and no roots.
The square identity at `t=1/25` and the independent finite-quotient reductions
give the final rank lower bound. QED.

The replay is
[`certify_r17_norm12_103b2_mw_glue.sage`](scripts/certify_r17_norm12_103b2_mw_glue.sage),
with artifact
[`elkies-k3-r17-norm12-103b2-mw-glue-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-103b2-mw-glue-v1.json).
The exact claim is saturation of the displayed cover-level rational span and
specialized rank at least 18.  Exact eclib saturation at each prime dividing
the index bound is run in an isolated process and proves that the displayed
specialized rank-18 subgroup is primitive.  An upper bound for the
specialized rank is not claimed.

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

- V. V. Nikulin,
  [*Integral symmetric bilinear forms and some of their applications*](https://www.mathnet.ru/eng/im1677),
  especially Proposition 1.4.1 for even overlattices and the indefinite
  uniqueness results used in Theorem H2.
- J. H. Conway and N. J. A. Sloane,
  [*Low-dimensional lattices. I. Quadratic forms of small determinant*](https://doi.org/10.1016/0022-314X(82)90084-1),
  for the Niemeier root systems and glue-word table used in the `24A1`
  marked-embedding counterexample of Theorem H0k.
- G. Chenevier,
  [*Statistics for Kneser p-neighbors*](https://arxiv.org/abs/2104.06846),
  for good-prime neighbour parametrization and the equidistribution theorem
  used in Theorem H3.
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
- D. Kubert,
  [*Universal bounds on the torsion of elliptic curves*](https://doi.org/10.1112/plms/s3-33.2.193),
  for the classical marked-point Tate normal form.
- M. Cvetic, D. Klevers, and H. Piragua,
  [*F-Theory Compactifications with Multiple U(1)-Factors: Constructing Elliptic Fibrations with Rational Sections*](https://arxiv.org/abs/1303.6970),
  for a global two-point `dP2` model and its birational Tate/Weierstrass map.
- H. Pasten and C. Salgado,
  [*Non-thin rank jumps for double elliptic K3 surfaces*](https://doi.org/10.1007/s00229-024-01554-2),
  *Manuscripta Mathematica* **175** (2024), 771--781, Theorem 1.1.
- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159),
  for the geometric relation between multisections and rank jumps.

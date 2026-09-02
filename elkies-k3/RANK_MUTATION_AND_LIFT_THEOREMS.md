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

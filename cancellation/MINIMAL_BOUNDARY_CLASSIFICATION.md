# Minimal-boundary Keller classification

This note states the central classification target for the two suspension
mechanisms in this repository and proves the degree-three reduction that is
available before the final chart theorems.  It is deliberately split into:

1. an intrinsic conjecture about arbitrary nonproper Keller maps;
2. the conditional suspension theorem already supplied by the log and
   one-boundary chart results;
3. an unconditional calculation showing that the two resulting branches
   have the same unique model in geometric degree three.

Work over an algebraically closed field `k` of characteristic zero.  Two
polynomial maps `F,G:A^3 -> A^3` are **left--right equivalent** if

\[
G=L\circ F\circ R
\]

for polynomial automorphisms `L,R` of the target and source.  Write

\[
\operatorname{gdeg}(F)
=[k(x,y,z):k(F_1,F_2,F_3)].
\]

This degree is invariant under polynomial left--right equivalence.

## 1. The intrinsic minimal-boundary package

For a dominant quasi-finite map `F:U=A^3 -> Y=A^3`, let

\[
\bar X_F=\operatorname{Norm}_Y k(U),\qquad
j_F:U\hookrightarrow\bar X_F,\qquad
\partial_F=(\bar X_F\setminus U)_{\mathrm{red}}.
\tag{1.1}
\]

The finite cover, the distinguished affine open, its boundary primes, their
valuations, ramification and residue degrees, critical Fitting ideals, and
scheme intersections form the intrinsic Zariski--Main package.  This package
is functorial under polynomial left--right equivalence by `S1`; none of its
markings may be chosen from a preferred formula for `F`.

The phrase **minimal-boundary hypotheses** in the conjecture below means:

1. the package selects one geometrically integral rational critical-boundary
   normalization `E`;
2. the affine links cut out by its height-one boundary primes are saturated:
   the two coordinate rings identify after inverting exactly the declared
   prime equations, including saturation of their zero and pole supports;
3. every positive-oriented link is boundary-monotone, so every affine target
   coordinate has nonnegative valuation at the corresponding source
   boundary;
4. the smooth completion of `E` has at most two places outside `E`;
5. the divisor ledger is minimal: no unrecorded divisorial valuation occurs
   on a normal graph model, and deleting the recorded boundary destroys the
   saturated link or the ledger.

The word “intrinsic” is essential.  Conditions imposed only after choosing
coordinates on a known weighted or cancellation formula would make the
classification circular.

Because a positive-dimensional affine curve is not proper, a smooth rational
`E` in (4) has either one or two places.  Over `k` these cases are

\[
E\simeq\mathbb A^1\qquad\text{or}\qquad E\simeq\mathbb G_m.
\tag{1.2}
\]

## 2. Central conjecture

### Minimal-boundary classification conjecture

Let `F:A^3 -> A^3` be a nonproper Keller map whose intrinsic Zariski--Main
package satisfies the minimal-boundary hypotheses of Section 1.  Then `F` is
polynomially left--right equivalent to at least one of the following
suspension types:

1. a **weighted tangent suspension**, whose marked plane core is
   \[
   (W,s)\longmapsto (s,Ws-H(W));
   \tag{2.1}
   \]
2. a **cancellation suspension**, whose plane core over the parameter `P` is
   \[
   (s,Q)\longmapsto
   \left(
     Q,\,
     C\int_0^s\{1-t(Q-Pt)^m\}^r\,dt
   \right),
   \qquad m,r\ge1.
   \tag{2.2}
   \]

The two alternatives can coincide in geometric degree three.  In every
noncubic clean case their intrinsic boundary packages distinguish them.

This is stronger than the proved marked-core theorem.  In particular, the
conjecture asserts that the markings needed to recognize (2.1)--(2.2) are
extracted from (1.1), rather than supplied as external coordinates.

## 3. What the repository already proves

The existing argument can be organized as the following implication chain.

### 3.1 Normal direction

A height-one saturated link between factorial affine threefold charts has a
rank-one unit lattice.  Its boundary equations therefore satisfy

\[
D=cA\qquad\text{or}\qquad D=cA^{-1}.
\tag{3.1}
\]

The exponent is exactly one.  Boundary monotonicity turns the positive sign
into a polynomial chart.  The negative sign is the reciprocal link.  An
exact relative Picard-rank audit removes hidden graph divisors, and
effectivity removes a hidden polynomial target ledger in the reciprocal
case.

### 3.2 Critical normalization

The number of places separates the two core types.  With the primitive
residue markings supplied:

- `E=A^1` and reduced ramification force the tangent normal form (2.1);
- `E=G_m`, a primitive unit `Y=Q-sP`, and
  `s|_E=Y^{-m}` force
  \[
  D=1-s(Q-sP)^m
  \]
  and the integral core (2.2).

In geometric degree three the first bullet no longer needs a separately
supplied primitive coordinate: the
[cubic marking-extraction theorem](CUBIC_MARKING_EXTRACTION.md) uses the
Abhyankar--Moh degree theorem to straighten the embedded `A^1` by a
triangular source shear preserving the core coordinate.  The analogous
two-place assertion is false at the level of bare plane cores: the same note
constructs the infinite reduced cubic atlas

\[
 D_b=1-s^2q^b,\qquad b\ge1\text{ odd},
 \tag{3.1a}
\]

whose `G_m` valuation vectors are not individually primitive for `b>=3`.
No polynomial Keller suspension of these defects is known.

### 3.3 Threefold lift

The reciprocal branch is complete after those height-one markings are
available.  Valuation straightening produces

\[
A=1+xy^m,\qquad P=AB,\qquad Q=y+xB.
\tag{3.2}
\]

The unsliced Hensel argument then forces

\[
B=y^{m+1}h_q(A)+A^{r+1}z
\tag{3.3}
\]

with the unique cancellation jet.  It also forces the full Stein degree and
the remaining rank-two flags; they are not additional assumptions.

In the positive branch, the plane core is known but an arbitrary positive
polynomial vertical chart has not been straightened to the explicit weighted
formulas.  The genuinely open implications are therefore:

1. **suspension and two-place marking extraction:** recover the
   coordinate-preserving core from the Zariski--Main package and, in the
   `G_m` branch, identify its quadratic coefficient as the primitive
   transverse conormal class and prove noncontraction; Theorem 3.3 of the
   cubic marking note then produces the primitive residue coordinate and its
   affine-linear lift automatically;
2. **positive-chart straightening:** identify every resulting positive
   polynomial Keller square with the weighted vertical chart.

In degree three the second item now has a positive straightening theorem
under explicit intrinsic saturation labels.  Once the ledger gives
`x=C/gamma`, it remains to extract

\[
 y=(W-\gamma)/C,\qquad
 z=(\gamma-1+\tfrac32xy)/(b x^2)
\tag{3.3a}
\]

as primitive intrinsic conormal classes, together with noncontraction,
primitive Jacobian content, and boundary Stein degree one.  Theorem 4.2 of
the [cubic marking note](CUBIC_MARKING_EXTRACTION.md) then produces the slice
and proves the weighted chart; cubic target polynomiality forces the
coefficient `-3/2` rather than assuming it.

Neither follows merely from rationality of `E`, its unit rank, or the
determinant ledger.  The countermodels in
[`ONE_BOUNDARY_CHART_CLASSIFICATION.md`](ONE_BOUNDARY_CHART_CLASSIFICATION.md)
show why those weaker statements are false.

Consequently the bare conjecture in Section 2 remains open.  The proved
one-boundary theorem gives the conjectured dichotomy after adding the
primitive boundary markings and, in the positive branch, weighted-chart
straightenability.

### 3.4 A finite-normalization frontend in degree three

There is now a second cubic route which does not assume a suspension square.
Let `B` be the rank-three normalization algebra of the canonical finite
cover.  The
[cubic normalization frontend](CUBIC_NORMALIZATION_FRONTEND.md) proves:

1. normality makes `B` flat outside a zero-dimensional intrinsic defect
   \[
   Z_{\rm flat}=V(\operatorname{Fitt}_3(B));
   \]
2. if this defect is empty, Deligne--Faddeev and Quillen--Suslin extract a
   global binary cubic, unique up to its Tschirnhausen basis;
3. if the resulting coefficient morphism is affine-linear of full rank and
   the distinguished affine source is the full simple-root locus, primitive
   resultant normalization is automatic;
4. the three cubic hyperplane orbits then force the tangent-nonosculating
   orbit, whose normalized source and map are the foundational model.

Thus, in this finite-flat affine-linear frontend, neither a
coordinate-preserving plane core nor the positive conormal quotients need be
supplied.  The remaining bare-package implication is sharpened to three
intrinsic assertions: point-flatness, affine-linearity modulo polynomial
Tschirnhausen gauge of the binary-cubic coefficient morphism, and absence of
an extra unramified boundary prime.  Affine Hartogs maximality rules out any
extra simple boundary supported only in codimension at least two.  The cubic
DVR degree budget also proves that the critical target divisor already has
exact decomposition `(e,f)=(2,1)+(1,1)`, so any unramified boundary prime
would have to map to a distinct second nonproperness divisor.
The first is a codimension-three condition and therefore cannot be read from
the height-one ledger alone.

### 3.5 The two cubic gateways

The branchwise and finite-normalization programs should not be run as a
single cumulative checklist.  They are alternative gateways.

#### Theorem 3.5 -- cubic gateway reduction

Let `F:A^3 -> A^3` have geometric degree three.  Either of the following
packages implies that `F` is polynomially left--right equivalent to the
foundational map:

1. **suspension gateway:** the intrinsic package extracts a
   coordinate-preserving one-boundary suspension, the two-place branch has
   the primitive quadratic conormal marking, and the positive branch has
   the quotient/Stein labels of the cubic positive-chart theorem;
2. **normalization gateway:** the canonical rank-three normalization has
   curvilinear closed collision fibers (or, more generally, empty
   point-flatness defect); curvilinearity follows, in particular, from the
   ramification/cotangent extension package of Proposition 1.13 in the
   cubic normalization frontend: `S_2` for the ramification support, `S_1`
   for its cotangent module, and codimension-one primitive generation.
   Equivalently, the two `Ext` obstructions of Proposition 1.14 vanish, or
   the two finite double-saturation quotients of Proposition 1.15 are zero.
   It also has no unramified boundary prime
   over a second target divisor, and its binary-cubic orbit has a full-rank
   tangent-hyperplane representative.  It is enough for that representative
   to lie in an invariant unipotent gauge orbit classified by the cubic
   gauge-straightening theorem.

#### Proof

For the first gateway, the cubic marking theorem supplies the weighted or
cancellation normal form, and Theorem 4.1 below identifies both with the
foundational map.

For the second, intrinsic curvilinearity gives flatness by the local
monogenicity theorem (or flatness is assumed directly).  Affine Hartogs
maximality and the absence of an unramified boundary prime identify the
Keller source with the full simple-root locus.  Deligne--Faddeev gives the
binary cubic, resultant normalization is automatic, and the tangent
hyperplane theorem gives the foundational map.  Invariant unipotent gauges
are removed by the gauge-straightening theorem.
QED

Once either gateway is crossed, the other structure exists automatically:
the foundational map has both its weighted/cancellation presentations and
its flat tangent-hyperplane cubic normalization, and these properties
transport under polynomial left--right equivalence.  Thus the open program
is to prove **one** gateway from the bare package.  Proving suspension labels
and coefficient gauge-linearity independently would duplicate work.

There is one important non-implication before a gateway is complete.
Height-one saturation in the current suspension theorem compares affine
models only after inverting their boundary equations.  It does not control
the finitely many closed points of the canonical normalization and therefore
does not eliminate its Fitting flatness defect.  A closed-point boundary
atlas with Cartier Cohen--Macaulay boundary would eliminate that defect by
Proposition 1.2 of the cubic normalization frontend, but it is a genuine
strengthening of the current height-one hypothesis.

There is a minimal intrinsic replacement for that stronger atlas.
Proposition 1.5 of the same frontend proves that point-flatness is equivalent
to **cubic fiber-minimality**: every scheme-theoretic canonical
normalization fiber has length three.  Any defect has length at least four.
This is a finite closed-point calculation in the intrinsic package.  It
cannot be inferred merely from “at most two places at infinity,” because
that reduced-place count does not see fiber multiplicity.

Nor can cubic algebra structure alone supply the missing implication.
Proposition 1.6 constructs the exact reflexive trace-module defect
`Fitt_3=(x,y,z)` and invokes the nonflat triple-cover correspondence to show
that such modules do occur in integral normal cubic covers.  What remains
to be proved is therefore genuinely Keller-boundary-specific: the
distinguished affine open and its scheme intersections must exclude the
length-four fiber.

The local problem is now finite and determinantal.  Proposition 1.7 writes
every defect trace module as the cokernel of an `(s+2)-by-s` matrix, with
fiber length `s+3` and defect ideal its maximal minors.  The first rung is a
single parameter column `(a,b,c)^T`; if the defect scheme is reduced, it is
formally the unique Koszul column `(x,y,z)^T`.  Higher rungs have maximal
minors in the square of the maximal ideal, so every reduced defect is
automatically this first rung.  Thus a reduced-defect argument need only
exclude this one local model, rather than arbitrary normal
non-Cohen--Macaulay cubic algebras.

That local model is now rigid as an algebra as well as a module.
Proposition 1.8 proves that its fiber is the square-zero algebra
`k plus k^3`, so it consists of one boundary point and the Keller map omits
the target point.  Over the critical target divisor, the generic ramified
boundary sheet and the generic affine sheet must therefore collide at that
point.  Corollary 1.9 turns reduced point-flatness into an intrinsic
collision-type theorem for precisely the scheme intersection already
recorded by the Zariski--Main package.  The flat foundational model
calibrates the distinction: its sheets meet in a curvilinear length-three
triple-root fiber, whereas the reduced defect is the noncurvilinear
square-zero length-four fiber.

Proposition 1.10 closes the entire point-flatness gap under one compact
intrinsic marking: every canonical collision fiber is curvilinear.  A
curvilinear fiber is monogenic; Nakayama lifts its generator to the local
rank-three algebra, which is then a free monic-cubic algebra.  This one
condition eliminates reduced, nonreduced, and higher determinantal defects
while retaining the foundational curvilinear triple-root collision.  Thus
the missing closed-point theorem should be stated as **intrinsic
curvilinearity**, not as global sheet disjointness or a full
Cartier-boundary atlas.

Proposition 1.11 makes this marking coordinate-free:
curvilinearity is equivalent to cyclicity of the relative cotangent module
on every collision fiber, or to the unit first Fitting ideal
`Fitt_1(Omega)`.  Relative differentials, their Fitting ideals, and completed
stratum maps already belong to the intrinsic package.  The issue is now to
derive cotangent cyclicity from saturation and monotonicity, not to add an
external root coordinate.

Equivalently, Proposition 1.12 asks for **closed-point conormal
saturation**: the primitive conormal class must generate the entire
nilradical of every collision fiber.  This is the exact bridge to the
existing branchwise markings.  The two-place theorem already extracts a
primitive quadratic conormal class generically, and the positive theorem
uses primitive quotient/Stein data; what neither currently proves is that
the chosen class stays a generator through every closed collision.  Proving
that extension simultaneously removes the normalization flatness defect
and supplies the local root marking.

Proposition 1.13 supplies the exact Hartogs theorem for this extension.
If the pure two-dimensional scheme-theoretic ramification support is `S_2`,
its rank-one full-support cotangent module is `S_1`, and the primitive class
generates in codimension one, then local cohomology forbids any closed-point
cokernel.  The class generates globally, so cotangent cyclicity and
point-flatness follow.  Proposition 1.14 sharpens this to a two-module
ledger: `Ext_A^2(T,A)` is the support-depth defect, and once it vanishes
`Ext_A^3(Omega_{B/A},A)` is the canonical dual of the primitive-generation
cokernel.  Proposition 1.15 makes both defects literal finite saturation
quotients.  If `C=T^[2]=Ext_A^1(Ext_A^1(T,A),A)` is the canonical `S_2`
hull, then `Ext_A^2(T,A)` is dual to `C/T`; after `C=T`,
`Ext_A^3(Omega_{B/A},A)` is dual to `Omega_{B/A}/T tau`.  The missing
closed-point theorem is exactly `C=T` and `Omega_{B/A}=T tau`.  This is
substantially weaker than assuming the entire canonical boundary is
Cartier and normal.

Propositions 1.16--1.17 compress this once more.  With `Z` the closed
collision locus and `P=H_Z^0(Omega_{B/A})`, there is an exact sequence

\[
 0\to P\to \Omega_{B/A}/T\tau\to C/T
 \to H_Z^1(\Omega_{B/A})\to0.
\]

Hence after `C=T`, the conormal defect is exactly `P`.  If
`F_1 -> F_0 -> Omega_{B/A}` has image `N` and
`I=Fitt_3(B)`, then

\[
 P=(N:I^\infty)/N.
\]

The entire normalization gap is therefore one canonical-bidual
surjectivity test and one module-saturation equality.

The extra-boundary condition likewise has a single intrinsic detector.
Proposition 2.1a writes the reduced nonproperness equation as

\[
 j_F=\delta_Fu_F,
\]

where `delta_F` is the unique cubic branch equation.  The factor `u_F` is a
unit exactly when no unramified boundary divisor remains.  Thus target
irreducibility of the nonproperness hypersurface closes this gap without
classifying compactifications.

These local generators cannot be globalized away.  Proposition 2.2 of the
frontend shows that a global monogenic cubic presentation would make its
derivative a constant unit on the full Keller étale locus, contradicting
geometric degree three.  The binary-root transition, which distinguishes
the one-place and two-place markings, is therefore essential rather than a
coordinate artifact.

The remaining coefficient-gauge problem is now graded.  Propositions
2.6--2.7 of the gauge note exclude every pair of nonzero monomial times and
show that a general multi-monomial two-shear recursion has only one scalar
obstruction, in the discriminant line, every fourth degree.  Proposition
2.8 proves that the first bilinear candidate always projects to zero on
that line; the exact quadratic recursion remains soluble through degree
eight on every coupled kernel-basis direction.  This closes the naive
first-obstruction attack negatively and makes extraction of the intrinsic
root line followed by Borel straightening the preferred route.
Propositions 2.9--2.10 make that route exact: a saturated projective flag
with scalar mixed coefficient gives `C_1=q`, after which the Keller map is
the pullback of the foundational tangent-hyperplane map along
`G=(C_0,C_2,C_3)`.  It remains to prove foundational base-change rigidity:
the minimal-boundary pullback with total simple-root space `A^3` forces
`G` to be a polynomial automorphism.  The ranked closure route through
these three certificates is recorded in
[`CUBIC_CLOSURE_ATTACKS.md`](CUBIC_CLOSURE_ATTACKS.md).

## 4. Cubic collapse inside the two branches

The following result is unconditional; it is a calculation inside the union
of the two classified suspension families.

### Theorem 4.1 -- unique cubic suspension model

Every weighted tangent suspension of geometric degree three and every
cancellation suspension of geometric degree three is polynomially
left--right equivalent to the foundational map

\[
F_0(x,y,z)=
\left(
\begin{aligned}
&u^3z+y^2u(4+3xy),\\
&y+3xu^2z+3xy^2(4+3xy),\\
&2x-3x^2y-x^3z
\end{aligned}
\right),
\qquad u=1+xy.
\tag{4.1}
\]

#### Weighted branch

After the intrinsic affine marking puts the double root at `0` and the other
distinguished root at `1`, an admissible cubic primitive has

\[
H(0)=H'(0)=H(1)=0.
\tag{4.2}
\]

Writing `H=a_3W^3+a_2W^2+a_1W+a_0`, equations (4.2) give

\[
a_0=a_1=0,\qquad a_2=-a_3.
\]

Hence

\[
H(W)=cW^2(1-W),\qquad c\in k^*.
\tag{4.3}
\]

There is no cubic weighted modulus.  This is the direct degree-three form of
the general dimension count `N-3`.

Let `G_{c,b}` be the weighted lift of (4.3), with nonzero vertical parameter
`b`, and let `G_{1,1}` be the normalized lift.  Directly from the weighted
formulas,

\[
G_{c,b}(x,y,z/b)
=\bigl(G_{1,1;1},\,cG_{1,1;2},\,G_{1,1;3}\bigr).
\tag{4.4}
\]

Moreover,

\[
F_0(x,y,z)
=\operatorname{diag}(1,2,2)\,
 G_{1,1}(x,y,-z/2).
\tag{4.5}
\]

All changes in (4.4)--(4.5) are polynomial automorphisms.

#### Cancellation branch

The generic degree of cancellation type `(m,r)` is

\[
N=r(m+1)+1.
\tag{4.6}
\]

For `N=3` and positive integers `m,r`, equation (4.6) has the unique solution

\[
(m,r)=(1,1).
\tag{4.7}
\]

Its parameter polynomial and unique normalized cancellation jet are

\[
M_{1,1}(q)=q-3,\qquad h(A)=3+9A.
\tag{4.8}
\]

Every additional term in `A^2k[A]` is removed by the polynomial source
shift in `z` from the cancellation construction, so (4.8) represents the
entire cubic cancellation branch up to source equivalence.

Put

\[
A=1+xy,\qquad
B=A^2z+y^2(3+9A),\qquad
P=AB,\qquad Q=y+xB,
\]

\[
R_C=C\int_0^{x/A}\{1-t(Q-Pt)\}\,dt.
\tag{4.9}
\]

Then exact polynomial identities give

\[
\boxed{
F_0
=\operatorname{diag}(1/3,1,2/C)
\circ(P,Q,R_C)\circ(x,y,3z).
}
\tag{4.10}
\]

Thus the first cancellation rung is not a second cubic class: it is the same
foundational map.

### Corollary 4.2 -- conditional minimal-boundary cubic classification

Assume the minimal-boundary classification conjecture.  If
`F:A^3 -> A^3` is a nonproper Keller map satisfying its intrinsic hypotheses
and `gdeg(F)=3`, then

\[
\boxed{F\ \text{is polynomially left--right equivalent to}\ F_0.}
\tag{4.11}
\]

The same conclusion is already a theorem without assuming the full
conjecture for any `F` to which the existing marked one-boundary dichotomy
applies and whose positive chart, if present, is known to be weighted
straightenable.

The distinction matters: Theorem 4.1 proves that the two candidate branches
collapse; it does not by itself extract those branches from the stated
intrinsic hypotheses.

## 5. The degree-three proof target

The first classification milestone can now be stated without quantifying over
all degrees:

> **Cubic extraction theorem (open).** Every geometric-degree-three
> nonproper Keller map satisfying the intrinsic minimal-boundary hypotheses
> supplies a coordinate-preserving suspension core; in the two-place case
> its ambient boundary algebra selects the primitive quadratic conormal
> coefficient (hence the cancellation marking),
> and in the positive case its polynomial chart straightens to a weighted
> chart.

Theorem 4.1 would then finish the classification immediately.  A proof should
therefore target the two extraction assertions, not repeat the cubic normal
form calculation.

Useful degree-three specializations are:

1. the critical incidence cover has `S_3` monodromy and no target-fixed deck
   transformation;
2. the one-place primitive plane-core marking is automatic by the cubic
   marking-extraction theorem;
3. in the two-place branch, (4.6) forces both the ramification exponent and
   tropical slope to be one;
4. the cancellation residue polynomial is linear, eliminating parameter and
   hidden-residue ambiguity;
5. both orientations land on the explicitly verified identities
   (4.5) and (4.10).
6. the canonical normalization is flat away from finitely many target
   points, and an empty flatness defect produces a global binary-cubic
   coefficient morphism;
7. affine-linearity of that morphism bypasses the branchwise extraction
   problem and forces the foundational map directly.
8. every nonlinear coefficient slice
   `C_1=q-3C_0h` with `q!=0` and `h` invariant under the translation
   locally nilpotent derivation is polynomially Tschirnhausen-gauge-equivalent
   to the tangent hyperplane; the opposite shear is symmetric.  These are
   exactly the polynomial automorphisms obtained from a single variable-time
   upper or lower unipotent shear.
9. the critical target divisor exhausts the cubic DVR degree as one
   ramified boundary sheet of index two plus one affine simple sheet; any
   extra simple boundary must lie over a separate unramified target divisor.

Items 1--4 compress the marking problem substantially, but none alone proves
that the `G_m` ambient function is intrinsically selected, that the
coordinate-preserving suspension exists, or that a positive vertical chart
is polynomially conjugate to the weighted one.  Items 6--9 give an
alternative frontend, but the minimal-boundary hypotheses have not yet been
shown to eliminate the point defect or force coefficient linearity.

## 6. Exact verification

Run

```bash
.venv/bin/python scripts/verify_minimal_boundary_cubic.py
```

The checker verifies (4.3)--(4.10), including the arbitrary nonzero weighted
parameters and cancellation scale.  It is a certificate for the branch
collapse, not a certificate for the open extraction theorem.

The structural inputs are:

- [`../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md`](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md);
- [`../verified/TANGENT_MAP_CORE.md`](../verified/TANGENT_MAP_CORE.md);
- [`LOG_GEOMETRY_OF_SUSPENSIONS.md`](LOG_GEOMETRY_OF_SUSPENSIONS.md);
- [`ONE_BOUNDARY_CHART_CLASSIFICATION.md`](ONE_BOUNDARY_CHART_CLASSIFICATION.md);
- [`CUBIC_MARKING_EXTRACTION.md`](CUBIC_MARKING_EXTRACTION.md);
- [`CONSTRUCTION.md`](CONSTRUCTION.md).

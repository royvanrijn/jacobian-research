# Ind-group constructibility and bounded orbit families

The complexity-filtered contact program has two logically separate parts.
The first is finite-type algebraic geometry: after fixing source degree,
target degree, parameter height, pole order, and stabilization dimension, the
space of left--right equivalences is a finite-type incidence scheme.  This
gives constructibility and a useful bounded-box algebraization theorem.  The
second is a coercivity problem: prove that the degree-five arc leaves every
fixed resource box after minimizing over target automorphisms.

The first part is established below.  It sharpens the exact list of remaining
steps in [complexity-filtered contact](COMPLEXITY_FILTERED_CONTACT.md).  The
current degree-five calculation already supplies the source-only coercivity
in every stabilization dimension, but target coercivity and a no-escape
statement at the base point remain open.

Work over a characteristic-zero field `k`.  The generic-uniformity theorem in
Section 3 additionally assumes that `k` is algebraically closed and
uncountable.

## 1. Paired finite-degree automorphism strata

For integers `N,D>=1`, let

\[
 \mathfrak A_{N,D}^{\pm}
 \subset \operatorname{End}_{\le D}(\mathbb A^N)^2
\]

be the coefficient scheme of pairs `(A,A^-)` satisfying

\[
 A\circ A^-=A^-\circ A=\operatorname{id}.          \tag{1.1}
\]

It is an affine finite-type scheme: both maps have finitely many coefficients
and (1.1) is a finite list of polynomial coefficient equations.  A point is a
polynomial automorphism together with an inverse, both of degree at most `D`.
The determinant-one locus is cut out by the additional equation
`det DA=1`.

This paired model is better suited to the present problem than a one-sided
degree filtration.  It is symmetric by construction, works over arbitrary
coefficient rings, and makes inversion a coordinate swap.  Over a field the
classical inverse-degree estimate

\[
 \deg A^{-1}\le (\deg A)^{N-1}                     \tag{1.2}
\]

shows that it defines the same ind-topology as the usual one-sided degree
filtration, up to reindexing.

The standard ind-variety results say more generally that `Aut(X)` is locally
closed in the ind-semigroup `End(X)` for an affine variety `X`.  For affine
space, however, the paired scheme above is enough and avoids importing that
machinery into any proof.

There is one important group-theoretic warning.  For constructibility the
ambient group should be `Aut`, or its closed determinant-one locus `Aut_1`,
not an unqualified bounded-degree piece of

\[
 \operatorname{SAut}(\mathbb A^N)
 =\langle\mathbb G_a\text{-subgroups}\rangle .
\]

No equality `SAut=Aut_1` is being used in dimension at least three, and a
degree bound on the resulting automorphism does not bound the number or
degrees of its additive-group factors.  If factorization complexity is
desired, the number and degree of factors must be added as separate resources.
The finite-jet interpolation theorem still gives `SAut` witnesses; enlarging
the minimization to `Aut_1` only makes a lower-bound problem stronger.

## 2. Bounded left--right incidence is constructible

Let `S` be a finite-type `k`-scheme and

\[
 \mathcal F:S\times\mathbb A^n\longrightarrow\mathbb A^n
\]

a polynomial family of uniformly bounded coordinate degree.  Fix a reference
point `s_0`, a stabilization dimension `r`, and put `N=n+r`.  For a degree
budget `D`, define the incidence scheme

\[
 \mathcal I_{r,D}\subset
 S\times(\mathfrak A_{N,D}^{\pm})^2               \tag{2.1}
\]

by

\[
 B\circ(\mathcal F_s\times\operatorname{id}_r)\circ A
 =\mathcal F_{s_0}\times\operatorname{id}_r.       \tag{2.2}
\]

Equation (2.2) is again a finite list of polynomial coefficient equations.
Thus `I_(r,D)` is finite type.  Chevalley's theorem gives:

### Proposition 2.1

The set

\[
 E_{r,D}=\{s\in S:\mathcal F_s
 \text{ is LR-equivalent to }\mathcal F_{s_0}
 \text{ in stabilization }r
 \text{ with two-sided degree }\le D\}             \tag{2.3}
\]

is constructible.

The same construction over `S x S` makes the bounded-degree equivalence
relation constructible.  This is the precise finite-type content of filtering
the automorphism ind-group by degree.

## 3. What pointwise equivalence gives generically

The countability of the resource indices now has a useful consequence.

### Theorem 3.1 -- generic boundedness after finite base change

Assume that `k` is algebraically closed and uncountable and that `S` is
irreducible.  Suppose that every closed fiber `F_s` is stably polynomially
left--right equivalent to `F_(s_0)`, with no a priori bound on automorphism
degree or stabilization dimension.  Then there are fixed integers `r,D` and
a dense open `U subset S` such that:

1. `E_(r,D)` contains `U`;
2. some irreducible component of `I_(r,D)` dominates `U`;
3. after a generically finite dominant base change `T -> U` and shrinking
   `T`, equation (2.2) is realized by one regular family of automorphism pairs
   of two-sided degree at most `D`.

### Proof

The countable union of the constructible sets `E_(r,D)` contains `S(k)`.
An irreducible finite-type variety over an uncountable algebraically closed
field is not a countable union of proper closed subvarieties.  Hence some
`E_(r,D)` is dense, and a dense constructible set contains a dense open `U`.

An irreducible component of its incidence scheme then dominates `U`.  Choose
a closed point of the generic fiber.  Its residue field is a finite extension
of `k(U)`.  Normalizing `U` in that extension and shrinking spreads the point
to the asserted regular family.  QED

This proves a substantial part of the desired pointwise-to-family bridge, but
not all of it.  It has three unavoidable qualifications:

- the equivalence may require a finite base change;
- the family need only exist generically;
- it can acquire poles when continued to the chosen base point `s_0`.

The last issue is not a technicality.  The closed incidence

\[
 \{(t,a):ta=1\}\ \cup\ \{(0,0)\}\subset\mathbb A^2 \tag{3.1}
\]

has a point over every `t`, but its only dominant branch is `a=1/t`.  A finite
base change does not make that branch regular at zero.  Consequently
constructibility alone does **not** turn pointwise equivalence into a based
regular automorphism arc.  One needs a no-escape lemma for the particular
orbit incidence, or one must record pole order as part of the complexity.

## 4. The full resource spectrum

For a rational function `f(t)=p(t)/q(t)` in reduced form set

\[
 h_t(f)=\max\{\deg p,\deg q\},
 \qquad
 \nu_0(f)=\max\{0,-\operatorname{ord}_0(f)\}.       \tag{4.1}
\]

For an automorphism pair, take the maximum over every coefficient of the map
and its displayed inverse.  A resource box is a tuple

\[
 Q=(r,D_x,D_t,\nu_0),                               \tag{4.2}
\]

recording stabilization dimension, the maximum two-sided coordinate degree
of the source and target pairs, their maximum parameter numerator/denominator
degree, and their pole order at the base point.  Thus `D_t` is the height
`h_t`, not merely the numerator degree.  On a general parameter curve, `D_t`
should be replaced by the degree of the pole divisor relative to a fixed
compactification.

For an arc `F_t` through `F_0`, let `K_m(F_t/F_0)` be the upward-closed set of
boxes `Q` for which there are rational source and target automorphism pairs in
that box satisfying

\[
 B(t)\circ(F_t\times\operatorname{id}_r)\circ A(t)
 \equiv F_0\times\operatorname{id}_r
 \pmod {t^{m+1}}.                                  \tag{4.3}
\]

For based regular contact impose `nu_0=0` and `A(0)=B(0)=id`.  Meromorphic
contact uses the `t`-adic valuation in (4.3), allowing cancellation of bounded
principal parts.

The set-valued spectrum

\[
 m\longmapsto K_m(F_t/F_0)                          \tag{4.4}
\]

is preferable to a single minimum.  A scalar minimum hides tradeoffs between
source degree and parameter height.  More seriously, minimizing over `r` at
each `m` allows the stabilization dimension to grow with the contact order,
whereas stable equivalence permits one **fixed** finite stabilization.  The
correct stable obstruction is therefore quantified as

\[
 \text{for every fixed }r,D_x,D_t,\nu_0,
 \text{ some finite contact order leaves that box}. \tag{4.5}
\]

### Theorem 4.1 -- bounded-box algebraization

Fix a resource box `Q`.  If `Q` belongs to `K_m(F_t/F_0)` for every `m`, then
there is one exact rational left--right equivalence satisfying (4.3) with
`equiv` replaced by equality and lying in the same box.  If `nu_0=0` and the
witnesses are based, the exact equivalence is regular and based at zero.

### Proof

Numerators and denominators of bounded degree form a finite-type coefficient
space after a normalization of scalar representatives.  Take its product
with the paired automorphism strata from Section 1.  There are only finitely
many denominator-valuation charts compatible with the fixed pole bound.  On
each chart the inverse equations and contact through order `m` cut out a
closed subset.  If witnesses exist at every order, one chart occurs at
arbitrarily high orders and hence at every lower order.  In that chart the
contact equations make

\[
 Z_0\supset Z_1\supset Z_2\supset\cdots             \tag{4.6}
\]

a descending chain of closed subsets in one Noetherian space.  If every
`Z_m` is nonempty, the chain stabilizes.  A point in the stable set makes all
coefficients of the error in (4.3) vanish, hence gives an exact equivalence.
QED

Thus escape from every fixed full resource box is not merely a necessary
condition for nonalgebraization; it is equivalent to nonexistence of an exact
bounded rational family.  This is the clean algebraization bridge that the
source-degree profile alone cannot provide.

Composing the arc with one fixed regular left--right gauge changes the boxes
by bounded degree multiplication and bounded height addition.  Therefore the
properties "some fixed box works for all orders" and "every fixed box is
eventually left" are intrinsic to the LR orbit germ.  Finer numerical growth
should be recorded only up to this bounded distortion; the exact coefficient
`34` is a gauge computation, while boundedness, unboundedness, and polynomial
growth exponent are plausible intrinsic outputs.

## 5. What is already stable for the degree-five arc

Let `G_t` be the determinant-normalized degree-five arc through `F_2`.  Its
canonical source trivializer satisfies

\[
 \deg_x V_m=\deg_x W_m=34m+1.                      \tag{5.1}
\]

This source-only lower bound already survives every identity stabilization.

### Proposition 5.1 -- stabilization does not change the canonical bound

For every fixed `r>=0`, the unique canonical source trivializer of

\[
 G_t\times\operatorname{id}_r
 \quad\text{relative to}\quad
 F_2\times\operatorname{id}_r
\]

is

\[
 \widehat\alpha_t\times\operatorname{id}_r.        \tag{5.2}
\]

Consequently the stabilized source-only canonical lower bound is still

\[
 b_m^{\mathrm{src},(r)}=34m+1.                     \tag{5.3}
\]

### Proof

The displayed product is a formal source trivializer.  Uniqueness in the
formal source-triviality theorem, applied in dimension `3+r`, says it is the
only one.  Its forward and inverse coefficient degrees are unchanged.  QED

Hence the existing computation proves, without any new hypothesis:

\[
 \boxed{
 \text{For every fixed stabilization dimension, there is no bounded-degree
 source algebraization in the chosen target gauge.}}                 \tag{5.4}
\]

This is stronger than the unstabilized statement presently highlighted in
the contact note.

It is not yet a left--right theorem.  A bounded-coordinate-degree formal
target gauge can have unbounded parameter height, and lifting that target
gauge through a noninvertible Keller map can create source coefficients of
unbounded degree.  Thus (5.1) by itself cannot distinguish high intrinsic LR
complexity from complexity transferred into the target or parameter.

## 6. The target-minimality problem

There are two viable routes.

### Route A: an intrinsic LR minimization

Work directly with the spectrum (4.4).  At order `m`, after fixing lower
coefficients, the new source and target coefficients enter through the linear
operator

\[
 (U,V)\longmapsto U\circ F_2+DF_2\,V.              \tag{6.1}
\]

The Catalan leading term producing (5.1) defines a class in the associated
graded cokernel of (6.1) with both the source and target degree filtrations.
The desired coercivity statement is:

\[
 \max\{\deg_x U_m,\deg_x V_m\}\ge c_r m-O_r(1)     \tag{6.2}
\]

for every solution in stabilization dimension `r`, or at least that the left
side is unbounded.  Proving (6.2) requires showing that the monomial

\[
 M^m=x^{20m}y^{8m}z^{6m}                            \tag{6.3}
\]

cannot be cancelled by the pullback of a lower-degree target vector field.
This is a concrete filtered-module or Gröbner problem.  It should be attacked
first for `r=0`, then with the extra identity variables retained symbolically.

If coordinate degree can be kept bounded, Theorem 4.1 says that parameter
height or pole order must escape instead.  This is why all three complexity
coordinates belong in the final invariant.

### Route B: an intrinsic target slice

Alternatively, impose canonical target conditions and prove that every LR
equivalence has a bounded-complexity representative in that slice.  The two
intrinsic boundary vertices and the marked-cover reconstruction suggest such
a slice, but using them would make this route depend on the normalization
boundary invariant.  It would prove target-gauge minimality efficiently, but
would not yet produce a conceptually independent invariant.

Route A is therefore the route to a genuinely different stable invariant.
Route B remains a useful intermediate theorem and a check on the computation.

## 7. The remaining no-escape statement

Theorem 3.1 turns pointwise stable equivalence into a bounded-degree regular
family only on a dense open after finite base change.  For an arc centered at
`F_2`, a true pointwise-to-based-family theorem still needs one of:

1. a valuative no-escape theorem for the bounded incidence `I_(r,D)`;
2. a proof that a dominating incidence component passes through a chosen
   equivalence of the special fiber;
3. a canonical slice with a regular section;
4. direct control showing that the minimal pole order is zero.

Without one of these inputs, "every nearby fiber lies in the orbit" and "the
arc is regularly trivial in the orbit" remain different statements.
Parameter pole order measures exactly that difference.

## 8. Resulting theorem ladder

The proposed stable invariant can now be developed in the following order.

1. **Finite-type strata -- proved here.**  Fixed two-sided degree and fixed
   stabilization give a finite-type LR incidence and a constructible orbit
   relation.
2. **Generic uniformity -- proved here.**  Pointwise stable equivalence over
   an uncountable base gives fixed `r,D` generically and a regular equivalence
   after generically finite base change, away from the base point.
3. **Full bounded-box algebraization -- proved here.**  Uniform bounds on
   coordinate degree, parameter height, pole order, and stabilization turn
   arbitrary finite contacts into one exact rational family.
4. **Stable source coercivity -- proved by the degree-five calculation.**
   The law `34m+1` survives every fixed identity stabilization in the chosen
   target gauge.
5. **Target coercivity -- open.**  Prove an associated-graded lower bound for
   (6.1), or construct an intrinsic target slice.
6. **No escape at the base -- open if a pointwise theorem is desired.**  Rule
   out poles in a bounded-degree dominating incidence branch.

Items 1--4 already define a rigorous asymptotic complexity framework.  Items
5--6 are exactly what separates that framework from a new proof of stable
degree-five moduli.  Once item 5 holds for every fixed stabilization
dimension, the escape spectrum (4.4) is a stable LR invariant genuinely
different in construction from boundary normalization.

## 9. References for the ind-group input

- J.-P. Furter and H. Kraft,
  [*On the geometry of the automorphism groups of affine varieties*](https://arxiv.org/abs/1809.04175),
  arXiv:1809.04175.  This supplies the general locally closed ind-group
  framework for `Aut(X)` inside `End(X)`.
- H. Bass, E. Connell, and D. Wright,
  [*The Jacobian conjecture: reduction of degree and formal expansion of the inverse*](https://doi.org/10.1090/S0273-0979-1982-15032-7),
  Bull. Amer. Math. Soc. 7 (1982), including the inverse-degree estimate
  attributed to Gabber.

Neither reference supplies the bounded-box contact theorem or the
pointwise-to-family formulation above; those are the elementary finite-type
consequences recorded in this note.

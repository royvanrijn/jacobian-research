# Intrinsic marked-line extraction

## 0. Purpose

The minimal-boundary program currently uses chart straightenability (`CS`)
as the last predicate in its suspension gateway.  That predicate already
asks for an isomorphism with one of the known weighted or cancellation
charts.  It is therefore too strong to serve as the missing intrinsic
classification theorem, and it does not leave room for the
root-engineered quadratic-gauge algebraization.

This note separates the open problem into two statements:

\[
\boxed{
 \text{canonical boundary data}
 \Longrightarrow
 \operatorname{MLE}
 \Longrightarrow
 \text{classification of polynomial algebraizations}.
}
\tag{0.1}
\]

The first arrow is the **intrinsic marked-line extraction problem**.  The
second is a controlled algebraization problem.  Neither arrow is included
in the definition of `MLE`.

Throughout, work over an algebraically closed field \(k\) of
characteristic zero.

## 1. Canonical input

Let

\[
 F:\mathbb A^3\longrightarrow\mathbb A^3
\]

be a nonproper Keller map, and let

\[
 \pi_F:\overline X_F\longrightarrow\mathbb A^3
\]

be its canonical finite normalization.  The input to extraction is:

1. the distinguished open
   \(\mathbb A^3\subset\overline X_F\);
2. the complete list of boundary primes, with their valuations,
   ramification indices, residue degrees, differents, and completed
   incidence maps;
3. the ordered list of divisorial target images of those primes;
4. a uniquely selected ramified boundary prime \(D_\Delta\);
5. its saturated primitive conormal class \(\tau\), including its
   specialization at every recorded collision;
6. the scheme structures on all intersections of target boundary images.

The bound relevant to the three known constructions is **at most two target
boundary images**.  It is not a bound of two on the number of source
boundary primes.  Several source primes may lie over the second target
image, and their full \((e,f)\)-multiset must be retained.

## 2. Marked controlled-incidence presentations

### Definition 2.1 -- marked controlled-incidence presentation

A **marked controlled-incidence presentation** of \(F\) consists of the
following data, functorially attached to the canonical input of Section 1.

1. **Intrinsic quotient flag.**  Normal affine quotient rings and rational
   quotient arrows give a commuting square
   \[
   \begin{array}{ccc}
   \mathbb A^3&\dashrightarrow&\mathbb A^3\\
   \downarrow&&\downarrow\\
   Z&\dashrightarrow&T,
   \end{array}
   \tag{2.1}
   \]
   where \(Z,T\) are normal affine surfaces.  The flag and its marked
   divisors are fixed by every automorphism of the canonical
   finite-normalization object.
2. **Primitive root and coefficient flag.**  On a declared birational
   chart there are rational functions \(p,s,q\) with
   \[
   k(\mathbb A^3)=k(p,s,q),
   \tag{2.2}
   \]
   where \(s\) is the primitive marked-root parameter induced by \(\tau\),
   and \(p\) is a preserved quotient parameter.  Changing \(s\) by an
   arbitrary rational coordinate is not allowed: its divisor and residue
   class are part of the marking.
3. **Controlled divisor.**  There is an irreducible divisor
   \(\Delta\) and an inverse equation
   \[
   E(s;y)=0
   \tag{2.3}
   \]
   for the finite cover such that
   \[
   \partial_sE=u\Delta^e,\qquad u\in k^*,\quad e\ge1.
   \tag{2.4}
   \]
   The normalization of \(\Delta=0\) is finite and birational over the
   selected reduced ramified target stratum.
4. **Rank-one oriented link.**  Every declared source/target localization
   link is a saturated rank-one affine modification.  Its unit lattice
   gives the orientation \(+1\) or \(-1\); the orientation is not chosen
   from a desired normal form.
5. **Complete contact monoid.**  The presentation includes the integral
   contact monoid of all boundary parameters, Rees generators, and
   reconstruction denominators.  It records the homomorphism to the
   source and target valuation lattices, not only the cone or a Hilbert
   basis of value vectors.
6. **Finite filtered presentation.**  The chart algebra and reconstruction
   module have a finite Khovanskii/SAGBI presentation for the contact
   filtration, with a saturated Rees module.  Every possible polar divisor
   of every reconstruction coordinate occurs in the declared ledger.
7. **Primitive noncontraction.**  The marked conormal generates after
   height-one saturation and at closed collisions, and its residue varies
   nontrivially along the normalized ramified stratum.

The presentation is required to reconstruct the original function-field
extension and the distinguished affine open.  It is not required to be
isomorphic to a weighted, cancellation, or quadratic-gauge chart.

### Definition 2.2 -- the extraction predicate

Write

\[
 \operatorname{MLE}(F)
\tag{2.5}
\]

when the canonical input of Section 1 admits a unique intrinsic marked
controlled-incidence presentation in the sense of Definition 2.1.

This predicate is strictly weaker than `CS`:

- it supplies a quotient flag, primitive marked root, controlled divisor,
  contact monoid, and a finite polynomiality certificate;
- it does not prescribe the degree or support of the horizontal incidence
  coordinate;
- it does not prescribe one of the three known reconstruction formulas;
- it does not assert polynomial left--right equivalence to a known map.

## 3. The two conjectural arrows

### Conjecture 3.1 -- intrinsic extraction

Assume that:

1. \(\pi_F\) is flat;
2. the canonical boundary has at most two divisorial target images and a
   complete prime ledger;
3. it has a unique intrinsically selected rational ramified stratum of
   puncture rank at most one;
4. every boundary link is a saturated rank-one affine modification;
5. the selected conormal is primitive and noncontracted; and
6. the associated contact/Rees data have the finite filtered presentation
   of Definition 2.1(5)--(6).

Then \(\operatorname{MLE}(F)\) holds.

The substantive assertion is the extraction of the quotient flag and
primitive marked-root coordinate.  Finite generation of a value semigroup
alone does not imply it.

### Conjecture 3.2 -- polynomial algebraization trichotomy

Let \(F\) satisfy \(\operatorname{MLE}(F)\), with critical normalization
\(\mathbb A^1\) or \(\mathbb G_m\) and at most two target boundary images.
Then:

\[
\begin{cases}
\operatorname{gdeg}(F)=3
 &\Longrightarrow
 F\text{ is stably left--right equivalent to the foundational cubic},\\
\operatorname{gdeg}(F)\ge4
 &\Longrightarrow
 F\text{ belongs to exactly one of the weighted tangent, cancellation,}\\
 &\hspace{3.8cm}\text{or root-engineered quadratic-gauge mechanisms.}
\end{cases}
\tag{3.1}
\]

The word “mechanism” is essential.  The assertion does not collapse the
stable moduli inside a family.

The expected intrinsic branch labels are:

\[
\begin{array}{c|c|c|c}
\text{mechanism}&\text{critical type}&
\text{Fitting-support rank}&
\text{boundary-contact index}\\ \hline
\text{weighted}&\mathbb A^1&\text{not needed}&\text{not needed}\\
\text{cancellation}&\mathbb G_m&1&mr(m+1)\\
\text{quadratic gauge}&\mathbb G_m&2&2.
\end{array}
\tag{3.2}
\]

For the foundational cubic, the last two rows both have support rank one
and contact index two, as required by the branch intersection.

## 4. Fixed-chart polynomial-\(Q\) reduction

The apparent non-fibre-affine enlargement lets the incidence coordinate
depend arbitrarily on the conjugate fibre variable \(Q\).  Keep

\[
 D=1-SQ+PS^2,\qquad B=Q+\beta(P,S),
\tag{4.1}
\]

and let \(C(P,S,Q)\) be any polynomial.  Rewrite it in the target fibre
coordinate:

\[
 \widehat C(P,S,B)
 =C(P,S,B-\beta(P,S)).
\tag{4.2}
\]

The change \((P,S,Q)\mapsto(P,S,B)\) has determinant one, and hence

\[
 \det\frac{\partial(P,B,C)}{\partial(P,S,Q)}
 =-\partial_S\widehat C,
\tag{4.3}
\]

where the derivative is taken at fixed \((P,B)\).  In these coordinates,

\[
 D=1-SB+S\beta(P,S)+PS^2.
\tag{4.4}
\]

Therefore the Keller identity

\[
 \det\frac{\partial(P,B,C)}{\partial(P,S,Q)}
 =-\lambda D
\tag{4.5}
\]

has the general polynomial solution

\[
 \widehat C
 =
 \lambda\int
 \left(1+S\beta(P,S)+PS^2\right)\,dS
 -\frac{\lambda}{2}BS^2
 +H(P,B).
\tag{4.6}
\]

The last term is exactly a triangular polynomial target shear.  Removing
it leaves the marked-line coordinate

\[
 \boxed{X=\frac{\lambda}{2}S^2.}
\tag{4.7}
\]

Thus arbitrary polynomial \(Q\)-dependence introduces no new mechanism on
the fixed reciprocal chart.  This is an all-degree statement in \(Q\), not
a degree-four truncation.  If the original target outputs are polynomial,
the shear by \(H(P,B)\) preserves polynomiality.

For example, in the quadratic ansatz

\[
 C=c_0(P,S)+c_1(P,S)Q+c_2(P,S)Q^2,
\tag{4.8}
\]

coefficient comparison gives

\[
 \partial_Sc_2=0,\qquad
 c_1-2c_2\beta=-\frac{\lambda}{2}S^2+f(P).
\tag{4.9}
\]

The target shear

\[
 C^\prime=C-c_2(P)B^2-f(P)B
\tag{4.10}
\]

recovers (4.7), agreeing with the general proof.

The exact regression is

```bash
.venv/bin/python scripts/verify_polynomial_in_q_incidence_reduction.py
```

The next genuinely open bounded ansatz must therefore alter the source
chart or mix the root with \(P\), rather than enlarge the second incidence
coordinate inside the fixed reciprocal chart.

## 5. Formal conormal integration

The primitive conormal data have two distinct parts.  The class itself is
transverse:

\[
 \tau\in H^0(D_\Delta,\mathcal N_\Delta).
\tag{5.1}
\]

Its residue coefficient is longitudinal:

\[
 \zeta_\tau\in k(E),
\tag{5.2}
\]

where \(E\) is the normalized ramified stratum.  Primitivity and
noncontraction give

\[
 \mathcal O(E)=k[\zeta_\tau]
 \quad\text{for }E\simeq\mathbb A^1,
\qquad
 [\zeta_\tau]\text{ generates }\mathcal O(E)^*/k^*
 \quad\text{for }E\simeq\mathbb G_m.
\tag{5.3}
\]

Thus the boundary root is already \(\zeta_\tau\), up to the orientation
and scalar normalization retained by the marked boundary link.  The formal
problem is to extend this residue root, with its prescribed first conormal
jet, into the marked chart.

### Proposition 5.1 -- unobstructed formal lifting

Let \(X\) be Noetherian, let \(D\subset X\) be an affine effective Cartier
divisor with ideal \(I\), and put

\[
 D_n=V(I^{n+1}).
\tag{5.4}
\]

Every function on \(D_n\) lifts to \(D_{n+1}\).  Consequently every
function on \(D\), together with any compatible prescribed finite jet,
admits a lift

\[
 \widehat s\in
 \Gamma(\widehat X_D,\mathcal O_{\widehat X_D})
 =
 \varprojlim_n\Gamma(D_n,\mathcal O_{D_n}).
\tag{5.5}
\]

If the initial function is a unit, it admits a formal unit lift.  At the
\((n+1)\)-st step the set of lifts, when nonempty, is a torsor under

\[
 H^0(D,I^{n+1}/I^{n+2}).
\tag{5.6}
\]

#### Proof

There is an exact sequence

\[
 0\longrightarrow I^{n+1}/I^{n+2}
 \longrightarrow\mathcal O_{D_{n+1}}
 \longrightarrow\mathcal O_{D_n}
 \longrightarrow0.
\tag{5.7}
\]

Because \(D\) is affine and \(I^{n+1}/I^{n+2}\) is quasi-coherent,

\[
 H^1(D,I^{n+1}/I^{n+2})=0.
\tag{5.8}
\]

The restriction map on global sections is therefore surjective.  Induction
gives (5.5), while its kernel gives (5.6).  A lift of a unit through a
nilpotent ideal is again a unit.  QED

If only the normalization of \(D\) is initially presented as affine, the
same conclusion applies once finite conductor descent has been checked:
a finite surjective image of an affine Noetherian scheme is affine.
Equivalently, one may work on the normalization and require at every order
that the lifted sections agree on the conductor equalizer.

Proposition 5.1 proves existence, not canonicity.  Even the primitive
conormal and a fixed first jet leave infinitely many higher choices.  For
example, on

\[
 X=\operatorname{Spec}k[x,y],\qquad D=(x),
\tag{5.9}
\]

the two formal parameters

\[
 x,\qquad \frac{x}{1-x}=x+x^2+x^3+\cdots
\tag{5.10}
\]

have the same primitive conormal class.  Likewise \(y\) and
\(y+x^2h(x,y)\) have the same restriction to \(D\) and the same first
normal jet.  The ambiguity is exactly the \(H^0\)-torsor in (5.6), not an
\(H^1\)-obstruction.

### Definition 5.2 -- formal incidence rigidity

Fix the oriented residue root, quotient flag, controlled divisor, and
allowed target-shear group.  Linearize the controlled-incidence and Keller
identities at Rees order \(n\):

\[
 L_n:
 H^0(D,I^n/I^{n+1})
 \longrightarrow \operatorname{Def}_n(F).
\tag{5.11}
\]

The marked presentation is **formally incidence-rigid** when, for every
\(n\ge2\),

1. the order-\(n\) obstruction lies in the image of \(L_n\); and
2. \(\ker L_n\) consists exactly of the order-\(n\) jets of declared target
   shears and marking-preserving gauges.

Under formal incidence rigidity, Proposition 5.1 gives a unique formal
marked-root **gauge class**.  The orientation and collision normalization
must then select a representative.  This is the precise uniqueness
statement that conormal primitivity alone cannot supply.

### 5.3 Formal effectivity and polynomiality

Finite generation of the boundary Rees algebra does not by itself
algebraize a formal function.  In (5.9), the Rees algebra of the principal
ideal \((x)\) is finitely generated, but the second expression in (5.10)
has a pole on the affine divisor \(x=1\).  More generally the completion
contains arbitrary power series that are neither rational nor polynomial.

The algebraization step must therefore be split into:

\[
 \boxed{
 \text{formal gauge class}
 \Longrightarrow
 \text{rational effectivity}
 \Longrightarrow
 \text{polynomial membership}.
 }
\tag{5.12}
\]

Here **rational effectivity** means that the selected formal class is the
completion of an element \(s\in k(\mathbb A^3)\).  Once such an \(s\) is
known, polar completeness excludes an undeclared divisor such as \(x=1\)
in (5.10), and saturation/Rees strictness upgrades divisorial
nonnegativity to membership in the actual reconstruction module.  Finite
Khovanskii/Rees generation makes these last tests finite; it does not
produce the rational representative.

This gives the corrected conditional extraction theorem:

> **Formal marked-root extraction criterion.**  Suppose the canonical
> boundary data determine an oriented primitive residue root and an affine
> Cartier formal neighborhood.  If the controlled-incidence equations are
> formally incidence-rigid, their unique formal gauge class has a rational
> representative, and the complete reconstruction ledger is polar-complete
> and Rees-strict, then the marked-root coordinate algebraizes.  Together
> with the intrinsic quotient flag, this proves \(\operatorname{MLE}(F)\).

The first new proof target is consequently not another unrestricted
coefficient search.  It is the homogeneous calculation

\[
 \boxed{
 \ker L_n
 =
 \{\text{target-shear and marking-preserving gauge jets}\}
 \quad(n\ge2).
 }
\tag{5.13}
\]

After (5.13), rational effectivity is the only genuinely nonformal step.

## 6. Relation to the surrounding program

The finite filtered-presentation requirement is developed in
[`FINITE_VALUATION_ALGEBRAIZATION.md`](FINITE_VALUATION_ALGEBRAIZATION.md).
The root-preserving fibre-affine classification is
[`../verified/INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md`](../verified/INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md).
The first root-only Möbius escape and its order-two pole obstruction are
recorded in
[`ESCAPE_THREE_SUSPENSION_FAMILIES.md`](ESCAPE_THREE_SUSPENSION_FAMILIES.md).

Together they leave one smallest uncontrolled chart class: a
\(P\)-dependent root change coupled to a non-fibre-affine determinant-one
source rechart.

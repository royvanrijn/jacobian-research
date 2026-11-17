# Puncture-rank frontier for Keller boundaries

This note audits the proposed assertion that a uniquely selected rational
ramified boundary of a nonproper Keller map has unit rank at most one.  It
separates a formulation issue from an all-degree direct-chart obstruction.

The outcome is:

1. the literal statement with `E=A^1` or `G_m` is not a statement about the
   canonical boundary prime used in this repository: that prime is a
   surface;
2. the natural corrected examples are
   `A^2` and `G_m x A^1`, whose unit ranks are zero and one;
3. finite normalization, tame generic ramification, rationality, and absence
   of a second boundary prime over the same target divisor do not by
   themselves turn Jelonek's uniruledness theorem into a puncture bound;
4. nevertheless, the most direct two-unit generalization of the reciprocal
   chart is impossible in every geometric degree.

The last item is an exact all-degree theorem.  A bounded integer-lattice
census is retained as a reproducible screen, but the coefficient obstruction
does not depend on that bound.

Work over an algebraically closed field `k` of characteristic zero.

## 1. The dimensional correction

For a dominant quasi-finite map

\[
 F:\mathbb A^3\longrightarrow\mathbb A^3
\]

the canonical finite normalization

\[
 \bar X_F=\operatorname{Norm}_{\mathbb A^3}k(\mathbb A^3)
\]

is three-dimensional.  Every prime component of
`\bar X_F - A^3` has codimension one, hence dimension two.  Thus a canonical
boundary prime `E` cannot literally be isomorphic to `A^1` or `G_m`.

In the two established suspension families, the relevant normalizations
have the form

\[
 E_{\rm wt}\simeq\mathbb A^1\times\mathbb A^1,\qquad
 E_{\rm rec}\simeq\mathbb G_m\times\mathbb A^1,
\]

or have the same unit lattices after the intrinsic affine-line parameter is
retained.  Consequently

\[
 \operatorname{rank}\mathcal O(E_{\rm wt})^*/k^*=0,
 \qquad
 \operatorname{rank}\mathcal O(E_{\rm rec})^*/k^*=1.
\]

A corrected puncture statement must therefore do one of the following:

- state only the unit-rank inequality for the boundary surface;
- canonically extract an `A^1`-fibration `E -> C` and state the conclusion
  for the relative curve `C`; or
- assume from the start that `E=C x A^1` and conclude
  `C=A^1` or `G_m`.

The last conclusion also needs normality (equivalently smoothness for a
curve here).  A geometrically integral rational affine curve of unit rank
zero need not itself be `A^1`; singular rational curves give immediate
counterexamples to that inference.

## 2. Why the suggested general inputs do not yet prove the bound

Two of the proposed hypotheses carry no extra force in the present setting.
For a generically finite polynomial self-map of affine space, normalization
in the finite function-field extension is finite because finitely generated
algebras over a field are excellent (in particular Nagata).  Generic DVR
ramification is tame automatically in characteristic zero.  The substantive
hypotheses are therefore the intrinsic selection, rationality, and the
no-additional-boundary condition; the last must specify whether it excludes
only boundary primes over the same target divisor or also transverse affine
divisors meeting the selected boundary.

Jelonek--Lasoń prove that the nonproperness hypersurface of a generically
finite polynomial map from affine space is covered by polynomial curves of
bounded degree.  For a curve, being covered by polynomial curves would force
a one-puncture normalization.  Here, however, the nonproperness set is a
surface and the selected canonical boundary is another surface finite over
it.  Uniruledness of the target surface does not bound the unit rank of a
curve factor in the source surface.

The generic DVR identity

\[
 \sum_{E_i\mid Z}e(E_i/Z)f(E_i/Z)=\operatorname{gdeg}(F)
\]

does not repair this gap.  The residue degree `f(E/Z)` is the generic degree
of one surface over another.  It is not the degree of a chosen punctured
curve over `A^1`.  In particular, a three-punctured relative curve need not
consume three residue sheets.

Likewise, tame log crepancy

\[
 K_{\bar X_F}+D_F=\pi^*(K_{\mathbb A^3}+B_F)
\]

is a codimension-one identity.  After adjunction it gives the expected
Riemann--Hurwitz identity on any genuinely extracted curve quotient, but it
does not construct that quotient or bound its punctures.

The missing general-theorem input is therefore precise:

> **Relative-curve extraction problem.**  Canonically extract from the
> selected boundary surface a curve quotient `E -> C` such that independent
> units of `E` descend to puncture units of `C`, and prove that the Keller
> reconstruction gives at most one independent puncture character.

Without this input, the proposed theorem is stronger than the listed tools.

## 3. Universal two-center reciprocal chart

The counterexample strategy in the question has a clean universal form.
Let `f in k[y]` have degree `d>=1`, let `r>=1`, and put

\[
 A=1+xf(y),\qquad
 B=A^{r+1}z+g(y,A),
\]

\[
 s=\frac{x}{A},\qquad P=AB,\qquad Q=y+xB,\qquad
 Y=Q-sP.
\]

Then identically

\[
 Y=y,\qquad
 D=1-sf(Y)=A^{-1}.
\]

The normal critical boundary of the coordinate-preserving core has ring

\[
 \mathcal O(E)
 =k[Y,P,f(Y)^{-1}].
\]

If `f` has `ell` distinct roots, then

\[
 \operatorname{rank}\mathcal O(E)^*/k^*=\ell.
\]

Thus two distinct roots give exactly the desired rank-two surface

\[
 \left(\mathbb P^1\setminus\{\alpha,\beta,\infty\}\right)
 \times\mathbb A^1.
\]

The chart determinant is independent of `f` and `g`.  First,

\[
 \det\frac{\partial(s,P,Q)}{\partial(x,y,B)}=-A^{-1}.
\]

Since `partial_z B=A^{r+1}`,

\[
 \boxed{
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}=-A^r=-D^{-r}.
 }
\tag{3.1}
\]

Define the core output

\[
 R=C\int_0^s\{1-tf(Q-Pt)\}^r\,dt.
\tag{3.2}
\]

Equivalently, this begins with the marked-root incidence equation

\[
 \boxed{
 \Psi_{f,r}(s;P,Q,R)
 =C\int_0^s\{1-tf(Q-Pt)\}^r\,dt-R=0.
 }
\tag{3.3}
\]

Here `s` is the marked root.  The normalization parameter below is
`Y`—it may be renamed `T` when writing the selected curve as
`\(\mathbb P^1_T\setminus\{0,1,\infty\}\)`.

Its root derivative is

\[
 \partial_s\Psi_{f,r}=CD^r.                            \tag{3.4}
\]

Thus (3.1) is precisely the divisor ledger

\[
 \operatorname{div}(J_\alpha)
 +r\operatorname{div}(D\circ\alpha)=0
 =F^*\operatorname{div}(J_\beta),                     \tag{3.5}
\]

with `beta=id`.  In particular, the integer-lattice determinant equation
is solved before any coefficient search.  Then `R_s=CD^r`, so (3.1) gives

\[
 \det\frac{\partial(P,Q,R)}{\partial(x,y,z)}=-C
\]

whenever (3.2) is polynomial.  Its generic inverse degree is

\[
 N=r(d+1)+1.
\tag{3.6}
\]

This is the smallest direct realization of
`D=1-sU_1^aU_2^b`: after affine normalization of the two roots,

\[
 f(y)=y^a(y-1)^b,\qquad a,b\ge1.
\]

On `D=0`, put `Y=Q-sP`.  Then

\[
 s=\frac1{Y^a(Y-1)^b},\qquad
 Q=Y+\frac{P}{Y^a(Y-1)^b},
\]

so the normalized controlled divisor has coordinate ring

\[
 k[Y,P,Y^{-1},(Y-1)^{-1}]
\]

and relative curve
\(\mathbb P^1\setminus\{0,1,\infty\}\).  Its puncture character lattice is

\[
 \Lambda=
 \{(v_0,v_1,v_\infty)\in\mathbb Z^3:
   v_0+v_1+v_\infty=0\},
\tag{3.7}
\]

with saturated basis

\[
 \operatorname{div}(Y)=(1,0,-1),\qquad
 \operatorname{div}(Y-1)=(0,1,-1).                   \tag{3.8}
\]

Hence the unit rank is two and the affine class group is zero.  The
degree-zero relation in (3.7) also shows why three independently variable
puncture characters are impossible: three punctures supply exactly two
independent characters.  The pair (3.8) is the maximal independent
reconstruction ledger and jointly records all three punctures.  The
controlled character is

\[
 \operatorname{div}(f)=(a,b,-a-b),                   \tag{3.9}
\]

primitive exactly when `gcd(a,b)=1`.

For a bounded ledger census, write two proposed character rows in the basis
(3.8).  A coefficient matrix of determinant zero is rejected for
insufficient unit rank; determinant of absolute value greater than one
leaves a finite nonsaturated class-lattice quotient; determinant
\(\pm1\) is a saturated basis.  With coefficients in `[-2,2]`, the exact
counts are respectively

\[
 129,\qquad392,\qquad104.
\tag{3.10}
\]

There is only one surviving lattice up to integral row change: (3.8).
For `1<=r,a,b<=4`, there are 44 primitive positive controlled ledgers
`(r,a,b)`.  All 44 pass the unit-rank and class-group screen.  Thus the
determinant and lattice ledgers pose no obstruction; polynomiality is the
entire issue.

## 4. The boundary-moment obstruction

Write

\[
 g_0(y)=g(y,0).
\]

In (3.2), substitute `t=sv`.  The identities above give

\[
 R=
 \frac{Cx}{A^{r+1}}
 \int_0^1
 \left[
 A-xv f\bigl(y+xB(1-v)\bigr)
 \right]^r\,dv.
\tag{4.1}
\]

If `R` is polynomial, the numerator in (4.1) is divisible by `A^{r+1}`.
Reducing its constant term modulo `A`, where

\[
 x=-f(y)^{-1},\qquad B=g_0(y),
\]

gives the necessary identity

\[
 \boxed{
 \mathcal M_{f,r}(y,g_0):=
 \int_0^1v^r
 f\left(
 y-\frac{g_0(y)}{f(y)}(1-v)
 \right)^r\,dv=0.
 }
\tag{4.2}
\]

After clearing denominators, (4.2) is a polynomial equation in
`k[y,g_0]`.  It is only the first cancellation equation, so failure of
(4.2) rules out every higher jet `g(y,A)` at once.

### Theorem 4.1 -- all-degree two-center reciprocal no-go

Let

\[
 f(y)=y^a(y-1)^b,\qquad a,b\ge1,
\]

and form (3.1)--(3.2), with any `r>=1`.  There is no polynomial `g(y,A)`
for which `R` is polynomial.  Consequently this entire direct two-center
reciprocal family contains no three-puncture Keller map in any geometric
degree.

#### Proof

Suppose first that (4.2) holds, and put

\[
 c(y)=\frac{g_0(y)}{f(y)}.
\]

At `y=0`, if `c` had a pole, the highest polar term of (4.2) would be

\[
 (-c)^{r(a+b)}
 \int_0^1v^r(1-v)^{r(a+b)}\,dv,
\]

whose beta integral is nonzero in characteristic zero.  Hence
`y^a|g_0`.  The same calculation at `y=1` gives
`(y-1)^b|g_0`.  Therefore `c in k[y]`.  If `deg(c)>1`, the identical leading
term at infinity, now with the leading coefficient of `c`, cannot vanish.
Thus

\[
 \boxed{c(y)=\lambda y+\mu.}                           \tag{4.3}
\]

Set

\[
 n=ar,\qquad m=br,\qquad p(y)=y^n(y-1)^m,
\]

and replace `v` by `1-t`.  The moment equation becomes

\[
 \int_0^1(1-t)^r p\bigl(y-tc(y)\bigr)\,dt=0.          \tag{4.4}
\]

If `lambda=0`, the coefficient of `y^(n+m)` in (4.4) is
`1/(r+1)`, a contradiction.  Assume `lambda!=0`, put

\[
 \eta=-\frac{\mu}{\lambda},\qquad z=y-\eta,
\]

and expand

\[
 p(\eta+z)=\sum_k p_kz^k.
\]

Since `c(y)=lambda*z`, equation (4.4) diagonalizes:

\[
 \sum_k p_k I_k(\lambda)z^k=0,\qquad
 I_k(\lambda)=
 \int_0^1(1-t)^r(1-\lambda t)^k\,dt.                  \tag{4.5}
\]

The constant coefficient gives

\[
 \frac{p(\eta)}{r+1}=0,
\]

so `eta=0` or `eta=1`.  If `eta=0`, then
`p(z)=z^n(z-1)^m` has nonzero coefficients in every degree from `n` through
`n+m`; in particular, (4.5) forces
`I_n(lambda)=I_(n+1)(lambda)=0`.  If `eta=1`, then
`p(1+z)=(1+z)^n z^m` instead forces
`I_m(lambda)=I_(m+1)(lambda)=0`.

Both conclusions are impossible.  Direct integration by parts gives, for
every `k>=0`,

\[
 \boxed{
 (r+k+2)I_{k+1}(\lambda)
 +(k+1)(\lambda-1)I_k(\lambda)=1.
 }
\tag{4.6}
\]

Two consecutive eigenvalues can therefore never vanish simultaneously.
This contradicts (4.4) and proves that the first necessary divisibility
condition already fails.  Hence no higher `A`-jet in `g(y,A)` can make `R`
polynomial. QED

### Bounded regression

The all-degree proof subsumes the earlier degree-four through degree-seven
factor screen.  The possible triples in that range are:

\[
\begin{array}{c|c}
r& (a,b)\\ \hline
1&(1,1),\\
1&(1,2),(2,1),\\
1&(1,3),(2,2),(3,1),\\
1&(1,4),(2,3),(3,2),(4,1),\\
2&(1,1).
\end{array}
\tag{4.7}
\]

For the first and smallest case no computer algebra is needed.  With
`r=a=b=1`, equation (4.2) is equivalent to

\[
 G^2-2y(y-1)(2y-1)G+6y^3(y-1)^3=0.
\tag{4.8}
\]

Its discriminant is

\[
 4y^2(y-1)^2(1+2y-2y^2),
\]

which is not a square in `k[y]`.  Hence (4.8) has no polynomial root.
The checker retains this calculation and the rational factorizations of all
eleven low-degree moments as regressions.  Absolute factorization is no
longer a proof dependency.

### Scope

The theorem allows:

- arbitrary positive multiplicities at the two boundary centers;
- every geometric degree in this chart;
- an arbitrary polynomial center `g(y,A)`, not merely a monomial or a
  bounded coefficient ansatz; and
- arbitrary higher `A`-jets and arbitrary `z`-terms, because the obstruction
  occurs already modulo `A`.

It does not exclude:

- a nonmultiplicative rational source chart;
- two independent reconstruction variables;
- a nontrivial target Jacobian ledger;
- a critical boundary not presented by `1-sf(Q-Ps)`; or
- a construction in which the two puncture units do not lift from two
  affine-linear boundary centers.

Those are exactly the escape routes already suggested by the
controlled-boundary programme.

## 5. What this says about the conjecture

The proposed rank bound should be split into two assertions.

### Corrected structural conjecture

Let `E` be the canonically selected normal ramified boundary **surface**.
Assume the canonical package also extracts an `A^1`-fibration

\[
 E\longrightarrow C
\]

onto a smooth rational affine curve, identifies
`\mathcal O(E)^*/k^*` with `\mathcal O(C)^*/k^*`, and supplies a single
primitive reconstruction character.  Under the saturated-link,
boundary-monotonicity, no-extra-boundary, and closed-point conductor
conditions, conjecture that

\[
 \operatorname{rank}\mathcal O(E)^*/k^*\le1.
\]

Then `C=A^1` or `G_m`; it is `C`, not the surface `E`, that has this
isomorphism type.

### Proven direct-chart theorem

For the direct reciprocal realization with two affine-linear centers,
the conjectural conclusion holds in every degree for the stronger reason
that polynomiality already fails at its first boundary moment.

This closes every coefficient search inside the same reciprocal chart.
A meaningful counterexample search must leave it, most economically by
introducing a target ledger or a second reconstruction variable.

## 6. What remains after the no-go theorem

The requested four-stage construction pipeline now terminates at stage
three for the direct chart:

1. the puncture lattice enumeration leaves the unique saturated rank-two
   lattice (3.8), with 44 primitive positive ledgers in the displayed
   coefficient box;
2. all those ledgers pass the unit-rank and class-group tests;
3. the three local conditions at `0`, `1`, and infinity reduce the numerator
   to `g_0=f(lambda*y+mu)`, and (4.6) rejects every coefficient choice;
4. consequently no candidate reaches polynomiality, so there is no
   polynomial map on which to test a constant determinant or complete
   collision.  The rational determinant is nevertheless the constant `-C`
   wherever the formulas are defined.

This is a no-go theorem for one complete all-degree ledger family, not a
classification of all low-complexity three-puncture diagrams.  The existing
[double-incidence core](../extended-geometry/THREE_PUNCTURE_DOUBLE_INCIDENCE_CORE.md)
uses two reconstruction variables and solves the determinant and
normalization halves independently.  Its
[nonlinear completion frontier](../extended-geometry/THREE_PUNCTURE_NONLINEAR_COMPLETION_FRONTIER.md)
proves two dimension-free rank-drop gates and rejects 80 coupled
degree-at-most-three `A^6` skeletons.  Those are the first remaining charts
after the present theorem; neither note claims a Keller map.

## 7. Exact verification

Run

```bash
.venv/bin/python scripts/verify_puncture_rank_frontier.py
```

The checker verifies:

1. the universal chart determinant;
2. the rank-two boundary ring for two distinct roots;
3. the bounded integer-lattice counts (3.10) and the 44 primitive positive
   controlled ledgers;
4. the degree list (4.7);
5. the hand quadratic and its nonsquare discriminant;
6. the consecutive-eigenvalue recurrence (4.6) for a broad exact regression
   range; and
7. the rational factorization data for every degree-four through
   degree-seven boundary moment.

The proof of Theorem 4.1 is the displayed symbolic argument, not an inference
from the finite regression range.  Singular is no longer required.

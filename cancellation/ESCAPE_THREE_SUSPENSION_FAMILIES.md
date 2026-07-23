# Escape from the three known suspension families

## 0. Scope and first result

The degree-four incidence theorem leaves four real gaps: non-fibre-affine
source recharts, root-changing charts, multi-boundary ledgers, and incidences
with more than one marked-line coefficient.  This note formulates the
counterexample search so that polynomiality is tested rather than inferred
from determinant cancellation, and closes the smallest root-changing
subcase.

The first exact result is negative:

> A non-affine Möbius change of the marked root, equipped with its canonical
> determinant-one lift of the quadratic reciprocal chart, produces a
> polynomial controlled divisor and generically a cubic horizontal
> coordinate.  Nevertheless its first incidence coordinate has an
> unavoidable order-two pole on the original reciprocal boundary.  Hence it
> does not give a polynomial Keller map.

This does not address a root change mixed with `P`, a non-fibre-affine chart,
or a target divisor ledger.  Those are the next bounded searches.

## 1. The equations that a coefficient search must impose

Let a rational source chart be

\[
 \alpha=(P,S,B):\mathbb A^3_{x,y,z}\dashrightarrow
 \mathbb A^3_{P,S,B}
\]

and put

\[
 C=Y(S,P,B)-B X(S,P,B).
\]

At fixed `(P,B)`,

\[
 D:=\partial_S C=Y_S-BX_S,
\]

and the chain rule gives

\[
 \boxed{\det D(P,B,C)=-D\,\det D(P,S,B).}              \tag{1.1}
\]

Thus a normalized search must solve the rational identity

\[
 D\,\det D\alpha=1                                    \tag{1.2}
\]

together with all of the following non-Jacobian conditions:

1. `P`, `B`, and `C` pull back to polynomials in `(x,y,z)`;
2. the graph of `alpha` is birational, certified by reconstruction
   identities or a saturated graph ideal;
3. the inverse equation
   \[
    E(T;P,B,C)=Y(T,P,B)-B X(T,P,B)-C
   \]
   has exact degree four and is generically separable;
4. every declared denominator divisor is recorded in the source/target
   ledger, including the divisor at infinity of a nonmonomial denominator;
5. the reconstruction open is saturated before components are compared.

Allowing arbitrary `B`-dependence in both `X` and `Y` introduces a gauge:
only their combination `C=Y-BX` occurs in (1.1).  Indeed any `C` admits the
presentation `X=0,Y=C`.  Consequently a finite search must either keep
`C` affine-linear in `B` with `X=X(S,P)`, or place a direct support bound on
`C` and forget the artificial pair `(X,Y)`.  Without this normalization,
"degree of the horizontal marked-line coordinate" is not invariant and the
coefficient system contains large tautological components.

## 2. The smallest root-changing chart

Start with the quadratic reciprocal chart

\[
 D=1-SQ+PS^2,\qquad
 \det D(P,S,Q)=D^{-1}.
\]

Every birational self-map of the marked-root line over the coefficient
field is Möbius:

\[
 T=\phi(S)=\frac{aS+b}{cS+d},\qquad
 \Delta=ad-bc\ne0.                                    \tag{2.1}
\]

Its determinant-one lift in the conjugate coordinate is

\[
 \widetilde Q=\frac{Q}{\phi'(S)}+H(P,T).               \tag{2.2}
\]

The unique shear which removes the old `PS^2` term from the controlled
divisor is

\[
 H=-\frac{PS}{\phi'(S)}.
\]

It gives the especially simple identities

\[
 \widetilde Q=\frac{Q-PS}{\phi'(S)},\qquad
 D=1-\kappa(T)\widetilde Q,                            \tag{2.3}
\]

where

\[
 \boxed{\kappa(T)=
 \frac{(dT-b)(a-cT)}{\Delta}.}                         \tag{2.4}
\]

For `cd!=0`, `kappa` has degree two.  The incidence equation
`X_T=lambda*kappa` therefore has a horizontal coordinate of exact degree
three.  This is a genuine case omitted by the root-preserving theorem: the
controlled divisor is polynomial and, generically, nonmonomial in the new
root.

## 3. Exact polynomiality obstruction

On the original affine source write

\[
 t=1+xy,\qquad S=x/t,\qquad Q-PS=y,\qquad P=tq.
\]

Equations (2.1)--(2.3) become

\[
 T=\frac{ax+bt}{cx+dt},\qquad
 \widetilde Q=
 \frac{y(cx+dt)^2}{\Delta t^2}.                        \tag{3.1}
\]

At the generic point of `t=0`, with `c!=0`, the new root is regular:

\[
 T\equiv a/c,\qquad P\equiv0.
\]

But

\[
 \left.t^2\widetilde Q\right|_{t=0}
 =\frac{c^2x^2y}{\Delta}\ne0.                         \tag{3.2}
\]

Every polynomial incidence shear `beta(P,T)` is regular at this valuation,
so it has zero `t^2` residue and cannot cancel (3.2).  Hence

\[
 \widetilde B=\widetilde Q+\beta(P,T)
\]

has an exact order-two pole.  The cubic horizontal incidence is a valid
rational constant-Jacobian suspension, but it does not extend to a
polynomial map of affine three-space.

If `c=0`, the root change is affine and stays in the already classified
quadratic chart.  If `c!=0` but `cd=0`, the same pole obstruction applies,
although `kappa` drops in degree.  Thus the entire Möbius root-only
rechart class is eliminated; in particular its generic cubic member is not
a fourth suspension family.

## 4. What the next finite search should be

The calculation isolates the necessary escape: the root change must mix
with `P` (or with an additional primitive coordinate) so that the new root
itself carries enough boundary valuation to cancel (3.2).  A useful next
ansatz is therefore not a larger polynomial `Y`; it is

\[
 T=\frac{U_0(S)+P\,U_1(S)}
         {V_0(S)+P\,V_1(S)},                           \tag{4.1}
\]

with `deg U_i,deg V_i<=2`, coupled to

\[
 R=A_0(S)P+A_1(S)Q+A_2(S),\qquad
 \widetilde Q=B_0(S)P+B_1(S)Q+B_2(S),                 \tag{4.2}
\]

initially with coefficient degrees at most two.  Impose:

1. the determinant of `(R,T,Q_tilde)` relative to `(P,S,Q)` is one;
2. the transformed `D` is affine in `Q_tilde` and polynomial in `(R,T)`;
3. its `Q_tilde` coefficient has degree at most three in `T`;
4. the incidence outputs pull back polynomially to `(x,y,z)`;
5. the graph and inverse equation saturations exclude denominator,
   degree-drop, and inseparable components.

The search should be stratified by denominator support.  A single
irreducible denominator tests root-changing/non-fibre-affine charts.
Two declared factors plus the extra three-term derivative factor from the
known two-boundary calculation test a three-boundary ledger.  A genuinely
nonmonomial irreducible denominator should be its own stratum, not replaced
by a monomial exponent vector.

Known families should be removed only after solutions are reconstructed.
Formula-level ideal quotients are unsafe because polynomial left--right
orbits and boundary degenerations are not single closed coefficient
components.  The reliable filters are the normalized critical curve
(`A^1` versus `G_m` or higher puncture rank), the complete boundary ledger,
the marked Fitting divisor, and the scheme-theoretic boundary-contact
index.

The old-boundary contacts of (4.1) are classified before coefficient
elimination in
[`FOURTH_SUSPENSION_VALUATION_FAN.md`](FOURTH_SUSPENSION_VALUATION_FAN.md).
For coefficient degrees at most two, only six effective-degree pairs
survive; the sole exceptional derivative-cancellation row fails its initial
residue-shape test unless the proposed root degenerates to a function of
`P`.

Imposing that `T` is itself a primitive root variable collapses those six
rows further; see
[`FOURTH_SUSPENSION_PRIMITIVE_ROOT_REDUCTION.md`](FOURTH_SUSPENSION_PRIMITIVE_ROOT_REDUCTION.md).
The one-prime affine and reciprocal-scaling cases fail respectively at the
controlled-divisor coefficient and the forced incidence residue.  A
unimodular nonmonomial Möbius denominator requires a nonpolynomial
third-order shear.  The smallest remaining search therefore has either a
determinant-supported multi-prime ledger or one additional primitive
reconstruction variable.

## Exact regression

Run

```bash
.venv/bin/python scripts/explore_root_changing_incidence_degree_four.py
```

The checker verifies the determinant-one lift, the transformed controlled
divisor, the cubic horizontal primitive, and the unavoidable order-two
source pole.

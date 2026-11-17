# Contracted-divisor Smith classification and the cubic E8 packet

> **Status.**  Generic local theorem and minimal cubic-jet normal-form
> theorem.  If a source boundary divisor `T` is contracted to a smooth
> target point and `(x,y)` has common order `h>=1` along `T`, then the
> logarithmic differential has generic Smith exponents
> `(h,delta-h)`, where `delta` is its determinant order along `T`.
> In particular, both exponents are positive: a contracted divisor is
> necessarily a positive-dimensional noncyclic `Fitt_1` packet.
>
> For the minimal cubic E8 attachment, transverse cusp divisibility forces
> the first two normal jets.  After removing the common contracted factor
> `t`, the logarithmic matrix is equivalent to
> `[[r,0],[t^2,r^2]]`.  Thus the unsaturated matrix is
> `t*[[r,0],[t^2,r^2]]`, its generic Smith form on `T` is `diag(t,t)`, and
> its saturated residual has determinant `r^3` and isolated
> `Fitt_1=(r,t^2)`.  The residual cokernel contains `R/(r^3)`, but its
> quotient is `R/(r^2,t^2)` of length four.  This completely classifies the
> minimal local module.  It does not yet compute the global degree of the
> universal rank-two quotient on `T` or decide whether that quotient glues
> to the terminal `A_6` cover.

The valuation arithmetic, forced cubic jets, matrix operations, Fitting
ideals, and quotient length are replayed by
[`verify_log_contracted_divisor_smith_classification.py`](../scripts/verify_log_contracted_divisor_smith_classification.py).

## 1. Generic theorem over the contracted divisor

Let `T=(t=0)` be a reduced source boundary component and let `K=k(T)` be
its function field.  Complete at its generic point:

\[
 R_T=K[[t]].                                      \tag{1.1}
\]

Suppose a morphism to a smooth target point has parameters

\[
 x=t^hX,qquad y=t^hY,qquad h\ge1,              \tag{1.2}
\]

where at least one of `X,Y` is a unit.  Use one tangential derivation
`partial_z` and the logarithmic normal derivation `t partial_t` on the
source.  Every entry of the logarithmic differential matrix is divisible by
`t^h`.  Conversely, characteristic zero gives

\[
 t\partial_t(t^hX)=t^h(hX+t\partial_tX),          \tag{1.3}
\]

so one entry has order exactly `h`.  Therefore the first Smith exponent is
`h`.

If

\[
 \delta=\operatorname{ord}_T(\det\Theta),         \tag{1.4}
\]

the elementary-divisor theorem over the DVR `K[[t]]` gives

\[
 \boxed{
 \operatorname{Smith}_{\eta_T}(\Theta)
 =\operatorname{diag}(t^h,t^{\delta-h}).}         \tag{1.5}
\]

Since all four entries have order at least `h`, necessarily
`delta>=2h`.  Hence both exponents in (1.5) are positive.  Equivalently,

\[
 \boxed{
 T\subset V(\operatorname{Fitt}_1(\operatorname{coker}\Theta)).} \tag{1.6}
\]

This proves that contraction to a target point and generic cyclicity are
mutually exclusive.  The isolated-`Fitt_1` positivity theorem can be applied
only after the common contracted factor has been removed.

## 2. Forced jets for cubic E8 inertia

Work in `R=k[[r,t]]`.  The affine E8 component is `E=(r=0)` and the
contracted boundary component is `T=(t=0)`.  Normalize the cusp
parametrization to

\[
 x(0,t)=t^3,qquad y(0,t)=t^5.                  \tag{2.1}
\]

The minimal contracted order is one, so after a boundary-preserving source
change write

\[
 x=t^3+a rt+b r^2t+O(r^3t),qquad a\ne0.        \tag{2.2}
\]

Impose cubic transverse contact:

\[
 y^3-x^5\in(r^3).                               \tag{2.3}
\]

The coefficients of `r` and `r^2` in (2.3) vanish successively.  They force

\[
\boxed{
 y=t^5+\frac53a rt^3
 +r^2\left(\frac53b t^3+\frac59a^2t\right)
 +c r^3t+O(r^4t).}                              \tag{2.4}
\]

No choice is involved in the two displayed transverse followers.  Direct
differentiation gives

\[
 J_{(r,t)}(x,y)
 =-\frac19r^2t\left(5a^3+O(r,t^2)\right).       \tag{2.5}
\]

Thus the logarithmic determinant has multiplicity three on `E` and two on
`T`, while `a!=0` makes the remaining factor a unit.

## 3. Saturation by the contracted divisor

At the SNC node use `(r partial_r,t partial_t)`.  Equation (2.2)--(2.5)
show that

\[
 \Theta=t\Phi,                                  \tag{3.1}
\]

where

\[
 \det\Phi=r^3\cdot\text{unit},qquad
 \operatorname{Fitt}_1(\operatorname{coker}\Phi)=(r,t^2). \tag{3.2}
\]

Indeed, the upper-left entry of `Phi` is `r` times a unit.  Use it to clear
the upper-right entry.  Subtract its unit multiple from the lower-left
entry; the remainder is `t^2` times a unit.  The determinant condition then
forces the lower-right entry to be `r^2` times a unit.  Scaling rows and
columns gives the normal form

\[
 \boxed{
 \Phi\sim
 \begin{pmatrix}r&0\\t^2&r^2\end{pmatrix},
 \qquad
 \Theta\sim
 t\begin{pmatrix}r&0\\t^2&r^2\end{pmatrix}.}    \tag{3.3}
\]

The equivalence uses invertible source and target basis changes over `R`.
Consequently

\[
\begin{aligned}
 \operatorname{Fitt}_0(\operatorname{coker}\Theta)&=(r^3t^2),\\
 \operatorname{Fitt}_1(\operatorname{coker}\Theta)&=t(r,t^2),\\
 \operatorname{Fitt}_0(\operatorname{coker}\Phi)&=(r^3),\\
 \operatorname{Fitt}_1(\operatorname{coker}\Phi)&=(r,t^2).
\end{aligned}                                    \tag{3.4}

This separates the two geometries that were previously mixed together:

- `t` is the common rank-two contracted packet;
- `Phi` is generically cyclic on `E` and has one isolated cusp defect.

## 4. Exact module filtration

Let

\[
 M=\operatorname{coker}\Theta,qquad
 M_{\rm sat}=\operatorname{coker}\Phi.           \tag{4.1}
\]

Factoring through the subbundle `F(-T)` gives the exact sequence

\[
 \boxed{
 0\longrightarrow M_{\rm sat}\longrightarrow M
 \longrightarrow F/F(-T)\longrightarrow0.}      \tag{4.2}

Locally, the final term is `(R/(t))^2`.  It is the entire
positive-dimensional split packet; it must not be recorded as a point
length.

For the saturated normal form in (3.3), the class of the first target basis
vector has annihilator `(r^3)`.  Killing that class leaves the second target
generator with relations `(r^2,t^2)`.  Therefore

\[
 \boxed{
 0\longrightarrow R/(r^3)\longrightarrow M_{\rm sat}
 \longrightarrow R/(r^2,t^2)\longrightarrow0.}  \tag{4.3}

In particular,

\[
 \ell R/(r,t^2)=2,qquad
 \ell R/(r^2,t^2)=4.                            \tag{4.4}

The isolated Fitting colength gives only the lower bound two; the actual
cyclic-submodule quotient in the integrable cubic normal form has length
four.  This distinction is invisible to the determinant and to generic
branch Smith data.

## 5. Global charge of the split factor

Globally, if `F=Omega_X^1(log D_X)`, the last term of (4.2) is

\[
 F\otimes O_T.                                   \tag{5.1}
\]

Its second Chern character is

\[
 \operatorname{ch}_2(F\otimes O_T)
 =c_1(F)\mathbin\cdot T-T^2
 =L_X\mathbin\cdot T-T^2.                      \tag{5.2}

For a rational SNC boundary component of valency `v`, logarithmic
adjunction gives `L_X.T=v-2`; hence

\[
 \boxed{
 \operatorname{ch}_2(F\otimes O_T)=v-2-T^2.}    \tag{5.3}

If `T` is contracted by the generically finite surface map, Hodge index
forces `T^2<0`.  Formula (5.3) is now the only global numerical datum still
missing from the cubic packet: its valency and self-intersection in the
compiled F2 graph.

This is substantially narrower than an arbitrary normalization defect.  To
finish the cubic row one must locate `T` in the compiled boundary tree,
read `(v,T^2)`, insert (5.3), and compare the resulting rank-two curve
charge plus the length-four quotient in (4.3) with the saturated global
budget.  The present theorem classifies the local module but does not make
that global attachment identification.

## 6. Exact saturated budget gate

There is a useful closed form for that final comparison.  Work on the
degree-six, one-component cubic E8 equality row, and assume `T` is the only
new contracted noncyclic packet.  Put

\[
 D'=D_{\log}-2T,qquad
 I=D'\mathbin\cdot T,qquad T^2=-n.              \tag{6.1}
\]

The unsaturated complete-chain residual is `u-1=2`.  Passing from `D_log`
to `D'` changes the half-square by

\[
 \frac12(D_{\log}^2-D'^2)=2I-2n.                \tag{6.2}
\]

Subtracting the rank-two curve charge (5.3), namely `v-2+n`, leaves total
saturated point budget

\[
 \boxed{
 P_{\rm sat}=4-v+2I-3n.}                        \tag{6.3}
\]

The forced local quotient (4.3) already consumes four.  Hence every other
isolated point correction has exact total

\[
 \boxed{
 P_{\rm other}=2I-v-3n\ge0.}                    \tag{6.4}
\]

This is the promised finite attachment gate.  If the saturated determinant
meets `T` only in the cubic affine row, then `I=3`.  Equation (6.4) becomes

\[
 P_{\rm other}=6-v-3n.                          \tag{6.5}

Since `v,n>=1`, the complete survivor list is

\[
\begin{array}{c|c|c}
v&n&P_{\rm other}\\ \hline
1&1&2\\
2&1&1\\
3&1&0.
\end{array}                                     \tag{6.6}

Thus an attachment with `I=3` is impossible if `T^2<=-2` or if `T` has
valency at least four.  The trivalent `(-1)` row is completely saturated;
the leaf and bivalent `(-1)` rows require exactly two and one additional
point units.  If another determinant component meets `T`, its multiplicity
enters only through `I` in (6.4).  No Laurent coefficients are needed to
apply this test once the attachment vertex is identified.

## Reproduction

```bash
.venv/bin/python scripts/verify_log_contracted_divisor_smith_classification.py
```

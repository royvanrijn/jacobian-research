# A homogeneous three-variable counterexample to GVC

## 1. Explicit statement

Work over a characteristic-zero field.  In variables \(x,y,t\), put

\[
 \rho=t^2+xy,
 \qquad
 A=\rho+x^2,
 \qquad
 C=y\rho^2-2xt^2\rho-x^3t^2.
 \tag{1.1}
\]

Define

\[
 \boxed{
 P=A C^2,
 \qquad
 \Delta=4\partial_x\partial_y+\partial_t^2,
 \qquad
 \Lambda=\Delta^6,
 \qquad
 Q=x^2.
 }
 \tag{1.2}
\]

Equivalently,

\[
 P=(t^2+xy+x^2)
 \left(
 yt^4+2xy^2t^2+x^2y^3
 -2xt^4-2x^2yt^2-x^3t^2
 \right)^2.
 \tag{1.3}
\]

The polynomial \(P\) is primitive, homogeneous of degree twelve, and has
23 expanded terms.  The operator \(\Lambda\) is homogeneous of order
twelve.  Its quadratic base symbol \(4XY+T^2\) has rank three and is not a
product of linear forms.

> **Theorem 1.1 — explicit GVC(3) counterexample.**  For every \(m\geq1\),
> \[
>  \boxed{\Lambda^m(P^m)=0,}
>  \tag{1.4}
> \]
> whereas
> \[
>  \boxed{\Lambda^m(QP^m)\ne0.}
>  \tag{1.5}
> \]
> More precisely,
> \[
> \boxed{
> \Delta^{6m+1}(x^2P^m)
> =2^{8m+1}(6m+1)!(2m)!
>   \frac{(12m+3)!!}{(4m+1)!!}
> \ne0.
> }
> \tag{1.6}
> \]

Since \(\Lambda^m=\Delta^{6m}\), equation (1.6) proves that the quadratic
polynomial \(\Lambda^m(x^2P^m)\) is nonzero for every \(m\).  Thus the
constant-coefficient Generalized Vanishing Conjecture fails already in
three variables.  Adjoining unused variables gives the same failure in
every dimension at least three.

This does not produce an ordinary-Laplacian counterexample: the operator is
\(\Delta^6\), not \(\Delta\).

## 2. Polynomiality and the homogeneous lift

The key identity is

\[
 \boxed{xC=\rho^3-t^2(\rho+x^2)^2.}
 \tag{2.1}
\]

Indeed, using \(xy=\rho-t^2\),

\[
\begin{aligned}
 xC
 &=xy\rho^2-2x^2t^2\rho-x^4t^2\\
 &=(\rho-t^2)\rho^2-2x^2t^2\rho-x^4t^2\\
 &=\rho^3-t^2(\rho+x^2)^2.
\end{aligned}
\]

On the affine quadric \(\rho=1\), this gives the Laurent presentation

\[
 C=\frac{1-t^2(1+x^2)^2}{x},
 \tag{2.2}
\]

and consequently

\[
 \boxed{
 P=x^{-2}(1+x^2)
       \left(1-t^2(1+x^2)^2\right)^2.
 }
 \tag{2.3}
\]

The expression in (2.3) is used only for coefficient extraction.  Formula
(1.2) proves directly that \(P\) is polynomial.  The square in (2.3) is the
first endpoint-contact power for which the winding \(x^{-2}\) can be
homogenized inside the same three variables.

## 3. All-order spherical moments

Over \(\mathbb C\), write

\[
 x=X+iY,\qquad y=X-iY,\qquad t=T.
 \tag{3.1}
\]

Then

\[
 \rho=X^2+Y^2+T^2,
 \qquad
 \Delta=\partial_X^2+\partial_Y^2+\partial_T^2.
 \tag{3.2}
\]

On the unit sphere \(\rho=1\), parameterize

\[
 x=\sqrt{1-t^2}\,e^{i\theta},
 \qquad
 y=\sqrt{1-t^2}\,e^{-i\theta}.
 \tag{3.3}
\]

Normalized surface measure is

\[
 \frac{dt}{2}\frac{d\theta}{2\pi},
 \qquad -1\leq t\leq1.
 \tag{3.4}
\]

Put \(u=x^2\).  Phase integration is constant-term extraction in \(x\).
Because the height dependence in (2.3) is even, for every \(m\geq1\),

\[
\begin{aligned}
 \int_{S^2}P^m\,d\sigma
 &=[u^m]K_m(u),\\
 \int_{S^2}x^2P^m\,d\sigma
 &=[u^{m-1}]K_m(u),
\end{aligned}
\tag{3.5}
\]

where

\[
 K_m(u)
 =(1+u)^m\int_0^1
   \left(1-v^2(1+u)^2\right)^{2m}\,dv.
 \tag{3.6}
\]

Let \(B=1+u\) and define

\[
 J_m(B)=\int_0^B(1-w^2)^{2m}\,dw.
 \tag{3.7}
\]

Changing variables \(w=vB\) in (3.6) gives

\[
 K_m(u)=B^{m-1}J_m(B).
 \tag{3.8}
\]

Now

\[
 J_m'(B)=(1-B^2)^{2m}
 \tag{3.9}
\]

has a zero of order \(2m\) at \(B=1\).  Therefore

\[
 J_m(1+u)=J_m(1)+O(u^{2m+1}).
 \tag{3.10}
\]

In every degree at most \(m\), equation (3.8) may consequently be replaced
by

\[
 K_m(u)=J_m(1)(1+u)^{m-1}+O(u^{m+1}).
 \tag{3.11}
\]

The degree-\(m\) coefficient is zero, while the degree-\((m-1)\)
coefficient is \(J_m(1)\).  Hence

\[
 \boxed{
 \int_{S^2}P^m\,d\sigma=0,
 \qquad
 \int_{S^2}x^2P^m\,d\sigma=C_m\ne0,
 }
 \tag{3.12}
\]

with

\[
\begin{aligned}
 C_m
 &=\int_0^1(1-w^2)^{2m}\,dw\\
 &=\frac{2^{2m}(2m)!}{(4m+1)!!}.
\end{aligned}
\tag{3.13}
\]

This is the same adjacent Taylor-coefficient mechanism as the Long/Hopf
counterexample, but with contact order two and an internal homogeneous
power lift.

## 4. From the Reynolds functional to differential identities

The general
[Reynolds--apolar transfer lemma](CUSP_PROFILE_SUSPENSION_THEOREM.md#2-reynolds--apolar-transfer)
says that, for every homogeneous polynomial \(F\) of degree \(2k\),

\[
 \boxed{
 \Delta^kF
 =2^kk!(2k+1)!!
   \mathcal R_\rho(F|_{\rho=1}).}
 \tag{4.1}
\]

This is an algebraic orthogonal-invariance identity: both sides are
invariant functionals on \(\operatorname{Sym}^{2k}\), and evaluation on
\(\rho^k\) fixes the constant.  On the real compact form,
\(\mathcal R_\rho\) is the normalized sphere average used in Section 3.
Thus the proof does not require expanded differentiation or Gaussian
integration.  The Gaussian/Wick calculation in the checker remains an
independent finite replay.

Since \(P^m\) has degree \(12m\), equation (3.12) and (4.1) give

\[
 \Delta^{6m}(P^m)=0.
 \tag{4.2}
\]

This is exactly (1.4).  Likewise, \(x^2P^m\) has degree \(12m+2\), so

\[
\begin{aligned}
 \Delta^{6m+1}(x^2P^m)
 &=2^{6m+1}(6m+1)!(12m+3)!!\,C_m\\
 &=2^{8m+1}(6m+1)!(2m)!
   \frac{(12m+3)!!}{(4m+1)!!}.
\end{aligned}
\tag{4.3}
\]

This proves (1.6).  The Reynolds functional and all displayed differential
identities are defined over the ground field, so the conclusion holds over
every characteristic-zero field.

## 5. Infinite homogeneous family

The construction is not isolated.  For every integer \(s\geq2\), put

\[
 P_s=(\rho+x^2)x^{s-2}C^s,
 \qquad
 \Lambda_s=\Delta^{3s},
 \qquad
 Q=x^2.
 \tag{5.1}
\]

Then \(P_s\) is homogeneous of degree \(6s\), and on \(\rho=1\),

\[
 P_s=x^{-2}(1+x^2)
       \left(1-t^2(1+x^2)^2\right)^s.
 \tag{5.2}
\]

The proof of Section 3 now uses

\[
 J_{s,m}'(B)=(1-B^2)^{sm},
 \tag{5.3}
\]

whose zero at \(B=1\) has order \(sm\geq m\).  Therefore

\[
 \Lambda_s^m(P_s^m)=0,
 \qquad
 \Lambda_s^m(x^2P_s^m)\ne0
 \quad(m\geq1),
 \tag{5.4}
\]

and the exact detector is

\[
 \boxed{
 \Delta^{3sm+1}(x^2P_s^m)
 =2^{4sm+1}(3sm+1)!(sm)!
   \frac{(6sm+3)!!}{(2sm+1)!!}.
 }
 \tag{5.5}
\]

The member \(s=2\) in Theorem 1.1 has the minimum degree in this internal
homogenization family.  No global minimum-degree claim is made.

## 6. Reproduction

Run

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
```

The dependency-free checker uses two independent exact replays through
\(m=6\): sparse differentiation by
\((4\partial_x\partial_y+\partial_t^2)^{6m}\), and direct complex-Gaussian
monomial moments.  It also verifies (2.1), homogeneity, primitivity, the
23-term expansion, the nonzero quadratic mixed output, and formula (1.6).
The all-order result is the proof in Sections 2--4, not the bounded replay.

The checked-in Lean package supplies an independent all-order formal proof
of Theorem 1.1.  It verifies the literal ternary definitions and
homogeneity, coefficientwise differential semantics, the algebraic
Reynolds/Laurent phase extraction, the endpoint coefficient calculation and
exact scalar, characteristic-zero coefficient base change, and
unused-variable padding.  Build it with:

```bash
make verify-gvc-lean
```

There are no `sorry`, `admit`, or explicit `axiom` declarations in that
package.  Primitivity and the 23-term expanded form are covered by the exact
checker above rather than by named Lean theorems.

## 7. Provenance and scope

The endpoint-contact coefficient mechanism comes from Christopher D.
Long's three-real-variable Gaussian counterexample and the repository's
Hopf-lift classification.  The new step here is to square the endpoint
factor and use the quadric identity (2.1) to obtain a homogeneous polynomial
without adding variables.  This note makes no priority claim before
independent literature and expert review.

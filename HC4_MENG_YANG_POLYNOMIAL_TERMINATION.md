# Meng--Yang polynomial-termination frontends

## Status

No `HC(4)` counterexample or general nontermination theorem is constructed
here.  This note begins the polynomial-termination analysis left by the
all-normal theorem `HC4MYGJ2` with six exact, degree-free statements and one
bounded calculation.

> **Theorem `HC4MYPT1` -- terminal coefficient cone gate.**  Let
> \(R\in K[x,y,p,q]\) have positive \(x\)-degree \(n\), write
> \(R=x^nU+O(x^{n-1})\), and suppose the Meng--Yang graph determinant is
> constant.  Put \(d=2n+6\).  Then
> \[
> (d-1)U\det H_U-(d+1)\nabla U^{\mathsf T}
> \operatorname{adj}(H_U)\nabla U=0.                 \tag{T.1}
> \]
> If \(U_m\) is the highest nonzero homogeneous part of \(U\) in
> \((y,p,q)\), then \(\det H_{U_m}=0\).  Hence, over a
> characteristic-zero field, the ternary Hesse theorem makes \(U_m\) a
> cone after scalar extension.  This is a necessary last-coefficient
> equation, not a proof that a polynomial branch cannot terminate.

> **Theorem `HC4MYPT2` -- pure weight-five one-channel obstruction.**  No
> polynomial graph
> \[
> R=\frac NL y^5f(xy),\qquad f\in K[z],              \tag{T.2}
> \]
> has constant Meng--Yang graph determinant when \(LN\ne0\).  This is an
> all-degree exclusion of the displayed one-function ansatz.  It does not
> exclude the coupled weight packet forced above a general plane trace.

> **Theorem `HC4MYPT3` -- coupled five-channel terminal split.**  Put
> \(z=xy\) and consider the exact five-function graph
>
> \[
> R=y^5f(z)+y^3p\,g(z)+y^3h(z)+y^2q\,j(z)
>     +yp^2\ell(z).                                \tag{T.3}
> \]
> Let \(n>0\) be the largest degree of the five polynomials and write
> \(A,B,C,D,E\) for their respective \(z^n\)-coefficients.  If the
> Meng--Yang graph determinant is constant, then
>
> \[
> D=0\qquad\hbox{or}\qquad B=E=0.                  \tag{T.4}
> \]
> Equivalently, if \(j\) has maximal degree, then both \(g\) and \(\ell\)
> have smaller degree.  This classifies the last-coefficient face of the
> coupled packet; it does not exclude either surviving branch or graph
> channels below (T.3).

> **Theorem `HC4MYPT4` -- coupled upper-Newton chamber gate.**  Assume
> \(f,g,j,\ell\ne0\) in (T.3), with respective degrees \(a,b,d,e\).  The
> complete weight-six five-function equation is independent of \(h\), and
> its upper Newton degree is
>
> \[
> \max\{5a+e+2d+22,\ 4a+2b+2d+22\}.              \tag{T.5}
> \]
> If \(d>0\), the strict chamber \(2b>a+e\) is empty in all polynomial
> degrees.  In the other strict chamber \(a+e>2b\), a balance must satisfy
> the explicit genus-one resonance (5.4) below; the wall has the amplitude
> equation (5.7).  If \(d=0\), either strict chamber first forces the
> constant \(j=3\), and the next Newton face excludes that exceptional
> value.  Thus every degree-zero \(j\)-balance lies on the wall.  This is an
> exact all-degree Newton gate, not an exclusion of the surviving resonance
> or wall strata.

> **Theorem `HC4MYPT5` -- arithmetic closure of the strict resonance.**
> The equation \(P_F(a,d)=0\) from (5.4) has no solution with integers
> \(a>0,d\ge0\).  Its only solution with \(a=0\) is \(d=2\).  Continuing
> that exceptional balance by one Newton face excludes every strict
> \(f^5\ell\)-balance except the ridge
>
> \[
> a=0,\qquad d=2,\qquad e=2b+2,\qquad b\ge1,       \tag{T.6}
> \]
>
> whose leading coefficients must satisfy
>
> \[
> 360AE+B^2D(b+5)^2=0.                            \tag{T.7}
> \]
>
> Thus the apparent positive-degree genus-one family is empty.  The
> exceptional ridge (T.6)--(T.7) is not excluded.

> **Theorem `HC4MYPT6` -- algebraic fifth-channel elimination.**  The
> complete coupled equation contains neither \(\ell'\) nor \(\ell''\) and
> is affine in \(\ell\):
>
> \[
> \mathcal E_5=\mathscr C\ell+\mathscr R,\qquad
> \mathscr C=2z(2Lz^5f-3Nz+2N)\mathscr Q.         \tag{T.8}
> \]
>
> Here \(\mathscr Q\ne0\) is independent of \(g,h,\ell\).  On every
> admissible wall balance (5.7), the coefficient that solves successively
> for the lower \(\ell\)-coefficients is
> \(-32L^4A^5D^2P_F(a,d)\ne0\).  Hence the wall expansion at infinity is
> unit-triangular in the fifth channel.  Polynomiality is reduced to the
> exact divisibility condition \(\mathscr C\mid\mathscr R\), together with
> its finite remainder equations; that divisibility is not proved here.

For the collision-containing plane-flat near miss of the
[graph-obstruction note](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md), exact
recursion gives nonzero coefficients through normal order five, and the
terminal bracket is nonzero at each of orders one through five.  Orders
three through five are new beyond the first-transverse obstruction.  This
finite prefix is an experiment/certificate, not an all-order nontermination
proof.

## 1. The semilinear normal equation

Use the notation

\[
 u=1+xy,\qquad
 A_0=u^3p+3xu^2q,
\]

and let \(\Phi_R\) be the five-variable Meng--Yang potential pulled back to
the graph \(r=R(x,y,p,q)\).  Direct expansion gives the exact differential
polynomial presentation

\[
 \Phi_R=F_0+cR+Lx^6R^2,                              \tag{1.1}
\]

where \(F_0\) is independent of \(R\) and

\[
 c=-2Lx^3A_0-Mx^3+N(2x-3x^2y).                     \tag{1.2}
\]

Only the \((x,x)\)-entry of \(\operatorname{Hess}\Phi_R\) contains
\(R_{xx}\).  Consequently

\[
 D_R=\mathcal A_R R_{xx}+\mathcal B_R,              \tag{1.3}
\]

with

\[
 \mathcal A_R=(c+2Lx^6R)
 \det\bigl((\Phi_R)_{ij}\bigr)_{i,j\in\{y,p,q\}},  \tag{1.4}
\]

and \(\mathcal B_R\) independent of \(R_{xx}\).  Thus the graph equation is
semilinear in the highest normal derivative.  Its leading coefficient
vanishes at \(x=0\); this is a regular-singular rather than an ordinary
Cauchy equation.  The Taylor/cofactor calculation in `HC4MYGJ2` packages its
indicial operator as

\[
 -4LN^3\partial_x^2(xR),                             \tag{1.5}
\]

and gives the exact recursion

\[
 R_k=
 \frac{[x^{k-1}](D_{R^{(<k)}}-C)}
 {4LN^3k(k+1)},
 \qquad
 R^{(<k)}=\sum_{i<k}x^iR_i.                         \tag{1.6}
\]

Equation (1.6) is the termination engine: no coefficient of \(R_k\) is a
search variable.

## 2. The backward terminal equation

Assume that a polynomial solution terminates at \(x\)-degree \(n>0\):

\[
 R=x^nU(y,p,q)+O(x^{n-1}),\qquad U\ne0.             \tag{2.1}
\]

The highest \(x\)-degree part of (1.1) is

\[
 Lx^dU^2,\qquad d=2n+6.                             \tag{2.2}
\]

For an arbitrary ternary polynomial \(W\),

\[
 \det\operatorname{Hess}_{x,y,p,q}(x^dW)
 =x^{4d-2}
 \det\begin{pmatrix}
 d(d-1)W&d\nabla W^{\mathsf T}\\
 d\nabla W&H_W
 \end{pmatrix}.                                    \tag{2.3}
\]

Set \(W=U^2\).  Since

\[
 \nabla(U^2)=2U\nabla U,
 \qquad H_{U^2}=2(UH_U+\nabla U\nabla U^{\mathsf T}),
\]

the rank-one determinant identity factors (2.3) as

\[
 8dU^4\left((d-1)U\det H_U
 -(d+1)\nabla U^{\mathsf T}
 \operatorname{adj}(H_U)\nabla U\right).          \tag{2.4}
\]

No other term of \(\Phi_R\) reaches the same \(x\)-degree.  Constancy of
\(D_R\) therefore proves (T.1).

Let \(U_m\) be homogeneous of degree \(m\ge2\).  Euler's identities give

\[
 \nabla U_m^{\mathsf T}\operatorname{adj}(H_{U_m})
 \nabla U_m
 =\frac{m}{m-1}U_m\det H_{U_m}.                    \tag{2.5}
\]

The bracket in (2.4) becomes

\[
 -\frac{d+2m-1}{m-1}U_m\det H_{U_m}.               \tag{2.6}
\]

Here \(d+2m-1=2n+2m+5\ne0\) in characteristic zero.  The cases
\(m=0,1\) already have zero Hessian.  This proves `HC4MYPT1`.

The gate is sharp enough to remove every terminal top monomial involving
all three tangential variables: the Hessian of
\(y^ap^bq^c\) is nonsingular when \(a,b,c>0\).  Cone terminal coefficients,
including pure powers, survive and require the next backward equation.

## 3. A polynomial balance at \(x=\infty\)

The first closed dominant-balance model is the pure weight-five graph

\[
 R=\frac NL y^5f(z),\qquad z=xy.                   \tag{3.1}
\]

Give \(x\) weight \(-1\) and \(y,p,q\) weight \(1\).  After differentiating
the full potential, set \(p=q=0\), replace \(x\) by \(z/y\), and take the
weight-six determinant coefficient.  Scaling \(R\) as in (3.1) removes
\(L,N\); \(M\) has lower weight.  The result is a semilinear equation

\[
 \mathcal E(z,f,f',f'')=0.                          \tag{3.2}
\]

The complete 143-term polynomial is generated by the checker rather than
printed here.  Its \(f''\)-coefficient factors compactly:

\[
\begin{aligned}
[f'']\mathcal E={}&-2z(z+1)^6(2z^5f-3z+2)\\
&\cdot(6z^6f+12z^5f+6z^4f-9z^2-12z+1)^2.          \tag{3.3}
\end{aligned}
\]

For a prospective polynomial balance \(f\sim az^m\), \(a\ne0\), the five
sectors of \(\mathcal E\), ordered by their degree in \(f,f',f''\), have
dominant data

| sector | \(z\)-degree | leading coefficient |
|---:|---:|---|
| 4 | \(4m+22\) | \(-144a^4(m^2-m-21)\) |
| 3 | \(3m+18\) | \(216a^3(3m^2-3m-64)\) |
| 2 | \(2m+14\) | \(-972a^2(m^2-m-23)\) |
| 1 | \(m+10\) | \(486a(m-6)(m+5)\) |
| 0 | \(6\) | \(2916\) |

For every integer \(m\ge0\), the quartic sector strictly dominates.  Its
resonance equation has discriminant \(85\), so it has no integral root.
The zero polynomial also fails because the sector-zero coefficient is
nonzero.  This proves `HC4MYPT2`.

This calculation is deliberately scoped to (3.1).  On the actual
collision-containing branch, determinant rank loss couples the weight-five
face to at least the four channels

\[
 y^3p\,g(z),\quad y^3h(z),\quad y^2q\,j(z),
 \quad yp^2\ell(z).                                \tag{3.4}
\]

Discarding those channels gives a false scalar recurrence.  The first exact
attack on their coupled system is the terminal split below, not an
extrapolation of (3.2).

## 4. The coupled five-channel terminal split

For (T.3), the coefficient of \(x^n\) is

\[
 U=A y^{n+5}+B y^{n+3}p+C y^{n+3}
      +D y^{n+2}q+E y^{n+1}p^2.                  \tag{4.1}
\]

Apply the universal terminal bracket (T.1), then restrict the tangential
variables to \(p=Py,\ q=Qy\).  Exact expansion gives

\[
\begin{aligned}
\mathcal T(U)|_{p=Py,q=Qy}=D^2y^{4n+6}\bigl[{}
 &\bigl(8AE n^2+14AE n-40AE+2B^2n+7B^2\bigr)y^2\\
 &+2BEP(4n^2+15n+8)y\\
 &+2E(n+2)(4n+11)(C+DQ+EP^2)\bigr].              \tag{4.2}
\end{aligned}
\]

This substitution is injective on \(K[y,p,q]\), since the exponents of
\(P,Q\), and then of \(y\), recover the original monomial.

If \(D=0\), this expression vanishes.  Suppose instead that \(D\ne0\).
The coefficient of \(Q\) is

\[
 2D^3E(n+2)(4n+11)y^{4n+6},                       \tag{4.3}
\]

so characteristic zero and \(n>0\) force \(E=0\).  Equation (4.2) then
reduces to

\[
 D^2B^2(2n+7)y^{4n+8},                            \tag{4.4}
\]

and forces \(B=0\).  Conversely \(D=0\), or \(B=E=0\), makes the complete
terminal bracket vanish identically.  This proves `HC4MYPT3`.

The split is useful but not yet a coupled nontermination theorem.  It leaves
two unequal-degree Newton branches: either \(j\) drops below the maximal
degree, or \(j\) is maximal while both \(g\) and \(\ell\) drop.  The channel
\(h\) is invisible at this terminal face, as is its leading coefficient
\(C\).  Lower graph-weight channels outside (T.3) can also enter subsequent
faces, so (T.4) is asserted only for the exact five-function graph.

## 5. The complete coupled upper Newton face

Differentiate the full potential before making the substitutions

\[
 x=z/y,\qquad p=Py,\qquad q=Qy.                   \tag{5.1}
\]

The coefficient of \(y^6\) in the Hessian determinant is a 2,348-term
differential polynomial in \(f,g,j,\ell\) and their first two derivatives.
The parameters \(M,P,Q\) and the complete function \(h\) cancel.  Denote this
equation by \(\mathcal E_5=0\).  Grouping its monomials by derivative count
leaves 396 Newton sectors.

Let

\[
 a=\deg f,\quad b=\deg g,\quad d=\deg j,\quad
 e=\deg\ell,                                      \tag{5.2}
\]

and let \(A,B,D,E\) be the corresponding nonzero leading coefficients.
Exact comparison of all 396 affine degree functions gives the upper envelope

\[
 \nu=\max\{5a+e+2d+22,\ 4a+2b+2d+22\}.           \tag{5.3}
\]

On the strict chamber \(a+e>2b\), and with \(d>0\), the unique coefficient
at degree \(\nu\) is

\[
 -32L^4A^5D^2E\,P_F(a,d),
\]

where

\[
\begin{aligned}
P_F(a,d)={}&2a^3-4a^2d+15a^2-2ad^2-34ad\\
            &+29a-30d+60.                         \tag{5.4}
\end{aligned}
\]

Thus \(P_F(a,d)=0\) is necessary.  For \(a=0\), it is exactly
\(-30(d-2)=0\).  For \(a>0\), its discriminant as a quadratic in \(d\) is

\[
 4(a+3)(a+5)(2a+1)(4a+15),                        \tag{5.5}
\]

so an integral balance must in particular make the four-factor product in
(5.5) a square.  Section 5.1 closes this integral-point condition
elementarily.

On the other strict chamber \(2b>a+e\), the upper coefficient is

\[
 -16L^4A^4B^2D^2\,P_G(a,b,d),
\]

where

\[
\begin{aligned}
P_G(a,b,d)={}&a^2-2ab^2+4abd-8ab-2ad^2+6ad-a\\
 &-b^2+4bd-18b-4d^2+6d-21.                       \tag{5.6}
\end{aligned}
\]

Put \(t=2b-a>0\).  The discriminant of (5.6), as a quadratic in \(d\), is

\[
-2\bigl(a^3+6a^2t+22a^2+at^2+44at+116a+60t+150\bigr)<0.
\]

Its leading coefficient is \(-2a-4<0\), so \(P_G<0\) for every real \(d\).
This excludes the complete strict \(g^2\)-dominant chamber in every
polynomial degree.

On the wall \(a+e=2b\), both sectors contribute and the upper equation is

\[
 2AE\,P_F(a,d)+B^2P_G(a,b,d)=0.                  \tag{5.7}
\]

Finally suppose \(d=0\), so \(D=j\) is the nonzero constant value.  The
three tied powers of \(j\) combine into the square \((D-3)^2\).  In the two
strict chambers the remaining factors are respectively

\[
 2a^3+15a^2+29a+60>0
\]

and \(P_G(a,b,0)<0\).  Hence either strict chamber forces \(j=3\).  On the
wall the residual amplitude equation remains and may instead cancel.

The exceptional value \(j=3\) can be continued exactly.  Substitution into
\(\mathcal E_5\) leaves 512 terms and 121 Newton sectors, with upper envelope

\[
 \max\{5a+e+20,\ 4a+2b+20\}.                     \tag{5.8}
\]

The strict \(f^5\ell\)-chamber coefficient contains

\[
 2a^3+19a^2+61a+90>0,                             \tag{5.9}
\]

while the strict \(g^2\)-chamber contains

\[
 Q_3(a,b)=a^2-2ab^2-12ab-9a-b^2-22b-31.          \tag{5.10}
\]

For \(t=2b-a>0\),

\[
4Q_3\left(a,\frac{a+t}{2}\right)
=-\bigl(2a^3+4a^2t+21a^2+2at^2+26at+80a
        +t^2+44t+124\bigr)<0.                    \tag{5.11}
\]

Thus neither strict chamber survives when \(d=0\).  This proves
`HC4MYPT4`.  The remaining strict resonance and wall are sharpened next.

### 5.1. Arithmetic closure of the strict resonance

Suppose \(a>0\) and \(P_F(a,d)=0\).  The values at \(d=0,1\) are positive,
so \(d\ge2\).  Reduction modulo \(a\) gives

\[
 P_F(a,d)\equiv-30(d-2)\pmod a.                  \tag{5.12}
\]

Consequently \(k=30(d-2)/a\) is a nonnegative integer.  Substitution into
(5.4) reduces the resonance to

\[
\begin{aligned}
Q_k(a)={}&a^2(k^2+60k-900)+a(630k-3150)\\
         &+450k+21150=0.                          \tag{5.13}
\end{aligned}
\]

For \(k\ge13\), all three displayed coefficients are positive; for the
quadratic coefficient use

\[
k^2+60k-900=49+(k-13)(k+73).
\]

It remains to inspect \(0\le k\le12\).  The discriminant of (5.13) in
\(a\) is

\[
 900(51-k)(25-k)(75-2k).                          \tag{5.14}
\]

After removing the square factor \(900\), its thirteen values are

\[
95625,87600,80017,72864,66129,59800,53865,48312,
43129,38304,33825,29680,25857.
\]

None is a square.  Thus no \(a>0\) integral resonance exists.  When
\(a=0\), equation (5.4) is \(-30(d-2)=0\), so \(d=2\).

At \(a=0,d=2\), the vanished top \(f^5j^2\ell\)-sector is followed by the
\(f^5j\ell\)-sector, whose coefficient is
\(5760L^4A^5DE\).  The competing \(g^2j^2\)-sector contains
\(P_G(0,b,2)=-(b+5)^2\).  Write \(e=2b+\delta\), where strict
\(f^5\ell\)-dominance gives \(\delta\ge1\).  The \(g^2j^2\)-sector dominates
when \(\delta=1\), the \(f^5j\ell\)-sector dominates when \(\delta\ge3\),
and both coefficients are nonzero.  At the only tie \(\delta=2\), their
sum is

\[
16L^4A^4D\bigl(360AE+B^2D(b+5)^2\bigr).          \tag{5.15}
\]

The terminal split `HC4MYPT3` excludes \(b=0\) on this tie because then
\(j\) and \(\ell\) both have maximal degree while \(D,E\ne0\).  This proves
`HC4MYPT5`: the strict chamber is reduced to (T.6)--(T.7), with \(b\ge1\).

## 6. Algebraic elimination of the fifth channel

The complete equation has a structural simplification stronger than its
Newton polygon.  Exact collection in \(\ell,\ell',\ell''\) gives

\[
 \mathcal E_5=\mathscr C\ell+\mathscr R,
 \qquad
 \mathscr C=2z(2Lz^5f-3Nz+2N)\mathscr Q,          \tag{6.1}
\]

with no derivatives of \(\ell\).  The expanded \(\mathscr Q\) has 311
terms, is nonzero, and depends only on \(z,L,N,f,j\) and their derivatives;
in particular it is independent of \(g,h,\ell\).  The remainder
\(\mathscr R=\mathcal E_5|_{\ell=0}\) has 1,895 terms.

The explicit factor \(z\) immediately gives a necessary finite-end
condition.  Writing \(f_i=f^{(i)}(0)\), \(g_0=g(0)\), and \(j_0=j(0)\), it is

\[
\begin{aligned}
0={}&80Lf_0j_0+492Lf_0-8Lf_1+4Ng_0^2-64Ng_0j_0\\
   &-356Ng_0+256Nj_0^2+2848Nj_0+7921N.            \tag{6.2}
\end{aligned}
\]

There is also a useful consequence at infinity.  On the wall
\(a+e=2b\), the next equation is linear in the coefficient of \(z^{e-1}\)
in \(\ell\), with multiplier

\[
 -32L^4A^5D^2P_F(a,d).                            \tag{6.3}
\]

This multiplier cannot vanish on an admissible wall.  Indeed,
`HC4MYPT5` removes every \(a>0\) zero of \(P_F\).  The remaining zero
\(a=0,d=2\) has \(P_G(0,b,2)=-(b+5)^2\ne0\), so the wall equation (5.7)
cannot hold.  Polynomial long division therefore determines every lower
\(\ell\)-coefficient uniquely from the other channels.  Equivalently, an
exact polynomial solution must satisfy

\[
 2z(2Lz^5f-3Nz+2N)\mathscr Q\mid\mathscr R.       \tag{6.4}
\]

This proves `HC4MYPT6`.  It turns the balanced wall from a growing
five-function coefficient ansatz into a four-function divisibility problem,
but it does not prove the divisibility impossible.  The same condition
applies to the exceptional ridge (T.6).

## 7. Exact bounded recurrence on the collision near miss

For the v2 parameters \((L,M,N)=(1,13,2)\), take the trace and target

\[
 T_*=\frac{51}{100}y^5-\frac{47}{10}y^3p
       -\frac{123}{20}y^2q,
 \qquad C=\frac{17165601}{25}.                      \tag{7.1}
\]

This trace, together with its forced first jet, contains both marked points
and makes the determinant equal to \(C\) on \(x=0\).  Applying (1.6) with
exact truncated polynomial arithmetic gives:

| \(k\) | terms in \(R_k\) | \(\deg R_k\) | extremal term |
|---:|---:|---:|---|
| 1 | 11 | 4 | \((1989/40)y^4\) |
| 2 | 17 | 7 | \((459/500)y^7\) |
| 3 | 24 | 8 | \((235467/4000)y^8\) |
| 4 | 29 | 9 | \(-(51583797/20000)y^9\) |
| 5 | 38 | 10 | \((1094341514607/10000000)y^{10}\) |

The \(k=2\) row is equivalent to the recorded nonzero first-transverse
coefficient.  The last three rows extend the exact prefix.  They show that the
unique formal branch above this plane-flat near miss has not terminated
before order five.  More strongly, substituting each complete \(R_k\), not
just its extremal monomial, into (T.1) gives a nonzero terminal bracket for
every \(1\le k\le5\).  Thus none of these five coefficients can be the last
one of a polynomial solution.

The emerging \(j\)- and \(\ell\)-channel coefficients

\[
 D_k=[y^{k+2}q]R_k,\qquad E_k=[y^{k+1}p^2]R_k
\]

are already both nonzero at every computed order:

| \(k\) | \(D_k\) | \(E_k\) |
|---:|---:|---:|
| 1 | \(-2727/200\) | \(348/5\) |
| 2 | \(243/100\) | \(-158457/50\) |
| 3 | \(-10145133/10000\) | \(53188533/400\) |
| 4 | \(2083052601/50000\) | \(-22284332193/4000\) |
| 5 | \(-880031929203/500000\) | \(46567165753557/200000\) |

This is evidence that the coupled face is persistent, not a proof that it
cannot eventually enter one of the branches in (T.4).  The higher corrections
are not asserted to retain the marked points.  The calculation therefore does
not prove
that a later coefficient cannot vanish or that no later terminal packet can
occur.  Moreover, each displayed extremal is a pure power and therefore
passes `HC4MYPT1`; a successful backward proof must use a lower face of the
terminal coefficient or the next terminal equation.

## 8. Reproduction and next target

Run

```bash
.venv/bin/python scripts/verify_hc4_meng_yang_polynomial_termination.py
```

The command checks the exact semilinear presentation, the abstract
bordered-Hessian factor (2.4), the complete determinant-defined equation
(3.2) and its dominant table, the coupled factorization (4.2), the complete
2,348-term coupled weight-six equation and all 396 of its Newton sectors,
the arithmetic reduction (5.12)--(5.15), the fifth-channel factorization
(6.1), its axis equation (6.2), and the rational recurrence (7.1) through
normal order five.  It requires only the locked Python/SymPy environment.

The positive-degree genus-one target and the free fifth-channel wall
recurrence are now closed.  The two coupled targets are sharper: disprove
the divisibility (6.4) on the exceptional ridge (T.6)--(T.7), and compute
the terminal remainder of the unit-triangular wall division.  In parallel,
the backward terminal calculation should extract the first cofactor equation
below (2.4).  The generic corank-one branch should force second-order
invariance of the penultimate coefficient; the pure-power branch visible in
Section 7 is the lower-rank exception that must be handled separately.

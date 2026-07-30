# The root-incidence derivative split for \(A_4\)

## 1. Scope

The cone construction has already isolated the determinant
\(4W^2K^3L\), and the subsequent ledger calculations show that a Keller
completion inside that model must alter at least two outputs by
source-dependent masks.  This note changes the starting point.  It returns
to the two-parameter generic \(A_4\) quartic and asks what the derivative
unit itself says in the rank-four root algebra.

The outcome is an exact partial construction:

- \(1/P'(T)\) has a compact four-term representative;
- the square discriminant can be split evenly between two incidence
  coordinates;
- the first new coordinate remains a primitive element, so the generic
  four-sheet field is not lost;
- one irreducible orientation pole survives polynomial pullback over the
  target coefficient ring.

Thus the quotient-algebra part of the proposed construction works.  It is
not yet an affine polynomial chart and does not give a Keller map.

## 2. Compact rank-four inverse

Use temporary coefficient variables \(A,B,C\) and put

\[
\begin{aligned}
P(T)={}&T^4-6ABT^2-8B^3T+B^2(9A^2-12CB),\\
R={}&4A^3B-3A^2C^2-6AB^2C+B^4+4BC^3.
\end{aligned}
\]

Then

\[
\operatorname{Disc}_T(P)=-110592B^8R.                \tag{2.1}
\]

In the free rank-four algebra

\[
\mathcal A=\mathbb Q(A,B,C)[T]/(P),
\]

define

\[
\begin{aligned}
n(T)={}&(B^2-AC)T^3+2B(A^2-BC)T^2\\
 &+B(3A^2C-7AB^2+4BC^2)T\\
 &-6B^2(A^3-2ABC+B^3).
\end{aligned}
\]

Exact reduction modulo \(P\) gives

\[
\boxed{\frac1{P'(T)}=\frac{n(T)}{48B^4R}.}            \tag{2.2}
\]

Equivalently, with

\[
N(T)=-2304B^4n(T),
\]

one has the unreduced discriminant identity

\[
\boxed{P'(T)N(T)\equiv\operatorname{Disc}(P)\pmod P.} \tag{2.3}
\]

The unreduced form is the useful one: it remembers that the denominator is
the square of one orientation equation.

## 3. Two-coordinate derivative-unit split

The polynomial \(N(T)\) has the explicit primitive

\[
\begin{aligned}
I(T)={}&576B^4(AC-B^2)T^4
       -1536B^5(A^2-BC)T^3\\
 &-1152B^5(3A^2C-7AB^2+4BC^2)T^2\\
 &+13824B^6(A^3-2ABC+B^3)T.
\end{aligned}                                       \tag{3.1}
\]

Thus \(I'(T)=N(T)\).  On the oriented chart

\[
\Omega^2=\operatorname{Disc}(P),
\]

introduce one incidence variable \(Q\) and set

\[
X=\frac{I(T)}{\Omega},\qquad Y=\frac Q\Omega.         \tag{3.2}
\]

Since \(\Omega\) is independent of \(T,Q\),

\[
\boxed{
\det\frac{\partial(X,Y)}{\partial(T,Q)}
=\frac{N(T)}{\Omega^2}
=\frac1{P'(T)}
}
\quad\text{in }\mathcal A.                           \tag{3.3}
\]

This is the desired algebraic cancellation: the two coordinates each
carry one copy of the orientation pole, and their combined Jacobian is the
derivative unit.

This split does not collapse the inverse cover.  Reduce \(I\) modulo \(P\)
and form the change-of-basis matrix from

\[
(1,T,T^2,T^3)
\quad\text{to}\quad
(1,I,I^2,I^3).
\]

Its determinant is a nonzero polynomial.  For example, at
\((A,B,C)=(0,1,0)\) one has

\[
I\equiv9216T\pmod P,
\]

so the determinant specializes to \(9216^6\ne0\).  Hence

\[
\mathbb Q(A,B,C)(I)=\mathbb Q(A,B,C)(T)              \tag{3.4}
\]

generically.  The root \(T\), and therefore the original quartic relation,
is rationally recoverable from \(I\).  This proves field recovery, not yet
a literal polynomial inverse equation in the new coordinates.

## 4. The surviving polynomiality obstruction

Now specialize to the Jensen--Ledet--Yui parameters

\[
\begin{aligned}
A_0={}&\alpha^3-\beta^3-9\beta^2-27\beta-54,\\
B_0={}&\alpha^3-3\alpha\beta^2+2\beta^3-9\alpha\beta
       +9\beta^2-27\alpha+27\beta+27,\\
C_0={}&\alpha^3-\beta^3+27.
\end{aligned}
\]

Put

\[
\rho=\beta^2+3\beta+9
\]

and

\[
\begin{aligned}
\sigma={}&2\alpha^3\beta+3\alpha^3
-3\alpha^2\beta^2-9\alpha^2\beta-27\alpha^2\\
&+\beta^4+6\beta^3+27\beta^2+54\beta+81.
\end{aligned}
\]

The orientation equation is the polynomial

\[
\boxed{\Omega=1728\rho B_0^4\sigma,}
\qquad
\operatorname{Disc}(P_{\alpha,\beta})=\Omega^2.      \tag{4.1}
\]

There is substantial cancellation in (2.2).  In fact

\[
R(A_0,B_0,C_0)=-27\rho^2\sigma^2,
\]

and all four coefficients of \(n(T)\) acquire a common factor
\(-3\rho\).  The reduced denominator of \(1/P'\) is consequently

\[
432\rho B_0^4\sigma^2.                               \tag{4.2}
\]

This cancellation is not enough to make (3.2) polynomial.  Modulo \(P\),
the \(T^3\)-coefficient of \(I\) is

\[
-1536B^5(A^2-BC).
\]

After specialization,

\[
A_0^2-B_0C_0=3\rho E,
\]

where

\[
\begin{aligned}
E={}&\alpha^4-\alpha^3\beta-6\alpha^3-\alpha\beta^3
+27\alpha\\
&+\beta^4+6\beta^3+27\beta^2+54\beta+81.
\end{aligned}
\]

Therefore the corresponding coefficient of \(I/\Omega\) is

\[
-\frac83\,\frac{B_0E}{\sigma}.                       \tag{4.3}
\]

The exact gcd calculation

\[
\gcd(\sigma,B_0E)=1
\]

shows that (4.3) has a genuine \(\sigma\)-pole.  Adding a
\(\mathbb Q[\alpha,\beta]\)-constant to the primitive changes only its
constant coefficient, so it cannot remove this obstruction.  The second
coordinate \(Q/\Omega\) retains the full orientation pole unless \(Q\) is
given its own source-dependent mask.

Hence the canonical triangular split proves precisely:

\[
\boxed{
\begin{gathered}
\text{the derivative unit admits a two-coordinate orientation split,}\\
\text{but target-only polynomial pullback fails already on }\sigma=0.
\end{gathered}}
                                                               \tag{4.4}
\]

This does not exclude nontriangular pairs linear in \(Q\), nor does it
exclude source charts in which two independent masks cancel the two poles.

## 5. Comparison with the ordinary root-incidence chart

The derivative-unit split must still be compared with the Jacobian of an
ordinary affine chart on the root incidence.  For the rational
\((U,V)\)-chart of the existing \(A_4\) construction, recall

\[
\alpha=\frac{N_1}{H},\qquad
\beta=\frac{N_2}{H},\qquad
\det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
=\frac{4K^3L}{H^3}.                                  \tag{5.1}
\]

Put

\[
\mathsf A=U^3-V^3-9V^2-27V-54.
\]

There are polynomials \(A_6,C_6\) defined by

\[
A_0(\alpha,\beta)=\frac{K^3A_6}{H^3},\qquad
C_0(\alpha,\beta)=\frac{K^3C_6}{H^3},
\]

while

\[
B_0(\alpha,\beta)=\frac{K^3L^2}{H^3}.
\]

The selected scaled root of \(P_{\alpha,\beta}\) on this chart is

\[
\boxed{
T_*=\frac{3\mathsf A K^3L}{H^3}.
}                                                     \tag{5.2}
\]

The checker verifies (5.2) by the denominator-free identity

\[
81\mathsf A^4-54A_6\mathsf A^2-24\mathsf A L^3
+9A_6^2-12C_6L^2=0.
\]

Define

\[
\Theta=
\frac{27\mathsf A^3-9A_6\mathsf A-2L^3}{2L}.
\]

Exact division shows that \(\Theta\in\mathbb Q[U,V]\), and direct
differentiation gives

\[
\boxed{
P'_{\alpha,\beta}(T_*)
=\frac{8K^9L^4\Theta}{H^9}.
}                                                     \tag{5.3}
\]

The square-discriminant orientation also pulls back completely.  With

\[
\rho_s=V^2+3V+9,\qquad
\sigma_s=\sigma(U,V),
\]

one has

\[
\boxed{
\Omega(\alpha,\beta)
=\frac{1728\rho_s\sigma_s\Theta K^{18}L^8}{H^{18}}.
}                                                     \tag{5.4}
\]

In particular,

\[
\frac{\Omega(\alpha,\beta)}
     {P'_{\alpha,\beta}(T_*)}
=\frac{216\rho_s\sigma_sK^9L^4}{H^9}.                \tag{5.5}
\]

This exposes the missing compatibility.  The ordinary suspension

\[
(U,V,Q)\longmapsto
\left(\alpha,\beta,\frac{Q}{P'_{\alpha,\beta}(T_*)}\right)
\]

has Jacobian

\[
\boxed{
\frac{4K^3L}{H^3}\,
\frac1{P'_{\alpha,\beta}(T_*)}
=\frac{H^6}{2\Theta K^6L^3},
}                                                     \tag{5.6}
\]

which is not constant.  The reciprocal unit actually required by the
ordinary \((U,V)\)-chart is

\[
\left(
\det\frac{\partial(\alpha,\beta)}{\partial(U,V)}
\right)^{-1}
=\frac{H^3}{4K^3L}.                                  \tag{5.7}
\]

Thus the current root-incidence parametrization returns exactly to the old
\(K^3L/H^3\) cone ledger.  The rank-four derivative split is not a
standalone escape from that ledger: a successful construction must also
replace the affine root-incidence chart or realize (5.7) through two
source-dependent polynomial masks.

## 6. The correct chart unit in the root basis

The comparison in Section 5 identifies the function that must actually be
realized:

\[
\mathcal J^{-1}
=\frac{H^3}{4K^3L}.
\]

Four-branch interpolation in the root basis gives an exact formula

\[
\boxed{
\mathcal J^{-1}
=\frac{
3B_0p_0+B_0p_1T+p_2T^2-p_3T^3
}{
72B_0^2\rho\sigma
}.
}                                                     \tag{6.1}
\]

Here

\[
\deg(p_0,p_1,p_2,p_3)=(9,6,6,3);
\]

their explicit integer coefficients are recorded in
[`verify_a4_chart_unit_rank_four.py`](../scripts/verify_a4_chart_unit_rank_four.py).
The checker derives the four conjugate values independently on the
\((s,t)\)-presentation and verifies every coefficient of (6.1) by exact
polynomial identities.

The denominator is structured:

\[
\boxed{B_0^2\rho\sigma.}                              \tag{6.2}
\]

Moreover \(B_0,\rho,\sigma\) are irreducible and pairwise coprime over
\(\mathbb Q[\alpha,\beta]\).  Thus (6.2) has three genuine target boundary
components, but it admits the formal two-mask grouping

\[
B_0^2\mid\rho\sigma.
\]

Write

\[
\mathcal N(T)=3B_0p_0+B_0p_1T+p_2T^2-p_3T^3.
\]

On the localized fourfold chart, introduce two mask variables \(Q_1,Q_2\)
and use

\[
\boxed{
X=\frac{Q_1}{B_0^2},\qquad
Y=\frac{\mathcal N(T)Q_2}{72\rho\sigma}.
}                                                     \tag{6.3}
\]

Then

\[
\det\frac{\partial(X,Y)}{\partial(Q_1,Q_2)}
=\mathcal J^{-1}.                                    \tag{6.4}
\]

Consequently the rational map

\[
(U,V,Q_1,Q_2)\longmapsto(\alpha,\beta,X,Y)
\]

has determinant one and retains the original degree-four \(A_4\) function
field: the two mask variables are recovered uniquely at the generic point.
This is the exact localized two-mask Keller suspension sought by the
rank-four calculation.

It is still not polynomial.  The first mask has a genuine \(B_0^2\)-pole.
The coefficient \(p_3\) is coprime to both \(\rho\) and \(\sigma\), so the
second mask also has genuine poles along both remaining components.
Equivalently, integrating \(\mathcal N(T)\) to make a triangular incidence
coordinate leaves a \(T^3\)-coefficient \(p_2/3\), with

\[
\gcd(p_2,B_0)=1.
\]

No target-only redistribution removes these poles.

### 6.1. Quick local viability screen

The reduced boundary arrangement is not simple normal crossings.  Exact
Gröbner calculations give

\[
\operatorname{Sing}(\sigma)=V(a^2,\rho)
\]

and, for every pair chosen from \(B_0,\rho,\sigma\), the nontransverse
intersection scheme is again

\[
V(a^2,\rho).
\]

All three components meet, with triple-intersection ideal

\[
(B_0,\rho,\sigma)=(a^3,\rho).                        \tag{6.5}
\]

This scheme consists of the two conjugate points

\[
a=0,\qquad \beta^2+3\beta+9=0.
\]

The local structure is explicit.  Put

\[
z=\rho,\qquad c=2\beta+3.
\]

Since \(c^2+27=4z\), the function \(c\) is a unit along \(z=0\).
Reduction by \(\beta^2+3\beta+9-z\) gives

\[
\begin{aligned}
B_0&=a^3-3az+cz,\\
\sigma&=ca^3-3a^2z+z^2.                              \tag{6.6}
\end{aligned}
\]

The coefficient \(p_3\) controlling the \(T^3\)-term of
\(\mathcal N(T)\) has local form

\[
p_3=4ca^2-16az+63a+4cz.                              \tag{6.7}
\]

At the common cluster \((a,z)=(0,0)\),

\[
\operatorname{ord}(B_0^2)=2,\qquad
\operatorname{ord}(\rho\sigma)=3,\qquad
\operatorname{ord}(p_3)=1.                           \tag{6.8}
\]

Thus the numerator supplies only one of the three local orders required by
the second mask.  The simplest two-normal-crossings reciprocal chart fails
this viability test.

This is not a theorem excluding every resolved affine modification.  Such a
modification could blow up the common cusp/tangency cluster and distribute
the resulting exceptional valuations.  It does show that the apparent
two-factor grouping in (6.2) hides a higher-contact three-component
boundary package; realizing it polynomially would require a nontrivial
resolution-aware construction.

### 6.2. Resolution and the full numerator

The coefficient order in (6.8) is only a first screen.  On the root
incidence, \(T\) also vanishes above the cluster, so the divisorial order of
the full numerator

\[
 {\cal N}(T)=3B_0p_0+B_0p_1T+p_2T^2-p_3T^3
\]

must be computed in the rank-four algebra.

Four ordinary point blowups give an embedded resolution.  Blow up
\((a,z)=(0,0)\), then the remaining common point in the \(a\)-chart.
There are then two different unresolved points: the \(B_0/\rho\) contact
and the corner through which the strict transform of \(\sigma\) passes.
Blowing up both gives exceptional valuations

\[
\begin{aligned}
 E_1&=(1,1),&E_2&=(1,2),&
 E_3&=(1,3),&F&=(2,3)
\end{aligned}
\]

in the coordinates \((a,z)\).  The exceptional chain is

\[
 E_1-F-E_2-E_3.
\]

The strict transform of \(\sigma\) meets \(F\), while those of \(B_0\) and
\(\rho\) meet \(E_3\) at two distinct points.  Thus the reduced total
boundary is simple normal crossings.  The complete exceptional parts of
the two denominator divisors are

\[
\begin{aligned}
\pi^*(B_0^2)
 &=2\widetilde B_0+2E_1+4E_2+6E_3+6F,\\
\pi^*(\rho\sigma)
 &=\widetilde\rho+\widetilde\sigma
   +3E_1+5E_2+6E_3+9F.                 \tag{6.10}
\end{aligned}
\]

To transform \({\cal N}\), form its characteristic polynomial in the
rank-four root algebra:

\[
 \operatorname{Res}_T(P_{\alpha,\beta}(T),X-{\cal N}(T)).
\]

Its Newton polygons give the following branchwise orders.  Multiplicity
notation records the number of roots with the displayed order.

| divisor | \(v(B_0^2)\) | \(v(\rho\sigma)\) | root orders of \({\cal N}\) | minimum residual pole of \({\cal J}^{-1}\) |
|---|---:|---:|---:|---:|
| \(E_1\) | 2 | 3 | \(6^1,4^3\) | \(0^1,1^3\) |
| \(E_2\) | 4 | 5 | \(8^1,7^3\) | \(1^1,2^3\) |
| \(E_3\) | 6 | 6 | \(9^4\) | \(3^4\) |
| \(F\) | 6 | 9 | \(15^1,12^3\) | \(0^1,3^3\) |

Here the residual pole is

\[
 \delta_D=\max\{v_D(B_0^2)+v_D(\rho\sigma)-v_D({\cal N}),0\}. \tag{6.11}
\]

The strict transforms add three further facts.  Along \(B_0\), every
target-normalized numerator order is \(3/2\); after the ramified root
normalization this leaves one integral mask order on every branch.  Along
\(\rho\), the orders are \(1,0,0,0\), so three branches retain one mask
order.  Along \(\sigma\), all four orders are one, so the full numerator
cancels \(\widetilde\sigma\) divisorially.  This last statement does not
contradict the coefficient gcd test: regularity on the normalized
incidence need not imply membership in the original nonnormal polynomial
root algebra.

The two-mask allocation is now an exact interval problem.  At a resolved
branch write

\[
 d_1=v(B_0^2),\qquad d_2=v(\rho\sigma),\qquad
 n=v({\cal N}).
\]

If \(x\) numerator orders are assigned to the first coordinate, the
remaining source-mask orders are

\[
 m_1=d_1-x,\qquad m_2=d_2-(n-x).                     \tag{6.12}
\]

For \(n<d_1+d_2\), every minimal allocation lies in

\[
 \max(0,n-d_2)\leq x\leq\min(n,d_1)
\]

and satisfies \(m_1+m_2=d_1+d_2-n=\delta_D\).  Thus the table solves the
divisorial allocation problem, but also proves that no allocation removes
all residual masks.  A convenient local choice assigns enough of
\({\cal N}\) to cancel \(d_2\) on every exceptional chart; all exceptional
residual orders then lie in the first mask.  The strict \(B_0\) order is
also forced into the first mask, while the three uncancelled \(\rho\)
branches are forced into the second.

This gives a first contraction screen.  The residual orders differ among
root branches over the same \(E_1,E_2\), and \(F\).  Consequently they
cannot be supplied by monomials pulled back from the resolved target.  A
polynomial contraction would have to construct root-dependent principal
Cartier masks on the normalized incidence, prove that they glue across the
four charts, and only then contract the exceptional chain.  No such
principal masks are constructed here.

### 6.3. The forced \(\rho\)-selector and polynomial-descent obstruction

The branch-dependent part of the second mask has a compact equation.
Reduction of the quartic modulo \(\rho\) gives

\[
 \boxed{
 P_{\alpha,\beta}(T)\equiv
 (T-3a^3)(T+a^3)^3\pmod{\rho}.}                     \tag{6.13}
\]

Thus

\[
 G=T+a^3                                               \tag{6.14}
\]

is the reduced equation of the ramified triple component over
\(\rho=0\), while \(T-3a^3\) gives the simple component.  The full
numerator distinguishes the two in the required direction:

\[
\begin{aligned}
 {\cal N}(3a^3)&\equiv0\pmod\rho,\\
 {\cal N}(-a^3)&\equiv144a^9(2\beta+3)\pmod\rho.
                                                               \tag{6.15}
\end{aligned}
\]

Consequently the chart unit is regular on the simple component and has its
uncancelled strict-\(\rho\) pole on the triple component.

The characteristic polynomial of \(G\) gives uniform exceptional orders

\[
\begin{array}{c|cccc}
 &E_1&E_2&E_3&F\\ \hline
 v(G)&1&2&3&3 .
\end{array}                                             \tag{6.16}
\]

On strict \(\rho\), its target-normalized root orders are
\((1/3,1/3,1/3,0)\); after normalization of the ramified component these
become one on the triple branch and zero on the simple branch.  Hence the
strict triple component is Cartier on the resolved charts.  Its local
equations are obtained from \(G\) by removing the exceptional factors:

\[
 \frac Ga,\qquad \frac G{a^2},\qquad
 \frac G{a^3},\qquad \frac G{u^3}                    \tag{6.17}
\]

on the \(E_1,E_2,E_3,F\) charts, respectively.

These quotients do not descend to the original polynomial root algebra.
Indeed, in

\[
 {\cal R}=\mathbb Q[\alpha,\beta,T]/(P_{\alpha,\beta}),
\]

the height-one prime of the triple component is

\[
 {\mathfrak p}_3=(\rho,T+a^3).
\]

It lies inside the cluster maximal ideal

\[
 {\mathfrak m}=(a,\rho,T).                            \tag{6.18}
\]

Every polynomial element vanishing on the forced triple component
therefore vanishes at the cluster.  Since \(a,\rho,T\) all have positive
order on every branch above the first exceptional divisor, such an element
has positive order on all four \(E_1\)-branches.

This contradicts the exact residual vector

\[
 \delta_{E_1}=(0,1,1,1).
\]

Indeed, if two regular polynomial masks had the required total divisor,
their product would vanish on \({\mathfrak p}_3\), hence lie in
\({\mathfrak m}\), but its order on the distinguished \(E_1\)-branch would
have to be zero.  Therefore

\[
\boxed{
\text{no two masks in the original polynomial root algebra realize the
exact resolved residual divisor.}}                    \tag{6.19}
\]

This is an all-degree local obstruction to direct polynomial descent.  It
does not exclude a genuinely new affine modification whose coordinate ring
adjoins the exceptional quotients in (6.17).

### 6.4. The forced affine-modification chain

The exceptional quotients can be followed exactly.  First adjoin

\[
 q=\frac{T+a^3}{a^3}.
\]

After eliminating \(T=a^3(q-1)\), this gives one affine hypersurface.
On \(a=0\), its equation is a unit times

\[
 (c^2+27)^4,\qquad c=2\beta+3.
\]

Thus \(a=0,\ c^2+27=0\), with \(q\) arbitrary, is a
codimension-one singular cylinder.  The direct affine-modification ring
fails \(R_1\), is nonnormal, and cannot be a polynomial ring.

The first integral correction is forced:

\[
 w=\frac{\rho}{a^3}.                                  \tag{6.20}
\]

Using \(c^2+27=4a^3w\), substituting
\(T=a^3(q-1)\), and dividing the exact total order \(a^{12}\) gives the
\(E_3\) strict-transform complete intersection.  It still has the two
conjugate singular points

\[
 a=0,\qquad q=1,\qquad 27w=c,\qquad c^2+27=0.         \tag{6.21}
\]

The other end of the fan is equally explicit.  On the \(F\)-chart put

\[
 a=u^2k,\qquad \rho=u^3k,\qquad
 T=u^3r-u^6k^3,
\]

so that \(G=T+a^3=u^3r\).  If \(F_F\) denotes the strict-transform
equation after division by its exact total order \(u^{12}\), then on the
exceptional divisor

\[
\boxed{
8F_F\big|_{u=0}
=
\bigl(2r-(27-3c)k\bigr)
\bigl(2r-(c-9)k\bigr)^3.}                            \tag{6.22}
\]

Hence the exceptional fiber is one simple line plus one triple line.  The
total affine surface is singular along the generic point of the triple
line, again violating \(R_1\).  Its next integral quotient is forced by the
Newton edge:

\[
 s=\frac{2r-(c-9)k}{u}.                               \tag{6.23}
\]

After adjoining \(s\) and dividing the exact new total order \(u^3\), the
resulting chart is still singular at

\[
 u=k=s=0,\qquad c^2+27=0.
\]

This is the adjacent fan center.  Repeating the process does not produce
one smooth affine chart; it reconstructs the neighboring charts of the
four-blowup resolution.

There are therefore two endpoints to the canonical construction:

- retaining one affine quotient chart leaves a nonnormal or singular
  surface;
- adjoining every chart gives the smooth resolution with exceptional
  chain
  \[
  E_1-F-E_2-E_3.
  \]

In that chain the self-intersections are \((-3,-1,-3,-1)\), and its
intersection matrix has leading principal minors

\[
 -3,\quad2,\quad-3,\quad1.
\]

It is negative definite and contracts back to the cluster, but the smooth
resolution contains complete exceptional curves and is not affine.
Moreover, proper birational descent over the normal root incidence gives
no new global regular functions, so the local Cartier quotients do not
become global polynomial masks on the full resolution.

This is the final conclusion of the present construction:

\[
\boxed{
\begin{gathered}
\text{the affine exceptional-quotient charts are singular or nonnormal,}\\
\text{while the smooth full-chain modification is nonaffine.}
\end{gathered}}
                                                               \tag{6.24}
\]

A polynomial counterexample cannot arise from the complete Rees/affine
modification dictated by (6.17).  A surviving model would have to remove a
codimension-one divisor from the resolved space, thereby introducing a new
boundary component which is absent from the current determinant ledger.

### 6.5. A corrected exceptional selector and exact rational masks

The deletion requirement can be used constructively.  On the first
exceptional chart put

\[
 z=ay,\qquad T=at.
\]

Modulo \(c^2+27=0\), the exceptional quartic factors exactly as

\[
 \bigl(2t-(27-3c)y\bigr)
 \bigl(2t-(c-9)y\bigr)^3.                            \tag{6.25}
\]

Thus \(H=2T-(27-3c)\rho\) selects the simple exceptional branch.  Its
order on that branch is still one too small at \(F\).  Expanding the simple
\(F\)-branch one order further gives the first correction

\[
\boxed{
 \widehat H
 =
 4c\bigl(2T-(27-3c)\rho\bigr)
 -27(c-9)a\rho.}                                    \tag{6.26}
\]

The characteristic polynomial of \(\widehat H\) has branch orders

\[
\begin{array}{c|cccc}
 &E_1&E_2&E_3&F\\ \hline
 v(\widehat H)
 &(2,1,1,1)&(3,2,2,2)&(3,3,3,3)&(6,3,3,3).
\end{array}                                          \tag{6.27}
\]

Since \(T\) and \(G=T+a^3\) both have uniform exceptional orders
\((1,2,3,3)\), subtraction of (6.27) gives

\[
\boxed{
 M=\frac{T(T+a^3)}{\widehat H}}                      \tag{6.28}
\]

with exceptional orders

\[
\begin{array}{c|cccc}
 &E_1&E_2&E_3&F\\ \hline
 v(M)
 &(0,1,1,1)&(1,2,2,2)&(3,3,3,3)&(0,3,3,3).
\end{array}
\]

These are exactly the residual orders in (6.11), including the previously
missing simple-\(F\) cancellation.  Along the strict transforms, \(T\) has
target-normalized order \(1/2\) on all four \(B_0\)-branches, \(G\) has
order \(1/3\) on the three ramified \(\rho\)-branches, and
\(\widehat H\) is a unit on \(B_0,\rho,\sigma\).  After root normalization,
(6.28) therefore has exactly the complete residual mask divisor.

There is also an exact two-factor allocation:

\[
\boxed{
 M_1=\frac{B_0T}{\widehat H},\qquad
 M_2=\frac{T+a^3}{B_0}.}                             \tag{6.29}
\]

On every exceptional branch, \(M_1\) has the full residual vector and
\(M_2\) has order zero.  Both are regular after deleting the strict
divisors \(B_0=0\) and \(\widehat H=0\).  Thus the earlier obstruction has
identified, rather than merely demanded, a new boundary divisor.

The norm

\[
 h(a,\beta)=\operatorname{Res}_T(P,\widehat H)
\]

is irreducible of total degree \(16\) over \(\mathbb Q\).  Its exceptional
multiplicities in chain order \(E_1,F,E_2,E_3\) are
\((5,15,9,12)\), so its strict transform meets \(F\) once and \(E_3\)
three times.  On the resolved coefficient plane, the four classes

\[
 F,\quad E_3,\quad\widetilde B_0,\quad\widetilde h
\]

form a unimodular basis of the relative Picard lattice.  Moreover,

\[
 F+E_3+\widetilde B_0+2\widetilde h
\]

has intersections \((1,1,2,6)\) with \(E_1,F,E_2,E_3\).  Hence this
four-divisor deletion passes the coarse relative unit, class-group, and
ampleness screens.

This is a candidate, not an affine-space theorem.  The class calculation
uses the norm divisor on the resolved coefficient plane.  The next
required step is to normalize the root incidence over this deletion,
compute its full divisor-class exact sequence, and identify its coordinate
ring.  Unimodularity downstairs does not by itself prove that the
normalized open is factorial or isomorphic to affine space.

### 6.6. Normalized-incidence class obstruction

The upstairs divisor count is decisive.  On both retained exceptional
divisors \(E_1\) and \(E_2\), the exceptional quartic has the factorization
in (6.25).  Centering the triple line gives coefficient orders

\[
 (1,1,1,0,0)
\]

in powers \(0,1,2,3,4\) of the shifted root coordinate.  Its lower Newton
edge has slope \(-1/3\).  Consequently the normalization has a simple
exceptional prime and a distinct triple exceptional prime above each of
\(E_1,E_2\).  The proposed open therefore retains at least four prime
exceptional curves.

The two horizontal deletions do not split enough to remove those four
class directions.  Along strict \(B_0\), all four roots have
target-normalized order \(1/2\).  After writing \(T=B_0^{1/2}U\), the
residual polynomial is

\[
 (U^2-3A_0)^2.
\]

The quadratic is irreducible over the function field of \(B_0\).  Indeed,
locally

\[
 B_0=a^3+(c-3a)z,\qquad A_0=a^3-(\beta+6)z,
\]

and on the normalization of \(B_0\)

\[
 A_0
 =a^3\frac{c-3a+\beta+6}{c-3a}.
\]

Its leading coefficient \(3(\beta+3)/c\) is a unit along \(\rho=0\), so
\(A_0\) has odd order three and is not a square.  The residue degree and
ramification degree are both two, exhausting the quartic degree.  Thus
strict \(B_0\) has one prime upstairs.  Likewise \(\widehat H\) is linear
in \(T\), with generic leading coefficient \(8c\), and its norm is
irreducible; hence it supplies one prime.

Pass to any smooth resolution of the normalized pullback.  Exceptional
curves over a normal surface point have negative-definite intersection
matrix and independent relative divisor classes.  Components deleted over
\(F,E_3\) remove at most their own class directions.  Any additional
resolution curve lying inside the deleted boundary adds one class and one
deleted component simultaneously.  After those cancellations, four
independent retained exceptional directions remain, while the strict
\(B_0,\widehat H\) divisors can kill at most two.  The divisor exact
sequence therefore gives

\[
\boxed{
 \operatorname{rank}\operatorname{Cl}(U)\geq2,
}                                                     \tag{6.30}
\]

for the normalized open \(U\) defined by the proposed four-divisor
deletion.  In particular,

\[
\boxed{
 U\not\cong\mathbb A^2,\qquad
 U\times\mathbb A^m\not\cong\mathbb A^{m+2}.
}                                                     \tag{6.31}
\]

Thus (6.28)--(6.29) solve the mask divisor exactly but do not produce an
affine-space source.  The corrected deletion candidate is closed.  Within
this resolved model, any further affine-space attempt must delete at least
two additional horizontal prime divisors, or replace the model so that the
simple/triple components over \(E_1,E_2\) are not retained.  Either choice
changes the determinant ledger again.

### 6.7. The normalized cluster and the index-two terminal obstruction

The normalized exceptional divisor can be computed completely.  In chain
order

\[
 S_1-F_s-S_2-Q-R_2-F_t-R_1,
\]

where \(S_i,R_i\) are the simple and triple components over \(E_i\), its
self-intersections are

\[
 (-3,-1,-3,-4,-1,-3,-1).                           \tag{6.32}
\]

The determinant of this chain is \(-4\), and its Smith form is

\[
 \operatorname{diag}(1,1,1,1,1,1,4).
\]

Blowing down, in order, \(R_1,R_2,F_t,F_s\), leaves

\[
 \begin{pmatrix}
 -2&1&0\\
 1&-2&1\\
 0&1&-2
 \end{pmatrix}.                                     \tag{6.33}
\]

Thus the cluster is an \(A_3\) rational double point and its local
discriminant group is cyclic of order four.  This also explains why
horizontal components cannot be assigned classes by looking at only one
factor of a principal divisor.  For example,

\[
 \operatorname{Norm}(T+a^3)=\rho\,q_{10},
\]

where \(q_{10}\) is irreducible of degree ten and also passes through the
cluster.

An explicit curvette identifies an odd generator.  Put

\[
 \chi=16a^2-8a\beta-12a+\rho
      =16a^2-4ac+\rho
\]

and

\[
\begin{aligned}
 L_\chi={}&(32a+9-4c)T-13311a^3\\
 &+(2241c-567)a^2+(162c-7290)a .
\end{aligned}
\]

On the first blowup, the strict transform of \(\chi\) meets \(E_1\) at
\(y=4c\).  The exceptional transform of \(L_\chi/a\) is

\[
 -4ct+162c+9t-7290.
\]

It vanishes on the simple root in (6.25) and takes the nonzero value
\(216(c-45)\) on the triple root.  Moreover,

\[
 \operatorname{Norm}(L_\chi)=\chi\,q_{14},
\]

with \(q_{14}\) irreducible of degree fourteen and coprime to \(\chi\).
The \(\chi\)-component is therefore a curvette of \(S_1\) and represents
an odd generator of the order-four discriminant group.

Let \(\Lambda^\vee\) be the full rank-seven local divisor lattice.  In the
exceptional coordinates of (6.32), it has basis

\[
 e_0,e_1,\ldots,e_5,\quad
 q=\left(\frac34,\frac54,\frac12,\frac14,
          \frac12,\frac14,\frac14\right).
                                                               \tag{6.34}
\]

The normalized exceptional valuations of \(B_0\) and \(\widehat H\) are

\[
\begin{aligned}
 b&=(1,3,2,3,6,3,3),\\
 h&=(2,6,3,3,6,3,3).
\end{aligned}
\]

Strict \(B_0\) has ramification multiplicity two, while strict
\(\widehat H\) has multiplicity one.  Their prime classes are consequently
\(-b/2\) and \(-h\).  In the basis (6.34), these are

\[
\begin{aligned}
 [D_B]&=(4,6,2,0,0,0,-6),\\
 [D_{\widehat H}]&=(7,9,3,0,0,0,-12).               \tag{6.35}
\end{aligned}
\]

Already the coordinate gcd of \([D_B]\) is two:

\[
 [D_B]=2(2,3,1,0,0,0,-3).                           \tag{6.36}
\]

Thus this irreducible boundary prime is nonprimitive and cannot occur in
any basis of \(\Lambda^\vee\).  Consistently, the gcd of all
\(2\times2\) minors of the matrix formed by the two classes in (6.35) is
also exactly two.

This closes the proposed boundary-enlargement repair, not just the original
four-divisor choice.  In the local divisor exact sequence, simultaneous
triviality of boundary units and of the relative class group requires the
boundary-class map to be an isomorphism.  Its columns would form a basis
and would include the forced strict \(B_0\) class, which is impossible.
Equivalently, every rank-seven boundary matrix containing this column has
even determinant.  Therefore

\[
\boxed{
\text{no resolved boundary enlargement retaining the exact masks
\(B_0,\widehat H\) passes the affine-space unit/class test.}
}                                                     \tag{6.37}
\]

The obstruction is the order-two quotient inside the \(A_3\)
discriminant lattice.  Adding simple- or triple-branch selectors cannot
remove it; the defect already belongs to the forced strict \(B_0\) prime.
A continuation must change the exact mask presentation so that \(B_0\)
is no longer a boundary prime.

### 6.8. A \(B_0\)-free curvette split

The odd curvette above provides exactly such a change.  The same total mask
(6.28) factors as

\[
\boxed{
 M_1^\chi=\frac{T L_\chi}{\widehat H},\qquad
 M_2^\chi=\frac{T+a^3}{L_\chi}.
}                                                     \tag{6.38}
\]

The characteristic root orders of \(L_\chi\) are

\[
\begin{array}{c|cccc}
 &E_1&E_2&E_3&F\\ \hline
 v(L_\chi)&(1,1,1,1)&(1,1,1,1)&(1,1,1,1)&(2,2,2,2).
\end{array}
\]

In the normalized seven-curve order (6.32), this is

\[
 \ell=(1,2,1,1,3,2,3).
\]

Using \(t=(1,3,2,3,6,3,3)\) for the normalized orders of both
\(T\) and \(T+a^3\), the two masks in (6.38) have exceptional orders

\[
\begin{aligned}
 v(M_1^\chi)&=t+\ell-h=(0,-1,0,1,3,2,3),\\
 v(M_2^\chi)&=t-\ell=(0,1,1,2,3,1,0).               \tag{6.39}
\end{aligned}
\]

Thus the only exceptional pole is the simple \(F\)-component \(F_s\).
The horizontal poles are the strict \(\widehat H\) divisor and the two
components

\[
 D_\chi,\qquad D_{14}
\]

of \(\operatorname{Norm}(L_\chi)=\chi q_{14}\).  Exact gcd calculations
show that neither denominator shares a horizontal component with its
numerator.  In particular, strict \(B_0\) is absent.

The two curvette classes are

\[
\begin{aligned}
 [D_\chi]&=(0,0,0,0,0,0,-1),\\
 [D_{14}]&=(8,13,5,2,3,1,-11)
\end{aligned}
\]

in the basis (6.34).  They meet \(S_1\) and \(R_1\), respectively.
The four forced classes

\[
 F_s,\quad D_{\widehat H},\quad D_\chi,\quad D_{14}
\]

span a primitive rank-four sublattice.  Among completions by three
exceptional primes, the unique unimodular completion is

\[
\boxed{
 F_s,\ R_2,\ F_t,\ R_1,\
 D_{\widehat H},\ D_\chi,\ D_{14}.
}                                                     \tag{6.40}
\]

This support also passes the positivity test.  The effective divisor

\[
 F_s+R_2+2F_t+6R_1
 +2D_{\widehat H}+D_\chi+5D_{14}
\]

has intersections

\[
 (2,1,1,7,1,1,1)
\]

with \(S_1,F_s,S_2,Q,R_2,F_t,R_1\).  Therefore the \(B_0\)-free
curvette split passes the full normalized local unit, class, and relative
ampleness screens.

This is the first resolved candidate to survive those three tests.  It is
not yet an affine-space identification or a polynomial Keller map.  The
remaining issue is now concrete: construct the affine complement of
(6.40), calculate its coordinate ring, decide whether it is
\(\mathbb A^2\), and express both quotients in (6.38) in those coordinates.

### 6.9. Positive-genus boundary obstruction

The horizontal boundary settles that question without a full presentation
of the coordinate ring.  Since \(\widehat H\) is linear in \(T\), its
prime divisor on the root incidence is birational to its coefficient-plane
norm curve

\[
 h(a,b)=\operatorname{Res}_T(P,\widehat H)=0.
\]

Section 6.5 already proves that \(h\) is irreducible of degree sixteen.
Exact normalization of its projective closure gives

\[
\boxed{g(\widetilde{V(h)})=13.}                     \tag{6.41}
\]

Thus \(D_{\widehat H}\) is a positive-genus boundary component of the
smooth resolved complement (6.40).

The other nonlinear horizontal prime independently fails the same test.
The selector \(L_\chi\) is also linear in \(T\), so \(D_{14}\) is
birational to its coefficient-plane norm component.  Exact division and
normalization give

\[
 \operatorname {Norm}(L_\chi)=\chi q_{14},\qquad
 \deg q_{14}=14,\qquad
 \boxed{g(\widetilde{V(q_{14})})=20.}                \tag{6.42}
\]

This is incompatible with an affine plane.  Indeed, suppose that the
complement were isomorphic to \(\mathbb A^2\).  Complete the resolved
surface smoothly and compare it with
\(\mathbb P^2\supset\mathbb A^2\).  A common resolution of the induced
birational map is obtained by blowing up boundary points.  The strict
transform of \(D_{\widehat H}\) still has genus thirteen.  Under the
birational morphism to \(\mathbb P^2\), it must either:

- map birationally onto the line at infinity, which would force genus zero;
  or
- be exceptional, whereas every exceptional component of a birational
  morphism between smooth surfaces is rational.

Both alternatives are impossible.  Therefore

\[
\boxed{
 U_\chi\not\cong\mathbb A^2.
}                                                     \tag{6.43}
\]

The finite-field behavior provides an independent diagnostic: at good
primes where \(\rho=0\) has no rational cluster point, the resolved open
has the same rational points as the root incidence with
\(\widehat H L_\chi\ne0\), and its counts already differ from \(p^2\).
Those counts are evidence only; either genus computation and the completion
argument prove (6.43) in characteristic zero.

Consequently the corrected total mask (6.28) is closed, not merely its
first factorization.  Every factorization retains the genuine strict
\(\widehat H\) pole, hence the genus-thirteen divisor at infinity; the
curvette split also retains the genus-twenty divisor \(D_{14}\).  A
polynomial affine-space construction must replace the corrected selector
itself by a principal mask whose horizontal prime components are rational,
or leave this root chart entirely.

### 6.10. First bounded genus-zero selector search

There is a smaller selector space than (6.26), but not a smaller horizontal
norm.  Give \(a,b,T\) ordinary degree one.  Let \(V_{d,r}\) be the vector
space of polynomials of total degree at most \(d\), root degree at most
\(r\), and normalized exceptional orders at least

\[
 h=(2,6,3,3,6,3,3)                                  \tag{6.44}
\]

in the seven-curve order
\((S_1,F_s,S_2,Q,R_2,F_t,R_1)\).  Equivalently, their four unsplit root
order vectors are at least those in (6.27).

The valuation ideals can be computed by exact truncated root jets.  Write
\(z=\rho\), \(c=2b+3\), so \(c^2+27=4z\).  On the simple components the
jets needed below the thresholds (6.44) are

\[
\begin{array}{c|c|c}
S_1&a=x,\ z=xy&
T=x(27-3c)y/2+O(x^2)\\
S_2&a=x,\ z=x^2y&
T=x^2(27-3c)y/2+O(x^3)\\
F_s&a=u^2k,\ z=u^3k&
\displaystyle
T=u^3\left(
\frac{(27-3c)k}{2}
+\frac{9(c+3)k^2}{8}u^2
\right)+O(u^6).
\end{array}                                          \tag{6.45}
\]

Reduction of every forbidden coefficient modulo \(c^2+27\) gives exact
rational linear algebra.  The result is

\[
\begin{aligned}
 V_{2,1}&=0,&
 V_{2,2}&=\langle T^2\rangle,\\
 V_{3,1}&=\langle q_0,q_1,q_2,q_3,q_4,q_5\rangle,&
 V_{3,2}&=V_{3,1}+\langle T^2,aT^2,bT^2\rangle,
                                                        \tag{6.46}
\end{aligned}
\]

where \(q_0=a^3,\ q_5=a^2T\), and, after harmless scalar
normalizations,

\[
\begin{aligned}
4q_1={}&4bT+81ab^2+243ab+729a\\
       &-72b^3-324b^2-972b-972,\\
4q_2={}&4b^2T+36T-243ab^2-729ab-2187a\\
       &+216b^3+972b^2+2916b+2916,\\
3q_3={}&3aT+4T-54ab^2-162ab-486a+12b^3-324,\\
q_4={}&abT-8T+27ab^2+81ab+243a-24b^3+648.
                                                        \tag{6.47}
\end{aligned}
\]

Each of \(q_1,q_2,q_3,q_4\) has the exact four branch vectors (6.27).
Thus the corrected degree-four selector was not minimal as a root-algebra
function.

The norm filtration supplies the useful bounded obstruction.  At infinity
give \(T\) weight three and \(a,b\) weight one.  The leading quartic has
weight twelve and nonzero constant term.  For

\[
 S=C(a,b)T^2+D(a,b)T+E(a,b)
\]

of ordinary total degree at most three, a nonzero quadratic-root
coefficient has

\[
 \deg\operatorname {Norm}(S)=24+4\deg C\geq24.
\]

If \(C=0\) and the quadratic part of \(D\) is nonzero, the norm has degree
twenty.  The quadratic parts of the \(T\)-coefficients in (6.46) are the
independent forms \(a^2,ab,b^2\), so they cannot cancel without deleting
the corresponding three basis directions.  What remains is

\[
 \langle q_0,q_1,q_3\rangle.
\]

An exact selector in this space has a nonzero linear part in its
\(T\)-coefficient: otherwise it is a multiple of \(q_0=a^3\), whose four
\(E_1\)-orders are all three.  Its norm therefore has degree exactly
sixteen.  Hence

\[
\boxed{
\begin{gathered}
\text{within total degree at most three and root degree at most two,}\\
\text{no exact selector of norm degree below sixteen exists.}
\end{gathered}}                                      \tag{6.48}
\]

In particular no selector below this horizontal degree can have both the
required valuations and genus zero.  The bound is sharp.  The norms of
\(q_1\) and \(q_3\) are irreducible of degree sixteen, but exact
normalization gives

\[
 g(\widetilde{V(\operatorname {Norm}(q_1))})=12,
 \qquad
 g(\widetilde{V(\operatorname {Norm}(q_3))})=14.      \tag{6.49}
\]

The displayed perturbation \(q_1+a^3\) remains irreducible of genus twelve.
These three genus computations are probes of the sharp degree-sixteen
parameter plane, not a classification of it.

The search order matters here.  Below degree sixteen there is no norm to
factor.  The three tested sharp-bound norms have one positive-genus
irreducible component, so they fail the rational-component gate.  No
primitive boundary-lattice completion, polynomial divisibility test, or
affine reconstruction claim is made for them.  Those later tests should be
run only after a rational point of the degree-sixteen parameter
stratification is found.

## 7. Structured next search

Both corrected-selector factorizations are now closed, and (6.48) gives a
sharp first lower bound for their replacement.  The next search should
start on the degree-sixteen norm stratum before doing any contraction:

1. stratify the projective parameter plane
   \(\mathbb P\langle q_0,q_1,q_3\rangle\) by norm factorization and
   normalization genus;
2. retain only parameters for which every indispensable norm component is
   rational;
3. only then compute the normalized discriminant-lattice span and
   relative-ampleness cone;
4. test polynomial divisibility and construct affine coordinates only for
   a candidate which passes those preceding gates;
5. if the degree-sixteen plane has no rational member, move to the
   degree-twenty root-linear stratum before the root-quadratic
   degree-\(\geq24\) stratum.

Random ambient shears, the uncorrected exceptional quotients, and the
coarse four-divisor contraction are no longer relevant to this branch.
Nor can extra selector primes repair (6.29): the terminal defect is already
present in the class of its forced strict \(B_0\) pole.  The curvette split
removes that defect but exposes independent genus-thirteen and
genus-twenty obstructions.
The next exact frontier is therefore a rational-horizontal selector search,
not polynomial contraction of either existing split.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_root_incidence_derivative_split.py
.venv/bin/python scripts/verify_a4_chart_unit_rank_four.py
.venv/bin/python scripts/verify_a4_two_mask_local_viability.py
.venv/bin/python scripts/verify_a4_affine_modification_obstruction.py
.venv/bin/python scripts/verify_a4_corrected_boundary_selector.py
Singular -q scripts/verify_a4_corrected_boundary_genus.sing
.venv/bin/python scripts/verify_a4_genus_zero_selector_search.py
Singular -q scripts/verify_a4_genus_zero_selector_search.sing
```

The checker verifies the compact inverse basis, the two-coordinate
Jacobian identity, generic root-field recovery, the square-discriminant
specialization, the residual \(\sigma\)-pole, the selected rational root,
and the complete comparison with the ordinary \((U,V)\)-chart ledger.  The
second checker verifies the exact rank-four expansion of the correct chart
unit, its three pairwise-coprime boundary factors, the localized two-mask
determinant-one suspension, the four-blowup resolution, the full
branchwise numerator transforms, and the residual divisor-allocation
intervals.  It also verifies the compact selector \(T+a^3\), its resolved
Cartier transforms, and the obstruction to descending the exact masks into
the original polynomial root algebra.  The third checker computes the
common nontransverse cluster, its local normal forms, and the
coefficientwise order deficit for the simple two-mask chart.  The fourth
checker follows the forced exceptional quotients, verifies the
codimension-one nonnormal loci and subsequent singular centers, and checks
the negative-definite full resolution which is smooth but nonaffine.
The fifth checker derives the corrected simple-branch selector, verifies
its full four-ray transforms, constructs the exact rational two-mask
allocation, and checks the unimodular relatively ample coarse deletion
set.  It then normalizes the retained \(E_1,E_2\) branches, counts the two
horizontal deleted primes, and proves the rank-two relative class
obstruction upstairs.  Finally it computes the normalized seven-curve
chain, contracts it to the \(A_3\) chain, constructs an explicit odd
curvette in its order-four discriminant lattice, and proves that the
forced strict \(B_0\) class has content two.  It then verifies the
\(B_0\)-free curvette factorization,
the unique seven-prime unimodular completion, and an effective boundary
divisor positive on every normalized exceptional curve.  The separate
Singular checker computes geometric genus thirteen for the irreducible
corrected-selector norm, verifies
\(\operatorname {Norm}(L_\chi)=\chi q_{14}\), computes genus twenty for
the degree-fourteen component, and rules out the surviving complement as
an affine plane.  The sixth Python checker computes the complete
total-degree-three valuation spaces in root degrees one and two, verifies
the exact order vectors of four new degree-three selectors, and proves the
sharp horizontal norm-degree floor sixteen.  The final Singular checker
verifies irreducibility and genera twelve and fourteen for two sharp-bound
norms, and genus twelve for one displayed \(a^3\)-perturbation.  It does
not classify the full degree-sixteen parameter plane.

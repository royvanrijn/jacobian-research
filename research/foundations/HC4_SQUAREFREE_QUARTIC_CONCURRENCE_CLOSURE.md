# Concurrence closure for the squarefree quartic denominator

## Status

This note proves `HC4NHM10`, the first exact closure after the squarefree
quartic synchronization frontend `HC4NHM9`.

Let

\[
P=L_1L_2L_3L_4,
\qquad
P^2\mid \det\operatorname{Hess}(h_5),
\tag{0.1}
\]

where the four lines are distinct.  On the clean generic-corank-one packet,
`HC4NHM9` supplies constant directions (v_i) satisfying

\[
D_{v_i}h_5\in
\begin{cases}
(L_i^2),&L_i(v_i)=0,\\
(L_i^3),&L_i(v_i)\ne0.
\end{cases}
\tag{0.2}
\]

The present note proves two exclusions.

1. If all four lines are concurrent, all sixteen flag patterns are empty.
2. If exactly three lines are concurrent and the flag on the fourth line is
   transverse, all eight remaining patterns on the concurrent triple are
   empty.

Thus the pencil arrangement is closed completely.  In the exactly-three-
concurrent arrangement only the eight patterns with a tangent fourth flag
remain at this stage.  The subsequent
[general-position closure](HC4_SQUAREFREE_QUARTIC_GENERAL_POSITION_CLOSURE.md)
proves `HC4NHM11` and closes the no-three-concurrent arrangement completely.
The later
[tangent-fourth closure](HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md)
proves `HC4NHM12` and eliminates these last eight rows.

Replay every determinant and polar-syzygy identity with

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_concurrence_closure.py
~~~

The theorem is an exclusion of part of the clean squarefree-denominator
packet.  It does not address positive-defect or lower-Smith components.

## 1. The all-concurrent arrangement

Let the four lines meet at (p=[0:0:1]), and write

\[
C=\operatorname{Hess}(h_5).
\]

Equation (0.2) implies (C(p)v_i=0) for every (i).  If (C(p)) has rank
two, all four directions coincide.  Then (D_vh_5) is divisible by the
product of four line squares, whose degree is eight, although (D_vh_5) has
degree four.  Hence (D_vh_5=0), (h_5) is a cone, and (det C=0).

Consequently every hypothetical nonzero determinant has

\[
\operatorname{rank}C(p)\leq1.
\tag{1.1}
\]

Also, (0.1) gives

\[
\operatorname{ord}_p(\det C)\geq8.
\tag{1.2}
\]

The rest of the pencil closure is a short local Hessian ladder.

### 1.1 Rank one with nonzero pure (z)-jet

Write a homogeneous quintic at (p) as

\[
h_5=A z^5+z^4F_1+z^3F_2+z^2F_3+zF_4+F_5,
\tag{1.3}
\]

where (F_j) is binary of degree (j).  If (A\ne0), rank one of the
constant Hessian says that the (F_1,F_2) part is the two-jet of a fifth
power.  A change (z\mapsto z+\ell(x,y)), which fixes (p) and preserves
the line pencil, therefore gives

\[
h_5=A z^5+z^2F_3+zF_4+F_5.
\tag{1.4}
\]

The order-two determinant face is

\[
20A z^7\det\operatorname{Hess}_{x,y}(F_3).
\]

Thus either (F_3=0), or the binary zero-Hessian theorem makes (F_3) a
cube.  In the latter case take (F_3=\alpha x^3).  Write

\[
F_4=\sum_{i=0}^4 f_i x^{4-i}y^i,
\qquad
F_5=\sum_{i=0}^5 g_i x^{5-i}y^i.
\]

The successive faces below order eight give

\[
f_2=f_3=f_4=0,
\qquad
g_2=\frac{3f_1^2}{4\alpha},
\qquad
g_3=g_4=g_5=0,
\tag{1.5}
\]

then (f_1=0), and finally (g_1=0).  The resulting quintic

\[
A z^5+\alpha x^3z^2+f_0x^4z+g_0x^5
\]

is independent of (y), so its Hessian determinant vanishes.

If (F_3=0) and (F_4\ne0), the first face instead makes (F_4) a fourth
power.  Taking (F_4=\beta x^4), the next two faces give

\[
g_2=g_3=g_4=g_5=0,
\qquad g_1=0,
\]

and again the quintic is a cone.  If (F_3=F_4=0), then

\[
\det\operatorname{Hess}(Az^5+F_5)
=20Az^3\det\operatorname{Hess}_{x,y}(F_5).
\]

Condition (1.2) makes the binary Hessian zero; hence (F_5) is a fifth
power and the result is again a cone.

### 1.2 Rank one with zero pure (z)-jet

Suppose (A=0) but (C(p)) has rank one.  A symmetric rank-one matrix with
zero (zz)-entry has zero (xz,yz)-entries.  After a binary rechart,

\[
h_5=\alpha x^2z^3+z^2F_3+zF_4+F_5,
\qquad \alpha\ne0.
\tag{1.6}
\]

The order-three face makes (F_3) divisible by (x^2), so another shift of
(z) removes it.  The next two faces give

\[
F_4=x^3(r_0x+r_1y),
\qquad
F_5=x^4(s_0x+s_1y).
\tag{1.7}
\]

The exact determinant is

\[
\begin{aligned}
4x^6(&-5\alpha r_1^2z^3-24\alpha r_1s_1xz^2
-24\alpha s_1^2x^2z+3r_0r_1^2x^2z\\
&+8r_0r_1s_1x^3+3r_1^3xyz-5r_1^2s_0x^3
+3r_1^2s_1x^2y).
\end{aligned}
\tag{1.8}
\]

Order at least eight forces (r_1=0), leaving

\[
\det C=-96\alpha s_1^2x^8z.
\tag{1.9}
\]

This is zero or has eighth-order tangent cone (x^8), not the square of
four distinct tangent lines.

### 1.3 Rank zero

Now (C(p)=0), so

\[
h_5=z^2F_3+zF_4+F_5.
\tag{1.10}
\]

The order-five determinant face is

\[
-4z^4F_3\det\operatorname{Hess}_{x,y}(F_3).
\tag{1.11}
\]

If (F_3\ne0), it is a cube; normalize (F_3=\alpha x^3).  The order-six
and order-seven faces give

\[
F_4=f_0x^4+f_1x^3y,
\qquad
F_5=g_0x^5+g_1x^4y+\frac{f_1^2}{4\alpha}x^3y^2,
\tag{1.12}
\]

and the complete determinant becomes

\[
\boxed{
\det C=-\frac{8}{\alpha}x^9(2\alpha g_1-f_0f_1)^2.
}
\tag{1.13}
\]

It is zero or has only one tangent direction.

It remains to consider (F_3=0).  The order-eight face is

\[
[\det C]_{\operatorname{ord}_p=8}
=-\frac43 zF_4\det\operatorname{Hess}_{x,y}(F_4).
\tag{1.14}
\]

This cannot be a nonzero scalar multiple of (P^2) for a squarefree binary
quartic (P).  Indeed, if (x^m\Vert F_4), then

\[
x^{2m-2}\Vert\det\operatorname{Hess}_{x,y}(F_4)
\qquad (m=1,2,3),
\]

with restrictions (-9y^4,-12c^2y^2,-9b^2), respectively, after normalizing
the remaining factor.  Thus a simple root of (F_4) occurs only once in the
product, while a repeated root occurs at least four times.  Neither is the
required multiplicity two.  For (m=4) the binary Hessian is zero.

Equations (1.9), (1.13), and (1.14) exhaust (1.1).  They prove that the
all-concurrent arrangement has no nonzero determinant satisfying (0.1).

## 2. Three concurrent lines and a transverse fourth flag

Normalize

\[
(L_1,L_2,L_3,L_4)=(x,y,x+y,z)
\tag{2.1}
\]

and assume the fourth flag is transverse.  Put (U_i=L_i^{m_i}S_{4-m_i}),
where (m_i=2) for a tangent flag and (m_i=3) for a transverse flag.
For the four patterns on the concurrent triple, exact linear algebra gives

\[
\begin{array}{c|cccc}
(m_1,m_2,m_3;m_4)&RRR;R&TRR;R&TTR;R&TTT;R\\ \hline
\dim\ker(\bigoplus U_i\to S_4)&1&2&4&6.
\end{array}
\tag{2.2}
\]

Crucially, every syzygy has fourth component zero.  If (h_5) were not a
cone, the polar map

\[
v\longmapsto D_vh_5
\tag{2.3}
\]

would be injective.  A relation among the four directions would therefore
give a nonzero polar syzygy.  By (2.2) its fourth coefficient is zero, so the
first three directions are dependent.

Scale that relation to (w_1+w_2+w_3=0), and write its polar components as
(q_i=D_{w_i}h_5).  Mixed partials require

\[
D_{w_1}q_2=D_{w_2}q_1.
\tag{2.4}
\]

The four possible patterns now close as follows.

### 2.1 Three transverse flags

The unique syzygy is

\[
\begin{aligned}
q_1&=-x^3(x+2y),\\
q_2&=y^3(2x+y),\\
q_3&=(x-y)(x+y)^3.
\end{aligned}
\tag{2.5}
\]

For (w_1=(a,b,c)) and (w_2=(d,e,f)), the four coefficients in (2.4) are

\[
4d+2e,qquad6d,qquad6b,qquad2a+4b.
\]

They force (a=b=d=e=0), making (w_1,w_2) dependent, contrary to
injectivity because (q_1,q_2) are independent.

### 2.2 One tangent flag

Up to the symmetry of the concurrent triple, take the tangent flag on
(x=0).  A general syzygy has two parameters (A,B).  Put

\[
w_1=(0,u,v),
\qquad
w_2=(p,q,r).
\]

The four mixed-partial coefficients are

\[
\begin{aligned}
&(A+3B)q+4Bp,\\
&3\bigl((A+3B)p+2(A+B)q\bigr),\\
&3\bigl(2(A+B)p-(3A+B)u\bigr),\\
&-4Au.
\end{aligned}
\tag{2.6}
\]

If (A=0), the first three equations give (p=q=u=0).  If (A\ne0), the
last equation gives (u=0), and the first three again give (p=q=0),
including the exceptional subcase (A+B=0).  This contradicts transversality
of (w_2) to (y=0).

### 2.3 Two tangent flags

Take the tangent flags on (x=0,y=0), and write

\[
w_1=(0,u,v),
\qquad
w_2=(p,0,r).
\]

A general syzygy has four coefficients (A,B,C,D).  Two coefficients of
(2.4) are (3Bp) and (-3Bu).  If (B\ne0), they contradict
transversality of (w_3=-w_1-w_2) to (x+y=0), which requires
(p+u\ne0).  If (B=0), the remaining equations include

\[
Dp=0,qquad Cu=0,
\tag{2.7}
\]

and

\[
2Au+3Cp+9Dp=0,
\qquad
2Ap+6Cp-9Cu+6Dp-3Du=0.
\tag{2.8}
\]

Separating (u=0) and (u\ne0) shows that (p+u\ne0) forces
(A=C=D=0).  The syzygy is zero, a contradiction.

### 2.4 Three tangent flags

Now tangency of (w_3=-w_1-w_2) gives (p=-u).  For the six-parameter
syzygy, the mixed-partial coefficients include

\[
-3Bu,quad -3Du,quad Dr-4Fu,quad -Bv-4Cu,
\tag{2.9}
\]

and

\[
2Au+Br+2Dr-3Eu-6Fu,
\tag{2.10}
\]

\[
-2Au-2Bv-8Cu-Dv-7Eu-2Fu.
\tag{2.11}
\]

If (u\ne0), these equations successively give
(B=D=F=C=E=A=0).  If (u=0), then (p=0).  Injectivity permits only a
relation supported on a repeated tangent pair, for example

\[
D_zh_5=\kappa x^2y^2.
\]

Thus

\[
h_5=\kappa x^2y^2z+F_5(x,y).
\tag{2.12}
\]

For a direction (v_4=(c,d,1)) transverse to (z=0), the coefficient of
(z) in (D_{v_4}h_5) is

\[
2\kappa(cxy^2+dx^2y).
\]

Divisibility by (z^3) forces (c=d=0), after which the constant term is
(\kappa x^2y^2\ne0).  Hence (2.12) cannot satisfy the fourth flag.  The
other two repeated pairs are equivalent by the symmetry of the concurrent
triple.

This proves all eight exclusions with transverse fourth flag.

## 3. Result and remaining squarefree rows

> **Theorem `HC4NHM10` -- Squarefree concurrence closure.**  In the clean
> generic-corank-one squarefree quartic-denominator packet of `HC4NHM9`, the
> all-four-concurrent line arrangement is empty for every tangent/transverse
> flag pattern.  In the exactly-three-concurrent arrangement, every pattern
> whose fourth, nonconcurrent line has a transverse kernel flag is empty.
> The exclusions occur already at the ternary Hessian/polar synchronization
> level, before the cleared Schur-gradient equations are imposed.

After the subsequent `HC4NHM11` general-position exclusion, the remaining
clean squarefree rows are exactly the eight exactly-three-concurrent patterns
with tangent fourth flag.  The all-concurrent cross-ratio is no longer a
parameter of the open problem.  The later theorem `HC4NHM12` closes those
eight rows.

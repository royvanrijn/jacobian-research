# General-position closure for the squarefree quartic denominator

## Status

This note proves `HC4NHM11`.  It continues the squarefree synchronization
frontend `HC4NHM9` and the concurrence closure `HC4NHM10` by excluding all
sixteen tangent/transverse flag patterns when the four denominator lines are
in general position.

Normalize the four lines to the projective frame

\[
(L_1,L_2,L_3,L_4)=(x,y,z,x+y+z).
\tag{0.1}
\]

For each line, `HC4NHM9` supplies a constant direction (v_i) and the polar
condition

\[
D_{v_i}h_5\in L_i^{m_i}S_{4-m_i},
\qquad
m_i=\begin{cases}2,&v_i\text{ tangent to }L_i,\\3,&v_i\text{ transverse.}\end{cases}
\tag{0.2}
\]

The exact projective direction atlases and all rank minors used below are
replayed by

~~~bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_general_position_closure.py
~~~

This is a clean generic-corank-one exclusion.  Positive-defect and
lower-Smith packets remain separate.

## 1. Zero or one tangent flag

For a line (L), put (U_L(m)=L^mS_{4-m}\subset S_4).  Direct coefficient
rank gives

\[
\ker\left(\bigoplus_{i=1}^4U_{L_i}(3)\longrightarrow S_4\right)=0
\tag{1.1}
\]

and

\[
\ker\left(U_{L_1}(2)\oplus\bigoplus_{i=2}^4U_{L_i}(3)
\longrightarrow S_4\right)=0.
\tag{1.2}
\]

The projective automorphism group of a projective frame acts as the full
permutation group on the four lines, so (1.2) covers every placement of the
single tangent flag.

If (h_5) is not a cone, its polar map

\[
V\longrightarrow S_4,
\qquad v\longmapsto D_vh_5,
\tag{1.3}
\]

is injective.  Four directions in the three-dimensional space (V) have a
nonzero linear relation.  Applying (1.3) would give a nonzero syzygy in
(1.1) or (1.2), a contradiction.  Thus every solution is a cone and has
zero Hessian determinant.

## 2. Exactly two tangent flags

By frame symmetry take the tangent lines to be (x=0,y=0), and retain the
transverse line (z=0).  The fourth line is not needed.  We must classify

\[
D_{v_1}h_5\in(x^2),
\qquad
D_{v_2}h_5\in(y^2),
\qquad
D_{v_3}h_5\in(z^3).
\tag{2.1}
\]

The tangent directions have four projective charts.

### 2.1 Both finite tangent charts

Write

\[
v_1=(0,a,1),
\qquad
v_2=(b,0,1),
\qquad
v_3=(c,d,1).
\tag{2.2}
\]

The direction determinant, up to sign, is

\[
F=ab-ac-bd.
\tag{2.3}
\]

The (30\)-by-(21) coefficient matrix of (2.1) has a maximal minor

\[
\text{unit}\cdot a^7b^3d^3F.
\tag{2.4}
\]

On the coordinate branches of (2.4), further exact maximal minors are,
up to nonzero constants,

\[
\begin{array}{c|c}
a=0&b^7d^5\\
b=0&a^7c^4d\\
b=d=0&a^7c^5\\
d=0&a^7b^3c^3\,a(b-c).
\end{array}
\tag{2.5}
\]

Away from (F=0), the only remaining rank-drop component is

\[
c=d=0,
\tag{2.6}
\]

where the complete nullspace is

\[
\left\langle x^2y^3,x^3y^2\right\rangle.
\tag{2.7}
\]

It is binary, hence a cone.

### 2.2 One tangent direction at infinity

For

\[
v_1=(0,1,0),
\qquad v_2=(b,0,1),
\qquad v_3=(c,d,1),
\]

the direction determinant is (c-b).  The covering maximal minors are

\[
b^3d^3(b-c),\quad c^4d\ (b=0),\quad c^5\ (b=d=0),
\quad b^3c^3(b-c)\ (d=0).
\tag{2.8}
\]

Again the sole independent-direction rank drop is (c=d=0), with nullspace
(2.7).  The chart with the boundary on (y=0) is symmetric; its covering
minors are

\[
a^3d^3(a-d),\quad d^5\ (a=0),\quad a^4c^3\ (d=0),
\tag{2.9}
\]

and its only rank drop is the same binary space.

If both tangent directions are at infinity, the generic nullspace is
(\langle z^5\rangle).  At (c=d=0) it enlarges to

\[
\left\langle z^5,x^2y^3,x^3y^2\right\rangle.
\tag{2.10}
\]

Every form in (2.10) is a cone.

### 2.3 The dependent-direction divisor

It remains to justify that (F=0) cannot hide a non-cone.  The only quartic
syzygy among

\[
x^2S_2,\qquad y^2S_2,\qquad z^3S_1
\]

is

\[
(-x^2y^2,x^2y^2,0).
\tag{2.11}
\]

If the polar map were injective, a dependent triple of directions would
therefore force (v_1=v_2=\partial_z) and

\[
D_zh_5=\kappa x^2y^2.
\]

Hence (h_5=\kappa x^2y^2z+F_5(x,y)).  For the transverse direction
(v_3=(c,d,1)), the coefficient of (z) in (D_{v_3}h_5) is

\[
2\kappa(cxy^2+dx^2y).
\]

Divisibility by (z^3) makes (c=d=0), after which the constant term
(\kappa x^2y^2) is nonzero.  Thus (kappa=0), making (h_5) a cone.

This completes all six exactly-two-tangent patterns.

## 3. Three or four tangent flags

Any three lines in the general-position arrangement are nonconcurrent.
Normalize three tangent lines to (x=0,y=0,z=0).  Their directions have the
eight projective charts

\[
\begin{aligned}
v_1&=(0,a,1)\ \text{or}\ (0,1,0),\\
v_2&=(b,0,1)\ \text{or}\ (1,0,0),\\
v_3&=(c,1,0)\ \text{or}\ (1,0,0).
\end{aligned}
\tag{3.1}
\]

On the main chart the direction determinant is

\[
F=ac+b.
\tag{3.2}
\]

Three maximal minors cover (F\ne0):

\[
a^7b^3F,
\qquad b^7\ \text{on }a=0,
\qquad a^7c^4\ \text{on }b=0.
\tag{3.3}
\]

On (F=0), a parameter-independent (20\)-minor is nonzero and the complete
nullspace is

\[
\left\langle (ac z-cy+x)^5\right\rangle.
\tag{3.4}
\]

The seven boundary charts give exactly the same dichotomy.  Two have
constant nonzero direction determinant and no solution.  Two have direction
determinants (a) and (c); on their zero loci the nullspaces are

\[
\langle y^5\rangle,
\qquad
\langle(x-bz)^5\rangle.
\tag{3.5}
\]

The three identically dependent charts have nullspaces

\[
\langle(y-az)^5\rangle,
\qquad
\langle z^5\rangle,
\qquad
\langle z^5\rangle.
\tag{3.6}
\]

The checker verifies a nonzero (20\)-minor on every dependent chart,
including all parameter boundary points.  Thus three nonconcurrent tangent
polar conditions force a pure fifth power.  This closes all patterns with
three or four tangent flags.

## 4. Result

> **Theorem `HC4NHM11` -- General-position squarefree closure.**  In the
> clean generic-corank-one squarefree quartic-denominator packet of
> `HC4NHM9`, the line arrangement with no three concurrent has no nonzero
> Hessian determinant for any of its sixteen tangent/transverse flag
> patterns.  The exclusion occurs at the ternary polar-synchronization level,
> before the cleared Schur-gradient equations.

Together, `HC4NHM10` and `HC4NHM11` reduce the clean squarefree partition
from forty-eight flag families to exactly eight: the exactly-three-concurrent
arrangement in which the fourth, nonconcurrent line has a tangent kernel
flag.  The subsequent theorem `HC4NHM12` closes those last eight rows in
[`HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md`](HC4_SQUAREFREE_QUARTIC_TANGENT_FOURTH_CLOSURE.md).

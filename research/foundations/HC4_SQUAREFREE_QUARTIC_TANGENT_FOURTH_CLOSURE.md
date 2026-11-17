# Squarefree quartic tangent-fourth closure

## Status

This note proves `HC4NHM12`.  It closes the eight rows left by the
squarefree frontend `HC4NHM9`, the concurrence theorem `HC4NHM10`, and the
general-position theorem `HC4NHM11`.  The result concerns the clean
generic-corank-one quartic-denominator packet.  Positive-defect and
lower-Smith packets are separate.

The exact replay is

```bash
.venv/bin/python scripts/verify_hc4_squarefree_quartic_tangent_fourth_closure.py
```

It uses SymPy over `QQ` and Singular 4.4.1 for exact ideal saturation.

## 1. The last arrangement

Normalize the exactly-three-concurrent arrangement to

\[
 (L_1,L_2,L_3,L_4)=(x,y,x+y,z).
\]

The fourth flag is tangent.  Up to permutation of the three concurrent
lines, the eight patterns have four representatives

\[
 \mathrm{RRR;T},\qquad \mathrm{TRR;T},\qquad
 \mathrm{TTR;T},\qquad \mathrm{TTT;T}.
\]

For a tangent flag the corresponding quartic polar belongs to
\((L_i^2)_4\); for a transverse flag it belongs to \((L_i^3)_4\).

Assume that \(h\) is a noncone ternary quintic satisfying the four flag
conditions, and write

\[
 p_i=D_{v_i}h.
\]

Because \(h\) is not a cone, its polar map

\[
 K^3\longrightarrow K[x,y,z]_4,\qquad v\longmapsto D_vh
\]

is injective.  The four directions are dependent.  Choose a nonzero
relation \(\sum_i\lambda_i v_i=0\), and put

\[
 w_i=\lambda_i v_i,\qquad q_i=\lambda_i p_i=D_{w_i}h.
\]

Then

\[
 \sum_iw_i=0,\qquad \sum_iq_i=0,
 \qquad D_{w_i}q_j=D_{w_j}q_i.                 \tag{1}
\]

Injectivity also gives

\[
 \operatorname{rank}\{q_i\}=\operatorname{rank}\{w_i\}. \tag{2}
\]

Zero relation coefficients cause no problem: the zero columns occur on
both sides of (2).

## 2. Exact rank-at-least-two exclusion

For each representative, form the vector space

\[
 \mathcal S=\left\{(q_1,q_2,q_3,q_4):
   \sum q_i=0,\ q_i\in(L_i^{e_i})_4\right\},
 \qquad
 e_i=\begin{cases}2&T,\\3&R.\end{cases}
\]

Its dimensions in the four rows are respectively

\[
 1,\quad3,\quad6,\quad9.
\]

Use coordinates \(s\) on \(\mathcal S\), eliminate \(w_4\) using
\(\sum w_i=0\), and let \(I\) be the coefficient ideal of (1), together
with all tangent incidence equations \(L_i(w_i)=0\).  Let \(Q_2\) be the
ideal of the \(2\)-minors of the coefficient matrix with columns \(q_i\),
and let \(W_2\) be the ideal of the \(2\)-minors of the matrix with columns
\(w_i\).

The exact rational computations give, in every row,

\[
 \boxed{(I:Q_2^\infty):W_2^\infty=(1).}        \tag{3}
\]

Thus (1) has no solution for which both families have rank at least two.
By (2), every relation arising from a noncone quintic must have rank one.

The checker constructs every syzygy and every mixed-partial coefficient
from scratch.  It clears rational denominators and asks Singular for the
four double saturations; each standard basis is exactly `[1]`.  No
finite-field inference or bounded parameter search is used.

## 3. Rank-one relations

By the repeated-direction lemma in `HC4NHM9`, a rank-one relation can be
supported only on a tangent--tangent pair.  Three active repeated
directions already make \(h\) a cone.  For two active columns the same
lemma gives

\[
 D_vh=\kappa L_i^2L_j^2,\qquad \kappa\ne0.     \tag{4}
\]

There are two normal forms.

### 3.1 Two lines in the concurrent pencil

Normalize the pair to \((x,y)\).  Its common tangent is \(v=\partial_z\),
so (4) integrates to

\[
 h=\kappa x^2y^2z+F_5(x,y).
\]

If the fourth line \(z\) is tangent with direction \((a,b,0)\), its polar
can be divisible by \(z^2\) only if the coefficient of \(z\) vanishes:

\[
 D_{(a,b)}(x^2y^2)=2xy(ay+bx)=0.
\]

Hence \(a=b=0\), impossible.  In the `TTR;T` row the unused pencil line
is instead transverse.  For \(x+y\), divisibility of the displayed binary
cubic by \((x+y)^3\) again forces \(a=b=0\), contradicting transversality.

### 3.2 One pencil line and the fourth line

Normalize the pair to \((x,z)\).  Its common tangent is
\(v=\partial_y\), hence

\[
 h=\kappa x^2z^2y+F_5(x,z).
\]

For a direction \((a,b,c)\) on the unused line \(y\), the coefficient of
\(y\) in its polar is

\[
 2\kappa xz(az+cx).
\]

Divisibility by \(y^2\), and a fortiori by \(y^3\), forces \(a=c=0\).
If the flag is tangent the direction is zero.  If it is transverse then
\(b\ne0\), but the constant term is
\(b\kappa x^2z^2\ne0\), again impossible.  Permuting the pencil lines
covers every remaining tangent pair.

Therefore rank one is also empty in all four representatives.

## 4. The theorem

> **Theorem `HC4NHM12` -- Tangent-fourth squarefree closure.**  In the
> clean generic-corank-one quintic Hessian--Schur packet with squarefree
> quartic minimal denominator, every exactly-three-concurrent flag pattern
> whose fourth line is tangent is empty.  Consequently all forty-eight
> split-squarefree line-arrangement/flag rows are empty.

Together with `HC4NHM7`, `HC4NHM8`, `HC4NHM10`, and `HC4NHM11`, this closes
every clean quartic-denominator packet whose denominator splits completely
into linear factors.  It does not close clean packets with irreducible
conic, cubic, or quartic denominator components, positive-defect components,
generic lower-Smith components, or the separate leading-Hessian rank-one/
rank-two synchronization frontiers.

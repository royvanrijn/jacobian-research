# A repeatable cubic index-four block with exact inverse reset

Let K be any field of characteristic zero. The following construction can be
iterated indefinitely while its inverse degrees remain bounded. It supplies
a direct control against inferring inverse growth from an indefinitely
extendible Engel algebra, growing coordinate weights, or nonzero products.

## The repeatable block

Write J=JH for the Jacobian of a homogeneous cubic map H on a polynomial
coordinate algebra K[x]. Assume that u,w,p are distinct coordinates and
\[
 J^4=0,\qquad H_u=0,\qquad dw\,J^3=0,\qquad H_p=u^2w .
\]
Here differentials are row vectors. Let V be a nonzero homogeneous quadratic
polynomial in the old coordinates with
\[
 dV\,J^3=0.
\]
Introduce three independent new coordinates z,y,t, retain every old
component, and set
\[
 H_z=uV,\qquad H_y=wV,\qquad
 H_t=uwz-u^2y+pV.                                      \tag{R1}
\]
Define the returned quadratic
\[
 V^+=pz-ut.                                            \tag{R2}
\]

**Theorem RepeatableIndexFourResetBlock.**
The enlarged map satisfies J_new^4=0 and
dV^+ J_new^3=0. The coordinates u,w,p keep the same hypotheses, so (R1)
can be repeated with V replaced by V^+. If the original Jacobian has
nilpotency index exactly four, every enlarged Jacobian also has index
exactly four.

For a direct proof put
\[
 s=dwJ,\quad r=dwJ^2,\quad
 a=dVJ,\quad b=dVJ^2,\quad c=dVJ^3=0 .
\]
We have du J=0 and
\[
 dpJ=2uw\,du+u^2dw,\quad dpJ^2=u^2s,\quad dpJ^3=u^2r.
\]
Ordinary block multiplication gives
\[
\begin{array}{c|c|c}
 &\hbox{third row}&\hbox{fourth row}\\ \hline
 z&ub&uc\\
 y&Vr+wb&wc\\
 t&uzr+pb&pc .
\end{array}                                             \tag{R3}
\]
For the y fourth row one also uses rJ=dwJ^3=0.
The old fourth rows are zero, proving the first assertion. Moreover
\[
 dV^+J_{\rm new}^3
 =z(u^2r)+p(ub)-u(uzr+pb)=0.                           \tag{R4}
\]
A stronger returned-row identity, obtained by one fewer multiplication, is
\[
 dV^+J_{\rm new}^2=-3u^2wV\,du.                         \tag{R4a}
\]
Thus its third row is zero because du J_new=0; this calculation does not
make its second row zero.
All old rows and their products with old columns are preserved. In
particular an old nonzero third row remains nonzero, and dw J_new^3=0.
The returned quadratic is nonzero because z,t are independent coordinates.

The identities are polynomial identities; no localization or division by
u,w,V occurs.

## The forced coefficient

In (R1), replace the coefficient of pV in H_t by a scalar lambda. Under
the same hypotheses the new last fourth row becomes exactly
\[
 (\lambda-1)u^2V r.                                    \tag{R5}
\]
Thus if r is nonzero, the coefficient ring is a domain and V is nonzero,
preservation of fourth-power zero is equivalent to lambda=1. This
particular attempted change of coefficient cannot evade the reset below.

For completeness, one wider identity explains both cancellations. Set
\[
 H_z=\alpha uV,\quad H_y=\beta wV,\quad
 H_t=\gamma uwz+\delta u^2y+\epsilon pV .
\]
Its new last fourth row is
\[
 u^2V(\delta\beta+\epsilon)r+
 u^2w(\gamma\alpha+\delta\beta)b+\epsilon pc.           \tag{R6}
\]
Consequently the sufficient coefficient equations
epsilon=gamma alpha=-delta beta return the invariant
gamma pz-ut. Its inverse reset is
gamma p_input Z-uT-delta u^3Y.
Equation (R6) does not assert necessity without an additional independence
hypothesis on its obstruction rows.

## Exact inverse reset and a uniform degree ceiling

Now assume in addition that F=Id-H is a polynomial automorphism. Denote
its inverse on the old input by G. Put
\[
 U=G_u=u_{\rm input},\quad W=G_w,\quad
 P=G_p=p_{\rm input}+U^2W,\quad v=V(G).
\]
For new target coordinates Z,Y,T, triangular substitution in (R1) gives
\[
 \begin{aligned}
 z_{\rm inv}&=Z+Uv,\\
 y_{\rm inv}&=Y+Wv,\\
 t_{\rm inv}&=T+UWZ-U^2Y+Pv .
 \end{aligned}                                        \tag{R7}
\]
In particular the enlarged map is again a polynomial automorphism, and
\[
 V^+(G_{\rm new})
   =Pz_{\rm inv}-Ut_{\rm inv}
   =p_{\rm input}Z-UT+U^3Y.                            \tag{R8}
\]
The right side is independent of both v and W. Its ordinary total degree
is exactly four. This is a reset, rather than a degree recurrence.

More precisely let
\[
 D=\deg G,\qquad B=\deg W\ge1,\qquad E=\deg V(G)\ge2 .
\]
The lower bound on E follows because G=x plus terms of degree at least
three and V is a nonzero homogeneous quadratic. Also deg P=B+2. In the
first block, the term Pv in t_inv has degree B+E+2 and contains no
new target variables; the other terms have degree at most B+2 or three.
Thus the first enlarged inverse has degree exactly max(D,B+E+2).

In every later block, its inherited v is the polynomial (R8) from the
preceding block and has degree four. Its top term is U^3 times that
block's independent Y coordinate. Hence Pv has degree B+6, and cannot
cancel with UWZ, U^2Y or T from the new block. All newly added coordinates
have degree at most B+6, with equality at t_inv. Therefore every iterate
has degree at most
\[
 \max\{D,B+E+2,B+6\},                                  \tag{R9}
\]
and after at least two blocks the full inverse degree equals this maximum.
No bound for arbitrary cubic index-four maps is asserted.

## An explicit infinite family whose inverse degree stays 207

Use the literal 58-variable map in
[the literal degree-207 map](CUBIC_INDEX_FOUR_SYZYGY58.md), with
\[
 u=x_0,\quad w=x_{47},\quad p=p_1=x_{48},\qquad
 V_0=uA-4uE+9p_1p_2 .
\]
Its stated positive weights give weights 1,29,31 to u,w,p, and weight90
to V_0. Direct multiplication of its literal Jacobian verifies
dw J^3=0, dV_0 J^3=0, and r=dw J^2 is nonzero. The full Jacobian has
fourth power zero and nonzero third power.

The old inverse has degree D=207, with deg G_w=29. Weighted homogeneity
gives deg V_0(G)<=90. Thus the first added block has inverse degree at
most121, and every later added block has inverse degree at most35.
The old degree207 coordinate remains present. We obtain, in dimensions
\[
 58+3k\quad(k\ge0),
\]
homogeneous cubic maps with Jacobian nilpotency index exactly four and
inverse degree exactly207.

The construction has compatible positive weights. Its kth returned
quadratic has weight90+32k, and the last coordinate in its kth block has
weight
\[
 121+32(k-1)\qquad(k\ge1).
\]
These weights tend to infinity. Formula (R8) nevertheless gives the
degree-four inverse value of each returned quadratic. In particular
unbounded weights in this family are not unbounded inverse degrees.

## Scope and classical context

The construction and reset are proved directly by block multiplication.
Umirbaev, Polarization algebras and their relations, MPIM14-60,
Lemmas7--8 and Corollary2, printed pages7--8, explain the polarization
algebra setting. The four-dimensional example on printed pages9--10
distinguishes the full algebra's nonnilpotence from the dimension of its
free one-generated algebra. That distinction is relevant context here;
no infinite inverse-degree conclusion is imported from it, and no
novelty claim about general Engel constructions is made.

The checker proves the displayed symbolic coefficient identities, then
reconstructs the entire literal 58-variable map and five iterates using
rational arithmetic. It checks the hypotheses, nonzero third powers,
zero fourth powers, each returned quadratic's annihilation, the growing
positive weights, two coefficient-perturbation failures, and exact
inverse resets at full-coordinate test inputs. The infinite statements
follow from (R3)--(R9), not from the finite number of iterates.

Run from a checkout, or specify its root explicitly:

    python3 research/jacobian/subtree-shuffling/verify_repeatable_index_four_reset.py \
      --repo-root . \
      --output /tmp/repeatable_index_four_reset_review.json

The input and mandatory source ledger accompany this proof. All
computations are internal exact checks; no external or formal review
is claimed.

Updated 20 November 2026.

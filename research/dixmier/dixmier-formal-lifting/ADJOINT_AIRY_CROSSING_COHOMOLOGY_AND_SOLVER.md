# Adjoint Airy cohomology across the crossing XQ=0

Let \(k\) have characteristic zero, put
\[
 M(t)=\begin{pmatrix}0&1\\t&0\end{pmatrix},\qquad
 V=\mathfrak{sl}_2(k[t]),\qquad D=\partial_t+\operatorname{ad}M(t).
\]
The one-variable operator \(D:V\to V\) is a polynomial isomorphism.
Its pullback along \(t=XQ\) nevertheless has three polynomial
degree-two cohomology classes. Thus one-variable invertibility does
not supply a primitive for every closed pulled-back matrix curvature.

## Exact one-variable solver

Write \(C=\left(\begin{smallmatrix}a&b\\c&-a\end{smallmatrix}\right)\).
For prescribed
\(R=\left(\begin{smallmatrix}u&v\\w&-u\end{smallmatrix}\right)\), the equation
\(DC=R\) is
\[
 a'+c-tb=u,\qquad b'-2a=v,\qquad c'+2ta=w.
\]
Eliminate \(a,c\):
\[
 Lb=w-u'-v''/2+tv,\qquad
 L=-\partial_t^3/2+2t\partial_t+1,                        \tag{1}
\]
\[
                 a=(b'-v)/2,\qquad c=u-a'+tb.            \tag{2}
\]
On \(t^n\), \(L\) has diagonal coefficient \(2n+1\) and only a
lower term of degree \(n-3\). It is therefore invertible on every
finite-degree polynomial space. Explicitly, for \(Lb=r\),
\[
 b_n=\frac{r_n+(n+3)(n+2)(n+1)b_{n+3}/2}{2n+1},          \tag{3}
\]
computed downward from \(\deg r\), with higher coefficients zero.
This is a finite exact solver. Equations (1)--(2) prove
\[
                        \deg_t D^{-1}R\leq\deg_tR+2.     \tag{4}
\]
In particular \(D\) has no nonzero polynomial horizontal traceless
endomorphism.

## The polynomial crossing complex

On \(k[X,Q]\), use the flat covariant differential
\[
 \nabla=d+\operatorname{ad}M(t)\,dt,\qquad t=XQ.
\]
Its coefficient derivations are
\[
 \nabla_X=\partial_X+Q\operatorname{ad}M(t),\qquad
 \nabla_Q=\partial_Q+X\operatorname{ad}M(t).
\]
The following cohomology is for polynomial forms with traceless
matrix coefficients:
\[
                 H^0_\nabla=H^1_\nabla=0,\qquad
                 H^2_\nabla\simeq\mathfrak{sl}_2(k).      \tag{5}
\]

Give \(X,dX\) weight \(1\) and \(Q,dQ\) weight \(-1\).
The connection preserves these integer weights. For
\(E=X\partial_X-Q\partial_Q\), \(dt(E)=0\), so the covariant Cartan
identity is
\[
                 \nabla\iota_E+\iota_E\nabla=\mathcal L_E.
\]
Every nonzero weight \(m\) is therefore contracted by \(\iota_E/m\).
Only weight zero contributes cohomology.

The zero-weight complex is exactly
\[
 V\longrightarrow V\oplus V\longrightarrow V,\qquad
 C\longmapsto(DC,DC),\qquad
 (A,B)\longmapsto(tD+1)(B-A),                             \tag{6}
\]
where a one-form is \(QA(t)dX+XB(t)dQ\) and a two-form is
\(R(t)dX\wedge dQ\). Since
\[
                              tD+1=D\circ t,
\]
\(D\) is bijective and \(D\circ t\) is injective. This proves
\(H^0=H^1=0\), and
\[
 H^2=V/D(tV)\xrightarrow[\sim]{D^{-1}}
                   V/tV=\mathfrak{sl}_2(k).              \tag{7}
\]
In particular the exact obstruction functional on an invariant
two-form is
\[
                    R(t)\longmapsto(D^{-1}R)(0).         \tag{8}
\]

For an arbitrary two-form, first keep its zero-weight coefficient
\(R_0(t)\). Its three obstruction coordinates are
\((D^{-1}R_0)(0)\). All other weights are automatically exact.

The classes may be represented by
\[
 D(E_{12})=-tH,\qquad D(E_{21})=H,\qquad
 D(H)=-2E_{12}+2tE_{21},\quad H=\operatorname{diag}(1,-1).
\]
Their obstruction vectors are the three constant basis matrices.
For example the closed constant curvature
\[
                              H\,dX\wedge dQ
\]
has no polynomial covariant primitive, since \(D^{-1}H=E_{21}\).

For full \(2\)-by-\(2\) matrices the trace summand is the ordinary
polynomial de Rham complex. Thus its \(H^0\) is \(kI_2\), \(H^1=0\),
and \(H^2\) is still the three-dimensional traceless group above.
These statements concern polynomial cohomology, not invariant
subspaces over a Weyl division ring.

## Polynomial primitives with a degree bound

Measure the ordinary degree of form coefficients, without assigning
degree to \(dX,dQ\). A closed polynomial traceless one-form of coefficient
degree at most \(d\) has a polynomial primitive of degree at most \(d+3\).
An arbitrary polynomial traceless two-form of degree at most \(d\)
has such a one-form primitive if and only if (8) vanishes; its
primitive can also be chosen with coefficient degree at most \(d+3\).

For a nonzero weight, the Cartan homotopy multiplies coefficients by
\(X\) or \(Q\) and raises degree by at most one.
A closed zero-weight one-form in (6) has \(A=B\); its primitive is
\(D^{-1}A\). Since its original coefficients contain an extra \(X\)
or \(Q\), (4) raises the ordinary position degree by at most three.

For a compatible zero-weight two-form, put \(U=D^{-1}R\).
Vanishing of (8) says exactly that \(U=tV\) for a polynomial \(V\).
Then
\[
                              XV(t)\,dQ                 \tag{9}
\]
is a primitive, because its differential is
\((tD+1)V\,dX\wedge dQ=R\,dX\wedge dQ\).
Using (4) and division by \(t\) gives the same bound \(d+3\).

As a concrete positive control,
\(R=E_{21}+tH=D(tE_{21})\) has primitive \(XE_{21}\,dQ\).
The constant \(H\) control above supplies a genuine negative
compatibility test.

## Adding the remaining affine source variables

Extend to the full source \(k[X,Q,P,Z]\), with the same connection
depending only on \(X,Q\). The polynomial de Rham complex in \(P,Z\)
contracts by their Euler field, and its contraction commutes with
the position connection. More explicitly, grade by total \(P,Z\)
degree plus the number of \(dP,dZ\) factors. Every positive grade is
contracted by the Euler contraction divided by that grade; grade
zero is the position complex (6). Thus the full source has the same
cohomology (5), and the same \(d+3\) primitive bound for closed
one-forms and compatible closed two-forms. The positive-grade
contraction raises coefficient degree by at most one.

To compute the full-source obstruction, restrict the closed two-form
to \(P=Z=0\) and discard \(dP,dZ\), extract its zero \(X,Q\)-weight
coefficient \(R_0(t)dX\wedge dQ\), and apply (8).

## The exact covariant linear solver for Long cochains

Let \(f=(R,T,D,S)\) be the actual canonical full polynomial Long frame,
with constant invertible Poisson matrix \(J_{ij}=\{f_i,f_j\}\), and put
\[
 \delta_i=\{f_i,-\},\qquad
 A_i=M(t)\{f_i,t\},\qquad \mathcal D_i=\delta_i+\operatorname{ad}A_i.
\]
The \(\delta_i\) form a global polynomial frame of commuting derivations.
Its polynomial dual forms are
\[
                    \theta^j=\sum_k(J^{-1})_{kj}\,df_k. \tag{10}
\]
The flat connection evaluated on \(\delta_i\) is \(A_i\), so the
cochain differential for the commuting \(\mathcal D_i\) identifies
with the covariant de Rham differential above.

In a matrix Moyal jet
\(F_i=f_iI_2+hA_i+\cdots\), changing the coefficient at \(h^m\),
\(m\geq2\), by \(C_i\) changes the commutator at \(h^{m+1}\) by
\[
                       \mathcal D_iC_j-\mathcal D_jC_i. \tag{11}
\]
The \(A_i\) are precisely the first coefficients of the admitted
polynomial-connection source embedding; they are not scalar
constant-sector shifts.

For a closed traceless defect cochain \(E_{ij}\), form
\[
                        \Omega=\sum_{i<j}E_{ij}\theta^i\wedge\theta^j.
\]
A polynomial correction exists exactly when its three crossing
coordinates vanish. If the maximum classical image degree is \(d_f\)
and the defect degree is at most \(d_E\), equations (9)--(10) give a
correction tuple of degree at most
\[
                             d_E+3d_f.                  \tag{12}
\]
Indeed two dual forms add at most \(2(d_f-1)\), the primitive solver
adds three, and evaluation on a frame derivation adds at most
\(d_f-1\). For the canonical Long frame \(d_f=23\), (12) is
\(d_E+69\). This is a finite linear algebraic correction bound at
one order, not a uniform all-orders spatial bound.

Closedness alone is insufficient: the cochain obtained by evaluating
\(H\,dX\wedge dQ\) on the \(\delta_i\) is closed and has nonzero
crossing obstruction. It is an explicit algebraic countercontrol,
not a claim that this defect arises from an actual lower quantum jet.

Although \(H^1=0\) gives covariant primitives for infinitesimal closed
corrections, an iteration of such primitives may have unbounded
degree. No bounded formal-gauge equivalence, finite matrix Long lift,
or scalar realization follows. Arbitrary noncompositional finite
matrix/graph target actions remain open.

The [exact checker](check_adjoint_airy_crossing.py) verifies the finite
descending inverse, the three obstruction representatives, closed
one-form and compatible two-form solvers, and the constant-curvature
negative control. Its [receipt](adjoint_airy_crossing_controls.json)
pins this note and checker. The all-polynomial cohomology, degree
bounds and Long comparison are the written arguments above, not an
extrapolation from finite controls. Formal verification, independent
whole-proof replay and external review are unclaimed.

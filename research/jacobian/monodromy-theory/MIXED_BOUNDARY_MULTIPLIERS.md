# Mixed multipliers cannot contain the critical boundary

A scalar factor supported on the critical hyperplane is impossible for
every primitive quadratic rank4 full-rank normalized input, even when
other factors accompany it. The logarithmic boundary correction is
retained throughout. Its polynomial residues, rather than its omission,
give the obstruction.

## 1. A plane-image boundary lemma with the mixed term

Work over algebraically closed characteristic-zero constants. Let
e be a positive integer, Q:k^3->k^3 homogeneous quadratic with image
in a constant plane P, ell a constant row, and J a nonzero homogeneous
polynomial. Define the rational matrix
\[
 B=C_e+\binom{\ell x}{Q}(0,d\log J),\qquad
 C_e=\begin{pmatrix}(e+1)\ell x&\ell\\eQ&DQ\end{pmatrix}.       \tag{1}
\]

**KMFI149.** If B is nilpotent identically, then ell annihilates P.
Consequently if the linear image span of Q has dimension at most1,
nilpotency forces ell=0.

Suppose instead ell is nonzero on P. By simultaneous constant
source/target coordinates put P={z=0} and the first coordinate x=ell*x.
The last output row of B is zero. Write J=z^b*J1 with z not dividing J1.
The first three rows and columns depend only on dlogJ in x,y,
so they are regular at the generic point of z0.
Restrict them there, setting j=J1(x,y,0), a nonzero binary homogeneous
polynomial of degree n>=0, and R=(R1,R2)=Q_high(x,y,0). The resulting
nilpotent matrix is
\[
 M=\begin{pmatrix}(e+1)x&1&0\\eR_1&R_{1x}&R_{1y}\\
 eR_2&R_{2x}&R_{2y}\end{pmatrix}
 +\begin{pmatrix}x\\R_1\\R_2\end{pmatrix}(0,j_x/j,j_y/j).      \tag{2}
\]
This step keeps every relevant mixed term; the possibly singular z
column does not enter this principal block, whose nilpotency follows
from the zero last row.

For the first matrix C in(2), Euler gives
C*(1,x,y)^t=(e+2)*(x,R1,R2)^t. The determinant lemma, also valid after
clearing denominators when C is singular, gives
\[
 \det M=\frac{e+2+xj_x/j+yj_y/j}{2}\,x\,\det D(R_1,R_2)
       =\frac{e+n+2}{2}\,x\,\det D(R_1,R_2).                \tag{3}
\]
The coefficient cannot vanish. Thus the binary Jacobian is zero,
and the two homogeneous quadratics are proportional, as in KMFI147.

First suppose R1 is nonzero. A simultaneous y->y-kx change removes
R2 and preserves x. The remaining nilpotent two-dimensional block,
with R=R1, is
\[
 \begin{pmatrix}
 (e+1)x&1+xj_x/j\\ eR&R_x+Rj_x/j
 \end{pmatrix}.
\]
Trace and determinant give
\[
 Rj_x/j=-(e+1)x-R_x,\qquad
                  e(xR_x-R)=(e+1)x^2.                   \tag{4}
\]
Write R=a*x^2+b*x*y+c*y^2. Then a=(e+1)/e and c=0.
If b!=0, the forced rational form
\[
 j_x/j=-\frac{(e+1+2a)x+by}{x(ax+by)}
\]
has residue -1 at x0, impossible for a polynomial logarithmic
derivative. If b=0, it is -(e+2)/x, again a negative residue.

If R1=0, R2=0 contradicts the trace immediately. Otherwise deleting
the zero R1 row leaves the block
\[
 \begin{pmatrix}
 (e+1)x&xj_y/j\\eR&R_y+Rj_y/j
 \end{pmatrix},\quad R=R2.
\]
Its trace and determinant give
\[
 R_y=(e+1)x/e,\qquad Rj_y/j=-(e+1)x-R_y.
\]
Thus R=x(ax+by) with b=(e+1)/e, and the residue at ax+by0
of the forced y-logarithmic derivative is -(e+1), impossible.
All cases contradict nilpotency, proving the assertion.

If image-span(Q) has dimension<=1 and ell!=0, choose a plane containing
that span on which ell is nonzero. The already proved assertion gives
a contradiction. This proves the final conclusion.

Polynomiality and positivity matter. With ell=(1,0,0),
Q=((e+1)x^2/e,0,0) and the rational J=x^(-(e+2)), matrix(1)
is nonzero and square-zero for every positive e. These rational controls
are not polynomial multipliers.

## 2. Boundary factors force matched charts in all four families

Let V be primitive homogeneous quadratic in rank4 with
detDV=delta*u^4 nonzero. Let the homogeneous scalar F=u^e*G, e>0, u not dividing G,
give a normalized polynomial input. No full-rank assumption is needed
in the following chart-matching conclusion.

**KMFI150.** Normalization forces the source and target boundary
hyperplanes to match. Consequently, under generic rankDX4, no such
F exists. Every full-rank multiplier for a primitive quadratic rank4
direction must be coprime to u.

By KMFI146, V_u=uL. If the charts differ, the high boundary
Q=V_high|u0 has image in a plane. Its exact nilpotent boundary matrix
is(1), with J=G(0,x,y,z), and L=a*u+ell*x.
Use the four exhaustive forms V=A*K of KMFI143.

For K_C,K_D the image span of Q has dimension1.
KMFI149 forces ell=0, so V_u=a*u^2, with a nonzero since detDV!=0.
This is a matched component, contradicting different charts.

For K_A,K_B, KMFI146 makes the first output row (a,b,0,0).
If b=0 the charts already match. Suppose b!=0.
Here ell=(b,0,0), and the boundary image plane has dimension2.
KMFI149 forces ell*Q=0, hence A12=A13=0.
Thus the first two rows of A form the block
\[
                    \begin{pmatrix}a&b&0&0\\c&d&0&0\end{pmatrix}. \tag{5}
\]
In particular the x output of Q is zero, and all of Q is independent
of z. Write its y component as
\[
 Q_y=h x^2+k y^2\quad(K_A),\qquad
 Q_y=h x^2+kxy\quad(K_B).                                  \tag{6}
\]
Here subscripts on Q denote output coordinates, not derivatives.

Take J's highest coefficient in z, denoted j(x,y), nonzero.
As z tends to infinity, its x,y logarithmic derivatives tend to
j_x/j,j_y/j, and its z logarithmic derivative tends to zero.
All other entries of(1) are independent of z. The rational limit
therefore exists and remains nilpotent.
Its x row and z column are zero, leaving the nilpotent block
\[
 \begin{pmatrix}
 (e+1)b x&b x\,j_y/j\\ eR&R_y+Rj_y/j
 \end{pmatrix},\qquad R=Q_y.                              \tag{7}
\]
Its trace and determinant force
\[
                 eR_y=(e+1)b x.                          \tag{8}
\]
For K_A, the left side is 2eky, contradicting b!=0.
For K_B, (8) gives k=(e+1)b/e, so k!=0.
The trace now requires
\[
       j_y/j=-\frac{(e+1)b+k}{hx+ky}.
\]
Its residue is -((e+1)b+k)/k=-(e+1), again impossible for polynomial j.
This excludes both unmatched chain cases.

All four forms therefore match whenever a positive boundary power
occurs. Under full rank the primitive matched class is excluded by
KMFI141, proving the asserted coprimality. No multiplier degree,
factor number, other-factor degree or multiplicity bound was used.

The full-rank hypothesis in this last exclusion is indispensable:
the normalized primitive quadratic rank3 family of KMFI148 has
F=u^e for every positive e and whole sweep Jacobian gamma^4.
It is compatible with, and shows the sharp boundary of, chart matching.

## 3. The remaining input problem and the full action objective

[Helper](mixed_boundary.py) · [checker](verify_mixed_boundary.py) ·
[model](mixed_boundary_model.json) · [receipt](mixed_boundary_receipt.json).

The proof removes every u factor, including inside mixed multipliers.
Together with KMFI145 and KMFI148, a remaining primitive quadratic
rank4 full-rank candidate must have
\[
 F=c\prod_i D_i^{m_i},\qquad m_i>0,\quad
 1\leq\deg D_i\leq3,\quad D_i\not\sim u,\quad F\notin k.
\]
The unmatched target hyperplane alpha is also forbidden by KMFI144.
Each factor still satisfies D_i(V)=c_i*u^(degD_i)*D_i.
Hence the whole scalar obeys F(V)=c'*u^(degF)*F.
Multiplicities and total degree remain unbounded; these necessary
conditions are not sufficient for normalization or full rank.
Arbitrary matrices with no boundary factor, higher primitive degrees,
nonhomogeneous fields, exact generic groups and absolute exceptional11
remain open. No new Keller occurrence group is supplied.

The previous goal turn was progress: KMFI146--148 established exact
mixed boundary matrices and the pure-power obstruction. This turn
extends the obstruction to every scalar containing u by retaining
and testing the mixed correction. Exact whole identities supplement
the unbounded written argument; no formal verification, independent
whole-proof replay or external review is claimed. The full requested
faithful-transitive action theory remains the active goal.

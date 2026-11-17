# The complex graded Long first-jet family at cap 24

At common spatial degree 24, the complete odd first-jet kernel has
dimension eight. Two directions admit uniformly bounded normalizations.
The full next equation excludes the new source \(X^2Z\) direction,
while retaining a coupled first-image-changing branch over
\(\mathbb Q(\sqrt{-2})\). This note classifies the field-valued jets
through commutator order \(h^3\), including their entire second
coefficient families. It proves no higher compatibility or full lift.

## Complete first-coefficient space

Use source order \((X,Q,P,Z)\), brackets
\(\{P,X\}=\{Z,Q\}=1\), weights \((-1,1,4,2)\), and \(h\)-weight 3.
Let \(f=(R,T,D,S)\) be the full canonical Long polynomial symbols from
[the exact input transcription](../dixmier-audit/verify_long_pbw.py).
Their weights are \(w=(-1,2,4,1)\), and they satisfy all six polynomial
Poisson frame equations. Work with Weyl symbols and the Moyal product.

Write
\[
                  F=f+h g^{(1)}+h^2g^{(2)}+O(h^3).
\]
The coefficient of \(h^2\) in the required commutators gives
\[
 \{g^{(1)}_i,f_j\}+\{f_i,g^{(1)}_j\}=0,\qquad i<j.             \tag{1}
\]
Every coefficient has ordinary source degree at most 24, and
\(g^{(1)}_i\) has weight \(w_i-3\). Its ENTIRE monomial basis is
\[
 X^aQ^bP^cZ^d,\qquad
 a=b+4c+2d+3-w_i\ge0,\quad a+b+c+d\le24.                      \tag{2}
\]
The component counts are \((91,127,155,114)\), totaling 487.

Put
\[
\begin{split}
 B&=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\qquad t=XQ,\qquad Y=RA_0-3t,\\
 H_4&=Y^4+4tY^2-\frac{16}{3}tY+10t^2,\\
 N&=9Y^5+45tY^3+45t^2Y-60tY^2+110tY-720t^2,\\
 H_6&=R^4D+N/25,\\
 M&=X^2Z,\\
 H_5&=9Y^5+45tY^3-60tY^2+135t^2Y-10tY-555t^2.
\end{split}                                                    \tag{3}
\]
These are commuting classical polynomial expressions, used as
Hamiltonians. The complete solution space of (1)--(2) is
\[
 g^{(1)}=\{G,f\},\qquad
 G\in\operatorname{span}_K
       \{t,Y,Y^2,R^2T,H_4,H_6,M,H_5\}.                         \tag{4}
\]
This first-kernel statement holds over every characteristic-zero field.

| Hamiltonian | Degrees of its \(R,T,D,S\) variations |
|---|---|
| \(t\) | \((-,15,23,11)\) |
| \(Y\) | \((-,8,16,4)\) |
| \(Y^2\) | \((-,15,23,11)\) |
| \(R^2T\) | \((-,-,18,6)\) |
| \(H_4\) | \((-,17,23,13)\) |
| \(H_6\) | \((12,17,23,20)\) |
| \(M\) | \((4,16,24,12)\) |
| \(H_5\) | \((-,24,23,20)\) |

A dash denotes zero. In particular,
\[
              \{H_6,R\}=R^4,\qquad \{M,R\}=-3X^4,
              \qquad \{H_5,R\}=0.                             \tag{5}
\]

The full matrix of (1) has 2666 support rows and 77128 nonzero
entries. The [witness](graded_long_cap24_complex_witness.json) specifies
a 479-square minor with determinant
\[
                         763154\pmod {1000003}.                \tag{6}
\]
All denominators are invertible at this prime. The eight exact
Hamiltonian variations in (4) are admitted kernel vectors, and a
specified eight-square coefficient minor has rational determinant
\(3440640\). Thus the matrix has rational rank exactly 479 and
kernel dimension eight. This proves completeness of (4).

## Two bounded normalizations and the full nonlinear next equation

There is the classical identity
\[
                              RS=-6t-3Y^2-8Y.
\]
The source \(t=XQ\) shear
\[
 P\mapsto P-bhQ,\qquad Z\mapsto Z-bhX,\qquad X,Q\ {\rm fixed}
\]
is linear in source variables. The target \(RS\) shear
\[
 T\mapsto T+bhR,\qquad D\mapsto D-bhS,\qquad R,S\ {\rm fixed}
\]
is linear in the four images. Both preserve any common spatial cap
and the grading. After only these two normalizations, the first
Hamiltonian has the complete form
\[
             G=aY+bR^2T+cH_4+dH_6+eM+fH_5.                   \tag{7}
\]
The new \(M=X^2Z\) parameter is retained in this equation. Its natural
source shear \(Q\mapsto Q+ehX^2,\ P\mapsto P-2ehXZ\) is nonlinear
in source variables; it is not used as a common-cap normalization.

Allow EVERY second coefficient of degree at most 24 and weight
\(w_i-6\). The complete component counts are
\[
                         (62,91,114,80),\qquad 347.
\]
Let \(\mathcal A\) be its linearized Poisson map
\[
            (\mathcal A U)_{ij}=\{U_i,f_j\}+\{f_i,U_j\}.
\]
The entire commutator-order \(h^3\) equation is
\[
 (\mathcal A g^{(2)})_{ij}
       =-\{g^{(1)}_i,g^{(1)}_j\}-\Pi^3(f_i,f_j)/24.            \tag{8}
\]
It includes the quadratic first-jet term. The right side has one
constant and 21 quadratic parameter coefficient vectors. Their union
with the coefficient columns has 2137 support rows.

**Field classification through \(h^3\).** Over every
characteristic-zero field \(K\), equation (8) has an admitted solution
if and only if
\[
                  e=0,\qquad d^2-250df+35625f^2=0.             \tag{9}
\]
The parameters \(a,b,c\) are unrestricted at this order.

For necessity, the witness retains three rational functionals
\(\ell_1,\ell_2,\ell_3\), supported respectively on 327, 63 and 62
specified rows. They annihilate all 347 columns of \(\mathcal A\)
and evaluate on the ENTIRE right side of (8) as
\[
                    d^2-250df+35625f^2,\qquad e^2,\qquad ef.   \tag{10}
\]
The replay checks each of these rational identities. Since \(K\) is
a field, \(e^2=0\) implies \(e=0\). No inference about nilpotent
coefficients over nonreduced base rings is made.

The quadratic condition in (9) is equivalently
\[
                         (d-125f)^2=-20000f^2.                 \tag{11}
\]
If \(-2\) is not a square in \(K\), it forces \(d=f=0\).
If \(K\) contains a root \(q\) with \(q^2=-2\), its two possibilities
are
\[
                         d=(125\pm100q)f.                      \tag{12}
\]
Neither conjugate branch is discarded. For sufficiency and the full
second-coefficient family, it remains to give exact points on these
branches.

## Exact polynomial families over \(\mathbb E=\mathbb Q(q)\)

Fix \(q^2=-2\), put \(z=f\), and use the branch
\[
              G=aY+bR^2T+cH_4+z\bigl(H_5+(125+100q)H_6\bigr).
                                                                    \tag{13}
\]
The Galois substitution \(q\mapsto-q\) gives the other branch.
The value \(z=0\) is included in this polynomial parametrization;
there is no division by \(z\).

The witness retains eleven admitted polynomial tuples
\[
 P_1,P_{aa},P_{ab},P_{ac},P_{bb},P_{bc},P_{cc},
                P_{az},P_{bz},P_{cz},P_{zz},
\]
and verifies the exact coefficient identities proving that
\[
\begin{split}
 \Phi(a,b,c,z)={}&P_1+a^2P_{aa}+abP_{ab}+acP_{ac}
                   +b^2P_{bb}+bcP_{bc}+c^2P_{cc}\\
                &+azP_{az}+bzP_{bz}+czP_{cz}+z^2P_{zz}
\end{split}                                                    \tag{14}
\]
solves (8) for every \(a,b,c,z\). In fact \(P_{bb}=0\).
Each polynomial coefficient is stored as a pair
\[
                  [r_0,r_1]\quad\hbox{meaning }r_0+qr_1,
                  \qquad r_0,r_1\in\mathbb Q.                 \tag{15}
\]
The eleven tuples have 797 coefficient-pair entries, containing 1224
nonzero rational components. The replay expands the branch
bilinearly with \(q^2=-2\), and checks the rational and \(q\) components
separately. No algebraic solver or numerical specialization is used.

The matrix \(\mathcal A\) has rank 344. A specified 344-square minor
has determinant \(381548\) modulo \(1000003\), and its three exact
independent homogeneous kernels are
\[
\begin{split}
 K_{R^3}&=(0,0,-3R^2,0),\\
 K_{R^4S}&=(0,R^4,-4R^3S,0),\\
 K_{X^3}&=\{X^3,f\}.
\end{split}                                                    \tag{16}
\]
The third has nonzero \(S\)-component, the second has nonzero
\(T\)-component and zero \(S\)-component, and the first has only
a nonzero \(D\)-component. Thus these are independent. The ENTIRE
admitted second-coefficient family over (13) is
\[
                g^{(2)}=\Phi(a,b,c,z)
                         +uK_{R^3}+vK_{R^4S}+\gamma K_{X^3},   \tag{17}
\]
with all three additional parameters retained.

Equations (14)--(17) also prove sufficiency in the field statement
(9). If \(f\ne0\), (11) supplies an element
\(q=(d/f-125)/100\) with \(q^2=-2\), so the corresponding embedding
of \(\mathbb E\) and (13) give the desired point. If \(f=0\), then
\(d=0\), and (14) at \(z=0\) has rational coefficients and descends
to every characteristic-zero field. This also handles the nonsquare
case without adjoining a root unnecessarily.

The first image actually changes on a nonzero complex branch:
\[
 \begin{split}
 R_1&=(125+100q)zR_0^4,\\
 R_2&=-\frac{4800}{7}(1+2q)z^2R_0^7.
 \end{split}                                                   \tag{18}
\]
All three homogeneous second kernels have zero \(R\)-component.
These are the first two coefficients only. Equation (18) does not
assert that all later first-image coefficients lie in a polynomial
ring generated by \(R_0\).

## Replay and scope

[The checker](check_graded_long_cap24_complex_family.py) uses Python
integers and Fraction, [the complete column engine](check_graded_long_cap23.py),
and [the independent polynomial curvature arithmetic](check_long_graph_curvature.py).
The canonical producer supplies input polynomials only.
[The retained receipt](graded_long_cap24_complex_controls.json) pins
repository-relative sources, the two rank minors, eight-kernel
independence, the three necessity functionals, all eleven exact
field-valued points and all three second kernels.

The classification is for characteristic-zero fields through
commutator order \(h^3\). It makes no nonreduced-base claim and no
higher compatibility, bounded full-lift, evenization or DC1/DC2 claim.
The surviving family at this order is not an assertion that it
extends at the next order.

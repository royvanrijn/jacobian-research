# Transverse double lines and conics in the rank-two quintic pencil

**Theorem QuinticRankTwoTransverseDoubleLineConicPencils.** Work over a
field of characteristic zero. Suppose a rank-two quintic Keller map
has cancelling quartic leading part \(R_4=z^2q\), where the binary
restriction of \(q\) to \(z=0\) has rank two. There is no Keller
completion with either of the following leading pairs:
\[
\begin{aligned}
 (P_5,Q_5)&=(z^5,\ b z^3q+\rho zq^2), &&\rho\ne0,\\
 (P_5,Q_5)&=(z(q^2+\delta z^4),\ z^3q+\nu z^5),
                                      &&(\delta,\nu)\ne(0,0).
\end{aligned}
\]
There is no assumption on \(R_3\) or on irreducibility of \(q\).

For the [degree-two repeated-line pencil](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md),
this leaves the two quintic spans \(\langle z^5,z^3q\rangle\) and
\(\langle zq^2,z^3q\rangle\) when \(q=xy\).
The same two spans are the ones outside the determinant proof for
a transverse irreducible conic. The tangent case is outside this
theorem. The use of the binary Jacobian invariant ring is the
same as for the square pencil with a split binary quadratic in
[Wang and Yan, §3.4.1](https://arxiv.org/pdf/2610.36495v1).

## Laurent modes and the cofactor

Embed the finitely generated coefficient field in \(\mathbb C\).
Make the binary restriction of \(q\) equal to \(xy\), and add
multiples of \(z\) to \(x,y\) to obtain \(q=xy+\gamma z^2\).
At \(z=1\) use \(r=q(x,y,1)\) and \(y=(r-\gamma)/x\).
The monomials \(1,x^i,y^j\), \(i,j\ge1\), form a degree-preserving
free basis over \(\mathbb C[q,z]\) before dehomogenization.
At \(z=1\) they give distinct Laurent modes in \(x\).
In particular the invariant ring of \(D=x\partial_x|_r\)
in the dehomogenized polynomial ring is \(\mathbb C[r]\).

Let \(H_i\) dehomogenize \((P_{5-i},Q_{5-i},R_{4-i})\), and let
\(E_i=\operatorname{diag}(5-i,5-i,4-i)H_i\). As in the
[reducible-square calculation](QUINTIC_KELLER_RANK_TWO_REDUCIBLE_SQUARE_CONIC_PENCIL.md#the-mixed-span-at-zero-parameters),
the homogeneous determinant of degree \(11-n\) at \(z=1\),
up to a common sign, is
\[
         J_n=\sum_{i+j+k=n}\det(E_i,\partial_rH_j,DH_k).
 \tag{1}
\]
Both derivations preserve Laurent modes. A component of original
degree \(d\) has no mode with absolute value greater than \(d\),
and the coefficient of its positive mode \(x^m\) has
degree at most \(\lfloor(d-m)/2\rfloor\) in \(r\).
Write \(\mathcal L=E_0\times\partial_rH_0\), so
the linear term in (1) at this step is \(D(\mathcal L\cdot H_n)\).
The coefficient identities below are verified by the
[exact Laurent-mode replay](verify_quintic_transverse_double_line_conic.py).

## The span containing \(z^5\)

Here \(H_0=(1,br+\rho r^2,r)\) and
\(\mathcal L=(br-3\rho r^2,-5,5(b+2\rho r))\).
For a positive mode in \(H_1\), let \((p,t,s)\) be its coefficient
triple. For mode two the bounds in \(r\) are affine, affine,
constant. The linear equation \(\mathcal L\cdot(p,t,s)=0\)
then forces \(p=0\) and \(t=(b+2\rho r)s\), by its
coefficients in degrees three and two and by \(\rho\ne0\).
For mode at least three it forces the triple to vanish.
Thus the mode-two coefficient is a scalar multiple of
\(V=(0,b+2\rho r,1)\). In \(J_2\) the \(x^4\) mode comes
only from two such mode-two terms of \(H_1\).
For coefficient \(A V\) it is \(20\rho A^2\).
Hence \(A=0\), and interchanging \(x,y\) eliminates the
negative mode of order two as well.

Mode one of \(H_1\) allows an affine \(s\). Solving the
same cofactor equation gives the two directions \(V\) and
\(W=(10/3,5br/3,r)\). The linear term in direction \(V\)
comes from translating \(q\) in the leading triple.
Indeed the homogeneous \(q\)-derivative of that triple is
\((0,z(bz^2+2\rho q),z^2)\). Any binary linear correction
to \(q\) can be induced by translating \(x,y\), so remove
the \(V\) direction in both modes of order one.
Let the remaining positive mode of \(H_1\) be \(\alpha x W\).

The \(x^2\) part of \(H_2\) has a constant coefficient triple.
For that triple \(J_2=0\) has the unique solution
\[
     x^2[ x^2 ]H_2
       =\alpha^2x^2(35/9,-5b^2/(18\rho),-5b/(18\rho)).
\]
Uniqueness follows from independence of the three coordinates
of \(\mathcal L\) over the constants. The only possible
mode-three terms of \(H_2\) satisfy the homogeneous cofactor
equation, since \(H_1\) now has only modes zero and one in
absolute value. That equation has zero constant kernel.
The \(x^3\) term in \(J_3\) therefore comes only from the
displayed highest modes of \(H_1,H_2\) and the leading triple.
It is
\[
           \frac{5\alpha^3}{9\rho}
                  (-5b+14r\rho)(b+2r\rho)x^3 .
\]
This forces \(\alpha=0\) and symmetrically kills the other
mode of order one. Hence \(H_1\) is radial, meaning independent
of \(x\) in the Laurent chart.

Now the linear equation for \(H_2\) has only nonradial
directions \(xV,yV\), by the bounds on \(H_2\).
In \(J_3\), any mode of order two from \(H_3\) must have
a constant coefficient in the kernel of \(\mathcal L\).
It therefore vanishes. The mode-two term of \(J_4\) is
then \(10\rho M^2x^2\) when \(M x V\) is the positive mode
of \(H_2\). Thus \(M=0\) and also the coefficient of \(y V\)
vanishes. Once \(H_1,H_2\) are radial, the bounds give only
constant positive-mode coefficient triples in \(H_3,H_4\).
The zero constant kernel makes those corrections radial,
contrary to the nonzero full determinant.

## Mixed spans

Now \(H_0=(r^2+\delta,r+\nu,r)\) and
\[
           \mathcal L=(r+5\nu,\ 3r^2-5\delta,\
                              -5r^2-10\nu r+5\delta).
\]
On bounds affine, affine, constant in \(r\), the kernel is
spanned by \(V=(2r,1,1)\) when \(\delta\ne15\nu^2\).
When \(\delta=15\nu^2\), add \(W=(r-5\nu,-1/3,0)\).
Indeed if \((p,t,s)=(a_1r+a_0,b_1r+b_0,s)\), then
\(b_1=0\), \(b_0=(5s-a_1)/3\), \(a_0=5\nu(2s-a_1)\), and
\((5\nu^2-\delta/3)(2s-a_1)=0\).
In the exceptional case the constant kernel is spanned
by \(K=V-2W=(10\nu,5/3,1)\). In the generic case it is zero.
With all three bounds affine the additional direction is
\(Z=(10\delta/3,5(r+2\nu)/3,r)\), since its third entry
accounts for the linear coefficient in \(s\).

For \(H_1\) the mode-two bounds are affine, affine, constant
and the mode-one bounds are all affine. Modes of order at
least three have zero kernel in the generic case. In the
exceptional case mode four has zero kernel since it lacks the
third entry of \(K\), but mode three can have coefficient
\(sK\). Its contribution to \(J_2\) at mode six is
\(-10s^2(r-3\nu)x^6\), from two copies of \(s x^3K\).
The linear term on \(H_2\) cannot reach mode six.
Thus \(s=0\), and interchange of \(x,y\) also eliminates
the opposite mode of order three.
In the generic case, the known term of \(J_2\) at \(x^4\),
from a mode \(A x^2V\) of \(H_1\), is
\(-4A^2(r+5\nu)x^4\). Hence \(A=0\). In the exceptional
case assume \(\nu\ne0\). For mode \(x^2(sV+tW)\) this term
is \(-4x^4 G(s,t)/3\), where
\[
                G(s,t)=15\nu s(s+t)+r(3s^2+3st+2t^2).
\]
As a polynomial in \(r\), \(G(s,t)=0\) implies
\(s(s+t)=0\) and then \(t=s=0\).
In both cases interchange \(x,y\) to remove the other
mode of order two.

The linear \(V\) terms of \(H_1\) are again removed by
translation of \(q\) in the leading triple. The remaining
positive mode is \(\alpha x Z\) in the generic case and
\(x(\alpha Z+tW)\) in the exceptional case. In either
case its contribution at \(x^2r^3\) in \(J_2\) is
\(-10\alpha^2/3\). The mode-two coefficients of \(H_2\)
are constant, so the cofactor cannot match the cubic in
\(r\). Thus \(\alpha=0\). In the exceptional case the
known term for \(t xW\) is \(-4r t^2 x^2/3\).
If the mode-two coefficients of \(H_2\) are \((p,t_2,s)\),
coefficient comparison at \(r^2,r,r^0\) in \(J_2=0\) gives
\[
      3t_2=5s,\qquad p-10\nu s=2t^2/3,\qquad
                       5\nu p-75\nu^2t_2+75\nu^2s=0.
\]
Substitution into the last equality gives \((10/3)\nu t^2=0\),
so \(t=0\). Interchanging \(x,y\) eliminates the other mode
of order one and \(H_1\) is radial.

For the generic case, \(H_2\) therefore has only the nonradial
directions \(xV,yV\), since its mode coefficients of order
at least two are constant. As in the first span, \(H_3\)
has no mode of order two. For a positive mode \(M xV\) in
\(H_2\), \(J_4\) has \(x^2\) term \(-2M^2(r+5\nu)x^2\).
It follows that \(H_2\) is radial. The remaining corrections
have zero nonradial kernel and hence the whole map is radial.

For the exceptional case, in \(H_2\) the mode-two kernel is
the constant direction \(K\), and the mode-one kernel is
spanned by \(V,W\). There is no mode three, since \(R_2\) of
degree two cannot supply the nonzero third coordinate of \(K\).
In \(J_4\), a mode \(s x^2 K\) of \(H_2\) contributes
\(-20s^2(r-3\nu)x^4/3\) at mode four. Thus \(H_2\)
has no modes of order two. Now \(H_3\) also has no modes of
order two: the known part of \(J_3\) has no such mode and
\(R_1\) of degree one cannot supply the third entry of \(K\).
It follows that mode two of \(J_4\) comes only from the
two drop-two columns. For \(x(sV+tW)\) it is
\(-2x^2G(s,t)/3\) and forces \(s=t=0\). By the same
argument in the other mode, \(H_2\) is radial.
Then \(H_3\) has at most the nonradial constant directions
\(xK,yK\). A positive mode \(s xK\) in \(H_3\)
gives the mode-two term of \(J_6\), namely
\(-10s^2(r-3\nu)x^2/3\). No contribution involving \(H_4\)
can reach this mode, since its degrees are at most one
and \(H_1,H_2\) are radial. Thus \(H_3\) is radial, and
the cofactor equation on \(H_4\) has zero nonradial kernel
because \(R_0\) is constant. This again contradicts the
Keller condition.

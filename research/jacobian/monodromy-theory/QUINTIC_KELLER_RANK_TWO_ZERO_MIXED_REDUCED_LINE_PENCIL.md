# The zero-parameter mixed pencil with a double line and a reduced line pair

**Theorem QuinticRankTwoZeroMixedReducedLinePencil.** Over a field of
characteristic zero there is no three-variable Keller map with ordinary
homogeneous leading parts
\[
                 P_5=z(xy)^2,\qquad Q_5=z^3xy,\qquad R_4=z^2xy .
\]
This is the mixed zero-parameter case in the
[transverse double-line pencil](QUINTIC_KELLER_RANK_TWO_TRANSVERSE_DOUBLE_LINE_CONIC_PENCILS.md).
The proof has no hypothesis on the cubic part of \(R\).

## Cofactor and mode bounds

Embed the finitely generated coefficient field in \(\mathbb C\).
Set \(q=xy\) and use source coordinates \((z,q,x)\) on \(xz\ne0\).
When \(z=1\) write \(r=q\), so \(y=r/x\).
The ordinary polynomial ring is free over \(\mathbb C[q,z]\)
with basis \(1,x^i,y^j\) for \(i,j\ge1\), by the relation \(xy=q\).
The basis elements have distinct Laurent modes in \(x\) at \(z=1\).
The invariant ring in the dehomogenized polynomial ring of
\(D=x\partial_x|_r\) is therefore \(\mathbb C[r]\).
This is the split binary quadratic invariant calculation in
[Wang and Yan, §3.4.1](https://arxiv.org/pdf/2610.36495v1).

Let \(H_i\) be \((P_{5-i},Q_{5-i},R_{4-i})\) on that chart,
with \(H_0=(r^2,r,r)\).
Write \(E_i=\operatorname{diag}(5-i,5-i,4-i)H_i\).
The determinant part of degree \(11-n\), restricted to \(z=1\)
and up to the common sign, is
\[
               J_n=\sum_{i+j+k=n}
                         \det(E_i,\partial_rH_j,DH_k).
 \tag{1}
\]
Indeed the coordinate change has determinant \(-x\).
At \(z=1\) use the Euler column
\(z\partial_z+2q\partial_q+x\partial_x\) for the first column,
and put the extra factor \(x\) into \(D\).
The nonzero-constant Jacobian condition requires \(J_n=0\)
for \(n<11\) and \(J_{11}\ne0\).

The linear part of (1) in \(H_n=(p,t,s)\) is
\[
                        r D(p+3rt-5rs).
 \tag{2}
\]
At positive mode \(x^m\), a component of original degree \(d\)
has coefficient of \(r\)-degree at most
\(\lfloor(d-m)/2\rfloor\), and has zero coefficient when \(m>d\).
The bounds on negative modes follow by swapping \(x,y\).
The derivatives in (1) preserve each Laurent mode. Constant
homogeneous parts of \(P,Q,R\) give zero derivative columns.

For mode coefficients in (2) define
\[
 V=(2r,1,1),\quad W=(r,-1/3,0),\quad
 Z=(0,5r/3,r),\quad K=V-2W=(0,5/3,1).
\]
On constant triples the kernel is the span of \(K\).
On affine, affine, constant triples it is the span of \(V,W\),
and on three affine entries it is the span of \(V,W,Z\).
These descriptions follow by equating coefficients of
\(p+3rt-5rs=0\).
In (1) the relevant two polynomial expressions in the
nonradial coefficients are
\[
          G(s,t)=3s^2+3st+2t^2,\qquad
          T(s,t)=s(2s+t)(3s-t).
\]
Their only common zero is \((0,0)\). Indeed the three factors
of \(T\) give \(s=0\), \(t=-2s\), or \(t=3s\), and in each
case \(G=0\) forces both scalars to vanish.
The mode coefficients below are identities from (1),
replayed by the
[exact coefficient checker](verify_quintic_zero_mixed_reduced_line_pencil.py).

## Eliminating the first correction

Equation \(J_1=0\) shows that the positive Laurent modes of
\(H_1\) at orders three, two and one respectively have coefficient
spaces \(\mathbb C K\), \(\langle V,W\rangle\) and
\(\langle V,W,Z\rangle\). A mode of larger order is zero because
its third entry is zero and the constant kernel requires a
nonzero third entry. If the mode-three term is \(s x^3K\),
the only contribution to mode six of \(J_2\) is quadratic in
that term, and is \(-10rs^2x^6\). Thus the third mode vanishes
in each direction.

If the positive mode-two term is \(x^2(sV+tW)\), then \(J_2\)
has mode-four term \(-4rG(s,t)x^4/3\).
The mode-six term of \(J_3\) is \(-2T(s,t)x^6/3\), since
\(H_2\) has modes of order at most three and no contribution
from \(H_0,H_1,H_2\) can reach order six. These two equations
kill the second mode in each direction.

Remove the mode-one terms in the \(V\) direction by translating
\(x,y\). In homogeneous coordinates the \(q\)-derivative of
the leading triple is \((2zq,z^3,z^2)\), and translating \(x,y\)
induces any binary linear correction to \(q\).
The possible remaining positive mode of \(H_1\) is now
\(x(aZ+tW)\). Its contribution to \([x^2r^3]J_2\) is
\(-10a^2/3\), while the mode-two coefficients in \(H_2\) are
constant, with no cubic in \(r\) in their contribution via (2).
Hence \(a=0\).

For the remaining \(txW\) in \(H_1\), the known mode-two term
in \(J_2\) is \(-4rt^2x^2/3\). The equation therefore gives
\(x^2(2t^2/3,5s/3,s)\) as the mode-two term of \(H_2\).
The mode-three term of \(H_2\) must be zero by (2), since
the known contribution has no mode three and \(R_2\) cannot
supply the third entry of \(K\) at that mode. In \(J_3\)
the \(x^3\) term is consequently
\[
                         \frac{5t}{9}(9rs-2t^2)x^3 .
\]
It forces \(t=0\). Interchanging \(x,y\) also eliminates the
other mode of order one. This makes \(H_1\) radial, meaning
independent of \(x\) at fixed \(r\).

## Eliminating the lower corrections

Now \(J_2=0\) is the homogeneous cofactor equation.
It gives at most the direction \(K\) for the second mode of
\(H_2\) and the directions \(V,W\) for its first mode.
There is no mode three because \(R_2\) has degree two.
For a positive term \(s x^2K\), the \(x^4\) term of \(J_4\)
is \(-20rs^2x^4/3\) from the two drop-two columns.
So \(H_2\) has no modes of order two.
Now \(J_3=0\) has no known terms of order two, and \(H_3\)
cannot support that mode homogeneously in (2), since
\(R_1\) is linear.

Write \(x(sV+tW)\) for the positive mode of \(H_2\).
In \(J_4\) its mode-two contribution is \(-2rG(s,t)x^2/3\),
with no other terms at that mode. The mode-three term of
\(J_6\) is \(-T(s,t)x^3/3\). Indeed \(H_0,H_1\) are radial,
\(H_2,H_3,H_4\) have modes of order at most one, and the
only summand of total drop six using three nonradial columns
is the one with three drop-two columns. Thus \(H_2\) is radial
by the common-zero calculation and interchange of \(x,y\).

With \(H_1,H_2\) radial, \(J_3=0\) permits only \(xK,yK\)
as the nonradial terms in \(H_3\). For positive mode
\(s xK\), the \(x^2\) term in \(J_6\) is
\(-10rs^2x^2/3\) from two drop-three columns.
No other term can reach that mode, because earlier drops
are radial and \(H_4\) is at most linear.
It follows that \(H_3\) is radial. Finally \(J_4=0\) makes
\(H_4\) radial, since \(R_0\) is constant and the kernel (2)
on constant triples with zero third entry is zero.
All components then belong to \(\mathbb C[q,z]\).
Their determinant vanishes, a contradiction.

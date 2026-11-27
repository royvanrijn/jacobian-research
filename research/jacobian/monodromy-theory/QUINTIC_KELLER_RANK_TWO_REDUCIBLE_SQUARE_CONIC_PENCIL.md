# A square of a reduced line pair in the rank-two conic pencil

**Theorem QuinticRankTwoReducibleSquareConicPencil.** In characteristic
zero there is no three-variable Keller map whose quintic leading span
has dimension two and is contained in
\[
                    \langle z^5,z^3q,zq^2\rangle,\qquad q=xy,
\]
and whose cancelling quartic has leading part \(q^2\).
In the [primitive degree-two repeated-line pencil](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md),
this excludes \(R_4=J^2\) when \(J\) is a product of distinct lines.

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. For the last assertion, neither factor of \(J\)
can be \(z\), since \(J\) is coprime to the other member \(z^2\).
If their restrictions to \(z=0\) were proportional, both factors
would belong to \(\mathbb C[z,L]\) for a linear \(L\), and \(J/z^2\)
would be a quadratic polynomial of \(L/z\). This contradicts
primitivity of the generator \(J/z^2\). The two factors and \(z\)
are therefore independent and can be made \(x,y,z\).

## The spans and the cofactor arguments

On \(z=1\) write \(r=xy\). Its polynomial Jacobian invariant
ring is \(\mathbb C[r]\), by the monomial eigenvalues of
\(D=x\partial_x-y\partial_y\), as in
[Wang and Yan, §3.4.1](https://arxiv.org/pdf/2610.36495v1).
There is a free basis \(1,x^i,y^j\) for \(i,j\ge1\) over
\(\mathbb C[q,z]\), by \(xy=q\) and the distinct Laurent modes
at fixed \(r\). The basis \(x^i y^j\), \(i=0,1\), used in the
generic cofactor proof is also available after a binary
linear change giving \(q\) a nonzero \(x^2\) coefficient.
Both bases respect homogeneous degrees.

If the span contains \(z^5\), choose leading parts
\(P_5=z^5\), \(Q_5=z^3q+\rho zq^2\) when the \(z^3q\) coefficient
is nonzero, and \(P_5=z^5\), \(Q_5=zq^2\) otherwise.
Both the [pure-fifth nonzero-coefficient proof](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_PENCIL_RESIDUAL.md#removing-the-nonzero-coefficient)
and the [pure-fifth zero-coefficient proof by Laurent modes](QUINTIC_KELLER_RANK_TWO_PURE_FIFTH_TRANSVERSE_DOUBLE_CONIC.md)
apply with \(r=xy\). Their cofactor calculations use only the invariant
ring, free basis, and translation of the binary linear term of \(q\).
Translations in \(x,y\) induce every binary linear correction to \(xy\).
Any multiple of \(z\) in the translating linear form can be absorbed
into the pieces in \(\mathbb C[q,z]\). In the Laurent proof this
makes its use of \(\gamma\ne0\) unnecessary, since the mode formulas
hold at \(r=xy\), including the geometric normalization of \(R_3\).

For \(P_5=z^5,Q_5=zq^2\), use the
[quartic-argument proof](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_QUARTIC_PENCIL.md)
with the same invariant ring and translations. Its other use of
irreducibility is to force a binary linear \(U\) to vanish when a
nonconstant quadratic polynomial of \(U\) is constant on \(r=0\).
For \(r=xy\), the quadratic relation makes \(U\) constant on
each axis, and those constants agree at the origin. A binary linear
\(U\) constant on both axes is zero. The polynomial division by \(r\)
and \(r^2\) in that proof does not require primality of \(r\).

If the span does not contain \(z^5\), its unique direction with zero
\(zq^2\) coefficient has a nonzero \(z^3q\) coefficient. Choose
\[
               P_5=z(q^2+\delta z^4),\qquad Q_5=z^3q+\nu z^5 .
\]
The [generic](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md)
and [nonzero exceptional](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_EXCEPTIONAL_PENCILS.md)
arguments apply at \(q=xy\) with the invariant ring and
translation above. The exceptional proof already treats a common
cofactor level consisting of two intersecting lines. In this case
the level is either of that form or is the irreducible curve
\(xy=r_*\) with \(r_*\ne0\). These arguments leave
\(\delta=\nu=0\).

## The mixed span at zero parameters

On \(xz\ne0\) use \((z,q,x)\), whose change of coordinates has
determinant \(-x\). At \(z=1\) use \(y=r/x\).
Write \(H_i\) for the dehomogenized triple
\((P_{5-i},Q_{5-i},R_{4-i})\), with \(H_0=(r^2,r,r^2)\).
If \(E_i=\operatorname{diag}(5-i,5-i,4-i)H_i\) and
\(D=x\partial_x|_r\), the degree-\((11-n)\) determinant at
\(z=1\), up to overall sign, is
\[
             J_n=\sum_{i+j+k=n}\det(E_i,\partial_rH_j,DH_k).
 \tag{1}
\]
This is the Euler-column formula with the factor \(x\) from the
coordinate change included in \(D\). In original binary variables
\(D=x\partial_x-y\partial_y\). The cofactor expression on
\(H_n=(p,t,s)\) is \(r^2D(6p-2rt-5s)\), from
\(E_0\times\partial_rH_0\).

The kernel of \(J_1\), modulo pieces in \(\mathbb C[q,z]\),
gives
\[
        (P_4,Q_4,R_3)=(q,3z^2,0)U_2+(5z/6,0,1)V_3
                       \pmod{\mathbb C[q,z]^3}
 \tag{2}
\]
with homogeneous \(U_2,V_3\) of degrees two and three.
Indeed the homogeneous cofactor condition on a nonradial basis
coefficient \((p,t,s)\) is \(6z^2p-2q t-5z^3s=0\).
It makes \(p-5zs/6\) divisible by \(q\) and \(t\) equal to
three times its quotient multiplied by \(z^2\).

In \(J_2\), reduce modulo \(r\) before replacing \(y\) by \(r/x\).
The two leading columns there are \(E_0=0\) and
\(\partial_rH_0=(0,1,0)\), so only the determinant of the first
and third entries of \(E_1,DH_1\) remains. If \(v\) is the
restriction of \(V_3\) at \(z=1\), and \(a,g\) are the constant
parts of the radial corrections to \(P_4,R_3\) on \(r=0\),
that determinant is, up to sign,
\[
                      \big(5v/6+4a-5g/2\big)Dv .
\]
On each axis \(D\) is a nonzero multiple of the coordinate
Euler derivation, whose polynomial kernel there is the scalars.
The product therefore makes \(v\) constant on each axis, with
the same constant at the intersection. Thus \(V_3=qL\) modulo
\(\mathbb C[q,z]\) for a binary linear \(L\). In the leading
triple \(H=(zq^2,z^3q,q^2)\), one has
\[
                (5z/6,0,1)qL
          = \tfrac12 H_q L-\tfrac16 z(q,3z^2,0)L .
\]
A translation removing the first term leaves only the first
direction of (2) modulo radial pieces. Consequently, with
\(V=(r,3,0)\), the nonradial part of \(H_1\) is \(U V\) where
\(U=A x^2+B y^2+\alpha x+\beta y\).

A term \(x^k r^j\) has Laurent mode \(k\).
Both \(\partial_r\) and \(D\) preserve the mode, and a term
of original ordinary degree \(d\) has no mode with
\(|k|>d\). Thus the only contribution to mode four of \(J_2\)
comes from the \(A x^2 V\) terms in two drop-one positions.
Its coefficient in (1) is \(24r^2 A^2\), so \(A=0\).
Interchanging \(x,y\) gives \(B=0\).

Now \(H_1\) has modes at most one in absolute value. In \(J_2\)
all known terms have mode at most two, so a mode three of \(H_2\)
would satisfy \(6p-2rt=0\) on constant mode coefficients
\(p,t\), with no mode three in \(R_2\). Thus \(H_2\) has no
mode three in either direction. At mode two its three
coefficients are constants by the degree bounds, and the known
part of \(J_2\) from \(\alpha x V\) is \(12r^2\alpha^2 x^2\).
The mode-two part of \(H_2\) is therefore
\[
                      x^2(5s/6-\alpha^2,0,s)
\]
for a scalar \(s\). The mode-three coefficient in \(J_3\)
is \(-15\alpha r(s+2\alpha^2)\); it only uses the displayed
highest modes of \(H_1,H_2\) and the leading triple, since
\(H_3\) has original degree at most two. The mode-four
coefficient in \(J_4\) is \(5s(s-12\alpha^2)/3\), since
\(H_1,H_3\) together cannot reach that mode. The two equations
are inconsistent if \(\alpha\ne0\). Swapping \(x,y\) gives
\(\beta=0\) as well.

With \(H_1\) now radial, \(J_2\) is the linear cofactor equation.
If \(W=(5/6,0,1)\), the nonradial part of \(H_2\) has the form
\(V U_1+W V_2\), with \(U_1\) binary linear and
\(V_2\) in the span of \(x^2,y^2,x,y\), by the degree bounds
in (2) after decreasing the degrees by one. The mode-four
coefficient in \(J_4\) is \(5s^2/3\), where \(s\) is the
coefficient of \(x^2\) in \(V_2\). Thus \(V_2\) has only
linear terms by the same argument in the \(y\) direction.
Now \(J_3\) has no known terms of mode two, and the homogeneous
cofactor equation \(6p-2rt=0\) on constant mode-two
coefficients of \(H_3\) eliminates that mode in both directions.
With \(x(uV+tW)\) the positive mode of \(H_2\), the mode-two
coefficient of \(J_4\) is therefore
\[
                       (72r^2u^2-60rtu+5t^2)/6 .
\]
It forces \(u=t=0\), so \(H_2\) is radial by also interchanging
\(x,y\).

With \(H_1,H_2\) radial, \(J_3\) is the cofactor equation.
Its only possible nonradial part of \(H_3\) is \(W(lx+my)\)
for scalars \(l,m\). The mode-two coefficient of \(J_6\)
is \(5l^2/6\) and comes from two drop-three columns and one
leading column. Indeed \(H_4\) has modes at most one and
all lower drops than three are radial. Hence \(l=0\), and
\(m=0\) by symmetry. With all earlier drops radial, the
cofactor condition \(J_4=0\) makes \(H_4\) radial because
\(\deg P_1,\deg Q_1\le1\) and \(R_0\) is constant.
This would make the entire map polynomial in \(q,z\), contrary
to the Keller condition.

The isolated coefficients in (1) are replayed with exact
rational arithmetic by
[verify_quintic_reducible_square_conic_pencil.py](verify_quintic_reducible_square_conic_pencil.py).

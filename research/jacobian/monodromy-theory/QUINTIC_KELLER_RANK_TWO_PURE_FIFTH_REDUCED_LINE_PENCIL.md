# The pure-fifth pencil with a reduced line pair

**Proposition KellerQuarticReducedLinePairCubicVanishing.** Let \(R\)
be a quartic component of a noninvertible three-variable Keller map
of degree at most five in characteristic zero. If \(R_4=z^2xy\)
in independent source directions, then \(z\mid R_3\).

**Theorem QuinticRankTwoPureFifthReducedLinePencil.** A
three-variable Keller map of characteristic zero with leading
forms \(P_5=z^5,\ Q_5=z^3xy,\ R_4=z^2xy\) is invertible.

## The cubic along the repeated line

Work over \(\mathbb C\) by embedding the finitely generated field
of coefficients. Assume \(R_3|_{z=0}\ne0\). The projective quartic
fibre has isolated singularities since \(R\) is a submersion and
the only repeated line at infinity is \(z=0\). Away from the two
triple intersections on that line, the
[double-line tests](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md#points-away-from-the-pencil-base-point)
at zeros of the cubic either extend the residue form on resolution
or give a nonzero period. The other intersection, of the two
simple lines, has a plane quadratic jet of rank two. At either
triple intersection a full jet of rank at least two gives the
[ordinary double-point extension](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve),
and a rank-one jet is a pure normal square and gives
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
Exactness of the affine fibre form excludes a holomorphic extension
on the whole projective resolution. Thus at least one of the two
triple intersections must have vanishing linear and quadratic jets.

Interchange \(x,y\) to take that point as the \(y\) direction.
Then \(R=a(x,z)y+b(x,z)\) with \(a_3=z^2x\) and \(b_4=0\).
The [cubic source-coefficient classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
excludes the repeated-parallel residual, whose leading coefficient
is a cube. Its other two cases allow affine coordinates with
leading directions \(z,x\), in which \(a=z^2x+\delta z-1\)
or \(a=z(zx-1)\), up to scaling. In the first case the
representative for \(b\) modulo \(a\) is affine in \(zx\).
In the second it is
\[
       (1-zx)(\lambda_0 x+\mu_0)+zx(\lambda_1 x+\mu_1),
\]
by the two restrictions at \(z=0\) and \(zx=1\).
Both representatives have degree at most three with every
cubic term divisible by \(z\). The quotient of \(b\) minus the
representative by \(a\) is constant, since \(\deg b\le3\).
Therefore all cubic terms of \(R\) are divisible by \(z\).
Translations used in the affine normalization do not affect this
conclusion because the first derivatives of \(z^2xy\) are
divisible by \(z\). This contradicts the standing assumption
and proves the proposition.

## The first coefficients in the pure-fifth span

For the theorem suppose the map is noninvertible. Write \(q=xy\),
dehomogenize at \(z=1\), set \(r=q(x,y,1)\), and use
\(y=r/x\) on \(x\ne0\).
Let \(H_i\) denote the dehomogenized piece
\((P_{5-i},Q_{5-i},R_{4-i})\) and put
\(E_i=\operatorname{diag}(5-i,5-i,4-i)H_i\),
\(D=x\partial_x|_r\).
The [Laurent determinant formula](QUINTIC_KELLER_RANK_TWO_TRANSVERSE_DOUBLE_LINE_CONIC_PENCILS.md#laurent-modes-and-the-cofactor)
gives, up to sign,
\[
       J_n=\sum_{i+j+k=n}\det(E_i,\partial_r H_j,DH_k)
\]
for the degree-\((11-n)\) equation. Here \(H_0=(1,r,r)\),
and the cofactor on \(H_n=(p,t,s)\) is \(D(rp-5t+5s)\).
The kernel of \(D\) in the polynomial ring in \(x,y\) is
\(\mathbb C[r]\), as in
[Wang–Yan, §3.4.1](https://arxiv.org/pdf/2610.36495v1).
We call the Laurent mode of \(x^i r^j\) its exponent \(i\)
in \(x\); the zero mode is radial. Nonradial
coefficients in \(\mathbb C[x,y]\) can be written uniquely on
the basis \(1,x^i,y^i\), \(i\ge1\), over \(\mathbb C[r]\).

Use the cofactor equation at \(n=1\) and \(z\mid R_3\).
At modes of order at least three the cofactor and degree bounds
force all nonradial terms to vanish. At mode two the leading
coefficients in \(r\) lie in the span of \(V=(0,1,1)\) and
\(W=(1,r/5,0)\), and the same holds at mode one: \(R_3\)
has binary degree at most two, so its coefficient at each of
these modes is scalar on \(x^i,y^i\).
Remove the linear terms in the \(V\) direction by translating
\(x,y\) in the leading triple. Remove the radial term \(zq\) of \(R_3\) by translating \(z\),
and remove the \(r\) terms of \(P_4,Q_4\) by subtracting
multiples of \(R\). More explicitly the \(z\) translation
changes \(R_3\) by a constant multiple of \(zq\),
\(Q_4\) by a multiple of \(z^2q\), and \(P_4\) by a multiple
of \(z^4\). In the resulting frame write
\[
 H_1=(a+t_0r^2,d+t_1r^2,g)+(A x^2+C y^2)V+
                      (B x+F y+B_2x^2+F_2y^2)W .
\]
At \(x^4r\), \(J_2=0\) forces \(B_2=0\), since the
coefficient after dividing by the mode is \(-2B_2^2/5\).
The \(y\) direction is symmetric. In the
[exact replay](verify_quintic_pure_fifth_reduced_line_pencil.py)
the now remaining compatibility equations include
\(-2 A t_0\) from order two, and
\(-5A^2t_1,\ 2B(50At_1-B^2)/25\) from order three.
It follows that \(B=0\) and, by symmetry, \(F=0\).
Any nonzero \(A\) or \(C\) forces \(t_0=t_1=0\).

The replay uses the basis \(x^i y^j\) in each unknown homogeneous
layer, with \(i+j\) bounded by the ordinary degree of the component.
Integrating known modes in the \(D\) direction divides a nonzero
mode coefficient by its integer exponent and retains the zero
mode as a separate compatibility condition. The linear
coefficient from the cofactor is imposed alongside it.
Computing a particular solution by the first independent
columns and rows and adding the ordered nullspace basis
parametrizes every solution; unmatched coefficients give the
compatibility equations quoted below. Thus the replay concerns
identities in arbitrary coefficients, with no numerical search or
assumption of generic parameters.

## Nonzero radial quadratic terms

If \(t_0\ne0\) or \(t_1\ne0\), the coefficients \(A,C\)
are both zero, so \(H_1\) is radial. The order-two linear
equation gives
\[
 H_2=(f_0+f_1r+5 f_4 y+5 f_5x,\
       f_2+f_3r+ f_4yr+f_5xr+f_7 y+f_8y^2+f_9x+f_{11}x^2,\
       f_6+f_{10}r+f_7 y+f_8 y^2+f_9 x+f_{11}x^2)
\]
for free coefficients, where constant radial terms can be
renamed harmlessly. If \(t_0\ne0\), the order-three compatibility
equations \(4f_5t_0,\ -2(10f_5t_1+f_9t_0),\ -2f_{11}t_0\)
force \(f_5,f_9,f_{11}=0\). If \(t_0=0,\ t_1\ne0\),
the same order-three equation gives \(f_5=0\), and the
order-four coefficients now give \(40 f_{11}t_1^2=0\),
then \(40 f_9t_1^2=0\). The negative modes vanish by
interchanging \(x,y\), so \(H_2\) is radial.
Then \(H_3\) has nonradial part \((g_5 y+g_6 x)V\)
by the cofactor and bounds. The order-four coefficient
\(-2g_6t_0\), or the order-five coefficient \(40g_6t_1^2\)
when \(t_0=0\), eliminates the positive mode and symmetry
eliminates the negative mode. The cofactor equation now makes
\(H_4\) radial: \(P_1,Q_1\) have degree one and \(R_0\)
is constant. This would make all components depend only on
\(q,z\) and contradict the Keller determinant.

## Reduction to a graph coordinate

Now \(t_0=t_1=0\). The quadratic part \((Ax^2+Cy^2)V\)
of \(H_1\) contributes zero to \(J_2\), and \(H_2\) has
the same form as above. If \(A\ne0\), order three gives
\(5A f_5=0\), and order four then gives \(-5A^2 f_4=0\).
By symmetry this holds also if \(C\ne0\).
When \(A=C=0\), order four gives
\(f_1 f_5=0,\ 5 f_{11} f_5=0,\ -f_1 f_{11}-10 f_5^2=0\);
these imply \(f_5=0\), and symmetry gives \(f_4=0\).
The replay then gives the following first entries:
\[
            (H_2)_1=f_0+f_1 r,\qquad
            (H_3)_1=g_0+g_1 r+ f_1(Ax^2+Cy^2),
\]
with a renaming of the constant terms. The lower part of \(P\)
has ordinary degree at most one. Consequently
\[
     P=F(z)+(f_1 z+g_1)xy+f_1(Ax^2+Cy^2)+\ell(x,y)
\]
where \(F\) is monic of degree five and \(\ell\) is affine linear.
If \(f_1\ne0\), the discriminant of \(P\) in \(x,y\)
has degree two in \(z\). By the
[quadratic-fibre argument](KELLER_QUARTIC_LINE_CONIC_STATIONARY_FIBRES.md#the-stationary-value-identity)
this is a nonzero square with a polynomial square root \(S\), and
\(1/S\) must be a rational derivative. But \(S\) is linear, so
its reciprocal has a nonzero residue. Hence \(f_1=0\).
If \(g_1\ne0\), the stationary value of \(P\) differs from \(F(z)\)
by a constant, while the stationary value identity forces a
nonconstant polynomial in it to have constant nonzero derivative.
This is again impossible. Thus \(P=F(z)+\ell(x,y)\),
and the submersion condition makes a binary coefficient of
\(\ell\) nonzero. Completing \(P\) by \(z\) and an independent
binary linear coordinate gives an inverse of degree at most five.
The [bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
proves invertibility, and inverse descent gives the same conclusion
over the original field.

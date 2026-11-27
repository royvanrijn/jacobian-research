# The two nonzero-parameter exceptional double-conic pencils

**Theorem QuinticKellerRankTwoDoubleConicExceptionalPencils.** Under the
hypotheses of the
[generic conic-argument exclusion](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md),
there is no Keller completion with either
\(\delta=-5\nu^2/9,\ \nu\ne0\) or \(\delta=0,\ \nu\ne0\).
Combined with that exclusion, the leading form
\[
                    (z(q^2+\delta z^4),\ z^3q+\nu z^5,\ q^2)
\]
can therefore occur only at \(\delta=\nu=0\).
The statements are valid in characteristic zero by the same coefficient
field embedding as in the generic case.

## The two linear kernels

Use \(r=q(x,y,1)\) and the coefficient basis \(x^i y^j\),
\(i=0,1\), over \(\mathbb C[r]\) from the generic proof.
For the first case rescale \(q,P,Q,R\) by nonzero constants
to normalize \(\nu=3,\ \delta=-5\).
The homogeneous linear operator in degree ten, after
dehomogenization and division by the common radial factor
\(r+5\), is
\[
                   6r p-2r(r-5)q_*-5(r+1)t .
 \tag{1}
\]
Here \(p,q_*,t\) stand for the dehomogenized coefficients from
\(P_4,Q_4,R_3\), respectively.
For positive basis degree two the bounds are affine, affine, constant.
The constant coefficient in (1) forces \(t=0\).
The remaining equation gives \(p=(r-5)q_*/3\),
so \(q_*\) is constant and the kernel is spanned by
\((r-5,3,0)\). At basis degrees at least three the smaller bounds
give zero kernel.

For the second case the common radial factor is \(r\), and the
divided operator is
\[
                    2(3r+5\nu)p-2r^2 q_*-5(r+2\nu)t .
 \tag{2}
\]
On the same affine, affine, constant bounds, write
\(p=a_1r+a_0,\ q_*=b_1 r+b_0,\ t=c_0\).
The coefficients of \(r^3,r^2,r^0,r^1\) successively give
\(b_1=0,\ b_0=3a_1,\ a_0=c_0,\ c_0=-10\nu a_1\).
Thus the kernel is spanned by \((r-10\nu,3,-10\nu)\).
Again the smaller bounds at higher basis degrees give zero kernel.

In both cases at basis degree one the third coefficient may
also be affine. Subtract a multiple of \((2r,1,2r)\), the
common \(q\)-translation solution, to remove its linear term.
The rest then lies in the kernel just computed. As in the
generic proof, translate the source to remove the
\(q\)-translation terms. The degree-one effects of the
translation on \(z\) are already forms in \(q,z\).
The surviving homogeneous correction \((P_4,Q_4,R_3)\)
therefore has the form
\[
                           V U_2+A_1(q,z),
 \tag{3}
\]
where \(U_2\) is an ordinary homogeneous quadratic,
\(A_1(q,z)\) is a vector of forms in \(q,z\), and respectively
\[
     V=(q-5z^2,3z^2,0)\quad\text{or}\quad
     V=(q-10\nu z^2,3z^2,-10\nu z).
\]
These are the homogeneous lifts of the dehomogenized kernel
vectors, with the third entry of degree one when nonzero.
In the notation of (1)–(2), both divided operators also have
zero kernel on constant triples: this is immediate from
the displayed rows, using \(\nu\ne0\) in (2).

## The quadratic term on the common cofactor fibre

Write \(H=(P_5,Q_5,R_4)\) and
\(D=\{q,\cdot\}_{x,y}\).
In the rational source coordinates \((z,q,x)\), use
\(F_z,F_q\) for derivatives at fixed values of the other
two coordinates. Then the determinant in the original
variables is
\[
                            \det(F_z,F_q,DF),
\]
since the determinant of \((z,q,x)\) is \(-q_y\) and
\(D=-q_y\partial_x|_{z,q}\).
These derivatives lower homogeneous ordinary degree by
one, two and zero respectively, so mixed terms in
degree \(11-i\) use a total drop of \(i\) from the three
top degrees.

For degree nine the total drop is two. The single drop-two
terms are the leading cofactor operator on \(P_3,Q_3,R_2\).
The other terms use two drop-one columns from (3).
Since \(D(VU_2+A_1)=V DU_2\), the terms involving \(U_{2,z}\)
or \(U_{2,q}\) in the other column vanish by repetition
of the vector \(V\). The contribution quadratic in \(U_2\)
is therefore exactly
\[
             \chi\,U_2 DU_2,\qquad
             \chi=\det(V_z,H_q,V)+\det(H_z,V_q,V).
 \tag{4}
\]
All remaining contributions are radial coefficients times
\(DU_2\), with coefficients polynomial in \(q,z\).
After dehomogenization the degree-nine equation thus says
that the leading cofactor expression on \(P_3,Q_3,R_2\),
plus \(\Theta(r)U_2+(\chi/2)U_2^2\), belongs to
\(\mathbb C[r]\), where the coefficients and \(U_2\) are
evaluated at \(z=1\). This uses \(D(\chi)=0\) and the kernel
criterion of the generic proof.

In the first case set \(r_*=-5\); in the second set \(r_*=0\).
All three undivided leading cofactor coefficients vanish
on \(r=r_*\). The coefficient \(\chi\) there is nonzero.
In the first case \(H_z=0\), \(H_q=(-10,1,-10)\),
\(V=(-10,3,0)\), \(V_z=(-10,6,0)\), giving \(\chi=300\).
In the second case \(H_z=(0,5\nu,0)\), \(H_q=(0,1,0)\),
\(V=(-10\nu,3,-10\nu)\), \(V_z=(-20\nu,6,-10\nu)\),
\(V_q=(1,0,0)\), giving \(\chi=150\nu^2\).
It follows that on the affine conic \(r=r_*\) a nonconstant
quadratic polynomial of \(U_2\) is constant.
Consequently \(U_2\) is constant there. This holds also if
the conic splits: it cannot be double or have parallel
lines, because either case would make \(q\) binary in
\(z\) and one other linear form. On two intersecting lines
the two constants must agree by evaluation at their
intersection. Since the conic is reduced and the dehomogenized
\(U_2\) has degree at most two, \(U_2\) is a polynomial in \(r\).
Thus all of (3) is in \(q,z\).

## The last two lower equations

With the drop-one terms now in \(q,z\), the degree-nine
equation is only the leading cofactor operator on \(P_3,Q_3,R_2\).
For positive basis degree the largest bounds are now affine,
affine, constant, and only at basis degree one. Thus the
same rows give a drop-two triple of the form
\(V U_1+A_2(q,z)\), with \(U_1\) binary linear after absorbing
its \(z\) term into \(A_2\).

In degree eight the possibilities for the drops are one
drop-three, or a drop-two with a drop-one. The latter terms
have exactly one possible normal derivative, from \(U_1\).
They give \(D(\Theta_1 U_1)\) for a radial coefficient.
Hence after dehomogenizing, the leading cofactor expression
on \(P_2,Q_2,R_1\) plus \(\Theta_1(r)U_1\) belongs to
\(\mathbb C[r]\). At positive basis degrees the coefficient
bounds for these new pieces are at most constant triples.
The divided leading operator is injective on such triples,
as checked in (1)–(2). Accordingly each degree-one basis
coefficient is its scalar coefficient in \(U_1\) times
the same inhomogeneous solution triple, while every higher
basis coefficient vanishes. In other words, the only
possible nonradial parts of \(P_2,Q_2,R_1\) are multiples
of \(U_1\) by radial coefficients.

In degree seven the total drop is four. The one drop-four
terms are the leading cofactors on \(P_1,Q_1,R_0\).
Terms with a drop-three and drop-one, or with a drop-two
and two drop-one columns, are linear in \(DU_1\) with
radial coefficients. Two drop-two columns have quadratic
part given by exactly (4) with \(U_1\) in place of \(U_2\);
terms using a radial part in either column are linear.
Thus reducing the dehomogenized bracket equation at \(r=r_*\)
once again leaves a quadratic polynomial of \(U_1\) with
coefficient \(\chi(r_*)/2\ne0\) constant on the conic.
The argument above forces \(U_1=0\) because it is binary linear.
The preceding injectivity removes the nonradial drop-three
pieces as well. The degree-seven equation is now the
leading cofactor operator on \(P_1,Q_1,R_0\), with zero
positive-degree kernel. Thus all components depend only
on \(q,z\), contradicting the Keller determinant.

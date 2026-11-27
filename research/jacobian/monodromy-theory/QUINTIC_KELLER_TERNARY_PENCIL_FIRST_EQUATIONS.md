# The first two determinant equations in the remaining independent quintic case

**Theorem QuinticKellerTernaryPencilFirstEquations.** Work over a
characteristic-zero field. Put \(p=z^2,\ q=x^2+yz\),
and fix the leading vector \(H_5=z(p^2,pq,q^2)\).
The degree-eleven determinant equation holds exactly for
quartic vectors of the form
\[
 H_4=(2pU+a q^2+b pq+c p^2,\ qU+pV+d p^2,\ 2qV+e p^2),
\]
where \(U,V\) are homogeneous quadratics and \(a,b,c,d,e\)
are scalars. Write \(u_i,v_i\) for the coefficients of \(U,V\)
in the order \(x^2,xy,xz,y^2,yz,z^2\), starting with index zero.
The degree-ten equation can then be solved for \(H_3\) exactly
when
\[
        u_1=u_3=v_3=0,\qquad
        a u_2=a v_1=a(u_0-u_4)=0 .
\]
For each such \(H_4\) the solution space for \(H_3\)
is affine of dimension ten. The later equations rule out
[a nonzero a](QUINTIC_KELLER_TERNARY_PENCIL_NONZERO_QUARTIC_EXCLUSION.md)
and force [v₁=0 and u₀=u₄](QUINTIC_KELLER_TERNARY_PENCIL_UNEQUAL_PARAMETERS.md)
when \(a=0\). The [lower determinant equations](QUINTIC_KELLER_TERNARY_PENCIL_TRANSVERSE_PARAMETER.md) then exclude the remaining cases.

This is the single independent leading form left by
[QuinticKellerLeadingPartResidual](QUINTIC_KELLER_LEADING_PART_RESIDUAL.md).
The equations here concern compatibility of the upper
homogeneous terms of the actual polynomial determinant.

## The quartic form

Set \(n=(q^2,-2pq,p^2)^{\mathsf T}\) and
\(D=-z\partial_x+2x\partial_y\).
The rows of the adjugate of \(JH_5\) are
\(5z^3 q^2 D,\ -10 z^5 q D,\ 5z^7 D\), interpreted as row
vectors of coefficients of the derivation. This follows
by taking the three cross products of the gradients of
\(z^5,z^3q,zq^2\). Thus the term from the vector \(H_j\)
and two leading columns is \(5z^3 D(n\cdot H_j)\).
The degree-eleven equation is \(D(n\cdot H_4)=0\).

In the polynomial ring the kernel of \(D\) is \(k[z,q]\).
Indeed after inverting \(z\), the coordinates \(z,q,x\)
identify the ring with \(k[z,z^{-1},q,x]\), where the derivation
is \(-z\partial_x\). In its invariant ring,
a negative power of \(z\) cannot occur in a polynomial of the
original variables. Otherwise reduction of \(z^m f\) modulo
\(z\), for the least \(m>0\) clearing the pole in the invariant
ring, would give a nonzero polynomial in \(q|_{z=0}=x^2\),
although \(z^m f\) vanishes there. This kernel computation is
the one used for the quartic square pencil in Wang and Yan,
[Section 3.4.1](https://arxiv.org/pdf/2610.36495v1).

By homogeneity, \(n\cdot H_4\) is a binary quartic in \(p,q\).
The five scalar parameters in the displayed \(H_4\) account
for its five coefficients. The kernel of contraction with
\(n\) consists exactly of \((2pU,qU+pV,2qV)\):
reduce modulo \(p\) and \(q\) and use coprimality to see that
\(p\) divides the first entry and \(q\) the third, then use
the equation for the middle entry. This proves the
claimed quartic form.

## The cubic compatibility

The degree-ten equation is the sum of the term
\(5z^3D(n\cdot H_3)\) and the three determinants using one
column of \(JH_5\) and two of \(JH_4\). In the basis of
degree-ten monomials this gives an explicit constant
linear operator on the thirty cubic coefficients and
an inhomogeneous column polynomial in the quartic parameters.

The [exact checker](verify_quintic_ternary_pencil_first_equations.py)
forms this operator over the rationals and verifies rank twenty.
Applying its left nullspace to the inhomogeneous column
gives the compatibility ideal. For necessity, the checker
verifies ideal membership of
\[
 u_1^3,\quad u_3^2,\quad v_3^2,\quad
 a u_2,\quad a^2v_1,\quad a(u_0-u_4)
\]
by polynomial division by a Gröbner basis. Their vanishing
over any field gives the stated equations. For sufficiency,
direct substitution into every left-nullspace equation kills
it on both branches, \(a=0\) and
\(a\ne0,\ u_2=v_1=0,\ u_0=u_4\), with \(u_1=u_3=v_3=0\)
on both. Thus the inhomogeneous column lies in the image
exactly under the stated conditions. The rank and
thirty-variable domain give the dimension of the cubic
solution space. The rational calculation applies over
any characteristic-zero extension.

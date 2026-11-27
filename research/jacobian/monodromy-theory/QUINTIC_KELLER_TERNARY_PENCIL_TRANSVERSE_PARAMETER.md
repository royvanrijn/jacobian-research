# Excluding polynomial completions of the ternary quintic pencil

**Theorem QuinticKellerTernaryPencilTransverseParameter.**
A characteristic-zero Keller completion of the leading form
\[
(z^5,z^3(x^2+yz),z(x^2+yz)^2)
\]
has \(u_2=0\) in the
[quartic parameterization](QUINTIC_KELLER_TERNARY_PENCIL_FIRST_EQUATIONS.md).
Here \(u_2\) is the coefficient of \(xz\) in \(U\). Together with
the [nonzero-\(a\) exclusion](QUINTIC_KELLER_TERNARY_PENCIL_NONZERO_QUARTIC_EXCLUSION.md)
and the [unequal-parameter exclusion](QUINTIC_KELLER_TERNARY_PENCIL_UNEQUAL_PARAMETERS.md),
this puts the normalized quartic in the form
\[
 H_4=(2p kq+b p q+c p^2,\ kq^2+p(vx^2+w yz)+d p^2,\ 
                  2q(vx^2+w yz)+e p^2).
\]

Work in the normalized equal-parameter sector \(a=0\),
\(U=kq+u_2xz,\ V=vx^2+w yz\), where \(k=u_0=u_4\).
If \(u_2\ne0\), use the source scaling
\((x,y,z)\mapsto(\lambda x,\lambda^2y,z)\) and the inverse
target diagonal weights \((1,\lambda^{-2},\lambda^{-4})\).
They preserve the degree-five form, preserve the normalized sector,
and send \(u_2\) to \(\lambda u_2\). Take \(u_2=1\).

The [exact checker](verify_quintic_ternary_transverse_parameter.py) uses
\(L(R)=5z^3D(n\cdot R)\) with \(D=-z\partial_x+2x\partial_y\)
and \(n=(q^2,-2pq,p^2)\). For determinant degrees ten, nine
and eight, the newly appearing term is \(L\) on a cubic, quadratic
or linear vector respectively. The other terms are obtained by
assigning the three determinant columns to the already known
homogeneous parts. The checker parameterizes the solutions of each
linear equation over the rationals, and retains its compatibility
equations by taking the left nullspace.

In degree nine, \(k(3b+10k)=0\) is a compatibility equation.
On substituting \(b=-10k/3\), another one gives \(k=0\), so \(k=0\)
in either case. With \(k=0\), the cubic parameters \(f_i\) of the
checker satisfy
\[
 b(5f_4+1)=0,\qquad f_1=1+bv-3bf_5 .
\]
If \(b\ne0\), substitute \(f_4=-1/5\) and the equation for \(f_1\).
The degree-eight compatibility includes \(-12b/25=0\), a
contradiction. Thus \(b=0\).

With \(k=b=0,\ f_1=1\), the degree-eight equations first force
\(f_4=0,\ f_8=4v(w-v)\) by square identities. The remaining equations
give \(c=d=e=0\) and \(g_1=2(f_5-v)\) in the quadratic parameterization;
the checker verifies these polynomial consequences by division
over \(\mathbb Q\). The resulting cubic, quadratic and linear vectors
contain every candidate in this normalized sector.

The final exact computation forms the coefficient equations of
determinant degrees seven, six and five and adds
\[
                              T\det JH_1-1=0 .
\]
The last equation imposes the nonzero constant term of the Keller
determinant. The checker clears rational denominators and computes
a standard basis of this ideal over \(\mathbb Q\) with Singular.
The remainder of \(1\) is zero. The system therefore has no solution
over any characteristic-zero extension. This rules out \(u_2\ne0\).

## The equal-parameter case with zero transverse coefficient

**Theorem QuinticKellerTernaryPencilCompletionExclusion.**
There is no characteristic-zero polynomial Keller map with degree-five
leading vector \((z^5,z^3q,zq^2)\), where \(q=x^2+yz\).

By the preceding results, the quartic in a candidate is in the displayed
form with \(U=kq\). Write \(\eta=w-v\),
\(\lambda=3b+10k\), \(M=f_4-2k\eta\), and \(N=f_8-4v\eta\).
The cubic parameterization in the same checker gives the
degree-nine compatibility equations
\[
                     -2\lambda M=0,\qquad 3\lambda f_5=0.
\]
If \(\lambda\ne0\), these yield \(M=f_5=0\). The degree-eight
compatibility then includes \(5N^2=0\) and \(-k\lambda f_9=0\).

If \(\lambda=0\), put \(b=-10k/3\) in the parameterization. The
degree-eight equations include \(5N^2=0\) and, with
\(S=-48ck-15f_1+220kv\),
\[
                    \frac{2MS}{15}=0,\qquad
                    \frac{2M(S-100M)}{5}=0 .
\]
They give \(M=N=0\). With these values the coefficient of \(x^7\)
in determinant degree seven is \(4k f_5^2\).
After \(f_5=0\), the coefficient of \(x^3z^4\) in the same degree
is \(5k f_9^2\). Thus for \(k\ne0\) we have \(f_5=f_9=0\)
in this case as well.

For \(k\ne0\) we now have \(M=N=f_5=f_9=0\) for arbitrary \(b\).
The checker substitutes these values and computes all solutions of the
degree-eight equation for \(H_1\). Their constant determinants
\(\det JH_1\) vanish identically.

For \(k=0\), the preceding equations already give \(f_4=0\) and
\(N=0\). If \(b\ne0\) they also give \(f_5=0\). If \(b=0\),
the degree-eight compatibility gives \(f_1 f_5=0\). Suppose
\(f_5\ne0\). Then \(f_1=0\), and the coefficient of \(x^6z\)
in determinant degree seven becomes \(-f_5g_1\), forcing \(g_1=0\).
The coefficient of \(x^6\) in determinant degree six is now
\(3f_5^3\), a contradiction. Thus \(f_5=0\) for all \(b\).

With \(k=0\) and \(f_4=N=f_5=0\), the checker parameterizes all
solutions of the degree-eight equation for \(H_1\). The determinant
equations of degrees seven and six, and \(T\det JH_1-1=0\), give a
rational ideal whose standard basis reduces \(1\) to zero. Thus
no constant-nonzero determinant is possible in the \(k=0\) case.

## Independent degree-five leading parts

**Theorem QuinticKellerIndependentLeadingExclusion.**
The three degree-five leading parts of a characteristic-zero polynomial
Keller map of total degree five in three variables are linearly dependent.

The [quintic leading-part reduction](QUINTIC_KELLER_LEADING_PART_RESIDUAL.md)
puts every independent leading vector into the form
\((z^5,z^3q,zq^2)\) over the algebraic closure by linear source and
target changes. These changes preserve the Keller condition. The
completion exclusion rules out this vector.

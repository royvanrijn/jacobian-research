# A cubic target component for three-variable quintic Keller maps

**Theorem QuinticKellerCubicTargetComponent.** Over a characteristic-zero
field, a three-variable polynomial map of ordinary degree at most five
with nonzero constant Jacobian is invertible if a nonconstant constant
linear combination of its coordinates has degree at most three.

The [degree-five inversion theorem](THREE_VARIABLE_KELLER_INVERSION_THROUGH_DEGREE_FIVE.md)
also covers the cases without such a component.

## The cubic reduction

We work first over \(\mathbb C\) and take the linear combination
as the third component \(R\). The [bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
covers degrees one and two. It also covers a cubic \(R\) that
extends to a polynomial automorphism of inverse degree at most
five. In either case the degree bound for the plane maps after
completion is at most 25.

For a cubic \(R\), use the argument of Proposition 3.16 in
[Wang–Yan, v1](https://arxiv.org/pdf/2610.36495v1), printed
pages 15–19. The statement there includes the quartic bound on
the other two components. The proof of this reduction uses
those components only to supply the following identities on
each smooth affine fibre \(S_t\) of \(R\):
\[
  dR\wedge\alpha_t=dx\wedge dy\wedge dz,\qquad
  \alpha_t=d\bigl(\kappa^{-1}P|_{S_t}\,d(Q|_{S_t})\bigr),
  \qquad \kappa=\det J(P,Q,R).
\]
Here \(dR\) is nowhere zero by the Keller condition. The first
identity denotes the residue volume form on the fibre. The
Jacobians give the second identity directly, regardless of the
degrees of \(P,Q\).

For the reduction, the boundary-period argument of Lemma 3.13
uses this exactness to exclude smooth and reduced nodal cubic
curves at infinity. The remaining leading cubics are binary,
cuspidal, or a conic with its tangent. Lemma 3.14 uses the same
exactness to exclude the respective exceptional affine fibre
forms by a residue on \(x=0\). The eliminations in Sections
3.6.1–3.6.4 use only the coefficients of \(R\) and the absence
of critical points. They exclude the cuspidal case, give the
explicit completion of inverse degree at most five in the
surviving tangent case, and give a completion or the scalar
\[
                        r=x+x^2y
\]
in the binary case. Thus the argument applies to a cubic
coordinate of the present degree-five map. The only case still
requiring exclusion is \(R=r\), up to affine source changes and
affine rescaling of \(R\), which preserve the degree bound.

## The exceptional scalar through degree five

**Proposition QuinticExceptionalCubicScalarExclusion.** If
\(P,Q\in\mathbb C[x,y,z]\) have ordinary degrees at most five, then
\(\det J(P,Q,x+x^2y)\) cannot be a nonzero constant.

Put \(B=\mathbb C[x,y]\) and
\(D=-x^2\partial_x+(1+2xy)\partial_y\). The determinant for the
exceptional scalar is \([P,Q]_D=D(P)Q_z-P_zD(Q)\).
Writing \(u=1/x\) gives \(\mathbb C(x,y)=\mathbb C(r,u)\)
and \(D=\partial_u\). In particular its rational kernel is
\(\mathbb C(r)\) and its polynomial kernel in \(B\) is
\(\mathbb C[r]\). For the latter assertion, a reduced rational
fraction \(A(r)/C(r)\) that is polynomial has \(C(r)\mid A(r)\),
where a Bézout identity for \(A,C\in\mathbb C[t]\) makes
\(A(r),C(r)\) coprime. The same argument shows that \(B\)
has no \(D\)-primitive of a nonzero constant. Otherwise
\(Dg=c\ne0\) would give \(g=c/x+A(r)/C(r)\), hence
\(C(r)\mid x\). A nonconstant \(C(r)\) has degree at least
three, so \(g-c/x\) would be polynomial, a contradiction.

We recall the two degree-independent restrictions in Wang–Yan's
Proposition 3.17. If \([P,Q]_D=\kappa\ne0\), neither \(P\)
nor \(Q\) can be independent of \(z\): the determinant would
express a unit as the product of a polynomial \(D\)-derivative
and a \(z\)-derivative. Neither leading \(z\)-coefficient can
be a nonzero constant. In fact, normalize such a \(P\) to be
monic of \(z\)-degree \(m\) and put
\(S=P^{-1/m}\in B((z^{-1}))\) with leading term \(z^{-1}\).
The chain rule gives
\[
 D(QP_zS)-\partial_z(QD(P)S)=-\kappa S.
\]
Taking \(z^{-1}\) coefficients gives a polynomial \(D\)-primitive
of \(-\kappa\). A pair affine in \(z\) is also impossible. Write
it as \(P=pz+a,\ Q=qz+b\), where \(p,q\) are nonzero. The
leading equation makes \(p/q\) a rational function of \(r\),
and \(p,q\) have no common zero since \(P_z,Q_z\) cannot both
vanish. Comparing coprime numerators and denominators gives
\(p,q\in\mathbb C[r]\). Then the constant equation is
\(D(qa-pb)=\kappa\), again impossible. These restrictions on
the determinant remain valid even when a target substitution
increases the ordinary degree.

Order the \(z\)-degrees as \(m\ge n\) and let \(p,q\) be the
leading \(z\)-coefficients. The leading determinant coefficient
gives
\[
                       p^n/q^m\in\mathbb C(r).
\tag{1}
\]
Both leading coefficients are nonconstant, so \(2\le m\le4\)
and \(1\le n\le m\). We use the following divisibility observation.
If \(f^i/g^j=A(r)/C(r)\) for coprime \(A,C\in\mathbb C[t]\),
any root \(c\) of \(A\) forces \(r-c\mid f\), and any root of
\(C\) forces \(r-c\mid g\). Indeed \(A(r),C(r)\) are coprime
and every polynomial \(r-c\) is squarefree of degree three,
since \(r_x,r_y\) have no common zero. This also applies
to \(c=0\), where \(r=x(1+xy)\).

If \(m=n\), (1) implies \(p/q\in\mathbb C(r)\), since
\(\mathbb C(r)\) is relatively algebraically closed in the rational
extension \(\mathbb C(r,u)\). If this ratio is nonconstant, one
of \(p,q\) contains a fibre factor \(r-c\). As \(m\ge2\) their
degrees are at most three, so that coefficient is a constant
multiple of \(r-c\). Equation (1) then forces the other coefficient
to lie in \(\mathbb C[r]\) as well. Both coefficients are linear
polynomials in \(r\). Nonproportional such coefficients have a
constant linear combination equal to a nonzero scalar, producing
a constant leading \(z\)-coefficient and contradicting the residue
restriction. Thus the ratio is constant. Cancel \(p\) by a constant
multiple of \(Q\) and reorder the components. This preserves the
ordinary degree bound and leaves unequal positive \(z\)-degrees.

Suppose \(m>n\) and \(m\ge3\). We have \(\deg p\le2\) and
\(\deg q\le4\), so \(p\) contains no fibre factor. If \(q\) does,
write \(q=(r-c)\ell\). This requires \(n\le2\) and
\(\deg\ell\le2-n\). Removing the fibre factor in (1) leaves
\(p^n/\ell^m\in\mathbb C(r)\), which must be constant, as neither
side contains a full fibre factor. For \(n=2\) this would make
\(p\) constant. For \(n=1\) nonconstant \(p\) requires
\(\deg p=m\deg\ell\ge3\), also impossible. If \(q\) has no fibre
factor, (1) is constant to begin with. Unique factorization then
gives \(\deg p\ge m/\gcd(m,n)\). This is at least three for
\(m=3\), and at least two for \(m=4\), where \(\deg p\le1\).
Both contradict \(\deg p\le5-m\).

It remains to take \(m=2,\ n=1\). Write
\(P=pz^2+cz+a,\ Q=qz+b\). If \(p/q^2\) is constant, subtract
the corresponding multiple of \(Q^2\) from \(P\) to get the
excluded affine pair. If it is nonconstant, one of \(p,q\)
contains a full fibre factor. In the case that \(p\) does,
\(\deg p\le3\) makes \(p\) a scalar multiple of that factor,
and the leading equation gives \(Dq=0\). Then the nonconstant
\(q\in\mathbb C[r]\) contains a full fibre factor as well.
In every nonconstant-ratio case we therefore have \(r-d\mid q\)
for some \(d\in\mathbb C\).

Reduce the bracket equation modulo \(r-d\). Since \(D(r)=0\),
both \(q\) and \(Dq\) vanish there, and the equation is
\[
                          -\overline{P_z}\,\overline{Db}=\kappa.
\]
On each irreducible component of this reduced fibre the equation
forces the coefficient of \(z\) in \(\overline{P_z}\) to vanish:
the other factor is independent of \(z\) and is a unit. Hence
\(r-d\mid p\). Now \(p\) is a nonzero scalar multiple of
\(r-d\) by its degree bound. The leading equation makes \(q\)
invariant under \(D\), so \(q\) is also a nonzero scalar multiple
of \(r-d\), using \(\deg q\le4\). Write \(p=\lambda q\).
The coefficient of \(z\) in the bracket gives
\(D(c-2\lambda b)=0\), so \(c=2\lambda b+L(r)\).
Finally the constant coefficient becomes
\[
                  D\bigl(qa-\lambda b^2-L(r)b\bigr)=\kappa.
\]
This contradicts the no-primitive restriction and completes the
exceptional scalar exclusion.

For a general characteristic-zero ground field, embed a finitely
generated subfield containing the coefficients of the map and
the chosen linear combination in \(\mathbb C\). The polynomial
inverse descends by uniqueness of the formal inverse after
translation to the origin, as in the bounded coordinate-completion
theorem.

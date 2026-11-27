# Two quintic spans for the tangent double-line scalar

**Proposition QuinticRankTwoTangentLineConicSpans.** Let \(F=(P,Q,R)\)
be a complex Keller map of ordinary degree five with
\[
             R=z+z^2(q+c),\qquad q=x^2+zy,\qquad p=z^2,
\]
where \(c\) is constant. Suppose the independent quintic leading
forms have the shape \(P_5=zU(p,q)\), \(Q_5=zV(p,q)\), with \(U,V\)
binary quadratics. Their span is one of
\[
                 \langle z^5,z^3q\rangle,\qquad
                 \langle z^3q,zq^2\rangle .
\]
The hypothesis on \(P_5,Q_5\) holds in the tangent irreducible
branch of the
[degree-two repeated-line pencil](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md),
using the
[tangent scalar normalization](KELLER_QUARTIC_LINE_CONIC_STATIONARY_FIBRES.md#tangent-conic).

## The linear equation for the first noninvariant terms

Use the derivations
\[
        \delta=z\partial_x-2x\partial_y,\qquad
        E=z\partial_z-x\partial_x-3y\partial_y,\qquad
        E_*=E-(z^{-2}+2cz^{-1})\partial_y .
\]
Both \(\delta,E_*\) annihilate \(R\), and the exact determinant
equation is \(z(E_*P\,\delta Q-\delta P\,E_*Q)=\det JF\).
Indeed the vector cross product of \(E_*\) and \(\delta\) is
\(z^{-1}\nabla R\).
After inverting \(z\), the coordinate \(u=x/z\) makes
\(\delta=\partial_u\) on \(\mathbb C[z,z^{-1},q,u]\).
Thus the polynomial kernel of \(\delta\) is \(\mathbb C[z,q]\).
To see the intersection with the original polynomial ring,
a polynomial in \(q,z,z^{-1}\) with a genuine \(z\)-denominator
would have reduced numerator at \(z=0\) a nonzero polynomial in
\(x^2\), which cannot cancel the denominator. Furthermore
\[
       \delta\mathbb C[x,y,z]\ \cap\ \mathbb C[z,q]
                               =z\mathbb C[z,q].
\tag{1}
\]
If \(\delta A=f(z,q)\), localizing gives \(A=xf/z+h(z,q)\).
The odd part under \(x\mapsto -x\) is \(xf/z\), so \(z\mid f\).
Conversely, \(\delta(xg)=zg\) for all \(g\in\mathbb C[z,q]\).

Define the invertible linear map \(L\) on binary quadratics
by multiplying the coefficients of \(p^2,pq,q^2\) by \(5,1,-3\),
respectively. Then \(E(P_5)=z(LU)(p,q)\) and likewise for \(Q_5\).
Let \(j\le4\) be the largest degree at which \(P_j,Q_j\) are
not both in \(\ker\delta\). Such a layer must exist by the
nonzero determinant. At degree \(5+j\) the bracket without
its exterior factor \(z\) gives
\[
                   (LU)\delta Q_j=(LV)\delta P_j .
\tag{2}
\]
Here \(E-E_*\) lowers homogeneous degree, and any bracket
between two invariants vanishes.

If \(LU,LV\) are coprime binary quadratics they remain coprime
after substitution of \(p,q\), because \(p,q\) are coprime.
Equation (2) either forces both derivatives to vanish by degree,
or, when \(j=4\), makes them a common scalar times \(LU,LV\).
The latter is also impossible by (1), since coprime substituted
forms cannot both be divisible by \(z\). Therefore \(LU,LV\)
have a common linear factor \(H(p,q)\). Independence rules out
a common quadratic. A constant target change makes \(LU=Hp\)
and \(LV=Hq\).

The equation becomes \(p\delta Q_j=q\delta P_j\).
On \(q=0\), after inverting \(z\), the class of \(Q_j\)
belongs to the kernel \(\mathbb C[z,z^{-1}]\) of \(\delta\).
Its ordinary homogeneous degree \(j\) makes this class a
multiple of \(z^j\). Since \(q=0\) is irreducible, the identity
holds without localizing on that surface.
Consequently there is a homogeneous \(W\) of degree \(j-2\) with
\[
                P_j=pW+I_j,\qquad Q_j=qW+J_j,\qquad
                I_j,J_j\in\mathbb C[z,q].
\tag{3}
\]
Indeed subtracting the invariant part of \(Q_j\) makes it
divisible by \(q\), and (2) then gives the assertion for \(P_j\).
Here \(\delta W\ne0\), so \(j=3\) or \(4\).

## Excluding a mixed linear factor

Suppose \(H\) is proportional to \(q-\lambda p\) with
\(\lambda\ne0\). Localize at \(z\) and reduce modulo this factor.
The resulting ring is \(\mathbb C[z,z^{-1},u]\), with
\(\delta=\partial_u\). In this ring \(p,q\) are both nonzero.
The class of \(W\) cannot be invariant: if it were, it would
equal a multiple of \(z^{j-2}\) by homogeneity, with any
difference divisible by the irreducible quadratic \(q-\lambda p\).
For \(j-2=1,2\) such a difference is zero or a constant
multiple of this invariant quadratic, contrary to \(\delta W\ne0\).

For \(j=4\), the bracket at degree eight integrates under
\(\delta\) to
\[
 zH(pQ_3-qP_3)+2pqW^2+
                     (qE(I_4)-pE(J_4))W\ \in\mathbb C[z,q].
\tag{4}
\]
In fact the first term comes from the quintic and cubic
layers, while \(E(p)=2p\) and \(E(q)=-2q\) give
\(\delta(2pqW^2+(qE(I_4)-pE(J_4))W)\) from the two quartic
terms. Corrections from \(E-E_*\) at degree eight use only
invariant terms. Reducing (4) modulo \(H\) makes a nonconstant
polynomial in \(u\) satisfy a quadratic with coefficients
in \(\mathbb C(z)\) and nonzero quadratic coefficient \(2pq\).
That is impossible by its highest \(u\)-degree.

For \(j=3\), the quartic terms \(P_4,Q_4\) are invariant.
The degree-seven equation is
\[
       zH(p\delta Q_2-q\delta P_2)+M\,\delta W=0,\qquad
                     M=qE(P_4)-pE(Q_4).
\]
The reduction of \(\delta W\) modulo \(H\) is nonzero, so \(H\mid M\).
Cancel \(H\). In the localized ring put \(m=\deg_u W\).
If either \(P_2\) or \(Q_2\) has \(u\)-degree larger than \(m\),
the equation after cancellation forces them to have the same
\(u\)-degree and ratio of their leading coefficients \(p/q\)
in the order \(P_2,Q_2\).
This cannot happen. For any ordinary homogeneous quadratic
of positive \(u\)-degree, its \(u\)-leading coefficient is a
constant times \(z^2\). Indeed insert \(x=zu,\ y=q/z-zu^2\)
and first replace any \(x^2\) by \(q-zy\).
In the resulting expressions the \(u\)-degrees of
\(z x,\ z y,\ xy,\ y^2\) are \(1,2,3,4\), with leading
coefficients \(z^2,-z^2,-z^2,z^2\), respectively, and all
other quadratic terms are invariant. It follows that \(P_2,Q_2\)
have \(u\)-degree at most \(m\).

Now take the degree-six equation, integrate under
\(\delta=\partial_u\), and reduce modulo \(H\). The contribution
from the quintic and linear layers is a multiple of \(H\)
and vanishes. From the two cubic terms, (3) gives \(2pqW^2\)
plus a linear term in \(W\), by the same calculation as in (4).
The quartic terms are invariant, so their bracket with the
quadratic terms integrates to
\(E(P_4)Q_2-E(Q_4)P_2\), of \(u\)-degree at most \(m\).
The only nonzero corrections from \(E-E_*\) at this degree
come from the quintic and cubic layers and integrate to
invariant multiples of the cubic terms, also of \(u\)-degree
at most \(m\). The integration constant is independent of \(u\).
On \(H=0\), the leading \(u\)-coefficient of the linear
polynomial \(W(x,y,z)\) after substitution is unchanged and
nonzero: its degree is one if its \(y\)-coefficient vanishes
and two otherwise. The term \(2pqW^2\) thus has degree \(2m\)
and cannot cancel. This proves that \(H\) is proportional to \(p\)
or to \(q\).

The inverse of \(L\) preserves each monomial line, so \(H=p\)
gives the span \(\langle z^5,z^3q\rangle\), and \(H=q\) gives
\(\langle z^3q,zq^2\rangle\). The proof transfers to characteristic
zero by embedding the finitely generated coefficient field.

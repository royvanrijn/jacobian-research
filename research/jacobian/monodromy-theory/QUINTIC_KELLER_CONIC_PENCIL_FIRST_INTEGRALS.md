# First integrals of a quintic Veronese conic pencil

**Theorem QuinticKellerConicPencilFirstIntegrals.** Let \(p,q\) be coprime independent homogeneous quadrics in three variables over \(\mathbb C\), and let \(\ell\ne0\) be linear. Suppose the pencil spanned by \(p,q\) contains a nonsingular projective conic and does not contain \(\ell^2\). There is no polynomial map with nonzero constant Jacobian and ordinary homogeneous leading part of degree five
\[
                     H=(\ell p^2,\ \ell pq,\ \ell q^2).
\]
The characteristic-zero version follows by embedding a finitely generated coefficient field in \(\mathbb C\).

This is the quadratic-argument pattern \((r,m,e)=(1,2,2)\) from
Wang and Yan's [homogeneous factorization lemma, printed pages 6–7](https://arxiv.org/pdf/2610.36495v1).
Their [Lemma 3.8 and Proposition 3.9, printed pages 11–12](https://arxiv.org/pdf/2610.36495v1)
use differential constants on the generic conic for the quartic leading form \(N(p,q)\). The linear factor in \(H\) changes the degree condition for first integrals. The following argument uses the pole orders of \(p/\ell^2\) on the generic conic and the defining quadratic equation of the Veronese image.

## Homogeneous first integrals of degrees at most nine

We prove that a nonzero homogeneous polynomial \(R\) of positive degree \(d\le9\) with
\(\det J(H_1,H_2,R)=0\) has \(d=5\) and belongs to the span of the components of \(H\).

Take \(\ell\) as the third linear coordinate and write \(P(x,y)=p(x,y,1),\ Q(x,y)=q(x,y,1)\). These two polynomials are algebraically independent. Indeed, if either were constant, the pencil would contain \(\ell^2\). If both were nonconstant and dependent, Wang and Yan's [Proposition 2.2(iii), printed pages 5–6](https://arxiv.org/pdf/2610.36495v1) writes them as univariate polynomials in a common polynomial. If both outer degrees are one, an affine relation produces \(\ell^2\) in the pencil. Otherwise, since \(\deg P,\deg Q\le2\), the common polynomial is affine linear. Homogenizing would then make \(p,q\) depend on only two independent linear forms, one of them \(\ell\), forcing every conic in the pencil to be singular.

Put \(r=Q/P\), \(w=t^5P^2\), and let \(R_*(x,y)=R(x,y,1)\). In the chart given by source coordinates \((tx,ty,t)\), the three functions in the determinant are \(w,rw,t^dR_*\). Its vanishing is equivalent to
\[
                       dw\wedge dr\wedge d(t^dR_*)=0.
\]
Thus the degree-zero ratio \((t^dR_*)^5/w^d=R_*^5/P^{2d}\) also has differential dependent on \(dw,dr\). The ratio and \(r\) belong to \(\mathbb C(x,y)\), while \(w_t\ne0\). It follows that
\[
                       d(R_*^5/P^{2d})\wedge dr=0
\]
in the binary function field. In characteristic zero this implies that \(R_*^5/P^{2d}\) is algebraic over \(\mathbb C(r)\).

For an indeterminate \(\lambda\), the generic projective conic \(q-\lambda p=0\) over \(\mathbb C(\lambda)\) is smooth by the nonsingular-member hypothesis. A common projective zero of \(p,q\) exists by Bézout and gives a rational point over \(\mathbb C(\lambda)\), so the generic conic is rational over that field. Its function field is identified with \(\mathbb C(x,y)\) over \(\mathbb C(r)\) by substituting \(\lambda=r\) in the chart \(\ell=1\). Therefore
\[
                        R_*^5/P^{2d}\in\mathbb C(r)^\times .
\]
On that conic \(P=p/\ell^2\) is nonconstant by the independence of \(P,Q\). It has a pole after extending constants to an algebraic closure. Its pole orders are at most four, since \(\ell=0\) meets a smooth conic with total multiplicity two and the denominator is \(\ell^2\). At a pole of order \(m\in\{1,2,3,4\}\), the valuation equality from the last display gives \(5\mid 2dm\). Hence \(5\mid d\), and \(d=5\) in the required range.

For \(d=5\), the same differential argument applied directly to \(t^5R_*/w\) gives \(R_*/P^2\in\mathbb C(r)\). Back in homogeneous coordinates this says
\[
                          R=\ell p^2 a(q/p),\qquad a\in\mathbb C(T).
\]
The rational function \(a\) cannot have a finite pole at a value \(\alpha\). Each irreducible factor of \(q-\alpha p\) is coprime to \(p\), so polynomiality of \(R\) would force \(q-\alpha p\) to divide \(\ell\), contradicting degrees. Thus \(a\) is a polynomial. If its degree were \(k>2\), the same coprimality at the divisors of \(p\) would force \(p^{k-2}\) to divide \(\ell\). Hence \(k\le2\), giving \(R\) in the asserted span.

## The Veronese equation

Suppose \(F\) had leading part \(H\) and \(\det JF=\kappa\ne0\). Define
\[
                         G=F_2^2-F_1F_3 ,
\]
so \(\deg G\le9\) and
\[
                         \det J(F_1,F_2,G)=-\kappa F_1 .
 \tag{1}
\]
In particular \(G\) is nonconstant. If \(G_d\) is its leading homogeneous term, the part of degree \(d+7\) of the determinant on the left of (1) is \(\det J(H_1,H_2,G_d)\). It is zero because the right side has degree five. The first-integral argument forces \(d=5\) and \(G_d\) to be a linear combination of the components of \(H\). Subtract the corresponding linear combination of \(F\) from \(G\), yielding \(G'\) of degree at most four. The right side of (1) changes only by a constant, so \(G'\) is also nonconstant and its positive-degree leading term must satisfy the same first-integral condition. This is impossible at degree at most four, proving the theorem.

## Scope of the remaining pencil types

The case of a pencil with a nonsingular member and \(\ell^2\) in the pencil remains outside this theorem. If every member is singular, the coprime pencil is binary. To see this, choose a rank-two member \(xy\) and write the other as
\(a x^2+bxy+c y^2+d xz+e yz+f z^2\) in a linear frame. A rank-two member exists because independent rank-one quadrics are squares of independent linears and their sum has rank two. The determinant of the pencil matrix gives \(f=de=0\) and \(ae^2+cd^2=0\), as in Wang and Yan's [Lemma 3.7, printed pages 10–11](https://arxiv.org/pdf/2610.36495v1). If \(d\ne0\), then \(c=e=0\) and both members are divisible by \(x\). If \(e\ne0\), they are both divisible by \(y\). Thus \(d=e=0\), giving two quadrics in the same two linear forms. The [independent binary theorem](QUINTIC_KELLER_INDEPENDENT_BINARY_LEADING_EXCLUSION.md) handles the case of a binary \(\ell\); the [external binary pencil theorem](QUINTIC_KELLER_EXTERNAL_BINARY_CONIC_PENCIL_EXCLUSION.md) handles the independent linear factor.

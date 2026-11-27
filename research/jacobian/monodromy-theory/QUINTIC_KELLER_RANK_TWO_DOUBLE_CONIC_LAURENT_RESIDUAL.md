# A Laurent fibre restriction in the mixed double-conic pencil

**Proposition QuinticKellerRankTwoDoubleConicLaurentResidual.** Work over
\(\mathbb C\). Suppose \(F=(P,Q,R)\) is a three-variable Keller map with
\(R_4=q^2\) for \(q=xy+\gamma z^2\), \(\gamma\ne0\). Assume that the
quintic span has a basis with
\[
             P_5=a zq^2+a_1 z^3q+a_0z^5,\qquad
             Q_5=b z^3q+b_0z^5,\qquad ab\ne0 .
\]
Then \(b_0+b\gamma=0\). After affine source transformations and rescaling
the second target component, \(Q\) has the form
\[
             Q=s^3xy+A(s)x+B(s)y+C(s),\qquad s=z-z_0,
 \tag{1}
\]
where \(A,B\) have degrees at most two and \(C\) has degree at most one.
After possibly interchanging \(x,y\), there are \(k\in\{1,2\}\),
\(v\ne0\), \(u\in\mathbb C\), and a polynomial \(\widetilde B\)
with \(\widetilde B(0)\ne0\) such that \(A(0)\ne0\) and
\[
 \begin{split}
 B(s)&=s^{3-k}\widetilde B(s),\qquad \deg\widetilde B\le k-1,\\
 s^k(Q-u)-v&=(s^k x+\widetilde B(s))(s^3 y+A(s)).
 \end{split}
 \tag{2}
\]
These are necessary forms. The proposition does not assert or exclude
their completion to Keller maps. With geometric factorization the
restriction extends to characteristic zero by coefficient-field embedding.

This concerns part of the degree-two pencil alternative for a double
conic. The degree-four pencil has a separate
[determinant exclusion](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_QUARTIC_PENCIL.md).

## The slice in the binary variables

Give \(x,y,z\) weights \(1,1,0\). Put \(w=xy,\ \rho=w^2\).
The first and third leading weight-four components are
\(a z\rho+P_4(x,y,0)\) and \(\rho\). Since \(Q_5\) has no
weight-four term, \(Q\) has weight-four part \(Q_4(x,y,0)\).
The top determinant coefficient forces
\(\{Q_4(x,y,0),\rho\}=0\), with the binary bracket in \(x,y\).
Thus \(Q_4(x,y,0)\) is proportional to \(\rho\) by homogeneity.
Subtract that multiple of \(R\) from \(Q\).

For the next weight, write the weight-three part of \(Q\) as \(zA_3+B_3\),
where \(A_3,B_3\) are binary cubics. The top determinant equation is
\[
        a\rho\{zA_3+B_3,\rho\}-A_3\{P_4(x,y,0),\rho\}=0 .
\]
Its \(z\) coefficient and then its constant coefficient give
\(A_3=B_3=0\). Indeed the polynomial kernel of the binary bracket
with \(\rho\) is \(\mathbb C[xy]\): in coordinates \(x,xy\) its
rational kernel is \(\mathbb C(xy)\), and a reduced univariate
rational fraction in \(xy\) is polynomial only if its denominator
is constant. This kernel has no nonzero binary cubic.

Let \(V\) denote the weight-two part of \(Q\), and put
\(T=P_4(x,y,0)/(a\rho)\). The top determinant equation now reads
\[
             \{V,\rho\}=\{T,\rho\}\,\partial_z V .
\]
With \(z+T\) as variable and \(x,y\) as coefficients, it says
the coefficients of \(V\) in \(z+T\) lie in \(\mathbb C(xy)\).
Binary homogeneity of degree two and weight zero of \(T\)
makes each such coefficient a constant multiple of \(xy\).
Therefore \(V=xy\,\Phi(z+T)\) for a cubic \(\Phi\), whose leading
coefficient is \(b\ne0\) from \(Q_5\).
If \(T\) had a pole along \(x=0\) or \(y=0\), then
\(xy\,\Phi(z+T)\) would have a pole there of order at least two.
Since those are the only possible poles of \(T\), it is polynomial
and homogeneous of degree zero, hence constant. We get
\[
                      Q=\Phi(z)xy+A(z)x+B(z)y+C(z).
\]
The degree bound and the specified quintic part give
\(\deg A,\deg B\le3\).

## A polynomial relation between the fibre parameters

On \(\Phi\ne0\), use \(X=\Phi x+B,\ Y=y+A/\Phi\). Then
\[
                       Q=XY+K(z),\qquad K=C-AB/\Phi.
\]
The change of variables to \((X,Y,z)\) has determinant \(\Phi\).
Thus the exact fibre form on \(Q=t\) is, up to sign,
\[
                               \frac{dX\wedge dz}{\Phi(z)X}
\]
on \(X\ne0,\ \Phi\ne0\), where \(Y=(t-K)/X\).
It has a regular primitive on the whole fibre from the two
other Keller components, by the
[exact residue identity of Wang and Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1).
For each root of \(\Phi\), take the product of a loop around that
root in the \(z\)-plane with a circle around \(X=0\).
Its period is a nonzero constant multiple of the corresponding
residue of \(dz/\Phi\), and is zero by exactness.
A simple root would have nonzero residue. Since \(\Phi\) is cubic,
it therefore has one triple root. Normalize it to \(s^3\) with
\(s=z-z_0\), by scaling \(Q\) and the corresponding fibre value.
The form \(dz/\Phi\) then has rational primitive \(L(z)=-1/(2s^2)\).

Choose the sign of the regular primitive \(\gamma_t\) so that
\(d\gamma_t=dX\wedge dz/(\Phi X)\). On the open fibre above,
\(\gamma_t\) has coefficients in
\(\mathbb C[z,\Phi^{-1},X,X^{-1},t]\) because the other two
Keller components are polynomials before substitution.
The closed form \(\gamma_t+L(z)dX/X\) has coefficient of \(X^{-1}\)
in its \(dX\) part independent of \(z\). Indeed the derivative
in \(z\) of this coefficient equals the \(X^{-1}\) coefficient
in the derivative in \(X\) of its \(dz\) part, which is zero.
Hence the coefficient is a polynomial \(H(t)\in\mathbb C[t]\).

For \(z\) with \(\Phi(z)\ne0\), set \(t=K(z)\).
The disc in \(X\) at fixed \(z\) with \(Y=0\) lies in the regular
fibre chart and includes \(X=0\). The regular primitive contributes
zero residue on this disc. It follows that \(H(K(z))=L(z)\)
as rational functions. This uses the actual polynomial primitive
to make \(H\) a polynomial.

Since \(L=-1/(2s^2)\), both \(K,H\) are nonconstant. For a polynomial
in \(K\), its pole set is that of \(K\), and pole multiplicities
are multiplied by the degree of the outer polynomial.
Consequently \(K\) has just one pole at \(s=0\), of order one or two.
If its order is two then \(H\) is linear and \(K\) is affine in
\(s^{-2}\). Thus in both cases \(K=u+v s^{-k}\), \(v\ne0\),
\(k\in\{1,2\}\). In the original normalization the degree-five
term \(C\) would give a pole at infinity of \(K\) if
\(b_0+b\gamma\ne0\), since \(AB/\Phi\) has growth at most cubic
there. Hence that degree-five coefficient vanishes.

Translate \(x,y\) by constants to remove their cubic coefficients
from \(A,B\), leaving \(K\) unchanged. The identity
\(AB=(C-u)s^3-vs^{3-k}\) shows that
\(\operatorname{ord}_{s=0}(AB)=3-k\). It also bounds
\(\deg C\le1\) since \(\deg A,\deg B\le2\).
When \(k=2\), one of \(A,B\) is a unit at zero by the
order-one product condition. When \(k=1\), if both vanished
they would both have simple zeros. Then \(Q_x=Q_y=0\)
at \(s=0\) and \(Q_s\) is a nonconstant affine polynomial
in \(x,y\), giving a critical point. Thus in either case we
may arrange \(A(0)\ne0\). Writing \(B=s^{3-k}\widetilde B\)
in the product formula \(Q=XY+K\) now gives (2).

## Exclusion when the conic meets the repeated pencil line transversely

**Corollary QuinticKellerRankTwoDoubleConicTransverseMixedExclusion.**
For the mixed conic-argument span above, a Keller completion is
impossible if the restriction of the irreducible conic \(q\) to
\(z=0\) has rank two. This holds in characteristic zero with geometric
factorization.

Diagonalize the rank-two binary restriction as \(xy\) and remove the
mixed \(xz,yz\) terms by adding multiples of \(z\) to \(x,y\).
Then \(q=xy+\gamma z^2\) with \(\gamma\ne0\) by irreducibility.
Scale the two quintic components and subtract a multiple of the
second from the first to get
\[
       P_5=z(q^2+\delta z^4),\qquad Q_5=z^3q+\nu z^5 .
\]
The proposition forces \(\nu=-\gamma\ne0\).
The [generic](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md)
and [exceptional](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_EXCEPTIONAL_PENCILS.md)
cofactor exclusions together require \(\delta=\nu=0\), a contradiction.

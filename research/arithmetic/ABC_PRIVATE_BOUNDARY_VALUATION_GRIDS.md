# Primitivity at the boundary of high-power components

A low-degree component cannot have a private affine zero of
a high-powered root. On primitive integer points, an arithmetic
Nullstellensatz then confines the prime support of that root.
This gives a determinant count with overlapping roots within
each letter and treats translated grids of exact valuations.

## PV1. Bounded releases with roots on the same letter

Fix \(m\ge2\). Take three nonzero signed terms
\[
       A_T \prod_{j\in I_T} x_j^{e_j},\qquad 1\le T\le3 ,
\]
with fixed integer coefficients and fixed positive integer exponents on active coordinates. Their sum is zero, the \(x_j\) are positive integers
at most \(H\), and there are at most \(m\) active coordinates.
Require that the supports of the moving variables on different
roles be disjoint. On each role designate a set of coordinates
of exponents at least \(E\). Every other active coordinate on
that role has support in the product of its designated roots.
There is no restriction on overlap within one role.

For \(E>4m\), the number of tuples is
\[
 O_{m,(e_j),E,\varepsilon}
       \left(H^{\,m^2(4m/E)^{1/(2(m-1))}+\varepsilon}\right).
                                                        \tag{PV1}
\]
The constant is uniform in the coefficients. With one active
coordinate the bound is constant, and with all active
coordinates on just one letter the divisor bound is available.
The exponent in (PV1) is explicit in \(m,E\); the other constants
are allowed to depend on \(m\) at each fixed budget.

The initial polynomial is nonzero and its geometric components
have dimension at most \(m-1\). In the latter one-letter case
any actual tuple bounds the fixed product by
\(H^{\sum_j e_j}\), so the coefficient uniformity of the divisor
bound follows. Otherwise work with the positive torus points on
its geometric components, using the following low-degree count.

### Private affine zeros force high degree

Let \(Z\) be a positive-dimensional geometrically integral
component or subvariety of projective degree \(d\), meeting the
open torus. Denote its affine part by \(Z_{\rm aff}\). Fix a
designated \(x_i\) of role \(T\) and put
\(P_{\ne T}=\prod_{j\notin I_T}x_j\).
Suppose \(Z_{\rm aff}\) contains a point with
\[
                     x_i=0,\qquad P_{\ne T}\ne0 .
\]
On the normalization of \(Z_{\rm aff}\), there is a prime
divisor with \(x_i\) vanishing and all other-role coordinates
units. Indeed this is the zero set of one nonzero regular
function on that open subset, and a component through the
chosen fibre has codimension one. All same-role exponents are
nonnegative. Thus the role-\(T\) monomial ratio \(F\) to either
other term, with signs chosen for \(1-F\), vanishes to order
at least \(E\) there.

Take a general linear section to a curve, as in the bounded
high-root release proof. A selected divisor meets this
section transversely on its smooth open part on the
normalization. One can avoid its intersections with the
other coordinate zero and pole divisors. The image curve
has degree \(d\), its normalization has genus at most
\((d-1)(d-2)/2\), and the orders of the three terms on the
selected branch have the asserted inequalities. The zeros
and poles of \(F\) and \(1-F\) lie in the union of the
coordinate hyperplane sections and the hyperplane at
infinity, with at most \(2md\) points in total. Riemann–Hurwitz
for \(F\) over \(0,1,\infty\) gives
\[
                        E\le\deg F\le 2g-2+2md<4m d^2 .
\]
Consequently, if \(d<D_0=\sqrt{E/(4m)}\), every designated
\(x_i\) has empty private affine zero set of this type.
[General sections and the three-value Riemann–Hurwitz argument](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

### A fixed integer supports every root on a low-degree component

Use the rational polynomial space of degree at most \(d\)
vanishing at all height-\(H\) integer tuples under consideration
on \(Z\). A maximal independent set of evaluation rows gives
an integral basis of this rational space, by minors of size
bounded solely in terms of \(m,d\). The coefficients are
\(H^{O_{m,d}(1)}\). Let \(Y\) be their common zero set.
It lies in \(Z_{\rm aff}\). For any projective point outside
\(Z\), choose a linear space through that point of complementary
dimension minus one disjoint from \(Z\), and project from a
hyperplane in that linear space avoiding the point. The image
of \(Z\) is a hypersurface of degree at most \(d\) and the
external point lies outside it. Pullback thus separates any
external point by a form of degree at most \(d\).
The complex polynomials of that degree vanishing on the
integer sample lie in the complexification of the rational
evaluation kernel. This also proves the containment for
components defined over a nonrational field.

For each designated \(x_i\), append \(x_i=0\) and
\(zP_{\ne T}-1=0\) to the equations of \(Y\), with new variable
\(z\). The system has no common complex zero. The arithmetic
Nullstellensatz supplies an integer combination of these
polynomials equal to a nonzero integer
\(B_i\) with \(|B_i|\le H^{O_{m,d}(1)}\). For an actual tuple
and a prime \(p\mid x_i\), primitivity across roles gives
\(p\nmid P_{\ne T}\). Evaluate the certificate modulo \(p\)
using \(z=P_{\ne T}^{-1}\) to get \(p\mid B_i\).
The given height bound now permits only \(H^{o(1)}\)
choices of each designated root by fixed-support enumeration.
The other roots use only the union of these same supports
on their role, so they also have \(H^{o(1)}\) choices.
No bound on the size of the original coefficients \(A_T\)
was used in producing the certificates.
[Krick–Pardo–Sombra, Theorem 1](https://www.maia.ub.edu/~sombra/publications/arith-nss/maino.pdf)

For \(d\ge D_0\) use the finite birational projection and
Salberger covering in the bounded-root release proof. At
dimension \(j\ge1\) the new components of lower dimension
are paid at \(H^{(j+1)/d^{1/j}+\varepsilon'}\), with their
degrees bounded at fixed parameters. The upper bound on
this exponent apart from \(\varepsilon'\) is
\(m D_0^{-1/(m-1)}\), and at most \(m\) steps occur.
A component with integer points contained in an intersection
of conjugates can either be passed through the same cover
or descended by the intersections, with bounded degrees.
At dimension zero the number of points is bounded by degree.
Choose the temporary epsilon smaller than the final allowance
divided by the bounded number of operations. This proves (PV1).
[Salberger, Theorem 0.10 and the projection argument](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

## PV2. Adaptive translated grids and growing numbers of grids

At height \(X\), divide the retained primes into disjoint
sectors on their original roles. In each sector allow an
arithmetic progression of valuations
\[
                         e=a+qt,\qquad t\ge0,
                      \qquad a,q\ge K .
\]
Both the choice of sectors and their orders may vary with the
row. On an ordinary divisibility sector one can just use
\(e=qt\) for positive \(t\). The exceptional prime powers
have combined radical at most \(X^f\). Use at most \(m\)
variables in total, charging a shifted sector two slots
and an ordinary divisibility sector one. Then at fixed
\(m,K\) the primitive count is
\[
           O_{m,K,\varepsilon}
               \left(X^{f+C_m K^{-1-1/(2(m-1))}+\varepsilon}\right).
                                                          \tag{PV2}
\]

For each shifted sector put its prime support in a squarefree
head \(R\) and its positive digits in an integer \(W\) supported
on \(R\). Its full product is \(R^a W^q\).
Choose the exceptional full coefficients via their radicals
at cost \(X^{f+o(1)}\). For a cutoff \(T_0\ge K\), a factor
raised to an order greater than \(T_0\) can be fixed at
\(X^{1/T_0+o(1)}\) cost, enumerating the exponent up to
\(\log_2 X\) when the root is nonunit. A fixed head fixes
the support of its digit, whose remaining numerical
possibilities cost a height subpower. Fixing a digit alone
preserves the high head. Drop unit variables. The remaining
exponents all belong to \([K,T_0]\) and the roots are bounded
by \(H=X^{1/K}\). In (PV1) take \(E=K\). Altogether the
exponent is at most
\[
              f+\frac{m}{T_0}
                 +\frac{m^2}{K}(4m/K)^{1/(2(m-1))}
                 +\varepsilon .
\]
The choice \(T_0=\lceil K^{1+1/(2(m-1))}\rceil\) gives (PV2).
A grid of step \(q\) and residue \(0\le r<q\) with quotient
at least one is represented by \(a=q+r\). Thus the count
applies also to shifts relatively prime to their steps.
It retains the exact original valuations.

For example one shifted grid per role uses \(m\le6\). Take
\(K=L^{-\theta}\) with fixed \(10/11<\theta<1\).
Then \(LK\to0\) and the loss in (PV2) is \(o(L)\).
Any fixed normalized exception cost below \(1/3\) now
gives a strict saving, without a purchased complementary
divisor. The same result allows a growing number of grids.
For
\[
       K=\left\lceil\frac{1}{L\log(1/L)}\right\rceil,\qquad
       m\log\log(1/L)=o(\log(1/L)),
\]
the explicit bound above with \(T_0=\lceil K^2\rceil\)
has both losses \(o(L)\). Indeed \(m/T_0=o(L)\) and the
negative power of \(K\) in the determinant term dominates
\(m^2\log(1/L)\). At every fixed \(L\), the values of
\(m,K,T_0\) are fixed before sending \(X\) to infinity.
Constants depending on this entire finite list are allowed.

The argument improves the bounded-grid residual estimate
for at least four disjoint ordinary root slots as well.
The three-single-root estimate of Browning–Verzobio is
stronger in its scope. A diffuse profile supported on prime
valuation indices can still avoid even these growing menus:
a grid of step at least \(K\) has \(O(1/(LK)+1)\) positions
in a fixed scaled band, while the prime degree indices
there number \(\asymp (L\log(1/L))^{-1}\).
Obtaining enough grid mass from the additive relation
and all three radical caps remains a separate extraction.
[Existing primitive Fermat and adaptive grid bounds](ABC_FINE_VALUATION_DIVISIBILITY.md)

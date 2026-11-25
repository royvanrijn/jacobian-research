# The joint curvature measure of ordinary binomial moments

For ordinary exponent rows, the coefficient arrays of every order have
one joint positive measure on the square. The measure determines the
multiset of unequal exponent rows. Its density on each curve also gives
an order comparison for inflation and identification of proportional
rows. At the end we give a collision in the relaxation that keeps the
joint positive measure and its support curves but drops the exact
factorial densities.

## A common measure for all orders

Let \(L(u,v)=\sum_i\log((a_i u+b_i v)!)\) for nonnegative
integers \(u,v\), and put
\[
 K(u,v)=L(u+2,v)-2L(u+1,v+1)+L(u,v+2).
\]
For an unequal row \(w=(a,b)\), set \(x=e^{-at}\), \(y=e^{-bt}\),
and push the measure
\[
 (x-y)^2\,\frac{e^{-t}}{t(1-e^{-t})}\,dt,\qquad t>0,
 \tag{1}
\]
to the point \((x,y)\). Let \(\mu\) be the sum over unequal rows.
The measure is finite because \(x-y=O(t)\) at zero and the
density decays exponentially at infinity.

**Theorem 1.** The whole coefficient curvature sequence is
\[
 K(u,v)=\int_{[0,1]^2}x^u y^v\,d\mu(x,y).                      \tag{2}
\]
In the notation \(P_n(z)=\sum_k\binom nk e^{L(n-k,k)}z^k\),
divide out the constant coefficient and scale the variable positively
to make both endpoints one. The resulting log multipliers are
\[
 J_{n,k}=\int
      \frac{x^{n-k}y^k-(1-k/n)x^n-(k/n)y^n}{(x-y)^2}\,d\mu.
 \tag{3}
\]
The logarithm of the ratio of the two normalization centers is
\[
 C=\int\frac{(y^{n+1}-x^{n+1})/(n+1)-(y^n-x^n)/n}{(x-y)^2}
                  \,d\mu .                                    \tag{4}
\]
The multiset of unequal rows is uniquely determined by the full
sequence \(K(u,v)\).

For the integral identity, integrating \(e^{-rt}\) first in
\(t>0\) and then in \(1\le r\le j\) gives
\(\log j=\int_0^\infty(e^{-t}-e^{-jt})\,dt/t\).
Summing for \(1\le j\le m\) gives
\[
 \log(m!)=\int_0^\infty
   \big(e^{-mt}-1+m(1-e^{-t})\big)\frac{e^{-t}}{t(1-e^{-t})}\,dt.
 \tag{5}
\]
The formula includes \(m=0\). All three linear combinations of
log factorials in (2)–(4) cancel the constant and the linear term
in their arguments, which proves those identities. In (4) the
quotient can have a first-order pole at the diagonal. The formula
integrates over the actual measures (1), and absolute convergence
follows directly from (5) before dividing by \((x-y)^2\).

In order \(n\), the three-term log difference at indices
\(k,k+1,k+2\) is \(K(n-2-k,k)\). Thus in the next order the two
neighboring log differences are the same integral with factors
\(x\) and \(y\). More generally, with forward differences
\(\Delta_u,\Delta_v\), the exact joint inequalities include
\[
 (-1)^{i+j}\Delta_u^i\Delta_v^j K(u,v)
       =\int x^u y^v(1-x)^i(1-y)^j\,d\mu\ge0 .
 \tag{6}
\]
Every matrix of integrals of products of real polynomials in
\(x,y\) is positive semidefinite. The orientation is retained by
splitting the measure into \(x>y\) and \(x<y\). On \(x>y\) the
shift in \(u\) reduces curvature less than the shift in \(v\);
on \(x<y\) the order is reversed.

For uniqueness, finite measures on the square agreeing on all
monomials agree on all continuous functions by polynomial
approximation, and hence coincide. Write each unequal row as
\(g(A,B)\) with \(g\) a positive integer and \((A,B)\) primitive.
The curves \((r^A,r^B)\), \(0<r<1\), of distinct primitive pairs
are disjoint, including the two private orientations. There is
no mass at the curve endpoints. With \(r=e^{-gt}\), the density
of \(\mu\) on a curve, with respect to \(dr\), is
\[
 \frac{(r^A-r^B)^2}{r(-\log r)}
      \sum_{g\ge1}\frac{m_g}{r^{-1/g}-1},
 \tag{7}
\]
where \(m_g\) is the multiplicity of \(g(A,B)\). Equality of
measures gives equality of the continuous densities in \(0<r<1\).
If two finite lists \(m_g\) differ, choose the largest \(g\)
with nonzero difference. After dividing by \(r^{1/g}\) and
letting \(r\downarrow0\), the difference of the sums in (7)
tends to that nonzero difference of multiplicities. This proves
uniqueness. Equal rows multiply a whole moment polynomial by a
single nonzero scalar.

## Inflation and identification along a primitive slope

**Theorem 2.** Fix an unequal primitive row \((A,B)\).
Inflating a row from \(g(A,B)\) to \(h(A,B)\), \(h>g\), or
identifying any two or more proportional rows into their sum,
increases the curvature measure of that block at every point
of its curve. Moreover, for the block by itself, write its
curvature at lattice point \((u,v)\) as \(K_{\rm old}(u,v)\)
and \(K_{\rm new}(u,v)\). Whenever \(A s+B t>0\) for a
nonnegative integer shift \((s,t)\), the operations give
\[
 \frac{K_{\rm new}(u+s,v+t)}{K_{\rm new}(u,v)}
       < \frac{K_{\rm old}(u+s,v+t)}{K_{\rm old}(u,v)} .
 \tag{8}
\]
For a whole profile containing the block, the normalized
multipliers \(e^{J_{n,k}}\) strictly decrease at every interior
index, simultaneously in all orders.

Write \(\phi_g(r)=(r^{-1/g}-1)^{-1}\) and \(z=-\log r>0\).
For \(h>g\), \(\phi_h>\phi_g\), and \(\phi_h/\phi_g\)
strictly decreases in \(r\), since the derivative in \(z\)
of its logarithm is
\[
 \frac{1/g}{1-e^{-z/g}}-\frac{1/h}{1-e^{-z/h}}>0 .
\]
This uses the strict increase of \(w/(1-e^{-w})\) for \(w>0\).
For a merger with \(G=\sum_j g_j\), the function \(\phi_g/g\)
is strictly increasing in \(g\) by the strict decrease of
\(w/(e^w-1)\). Hence \(\phi_G>\sum_j\phi_{g_j}\).
Also \(\sum_j\phi_{g_j}/\phi_G\) strictly increases with \(r\).
Thus in both operations the new density divided by the old one
is strictly decreasing in \(r\). To compare the two normalized
moments in (8), use the old probability law reweighted by
\(r^{Au+Bv}\) and its tilt by that decreasing density ratio.
The moment in the numerator of (8) adds the strictly increasing
factor \(r^{As+Bt}\). Its covariance with a strictly decreasing
factor is negative, by integrating the product of their signed
differences over two independent copies. This gives (8).
Finally, in (3) the numerator is strictly negative at every
interior index by weighted arithmetic-geometric mean. Increasing
just the block's measure therefore strictly decreases \(J_{n,k}\).

## A collision with the joint measure and the correct curves

The complete difference signs (6) and a common symmetric positive
measure on integer row curves still allow coprimality to fail.
Choose \(x=(\sqrt3-1)/2\), \(y=x^2\), so \(x+y=1/2\).
Use the two atoms \((x,y),(y,x)\) of mass \(\log4\) each.
These lie on the primitive curves for \((1,2)\) and \((2,1)\).
This measure gives \(K(0,0)=2\log4\) and
\(K(1,0)=K(0,1)=(\log4)/2\). Reconstruct the log multipliers by
their second differences and the two zero endpoints, or use (3).
In order two the interior multiplier is \(1/4\), and in order
three both interior multipliers are \(1/2\). The result is
\[
 1+\tfrac12 z+z^2,\qquad
 1+\tfrac32 z+\tfrac32 z^2+z^3
             =(1+z)(1+\tfrac12 z+z^2).
\]
The relative center in (4) is zero by symmetry, so the collision
persists with the linked centers. The atoms satisfy all the
moment positivity consequences of (2) but fail the densities
(7). Consequently a coefficient proof through this measure must
use additional inequalities from the actual densities. The row
operations above supply such inequalities for curvature, while
preservation of consecutive coprimality under those operations
remains a separate question.

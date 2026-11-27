# The five-cell standard-sector gap at every positive correlation

For the regular five-cell Gaussian simplex partition, the balanced
standard-sector interface form has a positive gap at each \(0<\rho<1\).
Together with [(S76)](SIMPLEX_INTERFACE_SYMMETRY.md#strict-supersolution-and-the-five-cell-reduction)
and the lift (S42), this gives a positive full interface gap off
geometric rotations for every ambient dimension \(d\ge4\), with the
gap constant independent of \(d\) for fixed \(\rho\). The proof uses
the [initial certificate](FIVE_CELL_STANDARD_GAP.md) through \(41/50\).
The comparison below gives a finite bridge from \(41/50\) to \(19/20\)
and an elementary estimate for higher correlations.

## An angular map to the \(A_3\) chamber

Use the scalar reduction (S49) on the full reference facet \(F\) in the
pair plane \(H\subset\mathbf1^\perp\subset\mathbb R^5\). Let
\(g_i=x_2-x_{i+2}\), \(1\le i\le3\). With respect to standard Gaussian
measure on \(H\), \(g\) has covariance \(C_0=I+J/2\), where \(J\) is the
three-by-three matrix of ones. Its Euclidean metric in the gap
coordinates is \(C_0^{-1}\), as follows from \(x_1=x_2=\sum_i g_i/5\)
and \(x_{i+2}=x_2-g_i\). Fold by the untied-coordinate permutations
to the chamber
\[
 F^+=\{g_1\ge g_2\ge g_3\ge0\}.
\]
Write \(\mu\) and \(\mu_+\) for the respective conditional Gaussian
probability measures on \(F\) and \(F^+\). Invariant functions have
the same ordinary mean and norm for these two probability measures.

For the reference chamber use \(H_*=\mathbf1^\perp\subset\mathbb R^4\),
with gaps \(g_i^*=y_1-y_{i+1}\) and the same ordering of the gaps.
This is one fundamental chamber for the permutation group \(S_4\).
Its gap metric is \((I+J)^{-1}\). The respective triples
\[
 \begin{split}
 U_0&=(1,1,0)\sqrt{5/6},&
 U_1&=(1,-1,0)/\sqrt2,&
 U_2&=(1,1,3)/\sqrt6,\\
 U_0^*&=(1,1,0),&
 U_1^*&=(1,-1,0)/\sqrt2,&
 U_2^*&=(1,1,2)/\sqrt2
 \end{split}
\]
are orthonormal in their gap metrics. On the sphere write
\[
 g=\cos\theta\,U_0+\sin\theta(\cos\varphi\,U_1+\sin\varphi\,U_2)
\]
and use the corresponding starred frame with angles \((\Theta,\psi)\)
in the target. These are polar coordinates at the right-angle vertex
of each spherical triangle. The first and last sorted gaps are
\(g_3=\sqrt{3/2}\sin\theta\sin\varphi\) and
\(g_1-g_2=\sqrt2\sin\theta\cos\varphi\). The middle gap is
\[
 g_2-g_3=\sqrt{5/6}(\cos\theta-d(\varphi)\sin\theta),
 \qquad d(\varphi)=a\cos\varphi+b\sin\varphi,\quad
 (a,b)=(\sqrt{3/5},2/\sqrt5).
\]
The target formulas are the same with the first-gap coefficient
\(\sqrt2\), the middle-gap coefficient \(1\), and
\(a_*=b_*=1/\sqrt2\). Write
\[
 W=1-\frac{d}{\sqrt{1+d^2}},\qquad
 W_*=1-\frac{d_*}{\sqrt{1+d_*^2}}.
\]
Thus \(0\le\varphi\le\pi/2\) and \(0\le u=1-\cos\theta\le W(\varphi)\).
The same description holds for \(\psi\) and \(u_*=1-\cos\Theta\).
Spherical area is \(du\,d\varphi\), so the triangle areas are
\(A_0=\int_0^{\pi/2}W\), \(A_*=\int_0^{\pi/2}W_*\).
Use \(S=A_*/A_0\). Define \(\psi=\psi(\varphi)\) by
\[
 \int_0^{\psi(\varphi)}W_*(s)\,ds=S\int_0^\varphi W(s)\,ds,\qquad
 c(\varphi)=\frac{W_*(\psi(\varphi))}{W(\varphi)},\qquad
 u_*=c(\varphi)u .
 \tag{A31}
\]
The azimuth map is a strictly increasing bijection with
\(\psi'=S/c\). Its primitive can be evaluated explicitly, since
\[
 \int_0^\varphi W(s)\,ds
 =\varphi-\arcsin\frac{a\sin\varphi-b\cos\varphi}{\sqrt{1+a^2+b^2}}
           -\arcsin\frac{b}{\sqrt{1+a^2+b^2}} .
\]
The Jacobian of (A31) in area coordinates is \(S\). Extend the map
radially without changing the radius, and call the extension \(T\).
It transports \(\mu_+\) to the conditional Gaussian probability
measure \(\mu_*\) on the \(A_3\) chamber.

We will use two properties of \(T\). It does not reduce angular
distances, and it does not reduce any of the three adjacent sorted
gaps. Both properties are homogeneous in the radius. We show them
on the unit sphere.

In orthonormal polar tangent frames the derivative of the map is
\(\left(\begin{smallmatrix}A&B\\0&C\end{smallmatrix}\right)\), with
\[
 X=A^2=\frac{c(2-u)}{2-cu},\qquad
 B=\frac{c'}{\sqrt{c(2-u)(2-cu)}},\qquad
 C^2=\frac{S^2}{X}.
 \tag{A32}
\]
For the derivative to expand lengths it suffices to check \(X>1\)
and \(B^2\le (X-1)(S^2/X-1)\), by the two principal minors of the
squared-length matrix minus the identity. In the checker, \(c>1\)
throughout each azimuth box. For a fixed azimuth,
\(X\) increases with \(u\), from \(c\) to \(c(2-W)/(2-W_*)\), and the
denominator in \(B^2\) decreases. The function
\((X-1)(S^2/X-1)=1+S^2-X-S^2/X\) is concave in \(X>0\).
The checker bounds \(B^2\) at the radial extreme and verifies the
inequality at the two enclosing endpoints for \(X\).
The inverse therefore contracts lengths on the target interior.
Geodesics between interior points of the reference triangle stay
inside it by spherical convexity. Applying the inverse to such a
geodesic, and then taking limits on the boundary, proves global
angular expansion.

Also
\[
 \frac{\sin\Theta}{\sin\theta}
 =\sqrt{\frac{c(2-cu)}{2-u}}
 \ \ge\ \sqrt{\frac{c(2-W_*)}{2-W}}\ \ge\ 1 .
 \tag{A33}
\]
The first inequality holds since the square of the ratio decreases
with \(u\). The second holds because \(\Theta\ge\theta\), with both
angles at most \(\pi/2\). For \(g_3^*\ge g_3\), a uniform bound
\(\psi'\ge t_0>0\) gives
\(\sin\psi\ge\min(1,t_0)\sin\varphi\) by integration and concavity of
sine. The checker verifies \(2\min(1,t_0)/\sqrt3>1\), and (A33)
then proves the first gap comparison.
For \(g_1^*-g_2^*\ge g_1-g_2\) the checker bounds the rightmost
radial factor in (A33) times \(\cos\psi/\cos\varphi\) from below.
On boxes near the azimuth right endpoint with the cosine vanishing it
uses \(\cos\psi/\cos\varphi\ge\min(1,t_{\rm suffix})\), where
\(t_{\rm suffix}\) is the lower bound for \(\psi'\) from that box
to \(\pi/2\). This follows by applying the same sine inequality to the
complementary angles.

For the middle gap, write \(u=qW\), \(q\in[0,1]\), so \(u_*=qW_*\).
After removal of the factor \(\sqrt{5/6}\) in the physical gap the
common function is
\[
 F(q,W)=1-qW-(1-W)r,\qquad
 r=\sqrt{\frac{q(2-qW)}{2-W}} .
\]
It is nonnegative by the chamber condition. For \(q>0\),
\[
 \partial_W F=-q+r-
     \frac{(1-W)r(1-q)}{(2-qW)(2-W)}\ge0 .
\]
Indeed \(r\ge\sqrt q\), so \((r-q)/r\ge1-\sqrt q\ge(1-q)/2\).
The fraction without \(r(1-q)\) is at most
\((1-W)/(2-W)^2\le1/4\).
At \(q=0\) the function is constant. Thus \(W_*\ge W\) and the larger
target coefficient prove the middle gap comparison.

The [checker](check_five_cell_A3_bridge.py) uses 256 adjacent rational
azimuth boxes, with coordinates in units of \(\pi/2\). It brackets the
exact monotone inverse in (A31) at every box endpoint on the 42-bit
dyadic grid using directed 192-bit Arb comparisons. Monotonicity gives
an enclosure on the complete box. It checks (A32) and the two outer
gap estimates just described on every box. These comparisons use
closed boxes and limiting endpoint values. In particular the minimum
enclosed last-gap ratio is greater than \(1.01475\), and the global
lower bound on \(\psi'\) is greater than \(0.92502\).

## Image comparison and Gaussian energy

Put \(\sigma=\sqrt{1-\rho^2}\),
\(p=\gamma_H(F)\), \(w=\nu_{12}(F)=\phi_0p\), and \(\ell=5w\).
By (FV1), \(p=\tfrac18+3\arcsin(1/3)/(4\pi)\).
The scaled local row is \(D_\rho=\sigma a_\rho/(\rho\ell)\).
On the full facet the kernel
\[
 B_\rho(x,y)=\frac{\sigma}{5}
             \left(2 M_\rho(x,y)+\sum_{l=3}^5 M_\rho(x,(2l)y)\right)
\]
has row \(D_\rho\) against \(\mu\), by (S49). Here \(M_\rho\) is the
ambient four-dimensional Mehler kernel. Thus
\[
 \frac{\rho\sigma}{w\ell}\,q_{\rm std,\rho}(h)
 =\frac12\iint B_\rho(x,y)(h(x)-h(y))^2\,d\mu(x)d\mu(y)
   -(1-\rho)\int D_\rho h^2\,d\mu.
 \tag{A34}
\]
On \(F^+\), the folded \(B_\rho\) is \(\frac16\sum_{\tau\in S_3}
 B_\rho(x,\tau y)\) against \(\mu_+\), with the same local row.

Since \(T\) preserves radii, the Gaussian exponential and prefactor
\(\sigma^{-3}\exp(-\rho^2(|x|^2+|y|^2)/(2\sigma^2))\) are common to
all image terms on the two sides, after taking out the factor
\(\sigma\) in \(B_\rho\). Write \(\xi=\rho/\sigma^2>0\). The native
exponents, with multiplicity, are \(x\cdot\tau y\) twice for each
\(\tau\in S_3\) and \(x\cdot(2l)\tau y\) once for each
\(3\le l\le5\), with the average over 30 terms.
The 6 self terms taken once together with the 18 adjacent terms are
exactly the group \(S_4\) acting on positions \(2,3,4,5\) and fixing
position 1. The remaining six terms are self terms. For every fixed
\(\tau\) the self dot dominates each adjacent dot, since their
difference is
\[
 (x_2-x_l)((\tau y)_2-(\tau y)_l)\ge0 .
\]
The average native exponential is therefore at least the average for
the 24-term subgroup.

Match each of these 24 terms to the permutation of the same ordered
positions of the target four-vector. To justify the dot comparison,
let \(u_i,v_i\) be the decreasing native coordinates on positions
\(2,3,4,5\), reindexed from 1 to 4. Use the convention that the
permuted dot product pairs \(u_k\) with \(v_{\pi(k)}\).
If \(\Delta_i(u)=u_i-u_{i+1}\), summation by parts gives the penalty:
\[
 x\cdot y-x\cdot \pi y
 =\sum_{i,j=1}^3
    \big[\min(i,j)-\#\{k\le i:\pi(k)\le j\}\big]
                           \Delta_i(u)\Delta_j(v).
 \tag{A35}
\]
The fixed coordinate 1 cancels. Every bracket is nonnegative. The
three gap comparisons show the same penalty is at least as large in
the target. The unpermuted native dot is at least the target dot by
angular expansion with the same radii. Hence every matched native dot
dominates its target dot. The resulting average Gaussian kernel
dominates the reference ordered-chamber kernel
\(\frac1{24}\sum_{\pi\in S_4}M^{(3)}_\rho(Tx,\pi Ty)\)
pointwise on \(F^+\times F^+\).

Against \(\mu_*\), that kernel is the ordinary Gaussian
Ornstein–Uhlenbeck operator on permutation-invariant extensions
from the chamber. Indeed the 24 chambers have equal Gaussian mass,
and partitioning the full Gaussian integral over them gives exactly
the averaged kernel. It has row one. For ordinary-mean-zero functions
the Hermite degree is at least two, since the fixed subspace of the
linear action on \(\mathbf1^\perp\) is zero. The degree-two invariant
space is spanned by
\(v=(|y|^2-3)/\sqrt6\), since a permutation-invariant quadratic form
on \(\mathbf1^\perp\) is scalar. All the other invariant degrees are at
least three. Pull back by \(T\), so \(v=(|x|^2-3)/\sqrt6\) on \(F^+\).
The square-difference energy in (A34), compared pointwise to this
reference, is at least
\[
 (1-\rho^3)\|h\|_{\mu_+}^2
       -(\rho^2-\rho^3)|\langle v,h\rangle_{\mu_+}|^2 .
 \tag{A36}
\]
Ordinary mean zero is preserved by \(T\).

## A scalar bound for the radial direction

Let \(L=2/(5p)\), so \(D_\rho\le L\) by (B7).
For \(\rho\ge19/20\) we have \(1+\rho>L\) by the explicit \(p\)
formula checked below. The degree-two lower bound
\((1-\rho^2)\|h\|^2\) on the reference energy already proves positivity
from (A34) there.

In the remaining interval \(41/50\le\rho\le19/20\), put
\(m=1+\rho+\rho^2-L>0\). Equations (A34)–(A36), after dividing by
\(1-\rho\), give the lower bound
\[
 \int W_\rho h^2\,d\mu_+
       -\rho^2|\langle v,h\rangle_{\mu_+}|^2,
 \qquad W_\rho=1+\rho+\rho^2-D_\rho\ \ge m .
\]
Cauchy–Schwarz applied to this single radial coefficient proves a
positive gap if
\[
       \rho^2\int v^2/W_\rho\,d\mu_+<1 .
 \tag{A37}
\]
The argument has used the ordinary Gaussian balance when invoking the
reference spectrum. The Cauchy–Schwarz bound is valid on that whole
balanced subspace.

For a useful upper bound on \(D_\rho\), use the cone derivative
(B7) as in (FV3). The probability \(D_\rho/L\) is at most the success
probability of the comparison with the smallest gap, \(g_3\).
Let \(s=\rho g_3/\sigma\ge0\). For \(N\) Rayleigh and independent of
a standard normal \(Z\), the marginal failure probability is
\[
 p_{\rm f}(s)=\Pr[Z\ge(N+\sqrt2\,s)/\sqrt3].
\]
Integration by parts in the Rayleigh integral gives
\(p_{\rm f}(0)=1/4\). Also
\(\mathbb E[Z1_{\{Z\ge N/\sqrt3\}}]=3\phi_0/4\),
by integrating first in \(Z\) and then using the Rayleigh density.
The Gaussian shift of the tail gives
\[
 p_{\rm f}(s)=\tfrac14 e^{-s^2/3}
 \mathbb E[e^{-\sqrt{2/3}\,s Z}\mid Z\ge N/\sqrt3]
 \ \ge\ \tfrac14 e^{-s^2/3-\sqrt{3/\pi}\,s}.
\]
The inequality is Jensen's inequality for the conditional expectation.
With \(2ab\le a^2+b^2\), the linear term in the exponent is at most
\((4/5)s^2+15/(16\pi)\). Thus
\[
 W_\rho\ge m+d\,e^{-k(\rho)g_3^2},\qquad
 d=\frac{L}{4}e^{-15/(16\pi)},\quad
 k(\rho)=\frac{17\,\rho^2}{15(1-\rho^2)} .
 \tag{A38}
\]

For \(x=e^{-k(\rho)g_3^2}\in[0,1]\) and \(m,d>0\), the secant
estimate with its quadratic correction is
\[
 \frac1{m+dx}\le
 \frac1m\left[1-\frac{d}{m+d}x-\frac{d^2}{(m+d)^2}x(1-x)\right] .
\]
Indeed the difference between the linear secant and the exact
reciprocal is \(d^2 x(1-x)/(m(m+d)(m+dx))\).
Put \(Z(t)=\int v^2 e^{-t g_3^2}d\mu_+\), so \(Z(0)=1\).
On a slab \([a,c]\), \(\rho^2/m(\rho)\) decreases because
its differentiated numerator has sign \(2+\rho-2L<0\).
The factor \(d/(m(\rho)+d)\) also decreases.
Moreover \(x_\rho\ge x_c\) and \(x_\rho(1-x_\rho)\ge
x_c(1-x_a)\) pointwise. Thus the left side of (A37) is at most
\[
 \frac{a^2}{m(a)}
 \left[1-\frac{d}{m(c)+d}Z(k(c))
 -\frac{d^2}{(m(c)+d)^2}
 \big(Z(k(c))-Z(k(c)+k(a))\big)\right].
 \tag{A39}
\]

These are exact tilted Gaussian polynomial moments. Use the ordered
coordinates \((g_1-g_2,g_2-g_3,g_3)\), called \(z\) here. Under the
unrestricted pair-plane Gaussian their covariance is
\(C=\left(\begin{smallmatrix}2&-1&0\\-1&2&-1\\0&-1&3/2\end{smallmatrix}\right)\),
and folding conditions on \(z\ge0\). In these coordinates
\[
 v^2=\frac{(z^\top C^{-1}z-3)^2}{6}.
\]
The tilt at \(t\) has covariance
\(C_t=(C^{-1}+2t e_3e_3^\top)^{-1}\); the moment integral of the
same quartic polynomial over the positive orthant with this
covariance is multiplied by
\(\sqrt{\det C_t/\det C}/\Pr_C(z\ge0)\).
The orthant formula (FV1) and integration by parts (FV7) evaluate
the moments by rational operations and balls.
The checker verifies \(Z(0)=1\) by enclosure and the expected
Gaussian mass \(\Pr_C(z\ge0)=p/6\) by enclosure as arithmetic
controls.

The contiguous exact slab endpoints are
\[
 .820,\ .825,\ .830,\ .835,\ .840,\ .850,\ .860,\ .870,\ .880,\
 .890,\ .900,\ .920,\ .940,\ .950 .
\]
The largest upper bound in (A39) is less than \(0.990793\), on
the first slab. The checker also verifies
\[
 m(41/50)(1-0.990793)>1/200,qquad 1+19/20-L>1/200.
\]
For the bridge, insert the Schur upper bound in (A37) and use
\(W_\rho\ge m(41/50)\). Above \(19/20\), use the degree-two
estimate and the second inequality. In the original normalization
this gives, on the whole range \(41/50\le\rho<1\),
\[
 q_{\rm std,\rho}(h)\ge
 \frac{\ell(1-\rho)}{200\rho\sqrt{1-\rho^2}}
          \|h\|_{L^2(\nu_{12})}^2 .
 \tag{A310}
\]
The initial certificate through \(41/50\) completes the gap
proof for all positive correlations.

If the initial global five-cell optimality interval ends below one
in the active space, the
[first-transition alternative](POSITIVE_CORRELATION_CONTINUATION.md#the-first-continuation-endpoint-has-a-tie-or-an-additional-zero-mode)
forces a nonsimplex tied maximizer there, because the gap excludes
an extra interface zero mode. The
[approximate reflection comparison](REFLECTION_SYMMETRY_OF_MAXIMIZERS.md#approximate-reflection-symmetries)
then forces a positive spanning-tree reflection defect on any
compact positive correlation band containing the endpoint.

Run the [checker](check_five_cell_A3_bridge.py) from the repository
root with Python 3.12.3 and python-flint 0.9.0:

~~~sh
.venv/bin/python research/plurality-stability/check_five_cell_A3_bridge.py
~~~

The [verification record](five_cell_A3_bridge_verification.json)
retains the internal replay and exact inputs. Independent mathematical
review and formal verification remain open.


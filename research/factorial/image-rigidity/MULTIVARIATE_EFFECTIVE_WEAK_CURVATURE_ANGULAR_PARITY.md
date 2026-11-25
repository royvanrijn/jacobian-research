# A weak-curvature band using even angular displacement

For the independent Gamma moments and midpoint scale \(\sigma\) in
the [weak-curvature proof](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md),
suppose \(N=n+1\ge3\) and
\[
                         \sigma\le(3\cdot10^7N^2)^{-1}.
                                                               \tag{1}
\]
Then both moment polynomials have simple roots and disjoint sets of
root arguments. In particular, for ordinary factorial moments
it is enough that \(S\le10^{-16}N^{-3}\).

The proof allows each individual complex root to move on the scale
of the same-order Hermite spacing. A holomorphic perturbation and
reflection give a smaller bound for its change in argument.

## A holomorphic interpolation

Fix either order \(j=n,N\). Use the real parameters
\(\rho_j\in[n/N,1]\), \(\sigma>0\), the centered logarithmic
cumulant function \(K_j\), and the real bounded \(h=\sigma Y_j\)
from the midpoint proof. Near zero write
\[
 K_j(z)=\frac{\rho_j z^2}{2}+e_{j,2}z^2+\sum_{k\ge3}e_{j,k}z^k.
\]
Its bounds give \(|e_{j,2}|\le\sigma^2/2\) and
\(|e_{j,k}|\le2\sigma^{k-2}/k\). For a complex parameter \(v\), put
\[
 K_{j,v}(z)=\frac{\rho_j z^2}2+v^2 e_{j,2}z^2
                     +\sum_{k\ge3}v^{k-2}e_{j,k}z^k,
\quad
 P_{j,v}(x)=\mathbb E_h\, j!\,[z^j]
          \exp\!\bigl(x(z+vh)+K_{j,v}(z+vh)\bigr).
\]
These expressions use the convergent series for
\(|v|\sigma|z+vh|<1\).
At \(v=0\) the function \(P_{j,v}\) is the Gaussian reference
polynomial \(G_j\). At \(v=1\) it is the centered exact moment \(F_j\),
by the homogeneous midpoint and uniform-tilt representation.

Write \(\alpha=|v|\sigma\). For \(|v|\sigma|z|\le1/2\),
\[
  |K_{j,v}(z)-\rho_j z^2/2|
                 \le \alpha^2|z|^2/2+(4/3)\alpha|z|^3.
                                                               \tag{2}
\]
Indeed the higher cumulants give the same geometric sum as in
the midpoint proof, with \(\sigma\) replaced by \(\alpha\).
Also \(|vh|\le \alpha N/2\). In particular if
\(\alpha\le(6\cdot10^6N^{3/2})^{-1}\), the midpoint relative
error estimate on \(|z|\le R=12\sqrt N,\ |x|\le3\sqrt N\),
with \(\sigma\) replaced by \(\alpha\), gives
\[
 \left|xvh+K_{j,v}(z+vh)-\rho_j z^2/2\right|
                   \le24\alpha R^3<1/2 .
\]
This estimate uses only absolute values in (2) and therefore
also holds for complex \(v\). The Taylor coefficients in \(z\)
of the exponential of this error minus one are bounded by
\(48\alpha R^3/R^k\) in degree \(k\). Uniformity in the bounded
real \(h\) ensures analyticity of \(P_{j,v}\) in \(v,x\)
throughout these domains.

## Individual root disks

For a Hermite root \(x_0\) of \(G_j\), use the disk radius
\(r=1/(80\sqrt N)\), one quarter of the same-order spacing.
On this circle put
\[
     \phi(t)=\sqrt N\,G_j(x_0+t/\sqrt N)/G_j'(x_0).
\]
It satisfies
\[
    \rho_j\phi''+(x_0/\sqrt N+t/N)\phi'-(j/N)\phi=0,\qquad
    \phi(0)=0,\quad\phi'(0)=1.
\]
On \(|t|\le1/80\) the two coefficients after solving for
\(\phi''\) are bounded by \(31/10\) and \(3/2\).
Integration along the straight segment from zero gives
\(|\phi|+|\phi'|\le e^{5|t|}<2\), by the integral
inequality for the resulting first-order system.
Then \(|\phi'(t)-1|\le(31/5)|t|\) and
\(|\phi(t)-t|\le(31/10)|t|^2\), by successive integration.
On the circle these imply \(|G_j'/G_j|<2/r=160\sqrt N\).
The differentiated Hermite equation now gives
\[
  |G_j^{(k)}/G_j|\le30(6\sqrt N)^k\qquad(0\le k\le j).
\]
The initial bounds follow from the logarithmic derivative and
the value one for \(k=0\). The normalized coefficients of the
inductive recurrence are bounded by \(3/4\) and \(1/24\).

As in the midpoint proof, extracting degree \(j\) in the
product with the Gaussian exponential expresses the relative
error as a sum of the preceding Taylor coefficient bounds
times the corresponding derivative ratios. After averaging
over \(h\), since \(R=12\sqrt N\), this gives on every circle
\[
 \bigl|P_{j,v}/G_j-1\bigr|
   \le48\cdot12^3\cdot31\,\alpha N^{3/2}<1/2
        \quad\text{for }\alpha\le(6\cdot10^6N^{3/2})^{-1}.
\]
Thus for every \(|v|\le L=(6\cdot10^6\sigma N^{3/2})^{-1}\)
there is one simple root \(x(v)\) in the disk about \(x_0\).
It is analytic throughout \(|v|<L\), by uniqueness and the
implicit function theorem.

## The even angular shift

The cumulants have real coefficients and satisfy
\(K_{j,-v}(-z)=K_{j,v}(z)\). The coefficient formula gives
\[
  P_{j,-v}(-x)=(-1)^jP_{j,v}(x),\qquad
  \overline{P_{j,v}(x)}=P_{j,\bar v}(\bar x).
\]
Since \(x_0\) is imaginary, the disk about it is invariant under
\(x\mapsto-\bar x\). Uniqueness of the root at real \(v\) gives
\(x(-v)=-\overline{x(v)}\). Therefore in the Taylor series for
\(x(v)-x_0\), the coefficients of odd degree are real and the
coefficients of even degree are imaginary. Cauchy's estimate
bounds each coefficient in degree \(k\) by \(r/L^k\).
Under (1) one has \(L\ge5\sqrt N\), so
\[
  |\Im x(1)-\Im x_0|
       \le\sum_{k\ge1}rL^{-2k}
       \le\frac1{1920N^{3/2}}.
\]
By the [quantitative Hermite gap](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_IMPROVEMENT.md#the-hermite-gap),
roots of the two reference orders have imaginary-coordinate
distance at least \(1/(200N^{3/2})\). The sum of the two
possible displacements is smaller than that gap. The physical
variable has argument given by the common scale \(\sigma\)
times the imaginary coordinate, up to the common base
argument and multiples of \(2\pi\). Here the entire ranges
of the disks are within one period since \(|x|<3\sqrt N\).
The exponential change of variable maps the \(j\) distinct
simple roots of \(F_j\) in these disks to \(j\) distinct
simple roots of the degree-\(j\) moment polynomial, so all
roots are accounted for.

For ordinary shapes,
\(\sigma^2\le W/n\le2S/n\). If \(S\le10^{-16}N^{-3}\),
then \(\sigma^2\le3\cdot10^{-16}N^{-4}\), which implies (1).

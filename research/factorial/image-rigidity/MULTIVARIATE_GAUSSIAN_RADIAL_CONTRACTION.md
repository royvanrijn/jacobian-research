# Log-radius contraction and physical center bounds for Gaussian factorial moments

The endpoint-normalized Gaussian deformation of a factorial
core contracts the full vector of log root radii. This gives
a physical center test for arbitrary diffuse tails, including
asymmetric cores and roots off the unit circle. A conjugate
pair uses the average of the two extreme radii rather than
the largest individual radius. The test has no simplicity
or contact-rank hypothesis.

## A majorization inequality for Gaussian deformation

Let \(P_0(z)=\sum_{k=0}^j a_k z^k\) have degree \(j\ge2\)
and nonzero constant coefficient, and put
\[
      P_W(z)=\sum_{k=0}^j a_k e^{-W k(j-k)/(2j)}z^k,\qquad W\ge0.
\]
Write \(x_1(W)\ge\cdots\ge x_j(W)\) for its log root radii,
counted with multiplicity, and let
\[
 \bar x=j^{-1}\log|a_0/a_j|,\quad
 D(W)=x_1(W)-x_j(W),\quad
 D^*(W)=2\operatorname{arsinh}\!\bigl(
                      e^{-W/4}\sinh(D(0)/2)\bigr).
\]
When \(D(0)>0\), set \(c(W)=D^*(W)/D(0)\).
When \(D(0)=0\) one can take \(c(W)=e^{-W/4}\).

**Theorem.** The mean \(\bar x\) is fixed, \(D(W)\le D^*(W)\),
and for every \(1\le k\le j\),
\[
      \sum_{i=1}^k(x_i(W)-\bar x)
                  \le c(W)\sum_{i=1}^k(x_i(0)-\bar x).         \tag{1}
\]
The corresponding lower partial sums satisfy the reverse
inequality. The diameter bound is attained in degree two,
for example by \(z^2+bz-1\) with \(b>0\).
The theorem allows complex coefficients and multiple roots.

In the common normalization, the polynomial equation is
\(\partial_W P_W=(D_z^2-jD_z)P_W/(2j)\),
where \(D_z=z\,d/dz\). This is the real contracting direction
of the multiplicative heat flow in
[Hall and Ho, Proposition 6.7 and Remark 6.8](https://link.springer.com/article/10.1007/s11005-025-01946-9).
At simple roots \(z_i=e^{x_i+i\theta_i}\), differentiating
the root equation and taking real parts gives
\[
 \dot x_i=-\frac1{2j}\sum_{l\ne i}
       \frac{\sinh(x_i-x_l)}
            {\cosh(x_i-x_l)-\cos(\theta_i-\theta_l)}.           \tag{2}
\]
In fact the complex log-root velocity is
\(-\frac1{2j}\sum_{l\ne i}(z_i+z_l)/(z_i-z_l)\);
this follows from the equation for \(P_W\) and
\(P_W''(z_i)/P_W'(z_i)=2\sum_{l\ne i}(z_i-z_l)^{-1}\).
The mean is fixed by the endpoint coefficients.

For \(d\ge0\), the nonnegative summand of (2) at
\(d=x_i-x_l\) is at least \(\tanh(d/2)\); a zero difference
contributes zero. Applying this at the extremal radii gives
\[
 \dot D\le-\frac1{2j}\sum_l
   \bigl(\tanh((x_1-x_l)/2)+\tanh((x_l-x_j)/2)\bigr)
       \le-\tfrac12\tanh(D/2).
\]
The last inequality uses
\(\tanh a+\tanh b\ge\tanh(a+b)\) for \(a,b\ge0\).
Integration of the inequality for \(\log\sinh(D/2)\)
gives \(D\le D^*\). If \(D\) vanishes it stays zero,
by the nonincrease of the diameter.

Next put \(S_k=\sum_{i=1}^k(x_i-\bar x)\).
The interactions within the top \(k\) cancel. Across
its boundary the differences are nonnegative, and
concavity of the hyperbolic tangent gives
\[
 \dot S_k\le-\frac1{2j}\frac{\tanh(D/2)}{D}
                               \sum_{i\le k<l}(x_i-x_l)
             =-\frac12\frac{\tanh(D/2)}{D} S_k .
\]
Here the cross sum is \(jS_k\). The function
\(\tanh(D/2)/D\) is decreasing in \(D>0\), and
\((\log D^*)'=-\tanh(D^*/2)/(2D^*)\).
This proves (1) by integration, including \(S_k=0\)
by the inequality for \(S_k/D^*\).
The log radii in sorted order and their sums are
continuous even at multiple roots. Away from such roots
the sorted sums are locally Lipschitz, and the estimates
apply almost everywhere, including ties. The coefficients
are polynomial in \(e^{-W/(2j)}\), and the discriminant
in this variable is nonzero at zero because the endpoint
binomial has simple roots. Thus only finitely many
multiple-root parameters require passage by continuity on
the real strength interval. The lower inequalities follow
from the fixed total sum. For \(z^2+bz-1\) the two radii
are reciprocal and their diameter is
\(2\operatorname{arsinh}(b e^{-W/4}/2)\), giving sharpness.

If only an upper bound \(\widehat D\ge D(0)\) is known,
one may replace \(D(0)\) by \(\widehat D\) in the factor
\(c(W)\) of (1). Indeed
\(\operatorname{arsinh}(e^{-W/4}\sinh y)\) is convex
for \(y\ge0\) and vanishes at zero, so its ratio to
\(y\) is nondecreasing. This allows direct use of the
[paired Gamma annuli](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md)
as initial data.

## The physical comparison for a conjugate pair

Fix adjacent orders \(n,N=n+1\) with \(n\ge2\), and a
finite ordinary core \(C\). In each order start from its
endpoint-normalized polynomial, with original rescaling
\(z=G_j^C t\). Thus \(\bar x_j=0\).
Define the nonnegative initial bounds
\[
   b_j=\tfrac12(x_{j,1}(0)+x_{j,2}(0)),\qquad
   a_j=-\tfrac12(x_{j,j-1}(0)+x_{j,j}(0)).
\]
One may instead use upper bounds on these quantities from
any initial annulus. Every nonreal root of the Gaussian
reference has log radius between \(-c_j(W)a_j\) and
\(c_j(W)b_j\). Indeed its conjugate has the same radius,
and (1) bounds the sum of the two largest or smallest
possible log radii. Real reference coefficients are used here.

For a diffuse tail use midpoint strength and bias
\[
 W=\sum_i d_i^2/m_i,\qquad \lambda=\sum_i d_i/m_i,\qquad
                   h=\sum_i d_i\in\mathbb Z .
\]
The positive physical Gaussian reference has relative
center
\[
 U(W,\lambda,h)=\log(G_N^C/G_n^C)
                           +h\log(N/n)-\lambda/(2nN).
\]
At a common physical root, its normalized log radii differ
by this amount, upper minus lower. The even reference is
strictly positive on the real axis at \(W>0\) by
[Gaussian averaging](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts).
Thus both roots in the comparison are nonreal.
A necessary condition for physical contact is
\[
 -c_N(W)a_N-c_n(W)b_n\ \le\ U(W,\lambda,h)\
                         \le c_N(W)b_N+c_n(W)a_n.             \tag{3}
\]
This improves any static radius bounds obtained by replacing
\(c_j(W)\) by one. It imposes no assumption on the angular
velocities.

For a bias interval given by two affine functions of \(W\),
the strict exclusions for the whole bias interval can be
checked at the two endpoints of a compact strength interval.
Indeed each \(c_j(W)\) is convex in \(W\), directly by
differentiating \(\operatorname{arsinh}(B e^{-W/4})\).
Thus the worst upper or lower radius bound after moving
the affine term for the bias to its side is attained at
an endpoint. In particular one can take \(-W\le\lambda\le W\)
and test every physical bias of ordinary diffuse rows at once.
The integer \(h\) stays fixed in each comparison.

## Exact ordinary transfer

Let \(\eta=\max_i |d_i|/m_i\) for the shared ordinary tail.
For bounded positive \(W\), the
[exact matched jet bounds](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#a-matched-expansion-with-a-second-order-ordinary-remainder)
give normalized coefficient errors \(O_{n,N}(\eta)\)
relative to the reference at \(W\), and relative center error
\(O_{n,N}(\eta)\) relative to \(U(W,\lambda,h)\).
The constants are uniform over the row counts, heights and gaps.
A strict violation of (3) on a compact strength and bias range
therefore gives coprimality of the exact ordinary pair for
an effective small enough \(\eta\), uniformly on that range.

For the use of the two-root averages, the transfer needs a
root matching with multiplicity. Both normalized polynomials
are monic and have endpoints one and coefficients bounded
by the binomial coefficients, so all their roots lie in a
fixed annulus. If their coefficient differences are \(O(\eta)\),
take disks of radius \(r=C_j\eta^{1/j}\) about the
reference roots, with a sufficiently large effective \(C_j\).
On the circle arcs bounding each component of their union
the reference polynomial has modulus at least \(r^j\).
The coefficient bound makes the error smaller than this.
Choose a slightly larger radius if necessary to avoid
tangencies. By Rouché each component contains the same
number of actual and reference roots, with multiplicities.
Each component has diameter \(O_j(r)\), so pairing roots
within each component changes log radii by \(O_j(\eta^{1/j})\).
An actual conjugate pair is therefore compared with two
reference radii even at a multiple or real reference root.
As the actual even moment is strictly positive on the real
axis, a physical common root necessarily enters this
two-root comparison. It follows that an exact common root
must satisfy (3) up to \(O_{n,N}(\eta^{1/N})\).
The absolute ordinary angular cutoff handles terminals above
its range when the test is applied to each remaining pair.

## A fully shared example at terminal fourteen

Use the one-row core \((10,11)\), \(n=13,N=14\),
and the tail parameters \((W,\lambda,h)=(10,6,-1)\).
All the diffuse bias and gap conditions permit these
parameters. We compare static Gamma bounds with (3).

Set \(A=10\) and \(g_j=\log(G_j^C/(jA))\).
The unit Gamma shape gives the initial endpoint-normalized
annulus of upper and lower log radii bounded by
\[
       b_j\le g_j-\log(1+1/(jA)),\qquad
       a_j\le\log(1+1/A)-g_j .
\]
In this display the right sides may be used instead of the
two-root bounds. Their sum for the diameter is less
than \(1/10\). Directly from the factorial endpoint product,
\[
      g_j=\frac1j\sum_{k=1}^j\log(1+k/(jA)),\qquad
      \frac{j+1}{2jA}-\frac1{2A^2}
                        \le g_j\le\frac{j+1}{2jA}.
\]
In particular \(U=g_N-g_n-6/(2nN)<0\), and
\[
       107/9100\ \le |U|\ \le 99/4550 .
\]
For the static Gamma interval the magnitude of its negative
endpoint, using the displayed upper bounds for \(a_N,b_n\),
is at least \(1503/18200\). It therefore contains \(U\).
Those same Gamma bounds give an upper bound
\(4619/47320\) for their sum. At \(W=10\) the contraction
factor of each order is below \(8/77\): use
\(c_j(W)\le e^{-W/4}\cosh(D(0)/2)\),
\(\cosh(1/2)\le8/7\), and \(e^{5/2}>11\).
Consequently the magnitude of the negative endpoint of
the contracted interval (3) is at most
\[
                    \frac{8}{77}\frac{4619}{47320}
                                      <\frac{107}{9100}.
\]
Thus the reference is radially separated there even though
the static annulus test does not certify separation. The
strict interval gap also excludes all sufficiently diffuse
ordinary tails in a neighborhood of these parameters.
The comparison is a radial exclusion and does not require
this parameter point to be an angular candidate.

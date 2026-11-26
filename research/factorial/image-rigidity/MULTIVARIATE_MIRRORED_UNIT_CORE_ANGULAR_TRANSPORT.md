# Angular transport near an exact mirrored unit-gap core

Take a distinguished shared mirrored pair
\((\ell,\ell+1),(\ell+1,\ell)\) with integer \(\ell\ge1\), and
any number of additional ordinary shared mirrored pairs
\((\ell_i+d_i,\ell_i),(\ell_i,\ell_i+d_i)\), where
\(\ell_i,d_i\ge1\). Let
\[
                         V_{\rm tail}=2\sum_i d_i^2/\ell_i .
\]
Let \(J\ge3\) be any integer at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If
\[
                V_{\rm tail}<\frac{1}{15625\pi^2 J^8\ell^3},
 \tag{1}
\]
then the profile has disjoint consecutive root arguments in
every degree and satisfies SFC(2). The pairs at terminal
\(N<J\) strictly interlace in symmetric angular order. One can
also take \(J\) beyond the
[eventual mirrored interlacing threshold](MULTIVARIATE_UNIFORM_EVENTUAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
to conclude strict interlacing in every degree. The bound in
(1) is polynomial in the threshold. It imposes smallness only
on the added rows.

The proof quantifies the exact consecutive Gegenbauer gap of the
core and uses angular transport under convolution to carry that
gap to the full profile.

## Angular transport under circular coefficient convolution

Fix a degree \(m\) and a monic reciprocal real polynomial \(Q\)
with simple unit-circle roots and nonzero coefficients. Define
\(T\) to be its binomial-basis convolution multiplier. Thus if
\(Q=\sum_k\binom mk q_k z^k\), then \(T\) multiplies the ordinary
coefficient of \(z^k\) by \(q_k\). We have \(q_m=q_0=1\) and
\(q_{m-k}=q_k\). Convolution preserves roots in the closed disk
and preserves roots on the circle. The disk statement is the
Grace–Szegő theorem with the second factor having circular roots.
The circle statement follows also by reversing a self-inversive
first factor. [Circular composition theorem](https://www.numdam.org/item/10.1016/j.crma.2014.10.017.pdf)

At each simple-circle input \(P\), the output is also simple.
Indeed, the set of simple-circle polynomials is the relative
interior of circle-rooted polynomials in the monic self-inversive
space with fixed endpoint phase. A simple circle root stays
on the circle under small perturbations in that space, by
uniqueness and reciprocal conjugation. A repeated root
\(\zeta\) can be split into \(\zeta r,\zeta/r\) for real
\(r\) near one, with the endpoint phase fixed, and so belongs
to the relative boundary. The invertible multiplier \(T\)
preserves this space and maps the circle-rooted set into
itself, hence maps its relative interior into itself.

Fix all but one input root \(a\). Write \(P=(z-a)R\),
\(A=T(zR)\), and \(B=T(R)\). With reversal taken in degree \(m\),
\(A^*=\xi B\) for a constant \(|\xi|=1\), by the
self-inversive identity for \(R\) and symmetry of the weights.
The output equation is \(A-aB=0\).
All roots of \(A\) lie in the closed disk by the convolution
theorem applied to \(zR\). Cancel its boundary factors in
\(A/A^*\) to obtain a finite Blaschke product times a unit
constant. At least one factor remains because \(A(0)=0\)
and \(A\) is monic. The boundary argument of this product
has strictly positive derivative, the sum of the positive
Poisson kernels of its interior roots. Thus each moving
output argument has nonnegative derivative in the argument
of \(a\). The canceled output roots are fixed.

It follows that the Jacobian of local output arguments in the
input arguments has nonnegative entries. Its row sums are
one, by rotation covariance. Its column sums are also one,
since \(T\) preserves both endpoint coefficients and hence
the product of the roots. In particular the largest absolute
instantaneous angular speed of any output is bounded by the
largest absolute speed of the inputs. More generally, suppose
two simple inputs have cyclically paired arguments with increasing
lifts \(\alpha_1<\cdots<\alpha_m<\alpha_1+2\pi\) and similarly
ordered lifts \(\beta_i\). Interpolate the lifts linearly. This
stays in the simple-circle chamber, including the circular gap at
the cut, because each gap is a convex combination of two positive
gaps. The output roots admit continuous lifts along the path.
Integrating the speed bound pairs the two output sets within the
input bound \(\max_i|\beta_i-\alpha_i|\). This comparison is
independent of where a lifted root crosses a chosen cut.

Suppose the input has circular roots
\(-e^{i\eta_i}\) with \(|\eta_i|\le\delta<\pi\).
Starting from the neutral factor \((1+z)^m\), move its roots
along \(-e^{it\eta_i}\), \(0\le t\le1\). For a simple final
input the roots are simple for \(t>0\), since the lifts of
the arguments do not differ by multiples of \(2\pi\).
The convolution output at \(t=0\) equals the simple polynomial
\(Q\). Integration of the speed bound along the path pairs
the final output roots with the roots of \(Q\), each with
angular displacement at most \(\delta\).

## A polynomial gap for the core alone

At degree \(m\), let \(\lambda_m=1+m\ell\). The exact core
polynomial has angles \(2\arccos x\), where \(x\) ranges over
the zeros of \(C_m^{\lambda_m}\). The Jacobi matrix of these
zeros has zero diagonal and positive off-diagonals
\[
 c_j(v)=\tfrac12\sqrt{1-u_j(v)},\qquad
 u_j(v)=\frac{v(v-1)}{(j+v)(j+v-1)} .
\]
For \(N=n+1\) and \(n\ge2\), the known consecutive comparison
uses the leading \(n\)-block and trailing \(n\)-block of
the upper Jacobi matrix. The leading-block comparison ratios
\(c_j(\lambda_n)/c_j(\lambda_N)\) decrease with \(j\) and
are at least one. The trailing-block ratios
\(c_{j+1}(\lambda_N)/c_j(\lambda_n)\) satisfy the same
properties. The ratio comparison uses diagonal congruence
of the two bipartite matrices and therefore remains valid
after factoring the minimum ratio out as a common scale.
[Consecutive Gegenbauer Jacobi comparison](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#consecutive-interlacing-for-a-single-mirrored-unit-gap)

Let \(\mu_i>0\) be a lower eigenvalue in decreasing order, and
write \(\nu_i,\nu_{i+1}\) for its corresponding upper neighbors.
If \(\rho_-=c_{n-1}(\lambda_n)/c_{n-1}(\lambda_N)\) and
\(\rho_+=c_n(\lambda_N)/c_{n-1}(\lambda_n)\), principal-block
interlacing therefore gives
\[
   \mu_i-\nu_{i+1}\ge(1-\rho_-^{-1})\mu_i,\qquad
   \nu_i-\mu_i\ge(\rho_+-1)\mu_i .
 \tag{2}
\]
For \(j=n-1\) and \(v\in[\lambda_n,\lambda_N]\),
\(v\le2n\ell\), \(j+v\le(5/2)n\ell\), and \(u_j(v)\ge4/25\).
Thus
\[
 u_j'(v)
  =j u_j(v)\left(\frac{1}{v(j+v)}
                   +\frac{1}{(v-1)(j+v-1)}\right)
  \ge\frac{2}{125n\ell^2}.
\]
Integrating over an interval of length \(\ell\) gives
\(1-\rho_-^{-1}\ge1/(125n\ell)\), using
\(1-\sqrt{1-s}\ge s/2\) for \(0\le s\le1\).
For the trailing comparison, write
\(u_j(v)=[v/(j+v)][(v-1)/(j+v-1)]\).
The first factor decreases by
\[
 \frac{\ell+1}{(\lambda_n+n-1)(\lambda_N+n)}
       \ \ge\ \frac{1}{6n^2\ell}
\]
between \((j,v)=(n-1,\lambda_n)\) and \((n,\lambda_N)\).
The second factor also decreases, and its lower-degree value
is at least \(1/2\). Thus
\(u_{n-1}(\lambda_n)-u_n(\lambda_N)\ge1/(12n^2\ell)\).
Since \(c_j(v)\le1/2\), the difference-of-squares formula
gives \(\rho_+-1\ge1/(24n^2\ell)\).

The smallest positive eigenvalue at degree \(n\) is at least
\(1/(2n^{3/2}\sqrt\ell)\). In fact the odd-even block has
a square lower-bidiagonal initial block with diagonal
\(c_1,c_3,\ldots\) and subdiagonal \(c_2,c_4,\ldots\).
Its inverse entries have modulus at most \(1/c_1\) because
the \(c_j\) increase with \(j\). The inverse has operator norm
at most \(n/c_1\), and inserting any extra row cannot reduce
the smallest singular value. Finally
\(c_1(\lambda_n)=1/\sqrt{2(\lambda_n+1)}
\ge1/(2\sqrt{n\ell})\).
Combining with (2) gives a real cross-gap at least
\(1/(250n^{7/2}\ell^{3/2})\). By reflection the same estimate
applies to negative eigenvalues. Gaps to a structural zero
are bounded by the even-degree smallest positive eigenvalue,
so the weaker bound with \(N\) in place of \(n\) holds
for these gaps too. The angle map has derivative of magnitude
at least two. Angles in opposite half-circles are at least
as far apart as their same-side reflections, so the core's
circular cross-gap is at least
\[
                             \frac{1}{125N^{7/2}\ell^{3/2}} .    \tag{3}
\]

## Transfer to the ordinary tail

Normalize the tail polynomial to be monic. Its roots lie on
the unit circle by the
[mirrored factor theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order),
and are simple unless the tail is empty.
If they are written \(-e^{i\eta_i}\) with principal
\(\eta_i\in[-\pi,\pi]\), then
\[
 \sum_i(1-\cos\eta_i)=m(1-a_{m,1}),\qquad
 -\log a_{m,1}\le\sum_i\frac{d_i^2}{\ell_i}=V_{\rm tail}/2 .
\]
Here \(m a_{m,1}\) is the monic coefficient at degree \(m-1\).
The last inequality follows by summing \(\log(1+x)\le x\)
in the exact row factor
\(\prod_{h=1}^{d_i}(m\ell_i+h)/(m(\ell_i+d_i)-d_i+h)\).
Using \(1-\cos\eta_i\ge2\eta_i^2/\pi^2\) bounds every
\(|\eta_i|\) by \((\pi/2)\sqrt{mV_{\rm tail}}\).

The full polynomial is the circular coefficient convolution
of the exact core and tail factors. For \(N<J\), angular
transport moves its two root sets from their core
references by a total at most \(\pi\sqrt{J V_{\rm tail}}\).
By (1) this is smaller than (3) throughout the interpolation
path in both degrees. Thus the core's strict cross-degree
interlacing holds for the full pair in these degrees.
The linear-to-quadratic pair also interlaces, by parity and
simplicity. The eventual theorem covers \(N\ge J\).
All factorial polynomials have rational coefficients, so the
result gives SFC(2) over every characteristic-zero field.

## Asymmetric ordinary tails

For the same distinguished unit-gap core, an ordinary tail may
also contain unmatched orientations. Put
\[
                  S_{\rm tail}=\sum_i
                      \frac{(a_i-b_i)^2}{\max(a_i,b_i)}
\]
over the unequal rows of the tail. If
\[
 S_{\rm tail}<
     \min\left\{\frac{10^{-16}}{4J^2},
                 \frac{1}{36\,000\,000\,J^7\ell^3}\right\},
 \tag{4}
\]
the exact full profile has disjoint consecutive root arguments
and satisfies SFC(2). This conclusion concerns angles even
when the tail and the full polynomial have roots off the circle.
The bound itself forces all unequal tail rows to be shared.

First obtain a small complex disk for the endpoint-normalized
tail roots. For \(2\le m<J\), apply the matched Gauss root-disk
construction at lower order \(m\) and upper order \(M=m+1\).
Its total reciprocal parameter \(X\) satisfies
\(X\le S_{\rm tail}/[m(1-S_{\rm tail})]\le2S_{\rm tail}/m\).
The first bound in (4) therefore implies
\[
             X\le10^{-16}M^{-3},\qquad M X\le3S_{\rm tail}.
\]
Here \(M\le J\) and \(M\le2m\) in the first inequality, and
\(M/m\le3/2\) in the second.
The root-disk proof gives the endpoint-normalized roots as
\(a=-e^{-\sqrt X\,z}\), with \(|z|\le3\sqrt M\), by the
radius-\(1/(80\sqrt M)\) disks around scaled Hermite roots.
Indeed the Hermite root bound is \(2\sqrt M\), and the
Gaussian variance in those disks is at most one. We use
only the complex disk assertion, which holds at
\(\sqrt X\le10^{-8}M^{-3/2}\).
[Matched Gauss root disks](MULTIVARIATE_MATCHED_GAUSS_UNIFORM_REMAINDER.md#a-uniform-complex-disk-for-the-matched-deformation)
Thus every root admits a log coordinate
\(a=-e^\omega\) with \(|\omega|\le\delta=6\sqrt{S_{\rm tail}}\).
The case of no unequal tail rows is the exact neutral factor.

Let \(Q\) now be the exact core polynomial, with binomial
coefficients \(q_k\), and \(\widetilde Q\) its symmetric
multi-affine polarization. At a root \(z\) of the full
convolution, one has
\[
                    \widetilde Q(-z/a_1,\ldots,-z/a_m)=0.
\]
This follows by writing the coefficients of the tail polynomial
in terms of the elementary symmetric polynomials of \(a_i^{-1}\).
The endpoint normalization of the tail only rescales the
entire argument by a positive number upon convolution.
The polarized arguments are \(z e^{-\omega_i}\), all in the
closed disk of center \(z\) and radius \(|z|(e^\delta-1)\).
By the Grace–Walsh–Szegő coincidence theorem the disk
contains a root \(\zeta\) of \(Q\).
[Coincidence theorem, Theorem 11](https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf)

The resulting circular angular distance from \(z\) to \(\zeta\)
is at most \(\arcsin(e^\delta-1)\le4\delta\), since (4)
gives \(e^\delta-1\le2\delta<1/2\).
For a pair with \(N<J\) the sum of these two angle errors
is at most \(48\sqrt{S_{\rm tail}}\). The second bound in
(4) makes this less than the exact core gap (3).
Thus there is no common root argument in these degrees.
The linear-to-quadratic case follows from positivity of the
quadratic moment for every real argument, and the absolute
eventual angular cutoff covers the remaining pairs.

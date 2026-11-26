# Eventual interlacing throughout subcritical mirrored curvature

Let the distinct fully shared ordinary exponent columns be
mirrored, omitting equal rows. Write \(W=\sum_i
d_i^2/(\ell_i+d_i/2)\) over both orientations, and denote the
order-\(j\) endpoint-normalized log-factorial symbol per degree
by \(L_j\). For every fixed \(0<W_+<4\) there is a terminal
threshold, uniform in the number of rows and their integer
exponents, with the following property. When \(W\le W_+\),
for every sufficiently large \(N=n+1\), the actual order-\(n\)
and order-\(N\) factorial polynomials strictly interlace in
root argument, with the circle cut at the positive real ray.
So do the two polynomials with the upper symbol frozen on
both grids:
\[
        \sum_{k=0}^m\binom mk e^{m L_N(k/m)}z^k,\qquad m=n,N.
\]
All four polynomials have simple roots on the unit circle in
endpoint normalization. This extends the
[uniform weak common-symbol theorem](MULTIVARIATE_UNIFORM_WEAK_COMMON_SYMBOL_INTERLACING.md)
to the full open subcritical range at large orders. The
threshold is unevaluated.

## A uniform angular gap in a compact subcritical band

Write \(G=L_\infty\) for the Stirling chord symbol and
\(\rho(q)=q(1-q)G''(q)\).
Each summand \(d_i^2q(1-q)/(a_i+(b_i-a_i)q)\) is strictly
concave on \([0,1]\): its second derivative is
\(-2d_i^2 a_i b_i/(a_i+(b_i-a_i)q)^3\).
The sum is symmetric under \(q\mapsto1-q\).
Thus \(\max\rho=\rho(1/2)=W/4\).
For \(0<W_-\le W\le W_+<4\), use the
[fixed-band angular comparison](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md)
with the subcritical margin \(1-W_+/4\).

For the mirrored profiles the chord-centered field \(w\) is
purely imaginary at every root. Indeed the endpoint center
and signed degree shift vanish, and the
[mirrored circle-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
puts the roots on the circle. Angles may be reflected to
\([0,\pi]\), since the circular distance between opposite
sides is at least the difference of their unsigned angular
coordinates. In this setting the fixed-band comparison gives
the following quantitative strengthening. There are a constant
\(g>0\) and a terminal threshold, uniform on this compact band,
such that every pair of adjacent-order actual root arguments
has circular distance greater than
\[
                              g n^{-5/3}.                    \tag{1}
\]

For the proof, suppose two physical arguments differ by at
most \(g n^{-5/3}\). For sufficiently large \(n\) their imaginary
fields lie in a common chart of the fixed-band cover.
The subcritical upper boundary is one-saddle and has a
zero-free neighborhood. The regular interior and the lower
boundary use the common phase difference \(D\) and amplitude
ratio \(A\), with
\[
       e^{jD(w_j)}=A(w_j)(1+O(1/j)),\qquad j=n,N .
\]
Dividing the two equations over the joining segment bounds
\(e^{D(w_n)}-1\) by \(O(n^{-1}+g n^{-2/3})\).
This contradicts the strictly positive regular chart margin.
The lower boundary uses its exact reflection factor in \(D\),
as in the fixed-band proof.

At a subcritical fold, the two Airy coordinates
\(z_j=j^{2/3}\lambda(w_j)\) are each within \(O(n^{-1/3})\)
of an Airy zero. The proposed physical field displacement
changes the scaled coordinates by at most
\[
 N^{2/3}\|\lambda'\|\,g n^{-5/3}
       +(N^{2/3}-n^{2/3})\|\lambda\|.
\]
The first term is \(O(g/n)\), and the second is the term
controlled by the fixed chart shrinkage in the
[fixed-band Airy label comparison](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#the-subcritical-fold-and-complex-order).
Its strict Airy spacing bound still forces the same label.
The same analytic order branch then has imaginary derivative
of fixed sign and magnitude at least a positive multiple
of \(|a_k|n^{-5/3}\). The exponentially small integer-order
errors cannot remove that displacement. Taking \(g\) smaller
than the resulting uniform lower bound proves (1).

For the deformation below, the same charts can use real
smaller exponents \(\ell_i\ge1\) with integer gaps \(d_i\).
Here are the needed extensions of the fixed-band argument.
The exact Rice residues are at \(q=k/j\) and recover the
Gamma coefficients for real smaller exponents, since the
rational binomial kernel and the integer order \(j\) are
unchanged. The logarithmic boundary heights of the stationary
map stay integer multiples of \(\pi\) because the gaps
stay integer. The compactness estimates still apply: at any
positive cutoff on the row relative gaps, the numbers of
retained rows and the integer gaps are bounded, and the
retained smaller exponents now range over bounded real
intervals. Subsequence compactness and interval meshes replace
the discrete choices of these exponents. The small relative
gap estimates are unchanged.

For clarity, the mirrored subcritical branch geometry itself
also has a direct real-exponent check. With the sum counting
each retained pair once and \(a_i=\ell_i/d_i+1/2\), the centered
map on \(q=1/2+iv\) has imaginary part
\[
       Y(v)=2\sum_i d_i\arctan(v/a_i)-2\arctan(2v)+\tau v.
\]
Here the diffuse curvature \(\tau\) is zero for an actual
finite row list and nonnegative in a profile limit.
The quantity
\[
       (v^2+\tfrac14)Y'(v)
            =2\sum_i d_i a_i\frac{v^2+1/4}{v^2+a_i^2}
                     +\tau(v^2+1/4)-1
\]
is strictly increasing in \(v^2\), starts at \(W/4-1<0\)
and has positive limit or tends to infinity. Thus the central vertical
branch has a simple fold in the open target strip.
The real-gap inversions give the proper two-sheet count as in
the fixed-band proof. No other interior fold can occur: the
component of the strip preimage containing the simple fold
already has degree two, and every such component is simply
connected by the maximum principle for the imaginary part
of the stationary map. Below the fold on the imaginary
axis the two saddles are reflected by \(q\mapsto1-\bar q\).
The exact Rice contour and its reflected integrand give
crossing weights of equal magnitude, and at least one is
nonzero by the same bounded-gap winding test. In the purely
diffuse limit the negative-real-ray winding test replaces
the bounded gap. Thus the active equal-height arc on the
mirrored field line has both terms. There is no exceptional
one-row class for a mirrored profile. These constructions
and the local Gamma bounds persist on the compact real
parameter meshes, so the regular and Airy margins used
in (1) are uniform there as well.

## Uniform freezing of the lower grid

The [full-mesh transfer theorem](MULTIVARIATE_MATCHED_LAPLACE_MESH_SHARPENING.md)
and the [symmetric compensation estimate](MULTIVARIATE_UNIFORM_WEAK_MIRRORED_LAPLACE_TRANSFER.md)
give, uniformly for \(0<W\le W_+\), angular motion at most
\[
                          C(W_+)\frac{1+\log n}{n^2}          \tag{2}
\]
from the actual lower polynomial to the upper-symbol
polynomial on the lower grid, once \(n\) is sufficiently large.
The paired roots remain simple and on the unit circle.
Here the smallness condition for the symmetric compensation is
obtained from large order, rather than from weak total curvature.

To check this uniformity, let \(L=\min_i\ell_i/d_i\).
The mirrored occurrence attaining \(L\) gives
\(L\ge L_{\min}=2/W_+-1/2>0\).
If \(L\le1600\), take \(n\ge4/L_{\min}\) in the full-mesh
estimate. Its constants are bounded uniformly on the compact
interval of \(L\). The rate gap satisfies
\[
 a_i-K_*\ge (\ell_i/d_i)/16
             \ge a_i/[16(1+1/(2L_{\min}))].
\]
In its jump bound the sum of inverse rate gaps and inverse
squared rate gaps is therefore bounded by a constant times
\(W\), since \(\sum_i1/a_i\le W\) over both orientations.
The weighted jump parameter is \(O(W/n)\), and the exponential
factor remains uniformly bounded as \(n\) increases.
The Rouché radius is \(O(W/n)\) in the scaled log coordinate,
giving \(O(W/n^2)\) angular motion in this regime.

If \(L\ge1600\), keep the improved Gegenbauer mesh and all
the constants of the symmetric compensation estimate.
Its second weighted jump moment \(T\) satisfies
\(K_1^2T\le648 W/n\) and \(b^2T\le2W(1+\log n)^2/n\).
Both tend to zero uniformly for \(W\le W_+\).
For all sufficiently large \(n\), the same numerical bounds
on the exponential and quadratic factors give
\[
           \rho_{\rm disk}/a_0
                  \le394752 W(1+\log n)/n<1,\qquad
           \rho_{\rm disk}/n
                  \le12336W(1+\log n)/(n^2\sqrt L).
\]
This proves (2) in the other regime. These estimates allow
real positive smaller exponents and gaps that remain integral.

For \(W\) bounded below, increase the terminal threshold so
(2) is less than half the right side of (1), using
\((1+\log n)n^{-1/3}\to0\).
The frozen pair then has disjoint arguments and the same
ordering across degrees as the actual pair.

## Setting the interlacing order

Choose a positive \(w_0\) below both \(W_+/2\) and the
absolute weak cap from the weak common-symbol theorem.
For \(W\le w_0\) that theorem fixes the strict interlacing
order. At large orders its quantitative comparison also
gives that order for the actual pair.
For \(w_0\le W\le W_+\), continuously increase all smaller
exponents by a common factor, keeping the gaps and mirror
pairing fixed, until the curvature reaches \(w_0\).
Along this path \(W\) decreases and remains in
\([w_0,W_+]\). Every row contribution to \(\rho\) decreases
pointwise as well. Thus the same uniform subcritical chart
and jump comparisons apply throughout.
Both pairs stay simple, have no common argument across
degrees, and retain the positive real cut since the
coefficients are positive. Their interlacing order at the
weak endpoint therefore holds at the original profile.

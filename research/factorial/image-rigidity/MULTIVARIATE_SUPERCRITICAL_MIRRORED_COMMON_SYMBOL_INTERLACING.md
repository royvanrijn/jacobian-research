# Eventual interlacing throughout supercritical mirrored curvature

Let the two distinct ordinary exponent columns be fully shared
and mirrored, omitting equal rows. Write
\[
                  W=\sum_i\frac{d_i^2}{\ell_i+d_i/2}
\]
over both orientations, with smaller exponents \(\ell_i\) and
gaps \(d_i\). Let \(L_j(s)\) be the endpoint-normalized
log-factorial symbol per degree at order \(j\).

**Theorem.** For each fixed compact band \(4<W_-\le W_+<\infty\),
there is a uniform terminal threshold for all profiles with
\(W\in[W_-,W_+]\). At every larger terminal \(N=n+1\) the actual
factorial polynomials at orders \(n,N\) strictly interlace in
argument, with the circle cut at the positive real ray.
The same is true of the arrays
\[
       \sum_{k=0}^n\binom nk
           e^{n[(1-\lambda)L_n(k/n)+\lambda L_N(k/n)]}z^k,
       \qquad
       \sum_{k=0}^N\binom Nk e^{N L_N(k/N)}z^k,
                         \qquad 0\le\lambda\le1 .
\]
Each polynomial has simple roots on the unit circle. Thus the
assertion includes strict interlacing at the common upper symbol
\(\lambda=1\). The threshold is unevaluated.

Away from the positive ray, the proof uses a regular saddle gap
and the nonlinear same-grid transfer. At the positive ray, the
two saddle contributions are conjugate for a mirrored profile.
Taking a single half-polynomial makes that reflection exact,
including the phase error from replacing the lower symbol.

## The mirrored contour hypotheses

Use the notation \(G=L_\infty\) for the Stirling chord and
\(\rho(q)=q(1-q)G''(q)\) as in the
[fixed-band angular comparison](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md).
Each row contribution to \(\rho\) is strictly concave and the
sum is symmetric, so \(\max\rho=W/4\). At supercritical
curvature the real entropy phase
\[
                   G(q)-q\log q-(1-q)\log(1-q)
\]
has its two equal global maxima strictly on either side of
\(1/2\). Write them as \(q_-<1/2<q_+=1-q_-\).

The regular and positive-boundary saddle estimates used below
also allow real smaller exponents \(\ell_i\ge1\) with integer
gaps. The
[Rice coefficient residues](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases)
still give the exact Gamma coefficients. The analytic row
estimates and Stirling bounds are unchanged. The integer-gap
logarithmic heights on the real cuts are unchanged as well.
For the compact profile approximation, bounded relative gap
cutoffs leave boundedly many retained rows and bounded integer
gaps. The retained real smaller exponents can be meshed over
bounded intervals. The small-relative-gap approximation and
finite strict contour tests then apply as in the
[regular cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md).
Coincident limiting real poles have the same sign of their
logarithmic coefficients, so they do not permit a new bounded
fiber at the cuts.

In this mirrored case the supercritical inverse-sheet divider
also follows directly. For \(q=1/2+iv\), \(v>0\), the
stationary height is
\[
             Y(v)=2\sum_i d_i\arctan(v/a_i)
                         -2\arctan(2v)+\tau v,
          \qquad a_i=\ell_i/d_i+\tfrac12 ,
\]
counting each retained pair once, with nonnegative diffuse
curvature \(\tau\) in a limit and \(\tau=0\) for a finite
profile. Then
\[
    (v^2+\tfrac14)Y'(v)
       =2\sum_i d_i a_i\frac{v^2+1/4}{v^2+a_i^2}
                              +\tau(v^2+1/4)-1
\]
is increasing and begins at \(W/4-1>0\).
Thus the central vertical line maps strictly above the target
strip \(-\pi<\Im w<0\). The two real-gap inverses belong to
opposite halves of the upper saddle plane, so the proper
two-sheet map has one unramified sheet in each half. The
real-gap count and the strict contour tests give the same
two-saddle regular witnesses in the real-smaller-exponent
case. In the interior, reflection in \(q=1/2\) makes the
crossing weights equal in magnitude on the imaginary field
line. At least one is nonzero by the winding test, or its
purely diffuse negative-real-ray version. There is no
one-row exception for a mirrored profile.

## Away from the positive axis

The actual mirrored polynomials have simple unit-circle roots
by the [Gegenbauer factorization](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md).
Their chord-centered fields \(w\) are purely imaginary.
On a fixed closed angular interval away from the positive
ray there is a uniform actual adjacent angular gap \(g/n\),
with \(g>0\) depending on the band and the interval.
For near-coincident roots the fields enter a common regular
chart of the fixed-band cover, including the negative-ray
boundary. The actual root equations there have the form
\[
                 e^{jD(w_j)}=A(w_j)(1+O(1/j)),\qquad j=n,N,
\]
with the same analytic functions and a positive lower bound
on \(|e^{D(w_n)}-1|\). If the fields differ by at most \(g/n\),
the ratio of the root equations would give
\(e^{D(w_n)}=1+O(g+1/n)\).
Choosing \(g\) sufficiently small contradicts that margin.
The comparison of circular angles reduces to the same
upper half by reflection.

For the change on the lower grid the
[uniform freezing argument](MULTIVARIATE_SUBCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md#uniform-freezing-of-the-lower-grid)
applies with any finite upper bound on \(W\).
Its needed rate parameter \(L=\min_i\ell_i/d_i\) has
a positive uniform lower bound even when the band exceeds
four. In fact \(L\ge\min(1,1/W_+)\): if \(L<1\), the pair
attaining \(L\) has \(d_i\ge1/L\), hence
\(W\ge2d_i/(L+1/2)\ge4/(3L)\).
The bounded-\(L\) and large-\(L\) estimates of that proof
therefore give the same uniform angular movement
\(O((1+\log n)/n^2)\) and preserve simple circle roots.
They also apply to every \(\lambda\in[0,1]\), since the
positive increment measure is just multiplied by \(\lambda\).
This movement is smaller than the actual gap away from the
positive ray.

## A single analytic half at the positive ray

For the order-\(m\) actual array write
\(c_{m,k}=\binom mk e^{mL_m(k/m)}\).
Use the exact left half
\[
 \mathcal H_m(\theta)=e^{-im\theta/2}
            \left(\sum_{k<m/2}c_{m,k}e^{ik\theta}
                      +\tfrac12{\bf1}_{2\mid m}
                              c_{m,m/2}e^{im\theta/2}\right).
\]
Define \(\mathcal H_n^\lambda\) in the same way for the
interpolated lower array. For real \(\theta\) the full
centered polynomial is \(2\Re\mathcal H_m(\theta)\), or
\(2\Re\mathcal H_n^\lambda(\theta)\), by exact reflection.

The [positive-boundary lattice expansion](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md#upper-boundary)
on a small complex \(\theta\)-disk applies to the left half
alone. Its retained interval is around \(q_-\), strictly
inside \((0,1/2)\). All other terms of this half have a
strict entropy loss at \(\theta=0\), including the optional
half of the center term. The loss persists on a smaller
complex disk, so the off-interval sum is exponentially small
relative to the retained packet. The same cotangent panel
deformation around the left interval gives the one-saddle
expansion
\[
        \mathcal H_m(\theta)
              =A(\theta)e^{m\Phi(\theta)}(1+r_m(\theta)),
                    \qquad \|r_m\|\le C/m,                 \tag{1}
\]
where \(A\) is nonzero, \(A(0)>0\), and \(\Phi(0)\) is real.
The normalized phase here is
\[
  \Phi(\theta)=
        G(q_-(\theta))-q_-(\theta)\log q_-(\theta)
        -(1-q_-(\theta))\log(1-q_-(\theta))
        +i\theta(q_-(\theta)-1/2),
\]
with \(q_-(\theta)\) the continued saddle.
All constants are uniform on the compact band. A positive
common normalization, if used in the lattice expansion, has
no effect on the phase.

The same leading \(A,\Phi\) in (1) work for
\(\mathcal H_n^\lambda\). On the saddle disk and its contour
neighborhood, write the analytic Gamma expansion as
\[
 L_j(q)=G(q)+A_0(q)/j+A_1(q)/j^2+O(j^{-3}),
\]
using the Stirling terms through inverse order and subtracting
the endpoint chords. The bounds are uniform there because the
row arguments divided by \(j\) have positive real parts, and
the reciprocal powers sum to a bound from \(W_+\).
The additional coefficient multiplier for the interpolated
array is therefore
\[
          \exp\big(n\lambda(L_N-L_n)(q)\big)
                  =1-\lambda A_0(q)/N+O(n^{-2}).
\]
It is \(1+O(n^{-1})\) on the retained lattice contour.
Applying the same local Gaussian expansion to this analytic
multiplier changes only the \(O(n^{-1})\) relative remainder
in (1). The strict entropy loss outside the left interval
also persists for the interpolated coefficients.

Put \(d(\theta)=-\Im\Phi(\theta)\) on real \(\theta\).
Then \(d(0)=0\) and
\(d'(0)=1/2-q_->0\), because the phase is stationary in \(q\).
This derivative stays uniformly positive on a smaller
positive interval. Define \(f_m=-\arg\mathcal H_m\), and
\(f_n^\lambda=-\arg\mathcal H_n^\lambda\), with each branch
zero at zero. Equation (1) and Cauchy's estimate for
\(\log(1+r_m)\) give
\[
 \begin{split}
      f_m(\theta)&=m d(\theta)-\Im\log A(\theta)
                                      +O(\theta/m),\\
      f_n^\lambda(\theta)&=n d(\theta)-\Im\log A(\theta)
                                      +O(\theta/n).
 \end{split}
\]
Indeed the remainder logarithms are real at zero, so the
imaginary errors are proportional to \(\theta\).
The derivatives of the phases are positive for all large
orders. After another fixed decrease of the interval,
\[
            0<f_N(\theta)-f_n^\lambda(\theta)<\pi
                            \qquad (\theta>0).
\]
A circle root has half phase in \(\pi/2+\pi\mathbb Z\).
Consequently there can be no common argument near the
positive ray, including at all interpolation parameters.
At the ray itself the coefficients are positive. Together
with the comparison away from the ray, this proves disjoint
arguments along the whole interpolation on the lower grid.

## Orienting the interlacing

For the actual ordinary pair, add 27 mirrored unit-gap
pairs, increasing their curvature continuously from zero by
decreasing their common smaller exponent from infinity to
one. At the infinite endpoint their normalized multipliers
tend to one. On this path \(W\) stays in a fixed compact
supercritical band, since the added terms increase \(\rho\)
pointwise and contribute at most 36 to \(W\).
The same real-exponent chart estimates apply along the path.
The actual roots remain simple on the circle with disjoint
arguments across degrees once the uniform threshold for that
enlarged band is met.

At the other endpoint the ordinary absolute curvature
\(S\) is at least 27. The mirror case of the
[absolute curvature 27 half-phase proof](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#paired-radius-and-phase-comparison)
identifies the strict interlacing order there. To make the two
halves reflected copies, transfer half of the even middle
coefficient from the first half to the second. Its coefficient
bound on the disk of radius \(6/5\) is
\((12e\,e^{-27/4}/5)^{m/2}\), by the same entropy bound.
The half-log norm \(3(D+1)\), with \(D=N e^{-S/2}\),
shows this correction is exponentially small relative to
either half. Cauchy's estimate on the smaller disk preserves
the strict derivative margins. With this symmetric half
denoted by \(H_m\), the functions on the circle
\[
     p_m(\theta)=m\theta-2\Im\log H_m(e^{i\theta})
\]
in that proof are strictly increasing from zero to \(m\pi\)
on \([0,\pi]\). The half-log derivative bounds give
\(|p_m'-m|<m\) and
\(|p_N-p_n-\theta|<\theta\) throughout the open interval.
The roots satisfy \(p_m(\theta)\in(2\mathbb Z+1)\pi\).
At each lower root the upper phase is strictly between
that level and the next one. Monotonicity and reflection
therefore give strict interlacing, with the structural root
at the negative ray in the odd degree.

The absence of common arguments along the added-row path
preserves that order at the original profile. The earlier
lower-grid interpolation then preserves it from the actual
pair to every \(\lambda\), in particular to the common
upper symbol.

# All-degree angular separation from a large mirrored pair with arbitrary other rows

Suppose an ordinary two-monomial profile contains
\((\ell,\ell+d),(\ell+d,\ell)\) with integers \(\ell\ge1\) and
\(d\ge7\ell\). All other ordinary rows are arbitrary, including
unmatched orientations and private rows. Then consecutive factorial
moment polynomials have disjoint root arguments in every degree and
satisfy SFC(2) over every characteristic-zero field.

The same angular conclusion holds if a common physical Gaussian
strength \(\tau\ge0\) multiplies the endpoint-normalized coefficient
at index \(k\) of order \(m\) by \(e^{-\tau k(m-k)/(2m)}\).
The proof controls the two oriented coefficient halves separately.
It needs no diffuse assumption, signed radial condition, or
individual-degree root checks.

## The two endpoint estimates

Let \(s_h=2d^2/(\ell+d)\), and let \(s_0\) be the absolute curvature
of the other rows plus \(\tau\). Put
\[
 S=s_h+s_0,\qquad q=\left(\frac{\ell}{\ell+d}\right)^d e^{-s_0/2},
 \qquad L_*=\frac{d^2}{\ell}+\frac{s_0}{2}.
\]
Also let \(V\) be the sum of \(d_i^2/\min(a_i,b_i)\) over shared
unequal rows, of \(\max(a_i,b_i)\) over private rows, and of \(\tau\).
Write \(W=1+V\). The bounds \(S\le V\le S(S+2)\) hold, as in the
[absolute curvature comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md),
and \(L_*\le V/2\).

In either orientation, let \(c_{m,k}\) be the endpoint-normalized
coefficient at index \(k\), including the binomial coefficient
and the Gaussian factor. For \(0<k\le m/2\), put \(w=k/m\).
The exact selected pair multiplier satisfies
\[
 \prod_{h=1}^{dk}\frac{m\ell+h}{m(\ell+d)-dk+h}
       \ \le\ \left(\frac{\ell+dw}{\ell+d}\right)^{dk}.
\]
The rest of the rows, including the Gaussian, supply
the upper bound \(e^{-s_0 k(1-w)/2}\) by the
[discrete curvature inequality](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate).
Using either \(\binom{m}{k}\le(em/k)^k\) or
\(\binom{m}{k}\le m^k/k!\), and \(\log(1+x)\le x\), gives
\[
 c_{m,k}^{1/k}\le \frac e w
       \left(\frac{\ell+dw}{\ell+d}\right)^d e^{-s_0(1-w)/2},
 \qquad
 c_{m,k}\le\frac{(mq)^k}{k!}e^{L_* k^2/m}.                    \tag{1}
\]
The first bound takes its maximum on a closed interval of
\(0<w\le1/2\) at one end. Its logarithmic derivative times \(w\)
is \(-1+ (d^2/\ell)w/(1+(d/\ell)w)+s_0w/2\), increasing in \(w\).
The value at \(w=1/2\) is at most \(2e(9/16)^7<1/10\).

## Short halves and large absolute curvature

Put \(N=n+1\) and assume \(n\ge2\). Use the oriented halves
\(F_m(z)=L_m(z)+z^m U_m(1/z)\), assigning the middle coefficient,
if any, to just one half. Each half starts at one.
For \(N\le30001\), the upper bound at \(w=1/m\) in (1) is at
most \(em[(1+7/m)/8]^7\). Its logarithmic derivative has the
sign of \(m+7-49\), so the maximum for \(2\le m\le30001\)
is at an endpoint. At \(m=2\) it is \(2e(9/16)^7<1/10\).
At \(m=30001\) use \((1+7/30001)^7<2\) and \(e<3\)
to obtain the same bound. Hence each nonconstant entry in
both halves is at most \(10^{-k}\).

If \(N>30001\) and \(S\ge2\log N+5\), the
[large-curvature geometric estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
with this \(S\) gives \(8^{-k}\) instead. In these two regions
let \(R_m=e^{1/(2m)}\). The geometric ratio on this disk is
less than \(1/7\), using \(e^{1/4}<4/3\) in the first and
\(e^{1/(2m)}<8/7\) in the second. Either half \(H\) satisfies
\[
          |H-1|<1/6,\qquad |zH'/H|<7/30<1/4 .
\]

In the remaining region \(N>30001,\ S<2\log N+5\), set \(D=Nq\).
As in the
[short endpoint estimate for a mirrored core](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#direct-bounds-in-the-short-endpoint-ranges),
\(L_*/n<1/32\), using \(d\le S\) and
\[
 \frac{L_*}{n}\le
 \frac{(2\log N+5)^2+(2\log N+5)/2}{N-1}.
\]
If \(D\le1/10\), (1) gives \(c_{m,k}\le(1/8)^k/k!\)
for \(k\le7\), and \(c_{m,k}\le10^{-k}\) for
\(8\le k\le m/2\) by the endpoint test at \(k=8\) and \(m/2\).
At \(R_m=e^{1/(2m)}<101/100\), the first seven coefficients
contribute at most \(u/(1-u)\) to each of \(|H-1|\) and
\(|zH'|\), where \(u=101/800\). This follows from the exponential
series for \(H-1\) and from \(u e^u\le u/(1-u)\) for the
derivative. The remaining geometric sums and their derivatives
at ratio \(101/1000\), starting at eight, are each below \(1/1000\).
Thus both norms are below \(3/20\), and \(|zH'/H|<1/4\).

In all three short regimes, \(R_m^m=\sqrt e>7/5\) gives
Rouché domination of the outer half at \(R_m\) and of the
inner half for the reciprocal. All roots satisfy
\(|\log|z||<1/(2m)\). At a root with argument
\(0<\theta\le\pi\), the root phase equation is
\[
 m\theta=(2j_m+1)\pi+
     \Im\{\log L_m(e^{u_m+i\theta})-\log U_m(e^{-u_m-i\theta})\},
 \qquad u_m=\log|z|.
\]
The logarithms are real along the positive radius. By the
derivative bounds, the phase correction has modulus less than
\(\theta/2\). At a common argument subtracting the two phase
equations would give a multiple of \(2\pi\) strictly between
zero and \(2\theta\le2\pi\). This rules out common arguments
in the three regions.

## Long halves and the change of order

Now \(D>1/10\) and \(S<2\log N+5\). The selected pair implies
\[
 Wq\le (s_h+s_0+1)^2
          \left(\frac{\ell}{\ell+d}\right)^d e^{-s_0/2}
      \le (2d+1)^2 8^{-d}\le225/8^7,
\]
since the first expression decreases in \(s_0\ge0\),
\(s_h\ge d\ge7\), and \((2d+1)^2 8^{-d}\) decreases for
\(d\ge7\). Since \(q\le8^{-7}\) and \(D>1/10\), \(n>200000\).
With \(K=\lceil32D\rceil+14\) this gives
\[
 D/n<1/720000,\qquad WD/n<1/9000,\qquad
 K\le182D,\qquad WK/n<182/9000,\qquad K+1<n/2 .                \tag{2}
\]
The first bound of (1) applied at \(k=K,m/2\) gives
\(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\).
At the first endpoint use
\((e/32)\exp(L_*K/n)<1/8\).

Consider a half with all its starting exponents positive.
Let \(B_m\) be its first coefficient, interpolated for real
\(m\in[n,N]\) with the same Gamma normalizations. The
[general first-coefficient derivative estimate](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#the-endpoint-cluster-region)
gives
\[
 \left|(\log B_m)'-1/m\right|
       \le \frac{5V+S(3+\log N)}{m^2}<\frac1n .               \tag{3}
\]
For the Gaussian factor this derivative gains just
\(-\tau/(2m^2)\), absorbed by replacing the ordinary \(V,S\) by
the augmented quantities. For the last inequality use
\(S<2\log N+5,\ \log N>12\), and \(V\le S(S+2)\).
The resulting numerator is at most \(40(\log N)^2\).
At \(m=n\), (1) bounds \(B_n\) by \(D e^{L_*/n}\).
Integrating (3) now gives \(B_m<101D/100\), since (2) and
\(D>1/10\) imply \(L_*/n<1/1800\).

The logarithmic interpolation of \(c_{m,k}k!/B_m^k\)
uses the usual finite differences. The Gaussian term adds
\(\tau k(k-1)/(2m)\), contributing \(\tau/m\) to the second
finite difference and \(-\tau/m^2\) to its derivative.
Thus the
[endpoint interpolation bounds](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
apply with \(\alpha=2\max(1,\max_i|d_i|/a_i)/n\),
where \(a_i\) is the starting exponent in this orientation,
and \(\beta=8eW/n\). The change in the second difference
and its derivative is absorbed by the increased \(W\).
Then (2) gives \(\alpha K<1/20\), \(\beta K<1/2\), and
\(\beta B_m<1/300\). These are exactly the inputs of the
[paired endpoint cluster estimate](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#a-paired-estimate-for-the-two-long-endpoint-sums),
with the same bound \(|(\log B_m)'|\le2/n\).
Its proof uses only those bounds, interpolation, and the
integer coefficient tails. It gives the analytic logarithm
estimates on the disk \(6/5\), for this orientation,
\[
 \|\log H_m\|\le M:=3(D+1),\qquad
                            \|\log H_N-\log H_n\|<1/55 .      \tag{4}
\]
The use of (3) replaces any need for monotonicity of
\(B_m/m\) in an asymmetric profile.

If the half has a private starting row of gap \(v\ge1\),
its normalized factorial factor is at most
\(f_{m,k}^v\), with \(f_{m,k}=k!/(m!)^{k/m}\).
Use \(\binom{m}{k}f_{m,k}\le e^k\), bound its remaining
\(v-1\) factors by curvature, and use the selected mirrored
pair separately. Since \(s_0\) includes the private row's
curvature \(v\), this gives
\[
 c_{m,k}^{1/k}\le e
    \left(\frac{\ell+dw}{\ell+d}\right)^d
                       e^{-(s_0-1)(1-w)/2}.
\]
At radius \(6/5\), the ratio is below
\((6/5)e(13/48)^7<1/2000\) on \(w\le1/6\), and below
\((6/5)e(9/16)^7<1/10\) for \(1/6<w\le1/2\).
The latter terms start beyond \(m/6\) with \(m>200000\).
The two geometric sums give \(\|H_m-1\|<1/200\), hence
\(\|\log H_m\|<1/190\). Both orders have the same private
starting row in this orientation, so these bounds give (4)
also in this case.

## The two radii do not affect the phase contradiction

By (4) on both halves, Rouché at radius
\(\exp(7(D+1)/m)\), then for the reciprocal, gives
\(|u_m|\le7(D+1)/m\). By (2) the moduli of the half arguments
are at most \(21/20\), since \(7(D+1)/n<77/720000\).
Use the notation \(J_m(u,\theta)\) for the difference of the
two half logarithms in the root phase equation. Cauchy's
bounds on the two nested disks, as in the
[absolute curvature phase comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#paired-radius-and-phase-comparison),
give, with \(\delta=1/55\),
\[
 |J_N-J_n|\le2\delta,\qquad
 |\partial_\theta(J_N-J_n)|\le14\delta,\qquad
 |\partial_u J_m|\le14M,\qquad
 |\partial_u\partial_\theta J_m|\le210M .
\]
At a common upper argument, the real parts of the two exact
root equations give
\[
  (n-14M)|u_N-u_n|\le7(D+1)/n+2\delta,\qquad
                               n|u_N-u_n|<1/25 .
\]
For the last inequality use
\(7(D+1)/n<77/720000\) and \(M/n<33/720000\).
The logarithms are real at \(\theta=0\) and real \(u\).
Integrating the derivative bounds therefore gives
\[
 |\Im(J_N(u_N,\theta)-J_n(u_n,\theta))|
 \le (14\delta+210M|u_N-u_n|)\theta<\theta/3 .
\]
Subtracting the two phase equations at the common argument
again contradicts the integer multiple of \(2\pi\).

Terminal two separates because the ordinary quadratic moment
is strictly positive on the real axis for these distinct
monomials. A Gaussian deformation preserves this positivity
by positive lognormal averaging. At \(\tau=0\), the absence
of a common complex root gives a rational Bézout identity
and thus ordinary coprimality over every characteristic-zero
field.

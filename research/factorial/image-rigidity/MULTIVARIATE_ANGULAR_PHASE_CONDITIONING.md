# Angular phase conditioning away from the flat-limit ray

The [matched tanh kernel transport](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md)
bounds smooth signed and unsigned kernel differences at a common
location. The [fractional grid transport](MULTIVARIATE_MATCHED_FRACTIONAL_LOG_TRANSPORT.md)
also bounds the assembled normalized log coefficients on one
continuous index interval. A direct stability estimate from these
coefficient norms would need control of the root phase denominator.
That denominator can be exponentially small relative to the
coefficient mass, even at angles a fixed distance from the negative
real axis.

## Phase derivatives for reciprocal references

Let \(G(z)=\sum_{k=0}^j a_k z^k\) be a real reciprocal polynomial
with \(a_k>0\). Let \(t=e^{i\varphi}\) be a simple root and write
\(v_k=k-j/2\). For a real log perturbation \(D_k\), use coefficients
\(a_k e^{sD_k}\) and follow the root near \(t\) as \(s\) varies.
The centered angular denominator is
\[
 H=\sum_k a_k v_k\sin(v_k\varphi),\qquad |H|=|G'(t)|.
\]
Indeed, at the root the derivative of the centered polynomial with
respect to \(\log z\) is \(t^{-j/2}tG'(t)=iH\).

For an even direction \(D_{j-k}=D_k\), the simple root stays on the
unit circle locally by reciprocal conjugation, and its argument has
derivative
\[
         \varphi'(0)
              =\frac{\sum_k a_kD_k\cos(v_k\varphi)}{H}.
\]
For an odd direction \(D_{j-k}=-D_k\), put
\(\gamma=\sum_k a_kD_k\sin(v_k\varphi)/H\).
Implicit differentiation in \(\log z\) gives
\[
 \frac{d}{ds}\log|t(s)|\bigg|_{s=0}=-\gamma,\qquad
 \varphi'(0)=0,\qquad
 \varphi''(0)
     =\frac{\sum_k a_k(D_k-\gamma v_k)^2\cos(v_k\varphi)}{H}.
\]
In the second differentiation the numerator before division by
\(iH\) is the sum of the centered squared direction. Its imaginary
part vanishes by reflection. These formulas hold at any simple
unit-circle root.

## An exponential phase sensitivity at moderate angles

Fix \(c=1/20\) and consider the
[Gaussian reference polynomials](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
\[
                  G_j(z)=\sum_{k=0}^j
                         \binom jk e^{-(c/j)k(j-k)}z^k.
\]
Their roots are simple and on the unit circle. Write them as
\(t_i=-e^{i\theta_i}\), with \(\theta_i\in(-\pi,\pi]\).
The first two Newton sums give
\[
 \begin{aligned}
 \frac1j\sum_i\cos\theta_i
     &= e^{-c(j-1)/j},\\
 \frac1j\sum_i\cos(2\theta_i)
     &= j e^{-2c(j-1)/j}-(j-1)e^{-2c(j-2)/j}.
 \end{aligned}
\]
The second expression tends to \((1-2c)e^{-2c}\).
For \(x_i=1-\cos\theta_i\), the mean and the second moment
therefore tend to
\[
       m=1-e^{-c},\qquad
       s=\frac32-2e^{-c}+\frac12(1-2c)e^{-2c}.
\]
The elementary exponential bounds give \(m\ge c-c^2/2\) and
\(s\le3c^2\). In particular \(m-s>1/40\).
For all sufficiently large \(j\) at least one root has
\[
                              1/40<x_i<1 .
\]
Otherwise \(x_i\le1/40\) or \(x_i\ge1\) for every \(i\), so
the mean of \(x_i\) is at most \(1/40\) plus its second moment.

Fix such a root \(t\). Its displacement from the negative real
axis is bounded below by \(\arccos(39/40)>0\) and above by
\(\pi/2\). Since \(G_j\) is monic, the arithmetic and geometric
mean inequality on squared chord lengths gives
\[
 |G_j'(t)|\le
 \left(\frac{2j}{j-1}\right)^{(j-1)/2}.
\]
Indeed the sum of the squared chord lengths to all roots is
\(2j(1-e^{-c(j-1)/j}\cos\theta_i)\le2j\). On the other hand
\(G_j(1)\ge2^j e^{-cj/4}\).

Take \(D_k=\cos((k-j/2)\arg t)\) at the interior coefficients
and set both endpoint entries to zero. This is an even log direction
of sup norm at most one. The Gaussian coefficient sequence is
unimodal, since
\[
 \log\frac{j-k}{k+1}\ \ge\ \frac{2(j-2k-1)}{j+1}
\]
on the increasing half. Each of the two closest adjacent central
coefficients is at least \(G_j(1)/(2(j+1))\).
For the central pair the sum of its two squared cosines is at least
\(1-|\cos(\arg t)|>1/40\). Therefore the phase derivative satisfies
\[
 |\varphi'(0)|\ \ge\
 \frac{G_j(1)}{80(j+1)|G_j'(t)|}
 \ \ge\
 \frac{\sqrt{2/e}}{80(j+1)}
         \exp\!\left(j\left(\frac{\log 2}{2}-\frac{c}{4}\right)\right).
\]
The perturbation preserves positivity, reciprocity and endpoint
normalization. Thus coefficient sup norms alone have no uniform
phase stability bound at these angles. An application of matched
transport must control the alternating numerator and denominator
together for the specific factorial directions.

## A common root with linked outer losses and strict log convexity

Use the same \(c=1/20\) at consecutive orders \(n,N=n+1\), with
Gaussian temperatures \(c/n\) and \(c/N\). The reflected outer losses
are \(U=c(n-1)/n\) and \(W=cn/(n+1)\), so
\[
                         W/U=\frac{n^2}{n^2-1}.
\]
For every sufficiently large \(n\), the Newton sum argument above
supplies an upper-order root \(t=-e^{i\theta}\) with
\(1/40<1-\cos\theta<1\). Applying the chord-square bound to all
\(n\) lower-order roots gives \(|G_n(t)|\le2^{n/2}\).

Among the two central adjacent indices of \(G_n\), choose an index
\(k\) satisfying
\(|\cos((k-n/2)\arg t)|\ge1/\sqrt{80}\). For sufficiently large
\(n\), both \(k\) and \(n-k\) are interior and differ from \(1\) and
\(n-1\). Write \(B(z)=z^k+z^{n-k}\), using only \(z^k\) if the
indices coincide. Both \(t^{-n/2}G_n(t)\) and \(t^{-n/2}B(t)\) are
real, and the latter has modulus at least \(1/\sqrt{80}\).
The real number
\[
                         \delta=-G_n(t)/B(t)
\]
therefore makes \(\widetilde G_n=G_n+\delta B\) vanish at \(t\).
The changed coefficients each have initial weight at least
\(G_n(1)/(2(n+1))\). Their relative changes are consequently bounded
by
\[
   2\sqrt{80}\,(n+1)
         \exp\!\left[-n\left(\frac{\log 2}{2}-\frac c4\right)\right].
\]
For sufficiently large \(n\) this preserves positivity. It also
preserves reciprocity, endpoint normalization and both outer
coefficients, hence the exact relation between \(U\) and \(W\).

Strict log convexity after removing the binomial factors holds
simultaneously at both orders. The second differences of the
original binomial-stripped logs are \(2c/n\) and \(2c/N\).
If \(\varepsilon_n\) is the sup norm of the change in the lower logs,
the change in each second difference has modulus at most
\(4\varepsilon_n\). Here
\(\varepsilon_n=O(n\exp[-n(\log 2/2-c/4)])\), so
\(2c/n-4\varepsilon_n>0\) for all sufficiently large \(n\).
This gives an exact common root in an envelope of arbitrarily
small polynomial log-coefficient width while keeping these
paired constraints. A transport argument for factorial profiles
must therefore retain additional joint information about the
coefficient directions.

## A bound for polynomial directions

Suppose all reference roots are simple and on the unit circle, with
angles \(\varphi_\ell\). At the root with angle \(\varphi_i\) put
\[
 L_i=1+\sum_{\ell\ne i}
         \bigl(1+|\cot((\varphi_i-\varphi_\ell)/2)|\bigr).
\]
For \(F(\varphi)=e^{-ij\varphi/2}G(e^{i\varphi})\), the common
root product cancels in the exact identity
\[
 \frac{F(\varphi_i+x)}{F'(\varphi_i)}
 =2\sin(x/2)\!\prod_{\ell\ne i}
     \bigl(\cos(x/2)+\cot((\varphi_i-\varphi_\ell)/2)\sin(x/2)\bigr).
\]
On \(|x|=1/L_i\), its modulus is at most \(e^{1/2}/L_i\).
To see this, put \(r=|x|\) and use
\(|2\sin(x/2)|\le re^{r/2}\),
\(|\cos(x/2)|\le e^{r/2}\), and
\(1+(r/2)|\cot\alpha|\le e^{(r/2)|\cot\alpha|}\).
Cauchy's estimate then gives
\[
        \left|\frac{F^{(m)}(\varphi_i)}{F'(\varphi_i)}\right|
                \le e^{1/2}m! L_i^{m-1},\qquad m\ge1.
\]
In particular, for a polynomial direction \(D(v)=\sum_m d_m v^m\),
\[
 \left|\frac{\sum_k a_kD(v_k)t_i^k}{t_i G'(t_i)}\right|
          \le e^{1/2}\sum_{m\ge1}|d_m|m! L_i^{m-1}.
\]
The same estimate controls the odd-direction phase formula by
applying it to the squared, radially recentered polynomial direction.
For fixed degree directions, the bound is polynomial whenever the
angular gaps have an inverse-polynomial lower bound. Applying it to
the full transported factorial directions still requires controlling
their higher weights relative to the angular gap.

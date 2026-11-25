# All-degree Curie–Weiss interlacing for \(s\ge24\)

For \(s>0\) let
\[
 B_j(z;s)=\sum_{k=0}^j\binom jk e^{-s k(j-k)/(2j)}z^k .
\]
For every \(s\ge24\) and \(j\ge1\), the roots of \(B_j\) and
\(B_{j+1}\) are simple and strictly interlace on the unit circle.
In the [Gaussian reference notation](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
these are the orders at the physical lower time ratio
\(\tau'/\tau=j/(j+1)\), with \(s=2j\tau\).
Outward monotonicity and the earlier inner comparison then give strict
interlacing at every time ratio in \([j/(j+1),1]\) when \(2j\tau\ge24\).

The [previous phase proof](MULTIVARIATE_DIFFUSE_LOW_TEMPERATURE_INTERLACING.md#the-phase-comparison)
gives the all-degree conclusion for \(s\ge64\).
The improvement below estimates the difference between the
integrand and its matching Gaussian also on the far horizontal
part of the contour.

## The saddle and the normalized half-integral

Put \(a=s/2\ge12\), and define
\[
 \begin{split}
 G(x,z)&=e^{-x^2/(4a)}(1+ze^{x-a}),\\
 I_j(z)&=\sqrt{\frac{j}{4\pi a}}\int_{-\infty}^a G(x,z)^j\,dx .
 \end{split}
\]
Expansion and integration over the full real line give \(B_j\).
The substitution \(x\mapsto2a-x\) on the right half yields
\(B_j(z)=I_j(z)+z^j I_j(1/z)\) for \(z\ne0\).

For \(|z|\le2\), put \(F(x,z)=-x^2/(4a)+\log(1+ze^{x-a})\),
using its logarithmic series in the disk \(|x|\le10^{-3}\).
The stationary equation is
\[
                 x=\frac{2az e^{x-a}}{1+ze^{x-a}}.
\]
On the disk the right side has modulus at most
\(4ae^{10^{-3}-a}/(1-2e^{10^{-3}-a})<10^{-3}\),
since \(ae^{-a}\) decreases and \(e>19/7\).
Rouché's theorem gives exactly one zero \(x_*(z)\) of the
stationary equation inside the disk, counted with multiplicity.
Thus the saddle is simple and depends analytically on \(z\) on
a neighborhood of \(|z|\le2\). It has conjugation symmetry.
Write \(F_*=F(x_*(z),z)\). Since
\(|ze^{x_*-a}|<1/50000\) and \(|x_*|<10^{-3}\),
the logarithmic series gives \(|F_*|<10^{-3}\).

Move the integration contour to the horizontal line through \(x_*\)
up to real part \(a\), followed by the vertical connector to \(a\).
This uses the entire function \(G(x,z)^j\) and Gaussian decay on
the left closing side. On the horizontal put \(x=x_*+t\), with
\(t\) real. We compare \(e^{-jF_*}I_j\) with the normalized integral
of \(g_0(t)^j=e^{-jt^2/(4a)}\) over the real line. This latter
integral is exactly one.

## Bounds on the four pieces

On the horizontal part with \(\Re x\le a/2\), use the
logarithmic series, now valid throughout that half-plane. There
\[
 \left|F_{xx}(x,z)+\frac1{2a}\right|
       \le\frac{2e^{-a/2}}{(1-2e^{-a/2})^2}
       \le\frac{41}{20}e^{-a/2}.
\]
Since \(ae^{-a/2}\) decreases for \(a\ge12\) and
\(24(41/20)(7/19)^6<1/8\), stationarity and integration give
\[
 F(x_*+t,z)-F_*=-t^2/(4a)+E(t),\qquad
 |E(t)|\le \frac18\,\frac{t^2}{4a}.
\]
The inequality \(|e^w-1|\le |w|e^{|w|}\) gives a normalized
error bounded by the full real Gaussian integral with integrand
\(\delta jt^2/(4a)\exp(-(1-\delta)jt^2/(4a))\), where
\(\delta=1/8\). Its value after normalization is
\(\delta/[2(1-\delta)^{3/2}]<2/25\).

For \(a/2\le\Re x\le a\) set \(g(t)=G(x_*+t,z)e^{-F_*}\).
The exact ratio is
\[
       g(t)=g_0(t)L(t)(1+ze^{x_*+t-a}),\qquad
       L(t)=\frac{e^{-x_*t/(2a)}}{1+ze^{x_*-a}}.
\]
Here \(|x_*t/(2a)|<10^{-3}\), hence \(|L-1|<1/500\)
by the exponential series and \(|ze^{x_*-a}|<1/50000\).
Writing \(u=\Re x\), we have \(t=u-\Re x_*\) and
\[
 g_0(t)\le (1001/1000)e^{-a/16},\qquad
 g_0(t)|ze^{x_*+t-a}|\le 2e^{1/2000}e^{-a/4}.
\]
For the second bound, the exponent
\(u-a-(u-\Re x_*)^2/(4a)\) increases up to \(u=a\),
where it is at most \(-a/4+1/2000\). Consequently
\[
 |g-g_0|\le D(a):=\frac{21}{10000}e^{-a/16}
                         +\frac{201}{100}e^{-a/4},\qquad
 \max(|g|,|g_0|)<\frac{29}{50}.
\]
The latter uses \(e^{-3/4}<19/40\) and
\(e^{-3}<(7/19)^3<51/1000\), at \(a=12\).
Use the difference of powers and the length \(a/2\).
The normalized contribution of this horizontal part is at most
\[
 \frac{\sqrt a}{4\sqrt\pi}\,
            j^{3/2}(29/50)^{j-1}D(a)
                   <\frac{23}{250}.
\]
Indeed \(j^{3/2}(29/50)^{j-1}\le7/4\), by squaring to check
\(j=1,2,3\) and using the decreasing ratio thereafter.
The product \(\sqrt aD(a)\) decreases for \(a\ge12\).
Finally \(\sqrt{12}/(4\sqrt\pi)<1/2\) and \(D(12)<13/125\).

The Gaussian completion beyond \(t=a-\Re x_*\) has normalized
mass at most \(e^{-b^2}/(2b\sqrt\pi)\), where
\(b=(a-\Re x_*)\sqrt{j/(4a)}>17/10\). This mass is less than
\(1/75\), using \(\sqrt\pi>5/3\) and
\(e^{289/100}>(19/7)^2(1+89/100)\).

On the connector, \(\Re x=a\) and \(|\Im x|<10^{-3}\), so
\[
 \left|G(x,z)e^{-F_*}\right|
      \le 3e^{10^{-3}+10^{-6}/(4a)}e^{-a/4}<4/25.
\]
The connector has length at most \(10^{-3}\).
After normalization its cost is less than \(10^{-3}\), since
\(\sqrt j(4/25)^j<1\).

## Strict alternating phases

We have proved
\[
 I_j(z)=e^{jF_*}(1+\eta_j(z)),\qquad
 |\eta_j(z)|<\frac{2}{25}+\frac{23}{250}+\frac{1}{75}
                         +\frac{1}{1000}<\frac{19}{100}.
\]
Set \(H_j=jF_*+\log(1+\eta_j)\), using the logarithmic power
series. This function is analytic on \(|z|\le2\) and has conjugation
symmetry. On this disk \(|\log(1+\eta_j)|\le19/81\). Cauchy's
estimate on radius-one disks gives for \(|z|\le1\)
\[
        |H_j'|<\frac{j}{1000}+\frac{19}{81}<\frac j4,\qquad
        |H_{j+1}'-H_j'|<\frac1{1000}+\frac{38}{81}<\frac12 .
\]
For \(z=e^{i\theta}\), the two terms of the exact decomposition
cancel just when \(T_j(\theta)=j\theta-2\Im H_j(e^{i\theta})\)
is an odd multiple of \(\pi\). The first derivative bound makes
\(T_j\) strictly increasing from zero to \(2j\pi\).
The difference \(\Im(H_{j+1}-H_j)(e^{i\theta})\) vanishes at
\(\theta=0,2\pi\), and its absolute value is strictly less than
\(\tfrac12\min(\theta,2\pi-\theta)\) in between. Thus
\(0<T_{j+1}(\theta)-T_j(\theta)<2\pi\) for \(0<\theta<2\pi\).
This places exactly one upper-order root before and after the
extreme lower-order roots and between each pair of successive
lower-order roots. The roots account for both degrees and are
simple.

The elementary numerical comparisons can use
\(e>163/60>19/7\), the Taylor sum through fourth order at
\(3/4\), and \(e^y\le1/(1-y)\) for \(0\le y<1\).
The bound is independent of \(j\) and concerns the Gaussian
reference family.


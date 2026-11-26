# Absolute signed imbalance for two private unit rows

Consider an ordinary two-monomial profile with one private unit
row on each side, arbitrary shared rows, and signed degree
difference \(\Delta=\sum_i(b_i-a_i)\).

**Theorem.** If \(|\Delta|\ge1647\), every ordinary consecutive
pair is coprime. The bound is absolute and uses no eventual
terminal cutoff.

The proof traps both sets of endpoint-normalized root radii
in order-dependent annuli. The radius control is uniform in the
dimension and in the shared curvature below the range already
covered by the angular theorem.

## Endpoint means and concentration

Write \(V=\sum_{\rm shared}d_i^2/\ell_i\) and
\(S=\sum_{\rm shared}d_i^2/(\ell_i+d_i)\), where \(d_i=|b_i-a_i|\)
and \(\ell_i=\min(a_i,b_i)\) on unequal rows. The
[private-unit curvature theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
covers \(S\ge15/8\). For \(0<S<15/8\), each \(d_i/\ell_i\le2\):
a larger ratio requires \(d_i\ge3\) and contributes
\(d_i(d_i/\ell_i)/(1+d_i/\ell_i)>2\) to \(S\).
It follows that \(V\le3S<6\).

The private unit factorials cancel the binomial coefficient.
In degree \(m\ge4\), normalize the endpoint coefficients to
one and write the polynomial as
\[
 P_m(z)=\sum_{k=0}^m \mathbb E_\mu T^k\,z^k,\qquad T>0 .
\]
Here \(\mu\) is the geometric positive moment law obtained by
the independent Gamma row integrals and the endpoint scaling.
Write \(f(x)=\log\mathbb E_\mu T^x\) on \(0\le x\le m\).
The endpoints are zero. The upper trigamma estimate gives
\(f''\le V/m\), so \(|f'|\le V/2\).
For the other endpoint law \(\mu^\vee\), obtained by
\(T^m\) tilt followed by inversion, the corresponding log
Mellin function is \(g(x)=f(m-x)\), with the same bounds.
[Geometric positive moment representation](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)

There is a stronger estimate at the left endpoints of these
two Mellin functions. Indeed
\(-f'(0)=\int_0^m (1-x/m)f''(x)\,dx\).
For each row, the lower trigamma estimate gives
\[
 -f_i'(0)\ge
 d_i^2\int_0^1\frac{1-s}{a_i+(b_i-a_i)s+1/m}\,ds
              \ge \frac{d_i^2}{2(\ell_i+2d_i/3+1/m)} .
\]
The density \(2(1-s)\) on \([0,1]\) has mean \(1/3\);
apply Jensen to the reciprocal denominator and use both
possible orientations. The denominator on the right is at
most \(31\ell_i/6\) when \(d_i/\ell_i\le2\) and \(m\ge4\).
Thus \(f'(0),g'(0)\le-6V/31\le-3V/16\).

For either law, put \(B=\log T\). Its mean lies in
\([-V/2,-3V/16]\) and its variance is at most \(V/m\).
Taylor's integral formula and \(f'',g''\le V/m\) also give
\[
     \mathbb E e^{s(B-\mathbb E B)}
                       \le e^{s^2 V/(2m)}\qquad(0\le s\le m).
 \tag{1}
\]
These estimates will be applied separately to the two laws,
so their probability consequences also hold for any mixture.

## Endpoint packets on the two radial boundaries

Suppose \(mV\ge1536\). Fix \(r=e^{12/m}\) or
\(r=e^{-12/m}\) and put \(X=B+\log r\).
Then \(|\log r|\le V/32\).
By Chebyshev's inequality the event
\(|B-\mathbb E B|\le V/16\) has probability at least \(3/4\).
On that event
\[
                           -19V/32\le X\le-V/16 .
 \tag{2}
\]
Also \(\mathbb E X\le-5V/32\). Taking \(s=m/8\) in
(1) gives
\[
              \mathbb P(X>-V/32)\le e^{-mV/128}.
 \tag{3}
\]
Since \(u e^{-u/128}\) is decreasing for \(u\ge1536\),
we have \(mV e^{-mV/128}<1/64\). At the endpoint this
follows from \((8/3)^{12}>1536\cdot64\) and \(e>8/3\).

For a reference angle \(\pi/m\le\theta\le\pi\), define
\[
                 E=\mathbb E(1-e^{X+i\theta})^{-1}.
\]
For either endpoint law, and for their mixtures, we claim
\[
                 \frac1{512(V+\theta)}
                           \le |E|\le\frac{67}{V+\theta}.
 \tag{4}
\]
For \(0<\theta\le\pi/2\), the kernel has modulus at
most \(1/\sin\theta\), by minimizing the distance
\(|1-u e^{i\theta}|\) over \(u\ge0\).
At the angles in question this is at most \(m/2\).
Thus the packet from the event in (3) has modulus at most
\(1/(128V)\). For the upper bound on \(\theta\le V\) and
\(\theta\le\pi/2\), use \(X\le-V/32\) on the remaining event.
The modulus of its kernel is at most both \(64/V\)
and \(2/\theta\), hence at most \(66/(V+\theta)\).
Adding the tail gives the stated upper bound. For
\(V\le\theta\le\pi/2\) the bound \(2/\theta\) alone
suffices. For \(\theta\ge\pi/2\), the kernel has modulus
at most one, and \(V+\theta<10\).

For the lower bound when \(\theta\le V,\ \theta\le\pi/2\),
the kernel has nonnegative real part whenever \(X\le0\).
On (2), this real part is at least
\[
  \frac{1-e^{X}}{(1-e^{X})^2+2e^X(1-\cos\theta)}
     \ \ge\ \frac{V/32}{V^2+\theta^2}\ \ge\ \frac1{64V}.
\]
Here \(1-e^{-V/16}\ge V/32\) since \(V<6\).
The possible negative contribution from \(X>0\)
costs at most \(1/(128V)\) by (3).
This leaves \(\operatorname{Re}E\ge1/(256V)\).
For \(V\le\theta\le\pi/2\), the imaginary part is positive
everywhere. On (2), use \(e^X>1/64\),
\(\sin\theta\ge\theta/2\), and the denominator at most
\(V^2+\theta^2\le2\theta^2\). This gives
\(\operatorname{Im}E\ge1/(512\theta)\).
For \(\theta\ge\pi/2\), the real part is nonnegative
everywhere and on (2) it is at least \(1/4\);
the resulting \(3/16\) lower bound suffices.
All inequalities were either pointwise or applied
separately to the two component laws, so they hold for the
mixtures as asserted.

At \(z=r e^{i\theta}\) the finite geometric sum splits
the normalized polynomial exactly as
\[
 P_m(z)=\mathbb E_\mu(1-zT)^{-1}
                   +z^m\mathbb E_{\mu^\vee}(1-z^{-1}T)^{-1}.
 \tag{5}
\]
The two moduli of the endpoint expectations are bounded
as in (4), allowing \(r\) and \(r^{-1}\).
Their ratio lies between \(1/(67\cdot512)\) and
\(67\cdot512\). Since \(e^{12}>67\cdot512\), (5) cannot
vanish at \(r=e^{\pm12/m}\).

Apply this radial boundary exclusion along the linear
interpolation of the moment law \(\mu\) to the reciprocal
mixture \((\mu+\mu^\vee)/2\). Its reversal is the reciprocal
interpolation of the two laws, so (4) continues to hold in
both terms of (5). Positive polynomial coefficients exclude
arguments of absolute value below \(\pi/m\), since there the
terms have arguments spanning less than a half-circle.
The reciprocal mixture has all simple roots on the unit
circle by the
[geometric argument boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
Continuity of the polynomial roots, together with the radial
boundary exclusion, puts all actual roots between the same
two boundaries.

When \(mV\le1536\), the consecutive coefficient ratios are
between \(e^{-V/2}\) and \(e^{V/2}\) by \(|f'|\le V/2\).
Eneström–Kakeya gives \(|\log|z||\le V/2\le768/m\).
Combining the two regimes gives the bound \(768/m\)
at every \(m\ge4\).

## Physical drift

Let \(g_m\) be the positive endpoint-normalizing argument
scale of the physical order-\(m\) polynomial. The
[linked center estimate](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
for two private units gives
\[
              |\log(g_N/g_n)-\Delta\log(N/n)|\le V/(nN).
\]
For an actual common physical root, the radial bound implies
\[
     |\Delta|\log(N/n)\le768(1/n+1/N)+V/(nN).
\]
For \(n\ge7\) use \(\log(N/n)>1/N\), \(N/n\le8/7\) and
\(V<6\). This gives
\(|\Delta|<768\cdot15/7+6/7<1647\).
The
[terminal-seven theorem and its predecessors](MULTIVARIATE_TERMINAL_SEVEN_ALL_PROFILES.md)
cover \(N\le7\).
At zero shared curvature the signed gap is zero.
This proves coprimality at the stated imbalances, also over
every characteristic-zero field by a rational Bezout identity.

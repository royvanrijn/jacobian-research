# Density of isolated pockets in a fixed Gaussian strip

Fix \(0<\rho<1\) and an ordinary stationary Gaussian winner
partition \(A\) in a fixed finite dimension. Fix a host label
\(j\), a unit direction \(v\), and positive constants \(D,\eta,b,M\).
Consider points \(p_n=s_nv+y_n\), where \(y_n\perp v\),
\(|y_n|\le M\) and \(s_n\ge2\). Suppose for each \(n\) a
measurable set \(S_n\) lies in \(B_D(p_n)\), is a subset of
the nonhost labels \(A_j^c\), and has Lebesgue measure at
least \(b\). Assume its essential distance from
\(A_j^c\setminus S_n\) is at least \(\eta\).
If the \(s_n\) have any fixed positive separation, then
\[
                              \sum_n\frac1{s_n}<\infty.        \tag{PS1}
\]
In particular, pockets of uniformly positive physical size
and fixed positive clearance cannot recur with positive
linear density along a fixed strip. Their shapes can vary
arbitrarily.

The result uses the horizontal derivative integrability for
stationary winner partitions. In coordinates \((s,y)\) along
\(v\) and its orthogonal complement, let
\(\sigma=\sqrt{1-\rho^2}\), \(H^s_\sigma\) be centered
one-dimensional Gaussian convolution of standard deviation
\(\sigma\), and \(T^y_\rho\) be the Gaussian noise operator
in the transverse block. For \(f=1_{A_j^c}\), (PT8) gives
\[
 \int_1^\infty\frac1s
       \int_{\mathbb R^{d-1}}
          \big|\partial_s(H^s_\sigma T^y_\rho f)(s,y)\big|^2
                         d\gamma_{d-1}(y)\,ds<\infty.          \tag{PS2}
\]
The case of a scalar host complement follows from the
coordinatewise bound, or \(f=1-1_{A_j}\).
[Horizontal dissipation in a fixed strip, PT8](GAUSSIAN_PERIODIC_INTERFACE_TAIL_OBSTRUCTION.md#horizontal-dissipation-in-a-fixed-strip)

## A uniform derivative cost for a pocket

There is a fixed \(c>0\), depending on the specified
geometric constants, the dimension and \(\rho\), such that
\[
 \int_{s_n-1}^{s_n+1}\int_{B_1^{d-1}}
          \big|\partial_s(H^s_\sigma T^y_\rho f)(s,y)\big|^2
                         d\gamma_{d-1}(y)\,ds\ge c .           \tag{PS3}
\]
It suffices to prove this bound over all bounded nonhost
indicators and all isolated pockets with the stated geometry.
No maximizing equation is needed for this local lower bound.

Translate only in the horizontal direction by each \(s_n\).
If the lower bound failed, use a sequence of such indicators
\(f_n\) with energy tending to zero on
\([-1,1]\times B_1^{d-1}\). The corresponding translated
pocket indicators \(g_n\) satisfy \(0\le g_n\le f_n\le1\)
and have support in a fixed compact ball and integral at
least \(b\). Pass to locally weak-* limits \(g,f_*\) and
a Hausdorff limit \(K\) of the closed essential supports
\(K_n\) of \(g_n\). The compactness of these nonempty
compact sets takes place in that same fixed ball. Then
\(0\le g\le f_*\), \(\int g\ge b\), and \(g\) is
supported on \(K\).

In the open collar \(0<\operatorname{dist}(z,K)<\eta/2\),
the limit \(f_*\) vanishes almost everywhere. For each
compact subset of this collar, Hausdorff convergence
places it outside \(K_n\) and at distance less than
\(\eta\) from \(K_n\) for large \(n\). The pocket indicator
vanishes there, and the remaining nonhost material
cannot enter by the separation assumption. Pass to the
limit and exhaust the collar by compact sets.

Gaussian smoothing and its horizontal derivative pass
to the locally weak-* limit uniformly on compact
observation sets. To verify this, first discard a fixed
source tail using uniform integrability of the Gaussian
derivative kernels for bounded observations, and then
use their compactness in \(L^1\) on the bounded source
set. Thus the smoothed horizontal derivative of \(f_*\)
vanishes on the window, and hence throughout the space
by real analyticity. The smoothing here is ordinary
isotropic Gaussian convolution in the source variables
evaluated at \((s,\rho y)\). For every horizontal
translation, the two smoothed translates consequently
agree. Fourier division by the nonvanishing Gaussian
multiplier on compact frequency tests gives equality of
the two original translates as tempered distributions.
Equivalently \(f_*\) is independent of \(s\) almost
everywhere, as can be seen by mollifying and passing to
local \(L^1\).

This invariance contradicts the collar. Indeed \(g\) has
positive integral on \(K\), so on a positive-measure set
of transverse coordinates \(y\), the section \(K_y\)
is nonempty and the constant value \(f_*(y)\) is
positive. For almost every such \(y\), the compact section
\(K_y\) has a right endpoint \(a_y\). Each point
\((s,y)\) with \(a_y<s<a_y+\eta/2\) lies outside \(K\)
and in its open collar. Fubini gives zero \(f_*\) there,
contradicting the constant positive value on that fiber.
This establishes (PS3).

Finally the fixed horizontal windows in (PS3) have
bounded overlap by the assumed separation of the \(s_n\).
On each window \(1/s\ge 1/(s_n+1)\). Summing (PS3)
after weighting by \(1/(s_n+1)\), and applying (PS2),
gives (PS1). For example, the number of these centers
in \([R,2R]\) is \(o(R)\), and the normalized counts
are summable across dyadic intervals. The estimate
allows a sequence of pockets at geometric distances.

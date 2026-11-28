# Cumulative Gaussian order of the spectator projection

For a balanced maximizing compact discrepancy from a centered
three-cell simplex, the full reference score lines give one
spectator projection for each cell. A cumulative order for
these projections suffices to exclude a nonsimplex tie
approached from the right by simplex-optimal correlations.
The order is weaker than requiring one sign change of the
transverse source. A nonsimplex right-contact maximizer must
violate it by a uniform amount on any compact correlation
interval, independently of the radius of the discrepancy.

## An exact cumulative test

Write \(0<r<1\), \(\sigma^2=1-r^2\), and let
\(\gamma_{\sigma^2}\) denote the centered Gaussian law of
variance \(\sigma^2\). For spectator \(k\) write \(t,v\) for
the coordinates along \(\tau_k=-w_k\) and normal to it, where
\(w_k\) is the unit center of the reference sector \(T_k\).
Set \(q_k(t)=\int 1_{A_k}(t,v)\,d\gamma_{\sigma^2}(v)\) and
\[
 C_k(b)=\int_{-\infty}^b(1-3q_k(t))\,d\gamma_{\sigma^2}(t)
       =\gamma_{\sigma^2}((-\infty,b])
          -3\gamma_{\sigma^2 I_2}(A_k\cap\{t\le b\}).
\]
Here \(\gamma_{\sigma^2 I_2}\) is the planar centered law
with covariance \(\sigma^2 I_2\). The equal-price lines
give a triple score tie at the origin, so
\(\gamma_{\sigma^2 I_2}(A_k)=1/3\) and \(C_k(+\infty)=0\).
For \(s\ne 0\), set \(z=rs/\sigma^2\). Integration by parts
in the Gaussian convolution gives
\[
 \frac{2h_k(s)e^{r^2s^2/(2\sigma^2)}}{z}
      =-\int_{\mathbb R} e^{zt} C_k(t)\,dt .                \tag{CO1}
\]
There is no boundary term. Indeed \(C_k\) is absolutely
continuous and
\[
 |C_k(t)|\le 2\min\{\gamma_{\sigma^2}((-\infty,t]),
                         \gamma_{\sigma^2}([t,\infty))\},  \tag{CO2}
\]
using the two tail representations of the primitive and
\(|1-3q_k|\le 2\).

The condition \(C_k\le0\) is precisely first-order
stochastic dominance by the unconditional tilted
longitudinal coordinate over that coordinate conditional
on \(A_k\) under the tilted planar law. If \(C_k\le0\) and
\(C_k\not\equiv0\), then (CO1) gives
\(\operatorname{sign}h_k(s)=\operatorname{sign}s\) for every
\(s\ne0\). In particular the cumulative condition on all
three lines excludes both the antipodal germ and the
forward excursion required of a nonsimplex right contact.
For a compact discrepancy \(C_k=C_k^T<0\) sufficiently far
out on both sides, where \(C_k^T\) is the reference
primitive. For the equality to the right use the zero
total integral. The strict reference sign follows because
\(1-3q_k^T\) goes from negative to positive in one
crossing and has zero tilted integral.
[Central germ and the right-contact arm test](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#genuine-triple-points-on-the-reference-lines)

The same identity counts real score zeros with
multiplicity. If a nonzero \(C_k\) has \(m\) essential sign
changes, then \(h_k\) has at most \(m+1\) real zeros counted
with multiplicity. Indeed (CO2) and the Gaussian tail
bound make \(e^{t^2/(2\sigma^2)} C_k(t)\) bounded.
Apply the Laplace-transform zero count to this source.
The factor \(z\) in (CO1) adds one zero at the origin.
[The Laplace-transform zero count](THREE_CELL_TRANSVERSE_SOURCE_OSCILLATION.md#the-count-on-a-full-reference-line)

More generally, at any maximizer with the hypotheses of
(SO2), let \(\eta\) be the sign selecting its active germ
half-rays. If \(\eta C_k\le0\) on all three lines, then
(CO1) gives \(h_k(\eta u)>0\) for all \(u>0\), because the
nonzero germ excludes \(C_k\equiv0\). This contradicts
(SO2).
[Uniform germ and excursion](THREE_CELL_SPECTATOR_SOURCE_ORDER.md#a-uniform-excursion-on-one-of-the-three-germ-arms)

## A uniform positive cumulative violation

Fix a compact interval \(J\subset(0,1)\). Under the
hypotheses of (SO2) at \(r\in J\), there are \(M,\delta>0\),
depending only on \(J\), and some \(k\) and \(|b|\le M\)
with
\[
                              C_k(b)\ge\delta.             \tag{CO3}
\]
For an antipodal germ this holds for every \(k\), and the
derivative comparison in (SO2) is then unnecessary.

In the reference orientation, choose \(s>0\) and \(k\)
from (SO2). In the antipodal orientation, take \(s=s_0\)
in (SO1) for each \(k\). Thus in both cases \(s\) belongs
to a fixed compact subinterval of \((0,\infty)\), and
\(h_k(s)\le-m\) for a uniform \(m>0\). Formula (CO1)
bounds \(\int e^{zt} C_k(t)\,dt\) below by \(2m/Z\) where
\(Z>0\) is a uniform upper bound on \(z\). The integrand
in absolute value is dominated uniformly by (CO2) times
\(e^{Z|t|}\), which has Gaussian tails since \(\sigma\le
1\). Choose \(M>0\) so the absolute integral of this
bound outside \([-M,M]\) is less than \(m/Z\). The
integral of \(e^{zt}C_k(t)\) over \([-M,M]\) is then at
least \(m/Z\), proving (CO3), for example with
\(\delta=m/(2MZe^{ZM})\). The uniform Lipschitz bound
\(2/(\min_{r\in J}\sqrt{1-r^2}\sqrt{2\pi})\) for \(C_k\)
also gives an interval of uniform positive length on
which \(C_k\ge\delta/2\).

For a compact discrepancy, the strict negative
reference tails give at least two sign changes of \(C_k\)
on every line selected by (CO3). Thus the necessary
disorder persists after taking a Gaussian-weighted
primitive. A transfer that establishes the cumulative
order on the three lines would suffice for the
compact-defect exclusion.

The violation also supplies an explicit penalty in the
[pair-wall stability inequality](THREE_CELL_PAIR_ARC_CENTRAL_EXCLUSION.md#the-transverse-source-in-the-pair-wall-inequality).
Use \(P_s=N(rs,\sigma^2)\) on the line chosen in (CO3).
With \(B=(-\infty,b]\) and \(p=P_0(B)\), the covariance
identity and Cauchy–Schwarz give
\[
 |C_k(b)|=3\left|\int (q_k-P_s q_k)(1_B-p)\,dP_0\right|
  \le 3e^{r^2s^2/(2\sigma^2)}
                    \sqrt{\operatorname{Var}_{P_s}q_k}.
\]
Here \(\int(dP_0/dP_s)^2\,dP_s=e^{r^2s^2/\sigma^2}\).
In (PA10) the loss in normalized normal slope from the
spectator variance is consequently at least
\((\pi\delta^2/36)e^{-r^2s^2/\sigma^2}\), uniformly
positive on bounded observation windows. This penalty
applies to the second variation on the actual regular
pair-wall portions of the line.

## The cumulative condition is strictly weaker

The [shell construction](THREE_CELL_TRANSVERSE_SOURCE_OSCILLATION.md#ordinary-compact-examples-with-many-changes)
gives ordinary balanced compact discrepancies with all
three exact score lines for which \(C_k(t)<0\) at every
finite \(t\), while \(1-3q_k\) has arbitrarily many sign
changes on every line. To see the cumulative sign,
first use the fractional tuple before purification.
Conditional on any radius \(R>0\), let \(X\) be the
coordinate \(t\) under uniform angle. On an unchanged
circle, spectator membership is exactly
\(\{X\le -R/2\}\), of angular mass \(1/3\). The radial
contribution to \(C_k(b)\) is
\[
 \mathbb P(X\le b)-\mathbb P(X\le b\mid X\le -R/2)\le0,
\]
strict when \(|b|<R\). On a modified circle the uniform
fractional assignment gives zero. Integrate against the
tilted radial law, using unchanged circles at arbitrarily
large radii for strictness. Weighted-projection
purification retains \(q_k\) almost everywhere and hence
\(C_k\) everywhere. The purified examples choose the
reference simplex under their scores and fail the winner
equation for their own labels. Whether actual compact
maximality forces the cumulative order remains open.

## Forward spectator patches that violate the cumulative order

The compact-discrepancy geometry and the three exact pair
lines also allow cumulative violations on every line,
even when the scores still choose the reference simplex.
Fix \(r\in(0,1)\) and take \(R\) large. Set
\(d=2\sigma^2/R\), \(a=R-d\), and \(b=R+d\). In the
coordinates of each spectator \(k\), assign the rectangle
\(a\le t\le b,\ |v|\le4\sigma\) entirely to label \(k\).
Keep the reference partition outside the three rectangles.
They are disjoint for large \(R\), since their centers
are on the three forward rays and their diameters stay
bounded. Each rectangle lies in the complement of its
assigned label in \(T\), because \(t>0\) there.

The rotations and reflections of the reference simplex
permute the rectangles and their assigned labels in the
same way as the reference sectors. Rotation gives equal
ordinary and tilted Gaussian masses. Reflection in the
line for \(i,j\) exchanges the two labels, so their scores
agree on the whole line. Central and exterior assignments
are unchanged.

For a fixed spectator \(k\), the rectangle gives
\(q_k(t)\ge 15/16\) throughout \([a,b]\) by Chebyshev's
inequality. Write \(P(x)=\gamma_{\sigma^2}([x,\infty))\).
Shifting the Gaussian tail integral by \(2d\) gives
\(P(b)/P(a)\le e^{-2ad/\sigma^2}<1/8\) for all
sufficiently large \(R\). Consequently
\[
 -C_k(a)=\int_a^\infty(1-3q_k)\,d\gamma_{\sigma^2}
 \le P(a)-\frac{45}{16}\,[P(a)-P(b)]<0 .                \tag{CO4}
\]
Thus \(C_k(a)>0\) for all three spectators at a location
that can be arbitrarily far out.

The total rectangle area is \(48\sigma d\), tending to
zero. The conditional kernel
\((2\pi\sigma^2)^{-1}\exp(-|y-rx|^2/(2\sigma^2))\)
and its first observer derivatives are bounded uniformly
in \(x,y\) at fixed \(r\). Thus the change in each pair
score and its gradient is uniformly bounded by a constant
times the rectangle area. The pair score change vanishes
on its reference line. Integrating the gradient from that
line and combining with the uniform score bound gives a
bound by the same area times a constant times
\(\min\{\operatorname{dist}(x,L_{ij}),1\}\). The reference
score margin has a strictly positive coefficient in this
distance factor on \(T_i\) for every \(j\ne i\). Hence
the scores still choose \(T\) off its interfaces for large
\(R\). These partitions fail the winner equation on the
rectangles. Since the violations in (CO4) move out into
Gaussian tails, they do not give the bounded-window
violation required of a maximizer by (CO3). In fact
\(C_k\to C_k^T\) uniformly as the tilted Gaussian mass
of the changed rectangles tends to zero. On each
fixed compact interval the reference primitive has
a strictly negative maximum.
[The saturated reference score margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

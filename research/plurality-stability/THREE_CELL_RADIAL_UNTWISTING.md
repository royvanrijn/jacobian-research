# Removing a radial twist from two adjacent Gaussian cells

For a planar partition, rotate each circle centered at the
Gaussian origin by an angle that may depend measurably on its
radius. This preserves the three cell masses exactly. For the
following class, a common angular orientation maximizes the
interaction of every pair of circles at every \(0<\rho<1\).

On the circle of radius \(s>0\), suppose labels \(1,2\) occupy the
adjacent arcs \((\theta(s),\theta(s)+a(s))\) and
\((\theta(s)-a(s),\theta(s))\), where \(0\le a(s)\le\pi\).
Label \(3\) occupies the complement. The widths of labels \(1,2\)
agree at each radius and can vary arbitrarily with radius.
If all three labels have positive Gaussian mass, any
nonconstant \(\theta\) on the radii with \(a>0\) gives strictly
lower noise stability than a common angle, modulo changes on
null sets. In particular, a balanced global maximizer in this
class cannot have bounded cells \(1,2\).

## Alignment on a pair of circles

Conditional on positive radii \(s,t\), the joint angular
density is proportional to
\[
             k_\kappa(\varphi-\psi)=\exp(\kappa\cos(\varphi-\psi)),
          \qquad \kappa=\frac{\rho st}{1-\rho^2}>0 .
\]
This follows by writing the planar joint Gaussian density in
polar coordinates. The density of the two radii is strictly
positive.

For circular arcs of lengths \(a,b\), let \(C_{a,b}(u)\) be the
angular length of the intersection of two such arcs whose centers
are separated by \(u\). It is an even, \(2\pi\)-periodic function.
When \(a+b\le 2\pi\), its value on \([0,\pi]\) equals
\(\min(a,b)\) up to \(d=|a-b|/2\), then decreases with slope
\(-1\) to zero at \(h=(a+b)/2\), and vanishes beyond \(h\).
For \(a+b>2\pi\), passing to the complements changes this
function by an additive constant and reduces to the first case.

Put \(a=a(s)\), \(b=a(t)\). The sum of the two matched angular
overlaps for labels \(1,2\), at relative rotation \(u\), is
\[
           H_{a,b}(u)=C_{a,b}(u+d)+C_{a,b}(u-d).            \tag{RU1}
\]
This is even and nonincreasing on \([0,\pi]\). The second
summand is constant until \(u-d\) exits the plateau in the
positive direction, and is then nonincreasing. For the first
summand, the support near zero ends at \(h\) and no periodic
copy can contribute for arguments between \(d\) and \(\pi+d\),
since \(\pi+d\le 2\pi-h\) follows from \(\max(a,b)\le\pi\).
It too is nonincreasing. If \(a,b>0\), the value at zero is
\(2\min(a,b)>0\), whereas the value at \(\pi\) is zero.
The matched overlap for label \(3\) is
\(C_{2\pi-2a,2\pi-2b}(u)\), which is even and nonincreasing
by the same arc-overlap formula.

To evaluate the angular Gaussian pairing, convolve the
matched overlap sum with \(k_\kappa\). The convolution of
an even nonincreasing circular function with this kernel is
again even and nonincreasing, strictly decreasing on
\((0,\pi)\) if the function is nonconstant. Indeed, decompose
it into a constant and a positive integral of centered arc
indicators by its superlevel sets. For an arc of half-length
\(0<\ell<\pi\), the derivative of the convolved indicator at
\(0<u<\pi\) is \(k_\kappa(u+\ell)-k_\kappa(u-\ell)<0\), since
\(\cos(u+\ell)-\cos(u-\ell)=-2\sin u\sin\ell<0\).
The nonconstant contribution in (RU1) gives strictness if
both widths are positive.

The angular pairing for each pair of radii is therefore
maximized at equal angles. When both widths are positive
it is strictly maximized there modulo \(2\pi\).
Integrating against the positive joint density of the radii
gives the stated comparison. Equality implies equality of
\(\theta(s)\) and \(\theta(t)\) modulo \(2\pi\) for almost every
pair of radii with positive widths, and hence a common angle
on that set by Fubini.

## Infinitesimal rotation within a semicircle

There is a local consequence when the essential range of the
angles on positive-width radii lies in a closed semicircle.
Choose real lifts in an interval of length \(\pi\) and multiply
the angles by \(1-\epsilon\). Every pair difference in
\((-\pi,\pi)\setminus\{0\}\) strictly improves to first order
by the derivative calculation above, and all other pairs
have zero first derivative. If a set of pairs of positive
measure has such a difference, the total first derivative
is strictly positive. This excludes stationarity under
mass-preserving rotations of the circles.

The only remaining nonconstant case in a closed semicircle
uses exactly its two antipodal endpoints on sets of positive
radial measure. Here the first derivative vanishes. For an
arc of half-length \(0<\ell<\pi\), the second derivative of
its convolved indicator at \(u=\pi\) is
\[
   k_\kappa'(\pi+\ell)-k_\kappa'(\pi-\ell)
        = 2\kappa\sin\ell\, e^{-\kappa\cos\ell}>0 .
\]
The nonconstant overlap sum gives strictly positive second
derivative for pairs of opposite angles. Pairs of equal
angles have no change under this variation. Thus the total
second derivative of noise stability is positive.

The angle speeds used here are bounded and measurable.
They can be approximated almost everywhere by uniformly
bounded smooth radial speeds supported away from zero and
infinity. The first and second derivatives of the stability
converge, by boundedness of the speeds and Gaussian domination
of the polar kernel times \(1+\kappa^2\). Explicitly the
angular factor is bounded by \(e^\kappa\), so the remaining
joint radial Gaussian is bounded by a polynomial times
\(\exp[-(s^2+t^2)/(2(1+\rho))]\). Each smooth approximation
is an exact mass-preserving variation. The strict
first-order or second-order sign therefore holds for a
smooth mass-preserving variation as well.

## Consequence for winding compact networks

With a common angle, reflection about its line exchanges
labels \(1,2\) and fixes label \(3\). The reflection reduction
for a balanced global maximizer gives a finite analytic
height over the full half-line and unbounded cells of all
three labels. Thus bounded cells \(1,2\) in the stated
angular class are excluded. The proof allows arbitrary
winding in \(\theta\) and does not use a Cartesian graph
representation for the starting partition.
[Reflection reduction, RF1](THREE_CELL_REFLECTION_REDUCTION.md)

There is also a neighborhood exclusion for each specific
template \(A\) in this class with genuinely varying angle.
Apply its untwisting map \(U\) to any other balanced partition
\(P\). The map is invertible and Gaussian measure preserving,
so it preserves the sum of the \(L^1(\gamma_2)\) distances
between label indicators. For any two partitions,
\[
     |S_\rho(A)-S_\rho(P)|
       \le 2\sum_i\|1_{A_i}-1_{P_i}\|_{L^1(\gamma_2)},
\]
by expanding the quadratic form and using the Markov
contraction on bounded functions. If
\(\Delta=S_\rho(UA)-S_\rho(A)>0\), the same map gives
\[
   S_\rho(UP)-S_\rho(P)
          \ge \Delta-4\sum_i\|1_{A_i}-1_{P_i}\|_{L^1(\gamma_2)}.
                                                               \tag{RU2}
\]
This excludes all sufficiently close competitors, including
ones lacking equal adjacent angular arcs. No size-independent
neighborhood or conclusion about arbitrary compact networks
follows from this comparison.

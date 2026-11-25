# Eventual coprimality for fixed subcritical equal-degree supports

Let \(a,b\in\mathbb Z_{>0}^r\) be distinct exponent columns,
with equal rows removed, and let \(d_i=b_i-a_i\),
\(x_i(q)=a_i+d_iq\). Assume \(\sum_i d_i=0\) and
\[
 \max_{0<x<1}\rho(x)<1,qquad
 \rho(x)=x(1-x)\sum_i\frac{d_i^2}{x_i(x)}.
\]

**Theorem.** For every such fixed ordinary fully shared profile,
the moment polynomials \(M_n(t),M_{n+1}(t)\) are coprime for
all sufficiently large \(n\). The threshold depends on the
profile and is unevaluated. The result holds over every
characteristic-zero field.

Use the upper Rice phase and the stationary target
\(w=-\operatorname{Log}t\) from the
[exact contour](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases).
The proof classifies the equal-modulus arcs on which both
regular upper saddles contribute to that contour. This
classification allows the [exact phase matches with small curvature](MULTIVARIATE_SUBCRITICAL_SADDLE_PHASE_MATCH.md)
without requiring a separate comparison of their amplitudes.

## The stationary covering and the equal-modulus curves

Let \(\mathcal T=\{w:-\pi<\Im w<0\}\). The
[stationary map](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md#connectedness-and-ramification)
\[
 X(q)=\sum_i d_i\log x_i(q)+\log(1-q)-\log q
\]
identifies \(U=X^{-1}(\mathcal T)\cap\mathbb H\) as a connected
proper double cover with one simple ramification point \(q_c\).
Write \(w_c=X(q_c)\) and let \(\tau\) be its sheet involution.
The domain \(U\) is conformally a disk. Indeed, compose \(X\)
with a conformal map from \(\mathcal T\) to the disk taking
\(w_c\) to zero. Off \(q_c\) the result is a double cover
of the punctured disk. A generator around zero exchanges
the two sheets by simple ramification, so the image of
every closed curve in \(U\setminus\{q_c\}\) has even winding
about zero. The composed map has a single-valued holomorphic
square root there. It extends through its zero of order
two, and is a proper map of degree one onto the disk.

Set
\[
 F(q)=\sum_i a_i\log x_i(q)-\log(1-q),\qquad
 \mathcal D(q)=F(\tau(q))-F(q).
\]
These functions use the upper logarithmic determinations.
The stationary phase differs from \(F\) by a common
constant at the two saddles. The global difference
\(\mathcal D\) satisfies \(\mathcal D\circ\tau=-\mathcal D\).
It vanishes at \(q_c\), and its derivative is
\((q-\tau(q))X'(q)\), using \(F'(q)=-qX'(q)\).
Thus it has a zero of order three at \(q_c\) and no other
critical point. Its real zero set projects to the
equal-modulus set \(E\) of the two saddle values.
Off \(w_c\), a local choice of sheet gives a difference
\(D(w)\) with \(D'(w)=q-\tau(q)\ne0\). Consequently
\(E\setminus\{w_c\}\) is a smooth embedded real analytic
one-dimensional manifold. Near the fold, its arms are
the images of the six local real-zero rays of a cubic.

A *regular active arc* means a connected component of
\(E\setminus\{w_c\}\) where both actual Gaussian coefficients
are nonzero. This property is constant on each component:
the [regular contour localization](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#regular-values-and-nonzero-coefficients)
gives fixed integer contour weights on a neighborhood
of every regular equal-height point, with nonzero local
Gaussian factors. On overlapping neighborhoods each saddle
is strictly dominant on one side of the equal-modulus arc,
so the leading asymptotic there fixes whether its weight
is nonzero in either representation.

## Where an active arc can end

First, an active point \(w_0\) is approached by zeros of
the individual moment polynomials. To see this, use
\(D(w)-D(w_0)\) as a local coordinate and normalize the
two-saddle expansion on its disk by one exponential
and one nonzero amplitude. If \(L(w)\) denotes the
other nonzero amplitude ratio, the resulting expression is
\[
 e^{nD(w)}\bigl(1+O(n^{-1})\bigr)
       +L(w)\bigl(1+O(n^{-1})\bigr)+O(e^{-cn})
\]
on a smaller disk, with \(c>0\). Put
\(z=n(D(w)-D(w_0))\) and choose the bounded solution
\(z_n\) of \(e^{nD(w_0)}e^{z_n}=-L(w_0)\)
whose imaginary part lies in \([-\pi,\pi]\).
On a circle of fixed small radius around \(z_n\),
the two leading terms with \(L(w_0)\) have modulus
bounded away from zero. The errors are \(O(n^{-1})\).
Rouché's theorem gives a zero at distance \(O(n^{-1})\)
from \(w_0\) for each sufficiently large order.
The [positive-real one-saddle estimate](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis)
gives an eventual zero-free disk about every finite
point of the boundary \(\Im w=0\), since the profile is
strictly subcritical. Therefore an active arc cannot
accumulate anywhere on that boundary.

Next, the whole equal-modulus set has bounded real
coordinate. Both nearest row branch points \(p_L<0<1<p_R\)
are finite because the nonzero gaps have zero sum.
As \(\Re w\to+\infty\), the two upper inverses tend
uniformly to \(0\) and \(p_R\) for \(-\pi\le\Im w\le0\),
taking boundary limits at the endpoints of the height
interval. For zero, write \(X=-\log q+B(q)\), with
\(B\) analytic and real on the real axis near zero.
The coordinate \(q e^{-B(q)}=e^{-w}\) has a real
analytic inverse with positive derivative at zero.
It maps small upper-half-plane values into the upper
half-plane and gives the first inverse. Near \(p_R\), put
\(k=\sum_{p_i=p_R}|d_i|\ge1\) and use
\[
 X=-k\log(p_R-q)+B_R(q)-i\pi,\qquad
 (p_R-q)e^{-B_R(q)/k}=e^{-(w+i\pi)/k}.
\]
Here \(B_R\) is analytic and real on the real interval
to the left of \(p_R\). The coordinate in \(p_R-q\)
has a real analytic inverse with positive derivative
at zero, and the exponential on the right is in the
lower half-plane for interior heights. These two inverses
exhaust the degree. The value of \(F\) stays bounded
near zero and its real part tends to negative infinity
near \(p_R\), uniformly in these angles. The analogous
claim as \(\Re w\to-\infty\) follows by reversing
the columns and using \(q\mapsto1-\bar q\); the phases
of the two reversed saddles are conjugate to the original
ones up to the common \(X\). Thus neither end of the
strip contains points of \(E\).

On the lower boundary \(\Im w=-\pi\), the two inverses
are real with \(q_L\in(p_L,0)\), \(q_R\in(1,p_R)\).
On both gaps \(X'>0\), since all \(x_i(q)>0\) and
\(-1/[q(1-q)]>0\), and the real coordinate of \(X\)
runs from negative to positive infinity. Local holomorphic
inversion lifts each branch into the upper half-plane
just above the boundary. There are no more upper
inverses. With the ordering \(D=F(q_R)-F(q_L)\),
the boundary imaginary difference is \(\pi\), from
the term \(-\log(1-q_R)\). The real difference
strictly decreases because \(D'=q_L-q_R<0\).
It tends to opposite infinities at the two ends of this
boundary, by the same endpoint asymptotics. Hence there
is exactly one boundary point \(w_-\) of equal height.
There \(D=i\pi\), and local analytic inversion and
\(\partial_{\Re w}\Re D\ne0\) give exactly one regular
equal-modulus arc entering the strip.

## Classification of active arcs

No regular component of the equal-modulus set can
be a circle. Its lift, or a connected component of
the lift, would be a simple closed curve in \(U\)
on which \(\Re\mathcal D=0\). Since \(U\) is simply
connected, it contains the interior of that curve.
The maximum principle for \(\Re\mathcal D\) would
contradict the nonconstancy of \(\mathcal D\).
The same reasoning excludes an open regular component
with both ends at the fold. Such an arc has an
injective lift on its interior, with both ends tending
to \(q_c\) by properness. The completed lift is a
simple closed curve in \(U\).

An active component has compact closure in the closed
strip by the horizontal bound. It cannot accumulate on
the upper boundary, and on the lower boundary it can
accumulate only at \(w_-\). Any regular interior limit
point of the component belongs to the component: locally
the equal-modulus set is one smooth arc, and activity
is constant there. Therefore the ends of an active
component, if it is not a circle, can accumulate only
at \(w_c\) or \(w_-\). More explicitly, the component
is an embedded open interval, and any cluster set of
one end is an intersection of nested compact connected
closures of its tails. It is nonempty and connected,
so in this two-point set it is a singleton.
Both ends cannot approach \(w_-\) because only one
local arm enters the strip there. This also prevents
two distinct active components from approaching \(w_-\).
With the exclusions of circles and fold-to-fold arcs,
there is at most one active component and it must
connect \(w_c\) to \(w_-\).

On such an arc choose a continuous lift. The imaginary
part of its stationary difference is a strict coordinate
along the arc, since locally \(D'\ne0\) and \(\Re D=0\).
Its endpoint limits are zero at the fold and either
\(\pi\) or \(-\pi\) at the lower boundary, according
to the choice of sheet. Thus \(e^{D(w)}\ne1\) at every
regular active point.

## Consecutive exclusion

The [actual regular expansions](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#regular-values-and-nonzero-coefficients)
give local eventual zero-freeness whenever the two
saddles are unequal in height or only one contributes.
On a regular active arc the consecutive comparison would
require \(e^{D(w)}=1\) at any limit of common roots, which
has just been excluded. The
[simple-fold estimate](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#bounded-scaled-fields-for-equal-total-degrees)
excludes a whole neighborhood of the ramification value.
The positive and
[negative real scale estimates](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#local-exclusion-around-every-negative-real-scale)
exclude limits on the real axis. All possible common roots
lie in the [fixed annulus](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-fixed-annulus-for-possible-common-roots).
Compactness of the annulus and complex conjugation now
exclude any sequence of common roots at unbounded orders.
The integer coefficients then give rational Bézout identities
which extend to every characteristic-zero field. \(\square\)

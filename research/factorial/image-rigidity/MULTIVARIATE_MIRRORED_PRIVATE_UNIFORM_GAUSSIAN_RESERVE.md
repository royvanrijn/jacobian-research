# Uniform angular separation near mirrored private-unit cores

Take the private unit pair \((1,0),(0,1)\) and any finite
multiset of mirrored ordinary exponent pairs. Call this the
base. Equal rows are immaterial. The
[linked phase theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
gives strict angular interlacing for this base with any
Gaussian strength \(W\ge0\). The
[explicit reserve theorem](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md)
already gives numerical per-order bounds on that separation.
Here is an independent compactness proof of a uniform gap in the
relevant low orders, together with its transfer to both asymmetric
ordinary anchors and an arbitrary-strength diffuse tail.

**Theorem.** Let \(J\ge4\) be at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
There is an effective \(g_J>0\) such that for every base
and every \(W\ge0\), the endpoint-normalized reference
polynomials at each \(3\le N=n+1<J\) have cross-order
distance between unit root directions at least \(g_J\).
The bound is independent of the dimension, gaps, row
counts and curvature of the base.

There are consequently effective \(\kappa_0,\eta_0>0\)
depending only on \(J\) with the following property.
Keep any such base and append both:

- arbitrary ordinary unequal anchor rows of total curvature
  at most \(\kappa_0\);
- arbitrary ordinary shared unequal tail rows with relative
  gaps at most \(\eta_0\), with any total curvature and
  signed bias.

The resulting ordinary family has disjoint consecutive root
arguments in every degree. Its moments are coprime over any
characteristic-zero field. Every anchored core in this family
has no positive Gaussian physical contact in the relevant
range \(3\le N<J\), even when the anchors are not mirrored.
There is no restriction on the number of anchors or tail rows.

## Closed phase constraints

For a normalized Gaussian polynomial of degree \(j\), use
its private-unit geometric representation, with the endpoint
law of the positive ratio denoted by \(T\). Its coefficients
are \(\mathbb E_j T^k\), with \(\mathbb E_j T^j=1\).
The Gaussian multiplier is the moment sequence of an
independent lognormal factor with log mean \(-W/2\) and
variance \(W/j\). Tilting this factor by its \(j\)-th power
inverts it in law. Thus the
[exchange identity and resolvent](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md#the-exact-resolvent-and-angle-boxes)
still hold, including at \(W=0\).

For \(0<\theta<\pi\), write
\[
 D_\theta(T)=1-2T\cos\theta+T^2,\quad
 A_j=\mathbb E_j[D_\theta(T)^{-1}],\quad
 B_j=\mathbb E_j[T D_\theta(T)^{-1}],\quad
 \mu_j(\theta)=B_j/A_j .
\]
We have \(0<\mu_j\le1\). The resolvent root equation in
this coordinate is
\[
          \cos(j\theta/2)
                          =\mu_j(\theta)\cos((j/2+1)\theta).    \tag{1}
\]
Reflection gives the roots on the negative half of the
circle. In odd degree the remaining root is \(-1\).
In increasing order the positive root angles of degree \(j\)
are in the closed boxes
\[
           I_{j,k}=\left[\frac{(2k-1)\pi}{j},
                         \frac{2k\pi}{j+1}\right],
                         \qquad 1\le k\le\lfloor j/2\rfloor.    \tag{2}
\]
The geometric case uses the right endpoints. Under the
folded midpoint law, \(A_j\) and \(B_j\) are proportional
with the same factor to the numerator and denominator,
respectively, of \(R_j\) in the linked phase theorem.
Thus \(R_j=1/\mu_j\).

The exact paired inequality (4) of that theorem gives
\[
 \mu_N(\theta)\le H_n(\theta,\mu_n(\theta)),\qquad
 H_n(\theta,u)=
 \begin{cases}
  0,&u=0,\\[2pt]
  \displaystyle\frac{1}{
    1+r(n/N)^{2+\sin(\theta/2)+\operatorname{arcosh}(1+r)}},
              &0<u\le1,\quad r=u^{-1}-1 .
 \end{cases}                                                   \tag{3}
\]
The derivation of its equations (2)–(4) applies also to
\(\theta>\pi/2\). Indeed the common reweighting is positive
for all \(0<\theta<\pi\); \(D=\cosh x-\cos\theta\) is
positive, \(p=(\cosh x-1)/D\) remains increasing in \(x\),
and the bound on \(s(1-p)\) uses only
\(1-\cos\theta=2\sin^2(\theta/2)>0\). The logarithmic
derivative bound and its integration are unchanged.
For the empty geometric case (3) gives \(\mu_N=\mu_n=1\).

The extension by zero in (3) is continuous and effective.
For \(u\le1/2\), \(r\ge1\),
\(\operatorname{arcosh}(1+r)\le\log(4r)\),
\(\log(N/n)\le1/2\), and \(n/N\ge2/3\). Hence the term
added to one in its denominator is at least
\((4/27)\sqrt r\). In particular \(H_n(\theta,u)\le10\sqrt u\)
in this range, uniformly in \(\theta\).

There is also a uniform Lipschitz bound
\[
          |\mu_j(\theta)-\mu_j(\varphi)|
                         \le 9J^4|\theta-\varphi|,
                \qquad \theta,\varphi\in[\pi/J,\pi].           \tag{4}
\]
Indeed \(T\le2\) has probability at least \(3/4\) by
\(\mathbb E_j T^j=1\) and \(j\ge2\), so \(A_j\ge1/12\).
Write \(\delta=\sin^2(\pi/J)\).
For \(T\le2\), \(D_\theta(T)\ge\delta\) on the indicated
interval, and the integrands differentiated in \(\theta\)
have magnitudes at most \(4/\delta^2\) and \(8/\delta^2\)
for \(A_j\) and \(B_j\), respectively.
For \(T\ge2\), \(D_\theta(T)\ge T^2/4\), so those
bounds still hold. Differentiation is justified by these
uniform bounds. Then \(|\mu_j'|\le144/\delta^2\) follows
from \(B_j/A_j\le1\). The inequality
\(\sin(\pi/J)\ge2/J\) gives (4).

## The compact comparison separates uniformly

Fix \(n,N<J\) and two positive root labels. Consider the
compact set of quadruples \((\theta,\varphi,u,v)\) with
angles in their boxes (2), \(0\le u,v\le1\), equations
(1) at the two corresponding angles and parameters, and
\[
                   v\le H_n(\theta,u)+9J^4|\theta-\varphi|.   \tag{5}
\]
It contains the data of every finite mirrored base and
every Gaussian strength, by (3)–(4). We show this
compact set misses \(\theta=\varphi\).

Different labels cannot meet because the boxes for a lower
label \(k\) and the upper label \(k+1\) are separated by
at least \(\pi/N\), and the upper label \(k-1\) has
upper endpoint below the lower label \(k\)'s left endpoint.
Take equal labels \(k\). If a common angle exceeds
\(\pi/2\), the upper bound of the lower box gives
\(4k>n+1\). Hence the left endpoint of the lower box
is at least the right endpoint of the upper box.
Equality would require both endpoints, so \(u=0,v=1\)
by (1). This contradicts (5).

For \(\theta=\varphi\le\pi/2\), \(u=0\) implies \(v=0\)
by (5). Both phases would then have to be left endpoints,
and those differ because \(N>n\). If \(v=0\) alone,
the upper left endpoint is already below the lower box.
Thus \(u,v>0\), and \(n\ge3\) since for \(n=2\) the
only possible angle in this range is the lower left
endpoint, which requires \(u=0\). If \(r=u^{-1}-1\),
the tangent addition formula applied to (1) gives
\[
       \frac1v=\cos\theta+
                      \frac{(1+\cos\theta)r}{r+2}.
\]
For \(r=0\) this contradicts \(v\le1\).
For \(r>0\), subtract one from the right hand side.
The strict scalar inequality (5) of the
[linked phase theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md#a-scalar-phase-comparison)
bounds that expression strictly below
\(r(n/N)^{2+\sin(\theta/2)+\operatorname{arcosh}(1+r)}\).
This contradicts (5) here, since (5) at equal angles
requires \(1/v-1\) at least that quantity.

The compact sets therefore have a positive minimum
\(|\theta-\varphi|\) whenever nonempty. A rational
positive lower bound is computable. For example subdivide
the finite boxes for \((\theta,\varphi,u,v)\) and use interval
enclosures of the continuous equations and (5) to exclude
\(|\theta-\varphi|\le\epsilon\) while decreasing \(\epsilon\).
All these functions have computable uniform enclosures
on the box, including at \(u=0\) by the bound after (3).
Since the compact set at equality is empty, sufficiently
small fixed \(\epsilon>0\) and sufficiently fine subdivisions
certify the exclusion. Equivalently the maximum of the
constraint violations has a positive minimum on the closed
strip being excluded.

Each positive root angle lies between \(\pi/J\) and
\(\pi-\pi/J\) on these orders. Thus opposite-half comparisons
and comparisons with the odd structural root \(-1\)
also have a positive separated chord. Taking the minimum
over labels and orders gives \(g_J\) for the unit-circle
distances, using \(2\sin(|\theta-\varphi|/2)\) on the
same positive half. This proves the uniform angular reserve.

## Anchors, arbitrary-strength tails and an exact marked family

For the exact ordinary profile, the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
already separates the arguments if total curvature is at
least \(27\), so assume the diffuse tail has \(\tau\le27\).
In each order compare the actual endpoint-normalized
polynomial to the reference on the mirrored base at \(W=\tau\).
The anchor contribution has log coefficient modulus at
most \(j\kappa_0/4\) if \(\kappa_0\le1/2\). In fact such
a bound excludes private anchor rows, and on each anchor
factorial chord the relative gap is at most \(\kappa_0\).
The [quadratic chord bound](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
then gives at most \(j/4\) times the row curvature.
For the diffuse tail, the
[uniform ordinary comparison](MULTIVARIATE_UNIFORM_TIGHT_CORE_DIFFUSE_STABILITY.md#transfer-to-the-appended-ordinary-rows)
has log error at most \(j\tau\eta_0/4\) if \(\eta_0\le1/2\).
These errors add regardless of the signed bias or degree gap.

The earlier explicit reserve also gives immediate choices
\(\kappa_0=1/(10^7J^6)\) and
\(\eta_0=1/(27\cdot 10^7J^6)\) for ordinary extensions:
the log coefficient error is then at most half its
\(1/(10^7 J^5)\) tolerance.


The normalized coefficients of both polynomials are at
most one, because the private unit pair cancels the
binomial coefficient and the normalized factorial and
Gaussian chords each have nonpositive logarithm.
The endpoints are one. With
\(\varepsilon=J(\kappa_0+27\eta_0)/4\), the coefficient
differences are at most \(\varepsilon\).
All roots of both polynomials have modulus between
\(1/2\) and \(2\) by the leading term bound and reversal.
If \(z\) is an actual order-\(j\) root, the reference
monic polynomial evaluated there has modulus at most
\(3^j\varepsilon\). Its factorization gives a reference
root within \(3\varepsilon^{1/j}\). The unit directions
are then within \(12\varepsilon^{1/j}\). Choosing the
effective constants so that \(\varepsilon\le1\) and
\(24\varepsilon^{1/(J-1)}<g_J\) preserves separation
at every \(3\le N<J\). The pair starting at one follows
from even real positivity. The absolute angular cutoff
covers the remaining orders. For a Gaussian reference on
the anchored core the only comparison error is the anchor
error. It is independent of \(W\) and gives the no-contact
claim on the entire positive strength interval below \(J\).

For explicit nonmirrored marked cores, keep any finite
base of the theorem. For each \(3\le N<J\) choose
a distinct prime \(p_N>J^2\) also exceeding every
factorial argument of the base at the low orders.
By the Chinese remainder theorem there are arbitrarily
large \(A\) with
\[
                          NA+1\equiv p_N\pmod{p_N^2}
                              \qquad (3\le N<J).
\]
Use the single increasing anchor \((A,A+1)\).
At the lower order \(n=N-1\) none of its literals
\(nA+h\) can be divisible by \(p_N\), since a common
divisor of that literal and \(NA+1\) divides
\(Nh-n\), a nonzero integer of magnitude below \(J^2\).
At \(N\), only the first literal has positive valuation,
and it has depth one. These are the
[first-block marked cores](MULTIVARIATE_FIRST_BLOCK_MARKED_GAUSSIAN_CONTACTS.md)
at all the relevant pairs. The anchor is unreflected
and its curvature tends to zero. The uniform reserve
excludes every positive reference clock in this subclass,
including the algebraic possibilities left open by the
local marked-prime tests alone.

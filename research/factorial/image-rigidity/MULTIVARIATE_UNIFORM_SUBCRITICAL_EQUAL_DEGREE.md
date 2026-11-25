# Uniform eventual coprimality on subcritical equal-degree bands

Remove equal rows from an ordinary two-monomial support and
assume both entries of each remaining row are positive.
Use \(W,\rho\) from the
[bounded supercritical comparison](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md).
Here the total degrees agree, so \(\Delta=\sum_i(b_i-a_i)=0\).

**Theorem.** Given \(0<W_-\le W_+<\infty\) and
\(0<\kappa<1\), there is one terminal threshold for all
fully shared ordinary equal-degree profiles satisfying
\[
 W_-\le W\le W_+,\qquad \max_{0<q<1}\rho(q)\le1-\kappa .
\]
At every larger terminal the consecutive moments are
coprime, even when the dimension, the exponent columns and
the number of retained rows vary with the terminal. The
threshold is unevaluated, and the result holds over every
characteristic-zero field.

## The limiting double covering

Suppose a counterexample sequence has \(n\to\infty\) and
pass to a retained and diffuse limit. The normalized
parameter \(w=-\operatorname{Log}t-C_0\) is the same at
both orders, and the common roots have bounded \(\Re w\)
by the annulus. Use the notation \(G_*,X_*\), and the
limiting phases from the supercritical comparison.
The limit satisfies \(\max q(1-q)G_*''(q)\le1-\kappa\),
by convergence on \([0,1]\). It has positive mass because
of \(W_-\). Its proper degree over
\(\mathcal T=\{-\pi<\Im w<0\}\) is two, by the
[stationary compactness theorem](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md).
There is at most one branch point there. In fact any
interior critical point of the limit persists as a
critical point of the nearby ordinary centered map by
Rouché's theorem for the derivative, with the image
still in \(\mathcal T\). Each ordinary subcritical map
has exactly one simple upper branch point over that
strip. A critical point of the limit must be simple,
also by its degree-two count.

There must be a branch point. Otherwise both branches
would be conformal over the strip. They can be labeled
near its lower boundary by the two real-gap inverses.
On each gap \(X_*'>0\) and the real coordinate of
\(X_*\) runs onto the real line at height \(-\pi\).
This is the real-gap argument of the supercritical
comparison and does not need its curvature hypothesis.
For large positive \(\Re w\), the local inverse of
the singular term \(-\log q\) near zero continues from
the left real gap to \(0<q<1\) as the height varies
from \(-\pi\) to zero. More explicitly, the inverse
of \(q e^{-B(q)}=e^{-w}\), where \(B\) is real analytic
at zero, maps a small upper half-plane into the
upper half-plane. Similarly the local inverse near
one continues from the right real gap to \(0<q<1\)
at large negative \(\Re w\), by the term \(\log(1-q)\).
On \((0,1)\), \(X_*\) decreases strictly from
\(+\infty\) to \(-\infty\) with nonzero derivative.
The unique real inverse there lifts into the upper
half-plane just below height zero. Continuation
along this real boundary branch therefore connects
the two purported global sheet labels, a contradiction.
Thus there is exactly one simple ramification point
and the covering is connected.

## Actual coefficients and the active arc

The same finite Rice homotopies as for the
[supercritical band](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#finite-contour-deformation-for-the-limiting-profile)
apply here at every regular parameter. The real cut
exclusions are independent of criticality. The
far-flow speed proof is also unchanged except
possibly when \(\tau=0\) and all retained gaps
have one sign with total absolute gap at most two.
In that exception the retained list is finite and
\(\Re G_*'(q)\) grows like a nonzero signed multiple
of \(\log|q|\) uniformly for large \(|q|\) in either
half-plane. Thus the speed still has a positive
lower bound outside a large disk. The finite
energy and finite homotopy arguments apply.

The homotopy first isolates the larger saddle if the
heights differ. If its coefficient is zero, flow the
remaining strict sublevel chain into a rectangle about
the smaller saddle or its own strict sublevel. At least
one coefficient is nonzero, by the
[fixed bounded gap winding obstruction](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#regular-values-and-nonzero-coefficients).
If at least one row is retained there is a bounded real
gap and that obstruction applies verbatim. If no row is
retained, the limiting mass is purely diffuse and the
real phase on the entire negative real ray is bounded
below, because its leading term at infinity is
\(\tau q^2/2\) with \(\tau>0\). A loop in a sublevel below
that ray's minimum has the same winding at zero and at
infinity, also a contradiction.

Thus the finite limiting homotopies transport to all late
actual profiles as in the supercritical proof, with nonzero
integer crossing weights and strict loss on all bounded
exterior chains. At each regular equal-height point they
give uniform local expansions of the two actual terms.
Their weights are constant along the regular arc: this
follows by comparing the leading expansions on either
side of the arc, where each saddle in turn has strictly
larger limiting height. A regular arc is active if both
weights are nonzero. At any point \(w_0\) of an active
arc, actual individual zeros tend to \(w_0\) along the
approximating sequence. Indeed let \(D_*\) be the limiting
stationary difference, with \(D_*'(w_0)\ne0\).
In a small disk first move the finite-profile center
by \(o(1)\) to solve \(D_n^*(w)=D_*(w_0)\).
Then scale \(D_n^*(w)-D_*(w_0)\) by \(n\) and apply
Rouché to the sum of the two exponentials as in the
[fixed subcritical active-arc proof](MULTIVARIATE_SHARED_SUBCRITICAL_EQUAL_DEGREE.md#where-an-active-arc-can-end).
The nonzero analytic amplitude ratios are uniformly
bounded above and below there, so the small circles in
the scaled field can be chosen uniformly.

We can use the topology of the
[fixed subcritical covering proof](MULTIVARIATE_SHARED_SUBCRITICAL_EQUAL_DEGREE.md#classification-of-active-arcs)
on the limiting centered map. Its connected double cover
with one simple ramification is a disk, by the holomorphic
square root of a disk-valued target coordinate centered at
the branch value. On it the stationary difference on the two
sheets is a global holomorphic function, with a cubic zero at
the ramification and derivative nonzero elsewhere. Therefore
no regular equal-modulus component is closed, and no such
arc has both ends at the fold: in either case its lifted
closed curve would violate the maximum principle on the
covering disk.

An active arc cannot escape horizontally, because its
points are limits of actual individual roots and those
have a uniform centered annulus. Nor can it reach the
upper boundary, where the
[uniform positive subcritical saddle](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md)
and the endpoint normalization give a zero-free
neighborhood at every fixed real scale. On the lower
boundary the inverses lie in the two real gaps. There
\(D'\) is \(q_L-q_R<-1\), so \(\Re D\) is strictly
decreasing from positive to negative infinity and has
one zero, with \(D=i\pi\). The two vertical contour
expansions and the local inverse argument of the
supercritical comparison apply unchanged. In particular
only one local equal-modulus arm enters the strip there.
At every regular interior accumulation point the
equal-modulus set is locally a single smooth arc, so
an active component cannot stop there.

Within the bounded closed strip, an end of an active
component consequently has its cluster set in the
two-point set of the fold and the lower boundary point.
The cluster set of an end is connected, by the nested
closures of connected tails, hence is a singleton.
The no-loop argument and the unique boundary arm
show that the only possible active component joins
the fold to that lower endpoint. Along a lift of this
component the imaginary difference is a strict
coordinate, with limits zero and \(\pi\) or \(-\pi\)
depending on orientation. Thus the multiplier
\(\exp D_*\) is never one at any regular active point.

For common roots tending to a regular active point,
compare the two actual cancellation equations. Because
\(\Delta=0\), both use the same centered parameter, so
their quotient tends to \(\exp D_*(w_0)\) on the
left and to one on the right. This is impossible.
Regular points with one active saddle, or two unequal
heights, are already excluded by the actual
Gaussian expansions. At the lower real equal point
the multiplier is \(-1\), which likewise separates
consecutive orders. This leaves the interior fold.

## Uniformity at the interior fold

The fold is handled by transport of the full
[physical cubic comparison](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md).
This transport requires keeping the actual contour
class, a common coordinate disk, and the adjacent
amplitude estimates. At the limiting fold, the
finite cubic trapping homotopy applies in the present
cut domain by the boundary and escape estimates above.
It gives a nonzero local integer Airy class: a zero
class would allow descent below the fixed real gap
barrier, or below the negative real ray in the purely
diffuse case. This finite contour stays in the actual
domains for every sufficiently late profile.
Its exterior chain has a common strict phase gap.

The actual simple critical points and their values
converge to the limiting ones. Their analytic cubic
coordinates exist on one common smaller disk, with
uniform bounds. For example use the same bounded
inverse of the linearized cubic coordinate from the
[local construction](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#localizing-the-exact-contour-at-a-fold),
after shifting to each nearby actual fold. The
implicit-function estimates are uniform by analytic
convergence in the local chart and the nonzero limiting
cubic coefficient and field derivative.
The endpoints can be continued to the three nearby
actual decay sectors without leaving the strict
sublevel, so the integer Airy class remains fixed.
The bounded analytic coefficients in the uniform
Stirling expansion, the inverse quadratic-form
estimates on a compact local chart, and the common
pole-free domain provide the amplitude estimates of
the physical cubic comparison uniformly for all late
profiles. In particular the local amplitude at the
fold stays bounded away from zero. The complex-order
Stirling remainder on that local chart also has the
uniform derivative bound used in its bounded-field
proof, by Cauchy's estimate on a common sector around
the positive order ray.

The large-scaled-field determinant estimate of that
proof, in the equal-degree case, now applies with a
single sufficiently large cutoff on \(|n^{2/3}\lambda|\),
where \(\lambda\) is the actual profile-dependent Airy
coordinate. The bounded-scaled-field comparison uses
the same nonzero Airy function for the entire sequence,
with simple zeros away from zero. Its local
differentiable interpolation keeps the profile fixed
while comparing \(n\) with \(n+1\), and the uniform
bounds just verified make the adjacent branch
displacement estimate uniform. Because the centered
parameter agrees at the two orders, its actual
\(\lambda\) agrees as well at a common root. This is
impossible throughout a fixed disk about the
limiting fold. Thus all possible limiting roots are
excluded and the sequence of counterexamples cannot
exist. Rational Bézout gives the field assertion.

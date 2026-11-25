# Effective entry into the shared quartic neighborhood

Use the fully shared ordinary profiles, curvature \(W,\rho\),
centered fields \(w_j\), and varying real center \(w_*\) of the
[moving quartic window](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Put \(\delta=\max_{[0,1]}\rho-1\).

**Theorem.** For every prescribed positive rational curvature band
\(W_-\le W\le W_+\) and rational \(b>0\), one can compute a
positive rational \(\delta_b\) and integer \(J_b\) such that,
for \(N\ge J_b\) and \(|\delta|\le\delta_b\), a common root
of the orders \(n=N-1,N\) must satisfy
\[
                         |w_j-w_*|<b\quad (j=n,N).
\]
This is uniform in the ordinary exponent rows and the degree
difference. The constants are unevaluated.

## Limit witnesses at criticality

Use the compact centered field rectangle from the
[effective shared real-boundary proof](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md).
The retained and diffuse profile nets enclose the centered
phases and their derivatives on fixed contour neighborhoods.
On \([0,1]\) the real curvature and its derivatives have
uniform computable bounds, as in the
[effective bounded quartic proof](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#effective-center-and-localization-on-the-real-lattice).
That proof computes the real quartic center by locating
the unique curvature maximum, then the nearby simple zero
of the third derivative of the real phase. It gives a
computable initial bound on \(|\delta|\) on which \(w_*\)
is a continuous computably enclosed function of the profile.
At an exactly critical profile that center is \(w_c=X(q_c)\),
where \(q_c\) is the curvature maximizer.

Consider an exactly critical limiting profile and a pair of
fields at equal heights in the closed strip, with at least
one at distance at least \(b/2\) from \(w_c\).
The [critical limiting graph](MULTIVARIATE_CRITICAL_REGULAR_COMPLEMENT.md#the-critical-limiting-graph)
has one equal-modulus point at each interior height, with
two simple branches, distinct real parts, and
\(e^{P_R-P_L}\ne1\). The same statements hold at the
lower endpoint, with \(e^{P_R-P_L}=-1\). The two actual
terms are present on the active graph, while each strict
dominance region has a nonzero leading term by the
[critical contour argument](MULTIVARIATE_CRITICAL_REGULAR_COMPLEMENT.md#actual-contours-and-the-adjacent-comparison).
These contours give the rational witnesses of the
[effective regular cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md)
using finite polygons with strict sublevel loss, simple Rouché
saddle disks, and a checked winding at the coefficient poles.
The curvature condition at criticality is used there to
guarantee existence of the finite contour. The certificate
itself uses strict local inequalities and works for nearby
profiles on either side.

For an interior pair, if both fields lie at active balances
they coincide by the height coordinate. Use a common
convex chart with positive lower bounds on the real drift
slope and on \(|e^{P_R-P_L}-1|\), where the effective
regular comparison applies. Otherwise at least one field
has a dominant-term witness. At the lower boundary the
two real-gap witnesses and the same common-chart
comparison apply as in the real-boundary proof.

At the upper boundary at least one of the two fields differs
from \(w_c\) by \(b/2\) or more. Its real entropy phase has
only one maximum, which is simple. Indeed the second
derivative is nonpositive and vanishes only at \(q_c\),
whose target is \(w_c\). A one-maximum lattice witness is
therefore available on a small complex neighborhood of that
field. No asymptotic comparison at the cusp itself is used.

## Effective transfer to small \(|\delta|\)

Let \(\mathcal C\) denote the compact critical family of
all these field pairs and profiles. Every point of
\(\mathcal C\) has an open witness of one of the kinds above,
with a computable local terminal threshold. The exact
critical compact set is used only to prove termination.

Dovetail over rational curvature radii \(\delta_b\) below
the initial center bound, profile and field meshes,
interval precisions, and rational witnesses. At each
stage require a witness on every cell that cannot be
strictly excluded from
\[
              |\delta|\le\delta_b,\qquad
              \max(|w_n-w_*|,|w_N-w_*|)\ge b .
\]
Each proposed exclusion is a finite interval bound for
\(\max\rho\) or \(w_*\). A sufficiently small curvature
radius makes this closed admissible family lie in the
open union of a finite set of witnesses for \(\mathcal C\).
Otherwise a sequence of admissible profiles and field
pairs outside that open union with \(\delta_b\downarrow0\)
would have a limit at criticality with distance at least
\(b\), contradicting the cover of \(\mathcal C\).
The slack between \(b/2\) and \(b\) also allows the
outward errors of the \(w_*\) enclosures. Refining the cells
then verifies a finite set of witnesses and exclusions.
Taking the largest local terminal from that stage gives
\(J_b\) and proves the assertion.

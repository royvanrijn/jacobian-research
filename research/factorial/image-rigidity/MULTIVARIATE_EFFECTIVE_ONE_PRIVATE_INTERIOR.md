# A finite effective contour cover for the one-private interior

Consider ordinary profiles with exactly one private row \((0,1)\).
Use the fields \(\xi_j\), degree difference \(\Delta\), and shared
endpoint curvature \(\sigma\) of the
[one-private regular comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md).
Fix a rational \(0<\epsilon<\pi/2\).

**Theorem.** There is an algorithm computing an integer \(N_0(\epsilon)\)
such that no consecutive common root occurs for \(N\ge N_0(\epsilon)\),
\(0\le\sigma\le5/4\), and
\[
                  \epsilon\le\Im\xi_j\le\pi-\epsilon .
\]
The result is uniform in the profile and every integer degree
difference. The algorithm searches finite contour certificates.
Their termination uses the exact qualitative activation theorem.
No numerical cutoff or complexity bound for this search is asserted.

## Compact phase family and root-bearing fields

Write \(u_i=d_i/a_i\), \(s_i=d_i^2/a_i=|d_i||u_i|\), and
\[
 B(q)=\sum_i d_i\log(1+u_i q),\quad X(q)=\log(1-q)+B(q),\quad
 \Phi_+(q,\xi)=\xi q+\int_0^q X(z)\,dz .
\]
In retained and diffuse limits add \(\tau q\) to \(B\); then
\(\sigma=\sum s_i+\tau\). The lower phase uses lower conjugate
logarithms and the additional term \(-2\pi i q\).
The real field bound \(|\Re\xi_j|\le13\) for actual roots is
proved in the
[positive-boundary entry](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_POSITIVE_ENTRY.md#a-single-effective-saddle-at-the-positive-real-fields),
so use a closed field rectangle \(K\) with this real cap and
the stated strip margins.

Here is an explicit compact net for the phase family. At a rational
cutoff \(\rho>0\), the rows \(|u_i|\ge\rho\) have bounded integer
gaps \(|d_i|\le(5/4)/\rho\), bounded exponents
\(a_i\le(5/4)/\rho^2\), and each contributes at least \(\rho\) mass.
Their finite lists can therefore be enumerated. Mesh the possible
remaining mass \(\tau\) for each list in
\([0,5/4-\sum_{\rm retained}s_i]\). On a fixed disk where
\(\rho |q|\le1/2\), the omitted-factor expansion gives
\[
 \left|\sum_{|u_i|<\rho}d_i\log(1+u_i q)-\tau q\right|
                         \le (5/4)\rho |q|^2.              \tag{1}
\]
Cauchy bounds on a larger disk give all corresponding
fixed-derivative errors. The same bound holds with an additional
true diffuse component by adding it to the omitted mass.
Integrate (1) along straight segments to control the phase.
The large factors and their cuts are retained exactly.
A bounded region sees no discarded pole once \(\rho\) is small enough.
Thus these nets are effective on every fixed compact contour
neighborhood that avoids the retained cuts. This is the
boundary-centered form of the
[retained and diffuse compactification](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md#retained-and-diffuse-limits).

The one-sheet theorem gives exactly one simple upper saddle
\(q(\xi)\) with \(X(q(\xi))=-\xi\). Its stationary value
\(D(\xi)=\Phi_+(q(\xi),\xi)\) has \(D'(\xi)=q(\xi)\).
The effective fiber rectangle and inverse bounds follow also from the
[explicit fiber localization](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md)
and [inverse estimates](MULTIVARIATE_EFFECTIVE_STATIONARY_INVERTIBILITY.md).
For this application include the private row in \(W\), so
\(2\le W\le2+2\sigma\le9/2\). The change from the midpoint map
to \(X(q)\) is the real constant
\(\log2+\sum_i d_i\log(a_i/m_i)\), with a computable bound
from the shared ratio bounds. Enlarge the real target cap by
this bound when applying the fiber estimates.

For later use \(D\) is never zero inside the strip.
Its derivative has positive imaginary part, so \(D\) is injective
by integrating the derivative on the straight segment between two
fields. It approaches zero at an accessible real boundary field:
for \(\sigma<1\), the inverse at \(\xi=0,q=0\) is simple,
and for \(\sigma=1\) it has the local form
\(X(q)=-c_2q^2+O(q^3)\), \(c_2>0\).
The latter follows exactly as in the moving cubic phase,
with the diffuse term affecting only the linear coefficient.
It gives an upper inverse tending to zero as
\(\xi=i t\) tends to zero through \(t>0\).
If \(\sigma>1\), the strictly decreasing real function
\(R(q)=(1-q)B'(q)\) starts above one and ends at zero.
Thus real \(X(q)\) has a single interior maximum.
For real \(\xi\) just above minus that maximum, the right
stationary value is negative, since its limiting integral
is that of \(X(q)-\max X\) from zero to the stationary point.
At \(\xi=0\) the analogous value is positive.
The stationary value increases with the real field, so
it vanishes at a real coexistence field with simple saddle
in \((0,1)\). The upper inverse continues into that point.
In each case an interior zero of \(D\) would give a local
inverse near zero, contradicting injectivity for a sequence
of images approaching zero from the boundary.

For \(\sigma<1\), the
[finite activation and sign comparison](MULTIVARIATE_UNIFORM_SUBCRITICAL_PRIVATE_ACTIVATION.md#connected-sign-regions-and-the-two-anchors)
show that a field can carry the limiting balance of saddle and
endpoint only when \(\Re D=0\) and \(\Im D>0\). On that sign region
\(\Re q>0\). On a horizontal line the sign region is an interval,
because \(\partial_x\Im D=\Im q>0\). There is therefore at most
one possible balance point at each height, since
\(\partial_x\Re D=\Re q>0\) on the interval.

For \(\sigma\ge1\) the upper saddle lies on the right throughout
the strip. Indeed on the imaginary \(q\) axis,
\[
 \Im X(iy)=-\arctan y+
             \sum |d_i|\arctan(|u_i|y)+\tau y\ \ge0,
\]
by subadditivity of the arctangent and \(\sigma\ge1\).
Thus a strictly negative target height cannot occur there.
The branch is on the right near the negative-real gap, or near
the inverse tending to \(1\) at large positive real field, and
cannot cross the axis. In this case \(\Re D\) is strictly
increasing on each horizontal line throughout.
For both ranges of \(\sigma\), the conclusions needed below are:
a limiting balance is a regular two-contribution point with
\(\Re q>0\), and two such fields at the same height coincide.
The formulas hold for retained and diffuse profiles by the same
arguments, with the sign classification for \(\sigma<1\)
given in the cited activation proof.

## What a finite contour certificate records

Fix one limiting profile and interior field. The domain
\(\Omega\) formed by the upper and lower half-planes and the
real right gap \((1,p_R)\) is simply connected: its complement
in the Riemann sphere is the union of the two removed real
rays joined at infinity, with one ray omitted if \(p_R=\infty\).
On this domain the two phase charts have the same real height
and derivative across the gap. Their imaginary constants differ
by \(2\pi i\), immaterial in the integer-order kernel.
The exact Rice integrand for an ordinary profile has no poles in
the corresponding domain. Indeed all coefficient poles are
on \([0,1]\), and the right shared Gamma poles are strictly
beyond the shared gap endpoint. The private unit row only adds
poles on the removed left ray. A compact arc in the limiting
domain is contained in the finite-profile domains on a sufficiently
close profile cell by (1).

Keep the central pinched left segment at \(\Re q=-c/n\),
\(c=1/2\), with \(|\Im q|\le Y_n=C\log n/n\).
The complementary contour is an open arc from its lower endpoint
to its upper endpoint, through \(\Omega\); it may be replaced
by any other such arc by analyticity. One does not have to
certify a homotopy numerically.
At the limiting two copies of zero the endpoint descent
directions are respectively those for \(\xi-2\pi i\) and \(\xi\).
They point into the appropriate half-planes and have strictly
negative linear real-phase increments.
The [finite descent argument](MULTIVARIATE_UNIFORM_SUBCRITICAL_PRIVATE_ACTIVATION.md#a-finite-deformation-with-two-endpoint-traces)
in that setting is described in its sections on the phases and
the two endpoint traces.
For clarity its geometric part applies here also for \(\sigma\ge1\).
The lower vertical velocity is at most \(-\epsilon\).
Any retained negative row gives an upper velocity bounded
positively away from zero at infinity, since its integer gap
contributes at least \(\pi-\arg q-o(1)\) and
\(\arg(1-q)=\arg q-\pi+o(1)\).
With only positive retained rows, real descent bounds the real
coordinate above and a positive retained row then compensates
\(\arg(1-q)\) at infinity.
With a diffuse term alone it controls large upper heights;
at bounded height and far left the entropy argument tends to zero.
In the flat case the real logarithmic derivative grows in modulus.
These estimates and the real-cut exclusions cited in that proof
prevent escape above a fixed stopping level.
They do not depend on \(\sigma<1\).
The one-sheet count used in the same proof holds for the
whole one-private family.

If \(\Im D\ne0\), the upper endpoint trace cannot meet the saddle,
by constancy of the imaginary phase, and the lower trace cannot
return upward. If \(\Im D=0,\Re D>0\), neither endpoint can
reach the higher saddle by descent.
In either case the finite descent reaches a strict sublevel
below both the boundary and the saddle, outside its local
quadratic rectangle, and the endpoint traces themselves enter
that sublevel. The only remaining case \(\Im D=0,\Re D<0\)
allows stopping strictly between \(\Re D\) and zero, with
no saddle retained. The case \(D=0\) was excluded above.
The resulting finite contour crossing construction yields
the following witnesses.

- Two short endpoint rays in the respective half-planes, with
  strict linear loss from zero, and bounded endpoint trace pieces
  following each ray, whose real phase is strictly below zero.
- Finitely many remaining paths strictly below both saddle and
  endpoint heights, or strictly below zero in the no-saddle
  stopping case.
- A finite oriented number \(k\) of short saddle segments, when
  needed, in a direction of strict quadratic descent. Both ends
  join the strict sublevel through short connectors.

Here the endpoint trace pieces can be kept as part of the
endpoint error, so they need only be below zero.
In the finite construction, each initial descent arc can be
replaced by a short rational ray and a small connector with the
same strict linear loss. Indeed its initial tangent has a
strictly negative directional real derivative, an open
half-plane condition. Away from the saddle and the two copies
of zero all paths are compact in open strict regions.
They can be approximated by finite rational polygons without
losing the inequalities. Replace the local descending saddle
arc by a short straight segment with rational direction through
its center. Its endpoints and small joins retain the strict
sublevel condition. A witness can therefore use rational
vertices except for the center of a saddle segment, which is
enclosed analytically. All saddle traversals use the same segment, with the induced
sign of their orientation. Its concatenation is an oriented
path with the two prescribed endpoints, permitting a finite
number of loops.

Such a witness and all its strict inequalities can be checked
on a profile and field cell with finite interval tests.
Use (1) for value and derivative errors; for a polygon segment
use a finite subdivision with a derivative bound to verify
strict loss along the whole segment.
Membership in the contour domain is checked by separation
from the relevant retained cuts and the distant omitted poles.
For the saddle, search for a rational small disk on which
Rouché against a linear function certifies its location.
For example take center \(q_0\) and radius \(d>0\), and check
\[
 |X(q_0)+\xi|+\tfrac12 M_2 d^2<|X'(q_0)|d
\]
with \(M_2\) an upper bound for \(|X''|\) on the disk.
The exact one simple saddle ensures termination of a refining
search at a center profile. The strict inequality persists
under the interval errors of a small profile and field cell.
It also allows uniform bounds for the moving saddle segment
and its joins. A quadratic descent check on that segment
can use the second derivative at the center and a third
derivative majorant. Thus finding a positive margin requires
only finite interval evaluations and strict inequalities.
The qualitative finite descent is used for existence of a
witness, not as a numerical test or assumed finite flow time.

## Certified expansion errors

The exact central left integral equals
\[
                   L(\xi)+O(1/n),\qquad L(\xi)=(1-e^\xi)^{-1},
                                                               \tag{2}
\]
with uniform computable constants on \(K\).
Indeed on \(|\Im z|\le C\log n,\ z=nq\), the normalized
Gamma ratios differ from one by
\(O((1+|\Im z|)^2/n)\), as in the
[negative boundary proof](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md#the-two-exact-vertical-integrals).
The remaining model integrand
\(\Gamma(-z)\Gamma(1+z)(-e^\xi)^z\) has exponential decay
both ways at rate at least \(\epsilon\) on the left line,
so the exact truncated integral differs from its Mellin
integral by \(O(1/n)\), taking \(C\) sufficiently large.
This uses no decay assumption on the full actual left contour.

To attach the certified endpoint rays to the true left segment,
truncate each ray at imaginary coordinate \(\pm Y_n\)
and join horizontally to the corresponding segment endpoint.
At the true endpoint \(-c/n\pm iY_n\) the linear real loss is
at least a constant times \(Y_n\), uniformly on \(K\).
At the ray endpoint the same holds by the certificate.
The intervening linear real phase is affine and its nonlinear
error is \(O(Y_n^2)\), so the connectors retain strict loss.
Near zero, reflect the pinched Gamma pair
\(\Gamma(-z)\Gamma(1+z)=-\pi/\sin\pi z\).
For \(|\Im q|\ge Y_n\) its upper and lower exponential
bounds give the two phases already used, with a prefactor
bounded by \(Cn\) on the small disk.
The other Gamma factors there have uniform positive real
arguments after scaling, so their errors are effective.
The endpoint rays beyond height \(Y_n\) thus contribute
\(O(n^{-2})\) as an absolute normalized error for sufficiently
large \(C\). The remainder of either endpoint trace has
fixed strict loss below zero and contributes exponentially
little.

On all other compact contour pieces the exact normalized
Rice asymptotic has prefactor \(n\) and relative error
\(O(1/n)\), as in the
[regular right contour expansion](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-actual-right-saddle-on-compact-field-disks)
with private row count one.
These error constants are computable on each certified chart:
the Gamma arguments divided by \(n\) and their row scales
are separated from zero and the negative ray, or are in
the positive-real right gap.
The explicit Stirling remainder is bounded by the first
omitted term times a power of \(\sec(\arg z/2)\).
The summed row error uses \(\sum_i1/a_i\le\sigma\).
The amplitude has uniform upper and positive lower bounds by
the same finite factor separations for retained rows and
\(\sum_i|u_i|\le\sigma\) for the remaining logarithms.
This gives rational bounds from the interval certificates.
[DLMF Stirling remainder bounds](https://dlmf.nist.gov/5.11.ii)

The strict-sublevel pieces are therefore exponentially
small relative to the relevant terms.
For a retained saddle, symmetric integration in its
straight local direction gives relative \(O(1/n)\)
after subtracting odd amplitude and cubic Taylor terms,
as in the [negative boundary proof](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md#the-right-line-and-its-uniform-error).
The Gaussian is even and the explicit quadratic-loss
majorant bounds the remaining second, fourth and sixth
moments. Thus each witness certifies either endpoint alone,
or an expansion of the form
\[
 I_n(\xi)=L(\xi)+O(1/n)
     + k \sqrt n\,A(\xi)e^{nD(\xi)}(1+O(1/n)),              \tag{3}
\]
where the nonzero Gaussian amplitude \(A\) and its reciprocal
and field derivative are bounded effectively on the chart.
If \(k=0\) the endpoint alone suffices.
All leading coefficient functions are fixed by the profile
inside the chart and independent of the order.
The same saddle orientation is used throughout each chart.

## An effective cover of the pairs of fields

Consider the compact set of the profile closure and pairs
\((\xi,\xi')\in K^2\) with equal imaginary parts.
At every point of this set one of two kinds of certificates
exists:

1. At least one field has an endpoint-only witness or a witness
   with strict dominance, namely \(k\ne0\) and \(\Re D\)
   strictly positive or strictly negative. Such an expansion
   gives an explicit terminal threshold for no root at that
   field, uniformly on its cell.
2. Both fields are at a regular balance point.
   The uniqueness at fixed height proved above gives
   \(\xi=\xi'\), and \(\Re q>0\) there.
   One common two-term chart works on a convex field neighborhood
   of the pair, with \(k\ne0\) and \(\Re q\ge\kappa>0\) on the
   segment between the two nearby fields.

There is no third case. In particular, at equality with
\(\Im D<0,\sigma<1\), the source activation theorem gives
zero saddle coefficient. At \(\Im D=0\) equality was excluded
by injectivity and the boundary limit zero.
The crossing integer of a balance witness agrees with the
source coefficient, since an expansion valid on a complex
field neighborhood identifies it on the side of strict saddle
dominance. The individual saddle in case 2 is simple and the
condition \(\Re q\ge\kappa\) is a strict interval test.

In case 2, at an exact root the modulus equation from (3) is
\[
     j\Re D(\xi_j)=-\tfrac12\log j+C(\xi_j)+O(1/j),
              \qquad C(\xi)=\log|L(\xi)/(kA(\xi))| .
\]
The error, \(|C|\), and its horizontal derivative bound \(C_1\)
are computable on the cell.
Since the physical drift is
\(\xi_N-\xi_n=\Delta\log(N/n)\in\mathbb R\),
the real difference of \(D\) has modulus at least
\(\kappa|\xi_N-\xi_n|\).
Subtracting the two modulus equations gives
\[
          (N\kappa-C_1)|\xi_N-\xi_n|
                       \le C_2(1+\log n)/n.
\]
For \(\Delta\ne0\), the displacement is at least \(1/N\),
which contradicts the effective bound for large \(N\).
For \(\Delta=0\), direct subtraction gives \(\Re D=O(1/n)\),
inconsistent with \(n\Re D=-\tfrac12\log n+O(1)\).
This supplies a computable threshold for each case-2 cell.

For a finite algorithm, refine the retained-row cutoffs, remaining
mass intervals, and field grids, using (1) on tested contours.
The paired-field grid has separate real coordinates and one common
imaginary coordinate.
Enumerate rational polygon witnesses and their rational strict
margins in parallel with that refinement.
Check that every pair cell admits one of the two certificates,
using interval enclosures, and take the largest resulting terminal
threshold. The search can dovetail over mesh levels and proposed
witness lists, without waiting on a coarse cell that has no
certificate.
At every limiting profile and field pair, the applicable witness
above is valid on an open neighborhood, with retained cuts
matched on its compact contour domain.
A cutoff avoiding the finitely many possible slopes near its
positive boundary ensures this matching persists on a profile
neighborhood. The remaining phase errors tend uniformly to zero
there by (1).
Compactness supplies a finite subcover of those neighborhoods.
Sufficiently fine finite cells and interval precision therefore
certify the entire family, proving termination.


## Consequence for the whole subcritical and near-critical margin

Take \(\delta>0\) small enough for the
[positive-boundary entry](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_POSITIVE_ENTRY.md)
combined with the
[effective thermal endpoint theorem](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md).
Take \(\epsilon>0\) at most the positive and negative boundary
widths from those proofs and the
[negative-boundary exclusion](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md).
Use the interior algorithm for this rational \(\epsilon\).
The largest of the three computed thresholds gives one
computable absolute terminal \(N_*\) excluding common roots
for every ordinary one-private-linear profile with
\(0\le\sigma\le1+\delta\), at all physical fields.
The equal imaginary coordinates at the two orders put them
in the same one of the three sectors.

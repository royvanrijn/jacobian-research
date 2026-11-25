# An effective regular contour cover for fully shared profiles

Use the notation \(W,\rho,w\) of the
[uniform shared band proofs](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md).
In particular the exact consecutive field displacement is
\(w_N-w_n=-\Delta\log(N/n)\).
Fix positive rational \(W_-\le W_+\), a positive rational real field cap \(M\),
a rational interior strip margin \(\epsilon>0\), a rational thermal
margin \(\kappa>0\), and a rational saddle separation \(d>0\).
Use
\[
 K=\{w:|\Re w|\le M,\ -\pi+\epsilon\le\Im w\le-\epsilon\}.
\]

**Theorem.** A finite search computes a terminal threshold for
consecutive coprimality on each of the following regimes, with
\(W_-\le W\le W_+\) and both physical fields in \(K\).

- Fully shared supercritical profiles with \(\max\rho\ge1+\kappa\).
- Fully shared subcritical profiles with \(\max\rho\le1-\kappa\),
  whose two upper saddles are at least \(d\) apart at both fields.

The bounds are uniform in all ordinary rows and degree differences.
This is a computable regular-contour threshold, with the finite
search and constants unevaluated. It leaves the real strip
boundaries, the critical thermal transition, and the specified
subcritical fold neighborhoods outside its scope.

## Finite contour witnesses from actual limiting contours

Use the effective retained-row and diffuse-mass nets and the
compensated chord centering from the
[transversality proof](MULTIVARIATE_EFFECTIVE_SHARED_TRANSVERSALITY.md#computing-a-positive-margin).
The effective fiber rectangle encloses every upper saddle on
the target rectangle, for ordinary and limiting profiles.
The row approximation controls phases and their derivatives on
each bounded contour neighborhood avoiding the retained cuts:
approximate the curvature transform by its retained sum plus
the diffuse constant, and integrate with the chord normalization.
A small enough retained cutoff keeps every omitted pole outside
the neighborhood. The errors in the phase and its derivatives are therefore
computable.

For a fixed limiting profile, put \(p_L=\max_{d_i>0}(-a_i/d_i)\),
\(p_R=\min_{d_i<0}(-a_i/d_i)\), with the infinite endpoints
when the respective list is empty. Use the domain formed by
both open half-planes and the two real gaps \((p_L,0)\),
\((1,p_R)\). Circles and polygon loops in it have winding
number zero around every point on the outer removed rays:
each ray connects those points to infinity without crossing
the loop. Their winding number is constant on \([0,1]\),
which they also avoid. Consequently a loop with winding one
at \(1/2\) gives the exact
[Rice residue formula](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases)
at the actual profile, as long as it stays in the domain.
A compact witness loop in the limiting domain retains these
properties throughout a sufficiently small profile cell.
This makes the analytic equality of integrals a winding test.
The numerical search does not have to approximate the flow
homotopy of a contour.

For the two stated support bands, the finite contour constructions
are those of the
[supercritical proof](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#finite-contour-deformation-for-the-limiting-profile)
and the
[subcritical proof](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md).
At a regular field they provide finitely many oriented local
saddle traversals and strict-sublevel exterior chains.
When saddle heights differ, isolate the higher one. If its
coefficient vanishes, continue the descent to the lower
one and its strict sublevel.
At least one coefficient is nonzero.
In the supercritical case this follows from activation on the
equal-height graph and its complementary regions.
In the subcritical case it follows from the bounded-gap winding
obstruction, or the negative-real-ray obstruction for a purely
diffuse limiting profile, in the cited proof.
Near regular equal-height fields the two saddle rectangles and
the exterior pieces retain a common strict loss outside the
local segments. If only one saddle has a nonzero coefficient at
equal height, the zero coefficient cancels identical traversals
and the other saddle gives a one-term witness.

All paths in such a finite construction, except for the saddle
segments, are compact inside the domain and in the relevant
open strict sublevel sets. Approximate them by rational polygons.
For each retained saddle use a short straight segment in a
rational descent direction through its exact continued location.
A rational direction sufficiently close to the pure local
quadratic direction has strict negative real Hessian along the
segment. Its ends can be joined in the exiting strict sublevel,
so the crossing integer is unchanged.
Repeated traversals at one saddle use the identical segment
with orientation signs. In particular zero total crossing
coefficient cancels exactly even if that height is larger.
For every field in one small cell the certificate can use the
same polygons, with short connectors to the moving saddle
segments. Its center loop remains in the correct winding class
by strict separation from \(1/2\) and from the cuts.

Each condition on a witness is an interval test with a positive
margin. Enclose the one or two simple saddles by rational disks
using Rouché against their linearizations. A sufficient test
on the circle of radius \(r\) about \(q_0\) uses
\[
  |X(q_0)-w|+\tfrac12 C_2 r^2<|X'(q_0)|r,\qquad
             C_2\ge\sup |X''|.
\]
This also controls the saddle motion on a parameter cell.
A lower Hessian bound in the chosen direction and an upper
Taylor derivative bound certify the saddle segment loss.
Subdivide each compact polygon to test its strict phase
height bound with the derivative majorants from the nets.
The center-loop winding can be computed with a rational test
loop whose vertices approximate the certified segment ends,
since the segments stay away from \(1/2\).
These checks succeed on a small cell around each limiting contour
witness. The exact geometric construction gives existence
of such cells, rather than a prescribed computational flow time.

## Uniform expansion and the common-root certificate

On each witness the normalized Rice kernel has the uniform
Stirling expansion of the
[exact contour](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases).
Its constants are computable on the polygon and saddle
neighborhoods. All Gamma arguments divided by the order
and their row scale stay in sectors a positive distance
from the negative real ray and away from zero.
At a real gap one uses the positive-argument reflection
formula as in the actual contour proof.
The summed relative Stirling error is \(O((1+W)/j)\)
because \(\sum_i1/m_i\le W\). The leading row amplitude and
its reciprocal are bounded uniformly on the witness neighborhoods.
Indeed \(x_i/m_i=1+r_i(2q-1)\), and on these neighborhoods each
logarithm and its first few derivatives are bounded by a computable
constant times \(|r_i|\). The bounds follow by integrating the
logarithmic derivative in \(r_i\) from zero, or from a compact
enclosure for the finitely many retained types. Also
\(\sum_i|r_i|=\frac12\sum_i |d_i|/m_i\le W/2\).
Thus the amplitude bounds require no identification of the diffuse
row amplitude. For each actual profile the same amplitude is
used at both integer orders.
The sector constants are computable by the explicit
[DLMF Stirling remainder bounds](https://dlmf.nist.gov/5.11.ii).
All endpoint and exterior pieces are at strict real-height loss.
On each oriented saddle segment, use symmetric Gaussian
integration with odd linear-amplitude and cubic-phase
terms removed. The quadratic loss bounds the fourth and
sixth moment errors. This gives the relative \(O(1/j)\)
error, with a computable constant, for each retained
nonzero saddle term.

If a witness has only a dominant term, its nonzero leading
amplitude and the strict phase losses give a computable
zero-free threshold in the cell.
If two terms contribute at a regular balance, write the
uniform local expansion in the notation of the
[effective consecutive comparison](MULTIVARIATE_EFFECTIVE_SHARED_TRANSVERSALITY.md#use-in-a-certified-regular-contour-chart):
\[
 I_j(w)=A_1(w)e^{jP_1(w)}+A_2(w)e^{jP_2(w)}+R_j(w),
 \qquad |R_j|\le \frac Ej
       \bigl(|e^{jP_1(w)}|+|e^{jP_2(w)}|\bigr).
\]
The scalar normalization is nonvanishing.
The \(A_i\) here include the nonzero crossing integers and
are the same functions at the two orders.
Analytic saddle bounds give a computable lower and upper
amplitude bound and a bound for the real logarithmic derivative
of \(-A_2/A_1\).

For \(\Delta\ne0\), a case-2 certificate requires a single
convex field chart containing the two field cells and their
joining segment. In it require
\(|\Re (P_1-P_2)'|\ge c>0\).
In the supercritical band this slope is nonzero at every
regular balance by the real-part divider between the saddles.
In the subcritical band, the
[classification of vertical fibers](MULTIVARIATE_EFFECTIVE_SHARED_TRANSVERSALITY.md#the-exact-exceptional-geometry)
allows a zero slope only at an exceptional one-row profile on a
horizontal equal-modulus ray. That ray is inactive by the
[active-arc proof](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-active-arc-and-vertical-fibers):
its regular component goes to unbounded fields, whereas activity
is constant along it and the limiting root scales are bounded.
Such a point has a one-term witness with the other crossing
coefficient zero. Therefore the slope test is required only
where the slope is strictly nonzero, and it persists
on a smaller certified chart. The horizontal displacement
identity and the real logarithm comparison in the
transversality proof give an explicit threshold since
\(|w_N-w_n|\ge1/N\) for nonzero integer degree difference.
For \(\Delta=0\) require a positive lower bound for
\(|e^{P_1-P_2}-1|\) on the chart.
In the supercritical proof the imaginary difference along
the interior equal-modulus graph is strictly between its
boundary values \(0\) and \(\pi\), up to orientation.
On the subcritical active regular arc the
[covering argument](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md)
gives the same open bounds between the fold and negative
real boundary. Thus the phase test is strict at any
interior active regular balance. The local divided zero
equations then give the other explicit threshold.
Both tests can be imposed on the same sufficiently small
chart without changing its contour.

## Termination and exclusions

Consider pairs of fields of equal imaginary part in \(K\)
and the compact family of profiles with the stated bands.
For subcritical profiles retain the closed separation
condition at both fields.
The stationary trace discriminant from the
[effective inverse proof](MULTIVARIATE_EFFECTIVE_STATIONARY_INVERTIBILITY.md#a-finite-interval-certificate-for-a-shared-slab)
tests the squared saddle difference even at inadmissible
nearby profiles, so cells that strictly violate the
separation bound can be rejected effectively.
The real curvature maximum is tested by a grid on \([0,1]\)
with its uniform derivative bound. Thus a cell outside the closed
hypotheses has a strict exclusion certificate.

For each admitted limit profile and field pair, either
at least one field has a dominant-term witness, or both
lie at regular active balances. In the latter case they
coincide: the supercritical graph has only one balance
at each height, and the height along the subcritical
active arc is strictly monotone by the
[vertical-fiber comparison](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md#the-active-arc-and-vertical-fibers).
This last uniqueness uses the exact active-arc topology;
an equal-modulus arc with a zero saddle coefficient
does not enter the pair comparison.
A common convex chart therefore exists at the diagonal
pair, with all the strict slope, phase and amplitude
bounds verified above.

Dovetail over the profile cutoff, remaining diffuse-mass
mesh, paired field grids with common imaginary coordinate,
and enumerated polygon witnesses and margins.
Verify that each pair cell is inadmissible, has a
dominant-term certificate, or is contained in a common
two-saddle chart with the two consecutive-order tests.
The strict estimates persist on open neighborhoods in
the compact profile and field family, and matched retained
cuts keep their contour neighborhoods valid.
A finite subcover therefore exists. The interval net errors
go to zero uniformly on each of its finitely many contour
neighborhoods, so the finite search terminates.
The largest terminal bound certified by those cells proves
the theorem.

A subcritical fold violates the positive separation premise
of this construction. No uniform estimate as any of \(\epsilon,\kappa,d\) tends to zero is asserted.
A separate Airy or quartic analysis is needed in the
corresponding shrinking windows.

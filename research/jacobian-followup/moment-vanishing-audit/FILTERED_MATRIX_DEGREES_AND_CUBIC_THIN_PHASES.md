# Original matrix degrees and cubic thin phases

15 November 2026. Internal written proofs with selected exact controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_graph_hulls_and_matrix_degrees.py),
[receipt](graph_hull_degree_verification.json) and
[source audit](graph_hull_degree_source_audit.json).

The [cross-programme draft, Section1](../cross-programme-connections/PROOF_DRAFT.tex)
identified a missing conversion in the finite real-curve marker theorem.
That deduction is audited and proved here in the original degree filtration.
For **every** nonzero thin pure phase of matrix degree D, a fixed bad
marker has matrix degree at most D-8 on2O or D-12 on2I. This is an
all-degree original-coordinate bound, not just a normalized-curve bound.

The normalized polynomial has degree at least3 on either full quotient.
Its cubic topology, combined with the sphere Jordan exclusion, raises
this to at least4 in the original T-invariant algebras. Consequently
the thin-phase degree gates are24/36 on the full models and32/48 in
that restricted branch. These are necessary gates, not examples or
global positive MZ theorems. Both exceptional decisions remain OPEN.

## 1. The minimum representative filtration is multiplicative in degree

Give

    B=C[x,y,z,w]/(xw-yz-1)

the quotient filtration by minimum matrix-entry representative total
degree. For a finite right subgroup Gamma<=SL2 give A=B^Gamma the
induced filtration. Then

    gr B=C[x,y,z,w]/(xw-yz),
    gr A=(gr B)^Gamma.                                      (1)

For the first identity, the highest part of every nonzero multiple
(xw-yz-1)b is(xw-yz)b_top. Thus the initial ideal is exactly the
principal ideal(xw-yz). Its quadric is irreducible, so gr B is a domain.

For the second identity, the right action preserves every filtered
piece. Apply finite-group Reynolds averaging to the exact sequence
of consecutive filtered pieces. Taking invariants is exact in
characteristic zero, giving(1) in each degree, compatibly with
multiplication. In particular gr A is a domain too.

For nonzero f,g in A their highest symbols have nonzero product.
Therefore

    deg_matrix(fg)=deg_matrix f+deg_matrix g.                 (2)

If h is nonconstant, H=deg_matrix h>0, and p has degree d, the
powers h^j have distinct degrees jH. Its leading power cannot cancel
against lower powers, and

    deg_matrix p(h)=dH.                                     (3)

These statements use the minimum quotient filtration. Raw representatives
can have larger degrees: det has representative degree2 but equals1
and has filtered degree0. Arbitrary invariant-generator degree is a
different filtration and is not being substituted into(2)-(3).

## 2. A fixed bad marker below the phase in original coordinates

For either full exceptional ring, the
[normalization theorem](COMPLEX_HEIGHT_IMAGES_AND_CAUCHY_RIGIDITY.md)
writes a nonconstant thin phase as P=p(h), with h a regular real
self-adjoint seed. The
[finite-curve theorem](FINITE_CURVE_MARKERS_AND_POSITIVE_LOOP_MODEL.md)
proves that if P is nonzero and pure, one fixed marker h^j,
1<=j<d=deg p, has infinitely many nonzero marked moments.

Put D=deg_matrix P and H=deg_matrix h. The exact
[minimum-degree theorem](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md)
gives H>=8 on2O and H>=12 on2I. Equations(2)-(3) turn the marker
conclusion into

    deg_matrix(h^j)=jH<dH=D,
    deg_matrix(h^j)<=D-H
                  <=D-8 on2O, or D-12 on2I.                 (4)

The selected j is fixed for the entire infinite nonzero tail. No sign,
cutoff or specific j is asserted. There is no new existence claim.
This strengthens the previous normalized-curve-only caveat: the bound
now holds in the original matrix filtration on the full quotient.

In particular choose **any** fixed linear basis of the original invariant
functions of filtered degree<=D-8 or D-12. At least one basis marker
has an infinite nonzero tail: otherwise their finitely many eventual
cutoffs would have a common maximum, forcing eventual vanishing for
the bad seed-power marker in their span. Thus locating the normalized
seed is not required to obtain a complete finite marker list for this
architecture at each D. No finite moment cutoff is inferred here.

## 3. Linear and quadratic paths cannot be nonzero pure phases

A complex polynomial of degree<=2 maps a real interval to a line
segment or an injective arc. For a quadratic p(t)=a t^2+b t+c,

    p(t)-p(s)=(t-s)(a(t+s)+b).

If a,b are independent over R, there is no distinct real collision.
If dependent, the image lies in an affine real line. Lower degrees
are immediate. These images have connected complement, so the
retained positive Cauchy graph theorem excludes nonzero purity.
Consequently

    d>=3 on BOTH full exceptional models.                   (5)

It follows from D=dH that

    D>=24 on2O, D>=36 on2I.                                 (6)

At equality d=3 and H is exactly8 or12; the fixed bad marker is h
or h^2. This does not exclude other image-interior pure phases in
smaller matrix degrees.

## 4. Every cubic loop has only exterior tails

The following classification is used for **topology only**. Real-affine
changes of the complex plane do not in general preserve analytic
moments or corner angles.

Write a complex cubic as a real plane curve. If the cubic and quadratic
coefficient vectors are dependent, either all coefficient vectors are
collinear or a real linear coordinate kills both leading terms and
leaves a nonzero linear parameter. The image is a line or an injective
arc and cannot separate the plane.

Otherwise choose real-affine coordinates so one coordinate is quadratic
and the other cubic. A real parameter translation completes the square
in the first coordinate. Subtracting a multiple of it from the second
removes its quadratic term; rescaling gives

    x=t^2, y=t^3-beta*t, beta real.                          (7)

Equality of x for distinct parameters forces s=-t. Equality of y then
requires t(t^2-beta)=0. Thus there is a distinct collision only for
beta>0, at t=+/-sqrt(beta). The interval between these parameters
traces a regular Jordan loop, wholly in x<=beta. All parameter tails
outside that interval have x>beta and lie in the unbounded component.
There are no other real double points.

If a compact parameter interval fails to contain both double-point
parameters, its image is an arc. If it contains the loop and any
tail, that tail borders the unbounded complement on both sides.
The [exterior-edge theorem](HAAR_LOOP_CORNERS_AND_SPHERE_SYMMETRY.md)
excludes positive pure mass there. Therefore a nonzero pure cubic
polynomial path would have to be **exactly the Jordan loop**, without
tails. The affine coordinate changes in(7) only identify these topological
features; the Cauchy argument is still applied in the original complex
coordinate to the original pure measure.

## 5. The T-invariant cubic branch is closed in all seed degrees

The [polyhedral sphere theorem](HAAR_LOOP_CORNERS_AND_SPHERE_SYMMETRY.md)
excludes every simple Jordan compact image for a nonzero original
T-invariant thin pure phase on either exceptional model. Combined
with Section4 this excludes **every cubic normalized polynomial**,
without any bound on its real seed degree. Hence

    d>=4 in BOTH full original T-invariant algebras,
    D>=32 on2O, D>=48 on2I.                                 (8)

At equality d=4 and H is8 or12, with a bad marker among h,h^2,h^3.
No purity at either endpoint is asserted. The
[quartic interior-branch interval model](POSITIVE_GRAPH_HULLS_AND_INTERIOR_BALAYAGE.md)
shows why normalized degree4 is a genuine geometric transition for
abstract positive laws; no actual Haar seed has been exhibited.

The full non-T-invariant cubic branch is not closed by the sphere
argument. A real-affine deformation of a cubic loop can change its
corner opening, so the previously excluded particular right-angle
cubic does not exclude every cubic on a three-dimensional core.
Keep the full24/36 gates distinct from the restricted32/48 gates.

These results concern thin images of the PHASE. A thin auxiliary Hopf
height does not imply a thin associated phase. The one-domain harmonic
sweeping criterion and the Gaussian radial fixed-degree uniformity
remain separate structural tools; neither decides the full exceptional
MZ questions. No expanding degree or support campaign is required by
any result here.

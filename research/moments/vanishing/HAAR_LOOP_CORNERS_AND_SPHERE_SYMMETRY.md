# Haar loop corners and sphere symmetry

15 November 2026. Internal written proofs with selected exact geometry
controls. No novelty, independent specialist replay, external review or
formal verification is claimed. See the [checker](verify_haar_loop_corners.py),
[receipt](haar_loop_corner_verification.json) and
[source audit](haar_loop_corner_source_audit.json).

The [finite-curve marker theorem](FINITE_CURVE_MARKERS_AND_POSITIVE_LOOP_MODEL.md)
made any nonzero thin-image pure phase an automatic MZ failure. The next
gate was its existence for the actual Haar law. This note closes a global
geometric part of that gate, in every degree:

* On an N-dimensional compact smooth core, a pure polynomial Jordan loop
  with endpoint order k needs interior opening alpha>=2*pi*k/N at its join.
* The abstract cubic loop of the preceding note has opening pi/2 and
  endpoint orders1. It cannot come from **any C^2 real seed on any compact
  smooth core of dimension1,2 or3** with positive smooth volume. In
  particular it has no such Haar realization on either full quotient.
* For noncyclic rotation-invariant real seeds on S^2, an ordinary join
  is excluded. A further square-root-band and index argument closes the
  outward-cusp case for every rotation group that is neither cyclic nor
  dihedral. Thus **every simple Jordan image** is excluded for nonzero
  original T-invariant thin pure phases on both exceptional models.

This does not exclude non-Jordan separating graphs with several loops
or interior branches, or images with interior. The sphere exclusion applies to
original T-invariant phases. The dimension/corner bound also applies to
arbitrary original phases on the full three-dimensional compact cores.
It must not be replaced there by the stronger sphere conclusion. Both
exceptional MZ decisions remain **OPEN**.

## 0. A pure thin-image graph has no exterior hairs

There is also a constraint for arbitrary finite image graphs, including
graphs with several loops. Let a positive probability measure nu have
compact finite piecewise analytic graph support K, with positive C^1
open-edge densities after finite subdivision. If every positive analytic
moment vanishes, its Cauchy transform is1/s on the entire unbounded
component U of C\K. In particular0 cannot belong to U.

No open edge can border U on **both** sides. At any regular point of such
an edge away from0, both boundary branches of the transform equal the
same analytic function1/s. Their jump is zero. The local Plemelj formula
then makes the positive edge density zero, a contradiction. The jump
normalization and regularity hypotheses are retained from the
[Cauchy graph proof](COMPLEX_HEIGHT_IMAGES_AND_CAUCHY_RIGIDITY.md).

Thus every edge exposed to the exterior must have a bounded complementary
face on its other side. An exterior dangling branch is excluded even if
the rest of the graph contains loops. A bridge lying **inside** an outer
loop need not border U and is not excluded by this argument. This applies
to all original thin phases on both full exceptional quotients, through
their positive regular image measures, without T-invariance.

## 1. Pure positive measure on a Jordan curve is harmonic measure

Let gamma be a Jordan curve, D its bounded domain, and nu a positive
probability measure of full support gamma, with

    integral_gamma z^m dnu=0 for every m>=1.                    (1)

Then 0 belongs to D and

    nu=omega_D(0, .),                                         (2)

the boundary image of uniform circle measure under the Riemann map
phi:disk->D normalized by phi(0)=0.

First, if0 lies in the unbounded component, the exterior Cauchy transform
has germ1/s by (1). Analytic continuation in that component would give
1/s at0, where the integral transform is holomorphic, a contradiction.
Thus0 belongs to the closed domain.

By Mergelyan approximation on the full compact set closure(D), (1) implies

    integral_gamma f dnu=f(0) for every f in A(D),              (3)

where A(D) means holomorphic on D and continuous on its closure.
The approximation input is classical; a primary source stating it is
[Bishop--Lazebnik, Theorem10.7 and its proof](https://www.math.stonybrook.edu/~bishop/papers/Runge_Apr1_arxiv.pdf).
Only the ordinary uniform approximation conclusion is used.

Caratheodory's Jordan theorem gives a homeomorphism
psi:closure(D)->closed_disk, conformal inside. If0 were on gamma, set
zeta=psi(0) and f=(1+bar(zeta)*psi)/2. Its real part is at most1, with
equality on gamma only at0. Equation(3) and positivity would force
nu=delta_0, contrary to full support. Hence0 lies in D.

Normalize psi(0)=0 and use f=psi^m in (3). The positive circle measure
rho=psi_*nu has every positive Fourier coefficient zero; reality gives
every negative coefficient zero too. Uniform density of trigonometric
polynomials identifies rho with uniform circle measure, proving(2).
The conformal existence and boundary extension input is retained from
[Beliaev, Theorems2.3.1 and2.6.6](https://www.math.stonybrook.edu/~bishop/classes/math627.S25/Beliaev_notes.pdf).

## 2. The harmonic mass bounds used below

Fix a bounded Jordan domain D and interior point0. For z0 on its boundary,
let theta(r) be the total length of

    D intersection {|z-z0|=r}.

The Ahlfors distortion bound gives, with constants depending on fixed
normalizations,

    omega_D(0,B(z0,r))
       <= C exp(-pi integral_r^R dt/theta(t)).                 (4)

See [Bishop, Corollaries53--54, printed143--144](https://www.math.stonybrook.edu/~bishop/classes/math626.F08/rmt.pdf).
Scaling reduces to the normalization in that source. The outer portion of
the integral only changes C. Since theta(t)<=2*pi*t, (4) gives the general
Beurling bound

    omega_D(0,B(z0,r))<=C r^(1/2).                            (5)

If the boundary has interior opening alpha at z0, meaning

    lim_(r->0) theta(r)/r=alpha,

then for every beta<pi/alpha,

    omega_D(0,B(z0,r))<=C_beta r^beta.                        (6)

For alpha=0, (6) holds for every finite beta>0. Indeed choose
alpha'>alpha with beta<pi/alpha'; then theta(t)<=alpha'*t
for all sufficiently small t and integrate (4). No precise corner
asymptotic or equality case is asserted.

At a regular analytic boundary arc harmonic measure has bounded local
density in an analytic regular parameter. To see this, compose the
inverse Riemann map with that parameter, then apply Schwarz reflection
after a disk-to-half-plane Mobius transformation. Its derivative stays
bounded on a shorter arc. This is a local assertion away from vertices.

## 3. A dimension bound at the endpoint of a real seed

Let M be a connected compact N-dimensional smooth manifold **without
boundary**, with a probability measure having smooth strictly positive
density. Let h:M->R be nonconstant and C^2, h(M)=[a,b]. Suppose a complex
polynomial p maps [a,b] once around a Jordan curve: only a,b are identified
and p is injective on(a,b). Set P=p(h), and suppose all positive pure
moments of P vanish.

The pushforward has full support the loop, so Section1 identifies it
with harmonic measure. Let z0=p(a)=p(b), and suppose its interior opening
alpha exists. At a choose an extremum point of h. In local coordinates,
dh=0 and Taylor's theorem gives

    0<=h-a<=C |x|^2.

Positive source density on the ball of radius c*sqrt(epsilon) gives

    measure{a<=h<=a+epsilon}>=c' epsilon^(N/2).                (7)

If the endpoint order of p at a is k, then

    |p(t)-z0|<=C|t-a|^k.

Equations(2),(6) consequently give the upper bound

    measure{a<=h<=a+epsilon}<=C_beta epsilon^(k*beta)          (8)

for every beta<pi/alpha. If alpha<2*pi*k/N, choose beta with
N/(2k)<beta<pi/alpha. Bounds(7)-(8) contradict each other as
epsilon->0. The same argument applies at b. Thus necessarily

    alpha>=2*pi*k_a/N and alpha>=2*pi*k_b/N.                  (9)

An inward cusp alpha=0 is excluded in every finite dimension. Equality
in(9) is left undecided. No upper bound on the degree of p or h is used.
The manifold-without-boundary and C^2 hypotheses matter for(7).

For the earlier abstract cubic

    p(t)=t^2-1/2+i(t^3-t), -1<=t<=1,

the outgoing rays at z0=1/2 are(-2,2) and(-2,-2). The domain is to
their left and has opening alpha=pi/2; both endpoint orders are1.
Equation(9) forces N>=4. Its abstract positive harmonic interval law
therefore cannot be h_*volume for a C^2 seed on a compact core of
dimension<=3. This is an impossibility result for this loop architecture,
not a rejection of the preceding abstract measure construction.

In particular it cannot be realized by a regular real seed on either
SU(2)/2O or SU(2)/2I, or on their original sphere quotients. The same
statement excludes any real affine reparameterization of the seed
interval for this cubic. It does not exclude another polynomial loop.

## 4. Interior seed values have bounded local bands on a polynomial loop

Now let M=S^2 and h be real algebraic, so its critical sets have finitely
many semialgebraic components. Keep the one-turn Jordan hypothesis on p
and assume P is pure. For every fixed v in(a,b),

    measure{|h-v|<=epsilon}<=C_v epsilon                     (10)

for all sufficiently small epsilon.

If p'(v)!=0, p parameterizes a regular analytic boundary arc, and the
bounded local density statement after(6) proves(10). If p has order
k>=2 at v, the image of[v-epsilon,v+epsilon] lies in a disk of radius
C epsilon^k centered at p(v). Bound(5) then gives an upper bound
C' epsilon^(k/2)<=C'' epsilon. Injectivity on(a,b) ensures these
parameter bands are measured by the same boundary pushforward; overlaps
with another parameter interval are not being silently discarded.

The local band proof retained in the
[real-height obstruction](HOPF_HEIGHT_CRITERION_AND_REAL_OBSTRUCTION.md)
therefore applies at **every interior critical value**. In dimension2,
a degenerate critical point gives bands of mass at least c epsilon^(2/3)
or c epsilon^(5/6), according to Hessian rank. A nondegenerate saddle
gives at least c epsilon log(1/epsilon). All contradict(10).
Thus every interior critical point is a nondegenerate extremum. There
are finitely many, by semialgebraicity. This argument does not assume p
is immersed at every interior parameter and does not require globally
bounded endpoint density.

## 5. Noncuspidal joins force isolated endpoint extrema on the sphere

If the opening at z0 is alpha<2*pi, choose beta>1/2 in(6). Polynomial
Lipschitz continuity at a,b gives

    measure{h in[a,a+epsilon] union[b-epsilon,b]}
                                         <=C epsilon^beta. (11)

Suppose an endpoint level were positive-dimensional. A semialgebraic
level set on the sphere then contains a smooth curve segment. Along
that segment h equals its global extremum and dh=0. In a fixed tubular
chart, Taylor's theorem gives h-a<=C u^2 (or b-h<=C u^2), where u is
the normal coordinate. A strip of width c*sqrt(epsilon) therefore has
mass at least c' sqrt(epsilon), contradicting(11). Consequently both
endpoint levels are finite sets of isolated strict extrema.

This includes ordinary corners, smooth joins of opening pi, and inward
cusps of opening0. Only the outward-cusp opening2*pi escapes this
particular endpoint argument. Polynomial arcs have a limiting opening
at a join, allowing the two tangent directions to coincide; no angular
oscillation is introduced by a polynomial parameter.

## 6. All critical points being isolated extrema forces cyclic symmetry

Combine Sections4--5. The real algebraic h has finitely many critical
points, and all are isolated extrema, possibly degenerate at a,b. Each
has gradient index+1 on an oriented surface. For completeness this also
holds at a degenerate isolated strict extremum: in a small coordinate
disk with no other critical point, take a sufficiently small regular
sublevel component around a minimum. It is compactly contained in the
disk. A hole would contain a maximum above the boundary level, hence
another critical point. Thus it is a disk, and its gradient points
outward on the boundary, with degree1. At a maximum use -h; multiplication
by-1 on the target circle does not change degree. The boundary index
equals the local zero index by invariance of degree.

The Poincare--Hopf theorem now gives

    number of critical points=sum indices=chi(S^2)=2.        (12)

See [Wright--Dixon, Theorem15](https://websites.umich.edu/~alexmw/PoincareHopf.pdf)
for this standard index theorem. Hence h has one unique minimum and
one unique maximum.

If a finite noncyclic group K<=SO(3) preserves h, it must fix both
points individually. A rotation group fixing a point is cyclic. This
contradicts the hypothesis. We have proved:

    A nonconstant real algebraic K-invariant seed on S^2,
    with K noncyclic,cannot have a nonzero pure polynomial
    one-turn Jordan-loop phase whose join has alpha<2*pi.    (13)

The proof also shows that **any** surviving Jordan-loop pure phase on
such a sphere must have a positive-dimensional extremal level. Otherwise
Section4 and the isolated-extremum index argument still give(12),
regardless of its opening. Therefore its join must be an outward cusp.
At any endpoint corresponding to a positive-dimensional extremal level,
the order k of p is necessarily1: if k>=2, the general Beurling bound(5)
gives O(epsilon^(k/2))=O(epsilon), contradicting the strip lower bound
c sqrt(epsilon). Thus at least one such endpoint has p'!=0.

## 7. The square-root bound closes cusps for polyhedral symmetry

The preceding cusp possibility can be closed for O and I, without any
corner hypothesis. Polynomial Lipschitz continuity and(5) give

    measure{h in[a,a+epsilon] union[b-epsilon,b]}
                                            <=C sqrt(epsilon). (14)

We first classify the nonisolated endpoint extrema allowed by(14).
At a minimum write H=h-a in local coordinates. If its Hessian has
rank2, it is isolated. If it has rank1, the analytic splitting lemma
gives H=u^2+g(v), with g nonnegative real analytic. A nonisolated
zero forces g to have accumulating zeros, so g is identically0.
The endpoint level is then a smooth critical curve, nondegenerate in
its normal direction. The analytic splitting lemma is the usual
parameter form of the Morse lemma; it follows here by solving H_u=0
with the implicit function theorem, then factoring the remaining
quadratic difference and taking its positive analytic square root.

If the Hessian has rank0, nonnegativity makes the cubic Taylor term
zero. Let H_4 be its nonnegative homogeneous quartic Taylor term,
possibly identically0. Under x=epsilon^(1/4)y, Fatou's lemma gives,
for each fixed ball B_L,

    liminf measure{H<=epsilon}/sqrt(epsilon)
           >=c*area{y in B_L:H_4(y)<1}.                       (15)

The upper bound(14) implies that {H_4<1} has finite area in the
whole tangent plane. In polar coordinates that area is

    (1/2) integral_0^(2*pi) H_4(cos t,sin t)^(-1/2) dt.       (16)

A zero of a nonnegative analytic angular polynomial has even order
at least2, so(16) diverges if H_4 has a zero direction. It also
diverges if H_4 is identically0. Consequently H_4 must be positive
definite. Taylor's theorem now makes H>=c|x|^4 near the point,
so the extremum is isolated. In particular a rank-zero point cannot
belong to a nonisolated endpoint level. The same argument applies
at a maximum by replacing H with b-h.

Thus every nonisolated endpoint component is a smooth compact
critical circle with nondegenerate normal Hessian. There are no
crossings, endpoints or degeneracies on these circles. Compactness
and semialgebraicity give finitely many circles and isolated points.
All isolated critical points, including the interior ones from
Section4, are extrema and have gradient index+1.

A critical circle contributes total index0. Here is an explicit
reduction to the isolated-zero theorem. A circle in S^2 has trivial
normal bundle. In a tubular coordinate (t,u), factoring the normal
quadratic and changing u gives h=c+/-u^2. Perturb h near the circle
by a small Morse function of t with one maximum and one minimum,
using a normal cutoff. There are two new nondegenerate zeros, with
opposite determinant signs, so their indices sum to0. A sufficiently
small perturbation introduces no zeros in the cutoff annulus, where
the old normal derivative is bounded away from0. Carry this out in
disjoint tubes around all the critical circles, leaving the original
isolated points unchanged. Poincare--Hopf therefore gives

    number of original isolated critical points=2.            (17)

Any rotation group preserving h preserves this unordered pair.
A finite subgroup of SO(3) preserving a two-point set is cyclic or
dihedral: if their sum is nonzero, that vector is fixed and the group
is cyclic; if they are antipodal, the group preserves their axis,
with a cyclic rotation subgroup and at most one coset of half-turns.
Neither O nor I is of these types. More generally no finite rotation
group that is neither cyclic nor dihedral can preserve such an h.

We have therefore excluded **all one-turn polynomial Jordan-loop
pure phases on the polyhedral sphere models, including outward
cusps**, in every degree. The example h=z^2 explains the dihedral
boundary: it has one critical minimum circle and two isolated maxima,
and a noncyclic dihedral rotation group preserves that pair. This
example is not asserted to give any pure polynomial phase.

## 8. Original exceptional scope and what remains

For an original left-T-invariant thin phase on either exceptional model,
the [real-curve normalization theorem](COMPLEX_HEIGHT_IMAGES_AND_CAUCHY_RIGIDITY.md)
supplies a regular self-adjoint seed in that same T-invariant algebra.
On T\SU(2)=S^2 it is real algebraic and invariant under O or I, both
neither cyclic nor dihedral. Section7 excludes this complete all-degree
one-turn Jordan subclass. No finite support window or degree bound is
imposed.

In fact the one-turn condition is automatic when the normalized image
is a Jordan curve. The normalization construction has
C(P,P^dagger)=C(h), so a rational expression in p(h),bar(p)(h)
recovers h outside finitely many exceptional values. Consequently p
is injective on the real interval away from finitely many collisions.
Compose p with a homeomorphism from its Jordan image to the circle
and lift to a continuous real angle on[a,b]. If that angle reversed
direction, an interval of angles would have two preimages, contrary
to generic injectivity. It is therefore strictly monotone. Surjectivity
and generic injectivity make its total variation exactly2*pi: the
endpoints are identified and the interior is injective. Thus **every
simple Jordan compact image** of a normalized original T-invariant
thin pure phase is excluded on both exceptional models.

For arbitrary full phases use the three-dimensional compact core and
the general bound(9). The index argument on S^2 does **not** apply to
arbitrary phases or their seeds. In particular a reflex Jordan corner
is not excluded on the full quotient merely by the sphere result.

The thin-image existence frontier still contains non-Jordan separating
polynomial graphs without exterior hairs, including several loops or a
single loop with interior branches, and full non-T-invariant
architectures compatible with(9). Simple Jordan images are excluded
in the original T-invariant subalgebras; they are not all excluded in
the full quotients. Existence is not asserted for any remaining case.
If a nonzero pure thin
phase is found, the finite-curve theorem already supplies a fixed bad
seed-power marker. The full image-interior marked problem remains open.

This is a pure-phase theorem using the UNSYMMETRIZED image. The even
Hopf-height problem uses a SYMMETRIZED measure with a different target
moment sequence. It cannot be replaced by (1)-(2), and a thin height
does not imply a thin Hopf phase. Retain the full S_K/all-marker problem,
all global triangle conditions and both exceptional open MZ decisions.

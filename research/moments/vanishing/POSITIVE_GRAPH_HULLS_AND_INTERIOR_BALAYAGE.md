# Positive graph hulls and interior harmonic sweeping

15 November 2026. Internal written proofs with selected exact controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_graph_hulls_and_matrix_degrees.py),
[receipt](graph_hull_degree_verification.json) and
[source audit](graph_hull_degree_source_audit.json).

Every nonzero pure thin-image phase on either full exceptional quotient
has a **single outer Jordan domain**. The positive measure may charge
additional branches and loops inside that domain. Purity is exactly a
harmonic sweeping balance on its outer boundary; it does not identify
the entire measure with boundary harmonic measure.

This excludes separate exterior lobes as well as the previously excluded
exterior hairs. A disk boundary-cap estimate also constrains every
real-seed parameter sent to a regular outer boundary arc. An explicit
quartic polynomial path with an interior branch admits a positive pure
interval law, showing that the remaining geometry is real and cannot
be excluded by positivity and curve normalization alone. No actual Haar
realization of that new law is asserted. Both exceptional MZ decisions
remain **OPEN**.

## 1. Filling a finite graph leaves a finite collection of Jordan lobes

Let K be a connected finite piecewise analytic embedded planar graph
with at least one edge (the trivial measure delta_0 is outside this claim).
Subdivide at its finitely many vertices and crossings. Let U be its
unbounded complementary component and let

    H=C\U=K union all bounded complementary components.       (1)

This is the polynomial hull of K in one complex variable. The term
"filled graph" below means the explicit set in(1); no higher-dimensional
hull theorem is needed.

Suppose no open edge of K borders U on both sides. Then H is a finite
union of closed Jordan domains H_j=closure(D_j), whose open interiors
are disjoint and whose pairwise intersections are at most single
vertices. Their incidence through these vertices is a tree.

Here are the finite planar details. The outer frontier B=boundary(U)
is a connected finite plane graph: its boundary walk runs along the
outer face of connected K. Every edge of B is incident with U. A bridge
of B would have U on both sides and is therefore excluded. Split B at
cut vertices until its pieces are 2-connected. Every face of each piece
is bounded by a cycle. If a piece were not itself a cycle, it would
have an edge inside its outer boundary cycle, which could not be
incident with U. Thus the pieces are cycles. Their bounded disks have
disjoint interiors; a cycle inside another would not be exposed to U.
The cut-vertex splitting gives their tree incidence. Filling these disks
gives exactly H, while any remaining edges of K lie inside the disks.

These are the usual finite plane-graph boundary facts. See
[Diestel, Sections4.1--4.2, especially Proposition4.2.6](https://www.math.uni-hamburg.de/home/diestel/books/graph.theory/preview/Ch4.pdf).
Piecewise analytic finite arcs can be straightened for this purely
topological argument. The disk boundaries in the original embedding
remain Jordan curves, with their original corners or cusps.

## 2. Positivity forces exactly one outer lobe

Let nu be a positive probability measure of full support K, with positive
C^1 open-edge densities after finite subdivision, and suppose

    integral z^m dnu=0 for every m>=1.                         (2)

The [exterior-edge theorem](HAAR_LOOP_CORNERS_AND_SPHERE_SYMMETRY.md)
makes its Cauchy transform1/s on U. Thus0 belongs to H and exterior
hairs are absent, so Section1 applies.

Mergelyan approximation on the full compact set H gives

    integral f dnu=f(0) for every f in A(H),                   (3)

where A(H) consists of continuous functions on H holomorphic on its
interior. The ordinary approximation conclusion is the one retained
from [Bishop--Lazebnik, Theorem10.7](https://www.math.stonybrook.edu/~bishop/papers/Runge_Apr1_arxiv.pdf).

We need a positive-real function localized to a lobe H_j. Let S_j be
its finite set of shared boundary vertices, together with0 if0 is on
its boundary. Choose a conformal homeomorphism psi_j:H_j->closed_disk.
For the finite points zeta_l=psi_j(S_j), define inside the disk

    T(w)=sum_l (zeta_l+w)/(zeta_l-w),
    F(w)=1/(1+T(w)).                                         (4)

For empty S_j take F=1. The real part of T is positive inside; on the
circle away from its poles T is purely imaginary. Thus F is continuous
on the closed disk, vanishes at the listed points, and has strictly
positive real part elsewhere. Extend F composed with psi_j by0 on
the other lobes. Vanishing at S_j makes this a function in A(H).

If0 were on the outer boundary, use a lobe containing it. Then f(0)=0,
whereas the real part of its integral is positive by full support on
an open boundary edge, contradicting(3). Hence0 belongs to one D_0.
For any other lobe repeat the construction: f(0)=0 and its integral
again has positive real part. There can be no other lobe. Consequently

    H=closure(D), D a bounded Jordan domain, 0 in D.           (5)

In particular a figure-eight with separate exterior lobes cannot be
the full-support image of a pure phase. Nested loops and branches
inside a single outer domain are not excluded by(5).

All original thin-image Haar phases satisfy the measure hypotheses,
by the retained real-seed normalization and regular arc-density proof.
Thus(5) applies to both **full** exceptional quotient rings, without
T-invariance or a degree bound.

## 3. Purity is an exact harmonic sweeping equation

Normalize psi:D->disk by psi(0)=0. Write nu=nu_boundary+nu_interior.
For z in D, let omega_z be harmonic measure on boundary(D), and let
omega_0 denote the measure at0. Then(2) is equivalent to

    omega_0=nu_boundary+integral_D omega_z dnu_interior(z).    (6)

The rightmost operation sweeps interior mass to the boundary, often
called balayage. To prove(6), approximate any real continuous circle
boundary function by real trigonometric polynomials. Its harmonic
extension is then uniformly approximated on the closed disk by real
parts of analytic polynomials. Compose with psi and use(3). Equality
for all continuous boundary test functions gives(6). Conversely(6)
applied to the real and imaginary parts of z^m proves(2).

The standard conformal/Dirichlet interpretation of harmonic measure
is described in [Beliaev, Section4.2](https://www.math.stonybrook.edu/~bishop/classes/math627.S25/Beliaev_notes.pdf).
The balance and converse above are supplied here.

In disk coordinates let m be normalized circle measure and rho the
interior pushforward. The swept density is

    kappa(xi)=integral_disk (1-|w|^2)/|xi-w|^2 drho(w).

Equation(6) says exactly

    0<=kappa<=1 almost everywhere,
    psi_*nu_boundary=(1-kappa)m.                             (7)

There is no cancellation of positive interior mass; it consumes part
of the available harmonic boundary density. Equation(7) is necessary
and sufficient for a proposed positive interior measure plus its forced
boundary remainder to represent evaluation at0. Full support imposes
the additional positivity on all actual graph edges. This is an
all-order period criterion, not a finite moment test.

When the whole support is just the Jordan boundary, rho=0 and the
earlier uniqueness result is recovered. With interior branches, rho
need not be zero. Treating(6) as nu=omega_0 would be incorrect.

## 4. A universal boundary-cap bound and a seed consequence

For any positive probability rho_total on the closed disk with all
positive analytic moments zero, any |xi|=1 and 0<r<=1,

    rho_total{|w-xi|<=r}<=4r.                                (8)

Indeed the positive harmonic function

    J(w)=Re((1+(1-r)bar(xi)w)/(1-(1-r)bar(xi)w))

has integral1 by its convergent power series and the moments. On the
cap its numerator in the Poisson expression is at least r and its
denominator at most4r^2. Thus J>=1/(4r), proving(8).

At any regular analytic arc point z0 of the outer boundary, psi
extends holomorphically by Schwarz reflection and has bounded local
derivative. Pulling(8) back gives

    nu(B(z0,r))<=C r,                                       (9)

including interior graph mass in that ball, not just boundary mass.

Suppose nu is the law of P=p(h) for a C^2 real seed h on a compact
positive smooth N-dimensional manifold without boundary, h(M)=[a,b]. If p(t0) is on such
a regular outer arc and p has order k at t0, then

    measure{|h-t0|<=epsilon}<=C epsilon^k.                    (10)

The parameter band is simply mapped into a ball of radius C epsilon^k;
no injectivity or absence of interior branches is assumed. On a core
of dimension<=3, k>=2 is impossible: at a regular source point the
band has mass>=c epsilon by submersion coordinates, while at a critical
point it has mass>=c epsilon^(N/2) by the extremum-ball argument
applied to |h-t0|. Both contradict(10) for k>=2 and N<=3.

On the sphere, (10) with k=1 also excludes every degenerate critical
point and every nondegenerate saddle mapped to a regular outer arc,
by the retained local band estimates. Such critical points must be
nondegenerate extrema. These statements constrain **all** original
thin phases and their interior branches. They do not control an
arbitrary interior graph point or a singular outer boundary vertex.

## 5. An irreducible polynomial path with an inward branch really works abstractly

Positivity and normalized polynomial-curve geometry alone do not exclude
the remaining one-domain graphs. Consider

    p(t)=t^4-(3/2)t^2+1/4 + i*(t-t^3).                       (11)

On[-1,1] this is a regular Jordan loop with only endpoints identified
at-1/4. To check this, put u=t^2,v=s^2. Equality of real parts implies
u=v or u+v=3/2. In the first case equality of imaginary parts gives
t=s or the endpoint pair. In the second case the imaginary parts
must have the same sign. Squaring them gives uv=1/4; the two roots
(3+/-sqrt(5))/4 cannot both lie in[0,1]. Also its derivative components
4t^3-3t and1-3t^2 have no common real zero.

The loop crosses the real axis only at-1/4 and1/4, so0 lies in its
bounded domain D. Its signed area is12/35. The outgoing join rays
are(-1,-2),(-1,2); the interior is the reflex sector containing the
positive real direction. Continuing the parameter slightly below-1
therefore puts the additional arc strictly **inside** D, with no
new intersection for a sufficiently small delta>0. Choose delta<1/4.
Its full image on I=[-1-delta,1] is one outer Jordan loop and a
single interior branch. The new leaf lies inside D.

Let psi:D->disk be normalized at0 and put an unnormalized positive
measure on the interior branch by

    d sigma(p(t))=(1-|psi(p(t))|)dt,
                       -1-delta<=t<-1.

Every point of the open branch has positive density. Its swept density
with respect to omega_0 is bounded above by

    integral_(branch) (1-|w|)*(1+|w|)/(1-|w|) dt
                                  <=2delta<1/2, w=psi(p(t)).

Thus beta=omega_0-sweep(sigma) is a positive boundary measure with
beta>=omega_0/2. The measure nu=beta+sigma is a probability measure
of full support the connected polynomial graph and satisfies(6),
hence every positive analytic moment is zero. Pulling it back along
the loop and branch parameters gives a positive full-support interval
law mu on I with

    integral_I p(t)^m dmu=0 for every m>=1.                  (12)

The duplicate join point has zero mass, so the pullback is well defined
without endpoint atoms. This example has the same irreducible normalized
curve architecture as the thin-phase normal form: its real and imaginary
coordinate degrees are4 and3, so the degree of C(t) over C(x(t),y(t))
divides both4 and3 and is1. Since t is integral over the cubic y(t),
the normalization is C[t]. It is not a reducible circle-plus-segment
example.

By the [finite-curve marker theorem](FINITE_CURVE_MARKERS_AND_POSITIVE_LOOP_MODEL.md),
some fixed marker among t,t^2,t^3 has infinitely many nonzero marked
moments for this same law. No particular marker, sign or cutoff is
claimed. No actual exceptional Haar seed or Haar realization of mu is
asserted. This closes a positivity-only attempt to remove interior
branches; an actual Haar exclusion must use more information about
the seed distribution.

## 6. The remaining exceptional question

For original T-invariant phases the outer boundary is Jordan but the
**whole support** must have additional interior graph structure, since
the [sphere theorem](HAAR_LOOP_CORNERS_AND_SPHERE_SYMMETRY.md) excludes
the boundary-only case. Separate exterior lobes, exterior hairs and
ordinary critical parameters sent to regular outer arcs are excluded
in their stated scopes. The unsolved gate is the actual Haar sweeping
balance(6) for a non-Jordan support in one filled domain.

For arbitrary full phases keep this same one-domain criterion and
the three-dimensional corner restrictions where applicable; do not
extend the stronger sphere index theorem to those phases. The full
image-interior/all-marker problem and the separate symmetrized Hopf
height problem remain in scope. Both exceptional MZ decisions are OPEN.

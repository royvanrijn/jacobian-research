# Complex height images and Cauchy-transform rigidity

15 November 2026. Internal written proofs with exact algebraic controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_complex_height_images.py),
[receipt](complex_height_image_verification.json) and
[source audit](complex_height_image_source_audit.json).

The remaining complex-height problem has an all-degree geometric
division. A height whose real and imaginary parts are dependent is
a complex polynomial in one **regular real-valued** seed. If its
symmetrized image does not separate the plane, the required interval
moments force it to be real. Noncyclic polyhedral symmetry then
excludes it by the earlier sphere-density argument. In particular,
all complex line-valued heights are excluded, regardless of degree.

What remains is an image with interior, or a unary polynomial image
whose symmetrization contains a closed loop. Neither possibility is
proved to give a pure Hopf phase. Both exceptional MZ questions
remain **OPEN**.

## 1. Reality and one-dimensional regular images

Let A be a finitely generated normal C-domain with unirational
spectrum and A^*=C^*. Suppose an antilinear involution dagger on A
restricts to complex conjugation on a Zariski-dense compact real
core M, generically totally real of maximal dimension. These
hypotheses hold for C[SL_2/Gamma] and its original left-torus-invariant
subalgebra, with a compatible unitary realization of Gamma.

For a nonconstant a in A the following are equivalent:

    a,a^dagger are algebraically dependent;
    a=p(h), where h in A is regular,h^dagger=h,
                            and p in C[T] is nonconstant.   (1)

Proof. Set A_0=C[a,a^dagger] and normalize its affine curve in its
fraction field. The normalized ring B lies in A: every element is
integral over A_0 and belongs to Frac A, while A is normal. The
curve is unirational because it is dominated by Spec A. Its smooth
projective completion therefore has genus zero, by restricting a
rational parametrization to a general line and applying
Riemann--Hurwitz. Its affine normalization is P^1 minus finitely many
points. It has only constant units since B subset A, so exactly one
point is removed. Thus B=C[h_0]. This is the same normalization
mechanism as the retained [unary descent proof](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md).

Dagger preserves B and has h_0^dagger=alpha*h_0+beta with
alpha*bar(alpha)=1. An affine change of h_0 makes its new parameter
h self-adjoint. Explicitly choose c with bar(c)*alpha=c; the involution
condition makes bar(c)*beta purely imaginary, and an appropriate
constant shift solves d-bar(d)=bar(c)*beta. Then h=c*h_0+d is
self-adjoint. Since a belongs to B, it is p(h). The converse is
immediate because a^dagger=bar(p)(h).

On the compact core, (1) is also equivalent to the real differential
of a having rank at most1 everywhere on its generic smooth locus.
Indeed rank at most1 makes da wedge da^dagger vanish on the real
tangent spaces. Their complex spans are the full complex tangent
spaces; Zariski density and the characteristic-zero differential
criterion give algebraic dependence. Conversely polynomial dependence
places the image in a real algebraic curve, so rank2 is impossible.
If they are independent, rank2 occurs somewhere and the image has
interior by the inverse function theorem.

For a=p(h) in (1), h(M) is a nondegenerate compact interval I, and

    a(M)=p(I).                                              (2)

The image is a finite graph of real analytic arcs, possibly with
overlaps or closed loops. A polynomial path need not be an arc.

## 2. Regularity of the positive image measure

Use Haar's positive smooth probability on the compact core. For a
nonconstant real regular h, the critical values on that core form
a finite set: the critical set is semialgebraic with finitely many
connected components, and h is constant on each such component.
On any compact interval of regular values, a submersion-chart
partition of unity gives a smooth positive pushforward density.
Every individual level has measure zero, since a nonconstant real
analytic function cannot be constant on an open set. Hence there
are no atoms.

Passing through the polynomial p changes the density by the ordinary
one-variable change of coordinates on each regular branch. After
subdividing at its finitely many critical values and arc intersections,
the resulting measure has a smooth positive density on each open
edge of the image graph. Density may diverge integrably at the finite
vertex set. These statements also hold for the symmetrized measure
and image. The required semialgebraic finiteness is classical; see
[Coste's notes](https://perso.univ-rennes1.fr/michel.coste/polyens/SAG.pdf).

## 3. A positive graph measure with the interval moments

Let mu be a compact positive probability measure supported on a
finite piecewise analytic graph K in C. Assume its open-edge
densities are C^1 and positive, apart from a finite subdivision;
finite vertex atoms are allowed. Suppose C\K is connected and
the analytic moments are those of Uniform[-1,1]:

    integral z^(2n) dmu=1/(2n+1),
    integral z^(2n+1) dmu=0.                               (3)

Then

    mu=Uniform[-1,1].                                      (4)

Here is a proof retaining the positivity requirement. Its Cauchy
transform C(s)=integral dmu(z)/(s-z) has the germ

    C(s)=(1/2)log((s+1)/(s-1)),
    C'(s)=R(s)=-1/(s^2-1).                                 (5)

Connectedness makes the rational derivative identity valid on all
of C\K. It also forces both poles+1,-1 of R to belong to K.

At an interior point of a smooth edge, away from these poles, R
has a local primitive. The two one-sided branches of C differ by
a constant. Connect the two sides by a path in C\K and close it
through the small disk. The residue theorem gives their difference
as pi*i times an integer: the residues of R are+1/2,-1/2.

For an oriented arc z=gamma(t) carrying mass w(t)dt, the Cauchy
jump formula is, up to its orientation sign,

    C_+-C_-=2*pi*i*w(t)/gamma'(t).                         (6)

The local formula follows by subtracting the density at the point;
the remainder has equal boundary limits and the constant part is
the elementary Cauchy integral. See also
[Tao, Proposition1.1](https://www.math.ucla.edu/~tao/247a.1.06f/notes4.pdf)
for the straight-line normalization. Thus (5)-(6) make
w(t)/gamma'(t) a real half-integer. Since w is positive, each
mass-carrying edge is horizontal and has constant linear density.

After subdividing overlaps, mu is a finite sum of constant-density
horizontal segments plus vertex atoms. Its differentiated Cauchy
transform is now rational. A segment's endpoints give simple poles;
an atom gives a double pole. Equality with R excludes all atoms,
and all density jumps except at-1,+1. The only remaining positive
segment is[-1,1], with density1/2. This proves (4).

The same proof has a second useful form. If instead the germ is
C(s)=1/s, all residues of C'(s)=-1/s^2 vanish. Every edge jump is
zero and there is no edge mass. The finite atomic part must be
delta_0. Thus for these graph measures,

    integral z^n dmu=0 for every n>=1
                             implies mu=delta_0.            (7)

Neither result is asserted for a graph that separates the plane,
nor for an image with interior. On different complementary domains
the exterior germ need not determine the other Cauchy branches.

## 4. Consequences for every exceptional complex height

Let tau be ANY original T-invariant regular function on either
exceptional quotient, not just a chosen finite coefficient family.
It lies in the [complete global triangle algebra](GLOBAL_TRIANGLE_ALGEBRAS_AND_HAAR_BOUNDARIES.md).
Let nu be its positive compact pushforward and form

    mu=(nu+(-1)_*nu)/2,
    K=tau(S^2) union -tau(S^2).                             (8)

If tau has the required even moments1/(2n+1), mu has (3).
If tau has one-dimensional image and C\K is connected, Sections1-3
apply and give mu=Uniform[-1,1]. Positivity implies nu<=2mu, hence
tau is real-valued on the compact sphere and its law has bounded
density. The proved [noncyclic rotation obstruction](HOPF_HEIGHT_CRITERION_AND_REAL_OBSTRUCTION.md)
then excludes it.

Consequently every possible exceptional flat-moment height must
belong to one of these two architectures:

* tau,tau^dagger are algebraically independent; its compact image
  has interior;
* tau=p(h) for a regular real self-adjoint h, but the symmetrized
  polynomial path p(I) union -p(I) separates the plane.

For a finite planar graph separation is equivalent to containing a
simple closed loop. This is a necessary alternative, not an
existence or sufficiency theorem for either case.

In particular **no compact complex line-valued height works**.
Its symmetrized support is a union of two parallel line segments,
with connected complement. The argument does not require the line
to be the real axis, to pass through zero, or to have a prescribed
polynomial degree. Thus every affine complex function of every
regular real polyhedral-invariant seed is excluded.

The symmetrization in (8) is essential. Individual paths and their
unions can have different topology. A polynomial image can actually
produce the allowed separating geometry: for

    p(t)=t^2-1+i*t, -1<=t<=1,

p(I) and -p(I) are the parabolic arcs x=y^2-1 and x=1-y^2.
They meet at+i,-i and enclose an area8/3. This is only a geometry
example, NOT a flat-moment height or an exceptional counterexample.
It prevents silently excluding every unary polynomial path.

## 5. A full-phase consequence and what is still missing

For any regular phase P on SL_2/Gamma with one-dimensional compact
image K and connected complement, the positive pushforward obeys
the hypotheses of (7). If all pure moments vanish, P must vanish
on the compact core and hence be zero in the coordinate ring.
This conclusion includes arbitrary original phases under the stated
image hypothesis; it does not require T-invariance.

By (1), every other nonzero pure phase on either full quotient must
therefore have an image with interior, or be a polynomial in a regular
real-valued seed whose UNSYMMETRIZED image graph separates the plane.
For this pure-phase statement the graph is P(M), not (8).

No all-marker safety assertion follows merely from either remaining
architecture. Known safe pure phases can have image with interior.
As an exact control on actual source functions,the two canonical
Bezout heights lie in the rank-two case:the original D,E derivative
bracket with their compact conjugates is nonzero in each model.
This does not restore their purity,which was already excluded by
their nonzero first Hopf moments.
The full S_K problem, mixed-weight regular sections and both
exceptional MZ outcomes remain open. The geometric split narrows
the unknowns without replacing the full objective by a smaller
finite-degree or support calculation.

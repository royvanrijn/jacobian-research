# Unramified boundary curves have residue degree one

HC4MYBR0. 17 November 2026. Actual complex polynomial plane Keller maps,
with the entire original finite normalization retained.

For EVERY omitted unramified boundary prime D of the canonical entire
etale completion, its map onto its reduced target image Gamma has
residue-field degree ONE, and

    normalization(D) = normalization(Gamma) = A1.

The induced map between these two normalizations is an isomorphism.
The finite map D->Gamma is birational. No common source coordinate,
smooth algebraic boundary, root chart, or generic-degree bound is needed.

The target curve Gamma is a component of the actual nonproperness set.
By the separately retained no-affine-line-component theorem of Chau,
Gamma is singular. Its polynomial normalization is an immersion exactly
when D is smooth. Thus, whenever the individual-coordinate/root-chart
criterion supplies boundary smoothness, the target normalization is a
polynomial immersion with at least one self-identification. Tangencies
and higher contact between distinct smooth branches remain possible.

Residue degree one is not a claim that Gamma is an embedded affine line,
that all its singularities are nodes, or that the Keller map is invertible.

## 1. Actual finite normalization and the two curve normalizations

Let K=(P,Q):X=A2->Y=A2 be an actual complex Keller map. Let
Z->Y be the FULL finite normalization in the source function field,
and U its ENTIRE etale locus. The canonical construction and the
[boundary topology theorem](CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md)
give

    X open in U open in Z,
    D an irreducible component of (U minus X)_red,
    D is its whole closure in the affine Z,
    nu_D:Dnu->D finite, Dnu=A1,
    nu_D is bijective and D is analytically unibranch.    (1)

In particular D is closed in Z, so restriction of the finite morphism
pi:Z->Y makes

    g:D->Gamma=pi(D)

finite and surjective. Give Gamma its reduced irreducible curve
structure. Its normalization

    nu_G:Gammanu->Gamma

is a smooth affine curve, finite over Gamma. The residue degree is

    f=[C(D):C(Gamma)]>=1.                               (2)

No smoothness of D or Gamma has been assumed.

The composite Dnu->Gamma is finite. It factors uniquely through Gammanu,
because Dnu is normal and the map is dominant. The resulting map

    h:Dnu=A1->Gammanu                                  (3)

is finite, not merely generically finite. Explicitly, the integral closure
of C[Gamma] in C(Gamma) embeds into C[t]: its elements are integral over
C[Gamma] and hence over the normal ring C[t]. Also C[t] is already finite
over C[Gamma], so the same module generators make it finite over the
intermediate ring C[Gammanu]. Thus h is finite and surjective, of degree f.

## 2. Select the correct component of the etale fibre product

Form

    T=U times_Y Gammanu.

The morphism T->Gammanu is etale. Both schemes in the fibre product are
affine, so T is affine; it is a smooth curve. It may be disconnected.
We do NOT assume that T is finite over Gammanu.

The two maps Dnu->U and h:Dnu->Gammanu give a natural morphism

    alpha:Dnu->T.                                      (4)

It is finite. Indeed h is finite and T is separated over Gammanu:
the graph of alpha is a closed immersion into Dnu times_(Gammanu) T,
whose projection to T is finite by base change.

Its image is one irreducible component T_D of T. The image is closed
and one-dimensional, so it is a whole curve component. To verify that
alpha is birational onto this component, restrict to the dense normal
open of Gamma. There nu_G is an isomorphism, hence T->pi^-1(Gamma)
is an isomorphism. At the generic point of D the other normalization
nu_D is likewise an isomorphism. Thus alpha induces the identical
function field C(D) on its image.

The curve T_D is smooth and normal. A finite birational map onto a
normal integral curve is an isomorphism: on affine rings its integral
source algebra lies in the same fraction field and must equal the
integrally closed target algebra. Consequently

    Dnu=A1 is isomorphic to T_D.                        (5)

The irreducible components of a smooth curve are disjoint open and
closed components. Thus restriction of T->Gammanu to T_D is etale.
Together with Section1, (3) is a FINITE ETALE surjective morphism.

This is the step that handles both cusps and disconnected fibre products.
The normalization is performed on the target curve before using the
etale ambient map. Smoothness is claimed for T_D and Dnu, not for D.
Other components of T can be nonfinite over Gammanu and play no part
in the finite-cover argument.

## 3. Euler characteristic forces the residue degree to be one

A finite etale morphism of smooth affine complex curves is a finite
topological covering. Euler characteristic is multiplicative, so

    1=chi(A1)=f*chi(Gammanu).

The Euler characteristic on the right is an integer and f is a positive
integer. Therefore f=1 and chi(Gammanu)=1. A finite etale map of degree
one is an isomorphism, yielding

    h:Dnu -> Gammanu is an isomorphism,
    Gammanu=A1,
    C(D)=C(Gamma).                                      (6)

The Euler and finite-cover facts here are the same classical topological
foundations already used in the canonical completion theorem. No affine
curve classification or common-pencil theorem is imported.

Since D->Gamma is finite with the identical function field, it is
finite birational. It can still identify distinct points of D. The
normalization diagram does not make that map an isomorphism when
Gamma is singular.

## 4. Target branches, cusps and immersion

Identify Dnu and Gammanu using (6), and choose a parameter t. The target
map is a finite polynomial normalization

    beta:A1_t->Gamma subset Y,
    beta(t)=(p(t),q(t)).                                (7)

For every target point gamma in Gamma, the bijective normalization
nu_D in (1) and the isomorphism h give

    # {d in D : pi(d)=gamma}
       = # nu_G^-1(gamma)
       = number of analytic branches of Gamma at gamma. (8)

All branches of the target curve are represented, with one actual
boundary point for each branch. This is a statement about this finite
birational curve map; it is not an identification of a branch count
with a ramification index.

At d=nu_D(t0), the differential of pi:U->Y is an isomorphism.
Since beta=pi composed with nu_D,

    d beta_(t0)=0 if and only if d nu_D_(t0)=0.           (9)

For an analytically unibranch reduced curve in a smooth surface,
the normalization differential is nonzero precisely at smooth points
of the curve. In a normalization parameter, nonzero differential means
some local function has order one; the completed local ring then contains
C[[t]] and is equal to it. Conversely a smooth curve has an isomorphic
local normalization. Thus

    beta is an immersion everywhere
       if and only if D is a smooth algebraic A1.        (10)

The ambient etale map identifies the germ of D at d with the corresponding
ONE analytic branch of Gamma, not necessarily with its whole multibranch
germ. Equations (8)--(10) therefore retain possible cuspidal branches
until smoothness of D is separately established.

If D is smooth, every branch of Gamma is smooth. Any singular target
point must then have at least two distinct normalization preimages.
The branches may have arbitrary tangency order; ordinary nodes are not
forced. The conditional whole-root-chart theorem proves D smooth when
an individual regular source coordinate exists. This residue theorem concerns the finite map of curves.
[The local finite-degree argument](../../jacobian/monodromy-theory/KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#local-degree-and-the-finite-branch-quotient)
now proves smoothness of every unramified omitted prime and ensures
that its target normalization is immersive. This uses Orevkov's
finite-boundary resolution and purity. The existence of a common
source coordinate is a separate question.

## 5. The actual target image is singular

The nonproperness set of K equals pi(Z minus X). For a point outside
that finite closed image, the entire finite pi fibre is in X; over that
open subset K is a restriction of the finite map pi and is proper.
Conversely, if z lies in Z minus X, choose source points tending to z.
They escape every compact subset of the original A2, since a bounded
subsequence would have a source limit, contradicting the open embedding
into the separated Z. Their target images tend to pi(z). This proves
the asserted equality with the nonproperness set.

The omitted boundary of the normal affine Z is pure of dimension one:
the same normal-Hartogs argument as in Section2 of the topology source
rules out isolated components because X is affine. Finite morphisms
preserve the dimension of these curve components. Hence Gamma is an
irreducible component of the actual nonproperness set.

The separately retained EXT-CHAU-COMPONENT states that a nonsingular
polynomial map of C2 has no irreducible exceptional-value component
isomorphic to A1. For Keller maps the critical-value set is empty, so
it applies to nonproperness components.

The exact primary source is Nguyen Van Chau, *Two remarks on non-zero
constant Jacobian polynomial maps of C2*, cached at
[printed pp41--43](../../jacobian/affine-complements/sources/global-boundary-topology/chau-2003-two-remarks.pdf).
We inspected Lemma5, Lemma6(ii) and the component Remark on printed p43.
The proof uses the polynomial dicritical parametrizations and excludes
a line component after the Abhyankar--Moh coordinate change; these
classical inputs remain adopted at their retained scope. It is the
component Remark, not the stronger-looking Corollary3 about an entire
exceptional curve, that is used.

If Gamma were smooth, its normalization in (6) would be an isomorphism
and Gamma would be A1, contradicting that theorem. Therefore Gamma is
singular. This last singularity corollary uses the explicit Chau input;
the residue-degree and normalization theorem in Sections1--4 does not.

## 6. Exact controls and limits

The cusp t->(t^2,t^3) is finite birational with normalization A1 but has
zero derivative at zero. With the identity ambient etale map, it shows
that the curve-normalization argument alone does not prove smoothness
or immersion. This comparison is not a marked-plane completion.

For k>=1 the polynomial normalization

    p=t^2-1, q=t*(t^2-1)^k

is an immersion. It identifies t=1 and t=-1, and has branch contact
order k there. It is finite because t^2=p+1, and is birational because
t=q/p^k. Thus even residue degree one PLUS immersion does not imply
an embedded affine line or ordinary double points.

A concrete etale ambient comparison for k=1 is

    pi(t,s)=(t^2-1+s, t*(t^2-1)+(3/2)*t*s),
    U={1-(3/2)*s !=0}, D={s=0}.

Its target image is the preceding node. In the fibre product with its
normalization parameter r, eliminating s gives

    (t-r)*(t^2+r*t-2r^2+2)=0.

The selected component t=r is the whole finite A1. The other component
loses its branch points on the indicated etale open and is not finite
over the r-line. This guards the choice of component in Section2.
Its nonconstant ambient Jacobian means this is not an actual whole-plane
Keller pair.

Finally t->t^m, m>1, has residue degree m but ramifies at zero.
Removing zero gives the familiar degree-m finite etale cover C*->C*.
These controls show why both etaleness and the whole A1 normalization
are retained in the Euler argument.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_unramified_boundary_residue_degree.py

The checker verifies exact normalization parametrizations, ambient etale
and disconnected-fibre-product controls, arbitrary-contact polynomial
immersions, curve differential orders and the finite-cover Euler equation.
It does not replace the unrestricted finite/birational/etale argument.
All original finite-normalization primes stay in scope. No common source
coordinate, target-coordinate embedding, generic-degree bound, trace
perfection or HC4/JC2 inversion is asserted. Internal replay is not
formal verification or external human review.

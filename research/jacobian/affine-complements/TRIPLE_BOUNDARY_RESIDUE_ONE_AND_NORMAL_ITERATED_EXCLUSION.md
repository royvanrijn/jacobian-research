# Triple-boundary residue one and the normal iterated exclusion

PF2RULED105, 16 November 2026. All fields and germs are complex. Retain
the ACTUAL finite normalization N of a plane Keller map and its WHOLE
marked source U=A2. No smoothness of N is assumed.

**All-degree theorem.** Every omitted prime with transverse index e<=3
has residue degree f=1. At a singular point on an e=3 omitted prime,
the surface singularity is A1, the map to the normalization of the image
is locally unramified, and for some integer a>=2 the local surface degree
is 2a and the chosen image branch has the single Puiseux pair (a,2a-1).

Consequently EVERY omitted residue degree is one whenever the generic
field degree d<=8. This is a residue theorem, not an automorphism theorem
for all degrees through eight. A nontrivial omitted residue cover in
ANY degree must have e>=4 and must ramify at the cyclic quotient point
of its boundary prime.

In particular the natural eight-sheet colored iterated family of
[PF104](FIBRE_TRANSITIVE_SHEAR_RECURSION_AND_SPIN_BUDGET_EXCLUSION.md)
is excluded for arbitrary NORMAL actual N: its ramified peripheral orbit
(e,f,a_shear)=(3,2,0) is impossible. This removes the smoothness hypothesis,
at every finite cable length and all the independent parameters in PF104,
with any extra omitted primes and arbitrary finite image singularities.
The conclusion applies to ANY actual marked infinity cover with an
omitted index-three longitude-joining orbit f>1, regardless of the knot
type, field degree, end homology or the specific colored family.

## 1. The primary boundary input and its attribution

Use Borisov, *On the Stein factorization of resolutions of two-dimensional
Keller maps*, Beitr. Algebra Geom.56(2015),299-312,
[DOI10.1007/s13366-014-0208-4](https://link.springer.com/article/10.1007/s13366-014-0208-4).
The cached author manuscript has a different title. Its complete main
Theorem(4)--(6) on printed p.3 and Theorems3.4--3.5 with their proof
passages on pp.11--13 were read, rendered and visually inspected.
The exact pins and scopes are in the
[manifest](sources/triple-residue-excess/manifest.json).

Besides the disjoint unibranch A1 boundary and cyclic quotient surface
points already imported in [PF2BATT1](../jc2-boundary-attachment/PROOF.md),
we use the following stronger part of the SAME primary theorem:

    If D has no finite cyclic quotient point, D^nu -> C_D^nu is an
    isomorphism. Otherwise, after affine coordinates, the normalized
    map is t -> t^f. Its possible finite ramification point is the
    point over the cyclic quotient point A_D; its local degree is f. (1)

The location in (1), not merely the abstract power-map form, matters.
It is the main Theorem(6), with A_D from (4)--(5), and the proof of
Theorem3.5. On the projective boundary, the other ramification point
is the common point over target infinity, which is removed in N.
The finite Stein normalization is the original N by uniqueness, so no
target blowup is silently substituted. In a relatively minimal source
model, nonempty affine endpoint chains contract to nonsmooth points;
Orevkov explicitly makes this reduction in the proof on p.8.

This normalized power-map input was ALREADY used in the HC4
[source-module layers note](../../hessian/hc4-triangular-graphs/SOURCE_MODULE_LAYERS_AND_TRACE_ZERO_KERNEL.md).
It is not a new geometric discovery here. In particular the smooth
f>1 exclusion, including PF104's smooth natural8 consequence, already
follows from (1). PF104's source lattice and full-shear transfer proofs
remain independent results; the smooth exclusion is not claimed as
a new consequence beyond the published result. The new local deduction below
handles the nonsmooth index-three point and removes that remaining
source hypothesis.

The e1 and e2 conclusions were already established in PF2BATT1 and
[PF2OMIT7](../jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md).
Their predecessor proofs, receipts and scopes are preserved.

## 2. Corank one gives locally primitive residue at EVERY index

Let h:(C2,0)->(C2,0) be finite, with one smooth reduced ramification
branch A of generic index e>=2. Suppose dh has rank one. Choose the
first source coordinate u to be a target coordinate, so h=(u,G(u,w)).
If A is transverse to the kernel, u|A is a local parameter and the
map A -> normalization(h(A)) has local degree one.

Otherwise A is tangent to the kernel, hence A={u=phi(w)} with
phi(w)=c*w^m+higher terms, c!=0, m>=2. Finiteness rules out phi=0
identically. The Jacobian divisor gives

    G_w=V(u,w)*(u-phi(w))^(e-1), V(0,0)!=0.

Subtract the analytic target function G(u,0) from G. Restricting its
integral in w to u=phi(w) gives first orders

    ord(u|A)=m,
    ord(G|A)=m*(e-1)+1.                             (2)

The second leading coefficient is

    V(0,0)*c^(e-1)*Integral_0^1 (1-s^m)^(e-1) ds,
    Integral = product_(j=1)^(e-1) (m*j)/(m*j+1)>0.  (3)

Higher analytic terms have strictly higher weighted order, so cannot
cancel (3). If the local normalized residue degree were q, q would
divide both orders in (2), whose gcd is one. Thus q=1. Also the local
surface degree is ord G(0,w)=m*(e-1)+1. This is an all-index statement,
not an extrapolation from folds or a bounded jet calculation.

## 3. A corank-two germ with Jacobian twice a smooth branch

Now let h=(u,v) be finite on C2, dh_0=0, and

    Jac(u,v)=epsilon*x^2, epsilon(0)!=0, A={x=0}.    (4)

Both components have quadratic initial forms U,V, and their homogeneous
Jacobian is a nonzero multiple of x^2. Write

    U=a*x^2+b*x*y+c*y^2, V=d*x^2+e*x*y+f*y^2.

Vanishing of the xy,y^2 Jacobian coefficients gives af=cd and bf=ce,
while ae-bd!=0. If c!=0, these equations force ae-bd=0; hence c=0.
Then f=0 as well. Thus U,V are x times two independent linear forms.
Target linear combinations put them in the form x^2,xy.

Choose the second target coordinate GENERICALLY. Its quadratic initial
form is nodal, and it has minimal order on A among target linear forms.
Analytic factorization at the two distinct tangent lines gives source
coordinates X,Y with

    v=X*Y, A={X=beta(Y)}, ord beta=b>=2.

Here beta is not zero identically: otherwise v|A=0, and (4) forces
u|A constant, contradicting finiteness. The coordinate change

    Y'=Y*(beta(Y)/Y^b)^(1/(b+1)), X'=X*Y/Y'

is analytic and invertible, preserves v, and makes A={X'=Y'^b}.
Relabel X',Y' as X,Y. Equation (4) becomes

    (X*d_X-Y*d_Y)u=unit*(X-Y^b)^2.                  (5)

Subtract a convergent target function A0(v) to remove all diagonal
monomials X^iY^i in u. The diagonal series converges by the Cauchy
coefficient bounds for the original analytic u. With weights(b,1),
the right side of (5) has first weight 2b and initial form
epsilon0*(X-Y^b)^2. The operator multiplies X^iY^j by i-j. Therefore
all nondiagonal weights below 2b vanish, and the first weighted part is

    epsilon0*(X^2/2+2*X*Y^b/(b-1)-Y^(2b)/(2b)).     (6)

On A this has nonzero coefficient

    epsilon0*(b+1)^2/(2*b*(b-1)).                   (7)

Thus in these actual analytic coordinates the two restricted target
functions have orders b+1 and 2b. Their normalized local residue degree
q_h must divide

    gcd(b+1,2b)=gcd(b+1,2)<=2.                     (8)

This is the critical bound. It does NOT assert q_h=1: the genuine
degree-eight finite polynomial countercontrol of PF104 has q_h=2.

For later degree information, set m=ord(v|A)=b+1. Since v=XY, its two
branches are the axes. Equation (6) gives orders 2 and 2b for u on
those axes. Additivity of intersection multiplicity yields

    deg_0 h=I(u,v)=2+2b=2m.                         (9)

The quadratic-jet/two-branch intersection mechanism is retained with
credit from the cached
[endpoint note](../../hessian/hc4-triangular-graphs/sources/alok-one-dicritical-source-smoothness.md).
The coefficient argument and weighted normalization above supply the
additional LOCAL RESIDUE bound (8); no broader one-dicritical conclusion
from that exploratory source is assumed.

## 4. Pull back a singular index-three point to its uniformization

Let p be a singular point on an actual omitted prime D with e=3.
Use the precise fixed-affine-target endpoint input of PF2OMIT7:

    pi1(link(N,p) minus link(D,p))=Z.                (10)

The imported computation is Orevkov Lemma2.1 and the end-chain group
presentation on p.9. The presentation deletes the relation at the
endpoint curvette and expresses every meridian in one generator, with
no final relation. THIS computation is independent of ell_p=e_D;
the final image-smoothness conclusion of Lemma5.2 is NOT imported
without its separate equality hypothesis. The full relevant primary
pages2,7,8,9 were freshly rendered and read.

The cyclic quotient uniformization rho:(C2,0)->(N,p) has degree n>1
and is quasi-etale. The reduced inverse image A of D is ONE smooth
branch. Indeed its complement covers the complement in (10), so its
fundamental group is a finite-index subgroup of Z, hence Z. Its first
homology gives one branch; the plane-branch knot criterion then makes
it smooth. This is the explicit imported topology already audited in
PF2OMIT7, not an assertion about arbitrary cyclic quotient boundary
curves. The quotient action is free away from zero, so A->D^nu has
local degree n.

For h=pi composed with rho, the different pulls back in codimension
one to twice A. Thus Jac(h)=unit*x^2 after straightening A. There is
no extra isolated zero of a holomorphic function. The two target
components are invariant under the nontrivial SMALL cyclic action;
there are no invariant linear forms, so dh_0=0. Section3 applies.

If q is the local degree D^nu->C_D^nu at p, composition of these curve
maps gives q_h=n*q. Bound (8) therefore forces

    n*q<=2, n>=2, q>=1  ==>  n=2 and q=1.           (11)

The small order-two action is analytically -I; its quotient is the A1
surface ac=b^2. This also recovers the order-two conclusion of the
earlier endpoint quadratic-character argument.

Here q_h=2, so m=b+1 is even. Put m=2a, a>=2. The normalized image
parameter divides the orders in (6)--(8) by two, giving the coprime
first orders (a,2a-1). Thus there is just that one Puiseux pair.
Local surface degrees multiply in the finite composition; (9) gives
2*ell_p=deg_0 h=2m, hence ell_p=m=2a. This proves all singular-point
claims without assuming that N was smooth.

## 5. Global residue one and complete normal consequences

For e=3, if f>1, the primary power-map/location statement (1) requires
the local degree at the finite cyclic quotient point to be f>1.
But (11) proves it is one. Contradiction. Hence f=1 for every actual
index-three prime at EVERY field degree. Combining with the retained
e1/e2 results proves the all-degree e<=3 theorem.

The full residue Euler identity, with every prime and finite excess
retained, is

    d-1=sum_D e_D*f_D+sum_(D,p) eta_p, eta_p>=0.     (12)

If d<=8 and f>1, the theorem forces e>=4, so a single prime already
costs e*f>=8>d-1, impossible. Thus every omitted f is one in degrees
through eight, extending PF2OMIT7's through-six residue statement.
The analogous consequence with ALL indices odd is f=1 through degree10,
since the first possible nontrivial residue cover would have e>=5 and
cost at least10. Neither statement excludes all maps in those degrees.

Under the actual marked radial comparison of PF98, each omitted A1
has one end circle, and its joint peripheral lattice has longitudinal
degree exactly f_D. Hence no omitted index-three prime can give an
orbit joining more than one meridian cycle, in ANY degree or knot type.
Normal shear remains an independent integer and is not forced to zero.
Ordinary, non-omitted preimage curves can still have f>1; they must not
be added to the omitted ledger.

PF104's natural8 family retains one ramified orbit (3,2,0) at every
finite length and all of its independent coprime1mod6 cable parameters.
It would have to be an omitted prime in any Keller realization. Its
f=2 now contradicts the all-degree theorem even for SINGULAR normal N.
Thus that entire colored family is closed without smoothness or any
end-homology hypothesis. Other representations of the same knots are
not automatically excluded. PF104's regular24 exclusion also persists.

## 6. An unbounded genuine finite counterfamily above index three

The local conclusion does not extend automatically to all indices.
For any integer s>=1 put p=2s, m=4s, b=4s-1, and define on C2_(X,Y)

    V=X*Y,
    U=sum_(j=0)^p [(-1)^j binom(p,j)/(p-m*j)]
                         *X^(p-j)*Y^(b*j).          (13)

No denominator is zero, since p-m*j=2s*(1-2j). Direct monomial
differentiation gives

    Jac(U,V)=(X-Y^b)^p.                             (14)

Both variables are integral over C[U,V]. For Y, multiply (13) by Y^p
and replace XY by V; the resulting polynomial has nonzero constant
leading coefficient at Y^(m*p). For X multiply by X^(b*p) instead;
the same argument gives a monic equation after constant rescaling.
The only point of the fibre U=V=0 is zero; the two axes of V=0 give
local degree p+b*p=m*p=8s^2. Thus this is an ACTUAL finite polynomial
map of that generic degree, not merely formal data.

All monomials in U and V have even total degree. They descend through
the free-away-from-zero quotient by (X,Y)->(-X,-Y), giving a finite
map from an A1 surface N_s to A2 of generic degree4s^2. Its smooth
critical branch A={X=Y^b} descends to one D=A1 through the singular
point. In the coordinate t=Y^2 on D,

    V|D=t^(2s), U|D=C_s*t^(s*(4s-1)),
    C_s=2^p*(p-1)!/(2p-1)!! !=0.                   (15)

For the coefficient, the partial-fraction identity

    sum_(j=0)^p (-1)^j binom(p,j)/(z+j)
                     = p!/[z(z+1)...(z+p)]

at z=-1/2 proves (15) for ALL s. Thus the generic transverse index is
e=2s+1 and the residue degree downstairs is f=s. At s=1 this is the
allowed singular e3,f1 germ of degree4; at s>=2 it gives arbitrarily
large residue degrees at higher odd indices, starting at e5,f2.

These normal finite maps have the power-map form and a single cyclic
quotient point, yet are NOT Keller-plane realizations. In coordinates
g=X-Y^b,Y the involution is still -I. Removing D gives the quotient
of {g!=0}, with coordinates g^2 invertible and g*Y, hence

    N_s minus D = C* times C, NOT A2.               (16)

Their global length ledger has d=e*f+eta, differing by exactly one
from (12), because the marked open has Euler characteristic zero.
They preserve both the singular e3 possibility and the essential limit
of any attempted all-index residue-one extension. The degree-four
instance is an abstract normal finite cover, not a reopening of the
published Domrina quartic Keller exclusion.

## 7. Verification and remaining obligation

[The checker](verify_triple_boundary_residue.py) verifies symbolic
weighted coefficient identities and gcd arguments, explicit rank-one
and rank-zero germs (including primitive higher-characteristic and
nonprimitive controls), quotient degrees, the genuine finite family,
global degree-budget endpoints, source pins and predecessor seals.
The [receipt](triple_boundary_residue_receipt.json) distinguishes these
finite controls from the written all-index/parameter proofs and the
published boundary inputs. No finite enumeration proves the theorem.

The all-degree global Wright objective remains OPEN: f=1 packets,
e>=4 cyclic quotient residue covers, other colored infinity actions,
and the actual Laurent-terminal finite-algebra implication remain.
The new result is not all-degree JC2, all eight-sheet exclusion, or
automatic absence of singular triple points. External review stays
deferred; formal and independent external verification remain false.

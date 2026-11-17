# The residual curve of a finite quartic Wright map

14 November 2026. PF2RULED33 gives a uniform global reduction, with
no pole-height or coefficient-degree cutoff. The normal coefficient
identity and the residual discriminant calculation hold in every
possible finite degree. For degree four they determine the entire
normalization of the residual curve from the boundary collision data.

Suppose there are n collision values, and the two smooth image branches
at the i-th value have intersection multiplicity m_i. Let o count the
odd m_i. Then the normalization of the residual source curve is:

* connected exactly when o>0;
* in the connected case, an open curve of compact genus o-1 with
  4n-2o+3 punctures;
* in the all-even case, two punctured rational curves.

For ordinary nodes, o=n: the genus is n-1 and there are 2n+3 punctures.
This constrains every pole height at once. It does not prove that the
required curve cannot occur on the actual Wright surface.

## Actual hypotheses and notation

Keep [PF2RULED32](GLOBAL_BOUNDARY_OBSTRUCTION.md). A finite map
pi=(P,Q):N->A2 comes from P,Q in the actual marked Wright ring, with
nonzero constant source Jacobian. Its source open U=N\R is the specified
polynomial plane. The omitted R is A1_t, and its finite birational image
normalization is gamma(t)=(p(t),q(t)). Choose a reduced irreducible
image equation H(p,q)=0. Write

    h=gcd(p',q'), a=p'/h, b=q'/h,
    p1=lambda*a, q1=lambda*b,
    W=a*b'-a'*b,
    2h*(p2*b-a*q2)+lambda^2*W=kappa != 0.

These are the actual all-height chart identities in
[PF2OMIT9](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md).
The boundary budget gives d>=3, and a compulsory collision already
costs at least four sheets. Thus every possible map here has d>=4.
In the omitted chart the source coordinates are X=1/u,
Y=-u^2+u^3*t. Define

    beta(t)=[u^2] H(P(u,t),Q(u,t)).

Since (a,b)=1, differentiating H(gamma)=0 gives a unique polynomial
rho such that

    (H_p,H_q)(gamma)=rho*(-b,a).

The residual divisor T is defined by

    div_N(H(P,Q))=2R+T.

It is reduced at its generic points and has generic degree d-2 over
the image curve. Let G=F^(-1)(C) be the plane source curve. Then
G=T intersect U; normalization commutes with this open restriction.

## A normal-coefficient identity in every finite degree

The Hessian chain rule along gamma gives

    Hess_H((a,b),(a,b))=-rho*W/h.

The second normal coefficient is consequently

    beta=-rho*(p2*b-a*q2)-rho*lambda^2*W/(2h)
         =-kappa*rho/(2h).                            (1)

This is an identity of rational functions, and beta itself is a
polynomial. Therefore h divides rho. No parameter value with h=0 has
been removed: the identity extends there polynomially. Higher normal
coefficients and all pole spaces are unrestricted in this calculation.

Let E_res be the generic etale algebra of the d-2 unramified sheets
above C, over k=C(t). Its discriminant squareclass is

    disc(E_res) = beta in k^*/k^*2.                   (2)

Here discriminant means the orientation quadratic algebra, including
the split case. The standard orientation interpretation is
Biesel--Gioia, *A New Discriminant Algebra Construction*, Documenta
Math.21(2016), Theorem6.8, with the characteristic-not-two description
in its introduction and Section6.
[Primary publisher PDF](https://ems.press/content/serial-article-files/26337?nt=1).
This standard interpretation is imported and pinned; the specialization
argument follows below.

The global sign cover has squareclass H: the only finite branch prime
is C and its generic monodromy is a transposition. Factorization in
the target polynomial UFD removes all even valuations, leaving H up
to a square and a nonzero complex constant. On the source,
H(P,Q)=u^2*(beta+O(u)); after dividing its square root by u, the sign
cover is unramified at generic R and has residue class beta.

Its residue orientation is exactly the orientation of E_res. Indeed
the generic inertia swaps the two ramified sheets. Its centralizer
preserves that pair and the other d-2 sheets. The decomposition group
contains this swap; fixing one ramified sheet removes the pair swap,
so the remaining sign is precisely the sign on the unramified sheets.
This proves (2). Equivalently, a tame quadratic block over the target
DVR has trace discriminant 4H/beta; multiplication with the unramified
block's discriminant gives the same squareclass identity.

The residual algebra lives over the boundary function field C(t).
It is not a degree-two intermediate field of the non-Galois surface
extension C(X,Y)/C(P,Q). Field primitivity does not exclude it.

## Complete quartic collision data

Now let d=4. There is one clean ordinary cusp parameter s, with h=t-s.
By PF2RULED32, n>=1. Each collision value has exactly two distinct
immersed parameters a_i,b_i. All these parameters are distinct and
different from s. The intersection multiplicity m_i is any positive
integer; no nodality assumption is made.

At the cusp rho has order three and h has order one. At either
parameter of a two-smooth-branch collision, rho has order m_i and h
is a unit. At smooth image points rho is a unit. Thus (1) gives the
complete finite zero divisor

    beta=c*(t-s)^2*product_i ((t-a_i)*(t-b_i))^m_i,
    c != 0.                                         (3)

These orders follow directly from the plane-curve gradient. At a
collision, the other local branch contributes its intersection order
m_i to the gradient on the chosen branch. For the ordinary cusp
(t^2,t^3), the gradient divided by the primitive tangent has order
three. Analytic coordinate and equation units preserve these orders.

Put

    S(t)=product_(m_i odd) (t-a_i)*(t-b_i).

All constants are squares over C. The cusp factor and the even
collision factors in (3) disappear from the squareclass. Since a
rank-two etale algebra is its own orientation algebra, the complete
normalization T^nu is

    w^2=S(t) over A1_t.                              (4)

More explicitly, it is the normalization in that quadratic algebra.
The displayed model is already smooth when S is nonconstant, because
S is squarefree. If S=1 it is the disjoint union of two affine lines.
It is finite over A1: normalization of the finite residual curve maps
finitely to the normalization of C, and a normal finite cover is
uniquely determined by its function algebra. Hence no extra finite
points are omitted from (4).

An independent local check of (4) is useful. At a fold, use
(z,y)->(z^2,y). Pulling back another smooth branch x=f(y), ord f=m,
gives the residual equation z^2=f(y). Its normalization has one point
of index two when m is odd, and two unramified points when m is even.
The same occurs at both parameters of each collision.

## Which points are missing from the plane source

At every collision value the original source fibre is empty. Thus
all points of (4) over a_i and b_i lie on R and must be deleted from
T^nu to obtain G^nu.

Over the cusp parameter there are two unramified points of T^nu.
Exactly one lies on R and is deleted; the other normalizes the unique
retained source cusp. This can be checked in the complete cusp block

    P=u-3v^2, Q=uv-2v^3,
    4P^3+27Q^2=u^2*(4u-9v^2).

The residual curve is parametrized by u=9t^2, v=-2t and maps to the
boundary normalization (-3t^2,-2t^3) with index one. The fourth, etale
sheet supplies the other point. All other finite parameters have two
retained source points. This identifies every finite deletion.

If o>0, S has degree 2o. The compact double cover has exactly those
2o simple branch points, with two unramified points at infinity.
The [Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
therefore gives genus o-1. This standard projective formula is imported
with its smoothness, separability and properness hypotheses; it is
applied to the smooth projective completion, not directly to G.

Finite deletions comprise 2o points over odd collision parameters,
4(n-o) over even ones, and one cusp point. Adding the two points at
infinity gives

    punctures(G^nu)=4n-2o+3,
    chi(G^nu)=2-2(o-1)-(4n-2o+3)=1-4n.               (5)

If o=0, each component is rational. Delete all 2n collision parameters
on both copies of A1 and the cusp point on one copy. Their compact
puncture counts are 2n+2 and 2n+1, with the same total Euler
characteristic 1-4n. In particular G is irreducible if and only if
at least one collision multiplicity is odd.

The only finite singularity of G is its retained ordinary cusp:
F is etale, and the other singular image values are omitted. That
cusp is unibranch, so passing from G^nu to G does not change the
Euler characteristic in (5).

## Coupling the abstract cover to the actual Wright surface

The source coordinate X has exactly the following finite poles on
T^nu:

| Deleted point | Pole order of X |
| --- | ---: |
| Deleted cusp point | 2 |
| Point over an odd collision parameter | m_i |
| Each of the two points over an even collision parameter | m_i/2 |

At an odd collision z=local_parameter^m_i and the base parameter has
order two. At an even collision each branch has z of order m_i/2.
The marked R-uniformizer differs from z by a unit. At the cusp the
explicit block gives order two. Since X=1/u, these are exact pole
orders. Y has twice the corresponding zero order.

The regular Wright function D=X^3Y+X takes value b_i at every deleted
point over a_i, and value a_i over b_i. At the deleted cusp point its
value is s. These opposite-parameter conditions retain the actual
source marking; an abstract double cover alone does not supply them.
The behavior of the source functions at infinity remains part of
the embedding problem.

There is also a global divisor-class constraint. The source volume
form has divisor R, so K_N is represented by R. The definition of T
gives

    T linearly equivalent to -2K_N.                  (6)

Thus the required curve must occur in this specific divisor class
on the fixed Wright surface, with all the poles, deletions and marked
values above. An arbitrary plane embedding of a punctured curve is
not sufficient.

For example, one nodal collision gives the conic w^2=t^2-1, whose
compactification is rational. The substitutions t=(z+z^-1)/2,
w=(z-z^-1)/2 identify it with G_m. Deleting its two branch points and
one cusp point gives a rational curve with five punctures, as (5)
predicts. This is a consistent abstract passport, not a finite Wright
map. It shows why the residual reduction is a global constraint,
not an exclusion by itself.

## The associated sign cover has trivial canonical sheaf

There is a second global consequence of (6). Let W be the normal
target sign double plane z^2=H, and let Z be the normalized pullback
of that sign cover to N. The map Z->N is a double cover branched
exactly along T: the twice-occurring R drops from the branch divisor.
Since T is linearly equivalent to -2K_N, this is the canonical
antibicanonical double-cover construction. Its canonical sheaf is
trivial. Equivalently, the source volume form divided by z pulls
back to a nowhere-vanishing canonical form: its R zero cancels the
R factor in z, and its pullback zero at the double-cover branch
cancels the remaining simple z factor.

In degree four, Z has one A2 singularity over the retained source
cusp and two A_(m_i-1) singularities for each collision, with A0
meaning a smooth point. These are the only finite singularities.
At a collision the local equation is

    eta^2=z^2-y^m, hence (z-eta)*(z+eta)=y^m.

At the deleted cusp point T is smooth, so Z is smooth there. At the
retained cusp the double cover of an ordinary cusp gives A2. All
these singularities are rational double points; the minimal resolution
has a nowhere-vanishing canonical form as well.

In every finite degree the map Z->W has degree d and is unramified
in codimension one. The sign extension is not contained in the
primitive degree-d field: its point stabilizer S_(d-1) contains odd
permutations. Their intersection is the target field. In the quartic
case the map therefore has degree four. Locally,
the sign base change kills the generic index-two ramification, which
proves the codimension-one assertion. No etale assertion is made over
the singular points of W. In fact W has A2 at its cusp and
A_(2m_i-1) at a collision, so such points matter.

These surface invariants are independent of pole height:

    chi(T)=2-2n,
    chi(Z)=2*chi(N)-chi(T)=2+2n,
    chi(minimal_resolution(Z))=4+2*sum_i m_i.

For the first equality, normalization identifies two pairs at each
even collision and no pairs at the unibranch singularities. Also
chi(N)=2 since N is the disjoint union of its plane open and R.
Resolving A_j replaces a point by a chain with Euler characteristic
j+1, adding j. This gives the last equality.

The remaining cover must therefore realize both the residual curve
passport and this specific canonical-trivial, finite sign-cover tower.
These numerical and local data are not by themselves a contradiction.
Canonical triviality here concerns an affine normal surface and its
quartic DuVal resolution; no projective K3 compactification is supplied.

## Next global gate and assurance

The next gate is an obstruction to embedding this residual curve
in the class -2K_N with the actual source-coordinate pole and value
data, or an obstruction to the associated degree-four sign-cover extension
with these canonical and rational-double-point data.
The discriminant reduction replaces the entire residual rank-two
cover with (4); it leaves the global embedding and finite algebra
conditions to be proved. No coefficient or tangency order is bounded.

Run the [checker](verify_residual_cover.py):

```sh
.venv/bin/python research/jacobian/affine-complements/verify_residual_cover.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The [receipt](residual_cover_receipt.json) pins the primary
[orientation and Hurwitz sources](sources/residual-cover/manifest.json)
and predecessor evidence. It checks the universal normal-coefficient
identity, trace discriminants, sign controls, complete cusp block,
all-order local parity models and symbolic genus/puncture formulas.
The unbounded scope and global normalizations are written proofs,
not conclusions inferred from finitely many controls.

No actual plane embedding or finite Keller map is constructed, no
colliding-image family is excluded in full, and no arbitrary-normalization
compatible completion is proved. The full goal remains active, without
formal verification, independent whole-proof implementation, external
review or novelty-priority claims.

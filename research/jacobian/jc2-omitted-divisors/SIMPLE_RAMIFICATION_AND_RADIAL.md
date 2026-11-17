# Simple boundary ramification and complete residual singleton radial fibres

14 November 2026. PF2OMIT7 refines [PF2OMIT6](RESIDUE_LEDGER.md).
It proves that a generic index-two omitted prime lies entirely in the
smooth locus of the original affine finite normalization and has residue
degree one. In the residual r7 excess-one row, the remaining exceptional
point is therefore a smooth-source ordinary cusp with fibre algebra
C[w]/(w^3). All three residual rows now have smooth affine N.

The same local calculation gives every omitted radial prime, its e,f,
its image and its fibre length, over a singleton zero in those rows.
The actual global affine image equations remain unknown. No unrestricted
singleton-zero properness or JC2 conclusion is made.

## 1. The boundary end-curve complement has cyclic fundamental group

Use the original affine target A2, its finite normalization N, and the
Keller source open U=A2. Retain the boundary structure in PF2BATT1.
In particular all omitted primes are disjoint and unibranch; every
possible finite surface singularity on one is a cyclic quotient.

At a point q on an omitted prime D, the end-chain computation in
[Orevkov, Lemma2.1 and the proof of Lemma5.2](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf)
gives

    pi_1(link(N,q) minus link(D,q)) = Z.                (S1)

This particular computation in the proof does not require the later
comparison ell_q=e_D. It uses the linear exceptional resolution chain,
with the strict boundary meeting its end. The relations express every
meridian in terms of one end meridian, with no remaining relation.
The finite Stein-model comparison is the one already written in
PF2OMIT6, with the affine target unchanged. The source is
[pinned](sources/residue_manifest.json); printed pages8--9 contain the
end-chain argument. This is an explicit imported topological input.

Suppose q is singular. Write its cyclic quotient uniformization

    rho:(C2,0)->(N,q)=(C2,0)/mu_r, r>1,
    (x,y)->(zeta*x,zeta^a*y), gcd(a,r)=1.

It is unramified away from the origin. Let A be the reduced inverse
image of D. The connected complement of this plane analytic curve
covers the complement in (S1), so its fundamental group is a finite-index
subgroup of Z, hence Z. The first homology shows that A is unibranch.
The plane-curve knot criterion then makes A smooth: a singular plane
branch has a nontrivial algebraic knot and a nonabelian complement group.
This is the same Milnor/Zariski knot input used in Orevkov's proof.

## 2. A singular boundary point cannot have generic index two

Let e=e_D. The two target coordinates pulled back by pi composed with
rho are invariant convergent functions P,Q on C2. Subtract their
constant terms. The nontrivial small cyclic action has no invariant
linear form, so P,Q belong to the square of the maximal ideal. Hence

    ord_0 Jac(P,Q)>=2.                                (S2)

On the other hand rho is unramified in codimension one, and D is the
only ramification germ of pi near q. Thus the Jacobian divisor of
pi composed with rho is exactly (e-1) times the smooth curve A. A
holomorphic unit cannot have an extra isolated zero, so

    ord_0 Jac(P,Q)=e-1.                               (S3)

Equations (S2),(S3) imply e>=3. In particular every e=2 boundary point
is smooth on N. The e=1 case was already etale by PF2BATT1. The proof
uses neither source finality nor any conclusion of Borisov Theorem3.8.
It is not a statement about an arbitrary finite normal surface cover:
the boundary end-chain condition (S1) is essential.

## 3. Simple ramification has locally primitive image parametrization

Now e_D=2. The surface and D are smooth at every finite point. The
Jacobian has a simple zero along D. The differential has rank one at
the point: rank zero would make every matrix entry vanish, forcing its
determinant into the square of the maximal ideal.

Choose source and target coordinates so that pi=(u,g(u,w)). At a fold,
D is transverse to the kernel and its image is immersive, with local
degree two. Otherwise g_ww(0)=0 and smoothness of the Jacobian divisor
forces g_uw(0)!=0. Write

    g_w=unit*(u-h(w)), ord h=m>=2.

The function h is not zero identically, since the map is finite.
After subtracting g_u(0)*u from the second target coordinate, the
restriction to D has

    ord_w u=m, ord_w g(h(w),w)=m+1,
    ell_q=ord_w g(0,w)=m+1.                           (S4)

Indeed g_u restricted to D has order one, and differentiating the
second restricted coordinate gives g_u*h'. The coprime first two
orders m,m+1 make the parametrization locally primitive. Consequently
the polynomial residue map D^nu=A1->C_D^nu=A1 is unramified at every
finite parameter, and therefore has degree one.

For an actual birational polynomial image parametrization (p(t),q(t)),
(S4) also gives the useful exact bookkeeping identity

    E_D=deg gcd(p',q'),
    ell_(t0)=2+ord_(t0) gcd(p',q').                    (S5)

These formulas require this global simple-ramification setup. For an
arbitrary parametrization they do not supply an ambient cover. For
example (t^2,t^5) has a first-degree derivative gcd but is not the cusp
of a length-three simple-boundary factor described here.

Together with PF2BATT1 and the residue budget of PF2OMIT6, this proves
that **every omitted residue degree is one for a plane Keller map of
geometric degree at most six**. Otherwise e>=3,f>=2 would cost at least
six in sum e*f, exceeding d-1<=5. This is only a residue statement,
not a small-degree automorphism theorem.

## 4. The residual excess-one point is an ordinary cusp

In residual row B, S has e=2 and E=1. It is now smooth on smooth N.
At its unique exceptional point ell=3, (S4) forces m=2. Thus the image
has an ordinary (2,3) cusp branch. Its image parametrization has a
unique simple common derivative zero t*, and

    det(gamma''(t*),gamma'''(t*)) != 0.                (S6)

The fibre algebra is C[w]/(w^3); the square-zero length-three algebra
left open in the preceding ledger cannot occur. A local representative
is pi=(u,w^3+u*w). All other omitted points in the three rows have the
cyclic representatives pi=(u,w^k), k=1,2,3, with their established
indices. Every affine N in the three residual rows is smooth.

At a singleton zero the complete original normalization fibre is:

| Residual singleton type | Omitted factors | Retained factor |
|---|---|---|
| A, five unramified points | C five times | C |
| A, one fold and three unramified points | C[w]/w^2, C three times | C |
| A, two folds and one unramified point | C[w]/w^2 twice, C | C |
| B | C[w]/w^3 at the cusp, C[w]/w^2 at a fold | C |
| C | C[w]/w^3 at the cyclic index-three point, C[w]/w^2 at a fold | C |

The rows B,C have the same Artin fibre algebras but different generic
ramification indices and different local branch maps. Those flags cannot
be recovered from the special fibre lengths alone.

## 5. The radial blowup ledger is explicit in every singleton type

Blow up the target at zero. For each omitted point, the radial ideal
in its regular local source ring is I=(u,w^k), with k its ambient
length. For the cusp (u,w^3+u*w) this is the same ideal as (u,w^3).
Its principalization has exceptional chain

    (-2,...,-2,-1), e=(1,2,...,k), f=(0,...,0,1).       (S7)

There are k components. The only radial horizontal is the last one;
all earlier components are contracted by the direction map. Directly
Q*e=-f and e^T*f=k. Equivalently the single compact Newton edge of
(u,w^k) has primitive normal (k,1), so its normalized blowup has one
radial prime, with (e_rad,f_rad)=(k,1). Its image is the whole target
exceptional P1. It is omitted; the actual source origin supplies a
separate retained radial prime (1,1).

The complete omitted radial lists are therefore

| Original singleton type | Omitted radial (e,f) | Total omitted degree |
|---|---|---:|
| A, five unramified | (1,1) five times | 5 |
| A, one fold | (2,1),(1,1),(1,1),(1,1) | 5 |
| A, two folds | (2,1),(2,1),(1,1) | 5 |
| B or C | (3,1),(2,1) | 5 |

No radial e=1 term is dropped. For the global completed source-volume
form, the principalization canonical vectors in (S7) are

    unramified k1: (1), fold k2: (2,3),
    cyclic k3: (3,4,5), cusp k3: (2,4,5).              (S8)

They follow from the local Jacobians 1,2w,3w^2,3w^2+u and the charts
u=a*w^j. Horizontals satisfy k_i=2e_i-1, verticals k_i>=2e_i,
as required by HC4MYRB1. Even the same radial (3,1) block retains a
different first vertical canonical order in the cusp and cyclic cases.
The strict old ramification curve meets the first exceptional component
in the fold/cyclic blocks, and the second component in the cusp block.
Including its coefficient e_old-1 makes the full canonical intersection
equalities agree with adjunction. Dropping that attachment would change
the canonical equations, as already warned in HC4MYRB1.

## 6. Special radial directions and the fixed-target scope

The normal finite pullback has explicit charts over the target blowup.
For a cyclic block v=w^k, its u-chart is w^k=u*tau, with tau=v/u;
its v-chart is v=w^k, u=sigma*v. The first chart has an A_(k-1)
singularity at tau=0 for k>=2; the second is smooth along the radial
prime. Every direction has exactly one point of this block, always
with ambient fibre algebra C[w]/w^k.

For the cusp v=w^3+u*w the two charts are

    w^3+u*w-u*tau=0,
    w^3+sigma*v*w-v=0.                                (S9)

In the first put z=tau-w to obtain u*z=w^3, an A2 singularity.
The second is smooth along the radial prime. Again the block has one
point of length three at every direction. Its projection is (u,tau)
with tau=z+w, rather than the cyclic projection (u,z).

At the special direction the strict old boundary and the new radial
boundary meet. Count their common point **once**, with length k;
summing their generic indices would double count it. The cubic cusp
point lies on indices2 and3 but still has length3. Different original
normalization points stay different blocks, even when tangent directions
coincide on the exceptional line.

These node surface germs occur **after a target blowup**. They do not
contradict PF2BATT1's prohibition on a collided-node projection in the
original N->A2. Its disjoint-boundary theorem cannot be transported to
the new finite normalization over Bl_0 A2 without a separate comparison.
Here the charts explicitly show the new intersecting radial/strict
boundary germs. The original prohibition and the local HC4MYRB2 node
calibrations retain their stated fixed-target scopes.

## Verification and remaining data

[verify_simple_radial.py](verify_simple_radial.py) checks invariant
linear/jet obstructions, the corank-one controls, cusp jets and fibre
algebras, complete principalization matrices and canonical orders,
normalized chart equations, and all singleton radial degree sums.
[simple_radial_ledger.json](simple_radial_ledger.json) keeps original
and radial indices, omission flags, images and special lengths separate.

The end-chain topology, cyclic quotient uniformization and plane-curve
knot criterion are explicit inputs, not consequences of bounded checks.
Actual global images C,B,H, their intersections and the exceptional
target value remain unconstructed. The research goal remains open.

# Nodal deformation and a nonpositive normal-degree obstruction

15 November2026. Claim PF2RULED59. Programme status: archived.
The full research goal remains active.

## 1. What is now excluded

For the ENTIRE surviving quartic Keller packet of PF2RULED46, let
h:P1->X be the normalization of its ramified image C, on the fixed
efficient target compactification with infinity boundary E. Then

    eta=deg N_h=-K_X.Cbar-2=a_n-1 <= 0.                    (1)

This holds for arbitrary finite tangencies and all source rulings and
pole orders. The positive-eta range is excluded outright. The older
PF2RULED46 bound a_n<=1+2T and PF2RULED47 inequality10*eta<=2*delta
remain true, but large contact cannot rescue ANY eta>0 quartic.

The argument is not intrinsically quartic. For a reduced polynomially
parametrized plane curve B, take a smooth projective completion X of A2
on which every image normalization h_i:P1->X is immersive and meets the
fixed boundary transversely once, away from its crossings and the other
strict curves. Set eta_i=deg N_(h_i). The following sufficient conditions
give an abelian Orevkov invariant:

* every eta_i>=2, with arbitrary numbers of branches at finite values;
* every eta_i>=1, and each INDIVIDUAL component has at most two
  normalization preimages over any value. Different components may meet.

Thus a nontrivial zero-excess Keller branch union as in PF2RULED46 must
have some eta_i<=1. If each branch-image component is at most double
valued, some eta_i<=0 is necessary. Use a simultaneous completion of the
WHOLE reduced branch union; separate infinity diagrams do not suffice.

In particular, with a sole index-two ramified prime, immersive residue-one
image and arbitrary unramified omitted primes, O(C) surjects onto S_d.
For ALL d this forces eta<=1; eta=1 requires a value with at least three
normalization preimages. When d<=5 the sheet budget forbids such a value,
so eta<=0. The quartic consequence (1) therefore closes a full contact
range rather than another bounded calculation. No result is asserted for
eta<=0, or for arbitrary non-zero-excess/nonimmersive branch packets.

## 2. Deform maps on the fixed compactification

For a birational immersion h:P1->X the exact tangent sequence is

    0 -> O(2) -> h*T_X -> N_h=O(eta) -> 0.                (2)

Dedieu--Sernesi, Section2.1, equations(2.1.2),(2.2.1)-(2.4), gives the
normal-sheaf deformation and obstruction spaces. The complete statements
and discussion are on author-PDF7-8. Here H1(O(2))=0 and H1(O(eta))=0
for eta>=1. Consequently H1(h*T_X)=0 and every normal section lifts to
an actual small deformation of h in the smooth morphism space. In
particular this is not a conjectural equigeneric smoothing of an
arbitrary divisor. We start with an immersion and check unobstructedness.

On P1 the evaluation map on k prescribed distinct points is surjective
for O(eta) whenever eta>=k-1: its cokernel is H1(O(eta-k)). For eta>=2
three-point separation, including jets, supplies the nodal-deformation
argument of Dedieu--Sernesi Theorem2.8.2 and its full proof on PDF10-11.
A branch can be moved off a multiple value while keeping two other
branches fixed; a branch can also be displaced transversely off another
branch's tangent. General nearby images have only ordinary nodes.

If the original component has only double values, eta>=1 suffices.
Immersivity is open on the compact source. The absence of three distinct
preimages over one value is also open near this immersion: any limiting
triple with distinct parameters would be an original triple, whereas
coalescing parameters are ruled out by uniform local injectivity of an
immersion on a compact source. Thus a sufficiently small family remains
at most double valued. Two-point evaluation now moves either branch
normally while fixing the other, breaking every tangency generically.
This is exactly the two-branch part of the cited proof; no three-point
separation is claimed for O(1).

For completeness, at a contact of order m use analytic coordinates
v=0 and v=u^m. A deformation with nonzero relative normal displacement
has branch difference u^m-epsilon+terms of higher deformation order.
Generic small epsilon splits the contact into m transverse nodes. The
discriminant in the deformation parameter is nonzero. There are only
finitely many contacts, so a general parameter avoids all their proper
bad loci simultaneously. Their multiplicities m are unrestricted.

For a union, deform the normalization maps independently. Under the
first sufficient condition, each component becomes nodal. Under the
second, it becomes nodal by the preceding two-point argument. Relative
positions of different components are general: normal sections have no
base point, so cross-component tangencies break. A triple involving two
branches of one component and one of another is avoided using that
component's two-point evaluation and the other component's independent
motion. A triple on three components is avoided by their independent
motions. Higher multiple values contain such triples. This gives a
nodal union under either sufficient condition.

## 3. Why this deformation preserves the whole affine end

The compactification X and its boundary E are held FIXED. Initially
Cbar_i.E_j is1 for the unique component met by h_i and0 for the others.
These intersection numbers do not change in a small family of maps.
No curve becomes a boundary component. Positivity of intersection and
openness of transversality imply that every nearby normalization still
meets E exactly once, transversely at a smooth point, with no intersection
with its other components. These intersection points remain separate
for different strict curves. The maps stay birational immersions near
the original maps, and the affine normalizations are still A1.

Near the compact boundary this is a smooth family of transverse pairs,
hence an ambient isotopy preserving E. Equivalently choose a regular
neighborhood of E and a common transverse boundary; the compact pairs
on that neighborhood are isotopic. This gives the SAME affine link at
infinity and a proper isotopy on the end of X minus E=A2. It avoids any
assumption that a perturbation of low polynomial coefficients must
preserve the infinity splice. That assumption would be unjustified.

Choose disjoint small balls about all original finite singular values.
The normalization disks deform analytically in these balls; their
number of analytic components is preserved. Outside smaller balls,
including the end just treated, there is a smooth proper isotopy.
Thus the family is local-component-preserving (LCP) in the precise
Neumann--Norbury Definition2, PDF2. It need not preserve singularity
types or the original complement fundamental group. Their Theorem1
and its complete proof on PDF4-5 give

    O(B) = O(B_t).                                        (3)

Only a small connected deformation germ is used. There is no appeal to
connectedness of the entire polynomial-parameter moduli space or to the
extra positivity hypothesis in their global moduli theorem.

## 4. Apply the nodal theorem after deforming

For every component of the nodal B_t, the normalization has genus0
and its class on X is unchanged. If r_i counts its SELF nodes, adjunction
gives

    Cbar_(i,t)^2-2*r_i = -K_X.Cbar_i-2 = eta_i > 0.        (4)

Mutual intersections of different components are not counted as
singularities of an individual component in this formula. The nodal
union meets E transversely. Nori Proposition3.27, printed331/PDF28,
therefore makes the kernel of

    pi1(X minus (Bbar_t union E)) -> pi1(X minus E)

abelian. The filled complement is A2, so the latter group is trivial.
The whole affine complement group of B_t is abelian. Its local nodal
groups are already abelian, hence O(B_t)=pi1(A2 minus B_t). Equation(3)
proves the sufficient abelianity statements in Section1.

For a zero-excess full Keller branch union, PF2RULED46 gives a surjection
O(B)->G onto a NONABELIAN finite monodromy group, for arbitrary indices
and numbers of ramified primes. Either sufficient condition is therefore
impossible. For a sole index-two ramified prime it gives O(C)->>S_d.
Every normalization preimage of C consumes two finite sheets, so a triple
value needs d>=6. The asserted all-degree and d<=5 bounds follow.
For the surviving quartic, only double values occur, so (1) follows.

## 5. Exact limits and a control against a stronger false conclusion

The result is about O(B), not the complement group of the original
tangential curve. Retain PF2RULED56's irreducible control

    G0=x^2-(x+y^2)^3,
    x(t)=(1-t^2)^3,       y(t)=t(1-t^2).                    (5)

It is immersive with sole contact3 at t=1,-1. In the polynomial target
coordinates u=x+y^2, y, its degrees are4,3, with coprime leading orders;
the infinity knot is T(3,4). Thus a_n=7 and eta=6>0. Nevertheless the
cubic w^3-3(x+y^2)w-2x has discriminant -108G0 and a triple fibre at the
origin with local S3 monodromy, as retained in PF56. Its global complement
monodromy contains S3. There is no claim that pi1(A2 minus G0=0) is
abelian. The local nonabelian image is killed in O(G0); this is why the
actual finite Keller map's factorization through O(C) is essential.

At eta=1, O(1) separates TWO points but not THREE. The unrestricted
multiple-value argument accordingly stops at eta>=2. We do not assert
that every eta=1 immersion with triple values has a nodal deformation.
At eta=0 even two-point evaluation fails. No conclusion is drawn there
or at negative eta. Nor does this theorem bound all target degrees,
exclude all iterated cablings, remove unramified primes, or apply to
the unramified image H as though H were a branch component of pi.

The old contact budgets are retained evidence. For the sole-branch
quartic, their eta>0 sectors are now closed by deformation. PF2RULED47's
selected-component centralizer bound remains useful when other branch
components fall outside the simultaneous nodal theorem. PF2RULED48's
unramified-sheet bound remains separate.

[The checker](verify_nodal_deformation.py) verifies the exact P1
evaluation controls, arbitrary-contact derivative identity, normal-degree
and nodal adjunction identities, and the negative control's polynomial
and discriminant data. The written proofs and imported deformation,
LCP and nodal theorems are not replaced by these finite checks.
[The manifest](sources/nodal-deformation/manifest.json) pins the fully
inspected primary pages and sealed predecessors. No full quartic or
higher-degree exclusion, formal verification, independent full reproof,
external review, novelty priority or research-goal completion is claimed.

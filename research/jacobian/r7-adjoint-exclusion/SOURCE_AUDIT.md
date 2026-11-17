# External geometry and the exact source-model bridge

14 November 2026. The source is Alexander Borisov, *On the Stein
factorization of resolutions of two-dimensional Keller maps*, Beitr.
Algebra Geom.56(2015),299--312, DOI10.1007/s13366-014-0208-4. The pinned
[author PDF](../jc2-boundary-attachment/sources/borisov-stein-author.pdf)
has a different title and is 15 pages. Page numbers below are its printed
page numbers, not a silent identification with the publisher pagination.

The full source was not independently reproved or formalized. We inspected
the relevant statements and proof passages. There is no communication
with the author and no external endorsement of the present deduction.

## Inputs retained from the earlier boundary audit

- Proposition3.3 and Theorem3.4, pp.11--13: the projective finite Stein
  surface W contains the whole affine source; all its boundary curves
  meet at the same point A over target infinity and nowhere else. W is
  singular at A; elsewhere only the named cyclic quotient points on
  type3 curves can be singular.
- The affine consequences PF2BATT1/PF2R7BND1 remain exactly as stated.
  In particular there can be generically unramified omitted primes.

## Additional primary-source inputs used here

1. Theorem3.5(1), p.13: every type1 residue degree f_i is one. The theorem
   concerns the target P2 model of that paper, as here. It is not applied
   after an arbitrary target blowup or nonlinear target reparametrization.
2. The source-model assertion in the proof of Theorem3.8, pp.14--15:
   for a specified type1 curve the source resolution Z can be chosen so
   that its strict transform is final, hence has square minus one. Our
   transport argument can use separate choices for separate curves.
   Final curves were
   defined on p.7, and minimal source resolutions on p.8. This assertion
   is an explicit external proof input; it is not a consequence of the
   term 'minimal resolution' for our Y. The published theorem's headline
   statement alone would not supply the needed self-intersection bound.

The second item is particularly important. Without that cited source
assertion, a type1 curve on Y could a priori have square below minus one,
and this fixed-infinity parabola reduction would not follow. The local
polynomial certificates do not check that source assertion. A future
review that invalidates or weakens it must reopen this source-relative
route, while preserving the adjoint-rank theorem and algebraic exclusions.
The subsequent [PF2R7EX2 proof](../jc2-r7-quartic-contact/PROOF.md) avoids
that assertion for both nonimmersed classes through cubic and quartic
conic-contact exclusions. Its scoped conclusion has separate evidence.

## Why the change of model is legitimate

Z resolves the same finite normalization W as our minimal smooth Y.
Every resolution of a normal complex surface factors through its minimal
resolution. Thus Z->Y is a birational morphism of smooth projective
surfaces, a composition of point blowdowns. Every removed curve is
contracted by g, because W->P2 is finite. On a smooth intermediate model
the pullback H descends and Q changes by twice the exceptional class.

If C is a surviving type1 curve and v is its nonnegative intersection
with a removed (-1)-curve, then

    Q_new.C_new = Q_old.C_old-2v.

Initially C is a smooth (-1)-curve with H.C=f_i=1, so Q.C=1. Finally Q
is nef. Every v must be zero. This proves that these particular curves
stay smooth with square minus one on Y. We do not transport arbitrary
SNC labels, assume all boundary curves remain smooth, or replace Y by
the singular Stein surface. We do not require nefness on intermediate Z.

## Cusp local covers

The two possible cusp monodromy actions are exactly those in PF2R7BND1.
Equal transpositions give a normal double cover of a plane cusp,
z^2=x^2-y^3, an A2 rational double point. Distinct intersecting
transpositions give the smooth cubic cusp cover. This follows by the
classification of finite analytic covers by monodromy off the branch
curve and uniqueness of normal extension; the standard explicit models
realize the two actions. The A2 resolution has two (-2)-curves; the
smooth cubic model has none. No collided two-branch node projection is
being reintroduced here.

The proof uses the earlier global boundary theorem to separate ramified
points over the nodal and A5 target values. Without that separation,
other source singularities and exceptional counts would have to be kept.

## Other foundations

Named standard inputs are surface Hodge index, adjunction and rationality,
Stein factorization of a pencil, curve Riemann--Hurwitz, minimal surface
resolution and smooth blowdown, and characteristic-zero function-field
degree and quadratic-extension theory. A degree-one finite map to a
normal curve is an isomorphism; a nonconstant polynomial on A2 over C is
surjective onto A1. The latter fact supplies the critical-point
contradiction in the double-line alternative.

The proof imports none of the unrelated generic-cover, Gorenstein,
ordinary-Laplacian or Dixmier claims in the literature review.

## Status of the audit

The model transport is proved here. Source finality and type1 residue
degree are cited primary-source inputs with the exact passage identified.
They are not newly independently established theorems of this repository.
The complete row exclusion is a written deduction relative to these
inputs and JC2POS1--3; finite certificates support only its algebraic
endpoints and numerical identities. Independent end-to-end external
review and formal verification remain absent.

# Experimental geometry: C05--C16

## What is the main theorem?

The programme describes the discriminant, image, omitted values, collision
strata, finite-field behavior, stable normal forms, and dicritical boundary
of the weighted marked-root maps introduced in C04. The central geometric
claim is that generic inverse pencils have controlled rational discriminants,
while special seeds are organized by root-contact partitions.

## Why is it interesting?

It turns a concrete nonproper Keller map into a test family where inverse
fibers, monodromy, omitted values, and divisors at infinity can be computed
uniformly. It also exposes which phenomena are consequences of the marked-root
cover and which depend on the chosen seed.

## What does it depend on?

The whole area depends on the normalized marked-root theorem and `S_n`
monodromy in [C04](../../verified/WEIGHTED_SEED_THEOREM.md). Individual claims
also use normalization, Zariski's Main Theorem, purity, Mason--Stothers,
formal Hensel factorization, valuation theory, and standard plane-curve
singularity formulas.

## What is established?

The repository contains written proofs for C05--C11, C13, C15, and C16, with
exact scripts checking their finite identities and representative models.
Completion and review information is maintained only in
[STATUS.md](../../STATUS.md).

## Independent checks

C12's completed-local model and C14's quartic geometry now have complete
internal proofs plus dependency-free clean-room checks of their decisive
algebra.  An unrelated-CAS reproduction of the full C12 tangent matrices and
the complete C14 radical/path suite would still strengthen them further.
Bounded degree, multiplicity, and finite-field runs elsewhere remain
regressions and do not prove uniform quantifiers.

## What is the likeliest failure point?

The highest-risk steps are global exhaustiveness in the discriminant and
dicritical compactifications, plus descent from completed local
factorizations to global scheme statements. A failure would most likely be a
missing boundary component, an unstated genericity hypothesis, or a
normalization/conductor issue—not one of the displayed polynomial identities.

## What review expertise is needed?

Birational and logarithmic geometry, normalization of finite covers, plane
curve singularities, arithmetic geometry for C13, and computational algebra
for the exact elimination certificates.

Suggested entry points:
[generic discriminants](GENERIC_DISCRIMINANT_CURVE.md),
[contact strata](COINCIDENT_ROOT_REBUILD.md), and
[dicritical compactification](DICRITICAL_COMPACTIFICATION.md).

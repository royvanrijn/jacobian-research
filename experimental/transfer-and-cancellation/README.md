# Transfer and cancellation: C17--C24

**Workflow state:** C17--C23 are frozen. C24 is the repository's only active
research branch.

## What is the main theorem?

C17--C23 describe completed square/cube transfer blocks and their global
factorization equalizers. C24 constructs the master cancellation family,
computes its generic inverse degree and monodromy, and begins its
classification under polynomial left--right equivalence.

## Why is it interesting?

The transfer blocks isolate the local algebra behind collisions of
factorization strata. The cancellation family supplies a broader source of
explicit noninjective Keller maps and asks whether visibly different inverse
resolvents can define equivalent polynomial maps.

## What does it depend on?

The local transfer theory uses formal Hensel factorization, relative Groebner
bases, symmetric-group invariant theory, and completed tensor products. C24
uses its explicit reconstruction formulas, finite normalization, DVR
ramification, and tame branch cycles. Its `(1,1)` member links back to
[C01](../../verified/FOUNDATIONAL_GEOMETRY.md).

## What is fully proved?

Internal proofs are recorded for C17, C18, and C21--C23. The construction,
degree, reconstruction, monodromy, and `r>=2` ramification obstruction in C24
are also proved internally. These results have not received a full
independent audit.

## What is only computationally tested?

The explicit `Z_3` and `Z_4` presentations (C19--C20), higher transfer-block
regressions, parameter-polynomial factorizations, and displayed small C24
maps are computational evidence. They do not replace the written uniform
arguments.

## What is the likeliest failure point?

For the frozen transfer theory, the risk is descent from completed Boolean
models to the global affine equalizer. For active C24 work, the risk is
mistaking a primitive-resolvent artifact for an invariant of the polynomial
map. Every proposed obstruction must be identified on the canonical finite
normalization.

## What review expertise is needed?

Finite covers and valuation theory, permutation monodromy, formal and
derived local algebra, invariant theory, and polynomial automorphism/stable
equivalence techniques.

## Active milestone

**Prove a coordinate-independent invariant showing when two
resolvent-cancellation constructions cannot be polynomially left--right
equivalent.**

The unresolved case is C24 with `r=1`, `m>1` versus a generic weighted seed.
See [the signature](RESOLVENT_RAMIFICATION_SIGNATURE.md) and the
[precise open problem](OPEN_PROBLEMS.md). No other branch in this directory
is active.

# Transfer and cancellation: C17--C24

**Workflow state:** C24 is the main active branch. C22 is reopened only for a
corrective audit; C18, C21, and the all-`k` portion of C23 are conditional on
that repair.

## What is the main theorem?

C17--C23 describe completed square/cube transfer blocks and their global
factorization equalizers. C24 constructs the master cancellation family,
computes its generic inverse degree and monodromy, and begins its
classification under polynomial left--right equivalence.  The
[generalized cancellation note](GENERALIZED_CANCELLATION_MECHANISM.md)
extracts the determinant argument for arbitrary one-variable input `f(y)`
and reduces polynomiality to a finite functional equation.

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

C17 and the affine-difference Wronskian argument remain internally proved.
The former all-`k` Boolean model in C22 is refuted by an explicit second-order
deformation, making the dependent C18/C21 conclusions conditional. The
construction, degree, reconstruction, monodromy, and `r>=2` ramification
obstruction in C24 remain internally proved.

## What is only computationally tested?

The explicit `Z_3` and `Z_4` presentations (C19--C20), transfer-block ranks
through the bounded cases, parameter-polynomial factorizations, and displayed
small C24 maps are computational evidence. They do not prove an all-`k`
transfer theorem.

## What is the likeliest failure point?

For transfer theory, the Boolean norm is now known to omit collision
directions, beginning with `X^2` at `k=2`. For active C24 work, the risk is
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
[precise open problem](OPEN_PROBLEMS.md).  A second C24-only target is the
common-root and polynomial-descent problem isolated in the
[generalized mechanism](GENERALIZED_CANCELLATION_MECHANISM.md).  No other
branch in this directory is active.

The C22 corrective target is stated in the
[conductor-ribbon counteraudit](C22_CONDUCTOR_RIBBON_AUDIT.md): construct a
direct all-base straightening basis in the factorization ring, without using
the refuted norm isomorphism.

# Transfer and cancellation: C17--C24

**Workflow state:** C24 is the sole active research branch.  The proposed
C22 transverse filtration is refuted at `k=2`; C18, C21, and all
transfer-dependent all-`k` claims have been archived.  C23's independent
affine-factorization-complement results remain in scope.

## What is the main theorem?

C17, C19, and C20 now supply only exact low-`k` local-algebra examples.
[C23](UNIVERSAL_FACTORIZATION_GEOMETRY.md) independently studies when
universal marked-factorization complements are affine spaces. C24
constructs the master cancellation family,
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

The retained low-`k` transfer examples use exact relative Groebner bases. C24
uses its explicit reconstruction formulas, finite normalization, DVR
ramification, and tame branch cycles. Its `(1,1)` member links back to
[C01](../../verified/FOUNDATIONAL_GEOMETRY.md).

## What is fully proved?

C17 remains an exact low-degree theorem.  The independent cubic C23
classification and `(2,3)` obstruction remain internally proved.  The
construction, degree, reconstruction, monodromy, and `r>=2` ramification
obstruction in C24 remain internally proved.  The new
[quadratic-remainder audit](QUADRATIC_REMAINDER_ALGEBRA.md) proves that the
proposed C22 associated graded is false.

## What is only computationally tested?

The explicit `Z_3` and `Z_4` presentations (C19--C20), parameter-polynomial
factorizations, and displayed small C24 maps are computational evidence.
They are retained as bounded examples, not evidence for an active all-`k`
programme.

## What is the likeliest failure point?

For transfer theory, the quadratic-remainder tangent cone is now known to
omit `X^3` at `k=2`, while the older Boolean norm omits `X^2`.  For active
C24 work, the risk is
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

The archived all-`k` transfer programme is preserved under
[archive/transfer-all-k](../../archive/transfer-all-k/).  It is not an active
dependency chain.  Only the exact low-degree examples in
[TRANSFER_BLOCKS.md](TRANSFER_BLOCKS.md) remain in current navigation.

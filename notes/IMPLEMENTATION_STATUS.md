# Implementation status

## Exact core

- Two independent implementations verify `det DF = -2`, the three-point
  rational collision, and coordinate degrees `(7,6,4)`.
- The dependency-free verifier uses standard-library rational arithmetic and
  its own sparse polynomial representation.
- The primitive cubic, rational reconstruction, and discriminant identity are
  checked symbolically.
- Exceptional fibers, the exact image, and both inclusions in the stated
  nonproperness set have executable checks.

## Extended three-dimensional analysis

- Finite-field distributions and collision refinements are implemented.
- The commuting inverse-Jacobian frame, low-degree symmetries, and tangent
  constraints are computed exactly.
- The three-minimum polynomial, its Palais--Smale escape curve, and an
  adversarial homotopy benchmark are implemented.
- A general weighted-seed constructor is implemented. The quartic-sheet model
  has exact verifiers for polynomiality, determinant, collision, inverse,
  discriminant, special fibers, nonproperness, singular locus, and image.
- A bounded scan of 627 weighted seeds is exploratory and explicitly separated
  from theorem-level claims.

## Explicit normal forms

- The Bass--Connell--Wright/Yagzhev construction produces a 95-dimensional
  cubic-homogeneous map with an exact transported collision.
- The Gorni--Zampieri/Drużkowski pairing produces a 510-dimensional
  cubic-linear map with an exact transported collision.
- Both large artifacts are generated locally and checked by dedicated
  verifiers; they are ignored because they are reproducible and large.

## Remaining audit work

- Produce a short structural determinant proof suitable for hand checking.
- Obtain an independent second-CAS verification and archival reproduction log.
- Complete the provenance record and separate priority from verification.
- Independently audit irreducibility, normalization, monodromy, and every
  boundary stratum used by the strongest global-geometry claims.
- Audit the stable-equivalence steps and external implication theorems against
  their primary sources.
- Generalize the proved quartic weighted geometry only after stating precise
  hypotheses on primitive zeros and ramification.

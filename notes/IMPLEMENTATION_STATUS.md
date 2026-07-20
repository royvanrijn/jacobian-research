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
  discriminant, special fibers, nonproperness, singular locus, image, and full
  geometric/arithmetic monodromy `S_4`.
- The universal inverse-pencil theorem now proves irreducibility, birational
  discriminant normalization, and geometric/arithmetic monodromy `S_n` for
  every characteristic-zero weighted seed; canonical and deformed degrees
  through eight have an exact regression audit.
- The canonical family `H_d=W^d(1-W)` now has a uniform image and
  nonproperness theorem: the inverse-degree-three and four members have one
  omitted curve each, all inverse degrees at least five are surjective, and
  the `C=0` boundary fibers and saturation factors are explicit.
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
- Independently audit normalization and every boundary stratum used by the
  strongest global-geometry claims; the quartic irreducibility and monodromy
  now have exact certificates.
- Audit the stable-equivalence steps and external implication theorems against
  their primary sources.
- Extend the image and nonproperness theorem to deformed seeds, keeping every
  additional primitive-zero boundary branch explicit.

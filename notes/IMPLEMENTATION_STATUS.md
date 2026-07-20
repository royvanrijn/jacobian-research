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
- Boundary-clean deformed seeds now have a uniform `C=0` fiber and
  nonproperness theorem. The complete one-extra-simple-zero family has an
  explicit omitted double-double curve, including its quadruple-root
  degeneration.
- Omitted values for arbitrary fixed seeds now have a complete exact
  multiplicity-partition classifier. For two additional simple zeros, one
  explicit coefficient equation detects the exceptional nonsurjective
  subfamily; generic two- and three-extra-root representatives are surjective.
- Repeated extra primitive roots now have an exact boundary theorem. Their
  saturation exponent, boundary trace, direct fibers, and escaping-branch
  rates are determined by the complete root-multiplicity profile, including
  nonsplit factors.
- Universal `S_n` monodromy now yields a good-reduction finite-field
  Chebotarev theorem: degree-`n` fiber sizes converge to the fixed-point law of
  `S_n`, with target counts `p_{n,j}q^3+O(q^(5/2))` and limiting image density
  `1-D_n/n!`.
- The generic degree-`n` discriminant theorem now identifies the curve as a
  rational dual with `n-2` cusps and `(n-2)(n-3)/2` nodes in every degree. A
  uniform contact-incidence dimension proof and tangent-chord normalization
  show that the good locus meets every admissible degree. Thus generic seeds
  are surjective for all inverse degrees at least five. Exact rational seeds
  through degree ten remain as independent regression certificates.
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
- Turn the dual-curve incidence proof into explicit saturated universal ideals
  and compare it archivally with classical projective duality, then make
  good-prime exclusions effective for individual seeds.

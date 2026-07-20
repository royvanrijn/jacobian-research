# Implementation status

This document uses three distinct kinds of support:

1. **Executable certificate:** an exact computation directly establishes the
   stated finite identity or example.
2. **Uniform written proof:** a conventional mathematical argument establishes
   a quantified theorem in every degree.  Scripts may check its algebraic
   lemmas, but do not constitute a formal proof assistant derivation.
3. **Regression test:** exact bounded-degree examples test an implementation
   and guard the uniform argument against algebraic mistakes; they do not
   establish the all-degree quantifier.

Unless explicitly called an executable certificate, the all-degree statements
below have status (2), generally accompanied by support of type (1) for key
identities and type (3) in selected degrees.

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
- A uniform written proof of the universal inverse-pencil theorem establishes
  irreducibility, birational discriminant normalization, and
  geometric/arithmetic monodromy `S_n` for every characteristic-zero weighted
  seed.  Canonical and deformed degrees through eight have exact regression
  tests; those tests do not prove the all-degree statement.
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
- A uniform written proof of the generic degree-`n` discriminant theorem
  identifies the curve as a rational dual with `n-2` cusps and
  `(n-2)(n-3)/2` nodes in every degree. A
  uniform contact-incidence dimension proof and tangent-chord normalization
  show that the good locus meets every admissible degree. Thus generic seeds
  are surjective for all inverse degrees at least five. Exact rational seeds
  through degree ten remain as independent exact regression tests.
- The higher-cusp, cusp-plus-branch, tritangent, and ordinary-bitangent strata
  now have reusable universal ideals with explicit Rabinowitsch saturation by
  every diagonal, cusp, degree, and weighted-admissibility factor.
- The incidence API now accepts every contact partition, quotients equal parts
  before elimination, and supports exact residual factors. It recovers the
  degree-five polynomial `F(R,P)`, identifies the degree-six main locus as a
  rational irreducible quartic surface, and proves that degree-seven
  nonsurjectivity first occurs in codimension two.
- The written uniform exceptional-seed proof identifies the nonsurjective locus as
  the union of all full-contact strata and proves
  `dim E_lambda=ell(lambda)-1` by a weighted-Vandermonde determinant.  The
  multiple-omission API separates common collision values from genuinely
  distinct omitted values.  Exact degree-six and degree-eight calculations
  find only the predicted common collision boundaries and no off-diagonal
  intersections.
- Merging partition parts defines an executable combinatorial collision-order
  calculation.  The written tangent-chord deformation proves every coarser stratum lies in the closure
  of every refining stratum, and Mason--Stothers uniformly rules out
  off-collision two-omission solutions for distinct maximal 2/3 partitions.
- The maximal root hypersurfaces `Phi_(2^a 3^b)` are uniformly irreducible by
  a written proof: primitive-linear for `b>=3`, nonsquare-quadratic for `a>=3`,
  and has seven exact endpoint-rank certificates.  Thus maximal 2/3 types
  index the actual irreducible components of the exceptional-locus closure.
  Their dimensions are `a+b-1`, the full exceptional codimension is
  `ceil(n/2)-2`, and component intersections are exactly the common collision
  coarsenings.
- The contact-atom theorem explains the occurrence of twos and threes as the
  indecomposable elements of the allowed multiplicity semigroup.  Its exact
  threshold-`r` generalization is implemented.  The excess identity sharpens
  Mason separation from maximal types to every pair of distinct full-contact
  partitions.
- The unique omitted-value theorem closes Mason's sole support-equality case:
  two distinct monic all-double polynomials are squares whose difference has
  degree at least two.  Every seed therefore omits at most one normalized
  pencil value, and the exact contact strata form a disjoint stratification
  of the nonsurjective locus.
- Every exceptional component now has an explicit smooth normalization: the
  admissible `Q^2R^3` quotient hypersurface is smooth even across collision
  diagonals, and its seed map is finite and generically degree one.  Collision
  fiber cardinalities are computed by allocations `2i+3j=m` at each root.
- A bounded scan of 627 weighted seeds is exploratory and explicitly separated
  from theorem-level claims.

## Explicit normal forms

- The Bass--Connell--Wright/Yagzhev construction produces a 95-dimensional
  cubic-homogeneous map with an exact transported collision.
- The Gorni--Zampieri/Drużkowski pairing produces a 510-dimensional
  cubic-linear map with an exact transported collision.
- Both large artifacts are generated locally, checked by dedicated verifiers,
  and retained as tracked reference certificates. Newly generated result files
  remain ignored by default.

## Remaining audit work

- The README now gives a short hand-checkable structural determinant proof
  through the `(t,r,c)` inverse chart.
- Obtain an independent second-CAS verification and archival reproduction log.
- Complete the provenance record and separate priority from verification.
- Independently audit normalization and every boundary stratum used by the
  strongest global-geometry claims; the quartic irreducibility and monodromy
  now have exact certificates.
- Independently re-audit the stable-equivalence steps and recorded external
  implications against their primary sources.
- Compare the dual-curve incidence proof archivally with classical projective
  duality, then make good-prime exclusions effective for individual seeds and
  eliminate the universal exceptional strata in selected low degrees.
- Refine the set-theoretic common-coarsening formula to scheme-theoretic
  intersections, including embedded components and collision multiplicities.
- Compute normalization ramification, conductor ideals, and local branch
  multiplicities at special points of each collision stratum; the generic
  geometric branch count is now the coefficient of
  `U^a V^b` in `product_rho sum_(2i+3j=m_rho) U^i V^j`.

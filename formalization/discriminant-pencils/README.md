# Lean formalization: polynomial tangent pencils

> **Status: partial formalization.** The Lean development is axiom-free and
> checks the results listed below, but it does **not** yet prove the paper's
> main theorem in full. In particular, it must not be cited or described as
> a completed machine-checked proof of the paper.

This package formalizes the algebraic core of
*Generic Discriminants of Polynomial Tangent Pencils*.

Build it with:

```sh
lake build
```

## Checked results

- `TangentPencil.lean`
  - the repeated-root parametrization
    `r ↦ (H'(r), r H'(r) - H(r))`;
  - the coordinate derivative identity
    `ν'_H(r) = H''(r) (1, r)`;
  - the tangent-direction determinant;
  - the cusp-jet determinant `2 H'''(r)^2`.
- `LocalSingularities.lean`
  - contact at least two with the tangent line;
  - exact second and third jets;
  - nonvanishing of the ordinary-cusp determinant;
  - transversality of two distinct unramified branches.
- `ContactStrata.lean`
  - the factorization family `H - ℓ = M_μ Q`;
  - divisibility by every marked contact factor;
  - the source-dimension calculation `n - 2` for `(4)`, `(3,2)`,
    and `(2,2,2)`.
- `ProjectiveParametrization.lean`
  - the three homogeneous coordinate forms from equation (2.1);
  - recovery of the affine parametrization on `U = 1`;
  - evaluation of homogenized forms at `U = 0`;
  - absence of base points for exact degree `n ≥ 2`;
  - the unique source point over the line at infinity.
- `InfinityLocal.lean`
  - fixed-degree reversal equals evaluation of the homogenized coordinates
    on the source chart `[1:q]`;
  - the `T` numerator is a unit in `k⟦q⟧`;
  - the two expansions in equation (2.2), including
    `ord(S/T) = 1` and `ord(Z/T) = n`;
  - a two-sided compositional inverse for `S/T`;
  - a formal series `G` with `Z/T = G(S/T)`, giving the formal-graph
    certificate used for smoothness at infinity.
- `CuspCount.lean`
  - a squarefree `H''` has `n - 2` distinct roots over an algebraically
    closed characteristic-zero field;
  - each such root is an ordinary cusp parameter.
- `CuspImages.lean`
  - the no-ramified-collision condition makes the tangent map injective on
    the cusp parameters;
  - consequently there are exactly `n - 2` distinct cusp image points.
- `BadContactBridge.lean`
  - a zero of `H''` gives divisibility by `(X-r)^3`;
  - an equal-image second parameter gives divisibility by `(X-u)^2`;
  - for `r ≠ u`, coprimality combines these into the bad `(3,2)` contact
    factor;
  - excluding `(3,2)` therefore supplies the no-collision hypothesis used
    by the distinct-cusp-image count;
  - simultaneous vanishing of `H''` and `H'''` gives the bad `(4)` factor;
  - three distinct equal-image parameters give the bad `(2,2,2)` factor.
- `AffineExhaustion.lean`
  - the three direct contact-factor exclusions imply simple ramification,
    absence of ramified collisions, and fibers of cardinality at most two;
  - every affine normalization fiber is consequently a singleton
    unramified branch, a singleton ordinary-cusp jet, or two unramified
    branches with transverse velocities;
  - the combined theorem `affineClassification_of_contactExclusions`
    includes the exact count of `n - 2` distinct cusp images.
- `NumericalCount.lean`
  - the final genus arithmetic
    `N = (n - 2)(n - 3)/2`, once the doubled genus formula is supplied.

All these files are theorem-checked and contain no `sorry` or `axiom`.

## What is not yet formalized

This is not yet a complete formalization of the paper's main theorem.
The following global plane-curve results used by the manuscript do not
currently exist in mathlib at the required level:

1. the dimension-of-image theorem needed to turn the three
   `n - 2`-dimensional contact families into proper closed subsets of the
   `n - 1`-dimensional coefficient space;
2. the theorem relating a base-point-free birational map from `P¹` given by
   degree-`n` forms to an integral image curve of degree `n`, together with
   identification of `P¹` as its normalization;
3. formal local-germ criteria identifying the checked fiber/jet
   classification with smooth points, ordinary plane cusps, and ordinary
   nodes.  The normalization-fiber exhaustion, transverse velocities, and
   formal graph/uniformizer calculation at infinity are checked, but they
   still need to be connected to mathlib's scheme predicates;
4. the arithmetic-genus/normalization delta formula for integral plane
   curves, plus the computation that ordinary nodes and cusps have delta
   invariant one.

No placeholders for these theorems are introduced here: in particular, the
package does not disguise them as axioms. A full proof requires developing
these foundations first.

## Recommended continuation point

The next coherent block is the global genericity theorem: define the three
incidence morphisms for `(4)`, `(3,2)`, and `(2,2,2)` on coefficient affine
space, then prove that the closures of their images are proper.  This needs
a reusable finite-type image-dimension result before any further
paper-specific Lean code will close the main theorem.

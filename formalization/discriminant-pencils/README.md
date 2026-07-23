# Lean formalization: polynomial tangent pencils

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
- `CuspCount.lean`
  - a squarefree `H''` has `n - 2` distinct roots over an algebraically
    closed characteristic-zero field;
  - each such root is an ordinary cusp parameter.
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
3. formal local-germ criteria identifying the checked nonzero jet
   determinants with ordinary plane cusps and ordinary nodes, including
   smoothness of the formal graph at infinity;
4. the arithmetic-genus/normalization delta formula for integral plane
   curves, plus the computation that ordinary nodes and cusps have delta
   invariant one.

No placeholders for these theorems are introduced here: in particular, the
package does not disguise them as axioms. A full proof requires developing
these foundations first.

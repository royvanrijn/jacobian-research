# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Every Nonzero Finite Étale Algebra Except Rank Two Is a Keller Fiber*

This manuscript remains under active development. The directory name is a
stable path retained from its earlier draft. The complete constructive
polynomial-presentation layer is formalized in `formal/finite-etale-keller`:
automatic translation choice, the actual arbitrary-degree map, its Jacobian
and effective degree bound, the literal scheme fiber, quotient transport, and
naturality. The manuscript and its verification matrix state the remaining
boundary precisely: generic function-field degree, monogenicity, the classical
rank-two input, and the number-theoretic applications.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

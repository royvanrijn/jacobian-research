# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Every Finite Étale Algebra Except Rank Two Is a Keller Fiber*

This manuscript remains under active development. The directory name is a
stable path retained from its earlier draft. Its scheme-reconstruction core,
existence of an admissible translation, quotient transport, and naturality are
formalized in `formal/finite-etale-keller`. The manuscript and its README state
the exact boundary between that Lean certificate, the coefficientwise
all-degree displayed-map proof, its exact regression checks, and the classical
rank-two input.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Finite Étale Algebras as Keller Fibers*

This manuscript remains under active development. The directory name is a
stable path retained from its earlier draft. Its scheme-reconstruction core is
now formalized as a natural represented-fiber theorem in
`formal/finite-etale-keller`; the manuscript and its README state the exact
boundary between the Lean certificate, the general displayed-map identities,
and the classical rank-two input.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

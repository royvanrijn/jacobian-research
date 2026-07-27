# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Prescribed Finite Étale Algebras as Full
  Fibers of Keller Maps with Symmetric Monodromy*

This manuscript remains under active development. The directory name is a
stable path retained from its earlier draft. The complete constructive
polynomial-presentation layer is formalized in `formal/finite-etale-keller`:
automatic translation choice, the actual arbitrary-degree map, its Jacobian
and effective degree bound, the literal scheme fiber, quotient transport, and
naturality. The same Lean project now proves monogenicity for abstract finite
étale algebras in characteristic zero and composes it with that construction.
The manuscript derives absolute and stable compositional atomicity from its
symmetric-monodromy theorem. It additionally proves the exact reduced
nonproperness locus and boundary-sheet ledger as ordinary mathematics; the
verification matrix marks these results explicitly as not Lean-formalized.
Hasse and fixed-quintic applications are maintained as separate verified
notes rather than part of this paper.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

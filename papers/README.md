# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`
- `three-pair-image-counterexample`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Every Finite Étale Algebra of Rank at Least
  Three Is a Full Keller Fiber*
- `quadratic-gauge-nonproperness` — *The Exact Nonproperness Locus of the
  Quadratic-Gauge Keller Maps*

This manuscript remains under active development. The directory name is a
stable path retained from its earlier draft. The complete constructive
polynomial-presentation layer in characteristic zero is formalized in
`formal/finite-etale-keller`: automatic translation choice, the actual
arbitrary-degree map, its Jacobian and effective degree bound, the literal
scheme fiber, quotient transport, and naturality. The same Lean project proves
monogenicity for abstract finite étale algebras in characteristic zero and
composes it with that construction.
The manuscript is focused on the prescribed-fiber theorem: supplied
separable polynomial presentations in characteristic different from two,
inverse reconstruction, literal fibers, function-field degree, monogenicity,
base change, and one explicit arithmetic example. Lean formalizes the
complete characteristic-zero specialization and abstract finite-étale
corollary. Symmetric monodromy and stable atomicity remain in their verified
companion notes. Exact nonproperness and boundary-sheet accounting have moved
to the separate geometric manuscript `quadratic-gauge-nonproperness`.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

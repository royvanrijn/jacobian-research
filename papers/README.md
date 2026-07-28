# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`
- `three-pair-image-counterexample`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `common-arithmetic-fibers` — *Over Characteristic Zero, Every Finite Étale
  Algebra of Rank at Least Three Is a Full Keller Fiber*
- `fixed-map-hasse-failures` — *Quantitative Hasse-Principle Failures in the
  Fibers of a Fixed Keller Map*

`common-arithmetic-fibers` remains under active development. Its directory
name is a stable path retained from an earlier draft. The complete constructive
polynomial-presentation layer in characteristic zero is formalized in
`formal/finite-etale-keller`: automatic translation choice, the actual
arbitrary-degree map, its Jacobian and effective degree bound, the literal
scheme fiber, quotient transport, and naturality. The same Lean project proves
monogenicity for abstract finite étale algebras in characteristic zero and
composes it with that construction.
The manuscript is focused on the prescribed-fiber theorem: its headline
characteristic-zero polynomial-presentation theorem matches the end-to-end
Lean theorem, while the broader characteristic-not-two supplied-presentation
result is a separate proposition not yet formalized at that generality.
Inverse reconstruction, literal fibers, function-field degree, monogenicity,
base change, and one explicit arithmetic example remain in scope. Symmetric
monodromy and stable atomicity remain in their verified companion notes.
Exact nonproperness and boundary-sheet accounting have moved to the separate
geometric manuscript `quadratic-gauge-nonproperness`.

`fixed-map-hasse-failures` is the focused arithmetic sequel. It fixes one
Jacobian-one map of geometric degree five, proves a uniform family of full
Hasse-failing fibers, and counts the constructed targets by height using
Selberg--Delange.

## Companion draft

- `quadratic-gauge-nonproperness` — *The Exact Nonproperness Locus of the
  Quadratic-Gauge Keller Maps*

This geometric companion is not part of the current two-paper publication
sequence and is excluded from `make verify-papers`. It remains in
`make clean-papers` and can be compiled directly.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

# Paper registry

The manuscript directories remain at stable paths. Their workflow status is:

## Finalized / frozen / preprint

- `gaussian-moments-two-variables`
- `sparse-minimality-gaussian-moments-dimension-three`
- `three-pair-image-counterexample`

These manuscripts are frozen preprints. Changes should be limited to
corrections that are deliberately carried into a new deposited version.

## Active

- `generalized-vanishing-two-variables` — *The Generalized Vanishing
  Conjecture: The Two-Variable Theorem and the First Failing Dimension*
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

`generalized-vanishing-two-variables` is the concise paper version of the
Hall-envelope proof.  It includes the Hall localization, shifted-ray
prime-dilation argument, common-threshold cutoff, and global moving-envelope
proof, followed by the homogeneous three-variable counterexample, its
connection to Long's Gaussian construction, the complete
winding--profile--radial family of failures, and the exact dimensional
classification of GVC.  It is an active internal draft and has not been
externally reviewed.

The companion `formal/gvc` package is currently a partial Lean audit, not a
complete certificate of the headline theorems.  It verifies the concrete
cusp identity and all-order endpoint-coefficient mechanism and makes the
remaining counterexample and binary-envelope bridges explicit; see its
README for the exact boundary.

## Companion draft

- `tschirnhaus-keller-non-descent` — *Arithmetic Descent and Geometric
  Non-Descent for Universal Keller Fibres* (prospectus)

- `quadratic-gauge-nonproperness` — *The Exact Nonproperness Locus of the
  Quadratic-Gauge Keller Maps*

This geometric companion is not part of the current two-paper publication
sequence and is excluded from `make verify-papers`. It remains in
`make clean-papers` and can be compiled directly.

`tschirnhaus-keller-non-descent` currently contains a theorem-and-scope
prospectus rather than a TeX manuscript.  It assembles the generic
non-descent theorem, the rank-five transition calculation, and the clean
decorated receiver while keeping the full marked stabilizer problem open.

## Parked

- `exact-real-chamber-spectra`
- `discriminant-pencils`

These manuscripts are retained as research archives, but are excluded from
the normal `make verify-papers` build and publication pipeline. They can still
be compiled directly from their source directories.

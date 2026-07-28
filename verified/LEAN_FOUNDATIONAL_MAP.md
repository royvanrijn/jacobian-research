# External Lean formalization of the foundational Keller map

Dean Cureton's [`deancureton/jacobian`](https://github.com/deancureton/jacobian)
is a separate Lean 4 formalization of the announced counterexample.  This
repository uses it as an external, independently implemented formal
certificate for the [Foundational Keller map](../STATUS.md), pinned at commit
[`0d4a9212d874226ad81ce5a926becddfa94e6a88`](https://github.com/deancureton/jacobian/commit/0d4a9212d874226ad81ce5a926becddfa94e6a88).

## Exact scope

The upstream file `Jacobian/Counterexample.lean` defines multivariate-polynomial
Jacobian matrices, determinants, and evaluation maps, then proves:

- `jacobianDet_F`: the displayed map `F` used here has determinant `-2`;
- `evalMap_F_p0`, `evalMap_F_p1`, and `evalMap_F_p2`: the three displayed
  rational points have the common image `(-1/4, 0, 0)` when `2 != 0`;
- `not_jacobianConjecture`: over any field with `2 != 0`, the unit-Jacobian
  implication to injectivity fails;
- `jacobianDet_G`: a diagonally rescaled form `G` has determinant `1`;
- `not_jacobianConjecture_all_char`: the analogous injectivity statement fails
  over every field, using a separate collision in characteristic two; and
- `not_jacobianConjecture_complex`: the specialization over `C`.

The all-characteristic theorem is a statement about the displayed algebraic
property over arbitrary fields.  The classical Jacobian conjecture itself is
the characteristic-zero statement.

This formalization covers the determinant and collision of the foundational
map.  It does not formalize this repository's cubic inverse geometry,
weighted families, or stable normal-form reductions.

## Reproduction

Run:

```text
make verify-lean-foundational
```

The target fetches the exact upstream commit into the ignored source checkout
`.cache/lean-foundational`, creates the generated worktree
`.cache/lean-foundational-stable`, and applies the checked
[`v4.32.1` stable-toolchain overlay](lean-foundational-v4.32.1.patch).  The
overlay changes only `lean-toolchain`, `lakefile.toml`, and
`lake-manifest.json`; the wrapper verifies the exact Git blob hash of each
overlaid file before downloading the stable Mathlib cache and running
`lake build`.  The audited Lean source remains byte-for-byte at the upstream
commit.

There is no Lean or Mathlib `v4.33.0` stable tag as of 28 July 2026, so the
repository uses the newest available stable patch release, Lean/Mathlib
`v4.32.1`, rather than retaining the upstream `v4.33.0-rc1` build.  It is
intentionally not part of `make verify`, because the first run downloads this
pinned stable toolchain and its Mathlib dependencies.

The stable-overlaid target was reproduced successfully on 28 July 2026.  The
upstream `#print axioms` checks report only `propext`, `Classical.choice`, and
`Quot.sound` for each of the three final theorems; in particular, the
reported dependency lists contain no `sorryAx`.

## Attribution and source boundary

The Lean definitions and proofs are Dean Cureton's work.  His README titles
the project “Levent Alpöge/Fable 5's counterexample to the Jacobian conjecture
in Lean 4,” preserving the construction attribution used by its author.  The
audited upstream commit contains no `LICENSE`, `COPYING`, or `NOTICE` file.
Accordingly, this repository links to and builds the upstream project but does
not copy or relicense its Lean source.  Attribution does not substitute for a
license; vendoring or adapting the code should wait for explicit permission.

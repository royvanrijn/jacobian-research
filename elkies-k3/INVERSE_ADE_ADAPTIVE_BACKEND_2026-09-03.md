# Effective rank-15 inverse-ADE backend

Date: 2026-09-03.

## Outcome

An adaptive implementation of the projective birth theorem is now available
in
[`scripts/benchmark_inverse_ade_adaptive_backend.sage`](scripts/benchmark_inverse_ade_adaptive_backend.sage).
It has three exact modes:

1. `expanded` enumerates the required scaled dual-coset shells once and hashes
   their nonzero projective reductions;
2. `orbit` compresses that projective support by the subgroup of `Aut(K)` that
   fixes every core class used by the graph glue;
3. `lazy` constructs no shell and makes exact affine-CVP queries for each
   candidate line, stopping at the first occupied graph cell.

The automatic selector uses a geometry-of-numbers representation estimate
only to choose a backend.  It never uses that estimate as a root decision.
Every accepted lazy rootless line has exhausted every old-root hyperplane and
every required affine graph cell, so its complete predicted physical root set
is empty before any child lattice is constructed.

## Terminal controls and the corrected recovery gate

Both historical terminal parents have only the sign automorphism.  Sign is
projectively trivial and does not fix the nonzero marked graph classes, so the
compatible subgroup has order one.  The automatic selector correctly chooses
lazy CVP rather than attempting shell or orbit expansion.

| control | prime | required graph cells | blind exact predictions | blind affine queries | search-time child constructions | final child constructions | reduction |
|---|---:|---:|---:|---:|---:|---:|---:|
| H3 `3A1 -> rootless` | 11 | 33 | 25,860 | 220,554 | 0 | 2 | 1,500x |
| Q80 `2A1 -> rootless` | 29 | 21 | 25,879 | 322,052 | 0 | 2 | 1,500x |

The backend input records contain only the parent `K`, bridge `C`, marked
graph glue `H`, prime `p`, and desired rootless ADE signature.  They contain
no marked target core, target Gram matrix, target-overlap fingerprint, or
historical line.  With the same fixed native-Kneser seed `1`, the blind lazy
searches return

```text
H3:  [1,2,4,3,8,6,6,5,7,8,8,3,0,7,4]
Q80: [1,10,17,0,2,13,22,8,18,20,27,28,17,24,1]
```

after 25,860 and 25,879 exact predictions respectively.  Each accepting query
predicts the complete empty root set before any child is constructed.  One
transported-glue child per selected line is then independently rootless.

Only after the corresponding target-free search finishes is its historical
line revealed as a second held-out predicate query and materialized once.
Those regressions also predict and obtain the complete empty root set.

The blind lines differ from the historical lines on both controls.  This is
not a weakened computational result; it proves that the former gate demanding
the exact historical projective coordinate was ill-posed.  The declared input
`K,C,H,p,ADE` specifies a solution set and contains no datum that distinguishes
one of two independently certified rootless lines.  The corrected extensional
gate is therefore blind recovery of a satisfying terminal line, while exact
historical-coordinate equality remains a reported, non-required diagnostic.

The exact nonsingular projective quadrics have respectively

```text
37,974,983,358,324
10,627,079,738,421,409,410
```

points.  Two post-prediction materializations (one blind discovery and one
held-out historical regression) are therefore far below the requested 5--10%
ceiling on each control, and reduce the old 3,000-materialization terminal
window by a factor of 1,500.

## Nonhistorical completion

A blind lazy-CVP search on the rank-15 terminal NS0024 parent, with fixed
repository seed `314159`, returns at `p=7`

```text
[1, 3, 5, 2, 0, 3, 1, 4, 3, 4, 4, 5, 1, 5, 0].
```

The search made 22,094 isotropic proposals, rejected 13,627 by old-root
hyperplanes, and made 8,467 exact birth-locus predictions using 112,710 affine
CVP queries.  It constructed no child during the search.  The accepting query
exhausted all 567 graph-cell/layer queries and predicted the complete empty
root set.  Only then was one child completion constructed.  Its root set is
again empty, and its reduced core Gram hash is

```text
sha256:787020c4b1ecca41d9f1fdd3853f0a8fd9f1c659147b6a36c3ee100b7dff9015.
```

The selected line differs from the historical NS0024 terminal line, which was
not consulted until after the search returned.

## Three-mode control and conjecture evidence

The existing diagonal half-glue control (`A1^3 + A1`, parent completion `D4`,
`p=5`) is replayed through all three backends.  All modes make identical exact
rootlessness decisions on all six isotropic lines.  Full expansion has 87
signed scaled-shell vectors and six projected birth points.  The compatible
automorphism subgroup compresses those six points to one projective orbit, so
`auto` selects `orbit` on this control.

For the three rank-15 parents actually measured here (H3, Q80, and NS0024),
`Aut(K)` has order two and gives no useful graph-compatible projective
compression.  Thus these examples are negative evidence for the strongest
universal automorphism-compression conjecture, while also demonstrating why
the lazy backend is necessary.  The finite-quadratic-module/Fourier or direct
shell-orbit-representative versions remain open.  Current orbit mode compresses
support after exact shell enumeration; it does not yet avoid that enumeration.

## Replay

```bash
sage -python elkies-k3/scripts/benchmark_inverse_ade_adaptive_backend.sage
sage -python elkies-k3/scripts/benchmark_inverse_ade_adaptive_backend.sage --check
```

The generated certificate is
[`../artifacts/generated-results/elkies-k3-inverse-ade-adaptive-backend-v1.json`](../artifacts/generated-results/elkies-k3-inverse-ade-adaptive-backend-v1.json).

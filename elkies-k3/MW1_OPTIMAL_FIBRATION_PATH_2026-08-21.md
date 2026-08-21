# An MW1 fibration on the discriminant-948 generic lattice (2026-08-21)

## Status

An exact elliptic-neighbor calculation has reached Mordell--Weil rank one on
the rank-19, discriminant-948 lattice used in the original search.  Rank zero
is impossible for that lattice, so rank one is its smallest possible
endpoint, not merely the smallest endpoint seen in a bounded search.

This is no longer the final lattice of the reconstructed rational surface.
The exact third section `S` raises that specialization to Picard rank 20 and
discriminant 43.  Consequently this note is retained as a correctly scoped
generic-lattice computation; the active reconstruction search uses
[`data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt`](data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt)
instead.

The currently preferred exact transition is

```text
canonical Kumar fibration: MW2, E7 + E8
  -- q=4, 2*2 --> MW1, E8 + D7 + A1.
```

This Kumar fibration is the symmetry-selected H2 comparison frame on the same
K3, not the corrected source polarization H3.  Its exact positive frame is
the height lattice `diag(4,237/2)` described in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).  The stable integral
NS identification with the recovered rank-17 lattice is proved there by the
local-genus calculation and Nikulin uniqueness.  What is not yet supplied is
an explicit integral transport or a chain of neighbor witnesses between the
rank-17 frame and this Kumar frame.  Thus the displayed `q=4` step is an exact
geometrically anchored improvement, but it must not be advertised as a
completed extension of the existing explicit `q=25,4,4,4` rank-17 chain.

## Exact neighbor and terminal data

In the basis of
[`data/fibrations/kumar_e7e8_mw2_frame_2.txt`](data/fibrations/kumar_e7e8_mw2_frame_2.txt),
the neighbor has

```text
q = 4,  (a,b) = (2,2),
v = (0,0,0,0,-1,-2,-1,0,0,0,0,0,0,0,0,-1,0).
```

The exact child is
[`data/fibrations/mw1_e8_d7_a1_frame.txt`](data/fibrations/mw1_e8_d7_a1_frame.txt).
Its invariants are

```text
frame determinant = 948
root rank/count/determinant = 16/326/8
root components = E8 + D7 + A1
root-plus-orthogonal index = 2
MW height Gram = [237/2].
```

The expected fiber presentation is

```text
II* + I3* + I2 + 3 I1.
```

The saturated glue gives the unique generator profile

```text
(E8 class, D7 class, A1 class; P.O) = (0, vector, 1; 58).
```

Indeed the local self-correction is `1+1/2=3/2`, and Shioda's
formula gives

```text
237/2 = 4 + 2*58 - 3/2.
```

There are no pairwise section gates at rank one, but the pole order `58` is
large.  Therefore this node is optimal for MW rank while the semistable
`A5+A4+2A3` MW2 node can still be the better equation-reconstruction chart.
Rank and reconstruction sparsity are different optimization objectives.

## Why MW0 cannot occur

Suppose an elliptic fibration on this K3 had MW rank zero, and let `R` be its
rank-17 ADE root lattice.  Shioda's determinant formula would give

```text
det(R) = |disc(NS)| * |MW_tors|^2 = 948 * |MW_tors|^2.
```

In particular `79` would divide `det(R)`.  This is impossible: an
irreducible rank-at-most-17 ADE lattice has determinant `n+1 <= 18` for
`A_n`, determinant `4` for `D_n`, and determinant `3,2,1` for `E6,E7,E8`.
No determinant of a direct sum of such factors is divisible by `79`.
Consequently every elliptic fibration on this NS has MW rank at least one.

[`scripts/verify_mw0_obstruction.py`](scripts/verify_mw0_obstruction.py)
independently enumerates all rank-17 ADE determinant products and checks the
prime obstruction.

## Lower-pole abstract alternative

The third genus-compatible Kumar frame has a sampled `q=9` MW1 neighbor with

```text
root system = E8 + A6 + A2,
height = 316/7,
profile = (A6 label 1 or 6, A2 label 0; P.O=21).
```

This is a substantially smaller section pole than `58`.  It is not selected
as the geometric path because the source involution and quotient labels pick
the second Kumar frame, not the third.  Nikulin uniqueness still proves that
the third stable lattice is abstractly the same NS; an explicit transport is
required before this lower-pole fibration can be used geometrically over the
current rational model.

## Bounded widened search

The discovery beam included the four canonical MW2 endpoints, alternate
representatives of the semistable branches, and all three genus-compatible
Kumar frames.  It searched proper-factor neighbors with `4 <= q <= 16`, six
restarted capped enumerations, and three retained representatives per root
invariant.  It retained 783 frames and found the MW1 systems

```text
E8+D7+A1, E7+D9, E7+D6+A3, E7+A9, and E8+A6+A2.
```

A second bounded round continued representatives of all five types.  A
separate descendant beam from the explicitly chained MW2 endpoints has not
yet found a rank-one child.  These negative observations are bounded search
results, not nonexistence theorems; only the determinant obstruction to MW0
is a proof of optimal rank.

## Reproduction

Replay the exact selected neighbor and terminal glue with

```bash
sage elkies-k3/scripts/verify_fibration_neighbor.sage \
  --parent elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_2.txt \
  --child elkies-k3/data/fibrations/mw1_e8_d7_a1_frame.txt \
  --q 4 --a 2 --b 2 \
  --v '0,0,0,0,-1,-2,-1,0,0,0,0,0,0,0,0,-1,0'

sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame elkies-k3/data/fibrations/mw1_e8_d7_a1_frame.txt \
  --name mw1-e8-d7-a1

sage elkies-k3/scripts/recover_mw1_e8_d7_a1_glue.sage
python3 elkies-k3/scripts/verify_mw0_obstruction.py
```

The bounded beam commands are recorded verbatim in
`artifacts/generated-results/elkies-k3-mw1-beam-round*/commands.txt`; the
input manifests alongside those directories pin every parent frame.

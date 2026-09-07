# Lower-q fibration path search (2026-08-20)

## Status

This is a **bounded exact lattice computation**, not an exhaustive
classification of elliptic fibrations and not a proof that any displayed path
is minimal.  Every retained frame, root system, height Gram, and neighbor
witness below is checked with exact integer or rational arithmetic.  The norm
shell traversal is capped and order-dependent.

The search answers the immediate question behind the old choice of `q=90`:
that choice was not minimal even inside a modest rerun.  The selected complete
replacement path is

```text
MW17
  -- q=25, 5*5 --> MW7, A3 + A1^7
  -- q=4,  2*2 --> MW4, D4 + A3 + A2^2 + A1^2
  -- q=4,  2*2 --> MW3, D6 + A5 + A3.
```

A subsequent widened beam has now continued an alternate child of the same
MW4 node to an exact MW2 frame.  See
[`MW2_FIBRATION_PATH_2026-08-21.md`](MW2_FIBRATION_PATH_2026-08-21.md) for the
new preferred `q=25,4,4,4` path and its exact witnesses.

The terminal frame is
[`data/fibrations/mw3_d6_a5_a3_frame.txt`](data/fibrations/mw3_d6_a5_a3_frame.txt).
The two intermediate frames are
[`data/fibrations/q25_mw7_frame.txt`](data/fibrations/q25_mw7_frame.txt) and
[`data/fibrations/q25_mw4_frame.txt`](data/fibrations/q25_mw4_frame.txt).

## Exact witnesses and invariants

Coordinates in each row refer to the positive frame Gram in the preceding
row.  With `NS = U + (-F)`, the isotropic class is `f=(a,b,v)` and satisfies
`a*b=q(v)=v*F*v/2`.

| step | `q=a*b` | `v` | root data `(rank,count,det)` | ADE |
|---|---:|---|---|---|
| MW17 -> MW7 | `25=5*5` | `(-1,0,-4,3,0,0,0,0,0,-1,1,0,0,0,-3,0,0)` | `(10,26,512)` | `A3 + A1^7` |
| MW7 -> MW4 | `4=2*2` | `(-1,-2,1,0,1,1,2,-3,0,-2,0,1,0,0,-1,0,0)` | `(13,52,576)` | `D4 + A3 + A2^2 + A1^2` |
| MW4 -> MW3 | `4=2*2` | `(-1,-1,2,0,-2,0,-1,0,-1,0,0,-1,-1,1,1,0,0)` | `(14,102,96)` | `D6 + A5 + A3` |

The exact saturated Mordell--Weil height lattices, in reduced bases, are

```text
MW7: (1/4) *
  [ 4  0 -2  2  0  2  0]
  [ 0  4  0  2  0  0  0]
  [-2  0  5  0 -1  1 -2]
  [ 2  2  0  6  0  2  0]
  [ 0  0 -1  0  5  1  2]
  [ 2  0  1  2  1  7 -2]
  [ 0  0 -2  0  2 -2  8]
det = 237/128

MW4: (1/12) *
  [8 2  2  2]
  [2 13 6  0]
  [2 6 22 -8]
  [2 0 -8 22]
det = 79/48

MW3: (1/12) *
  [23  8  1]
  [ 8 26 -8]
  [ 1 -8 35]
det = 79/8.
```

The expected semistable/additive fiber presentations are respectively

```text
I4 + 7 I2 + 6 I1,
I0* + I4 + 2 I3 + 2 I2 + 4 I1,
I2* + I6 + I4 + 6 I1.
```

For the terminal node, placing `I2*` at infinity suggests the short
Weierstrass bounds `deg(A)<=6`, `deg(B)<=9`.  This makes the node a serious
reconstruction candidate: it has only three reducible fibers and a relatively
small height Gram.

After flipping the sign of the third reduced-basis vector, use

```text
12*Gram = [23  8 -1]
          [ 8 26  8]
          [-1  8 35].
```

A canonical profile choice `(D6 class, I6 label, I4 label; P.O)` is

```text
P1 = (0, 2, 1; 0)
P2 = (v, 5, 0; 0)
P3 = (v, 2, 1; 1),
```

where `v` is the D6 vector discriminant class of inverse-Cartan norm one.
All three section intersections are one.  Exact Shioda correction replay in
[`scripts/verify_mw3_d6_a5_a3_profiles.py`](scripts/verify_mw3_d6_a5_a3_profiles.py)
recovers the signed Gram.  Thus the reconstruction uses two polynomial
sections, one simple-pole section, and minimal pair gates.  On the presently
known profile data this is materially cleaner than both the active A10 branch
and the A13 branch.

A further proper-factor scan through `q=12` found no MW2 child.  It did find a
`q=4` MW3 neighbor with ADE `D7+A5+A1^2`, expected fibers
`I3*+I6+2 I2+5 I1`, and an isometric reduced height Gram
`(1/12)*[23,-8,-1; -8,26,-8; -1,-8,35]`.  The selected `D6+A5+A3` node keeps
the same degree bounds with one fewer reducible fiber, so it remains the
preferred node on this path.

## Other lower-q entry nodes

The completed part of the first-stage rerun covered every composite
`4 <= q <= 46`, using one deterministic vector followed by six permuted-basis
restarts capped at 750 new sign-pairs each.  It found direct MW7 nodes already
at

| `q` | root data | ADE |
|---:|---|---|
| 16 | `(10,26,432)` | `A2^3 + A1^4` |
| 24 | `(10,28,384)` | `A3 + A2 + A1^5` |
| 25 | `(10,26,512)` | `A3 + A1^7` |
| 28 | `(10,28,324)` | `A2^4 + A1^2` |
| 32 | `(10,26,432)` | `A2^3 + A1^4` |

The `q=16` entry also has a sampled `q=4` MW4 child, with root data
`(13,52,480)` and ADE `A4+A3^2+A2+A1`.  Thus the lower-q phenomenon is not
isolated to the selected `q=25` path.

## Continuation of the old MW3 node

The old `A10+A2+A1^2` MW3 frame has a proper-factor `q=4` neighbor with

```text
ADE = A13 + A1
root data = (14,184,28)
expected fibers = I14 + I2 + 8 I1
height Gram = (1/14) * [29 -5 -2; -5 60 -18; -2 -18 60]
det = 237/7.
```

Its exact frame is
[`data/fibrations/mw3_a13_a1_frame.txt`](data/fibrations/mw3_a13_a1_frame.txt).
A bounded proper-factor continuation scan through `q=16` found many MW3
neighbors but no MW2 node.  A supplied component-profile calculation gives a
convenient signed basis with profiles/O
`(6,1)/1`, `(5,1)/2`, `(12,0)/1` and pairwise section intersections
`2,3,3`; the unflipped third profile `(2,0)` changes these to `2,3,5`.
Those profile claims are not yet independently checked by a repository
verifier.  The high-pole second section and nontrivial pair gates mean this
branch is not automatically simpler than the active A10 branch.

The two retained `q=8`, `(a,b)=(2,4)` continuations from this frame have now
been audited in the effective chamber. Their 39- and 40-reflection reductions
coincide at the primitive nef class

```text
(2,2,-32,23,22,15,-8,-4,-2,-23,6,-4,11,11,-18,17,-15,15,23).
```

The class has old degree two, `D.O=0`, and passes exact section and bisection
wall proofs. Both children nevertheless have `A10+A2+2A1` roots and the old
A10 MW Gram up to a basis sign, so they are genuine nef return loops rather
than reductions to a better frame. Replay the chamber certificate with
[`scripts/analyze_a13_q8_neighbors.sage`](scripts/analyze_a13_q8_neighbors.sage).

## Reproduction

The generalized discovery driver now accepts an arbitrary frame, suppresses
unit-factor rediscoveries, can cap retained representatives per root invariant,
and can write exact frame files:

```bash
sage elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame elkies-k3/data/lattice/rank17_gram.txt \
  --min-qnorm 25 --max-qnorm 25 \
  --enum-baseline-cap 1 --enum-restarts 7 --enum-cap 750 \
  --enum-seed 1729 --proper-factors-only --report 3 \
  --out artifacts/generated-results/elkies-k3-lower-q-path-search/first-stage/q25.txt \
  --frames-dir artifacts/generated-results/elkies-k3-lower-q-path-search/first-stage/q25-frames
```

The three exact neighbor replays use
[`scripts/verify_fibration_neighbor.sage`](scripts/verify_fibration_neighbor.sage).
The root components, saturated glue, and height lattices are replayed by
[`scripts/analyze_mw3_branch.sage`](scripts/analyze_mw3_branch.sage), with
full outputs in `artifacts/generated-results/elkies-k3-lower-q-path-search/`.

No search recorded here is exhaustive over a full norm shell.  In particular,
failure to find an MW2 child is only a bounded negative computation.

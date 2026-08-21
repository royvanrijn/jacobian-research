# Optimal MW1 path on the discriminant-43 K3 (2026-08-21)

## Result

The corrected Picard-rank-20 Neron--Severi lattice admits a two-neighbor
path from the reconstructed rational fibration to Mordell--Weil rank one:

```text
E6 + D4 + 2A2 + A1, MW3
  -- q=8, (a,b)=(2,4) --> A7 + A4 + A3 + A2, MW2
  -- q=9, (a,b)=(3,3) --> A12 + A3 + A2, MW1.
```

MW rank zero is impossible on this discriminant-43 lattice, so the endpoint
is optimal rather than merely the best result in the bounded search.

This supersedes the discovered three-step route through two MW5 frames as
the preferred lattice path.  It also supersedes the discriminant-948 MW1
nodes as targets for the reconstructed rational specialization: those nodes
remain valid only on the earlier generic rank-19 lattice.

## Exact neighbor witnesses

Write the Neron--Severi lattice at each stage as `U + (-M)`, where `M` is the
positive rank-18 frame.  A witness `(q,a,b,v)` means

```text
a*b = q,   v^t M v = 2q,
```

so `(a,b,v)` is a primitive isotropic vector in `U + (-M)`.

The first transition starts from
[`data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt`](data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt)
and has

```text
q = 8,  (a,b) = (2,4),
v = (-1,0,0,0,0,0,0,0,0,-2,0,0,0,0,-1,0,0,-1).
```

Its exact child is
[`data/fibrations/picard20_mw2_a7_a4_a3_a2_frame.txt`](data/fibrations/picard20_mw2_a7_a4_a3_a2_frame.txt),
with

```text
frame det = 43,
root rank/count/det = 16/94/480,
ADE = A7 + A4 + A3 + A2,
MW rank = 2,
MW height Gram = (1/120) * [34  6]
                               [ 6 39].
```

The second transition has

```text
q = 9,  (a,b) = (3,3),
v = (0,-1,-2,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0).
```

Its exact child is
[`data/fibrations/picard20_mw1_a12_a3_a2_frame.txt`](data/fibrations/picard20_mw1_a12_a3_a2_frame.txt),
with

```text
frame det = 43,
root rank/count/det = 17/174/156,
ADE = A12 + A3 + A2,
MW rank = 1,
MW height = 43/156,
MW torsion = 0.
```

Composing the two changes of basis gives the full integral NS transport
[`data/fibrations/picard20_mw1_a12_a3_a2_ns_transport.txt`](data/fibrations/picard20_mw1_a12_a3_a2_ns_transport.txt).
It has determinant one and carries the initial NS Gram literally to the
terminal NS Gram.  Its first row is the terminal fiber in the initial basis:

```text
(a,b) = (26,48),
v = (-12,0,0,0,-1,-2,0,0,-1,-25,0,0,0,0,-12,0,3,-13),
q = a*b = 1248.
```

Thus the chain can formally be compressed to one `q=1248` splitting, but
the two small steps `q=8,9` are far better for geometric reconstruction.

## Reconstruction profile

The terminal root data admit an all-multiplicative presentation.  Up to
inversion of the cyclic component groups, its primitive free generator has
profile

```text
(A12 label, A3 label, A2 label; P.O)
  = (3 or 10, 1 or 3, 1 or 2; 0).
```

The local self-correction is

```text
30/13 + 3/4 + 2/3 = 581/156,
```

and Shioda's formula gives

```text
43/156 = 4 - 581/156.
```

Thus the sole free section is polynomial in a compatible Weierstrass chart
with this presentation, there are no pairwise free-section gates, and the
expected fiber list is

```text
I13 + I4 + I3 + 4 I1.
```

This makes the MW1 endpoint a better reconstruction target than the earlier
MW1 `A9+D6+A2` endpoint: the selected root system has a semistable candidate
presentation, trivial torsion, and a generator disjoint from the zero
section.  The lattice alone does not distinguish an `I3` fiber from the
additive `IV` realization of `A2`, or prove that all residual fibers are
`I1`; that requires the terminal Weierstrass model.

## Optimality

If an elliptic fibration on this K3 had MW rank zero, its reducible-fiber
root lattice would have rank 18 and Shioda's determinant formula would force
its determinant to be divisible by `43`.  No rank-at-most-18 irreducible ADE
factor has determinant divisible by `43`: `A_n` has determinant `n+1 <= 19`,
`D_n` has determinant `4`, and `E6,E7,E8` have determinants `3,2,1`.
Therefore MW0 cannot occur and MW1 is the sharp lower bound.

## Reproduction and status

Replay the exact NS extension, the complete neighbor path, terminal glue, and
the independent ADE determinant obstruction with

```bash
sage elkies-k3/scripts/verify_picard20_ns_extension.sage
sage elkies-k3/scripts/verify_picard20_mw1_path.sage
sage elkies-k3/scripts/recover_picard20_mw1_a12_a3_a2_glue.sage
python3 elkies-k3/scripts/verify_mw0_obstruction.py
```

The discovery searches are retained in
`artifacts/generated-results/elkies-k3-picard20-neighbor-fast*` and
`artifacts/generated-results/elkies-k3-picard20-mw2-beam-round1/`.  They are
bounded searches; their failure to find other endpoints is not a
nonexistence theorem.  The pinned matrices and replay scripts above are exact
integral lattice certificates.

For comparison, a wider direct scan through `q=20` tested `387035` isotropic
presentations and retained `519` frame types.  It found no MW1 child.  Its
strongest new entrance was a `q=4` MW2 frame with roots `E7+D6+A3` and
root rank/count/determinant `16/198/32`.  A focused proper-factor child scan
of that frame through `q=16` tested another `38737` presentations and retained
`301` frame types, again with no MW1.  These are bounded negative results,
not proofs; they only support selecting the exact `q=8,9` route above.

The two discovery commands were

```bash
sage elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame elkies-k3/data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt \
  --min-qnorm 4 --max-qnorm 16 \
  --enum-baseline-cap 250 --enum-restarts 4 --enum-cap 250 \
  --enum-seed 4302026 --proper-factors-only --one-factor-order \
  --per-root-data-cap 3 --quiet-candidates --root-method pari --report 120 \
  --out artifacts/generated-results/elkies-k3-picard20-neighbor-fast.txt \
  --frames-dir artifacts/generated-results/elkies-k3-picard20-neighbor-fast-frames

python3 elkies-k3/scripts/run_fibration_beam.py \
  --input artifacts/generated-results/elkies-k3-picard20-mw2-beam-input.tsv \
  --out-dir artifacts/generated-results/elkies-k3-picard20-mw2-beam-round1 \
  --qmin 4 --qmax 25 --enum-baseline-cap 600 \
  --enum-restarts 6 --enum-cap 600 --enum-seed 4320260821 \
  --per-root-data-cap 4 --report 160 --workers 3

sage elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame elkies-k3/data/fibrations/picard20_e6_d4_a2a2_a1_mw3_frame.txt \
  --min-qnorm 2 --max-qnorm 20 \
  --enum-baseline-cap 1000 --enum-restarts 6 --enum-cap 1000 \
  --enum-seed 20260825 --one-factor-order --per-root-data-cap 5 \
  --quiet-candidates --root-method pari --report 200 \
  --out artifacts/generated-results/elkies-k3-picard20-neighbor-round1.txt \
  --frames-dir artifacts/generated-results/elkies-k3-picard20-neighbor-round1-frames

sage elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame artifacts/generated-results/elkies-k3-picard20-neighbor-round1-frames/frame-001.txt \
  --min-qnorm 4 --max-qnorm 16 \
  --enum-baseline-cap 600 --enum-restarts 6 --enum-cap 600 \
  --enum-seed 4320260822 --proper-factors-only --one-factor-order \
  --per-root-data-cap 5 --quiet-candidates --root-method pari --report 120 \
  --out artifacts/generated-results/elkies-k3-picard20-dense-mw2-child.txt \
  --frames-dir artifacts/generated-results/elkies-k3-picard20-dense-mw2-child-frames
```

This result does not yet execute the two inverse neighbors on the explicit
Weierstrass equation.  Turning the lattice path into rational parameter
changes and a terminal Weierstrass model remains a separate geometric step.

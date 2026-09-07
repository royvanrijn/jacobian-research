# Exact reverse transport from pinned R17 to H3

Status: exact integral lattice/marking certificate, 2026-08-23.

This note retraces the selected H3 corridor in the direction in which the
original reconstruction problem was posed:

```text
pinned recovered R17
  ~= H3-13 rootless/MW17
  <-q6-- A1/MW16
  <-q4-- 2A1/MW15
  <-q4-- 3A1/MW14
  <-q4-- 4A1/MW13
  <-q4-- 5A1/MW12
  <-q4-- A3+2A2/MW10
  <-q4-- 3A3/MW8
  <-q4-- 2A5/MW7
  <-q8-- A11/MW6
  <-q6-- D12/MW5
  <-q24- D13/MW4
  <-q8-- E8+E6/MW3
  <-q6-- H3 E7+E8/MW2.
```

The symbol `~=` in the first line is now an explicit determinant-one
positive-frame isometry, not merely equality of determinant, genus, ADE type,
or theta counts. The remaining arrows are inverted determinant-one
19-dimensional Neron--Severi basis changes. Thus the composition gives an
explicit integral matrix `R` satisfying

```text
R * (U + -rank17_gram) * R^t = U + -H3_frame.
```

This closes the previously missing identification between the unnamed
rootless endpoint of the forward H3 corridor and the pinned recovered
17-by-17 Mordell--Weil lattice.

## 1. The objects that had been conflated

There are three distinct layers.

1. [`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt) is the
   recovered positive-definite rank-17 Mordell--Weil Gram `P`.
2. The old forward replay ended at a rootless positive frame `F` stored only
   as `final_frame` in
   [`../artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json`](../artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json).
   It checked `det(F)=948` and rootlessness but did not compare `F` with `P`.
3. At q8, the first dominant D13 hit used by the pinned lattice corridor and
   the component-nef D13 representative used by the equation compiler are
   different marked elliptic fibrations. Their stored 17-by-17 positive
   frames are not equal.

For the endpoint, the pinned matrix
[`data/fibrations/h3_rootless_mw17_to_pinned_rank17_isometry.txt`](data/fibrations/h3_rootless_mw17_to_pinned_rank17_isometry.txt)
is a matrix `C` with

```text
det(C)=1,                 C^t * P * C = F.
```

This is the missing positive-frame certificate. It was found by exact PARI
`qfisom` and is replayed by direct matrix multiplication, so future checks do
not depend on rerunning an isometry search or on PARI returning the same
representative.

For q8, the new artifact retains an exact full-NS bridge from the dominant
D13 marking to the component-nef D13 marking. That bridge does **not** fix
the standard `U`: it moves the fiber and isotropic mate. Consequently the
two D13 records may be compared only after transporting the fiber, zero,
components, and section lifts through that bridge. Equality of the abstract
`D13/MW4` label is not enough.

## 2. What is retained

The forward artifact compressed each selected neighbor mainly to a witness,
child invariants, height Gram, and one final composite. The reverse ledger
[`../artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json`](../artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json)
retains all fourteen stages. For every stage it stores:

- the complete 17-by-17 positive frame;
- the incoming 19-by-19 neighbor basis and its integral inverse;
- the selected `q`, factor order, old-fiber degree, orbit, witness, and full
  isotropic fiber class when available;
- the stage basis in H3 NS coordinates and the H3 basis in stage coordinates;
- the stage basis in pinned-R17 NS coordinates and the pinned-R17 basis in
  stage coordinates.

The apparently redundant matrices are intentional. They make a later
divisor, component, zero section, or Mordell--Weil lift transportable in
either direction without reconstructing discarded root adaptations or
guessing a convention.

The ledger also stores both q8 representatives and their two inverse NS
bridges. This keeps the equation/compiler marking separate from the selected
lattice-corridor marking.

## 3. Replay

From the repository root run

```bash
sage -python elkies-k3/scripts/verify_rank17_to_h3_reverse_transport.sage
```

Expected terminal line:

```text
R17H3|endpoint_positive_isometry=PASS|det=1|stages=14|reverse_ns_transport=PASS|artifact=artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json|status=PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT
```

Pinned hashes for this certificate are

```text
99d5895ca61caf9f340d297cda9b036222ab0baad0b0a94cdf0fdaa36ba8e358  elkies-k3/data/fibrations/h3_rootless_mw17_to_pinned_rank17_isometry.txt
db9518ee9ba5ffb520898242cbff06894900ea1fea2908476e40433a212af4d2  artifacts/generated-results/elkies-k3-rank17-to-h3-reverse-transport.json
```

The checker also pins and verifies the hashes of the H3 q6/q8 entrance
artifact, the old eleven-step suffix artifact, the H3 frame, the recovered
rank17 Gram, and every selected neighbor-search artifact consumed by the
suffix.

## 4. Exact boundary

This proves the following lattice statement:

> The selected fourteen-stage marking ledger gives a lossless integral
> transport from the pinned recovered rootless rank-17 fibration back to the
> H3 `E7+E8/MW2` Neron--Severi frame.

It does not execute the eleven characteristic-zero equation pencils after
D13, transport explicit rational section functions through them, or identify
curve 273 as a specialization. It also does not prove that the selected
corridor is shortest or cheapest.

In particular, this reverse certificate does not make the active q24
equation problem disappear. It changes how that problem should be handled:
all future equation work should carry the full marked NS state alongside the
Riemann--Roch data, and any passage between the component-nef q8 equation
marking and the dominant D13 lattice marking must apply the recorded bridge.

## 5. Information-preserving protocol for later hops

For each future equation-level neighbor, retain at least:

1. the parent and child full NS Grams and both inverse integral basis changes;
2. fiber, isotropic mate, zero section, simple and affine components, and all
   section lifts in both parent and child coordinates;
3. every Weyl reflection/root adaptation as an ordered matrix or word;
4. divisor functions and their pole modules before quotienting or choosing a
   convenient zero;
5. translations, 2-cover multipliers, and cleared denominators as explicit
   maps, not prose conventions;
6. a qualified identifier for every representative (`dominant`,
   `component-nef`, `equation-zero`, and so on).

ADE type, Mordell--Weil rank, determinant, height Gram, and a successful
bounded search are useful summaries, but none of them is a replacement for
this marked transport data.
<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 2c6a2a36699933ab -->

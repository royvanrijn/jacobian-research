# Exact profile audit of the semistable MW2 endpoints (2026-08-21)

## Status

This is an **exact frame-glue and Shioda-height computation** for the three
all-multiplicative MW2 endpoints found by the bounded beam.  The discovery
beam was not exhaustive, but the component profiles, section intersections,
and basis optimizations below are exact for the three pinned frames.

The audit gives a clear semistable reconstruction order:

```text
A5+A4+2A3  >  A6+A4+A3+A2  >>  A9+A5+A1.
```

The first endpoint has an MW basis consisting of two polynomial sections.

## Optimized profiles

Component labels are ordered as the displayed `A` factors and are defined up
to independent diagram inversion.  Equal factors may also be swapped.

### 1. `A5+A4+2A3`

Expected fibers:

```text
I6 + I5 + 2 I4 + 5 I1.
```

A convenient optimized basis has

```text
60*Gram = [67  40]
          [40 130],       det(Gram)=79/40,

P1 = (2,4,0,1; P1.O=0)
P2 = (5,0,2,0; P2.O=0)
P1.P2 = 1.
```

Thus both generators are polynomial sections and their only pair gate is
minimal.  This basis is obtained from the initially reduced signed basis by
replacing its second vector by `P1-P2`.

### 2. `A6+A4+A3+A2`

Expected fibers:

```text
I7 + I5 + I4 + I3 + 5 I1.
```

An optimized signed basis has

```text
105*Gram = [81  39]
           [39 326],      det(Gram)=79/35,

P1 = (5,4,2,0; P1.O=0)
P2 = (5,1,0,2; P2.O=1)
P1.P2 = 1.
```

This gives one polynomial section, one simple-pole section, and one minimal
pair gate.

### 3. `A9+A5+A1`

Expected fibers:

```text
I10 + I6 + I2 + 6 I1.
```

The better of the two sign orientations has

```text
30*Gram = [117  24]
          [ 24 248],      det(Gram)=158/5,

P1 = (8,0,1; P1.O=1)
P2 = (1,5,0; P2.O=3)
P1.P2 = 5.
```

The large second pole and pair intersection make this endpoint substantially
less attractive despite its having only three reducible fibers.

## Exact basis optimality

For each frame, the verifier exhausts every MW vector satisfying

```text
P.O <= max(P1.O,P2.O)
```

in the displayed basis.  This is a complete finite enumeration because
Shioda's formula gives `height(P) <= 4+2(P.O)`.  It then checks every
unimodular pair and minimizes lexicographically

```text
(maximum pole, total poles, pair intersection).
```

The exact optima are

| endpoint | optimum | section interpretation |
|---|---|---|
| `A5+A4+2A3` | `(0,0,1)` | polynomial + polynomial; pair `1` |
| `A6+A4+A3+A2` | `(1,1,1)` | polynomial + simple pole; pair `1` |
| `A9+A5+A1` | `(3,4,5)` | poles `1,3`; pair `5` |

Consequently the two-polynomial basis for `A5+A4+2A3` is not merely an LLL
artifact: no section basis can improve any entry of its optimization tuple.

## Comparison with the reconstructed additive endpoint

The reconstructed `E6+D4+2A2+A1` endpoint has one polynomial section, one
simple-pole section, pair intersection `2`, and short-Weierstrass bounds
`deg(A)<=5`, `deg(B)<=8`.  A semistable K3 uses the generic bounds
`deg(A)<=8`, `deg(B)<=12`.

The `A5+A4+2A3` endpoint therefore trades seven additional ambient
Weierstrass coefficients for eliminating both section denominators and
reducing the pair gate to one.  On section-system complexity it is the most
promising semistable candidate and a credible alternative if its
reconstruction curve is more rational than the additive MW2 curve.

## Certificates and reproduction

The pinned frames are

- [`data/fibrations/mw2_a5_a4_a3a3_frame.txt`](data/fibrations/mw2_a5_a4_a3a3_frame.txt),
- [`data/fibrations/mw2_a6_a4_a3_a2_frame.txt`](data/fibrations/mw2_a6_a4_a3_a2_frame.txt), and
- [`data/fibrations/mw2_a9_a5_a1_frame.txt`](data/fibrations/mw2_a9_a5_a1_frame.txt).

Their SHA-256 hashes, in that order, are

```text
b15a6f98c8d03b1236768e387d4e2dd2b1b77a9b38fd1ca5517d7e40dbdb62b9
48eb8566bacc673302cdb85304812f5a865289c45e227f57bb9d18696bd8f5bf
475900679877885272bdcfd0ad2f08c7b7c9ad243cf2c33a72cb82f31f8e8616
```

Run

```bash
sage elkies-k3/scripts/recover_semistable_mw2_glue.sage
```

to recover the root decomposition, saturated MW lattice, both sign
orientations, optimized profiles, and complete bounded-by-height basis
certificate.  Each incoming neighbor also passes
[`scripts/verify_fibration_neighbor.sage`](scripts/verify_fibration_neighbor.sage)
with the witness recorded in the beam hit table.

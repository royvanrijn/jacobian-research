# Curve273 construction: H3 source and rational base

[MATH_STATUS.json](../../MATH_STATUS.json) is authoritative. This concise note
is the canonical source for `EC-K3-H3-SOURCE` and `EC-K3-H3-PTS`; the complete
dated investigation is [archived](../../archive/elliptic-curves/notes/ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.before-2026-09-13.md.txt).

## Exact H3 source family

The characteristic-zero normalization of the degree-21 H21/H92 component,
composed with the pinned H92 short-Weierstrass formulas over

```text
Y^2 = -27*X^6 + 198*X^4 - 171*X^2 + 576,
```

gives an exact H3 `E7+E8/MW2` source family. After inverting the published
linear-fractional `X` map and writing `Y0` for the normalized double-cover
coordinate, the H92 chart is

```text
r = (a + Y0)/2,  s = 2/(Y0 - a).
```

Exact function-field substitution gives the degree-21 component equation and

```text
v^2 = u^3 + (A1*tau^3 + A*tau^4)*u
          + (B1*tau^5 + B*tau^6 + B2*tau^7),
```

with `E7` and `E8` fibres and two individually rational Mordell--Weil
directions of height Gram `[[21/2,3],[3,46]]`. The published point
`(13/7,12048/343)` specializes to the pinned H92 point
`(-3621005/690947,158286/143585)` and its five exact coefficients.

The retained checker is
[`export_h3_level474_source_family.sage`](../../elkies-k3/scripts/export_h3_level474_source_family.sage);
the status-recorded artifact digest is
`3304f994a410cdfd1bd8bedc414efabda9e863b02ba6e22e773f3c2503e53e03`.

This proves the H3 source family. It does not prove the first q6 equation, a
downstream rootless MW17 family, or a specialization of such a descendant to
Curve273. The rootless-chart degree-24 `j`-recognition equation has no
rational root for Curve273; a shared construction would need a different
family, an isogeny-level explanation, or a construction certificate.

<!-- status-consumer: EC-K3-H3-SOURCE a4bb40c9c9d0ff09 -->

## All rational points on the H3 base

For the Q-isomorphic model

```text
C: y^2 = -3*x^6 + 22*x^4 - 19*x^2 + 64,
```

the current Sage 10.9 certificate uses its two degree-two elliptic quotients

```text
y^2=x^3+22*x^2+57*x+576,
y^2=x^3-19*x^2+1408*x-12288.
```

Both are exactly rank one, torsion-free, and have conductor `474`. The pinned
Bianchi--Padurariu bielliptic quadratic-Chabauty implementation at `p=11,41`
finds the three expected rational automorphism orbits and leaves 71,776
finite coordinate classes modulo `11^4*41^4`. The exact quotient
Mordell--Weil sieve eliminates every one by `q=1987`. Since `-3` is not a
square, there are no rational points at infinity. Hence, in published
coordinates,

```text
(X,Y) = (0,+-24), (+-1,+-24), (+-13/7,+-12048/343).
```

The generated [current certificate](../../artifacts/generated-results/elkies-k3-h3-level474-rational-points-qc-sage109.json)
has SHA-256 `ab4e535594984bbcf756fa6cdfddbc36b98cc944e44edef594ef70fcfe32ad5a`.
Replay it from `research/` with

```sh
sage elkies-k3/scripts/certify_h3_level474_rational_points_qc.sage --check
```

The checker source SHA-256 is
`92bf527506840bac0d11559e2ab4e2b22924137ea31d2930dd6d4fefb0c71f5e`.
It pins the upstream quadratic-Chabauty revision and its two Sage-10.9
point-at-infinity compatibility edits; it is an independent current proof,
not a replay of the historic Magma calculation. For an offline fresh checkout,
pass the single upstream file with `--upstream-source PATH`; its original
SHA-256 is checked before the compatibility patch is applied. Without that
argument, the checker first uses its hash-checked local cache and otherwise
downloads the exact pinned revision. It never accepts a changed upstream file.

The retained historical program
[`prove_h3_level474_rational_points.m`](../../elkies-k3/scripts/prove_h3_level474_rational_points.m)
records the complementary two-cover/elliptic-Chabauty route, with source digest
`6efd86a474664c0fda2b5487f4d04e70d8926762860969ee1b7883f1b1db19c5`.
Its original output remains unavailable. The earlier Sage
[`bounded quotient successor`](../../elkies-k3/scripts/sieve_h3_level474_rational_points_sage109_replay.sage)
remains a finite cross-check through `|n| <= 1000000`; neither historical
provenance limitation weakens the current certificate.

<!-- status-consumer: EC-K3-H3-PTS a98696b3defadd8f -->

## Certified H3 lattice transport

The selected H3 `E7+E8/MW2` corridor passes through `D13/MW4` and eleven
stored primitive integral `U`-neighbours with q-sequence
`(24,6,8,4,4,4,4,4,4,4,6)`, ending in a rootless MW17 frame. Exact chamber
checks prove the displayed factor presentations nef. A pinned
determinant-one positive-frame isometry identifies that endpoint with the
recovered rank-17 Gram, and inverting the full composite gives the explicit
integral reverse Neron--Severi transport. The fourteen-stage ledger separates
the lattice-corridor `D13` marking from the unequal component-nef `D13`
marking used by the q8 equation route; their shared label does not make the
marked fibrations interchangeable.

This is a lattice/marking certificate only. The equation-level pencils after
`D13/MW4`, transported section functions, and specialization to Curve273 are
open. The exact verifier is
[`verify_rank17_to_h3_reverse_transport.sage`](../../elkies-k3/scripts/verify_rank17_to_h3_reverse_transport.sage).

<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 2c6a2a36699933ab -->

## Replay and scope boundary

The historical H3 source artifact and its three generated input records remain
unavailable, so its frozen hash lineage is source-level evidence rather than a
portable historical replay package. A separately named Sage 10.9 successor
reconstructs the exact factor from the pinned public H21 input, normalizes it,
checks the non-CM anchor, and reruns the same source-family composition:
[Sage 10.9 successor checker](../../elkies-k3/scripts/export_h3_level474_source_family_sage109_replay.sage)
and its [generated artifact](../../artifacts/generated-results/elkies-k3-h3-level474-source-family-sage109-replay.json).
It has its own hashes and does not retag or replace the frozen historical
artifact.

The historical point proof requires the stated Magma two-cover and
elliptic-Chabauty runtime, and its original result is unavailable. The current
Sage 10.9 certificate above is the retained independent proof of the rational
base; it does not recreate or retag the historical Magma result.

For Curve273 itself, the independent rank lower bound is documented in
[ICARM_CURVE273_RANK30.md](ICARM_CURVE273_RANK30.md). It proves only
`rank E(Q) >= 30`; exact rank and a 31st independent point remain open.

The [H3 backtrack](../../elkies-k3/KUMAR_E7E8_BACKTRACK.md) and archived diary
retain the route analysis, failed approaches, bounded calculations, and
historical commands. They do not authorize a new construction campaign.

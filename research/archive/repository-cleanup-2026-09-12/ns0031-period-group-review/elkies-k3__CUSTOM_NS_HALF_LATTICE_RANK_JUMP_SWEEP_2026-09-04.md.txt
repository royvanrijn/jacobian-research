# Custom Picard-19 NS half-lattice rank-jump sweep

Date: 2026-09-04.

Status: **complete finite lattice sweep; arithmetic rank jumps not claimed**.

## Result

The sweep isolates the repository's custom rank-19 Néron--Severi labels
`NS0002` through `NS0048`; the published determinant-948 `NS0001` control is
not part of the target batch.  Of the 47 custom NS classes, 32 have at least
one attached rootless rank-17 frame whose integral Gram matrix is the
Mordell--Weil height lattice.  The remaining 15 have maximum catalogued MW
rank 15 or 16, so their rooted frame Gram cannot be substituted for a
Mordell--Weil Gram.

The exact eligible batch is

```text
32 custom NS classes
136 rootless MW17 frame classes
131,072 parity classes per frame
17,825,792 parity classes in total
```

For every parity class `c in M/2M`, the sweep computes its minimum integral
norm `mu_2(c)`.  The maximum spectrum across frames is

```text
max_c mu_2(c) = 12    on 135 frames
max_c mu_2(c) = 14    on   1 frame
```

The unique outlier is determinant-1494 `NS0021`, frame
`K3-2cacfeb9277cda96-F006`.  Its complete minimum-norm histogram is

```text
norm       0      4       6       8       10    12   14
classes    1   1045   21158   63577   44377   913    1
```

The unique norm-14 class has mask `116966 = 0x1c8e6` and one stored shortest
representative

```text
(2,-1,-3,0,0,1,1,1,2,2,-2,-1,0,2,-3,-1,1).
```

Corollary S4 of
[`SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md`](SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md)
therefore gives the exact discrete bound

```text
rho_2(M)^2 = 14/4 = 7/2,
rho(M) >= rho_2(M) = sqrt(7/2)
```

for that frame.  The published R17 reference has `max mu_2=12` and
`rho_2=sqrt(3)`.  Thus the custom frame has a strictly deeper midpoint hole
than the reference lattice.

## What the theorem test says about a search

The result is an old-subgroup exclusion certificate, not a positive
rank-jump signal.  It nevertheless gives a sharp, cheap future schedule:
if `NS0021` passes its arithmetic marking and equation gates, its unique
norm-14 class is the first half-lattice chart to compile.  It is both deeper
than every other custom class in this sweep and a singleton stratum.

Among the arithmetic-`UNKNOWN` custom NS rows, the next cheapest complete
deepest strata are

| NS | frame | determinant | `max mu_2` | deepest classes |
|---|---|---:|---:|---:|
| `NS0021` | `K3-2cacfeb9277cda96-F006` | 1494 | 14 | 1 |
| `NS0005` | `K3-29b863dcdefc0eff-F001` | 992 | 12 | 50 |
| `NS0022` | `K3-7a0bdbb6d96eb93d-F002` | 1018 | 12 | 86 |
| `NS0002` | `K3-333d8fc7f9c5b179-F001` | 1028 | 12 | 92 |
| `NS0011` | `K3-32f4cd1fe26b1dd3-F001` | 1056 | 12 | 100 |

This ordering is only the cost of exhausting the deepest S4 layer.  Neither
depth nor stratum size is treated as a success probability.  The other new
theorems remain correctly gated:

- S3 becomes usable only after an equation and a marked section center are
  available;
- S1 requires a specialization with preserved generic rank, certified new
  points, and a genuine residual Selmer calculation for its upper side;
- S5 requires actual exceptional-point Kummer classes in a declared number
  field and is explanatory, not prospective.

## Arithmetic gate

No custom row in this batch is currently authorized for equation-level work:
30 are `UNKNOWN`, while `NS0024` and `NS0031` are arithmetically excluded over
`QQ`.  In particular, their geometric spectra remain regression data only.
For `NS0021`, the next certified step is to embed its displayed split order
in `M_2(QQ)`, derive the exact congruence conditions and signature, and
compute the stable marking subgroup.  A rational noncuspidal non-CM lift is
still required before equation compilation.

## Certificate and replay

The machine-readable certificate is
[`../artifacts/generated-results/elkies-k3-custom-ns-rootless-half-lattice-sweep-v1.json`](../artifacts/generated-results/elkies-k3-custom-ns-rootless-half-lattice-sweep-v1.json),
with SHA-256
`c60fc5b3a369bd058bbb4fcfb93e0a395280c5ea068d4ed20101a51d459d3ddd`.
Every one of the `2^17` returned norms per frame is recomputed integrally.
Every deepest class and a deterministic stride through all other classes is
repeated with 256-bit MPFR Gram--Schmidt arithmetic.  The runner checkpoints
after each complete frame.

```bash
sage -python elkies-k3/scripts/sweep_custom_ns_half_lattice_depths.sage
sage -python elkies-k3/scripts/sweep_custom_ns_half_lattice_depths.sage --check
```

The sweep constructs no K3 equation or rational marking, performs no
specialization, finds no point, and computes no Selmer group.  It therefore
does not prove a rank jump or turn a bounded search miss into point absence.

<!-- status-consumer: EC-K3-CUSTOM-NS-HALF-LATTICE-SWEEP 9dc0e4d23f677392 -->

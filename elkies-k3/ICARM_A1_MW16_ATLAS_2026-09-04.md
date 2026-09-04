# ICARM norm-eight low-root atlas from the 11952 rootless frame

Date: 2026-09-04

Status: exact complete target screen and exact root-rank stratification for all
63,917 minimum-norm-eight, old-degree-two classes on
`norm12-orbit-11952`.  Curve 302 misses every stratum.  Curve 398 supplies the
requested positive control at priorities 16,875 and 63,669.

## Result

The curve-398 pipeline was made target-generic and replayed on every class.
An exact audit of the already-pinned minimum shells and singular-pencil
certificate then corrects the old uniform `A1/MW16` label.  If a parity class
has `m` minimum representatives up to sign, its pencil has `m` split members;
these exhaust the even part of the discriminant.  A complete norm-twelve shell
replay also gives every class a degree-one section.  Thus the actual strata are:

| classes | fibres | root lattice | geometric MW rank at rho=19 | curve 302 |
|---:|---|---|---:|---|
| 1,266 | `I2+22I1` | `A1` | 16 | exact miss |
| 8,410 | `2I2+20I1` | `2A1` | 15 | exact miss |
| 20,348 | `3I2+18I1` | `3A1` | 14 | exact miss |
| 21,405 | `4I2+16I1` | `4A1` | 13 | exact miss |
| 9,861 | `5I2+14I1` | `5A1` | 12 | exact miss |
| 2,280 | `6I2+12I1` | `6A1` | 11 | exact miss |
| 331 | `7I2+10I1` | `7A1` | 10 | exact miss |
| 16 | `8I2+8I1` | `8A1` | 9 | exact miss |

The sum is 63,917.  In particular, the run does not merely test the genuine
1,266-class `A1/MW16` stratum: it also completes the 8,410-class
`2A1/MW15` stratum and all subsequent strata through `8A1/MW9` inside this
norm-eight atlas.  This uses the full class ordering and does not reconstruct
cores from curve 302's known points.

The eleven-target screen comprises 703,087 target--class pairs:

| ICARM curve | rank lower bound | complete norm-eight outcome | hit priorities |
|---:|---:|---|---|
| 302 | 31 | miss in every stratum | — |
| 273 | 30 | miss in every stratum | — |
| 542 | 26 | hit in `A1/MW16` | 30,486 |
| 548 | 24 | three `A1/MW16` hits | 31,627; 54,835; 63,647 |
| 398 | 30 | two `A1/MW16` presentations | 16,875; 63,669 |
| 399 | 29 | miss in every stratum | — |
| 400 | 28 | two `A1/MW16` hits | 53,042; 62,992 |
| 403 | 28 | miss in every stratum | — |
| 401 | 27 | hit in `A1/MW16` | 57,487 |
| 402 | 27 | miss in every stratum | — |
| 10 | 24 | miss in every stratum | — |

All nine hit presentations lie in the genuine `m=1` stratum and compile with
`I2+22I1`, generic Mordell--Weil rank 16, and saturated height determinant
474.  The target parameter is rational and the specialized curve is
isomorphic over `QQ`, not merely a quadratic twist.  A later exact base-change
audit identifies the repeated labels on curves 398, 400, and 548 as
presentations of one fibration for each target.  Thus the nine labels give
exactly five fibration classes.  All nine coordinate presentations remain
useful for bounded parameter searches because their affine base changes do not
preserve projective-height boxes; their search responses are nested within the
five fibrations, not independent observations.  See
[`../elliptic-curves/notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md`](../elliptic-curves/notes/ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md).
<!-- status-consumer: EC-K3-ICARM-MW16-BLIND-LADDER c5b0b57ee01c5c23 -->

For curve 10, one class survives the modular chain, but its exact comparison
polynomial is irreducible of degree 24.  The other negative targets are
eliminated completely by modular no-root witnesses.  Hence curve 302's miss is
exact over `QQ` in each displayed stratum.

## Independent fixed-corridor control

The previously compiled fixed-corridor `2A1/MW15` pencil is separate from the
norm-eight atlas.  It excludes curve 302 modulo 1019 and curve 273 modulo
1009.  This is retained as an equation-level regression, not counted as an
additional complete layer.

## Operational supersession

This target-first atlas is calibration and regression material.  Its ICARM
equations, parameters, rank lower bounds, and target `j`-invariants do not
enter the target-free A1/MW16 parameter experiment, which samples new
parameters directly and measures exact specialization quotient gains.

## Certificates and replay

The root-stratum certificate is
[`../artifacts/generated-results/elkies-k3-icarm-11952-norm8-low-root-strata-v1.json`](../artifacts/generated-results/elkies-k3-icarm-11952-norm8-low-root-strata-v1.json),
and the target-by-stratum atlas is
[`../artifacts/generated-results/elkies-k3-icarm-11952-norm8-low-root-atlas-v2.json`](../artifacts/generated-results/elkies-k3-icarm-11952-norm8-low-root-atlas-v2.json).
The historical unstratified target ledger is preserved as a pinned input so
the corrected theorem does not discard any exclusion or hit evidence.

```bash
sage -python elkies-k3/scripts/certify_icarm_norm8_low_root_strata.sage
python3 elkies-k3/scripts/build_icarm_norm8_low_root_atlas.py

sage -python elkies-k3/scripts/certify_icarm_norm8_low_root_strata.sage --check
python3 elkies-k3/scripts/build_icarm_norm8_low_root_atlas.py --check
```

The first command streams the exact R17 shell through norm 12.  It verifies
all 63,917 priority rows, the complete split-member/even-discriminant
histograms, degree-one sections for every class, and one primitive-`U`/root
frame in each of the eight strata.  The second joins those strata to every
modular witness, exact survivor factorization, and compiled hit.

The original screen and independent fixed-corridor regression remain
replayable with:

```bash
python3 elkies-k3/scripts/run_icarm_norm8_a1_atlas.py \
  --curve-ids 302,273,542,548,399,400,403,401,402,10 --resume
python3 elkies-k3/scripts/build_icarm_a1_mw16_atlas.py --check
sage -python elkies-k3/scripts/screen_icarm_fixed_mw15_fibrations.sage --check
```

## Proof boundary

“Complete” here means the minimum-norm-eight, old-degree-two residual-chord
translation atlas on source chart 11952.  It is not a classification of all
elliptic fibrations on the determinant-948 K3: `A2/MW15`, other old degrees,
other trace norms, and other rootless source charts remain outside the
theorem.  Therefore this exact miss does not prove that curve 302 has no K3
parent.  ICARM ranks and specialization jumps remain lower bounds; no exact
target-rank claim is made.

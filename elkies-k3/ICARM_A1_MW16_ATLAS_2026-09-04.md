# ICARM A1/MW16 atlas from the 11952 rootless frame

Date: 2026-09-04

Status: exact complete screen for the declared 63,917-class A1/MW16 layer on
`norm12-orbit-11952`; nine target hits compiled exactly.  The subsequent MW15
test is exact but presently bounded to one fixed-corridor `2A1/MW15` pencil.

## Result

The curve-398 forensic pipeline was made target-generic and replayed in the
requested order:

```text
X948 rootless frame
  -> complete minimum-norm-eight translation class
  -> old-degree-two residual-chord A1/MW16 pencil
  -> projective modular target-j screen
  -> exact QQ factorization and twist-sensitive isomorphism
  -> exact equation and saturated generic-MW compile
```

There are 63,917 classes per target.  Eleven targets therefore give 703,087
target--fibration pairs.  The outcomes are:

| ICARM curve | rank lower bound | exact A1 outcome | compiled priority ranks | specialization jump lower bound |
|---:|---:|---|---|---:|
| 302 | 31 | complete-layer miss | — | — |
| 273 | 30 | complete-layer miss | — | — |
| 542 | 26 | hit | 30,486 | +10 |
| 548 | 24 | three hits | 31,627; 54,835; 63,647 | +8 |
| 398 | 30 | two hits | 16,875; 63,669 | +14 |
| 399 | 29 | complete-layer miss | — | — |
| 400 | 28 | two hits | 53,042; 62,992 | +12 |
| 403 | 28 | complete-layer miss | — | — |
| 401 | 27 | hit | 57,487 | +11 |
| 402 | 27 | complete-layer miss | — | — |
| 10 | 24 | complete-layer miss | — | — |

Every one of the nine hits compiles to a polynomial K3 pencil with
`I2+22I1`, generic Mordell--Weil rank 16, and saturated height determinant
474.  Each target parameter is rational and the specialized curve is
isomorphic over `QQ`, not merely a quadratic twist.  Curve 398's second
survivor also passes the unified compiler, so both of its norm-eight hits now
have independent equation and generic-lattice certificates.

For curve 10, one class survived the declared modular chain; its exact
comparison polynomial is irreducible of degree 24.  The other five negative
targets were eliminated completely by modular no-root witnesses.  Thus all
six negative outcomes are exact over `QQ` within this committed A1 layer.

## First MW15 extension

Because both priority targets 302 and 273 miss the complete A1 layer, the
screen was extended immediately to the exact source-identified fixed-corridor
`2A1/MW15` Jacobian.  The same projective target-j gate excludes curve 273
modulo 1009 and curve 302 modulo 1019.  These are exact no-rational-parameter
certificates for that pencil.

This MW15 result is deliberately not labelled an atlas: it contains one
certified `2A1/MW15` fibration and does not exhaust all `A2` or `2A1` frames.
The next mathematically honest step is a complete, explicitly bounded
enumeration of the next neighbor layer before any broader negative statement.

## Certificates and replay

The compact atlas certificate is
[`../artifacts/generated-results/elkies-k3-icarm-11952-norm8-a1-mw16-atlas-v1.json`](../artifacts/generated-results/elkies-k3-icarm-11952-norm8-a1-mw16-atlas-v1.json).
It hash-pins every modular, exact-factorization, and compiled input.  The
bounded MW15 certificate is
[`../artifacts/generated-results/elkies-k3-icarm-curve302-273-fixed-2a1-mw15-screen-v1.json`](../artifacts/generated-results/elkies-k3-icarm-curve302-273-fixed-2a1-mw15-screen-v1.json).

```bash
python3 elkies-k3/scripts/run_icarm_norm8_a1_atlas.py \
  --curve-ids 302,273,542,548,399,400,403,401,402,10 --resume

sage -python elkies-k3/scripts/compile_icarm_norm8_a1_survivors.sage \
  --curve-id 542 \
  --exact-survivors artifacts/generated-results/elkies-k3-icarm-curve542-11952-norm8-a1-exact-survivors-v1.json \
  --output artifacts/generated-results/elkies-k3-icarm-curve542-11952-norm8-a1-compiled-survivors-v1.json --check

python3 elkies-k3/scripts/build_icarm_a1_mw16_atlas.py --check

sage -python elkies-k3/scripts/screen_icarm_fixed_mw15_fibrations.sage --check
```

The target-generic compiler command is identical for curves 400, 401, 542,
and 548 after substituting the curve number.  Curve 398 uses its frozen exact
survivor ledger as the compiler input; the aggregate atlas pins the resulting
unified two-fibration compile certificate.

## Proof boundary

“Complete-layer miss” means no rational target parameter occurs among the
63,917 minimum-norm-eight, old-degree-two A1/MW16 classes on source chart
11952.  It says nothing about another rootless chart, a higher-norm divisor,
or an `A2/2A1` neighbor.  ICARM ranks and the displayed specialization jumps
are lower bounds from the pinned 573-curve snapshot; no exact target-rank
claim is made.

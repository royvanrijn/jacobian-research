# Curve 302: preserved strict-residual bootstrap chain

The canonical machine-checkable record is
[`curve302_residual_strict_bootstrap_chain_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_residual_strict_bootstrap_chain_v1.json).
It deliberately retains the adaptive path rather than replacing it with a
single retrospective account.

| Certified transition | strict direction | frozen degree-two shell |
| --- | --- | --- |
| 24 to 25 | 06 | rational norm 10, orbit 117420 |
| 25 to 26 | 04 | rational norm 10, orbit 64677 |
| 26 to 27 | 02 | rational norm 10, orbit 58145 |
| 27 to 28 | 01 | rational norm 10, orbit 106210 |
| 28 to 29 | 05 | rational norm 10, orbit 4761 |
| 29 to 30 | 03 | genus-one nearby norm 8, orbit 17845 |
| 30 to 31 | 07 | genus-one nearby norm 8, orbit 114326 |

The manifest also pins the prior certified calibration states (M_{17}),
(M_{19}), (M_{22}), and (M_{24}), every arm's source-hashed protocol,
centre/map/worker/replay checkpoints, exact new independent witness, and the
modulo 2, 3, and 5 certificates.  Every arm after (M_{24}) was
retrospectively calibrated but execution-blind; this is not evidence that the
selector was prospective, nor an exact-rank upper bound.

Replay with:

```bash
python3 elliptic-curves/cas/audit_curve302_residual_strict_bootstrap_chain.py --check
```

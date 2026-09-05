# Status source-hash reconciliation

Date: 2026-08-25.

## Scope

The repository status audit found eight `artifact_hash` values that no longer
matched their committed checker source. This maintenance change updates only
those source locks in `MATH_STATUS.json`. It does not refresh a generated
artifact, rerun Sage, alter a theorem scope, or promote a mathematical claim.

The audit and digest command were run from the repository root with the system
`python3` and GNU `sha256sum`:

```bash
python3 scripts/audit_status.py
sha256sum \
  elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_promoted_route.sage \
  elliptic-curves/cas/analyze_icarm_construction_fingerprints.py \
  elkies-k3/scripts/certify_h92_q24_a11_q8_equation_marking_qq.sage \
  elkies-k3/scripts/audit_h92_a5a5_q6o1307_physical_nef.sage \
  elkies-k3/scripts/certify_h92_a5a5_direct_physical_q10_promoted_route.sage \
  elkies-k3/scripts/audit_h92_a5a5_q4o230_effective_return_zero.sage \
  elkies-k3/scripts/certify_h92_first_q8_q4o11_promoted_route.sage \
  elkies-k3/scripts/certify_h92_d13_q4o11_promoted_route.sage
```

## Changed locks

| Status entry | Previous SHA-256 | Current SHA-256 |
|---|---|---|
| `EC-K3-H3-A11-R17-PHYSICAL-Q4O208-PROMOTED-ROUTE` | `9017715bfe81cd73c691db0d95fde7703f973cd9ebc574421049820b1ebd9430` | `ded9dd735c9ff1ed2c41f1907fca6d491ed2ebf7c966e4dc706ed3e34cb244c4` |
| `EC-ICARM-CONSTRUCTION-FINGERPRINTS` | `4da9e79e4142d1daf56554df982fde59a23ba8a96f74505261cc5f19abca393d` | `fdaa2de2aa8c312e131dd55760855bc1d701a851d6b2d506d04a6c93359b6b42` |
| `EC-K3-H3-A11-Q8-QQ-2A5` | `b9d46929f37de95c68910f6246ee07b32d7bdb9476b79e6629fda0cc24cb23f9` | `9a22977e01e0eb3c88b893313a7a2642fdc54b18b0cc3704f920f854875a20bc` |
| `EC-K3-H3-A11-R17-Q6O1307-PROMOTED-LATTICE-ROUTE` | `78be49e31776b18ed36a95d7539bf7d7904cf4c62175be8718ed497613262909` | `be693f7642be184b1a04305bb6eeec336431cdb211791bf1b5212b77f073d7f3` |
| `EC-K3-H3-A11-R17-PHYSICAL-Q10-PROMOTED-ROUTE` | `e8edd71839f682a280196b4c2071f88ab95b41eabffe7b8c349ec16aeef47de8` | `79b70f95f8c610db01c3e0314658d47d64fc951f0be9214ba87da638ad275aad` |
| `EC-K3-H3-A11-R17-Q4O230-Q6O1315-PROMOTED-LATTICE-ROUTE` | `38e70792d8a26433ffa49d7a9793a553438bb58ea06b197625d2a7a9d0d8c909` | `b8e5c8ba7876c00976b0a7814f11bf8b2e07c36c4fe1b1c133db5f89eaa60678` |
| `EC-K3-H3-FIRST-Q8-Q4O11-PROMOTED-LATTICE-ROUTE` | `23521a55f96bd0b10b778138cf28ea833b4f602b7ade710ea01e06f67e503ed0` | `860b3a407cd7e456b2c5ecc1b110c4c1f184fc03e840c6c4359165f0eef2b21b` |
| `EC-K3-H3-D13-R17-Q4O11-PROMOTED-LATTICE-ROUTE` | `d4eae1e0d0b0e9534ab9d90c84ab36118a39c96556bc6e393ff4769002dd8afd` | `a3e241c3ea8c8d88b2c2957556f0d3683ef2f56dd3513a5d210b9171c2f5783a` |

The seven K3 checker files above last entered committed history in `ea2ccd2`;
the ICARM checker last entered it in `ada9371`. The status audit is the
verification that the recorded locks now match those exact source bytes. A
future mathematical replay should still use the commands in the canonical
notes and record any regenerated-output hash separately.

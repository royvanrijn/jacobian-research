# Exact replay commands

The reusable superelliptic leading-block reducer is documented in
[`../SUPERELLIPTIC_DERHAM_ENGINE.md`](../SUPERELLIPTIC_DERHAM_ENGINE.md).
Its fast regression is:

```bash
python3 plane-jc/cas/test_superelliptic_derham.py
```

The source-aware compiler IR and the explicit `(72,108)` tail-basis
certificate are documented in
[`../NEWTON_DERHAM_COMPILER.md`](../NEWTON_DERHAM_COMPILER.md).  Run:

```bash
python3 plane-jc/cas/test_newton_derham_compiler.py
```

This also verifies that the incomplete `(96,144)` and `(75,125)` table rows
are rejected until their Laurent bands are derived; the former additionally
has a source-level lower-side conflict to resolve.

The exact forced F2 `j=1` skeleton and its machine-readable residual
obligations are tested separately:

```bash
python3 plane-jc/cas/test_f2_75_125_frontend.py
python3 plane-jc/cas/f2_75_125_frontend.py
```

The second command emits JSON.  Its `frontend_complete` field is intentionally
false until the lower Laurent boundary has been classified exhaustively.

To compile both audited Proposition 4.3 polygons through lattice supports,
the Laurent chart, all upper bracket layers, and the genus-three first block,
run:

```bash
python3 plane-jc/cas/test_laurent_band_frontend.py
```

The exact 90 MB Zenodo attachment and its extracted source snapshot are pinned
in the repository at
`plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip`
and `plane-jc/external/zenodo-21479814/bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/`.
First check the attachment:

```bash
md5 plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip
shasum -a 256 plane-jc/external/zenodo-21479814/jc2-72-108-exact-certificates-v1.0.1.zip
```

Expected values:

```text
MD5    91255150c8c689b26dc6fb61f9d80aec
SHA256 f7f0876de12d35badbed2be6a773d4a9dada50aff778c080126aab541deefcde
```

After extracting the outer archive and its
`release_bundle/jc2_72_108_exact_replay_v1.0.1.zip`, create a Python
environment from the archived `requirements.txt` and run from
`release_bundle/exact_replay`:

```bash
PYTHON=/absolute/path/to/venv/bin/python ./verify_all.sh
```

Expected final marker:

```text
JC2_72_108_EXACT_REPLAY_PASS
```

For the independent hard-certificate check, run from the repository root:

```bash
/absolute/path/to/venv/bin/python \
  plane-jc/cas/verify_h_certificate_independent.py \
  /absolute/path/to/exact_replay/hne0_polred.pkl \
  /absolute/path/to/exact_replay/hard/h_certificate_exact.txt
```

This checker shares only the serialized four-polynomial input and the text
certificate with the primary computation.  It does not import the generating
code or a CAS.

To verify that the hard membership identity and the specialized \(h=0\)
identity compose to a single unit certificate for the seven pre-division
Case-1 equations, run with the archive's FLINT-enabled Python environment:

```bash
/absolute/path/to/archive/venv/bin/python \
  plane-jc/cas/verify_case1_unit_composition.py \
  /absolute/path/to/exact_replay
```

This checker does not expand the combined multipliers.  It verifies their
factored straight-line representation: all six elimination lifts, the
invertible degree-five descent, and the lift of the \(h=0\) identity to
`1 = sum(A_i G_i) + B h`.

## Critical pinned inputs

| Artifact in `exact_replay/` | SHA-256 |
| --- | --- |
| `firstblock_Q_exact.sing` | `90fb933bf4ae75accddecb69993957db8e289b44c1ae2285dc91b9130f30c062` |
| `firstblock_Q_exact.out` | `2e965d03b39d87531228943cd634de4438d3311fb7a0660dd6dc43a768ae05cb` |
| `case2_compact4_exact.sing` | `38694ed8e9e3b9256b380edf882abefff0ad944113cf0d98e4a951e8f4e31030` |
| `case2_exact_certificate.pkl` | `cfbc3c39d7a28013671144f43ef76f0498542eaf6d562dd624bba3311194e4aa` |
| `case1_residuals_exact.txt` | `f026228c422a213eed853f684b0e2fd98cbc51338ee11c3254a77bf053957c2b` |
| `case1_branch1_after_w_eqs.txt` | `dd6d122161388b8a1961c7e52e92b990e62d0e95c15acc5db56423a1cdf8ea44` |
| `case1_branch2_after_w_eqs.txt` | `0f7c9a8a01a18d725e2a9dc663eec0e4ea018ecb1b1959de1ae7be0214218692` |
| `case1_branch1_after_w.pkl` | `368a1dafdb6d26708b85d652a437c848a7676ba1419e4575ae74186f022621b9` |
| `h0_branch1_exact_certificate.pkl` | `664de005e99bc6a0e61ba479ba64ed57ddbfa9c5399b955cbb12237ac70f8186` |
| `h0_exact_certificate.pkl` | `d7844c2f8edea62be4e3a7fe8a160dc4b7b70efd1262993ec6d2ed7e78722a1a` |
| `hne0_deg35.pkl` | `082471d05a2a7ceebca9fd3a615d8fb6fddaee8ff80afae161a043f71edcb575` |
| `hne0_polred.pkl` | `5a6e423d74ef09fc9c7a7282c500bda566018d7e56a93124665796bbe417cedf` |
| `fixed_matrix_p71.npz` | `3fc0a958672d361a343e8fbeae77c1012f8a5bb4b8e1aa3fb9acdb670c8726dd` |
| `pivot_scalar_rows_p71.npy` | `be4215d7f0303e54890494f7ee03936eff0f6f3e416fe274d7c8cf579938190d` |
| `hard/h_certificate_exact.txt` | `0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a` |
| `verify_all.sh` | `f40bacfa84d915b375b4158109d5120b2d964a8d8012be770750d310abfb5837` |

The exact base rings, variables, orders, localizations, and identity forms are
listed in [the reproduction note](../PAIR_72_108_REPRODUCTION.md).

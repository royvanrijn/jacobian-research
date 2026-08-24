# Stable elliptic-curve commands

These are the small, stable entry points for the programme. Detailed commands
and claim boundaries are in [`../REPRODUCE.md`](../REPRODUCE.md).

## Pinned generation/replay pairs

| Generate | Verify | Result |
| --- | --- | --- |
| `run_crt_lattice_calibration.py` | `verify_crt_lattice_calibration.py` | exact low-rank CRT calibration |
| `run_fermigier_crt_seed.py` | `verify_fermigier_crt_seed.py` | exact local seed, not a target |
| `run_fermigier_rank_certificates.py` | `verify_fermigier_rank_certificates.py` | generic rank lower bound and E22 rank-at-least-22 |
| `run_fermigier_rank20_near_miss.py` | `verify_fermigier_rank20_near_miss.py` | sub-cutoff rank-at-least-20 near miss |
| `run_kihara_rank14.py` | `verify_kihara_rank14.py` | Kihara rank-at-least-14 baseline |
| `run_e29_independence.py` | `verify_e29_independence.py` | public rank-at-least-29 baseline |

`verify_family_data.py` and `verify_benchmarks.py` check the normalized family
metadata. `verify_k3_chain_ledger.py` checks the cross-programme K3 ledger.
`verify_icarm_curve273_rank30_sage.py` is the independent Sage implementation
for the curve-273 certificate.

## Active searches

`search_fermigier_denominator_offsets.py` is the stable denominator-aware
low-conductor search. `evaluate_fermigier_specialization.py` is the common
specialization evaluator. Raw output belongs in
`artifacts/local/elliptic-curves/`.

The much larger set of completed search entry points is preserved under
[`archive/elliptic-curves/cas/`](../../archive/elliptic-curves/cas/) and indexed
by the archive manifest.

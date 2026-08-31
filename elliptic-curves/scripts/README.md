# Stable elliptic-curve commands

These are the small, stable entry points for the programme. Detailed commands
and claim boundaries are in [`../REPRODUCE.md`](../REPRODUCE.md).

## Pinned generation/replay pairs

| Generate | Verify | Result |
| --- | --- | --- |
| `recover_conductor_parameter.py SPEC --output RESULT` | same command with `--check` | two-chart local fingerprint and bounded CRT/Gauss family recognition |
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

Current high-rank and conductor replays are catalogued in
[`../REPRODUCE.md`](../REPRODUCE.md), including curve 302 at rank at least 31,
the full curve-285/286 local conductor replays, the Elkies rank-25--28
positive controls, and the fail-closed residual-descent entry points.

The generic input format and the curve-282 Fermigier replay for
`recover_conductor_parameter.py` are documented in
[`../notes/CONDUCTOR_PARAMETER_RECOVERY.md`](../notes/CONDUCTOR_PARAMETER_RECOVERY.md).

`analyze_elkies_bisection_visibility_and_record_curves.py` converts the split
bisection classes into canonical visible/complementary quotient bases and
performs exact `j`-recognition for the 2024 rank-29 curve and ICARM 273, 302,
and 398--400, with the published rank-28 fibre as a positive control.  Its
`--check` mode replays the pinned artifact.

## Active searches

The compact published Elkies `t` chart is the active high-rank search surface.
Its exact positive controls and residual 2-Selmer gate must be applied before
any expensive point search. The retained conductor-first work starts from the
four exact rank-19/20 near-miss descent inputs. The older
`search_fermigier_denominator_offsets.py` and
`evaluate_fermigier_specialization.py` remain stable calibration tools. Raw
output belongs in `artifacts/local/elliptic-curves/`.

The much larger set of completed search entry points is preserved under
[`archive/elliptic-curves/cas/`](../../archive/elliptic-curves/cas/) and indexed
by the archive manifest.

# Active CAS modules

This directory is intentionally narrower than the computational archive. It
contains exact status checkers and their shared arithmetic, plus code for the
current rank-32, low-conductor, residual-Selmer, and K3-construction gates.
Stable user-facing commands are listed in [`../scripts/`](../scripts/) and
[`../REPRODUCE.md`](../REPRODUCE.md).

## Start here

- `check_icarm_curve302_rank31_pinned.py`: deterministic rank-at-least-31
  replay and compressed-artifact hash check.
- `verify_icarm_curve273_rank30.py`: independent rank-at-least-30 replay.
- `analyze_icarm_7fff_zip_sequence.py`: exact independence replay for the
  public curves 281, 282, 285, and 286; the conductor field is still imported
  public data.
- `verify_icarm_curve245_rank20.py`: fully local low-conductor rank-at-least-20
  certificate.
- `certify_mestre_dsquare_rank19_frontiers.py`: exact rank-at-least-19
  low-conductor frontiers and conditional fixed-fibre diagnostics.
- `certify_nagao_rank20_t5081.py`: exact Nagao rank-at-least-20 certificate.
- `newfamily/certify_rank_t83_6.py`: exact-rank-14 Sage/PARI replay.

## Shared arithmetic

The principal reusable modules are `fermigier_mestre.py`, `mestre_root_tuples.py`,
`nagao_1994.py`, `nagao_linear_sections.py`, `multiple_root_lifting.py`,
`crt_lattice.py`, `finite_quotient_escape.py`,
`mod2_reduction_independence.py`, `mod_l_reduction_independence.py`, and
`pari_bridge.py`. Prefer extending one of these over creating another copy of
the same arithmetic.

## Current research code

- Files containing `bnf_free`, `residual_selmer`, or `curve273` implement the
  unfinished residual 2-Selmer chain. Intermediate success is not a rank
  theorem.
- The retained `fermigier_rank20`, `mixed_small_prime`, and
  `six_root_low_conductor` drivers support the open low-conductor gate.
- `newfamily/` has its own [workflow index](newfamily/README.md).
- The exact exceptional-transport and Mestre two-section checkers remain
  active because they have entries in `MATH_STATUS.json`; their old surrounding
  search campaigns are archived.

## Archive boundary

Superseded versions, completed bounded searches, negative scans, and their
tests/artifacts are indexed under
[`archive/elliptic-curves/`](../../archive/elliptic-curves/). Do not move an
archived bounded result back into the active tree merely because it is
interesting; promote only a compact reproducible result with the correct
evidence label and a canonical note.

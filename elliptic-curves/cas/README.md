# Active CAS modules

This directory is intentionally narrower than the computational archive. It
contains exact status checkers and their shared arithmetic, plus code for the
current rank-32, low-conductor, residual-Selmer, and K3-construction gates.
Stable user-facing commands are listed in [`../scripts/`](../scripts/) and
[`../REPRODUCE.md`](../REPRODUCE.md).

## Start here

- `check_icarm_curve302_rank31_pinned.py`: deterministic rank-at-least-31
  replay and compressed-artifact hash check.
- `verify_icarm_curve356_rank29.py`: independent rank-at-least-29, minimal
  model, conductor, and local-reduction replay for the new rank-29 size record.
- `analyze_icarm_curve356_lineage.py`: hash-pinned curve-351/356 ordered
  denominator and numerical height-Gram comparison; this does not identify a
  family.
- `verify_icarm_curve273_rank30.py`: independent rank-at-least-30 replay.
- `analyze_icarm_7fff_zip_sequence.py`: exact independence replay for the
  public curves 281, 282, 285, and 286, plus repository-local global
  minimality and complete local conductor reconstruction for 285 and 286.
- `verify_icarm_curve245_rank20.py`: fully local low-conductor rank-at-least-20
  certificate.
- `verify_icarm_curve394_rank21.py`: compact-`t=3/8` specialization identity,
  exact generic-plus-public rank-at-least-21 certificate, and complete local
  conductor replay for the current ICARM rank-21 conductor anchor.
- `certify_mestre_dsquare_rank19_frontiers.py`: exact rank-at-least-19
  low-conductor frontiers and conditional fixed-fibre diagnostics.
- `run_conductor_first_pari_diagnostic.py`: bounded PARI `ellrank` supervisor
  for the four conductor-first near misses; provisional upper endpoints stay
  GRH-conditional, while returned points receive independent exact mod-2
  certification.
- `build_conductor_first_s_class_envelopes.py`: exact four-target comparison
  of cubic-field discriminants and materialized Bach/ERH factor-base sizes;
  this orders the BNF-free relation collectors but is not a class-group or
  Selmer computation.
- `run_fermigier_rank20_fixedfb_quadratic_specialq.py`: reusable
  quadratic-in-theta special-q collector for a declared cubic/S-prime set;
  its optional sparse-hypergraph engine cancels multi-large-prime columns
  while retaining exact principal-generator witnesses, and its bounded
  adaptive mode can reuse residual degree-one ideals as new special ideals.
- `run_fermigier_rank20_minkowski_specialq.py`: full-ideal Minkowski
  special-q collector.  In addition to its historical one-large-prime graph,
  the opt-in `--norm-factor-mode exact --large-prime-merge-mode
  sparse-hypergraph` path retains fully factored multi-residual supports and
  exact dependency witnesses.  It supplies declared discriminant factors to
  PARI before maximal-order construction.  This is relation collection, not
  class-group completion or a Selmer calculation.
- `certify_nagao_rank20_t5081.py`: exact Nagao rank-at-least-20 certificate.
- `newfamily/certify_rank_t83_6.py`: exact-rank-14 Sage/PARI replay.
- `elkies_residual_selmer_gate.py`: fail-closed rank-32 residual-dimension
  policy; signatures and incomplete descent cannot authorize point search.
- `run_elkies_2026_rank28_residual_selmer.py`: resource-bounded genuine PARI
  or Selmer-only eclib 2-descent on the public rank-28 positive control. Its
  `pari-factored` backend consumes the separately proved complete
  2-division-discriminant factorization and records a strict memory/wall cap.
- `build_elkies_2026_rank28_bad_place_ledger.py`: proves the complete
  2-division-discriminant factorization and computes the generic-seventeen
  Kummer images at every bad finite place, at 2, and at infinity. The complete
  known-image ledger is exact input to descent, not a Selmer upper bound.
- `run_elkies_2026_rank28_s_class_pari.py`: isolates factor-supplied maximal-
  order certification, class-group relation generation, and `bnfcertify` in
  an owned resource-bounded worker. Its `--field-model polredabs` path records
  an exact reduced polynomial and generator map before the same certified
  maximal-order calculation. A completed class quotient would still precede
  the separate norm and local-solubility gates.
- `run_fermigier_rank20_minkowski_specialq.py --elkies-rank28`: feeds the same
  proved factor support into the exact BNF-free principal-relation collector.
  Its factor-base-1000 paired-cover pilot is explicitly uncertified and does
  not authorize search.
- `build_elkies_2026_rank28_local_coverage.py`: recomputes the bad-place
  signatures of the generic 17 and public complement 11. The eleven globally
  independent exceptional directions add zero local-signature rank, an exact
  positive control showing that signature rank is not a Mordell--Weil or
  Selmer quotient dimension.
- `audit_bnf_free_local_kummer_coverage.py`: compares known-point projections
  with exact odd-prime and real local-dimension bounds. Equality certifies a
  full local image at that place; every strict inequality and the two-adic
  place remain unresolved.
- `run_bnf_free_two_cover_local_supervisor.py`: replaces monolithic bounded
  cover audits with owned, resumable `(cover,p)` workers. Each cache block is
  input- and limit-bound; smooth lifts and obstructions are retained, while
  timeouts and state caps stay mathematically inconclusive.
- `build_elkies_2026_rank28_public_selmer_controls.py`: converts the certified
  public complement eleven into exact classes `X(Q)-theta`, square norms, and
  rational two-cover witnesses. Together with the finite-reduction quotient
  certificate this proves a residual 2-Selmer lower bound of 11 and validates
  the cover layer on genuine classes; it supplies no upper bound.
- `build_elkies_2026_rank28_relative_descent_magma.py`: replays the certified
  generic Kummer image and emits an unconditional basis-level Selmer job whose
  rejection gate precedes all residual-cover construction.
- `parse_elkies_2026_rank28_relative_descent.py`: accepts only a complete log
  from the exact generated Magma source and emits the common fail-closed gate.

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
- For Elkies compact-`t` candidates, residual dimension below 15 rejects rank
  32. Only a completed unconditional global/local descent of dimension at
  least 15 may unlock a same-curve two-cover or expensive point search.
- The retained `fermigier_rank20`, `mixed_small_prime`, and
  `six_root_low_conductor` drivers support improved-conductor and residual
  descent work now that the original rank-21/cutoff branch is closed.
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

# Active elliptic-curve artifacts

This directory is the complete active artifact surface for
[`elliptic-curves/`](../../../elliptic-curves/README.md). Exploratory outputs and superseded results live in the
[`archive`](../../../archive/elliptic-curves/README.md); large or resumable
runs belong under the ignored `artifacts/local/elliptic-curves/` directory.

[`CATALOG.tsv`](CATALOG.tsv) is the index. Every row gives the exact SHA-256,
an evidence label, the governing `MATH_STATUS.json` identifier when one
exists, and a one-line scope statement. The evidence labels are deliberately
strict:

- `theorem-certificate`, `exact-rank-certificate`, and
  `exact-lower-bound-certificate` are proof-bearing;
- `conditional-bound` states its hypothesis;
- `exact-computation` proves only the calculation described;
- `bounded-experiment` never promotes a negative search to a theorem;
- `partial-reproduction` records precisely what remains dependent on a public
  source or missing local replay;
- `source-transcription`, `reproducibility-fixture`, and `search-plan` carry no
  mathematical conclusion by themselves.

Important distinctions made explicit by the catalogue:

- `latent_lattice_calibration_truth_v1.json` contains exact withheld control
  embeddings; it is not selector input. `latent_lattice_calibration_v2.json`
  is the active corrected-semantics gate artifact.  It still rejects the
  bounded selector before any target search.  The superseded `v1` selector
  bytes are retained only for provenance and are not a nonexistence theorem
  for a common wgxli lattice.
  `latent_lattice_finite_calibration_v1.json` adds disjoint development and
  held-out finite-code ensembles.  It passes R17 proposal recall but fails the
  Fermigier control and finite-profile selection, so it also forbids target
  use.

- `icarm_273_282_302_family_discovery_v1.json` screens 2,334 generated
  one-parameter families exactly.  It rediscovers curve 282 in both the
  canonical Fermigier coordinate `u=11671/42` and the generated six-root
  Mestre coordinate `T=11671/21`, verifies `Q`-isomorphism in both models, and
  finds no match for curves 273 or 302 in this bounded construction space.
- `icarm_curve282_conductor_parameter_recovery_v1.json` verifies global
  minimality, selected local Tate data, two-chart discriminant-root profiles,
  bounded CRT/Gauss recovery of `u=11671/42`, and the exact Fermigier/target
  `j`-identity.  Its exact leaf residues are replay inputs; the valuations
  alone determine only coarser p-adic balls.
- `record_first17_subgroups_v1.json` exactly computes the first-seventeen
  coordinate, quotient, finite-Kummer, and bad-component codes for curves 273
  and 302, while its canonical-height/theta layer is numerical at 100 digits.
  Its saturation index one is only inside each displayed subgroup, not in the
  full Mordell--Weil group.
- `icarm_curve302_rank31_v1.json.gz` proves rank at least 31, not exact rank 31.
- The ICARM 285/286 analysis exactly proves independence of 21 displayed
  points and now independently replays global minimality and every local
  conductor exponent.
- `icarm_curve394_rank21_v1.json` specializes the compact Elkies R17 family at
  `t=3/8`, proves a generic-17 plus public-4 basis independent, and replays the
  exact conductor locally.  It proves rank at least 21 at
  `log(N)=166.252098...`, not exact rank 21.
- `conductor_first_near_miss_descent_targets_v1.json` pins exact known-subgroup
  Kummer inputs for four fixed fibres; it does not claim a complete Selmer
  group or rank upper bound.
- `conductor_first_family_anchor_pilot_v1.json` is the closed 27-fibre
  Fermigier/Mestre anchor-neighborhood ledger. Nine sieve survivors have
  exact global/local Tate data and full-dimensional known mod-2 Kummer
  subgroups. The `u=481` fibre has an exact rank-at-least-14 point certificate;
  no complete Selmer or residual-cover classification is claimed.
- The Elkies compact-`t` artifacts certify the rank-25--28 positive controls
  and the complete height-10000 three-block calibration. The PARI and eclib
  rank-28 residual-descent artifacts are strict timeouts with no Selmer bound;
  both explicitly forbid point search.
- `elkies_2026_bisection_specialization_controls_v1.json` is the complete
  195,600-test evaluation of the 39,120 equation-level bisections at the four
  rank-25--28 controls and ICARM curve 394. It finds split counts
  `6,3,2,1,25` and known-complement class-span dimensions `5,3,2,1,4`, with
  no finite-quotient escape. Full-rank relation blocks verified by exact group
  addition prove generated-subgroup ranks `25,26,27,28,21`; these do not give
  upper bounds for the full curves.
- `elkies_2026_bisection_visibility_record_curves_v1.json` row-reduces those
  classes into deterministic visible and complementary quotient bases.  It
  records the ten-dimensional rank-28 target packet and the exact mechanism
  boundary that translated trace shells repeat existing bisection classes.
  It also proves, through degree-24 irreducibility witnesses, that the 2024
  rank-29 curve and ICARM 273, 302, and 398--400 are not rational fibres of the
  published `R17` fibration, including after quadratic twisting.  The rank-28
  control recovers `5471*t+9529`; other fibrations, families, and isogeny
  constructions remain open.
- `elkies_2026_rank28_bad_place_kummer_ledger_v1.json` proves the complete
  factorization of the rank-28 2-division cubic discriminant and contains all
  thirteen finite/2-adic/real local blocks for the generic seventeen points.
  Their combined coordinate rank is 15; this is a known-image calculation,
  not an ambient `K(S,2)`, local-solubility, or Selmer certificate.
- `elkies_2026_rank28_residual_2selmer_pari_factored_8g_v1.json` is the first
  PARI descent supplied with that factor certificate. Its strict 600-second
  run reached 5,698,514,944 bytes peak observed RSS but returned no Selmer
  dimension, so it remains search-forbidden. The longer supervised artifact
  has suffix `_8g_30min_v1`; its strict 1,800-second run reached
  6,040,723,456 bytes peak observed RSS and likewise returned no dimension.
- `elkies_2026_rank28_s_class_pari_v1.json` isolates the next PARI stage. The
  exact factor-supplied maximal-order setup completes, but the strict
  120-second run stops inside class-group relation generation before
  `bnfcertify`; it is not an `S`-class or Selmer bound.
- `elkies_2026_rank28_s_class_pari_polredabs_v1.json` repeats that envelope on
  PARI's exact reduced cubic with an explicit original-generator map. The
  polynomial order index drops by 27, but the run reaches the same 153-request
  random-relation plateau and returns no class-group bound. It remains
  search-forbidden.
- `elkies_2026_rank28_bnf_free_s_class_pilot_v1.json` audits 172 exact
  canonical principal rows after a paired-special-ideal pilot. Its displayed
  factor-base quotient dimension 141 is explicitly uncertified because bound
  1,000 is below the 1,202,640 Bach/ERH generation threshold; the pilot found
  no noncanonical relation and forbids point search.
- `elkies_2026_rank28_generic17_local_signature_v1.json` recomputes all 53
  bad-place coordinates for the generic 17 and public complement 11. The
  generic and full rank-28 local-signature ranks are both 15, with zero
  incremental rank from every exceptional direction; this is exact evidence
  that known Kummer signatures are not a Mordell--Weil quotient or Selmer
  bound.
- `elkies_2026_rank28_generic17_local_coverage_v1.json` certifies full known-
  point coverage at four of eleven odd bad primes and at the real place. The
  other seven odd places and the two-adic place remain unresolved.
- `elkies_2026_rank28_norm_one_local_pilot12_v1.json` is a bounded selected-
  place audit of 12 of 49 norm-one cover candidates. Of 84 cover/place tasks,
  60 have certified local points and 24 remain inconclusive; none has a
  certified local obstruction. This is not class-group completion, everywhere
  local solubility, Selmer membership, or search authorization.
- `elkies_2026_rank28_public11_selmer_candidates_v1.json` records the exact
  cubic Kummer classes `X(Q)-theta` of the eleven certified public complement
  points. Their norms are displayed squares and their quotient independence
  proves a residual 2-Selmer lower bound of 11.
- `elkies_2026_rank28_public11_two_cover_controls_v1.json` materializes those
  eleven intersections of quadrics and verifies `[1:0:0:1]` on each one.
  `elkies_2026_rank28_public11_global_cover_witness_audit_v1.json` rechecks the
  witnesses and records local solubility at every place. These are positive
  controls, not an ambient Selmer enumeration or upper bound.
- `newfamily_rank14_t83_6_v1.json` proves only rank at least 14; the separate
  `newfamily_rank14_t83_6_pari_exact_rank_v1.json` supplies the PARI interval
  `[14,14]` used for the exact-rank statement.
- Bounded Fermigier searches are indexed as experiments even when every
  calculation inside the declared box is exact.

Run the catalogue/archive integrity check with:

```sh
python3 elliptic-curves/scripts/audit_artifact_catalog.py
```

The pre-cleanup bytes and every provenance-only refresh are recorded under
[`archive/elliptic-curves/`](../../../archive/elliptic-curves/README.md).

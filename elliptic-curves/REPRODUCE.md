# Reproducing the elliptic-curve programme

Run commands from the repository root. This catalogue covers the active
certificates and current research gates. The pre-cleanup catalogue, including
commands for every bounded historical scan, is preserved as
[`REPRODUCE_2026-08-24.txt`](../archive/elliptic-curves/REPRODUCE_2026-08-24.txt).

## Environment

The dependency-free checks use the repository virtual environment:

```sh
.venv/bin/python --version
```

Some certificates additionally require PARI/GP, Sage/eclib, Singular, or
Magma. Those requirements are stated beside their commands. Raw output and
restart state belong under the ignored `artifacts/local/elliptic-curves/`
tree; do not overwrite pinned generated results during an exploratory run.

## Standard checks

Compile the active Python surface, validate links/status/layout, and run the
current elliptic-curve regression suite:

```sh
make check
make verify-elliptic-curves
```

Audit the evidence labels, JSON/gzip readability, generators, and coverage of
the compact current artifact directory:

```sh
.venv/bin/python elliptic-curves/scripts/audit_artifact_catalog.py
```

## Primary record certificates

### ICARM curve 356: rank at least 29 and the new rank-29 size record

Replay point membership, trivial torsion, the exact finite-quotient
independence certificate, global minimality, every local reduction, conductor,
root number, and the complete reported prime factorizations:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve356_rank29.py \
  --verify-primality
```

Replay the hash-pinned public-source comparison with curve 351 and the
80-digit PARI height-Gram fingerprint of the first seventeen displayed
points:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_icarm_curve356_lineage.py
```

The first command proves `rank E(Q) >= 29`. The second is numerical evidence
for a common ordered 17-section template, not a family-recognition theorem.
See
[`ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md`](notes/ICARM_CURVE356_RANK29_AND_CONSTRUCTION.md).

Replay the later complete same-submitter sweep, its three bounded residual
components, and the five-fibre 351/356/376/377/385 interpolation export:

<!-- status-consumer: EC-ICARM-WGXLI-R17-LINEAGE 90790392f558f0a0 -->

```sh
.venv/bin/python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_lineages.py

.venv/bin/python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_lineages.py \
  --write-artifact \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_lineage_v1.json
```

This pins three additional candidate fibres and canonical short-model data;
it does not prove a common family or rootless-K3 realization.  See
[`ICARM_WGXLI_RANK17_LINEAGE.md`](notes/ICARM_WGXLI_RANK17_LINEAGE.md).

Calibrate the basis-independent latent-lattice selector against the four
published R17 positive controls and the exact Fermigier rank-12 negative
control before using it on the wgxli cluster:

```sh
PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/build_latent_lattice_calibration_truth.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_method.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_finite_aware_latent_lattice.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_shape.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_relation_consensus.py --check

sage -python \
  elliptic-curves/cas/calibrate_latent_lattice_hypergraph_matcher.sage --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_metric_relation_search.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_partial_replay.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 \
  elliptic-curves/cas/calibrate_latent_lattice_star_component.py --check
```

The active `v2` artifact uses corrected unit/scaled relation semantics and
record indexing.  Its bounded selector still fails the control gate, despite
later high-recall diagnostics recovering all four R17 spaces in their
proposal ledgers and improving ICARM 245 to 11/12 before extension.  Its status is therefore
`FAIL_CALIBRATION_TARGET_GATE_CLOSED`: it makes no claim about the existence
of a common lattice on curves 351, 356, 376, 377, and 385.  See
[`LATENT_LATTICE_CALIBRATION.md`](notes/LATENT_LATTICE_CALIBRATION.md).

The finite-aware replay uses disjoint development and held-out good-reduction
quotient ensembles.  It recovers R17 in all four bounded proposal ledgers but
reaches only 11/12, 11/12, and 8/12 on the three Fermigier controls; finite
profile matching leaves the rank-25 truth at rank 188.  Its status is
`FAIL_FINITE_PROPOSAL_RECALL`, so it likewise leaves the target gate closed.

The exact hypergraph validator recovers primitive rectangular `17 x r` maps
when supervised ray injections are supplied, and full-cloud replay sees
238--304 training-core rays in every R17 held-out fibre.  The blind bounded
metric/relation search nevertheless reaches only 49 replayed rays against a
100-ray gate on rank 25.  Its status is
`FAIL_BLIND_R17_RECOVERY_GATE_CLOSED`; no wgxli curve is loaded.

The proper-subspace replay calibrator saturates a supervised rank-16 R17
relation path before lifting it.  Its exact lift replays 194 source rays and
318 ternary relations and has primitive target image; the attached finite
signature uses six disjoint good-reduction quotient blocks.  A separate
400-audit oracle-center beam still loses this component and accepts no full
embedding.  Its status is `PASS_EXACT_PARTIAL_REPLAY_SELECTOR_FAIL`, so the
new invariant is a validator but not yet an authorized blind selector.  This
command requires system `python3` with NumPy and PARI/GP; the repository
`.venv` does not currently provide NumPy.

The center-star component experiment assigns several incident ternary
relations jointly and enforces exact mod-2/mod-3 rank compatibility before
partial replay.  In its oracle-center rank-25 box, 4,152 finite-incompatible
branches are rejected, but the best of 476 exact candidates intersects the
rank-11 visible truth star in dimension only 9.  The more expensive pair/triple
finite-matroid diagnostic gives the same result.  Its status is
`FAIL_STAR_COMPONENT_RECALL_GATE_CLOSED`; it is a bounded proposal failure,
not evidence against R17.

The shape replay adds exact cross-bound rational intersections and primitive
closure.  It tests the top 200 height-28 enclosures against the top 200
height-29 enclosures on ICARM 245 and puts the true rank-12 Fermigier subgroup
at rank 65 among 2,939 exact candidates.  All four R17 truth spaces occur in
the finite-seeded ledgers, and the rank-25 held-out experiment succeeds, but
symmetric leave-one-out selection succeeds on only one of four R17 fibres.
Its status is `PASS_PROPOSAL_CALIBRATION_SELECTOR_FAIL`; target use remains
forbidden.  The full `--check` regenerates the control clouds and can take
several minutes.

The relation-consensus replay is a supervised signal audit.  In each R17
leave-one-out split, the exact coefficient rays supported by two training
fibres and visible on the held-out fibre span rank 17.  This proves that a
full-rank exact common-relation signal survives the numerical cutoffs; it does
not recover the unknown integral alignments and is not a blind selector.

Eliminate the section and surface coefficients from the literal five-fibre
rootless-K3 ansatz, then exhaust every normalized projective first-jet chart
over `GF(17)`, `GF(53)`, and `GF(67)`:

<!-- status-consumer: EC-ICARM-WGXLI-R17-FIRST-JET 11b13e24c5e42a14 -->

```sh
sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 17 --jobs 4 --threads 1 --pair-timeout 60 --check \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod17_v2.json

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 53 --jobs 14 --threads 1 --pair-timeout 60 --check \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod53_v2.json

sage -python \
  elliptic-curves/cas/eliminate_icarm_wgxli_rank17_first_jet.sage \
  --prime 67 --jobs 32 --threads 1 --pair-timeout 180 \
  --reuse-unit-outputs --check \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_first_jet_mod67_v1.json
```

All 210, 2550, and 4160 fixed-parameter ideals respectively have Groebner
basis `[1]`, with no timeout. The compact certified published-R17 model and
all seventeen sections give nondegenerate positive controls at 53 and 67 in
the finite and both infinity orientations. This exact necessary-condition
obstruction does not exclude bad or colliding parameter reduction at all
three primes, a changed Mordell--Weil basis, or a different family shape. See
[`ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md`](notes/ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md).

Replay the bounded relative-sign and fingerprint-permutation search, the
one-elementary-mutation search, and the artifact-bound formal rejection:

<!-- status-consumer: EC-ICARM-WGXLI-R17-BOUNDED-REBASING 6e0c7b116b5b25c3 -->

```sh
sage -python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_rebasing.sage --check

.venv/bin/python \
  elliptic-curves/cas/construct_icarm_wgxli_rank17_bounded_mutation.py --check

.venv/bin/python \
  elliptic-curves/cas/certify_icarm_wgxli_rank17_bounded_rejection.py --check
```

The declared search has one retained exact mutation, `P4 -> P4-P1`; its
complete mod-17 and mod-53 projective first-jet charts are both empty. This is
a finite rejection inside the recorded sign, permutation, one-shear,
coefficient, anchor, and height bounds, not an unrestricted `GL(17,Z)`
search. See
[`ICARM_WGXLI_RANK17_BOUNDED_REBASING.md`](notes/ICARM_WGXLI_RANK17_BOUNDED_REBASING.md).

### ICARM curve 302: rank at least 31

The fast checker verifies both pinned hashes and recomputes the complete exact
certificate:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py
```

To generate an unpinned plain-JSON replay with optional primality checks:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve302_rank31.py \
  --output artifacts/local/elliptic-curves/icarm_curve302_rank31_v1.json \
  --verify-primality
```

This proves `rank E(Q) >= 31`, not an unconditional exact rank. See
[`ICARM_CURVE302_RANK31.md`](notes/ICARM_CURVE302_RANK31.md).

### ICARM curve 273: rank at least 30

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve273_rank30.py --check
```

The independent Sage replay is:

```sh
sage -python elliptic-curves/scripts/verify_icarm_curve273_rank30_sage.py
```

See [`ICARM_CURVE273_RANK30.md`](notes/ICARM_CURVE273_RANK30.md).

### Comparative height lattices: ranks 28--31

PARI/GP is required.  Compute the 100-digit canonical height matrices, LLL
transforms and reduced Grams for the public rank-28, rank-29, curve-273 and
curve-302 point lists:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/compare_record_height_lattices.py --digits 100
```

The bounded additive short-vector search is run separately for each declared
height cutoff.  For example, the rank-29 control is:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/search_record_rank17_core.py \
  rank29 --bound 60 --additive-pair-limit 1507 \
  --pool 300 --trials 2000 --seed 29017 \
  --output artifacts/local/elliptic-curves/rank29-r17-additive-ransac.txt
```

The matching commands for `rank28`, `curve273`, and `curve302` use respectively
`(--bound, --additive-pair-limit, --trials, --seed)=(60,2423,800,28017)`,
`(65,2500,800,27317)`, and `(70,3000,800,30217)`.  With those four local
search files present, exactly saturate the candidate spaces and generate their
100-digit core Grams and approximate 1,311-vector profiles:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/analyze_record_rank17_candidates.py --digits 100
```

Calibrate the forced rank-17 fingerprint against ICARM curve 245, whose exact
Fermigier--Mestre rank-12 parent is already reconstructed.  This also fits
exact unimodular R17 shell bases and evaluates the out-of-sample integral-point
enrichment of each selected core.  It exactly transports the thirteen known
generic curve-245 points into the public basis, verifies their rank-12 span and
relation, measures its intersection with the forced rank-17 control, and
replays all 1,311 R17 minimal lines through every fitted basis:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_record_rank17_fingerprint.py \
  --digits 100 --restarts 64
```

The negative control fails for both the optimized R17 basis-entry score and
direct pairwise GL(17,Z)-plus-scale fitting.  The fitted full shells are also
highly dispersed.  Integrality enrichment survives as evidence for structured
cores, but the exact curve-245 replay shows that it need not recover the true
generic subgroup.  None of these calculations is a specialization
certificate.

Run the stronger full-shell coordinate descent on all four records and the
known negative control:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/search_record_rank17_shell_embedding.py \
  rank28 rank29 curve273 curve302 curve245-negative-control \
  --restarts 16 --random-steps 80 --maximum-sweeps 25
```

This minimizes the variation of all 1,311 mapped R17 minimal-vector heights.
It also fails calibration: the curve-245 negative control scores better than
the known R17-positive rank-29 specialization, so the values for curves 273
and 302 are diagnostic only.

Extend the exact bounded Mestre fingerprint census to curve 302:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_icarm_construction_fingerprints.py \
  --include-curve302
```

These are numerical/bounded provenance calculations, not a K3 specialization
certificate.  See
[`RECORD_CURVES_28_29_273_302_HEIGHT_LATTICES.md`](notes/RECORD_CURVES_28_29_273_302_HEIGHT_LATTICES.md).

### ICARM 273/282/302: generated family discovery

Screen the declared polynomial and generated six-root Mestre construction
space, factor every modular survivor over `QQ`, and verify rational
isomorphism for every exact parameter:

```sh
.venv/bin/python elliptic-curves/scripts/discover_record_families.py \
  elliptic-curves/data/family-discovery/icarm_273_282_302.json \
  --output artifacts/generated-results/elliptic-curves/icarm_273_282_302_family_discovery_v1.json \
  --check
```

The bounded search tests 2,334 distinct families.  It rediscovers curve 282
at Fermigier parameter `u=11671/42` and generated six-root parameter
`T=11671/21`, with exact rational-isomorphism scales `882` and `147`.  It
finds no match for curves 273 or 302 in this declared space; that is not a
nonexistence theorem for other constructions.  See
[`GENERATED_FAMILY_DISCOVERY.md`](notes/GENERATED_FAMILY_DISCOVERY.md).

### ICARM curve 282: local-conductor parameter recovery

Replay the generic two-chart discriminant-root, CRT, Gauss-reduction, and
exact `j`-recognition pipeline on the known Fermigier fibre:

```sh
python3 elliptic-curves/scripts/recover_conductor_parameter.py \
  elliptic-curves/data/conductor-engineering/icarm_curve282_fermigier.json \
  --output artifacts/generated-results/elliptic-curves/icarm_curve282_conductor_parameter_recovery_v1.json \
  --check
```

This recovers `u=11671/42` from the five declared local branches and verifies
the target's global minimality and exact family `j`-match.  The branch
residues are replay inputs: the discriminant valuations alone select coarser
p-adic balls.  See
[`CONDUCTOR_PARAMETER_RECOVERY.md`](notes/CONDUCTOR_PARAMETER_RECOVERY.md).

### ICARM curves 285 and 286: low-conductor candidates

```sh
.venv/bin/python elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py --check
```

The command exactly replays equation membership, finite-reduction
independence, torsion, invariants, and pairwise `j` comparisons for curves
281, 282, 285, and 286. For curves 285 and 286 it proves the 21-point rank
lower bound, verifies that the displayed model is global minimal, and
reconstructs the exact conductor from every local Tate datum. See
[`ICARM_7FFF_ZIP_SEQUENCE.md`](notes/ICARM_7FFF_ZIP_SEQUENCE.md).

### Elkies 2026 compact-t rank-32 search gates

Replay the compact rank-17 model, the exact q12 coordinate match, and the four
rank-25--28 public positive controls:

```sh
SAGE=/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python
$SAGE elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
$SAGE elkies-k3/scripts/match_h92_q12o5867_to_elkies_2026_qq.sage
python3 elliptic-curves/scripts/verify_elkies_2026_high_rank_calibrations.py
```

The control checker proves combined quotient gains `8,9,10,11` over the fixed
generic rank 17. Run the complete height-10000 three-ensemble calibration with:

```sh
python3 elkies-k3/scripts/calibrate_elkies_2026_positive_controls_nagao.py
```

Its weakest-block ranks strongly recover all four controls, but remain
heuristic. Attempt the actual rank-28 residual 2-Selmer computation with:

```sh
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --timeout 300 --overwrite
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend eclib --timeout 300 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_eclib_v1.json \
  --overwrite
```

Both pinned runs are exact-backend timeouts, not Selmer bounds. Each emits a
fail-closed artifact: residual dimension below 15 would exactly reject rank
32, while only a completed unconditional dimension at least 15 authorizes
two-cover solving or expensive point search. The authorization is bound to the
same global minimal model.

The rank-28 2-division cubic discriminant has a complete, proved
factorization. Rebuild its factor proof and the Kummer images of the generic
seventeen points at all bad finite places, at 2, and at infinity with:

```sh
python3 elliptic-curves/scripts/specialize_q12o5867_candidate.py \
  --a -9529 --b 5471 --overwrite
python3 elliptic-curves/cas/build_elkies_2026_rank28_bad_place_ledger.py \
  --overwrite
```

All thirteen local blocks complete, and their combined known-image coordinate
matrix has rank 15. This exact ledger is not an ambient `K(S,2)`, class-group,
local-solubility, or Selmer upper-bound certificate. Use its factor table to
avoid repeating the hard discriminant factorization inside PARI:

```sh
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --backend pari-factored --timeout 600 \
  --pari-stack-bytes 8000000000 --rss-limit-bytes 12000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_pari_factored_8g_v1.json \
  --overwrite
```

The pinned 600-second run enters the descent immediately and reaches
5,698,514,944 bytes peak observed RSS, but still returns no Selmer dimension.
The pinned 1,800-second run with the same 8 GB stack reaches 6,040,723,456
bytes and likewise returns no dimension. Both remain fail-closed. Reproduce
the longer resource envelope by changing `--timeout` to `1800` and the output
suffix to `_8g_30min_v1`; never overwrite a completed Selmer result with a
timeout.

Isolate the preceding cubic `S`-class computation with a stage-aware,
resource-bounded worker:

```sh
python3 elliptic-curves/cas/run_elkies_2026_rank28_s_class_pari.py \
  --timeout 120 --c1 0.01 --c2 4 --nrpid 20 \
  --pari-stack-bytes 2000000000 --rss-limit-bytes 3000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_s_class_pari_v1.json \
  --overwrite
```

The pinned run completes factor-supplied `nfinit` and `nfcertify`, then
reaches the strict limit inside `bnfinit` at 265,261,056 bytes peak observed
RSS. It never reaches `bnfcertify`, so it supplies no class-group or Selmer
bound. Its debug tail identifies the random-relation plateau without treating
PARI progress messages as a certificate.

Replay the exact absolute polynomial reduction and matched resource envelope:

```sh
python3 elliptic-curves/cas/run_elkies_2026_rank28_s_class_pari.py \
  --field-model polredabs --timeout 120 \
  --c1 0.01 --c2 4 --nrpid 20 \
  --pari-stack-bytes 2000000000 --rss-limit-bytes 3000000000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_s_class_pari_polredabs_v1.json \
  --overwrite
```

`polredabs` gives the depressed cubic

```text
x^3
- 35676022072134269484503481261046298223875964999429256003*x
- 81734190921553911625559669772737848345984148653181341176726216553622238508296306498
```

with exact original-generator map `theta=-3*x+1`. Its polynomial
discriminant is the original one divided by `3^6`, its defining-order index is
smaller by 27, and factor-supplied `nfinit` plus `nfcertify` proves the same
maximal cubic field. The matched run nevertheless reaches the same 243-ideal,
153-request random-relation plateau and times out inside `bnfinit` at
264,839,168 bytes peak observed RSS. Polynomial reduction is therefore kept
as an exact input optimization but not promoted to a class-group or Selmer
bound.

Replay the complementary BNF-free paired-special-ideal pilot and its exact
principal-row audit with:

```sh
mkdir -p artifacts/local/elliptic-curves/elkies-rank28-bnf-free
$SAGE elliptic-curves/cas/run_fermigier_rank20_minkowski_specialq.py \
  --elkies-rank28 --factor-base-bound 1000 \
  --special-q-min 1009 --special-q-max 5000 --max-special-q 20 \
  --special-ideal-mode cycle-pairs --pair-cycle-length 20 \
  --trial-prime-bound 1000 --lattice-combination-bound 2 \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_v1.json
$SAGE elliptic-curves/cas/augment_bnf_free_canonical_principal_relations.py \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_v1.json \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_canonical_v1.json
$SAGE elliptic-curves/cas/audit_bnf_free_s_class_quotient.py \
  --relation-ledger artifacts/local/elliptic-curves/elkies-rank28-bnf-free/minkowski_fb1000_paircycle20_canonical_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_bnf_free_s_class_pilot_v1.json
```

The collector retains 10,288 sampled generators; augmentation adds 172 exact
canonical generators and principal rows. No noncanonical relation closes.
The resulting 327-column, 26-`S`-column quotient model has relation rank 172
and displayed dimension 141, but factor-base bound 1,000 is below the
1,202,640 Bach/ERH generation bound. The audit therefore reports
`UNCERTIFIED_FACTOR_BASE`, not an `S`-class or Selmer upper bound, and forbids
point search.

Rebuild the exact rank-28 local-signature positive control and audit where
the generic seventeen already span the full local Kummer image:

```sh
$SAGE elliptic-curves/cas/build_elkies_2026_rank28_local_coverage.py \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --overwrite
$SAGE elliptic-curves/cas/audit_bnf_free_local_kummer_coverage.py \
  --signature-map artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_coverage_v1.json
```

The first command independently recomputes the local rows of the generic 17
and public complement 11 and checks the generic rows against the bad-place
ledger. The full 28-point local-signature rank is still 15: all eleven
globally independent exceptional directions add zero rank in every bad-place
block. This is a positive control proving that the signature cannot stand in
for the Mordell--Weil quotient. The coverage audit certifies the full local
Kummer image only at odd primes `3`, `19`, `20650099`, and
`315574902691581877528345013999136728634663121`, plus the real place. Seven
odd places and the two-adic place remain unresolved.

Generate a bounded exact norm-one tranche, materialize its intersections of
quadrics, and audit selected finite places with resumable per-cover/per-place
workers:

```sh
$SAGE elliptic-curves/cas/generate_q12o5867_norm_one_cover_candidates.py \
  --signature-map artifacts/generated-results/elliptic-curves/elkies_2026_rank28_generic17_local_signature_v1.json \
  --coefficient-bound 2 \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_v1.json
$SAGE elliptic-curves/cas/build_bnf_free_two_covers.py \
  --candidates artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_v1.json \
  --output artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_covers_v1.json
python3 elliptic-curves/cas/run_bnf_free_two_cover_local_supervisor.py \
  --covers artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_covers_v1.json \
  --primes 3,5,7,11,13,17,19 --max-covers 12 \
  --max-enumeration-prime 19 --max-lift-precision 6 --max-lift-states 5000 \
  --timeout-per-place 2 \
  --cache-dir artifacts/local/elliptic-curves/elkies-rank28-bnf-free/norm_one_b2_local_pilot_blocks_v1 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_norm_one_local_pilot12_v1.json
```

The pinned pilot has 84 cover/place tasks: 60 smooth-mod-`p` local witnesses,
19 state-capped singular lift trees, five supervised timeouts, and no proved
local obstruction. It is deliberately truncated to 12 of 49 norm-one
candidates and seven finite primes. Use `--retry-incomplete` to revisit only
timed-out cache blocks. No survivor is a locally soluble global class, no
Selmer dimension is obtained, and no point search is authorized.

Calibrate the same cover layer on the eleven *genuine* residual classes from
the public rank-28 complement:

```sh
python3 elliptic-curves/cas/build_elkies_2026_rank28_public_selmer_controls.py \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_selmer_candidates_v1.json \
  --overwrite
$SAGE elliptic-curves/cas/build_bnf_free_two_covers.py \
  --candidates artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_selmer_candidates_v1.json \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_two_cover_controls_v1.json
$SAGE elliptic-curves/cas/audit_bnf_free_two_cover_reduction.py \
  --covers artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_two_cover_controls_v1.json \
  --primes 2,3,5,7,11,13,17,19,48463,20650099,315574902691581877528345013999136728634663121,376018840263193489397987439236873583997122096511452343225772113000611087671413 \
  --max-enumeration-prime 2 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_public11_global_cover_witness_audit_v1.json
```

For every public complement point `Q`, the first command constructs
`alpha=X(Q)-theta`, checks `Norm(alpha)=Z(Q)^2`, and records the cover point
`[1:0:0:1]`. The cover builder and local audit recheck that witness exactly;
it proves a rational point, hence local solubility at every place, on all
eleven covers. The existing finite-reduction certificate proves these classes
independent modulo the generic seventeen. Thus the residual 2-Selmer quotient
has certified dimension at least 11 on this positive-control fibre. This is a
lower bound only: it still lacks the complete ambient class computation and
does not pass the rank-32 threshold 15.

For an unconditional basis-level external computation, generate the exact
rank-28 Magma job with:

```sh
python3 elliptic-curves/cas/build_elkies_2026_rank28_relative_descent_magma.py \
  --output artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.m
```

It computes `TwoSelmerGroup(E : Bound := -1)` before any cover construction,
records both the dimension modulo the generic 17 and the equivalent dimension
beyond the full known rank-28 subgroup, and exits below residual dimension 15
(equivalently below four unexplained directions). Only its passing branch
calls `TwoDescent`, with all 28 certified points removed, so it materializes
only genuinely unexplained covers. It contains no point search. After an
external Magma run, parse the complete transcript with:

```sh
python3 elliptic-curves/cas/parse_elkies_2026_rank28_relative_descent.py \
  --program artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.m \
  --log artifacts/local/elliptic-curves/elkies_2026_rank28_relative_2selmer.log \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_rank28_residual_2selmer_magma_v1.json
```

The parser rejects partial or source-mismatched logs. Magma is not available
on the current host. See
[`ELKIES_2026_R17_PAPER_IMPACT_2026-08-27.md`](../elkies-k3/ELKIES_2026_R17_PAPER_IMPACT_2026-08-27.md).

## Bisection specialization controls

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-SPECIALIZATION-CONTROLS 04f49e48e1c1dd88 -->

Replay all 195,600 exact square tests obtained by evaluating the complete
39,120-record rootless bisection atlas at the rank-25--28 controls and at
ICARM curve 394:

```sh
.venv/bin/python \
  elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py \
  --check
```

Every split point is constructed in both square-root branches, verified on
the source and global minimal fibres, and checked against its stored trace.
The exact split counts are `6,3,2,1,25`; their finite-quotient class spans
inside the known public complements have dimensions `5,3,2,1,4`. No point
escapes those spans. Thus the rank-28 fibre exposes one of eleven known
exceptional directions, whereas the `t=3/8` splits recover all four known
directions beyond R17. This leaves the unconditional rank lower bounds
unchanged. Exact relation blocks prove that the displayed generated subgroup
ranks remain `25,26,27,28,21`; they do not bound the full curve ranks. See
[`ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md`](notes/ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md).

Resolve the visible and invisible quotient spaces and test the 2024 rank-29
curve and ICARM 273, 302, and 398--400 against the published fibration:

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES 1c39220ee5fedc77 -->

```sh
.venv/bin/python \
  elliptic-curves/scripts/analyze_elkies_bisection_visibility_and_record_curves.py \
  --check
```

The rank-28 visible span has dimension one and its canonical complement has
dimension ten.  Since the 39,120 atlas records already exhaust bisections
modulo section translation and sign, a higher translated trace shell cannot
enlarge this span.  Exact irreducible degree-24 recognition equations exclude
all six targets from being rational fibres of this published fibration,
including after quadratic twisting.  The rank-28 control recovers the exact
linear factor `5471*t+9529`.  The result does not exclude another fibration,
other families, or isogenies.  See
[`ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md`](notes/ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md).

## Low-conductor exact baselines

### Conductor-first descent targets

Replay the current rank-at-least-21 conductor anchor first:

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve394_rank21.py --check
```

This specializes the compact Elkies R17 family at `t=3/8`, independently
certifies the generic seventeen plus four public directions, and reconstructs
all local conductor exponents.  It proves rank at least 21 at
`log(N)=166.252098...`; no exact-rank statement is made.

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_targets.py --check

.venv/bin/python -m unittest \
  elliptic-curves/tests/test_conductor_first_near_miss_descent.py
```

This pins exact full-dimensional mod-2 known subgroups for ICARM 245, the
Fermigier rank-20 near miss, and both split-infinity rank-19 fibres.  It is a
descent-input certificate, not a complete Selmer result.  Compare the exact
BNF-free factor-base envelopes with

```sh
PYTHONPATH=elliptic-curves/cas \
  /home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/build_conductor_first_s_class_envelopes.py --check
```

This replay orders the custom collectors by materialized Bach/ERH factor-base
size: ICARM 245, family 2, Fermigier, then family 3.  Bach generation is
explicitly ERH-conditional, and the artifact supplies no class-group, Selmer,
or rank bound.  Generate the cheap unconditional dimension-only Magma job
with

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --target icarm-245 \
  --mode selmer-dimension \
  --output artifacts/local/elliptic-curves/icarm245-2selmer-dimension.m
```

If the residual dimension is positive, construct only the quotient covers:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --target icarm-245 \
  --mode relative-covers \
  --output artifacts/local/elliptic-curves/icarm245-relative-2selmer.m
```

An independent unconditional 3-Selmer rank-bound job is generated by replacing
the mode with `three-selmer-dimension`.  None of these modes enables GRH class
group bounds.

For a resource-bounded local PARI diagnostic, with fail-closed model/point
checks and no unconditional use of its provisional upper endpoint, run:

```sh
.venv/bin/python elliptic-curves/cas/run_conductor_first_pari_diagnostic.py \
  --target icarm-245 \
  --output artifacts/local/elliptic-curves/icarm245-pari-diagnostic.json
```

The JSON records wall/RSS limits and raw protocol output.  PARI's cubic-field
BNF is GRH-conditional unless separately certified; only returned rational
points that pass a fresh full mod-2 certificate raise the unconditional lower
bound.

The two exact large-prime collection geometries are implemented in
`run_fermigier_rank20_fixedfb_quadratic_specialq.py` and
`run_fermigier_rank20_minkowski_specialq.py`.  For the latter, the relevant
opt-in switches are

```text
--special-ideal-mode cycle-pairs
--norm-factor-mode exact
--large-prime-merge-mode sparse-hypergraph
--max-residual-primes <cap>
```

The ledger records every retained principal generator plus the exact
large-prime incidence rank and nullity.  This path has not completed a
full-Bach ICARM 245 exact-factor pilot: the pre-factor-hint run reached its
600-second wall limit before collection, and the next setup run was stopped
by the user before the new factor-hint path was replayed.  Do not treat either
terminal as a class-group, Selmer, or rank result.  The bounded parameters and
outcomes are recorded in
[`CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md`](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md).

The independent eclib route can be replayed directly on any integral global
minimal model, for example ICARM 245:

```sh
printf '%s\n' \
  '[1,-1,1,-25880411472355347134118026792,1606663697747901005185875883284420820193259]' \
  | mwrank -q -v 1 -s -x 22
```

On eclib 20231211/20231212 this exits without a bound at the native-integer
conversion in quartic enumeration (`lower bound on c too large`).  The same
failure has been replayed on all four fixed targets; it is a backend resource
failure, not Selmer evidence.

The repository pins the experimental widening patch used to test whether that
conversion is the only ICARM 245 obstruction:

```sh
git clone --depth 1 --branch v20231212 \
  https://github.com/JohnCremona/eclib.git /tmp/eclib-20231212-bigint-c

git -C /tmp/eclib-20231212-bigint-c apply \
  "$PWD/elliptic-curves/cas/eclib-20231212-bigint-quartic-c.patch"

cd /tmp/eclib-20231212-bigint-c
./autogen.sh
./configure --disable-allprogs --enable-mpfp
make -j2

printf '%s\n' \
  '[1,-1,1,-25880411472355347134118026792,1606663697747901005185875883284420820193259]' \
  | timeout --signal=INT --kill-after=5s 30s \
      ./progs/mwrank -q -v 1 -p 256 -s -x 22
```

This patch only changes `c`, `cmin`, and `cmax` to NTL integers and preserves
the exhaustive enumeration.  It passes the former conversion terminal, then
times out in the first Type-3 slice without a Selmer or rank bound.  It is a
backend diagnostic, not a practical descent implementation or mathematical
evidence.

Apply the same gates to the closed 27-fibre anchor-neighborhood pilot:

```sh
.venv/bin/python elliptic-curves/cas/build_conductor_first_family_pilot.py --check

.venv/bin/python elliptic-curves/cas/build_conductor_first_near_miss_magma.py \
  --manifest artifacts/generated-results/elliptic-curves/conductor_first_family_anchor_pilot_v1.json \
  --target family2-u481 \
  --mode selmer-dimension \
  --output artifacts/local/elliptic-curves/family2-u481-2selmer-dimension.m
```

The pilot runs no standalone point search.  Its `family2-u481` certificate
pins two exact points returned incidentally by the provisional descent run;
these raise the unconditional lower bound from 12 to 14 without importing the
uncertified Selmer upper bound.  It leaves nine fibres in the residual-Selmer
queue and the other eighteen at the exact-Tate gate.

See
[`CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md`](notes/CONDUCTOR_FIRST_NEAR_MISS_DESCENT.md)
for all four target identifiers and the strict family-stage queue.
<!-- status-consumer: EC-CF-NEARMISS-DESCENT-INPUTS 25c9f212e5162216 -->

### ICARM curve 245

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve245_rank20.py --check

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/explicit_formula_icarm_curve245_delta22.py --check
```

The first command proves `rank E(Q) >= 20` and independently reconstructs the
exact conductor. The second is a GRH-conditional fixed-fibre upper-bound
diagnostic and does not alter the unconditional status.

### Fermigier rank-20 near miss and `E22`

```sh
.venv/bin/python elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py
.venv/bin/python elliptic-curves/scripts/verify_fermigier_rank_certificates.py
.venv/bin/python elliptic-curves/scripts/verify_benchmarks.py
```

These commands respectively replay the sub-cutoff rank-at-least-20 near miss,
the exact generic-rank/E22 independence certificates, and the family/model
normalization. The literal parameter factor-two discrepancy in the printed
Fermigier source remains open.

### Mestre frontiers

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_dsquare_four_u197.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_mestre_dsquare_rank19_frontiers.py --check
```

The first is an exact rank-at-least-17 certificate. The second checks the two
rank-at-least-19 frontiers and their conditional explicit-formula diagnostics.

## Family and structural certificates

### New six-root family: exact rank 14 at `T=83/6`

This requires SageMath with eclib and PARI:

```sh
sage -python elliptic-curves/cas/newfamily/certify_rank_t83_6.py \
  --efforts 0 \
  --output artifacts/local/elliptic-curves/newfamily/rank_bounds_t83_6.json
```

The exact subgroup rank is 14 and PARI returns the unconditional interval
`[14,14]`. The pinned exact-rank and lower-bound artifacts are deliberately
separate. See [`NEWFAMILY_RANK14_T83_6.md`](notes/NEWFAMILY_RANK14_T83_6.md)
and [`cas/newfamily/README.md`](cas/newfamily/README.md).

### Kihara and Elkies--Klagsbrun baselines

```sh
.venv/bin/python elliptic-curves/scripts/verify_kihara_rank14.py
.venv/bin/python elliptic-curves/scripts/verify_e29_independence.py
```

These replay unconditional rank lower bounds 14 and 29. The public exact-rank
29 statement is conditional and is not used by the second command.

### Exact Nagao certificates

The retained exact/status entry points are:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_section7_picard_bound.py \
  --output artifacts/local/elliptic-curves/nagao_section7_picard_bound.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank17_frontier.py \
  --output artifacts/local/elliptic-curves/nagao_rank17_frontier.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank20_t5081.py \
  --output artifacts/local/elliptic-curves/nagao_rank20_t5081.json
```

The input discovery artifacts for these proofs are preserved and hash-pinned.
The much larger negative search history is archived.

### Fermigier exceptional transport and Mestre two-section geometry

The exact status checkers remain in `elliptic-curves/cas/`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_fermigier_exceptional_transport.py \
  --output artifacts/local/elliptic-curves/fermigier_exceptional_transport.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_fermigier_two_section_generic_rank13.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_diameter235_eight_companion_component.py
```

The first command is an exact transport classification; the separate bounded
point searches retain their bounded label even though their enumeration is
exact.

## Current open computations

- The curve-273 residual 2-Selmer pipeline is under `cas/` with `bnf_free`,
  `residual_selmer`, and `curve273` in the filenames. Its intermediate local
  artifacts remain ignored until a complete certificate exists.
- The current low-conductor searches are the retained Fermigier rank-20,
  denominator-offset, mixed-small-prime, and six-root drivers.
- The H3/rootless-MW17 equation transport lives primarily in `elkies-k3/` and
  has its own reproduction catalogue.

No command in this section turns a partial Selmer calculation, timeout, score,
or bounded negative search into a rank theorem.

## Historical searches

The archived command snapshot and manifest are:

- [`REPRODUCE_2026-08-24.txt`](../archive/elliptic-curves/REPRODUCE_2026-08-24.txt);
- [`MANIFEST.tsv`](../archive/elliptic-curves/MANIFEST.tsv);
- [`archive/elliptic-curves/README.md`](../archive/elliptic-curves/README.md).

The manifest maps every old path to its archive path and records its SHA-256.
Use the historical Git revision named in the archive README when an old script
must be run exactly in its former directory layout.

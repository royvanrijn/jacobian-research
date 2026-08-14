# Reproducing the elliptic-curve track

Run commands from the repository root.  The exact Python layer uses only the
standard library.  Conductor and optional descent calculations require a
working `gp` executable from PARI/GP.

## Environment

```sh
.venv/bin/python --version
gp --version
```

The initial development environment used Python 3.14.6 and PARI/GP 2.17.4.
Later versions are acceptable, but generated JSON records the versions it
actually used.  Sage is not required.

## Fast checks

Compile the modules:

```sh
.venv/bin/python -m py_compile elliptic-curves/cas/*.py
```

Run the exact regression tests:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python -m unittest discover -s elliptic-curves/tests -v
```

The suite includes the Fermigier polynomial and discriminant identities, its
13 visible quartic points and exact Jacobian images, complete multiple-root
lift profiles, fixed-divisor classes, both local root sets used by the simple
pilot, the K3 fixture, the singular-root trap, and the beam-width
counterexample, the arbitrary six-root constructor, and the Nagao replay and
search regressions.  Tests requiring PARI/GP are skipped if `gp` is
unavailable.

## Pinned verified computations

Reproduce the Fermigier family normalization, published minimal model,
conductor, and historical score table:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_fermigier_benchmark.py
```

This writes
[`elliptic_fermigier_benchmark.json`](../artifacts/generated-results/elliptic_fermigier_benchmark.json).
It checks 13 visible quartic points and their 13 Jacobian images exactly, but
cites rather than independently reproves Fermigier's rank-at-least-22 theorem.

Replay all 22 rational points printed by Fermigier:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_fermigier_rank22_points.py
```

This writes
[`elliptic_fermigier_rank22_points.json`](../artifacts/generated-results/elliptic_fermigier_rank22_points.json).
The verifier checks every coordinate on the published curve exactly, then at
96- and 192-digit precision obtains numerical matrix rank 22, positive smallest
eigenvalue `0.064621...`, and determinant `1.299202272...e22`, reproducing the
paper's rounded determinant.  Independently of that height calculation, exact
images of the same twenty-two points in finite quotients
`E(F_p)/2E(F_p)` have combined binary rank 22, and the prime 31 proves trivial
rational 2-torsion.  This is an unconditional exact rank-at-least-22
certificate.  It does not claim an upper bound for the rank, and the curve is
not a target hit: the verifier replays the exact conductor and proves
`ln N>182.72` by a rational upper bound for `e` and a single exact
integer-power inequality.

Replay the record quartic through height one million and search all 28
genus-one slices through its fourteen search-relative accidental abscissas:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_rank22_accidental_slices.py
```

This writes
[`elliptic_fermigier_rank22_accidental_slices.json`](../artifacts/generated-results/elliptic_fermigier_rank22_accidental_slices.json).
The default pass uses slice height 200,000.  It finds five completed
subtarget-conductor parameters but no rank-21 certificate; the strongest
specialized extension, `T=3115/3`, remains stable numerical rank 15 through
quartic height one million.  This is a bounded constructive search, not a
rank upper bound.

Search the ten extra genus-one slices coming from published preimages absent
from the height-one-million source set:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_rank22_missing_preimage_slices.py
```

[`elliptic_fermigier_rank22_missing_preimage_slices.json`](../artifacts/generated-results/elliptic_fermigier_rank22_missing_preimage_slices.json)
records four new decontaminated parameters.  Two completed conductors are
above the target and two calls time out without retry, so no specialized
point search is triggered.

Search the necessary pairwise product curves for two simultaneous published
directions.  The first command is the cheap height-5,000 checkpoint; the
second reaches height 50,000 and therefore includes the record calibration:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_published_pair_fiber_products.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_published_pair_fiber_products_h50000.py
```

The height-50,000 artifact completes all 220 pairs and recovers the record
fiber in every one.  Three extra product-square points fail the two individual
square tests, leaving zero new double-forced fiber and no conductor/rank call.
This is a bounded quotient search, not a proof about larger heights.

Replay the disjoint held-forward-score tranche from the full height-50,000
multiple-root population:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_record_residue_deep_tranche.py \
  --conductor-timeout 15 --point-timeout 45 --saturation-timeout 60
```

[`elliptic_fermigier_record_residue_deep_tranche.json`](../artifacts/generated-results/elliptic_fermigier_record_residue_deep_tranche.json)
records all 23,769 exact population members, the 48-fiber leakage-controlled
selection, 28 completed conductor calls, and the two subtarget fibers
`3206/265` and `1925/157`.  Their maximum stable numerical rank is 14; no
finite-reduction certificate was triggered.

Reproduce the smaller end-to-end K3 fixture:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_ek_k3_fixture.py
```

This writes
[`elliptic_ek_k3_crt_fixture.json`](../artifacts/generated-results/elliptic_ek_k3_crt_fixture.json).
Expected terminal summary:

```text
t=-1468/21 log(N)=148.626486493304... PARI rank bounds=[10,10]
target_hit=false
```

The JSON separates exact rational checks from the PARI/GP minimal-model,
local-reduction, conductor, and 2-descent computations.

### Rank-29 record replay and bounded rank-30 searches

Check the 2024 Elkies--Klagsbrun model and all 29 public points, transport them
to an integral short model, and build the exact finite-reduction certificate:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_elkies_klagsbrun_rank29.py
```

[`elliptic_elkies_klagsbrun_rank29_certificate.json`](../artifacts/generated-results/elliptic_elkies_klagsbrun_rank29_certificate.json)
contains a full-column-rank binary reduction matrix and the separate
`p=67` rational-2-torsion certificate.  It proves `rank >= 29`
unconditionally; the public conditional exact-rank statement is not used.

Run the complete shallow direct-chart tier:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30.py
```

Run the separately pinned deeper version of the same 1,647-chart manifest:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30.py \
  --output artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_search_deep.json \
  --x-pair-height 50000 \
  --x-offset-height 10000000 \
  --slope-offset-height 50000 \
  --slope-pair-height 50000 \
  --chart-timeout 8
```

The default and deep artifacts record 1,647 completed charts, zero timeouts,
and zero nonpublic images.  Search the alternate degree-two coordinates based
at a pinned tranche of exact weight-two/three subgroup sums:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30_alternate_covers.py
```

[`elliptic_elkies_klagsbrun_rank30_alternate_covers.json`](../artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_alternate_covers.json)
records construction of all 4,060 subset-sum covers, deterministic retention
of 64, and 448 completed offset/cross-ratio charts with zero timeout and no
nonpublic image.  These finite searches are negative evidence, not a rank
upper bound.

Run the disjoint signed higher-weight cover tranche:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30_higher_weight_covers.py
```

[`elliptic_elkies_klagsbrun_rank30_higher_weight_covers.json`](../artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_higher_weight_covers.json)
records 2,000 scored signed representatives across five weight bands and 50
completed charts on ten retained covers.  Its 57 exact images all replay in
the public rank-29 subgroup; this is another bounded negative tranche, not an
upper bound.

Finally, replay the strictly capped local PARI descent diagnostic:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/tools/probe_elkies_klagsbrun_rank29_descent.py
```

The pinned run accepted all 29 points and reached the strict 120-second wall
cap without returning an `ellrank` interval.  Its artifact records the timeout
and resource observations; it makes no rank claim.

### Rank-18 source recovery and geometric-rank audit

Replay the official-source availability audit and the Kumar--Kuwata
finite-fibre Galois/lattice computation:

```sh
.venv/bin/python elliptic-curves/cas/audit_elkies_rank18_sources.py \
  --output artifacts/generated-results/elliptic_elkies_rank18_source_audit.json

.venv/bin/python elliptic-curves/cas/audit_kumar_kuwata_f6_galois.py \
  --output artifacts/generated-results/elliptic_kumar_kuwata_f6_galois.json
```

The first command downloads the two pinned official arXiv source bundles if
local copies are not supplied.  The second requires `gp` and records fixed
rank 5 plus quadratic-character ranks `3,2,2,1,1,0,0` for the published
geometric rank-18 basis.  Its status is an exact finite-fibre computational
audit, not a symbolic function-field proof and not a specialization search.
See
[`ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md`](ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md)
for source hashes and archive limitations.

## Automatic Fermigier local discovery and expanded frontier

Discover and classify every compressed local ball over the declared prime
range, select clean split-multiplicative groups by the recorded objective, and
run the bounded CRT search:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_discovered_local_conditions.py
```

[`local_condition_discovery.py`](cas/local_condition_discovery.py) is the
exact library used by this driver.  The output
[`elliptic_fermigier_discovered_local_conditions.json`](../artifacts/generated-results/elliptic_fermigier_discovered_local_conditions.json)
records 44 primes scanned, 188 classified balls, 16 eligible clean split
groups, and the automatic selection at 7, 11, 13, 17 and 19.  The new
13-condition is `T=+2,-2 mod 13`, forcing `v_13(H)>=3`.  Best completed
conductor record is `T=154/103`, `ln N=162.234032455648...`; no rank is
claimed.

Exhaust the original five-group union through projective height 50,000:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/exhaustive_multiple_root_height.py
```

The artifact
[`elliptic_fermigier_multiple_root_height_h50000.json`](../artifacts/generated-results/elliptic_fermigier_multiple_root_height_h50000.json)
records all 23,769 primitive nonsingular sign-orbits, exact local checks, the
leakage-free flow `23769 -> 256 -> 32 -> 12`, and five completed conductors out
of six requests.  Best completed `ln N` is `192.051614237934...`; this is a
bounded negative result with no rank computation.

Compute conductors for the first 24 pinned short-height CRT candidates and
then run their uniform point-pool triage:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/screen_multiple_root_frontier.py \
  --height-count 24 \
  --score-count 0 \
  --timeout 30 \
  --output artifacts/generated-results/elliptic_fermigier_multiple_root_frontier.json
```

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/batch_rank_triage.py
```

The first command writes
[`elliptic_fermigier_multiple_root_frontier.json`](../artifacts/generated-results/elliptic_fermigier_multiple_root_frontier.json):
21 calls complete, with `644/87`, `847/184`, `70/223`, and `1057/218` below
the threshold.  The second writes
[`elliptic_fermigier_batch_rank_triage.json`](../artifacts/generated-results/elliptic_fermigier_batch_rank_triage.json).
At height 50,000 the four low-conductor additions remain at stable numerical
rank 12; only `T=1666/9` is escalated and retains numerical rank 16.  Exact
point membership is checked, but numerical height rank is not an independence
certificate.

## Bounded general six-root survey

Enumerate all affine-normalized integer six-root configurations through
maximum root 14, apply the exact Mestre quartic obstruction, and make one
effort-zero specialization probe for each nonsingular survivor:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/survey_mestre_root_tuples.py
```

[`mestre_root_tuples.py`](cas/mestre_root_tuples.py) supplies the generic exact
constructor.  The output
[`elliptic_mestre_root_tuple_survey.json`](../artifacts/generated-results/elliptic_mestre_root_tuple_survey.json)
records `1023 -> 68 -> 59`; 57 survivors are reflection-symmetric.  The two
nonreflection probes return `[9,9]` and `[4,4]`, and the maximum lower bound
among all 59 probes is 9.  These are bounded software probes, not generic-rank
claims.

## Algebraic Mestre family design, two rank-13 families, and rank 17

Extract Kihara's centered six-root obstruction and its three extra-square
identities:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/derive_kihara_rank14_identities.py
```

Derive the affine-normalized linear-section equations, verify Fermigier's
two-parameter locus, and replay the isolated six-section fiber of the new
root tuple:

```sh
Singular -q elliptic-curves/cas/mestre_affine_section_elimination.sing
Singular -q elliptic-curves/cas/verify_fermigier_affine_section_component.sing
Singular -q elliptic-curves/cas/verify_fermigier_affine_section_jacobian.sing
Singular -q elliptic-curves/cas/verify_mestre_02393128133175_moduli_fiber.sing
```

The universal ansatz is `x=x0+x1*T` with a cubic ordinate.  After the
leading-square condition, exact triangular elimination leaves three residual
equations.  The last command checks that the normalized roots
`(0,23,93,128,133,175)` have exactly the displayed reduced six-point
nonvisible affine-section fiber on the declared open chart.

Replay the generic-rank certificate and degree-40 discriminant geometry with:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python \
  elliptic-curves/cas/verify_mestre_rank13_02393128133175.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python -m unittest \
  elliptic-curves/tests/test_mestre_rank13_02393128133175.py
```

At `u=1`, the exact mod-3 dimensions are 11 for the twelve paired-root
images, 12 after one chosen nonvisible companion, and 13 after split infinity.
The specialization proves generic rank at least 13.  The verifier also proves
that the primitive base-changed discriminant frontier is irreducible,
squarefree, and degree 40; it does not identify every specialized conductor
or claim a rank-17 specialization.

Replay the second D-square family and its conductor-qualified specialization:

```sh
Singular -q elliptic-curves/cas/verify_mestre_02595143168205_moduli_fiber.sing
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/verify_mestre_02595143168205_rank13_section.py
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/verify_mestre_02595143168205_discriminants.py
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/verify_mestre_dsquare_four_u197.py
```

For centers `(0,25,95,143,168,205)`, the first two commands prove generic
rank at least 13 after `T=(39146-u^2)/(2u)`.  At `u=197`, the last command
checks 17 independent points, proves the strict cutoff exactly, and replays
the minimal model and conductor with PARI/GP.  The pinned artifact is
[`elliptic_mestre_dsquare_four_u197_rank17.json`](../artifacts/generated-results/elliptic_mestre_dsquare_four_u197_rank17.json),
SHA-256 `f1235d845653219c53d906a06042d4904686feeb42c379ed7f3d83e01d7f0563`.
The bounded discovery command is:

```sh
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/search_mestre_dsquare_four.py --workers 8
```

The archival discovery replay expects the locally installed `ratpoints`
bundle at `tmp/ratpoints/root/usr/bin/ratpoints`; the tracked certificate does
not.  Its 102 capped conductor calls and finite ratpoints box are not negative
upper-bound evidence.

## Nagao 1994 replay and rank-13 searches

Replay the exact rank-13 base change and Nagao's printed rank-21 curve:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_1994.py
```

This writes
[`elliptic_nagao_1994.json`](../artifacts/generated-results/elliptic_nagao_1994.json).
The verifier checks the 13 and 21 printed points exactly, matches the printed
models, and reproduces `ln N=165.406045732331...` for rank-13 base-change
`u=1` and `ln N=196.679545735892...` for the rank-21 curve.  The independence
theorems are cited from Nagao, not reproved by the numerical determinant.

Classify every polynomial section in the restricted ansatz
`x=m*T+n`, `deg(y)<=3`, and verify the generic group-law relations exactly:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_linear_sections.py
```

This writes
[`elliptic_nagao_linear_sections.json`](../artifacts/generated-results/elliptic_nagao_linear_sections.json).
It recovers the six slope/intercept pairs
`(+/-1/15,703/15)`, `(+/-7/15,928/15)`, and
`(+/-5/3,3628/15)`.  The `+1/15` section is Nagao's known extra section; exact
relations express the other five in the pinned generic basis.  This is an
exact classification only within the stated polynomial ansatz, not among all
rational sections.

Run the staged rare-event search that selects first by unexpected exact point
yield, then by stable numerical height rank, and computes conductors only
after selection:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank13_rank_gain.py
```

The output
[`elliptic_nagao_rank13_rank_gain_search.json`](../artifacts/generated-results/elliptic_nagao_rank13_rank_gain_search.json)
records the complete declared 9,196-parameter population and all stage
settings.  Its maximum stable numerical rank is 17, attained by `u=135/2`,
`42`, and `471/11`, with respectively 66, 33, and 32 unexpected abscissas at
height `10^6`; it finds no numerical rank 18.  The five dependent generic
companions are excluded from unexpected-point counts.  This artifact is a
bounded search, not an exact rank certificate.

The separately pinned mutation run is:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank13_rank_gain.py \
  --farey-denominator 1 \
  --mutation-denominator 64 \
  --mutation-numerator-radius 32 \
  --screen-height 5000 \
  --screen-keep 512 \
  --denominator-diversity 2 \
  --rank-height 50000 \
  --height-keep 256 \
  --rank-keep 20 \
  --final-height 300000 \
  --escalation-height 1000000 \
  --batch-size 64 \
  --batch-timeout 30 \
  --output artifacts/generated-results/elliptic_nagao_rank13_rank_gain_mutations.json
```

[`elliptic_nagao_rank13_rank_gain_mutations.json`](../artifacts/generated-results/elliptic_nagao_rank13_rank_gain_mutations.json)
records all 5,133 declared parameters, again reaches numerical rank 17 but not
18, and adds `u=74` with 20 unexpected height-`10^6` abscissas to the
frontier.

Replay the later leakage-controlled rare-event calibration gate:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank21_rare_event_model.py
```

[`elliptic_nagao_rank21_rare_event_model.json`](../artifacts/generated-results/elliptic_nagao_rank21_rare_event_model.json)
stores leave-one-positive-out scores for four exact rank-at-least-18/19 fibers,
three-fold cross-fitted controls, and a held-out published rank-21 calibration
fiber.  Neither predeclared model clears the recovery gate, so the script
correctly launches no broad candidate, conductor, or point search.  This is a
negative model-validation result, not evidence against the target curves.

Certify the four rank-17 specializations and replay their conductors directly
from the rational short models:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank17_frontier.py
```

[`elliptic_nagao_rank17_frontier_certificate.json`](../artifacts/generated-results/elliptic_nagao_rank17_frontier_certificate.json)
contains exact `E(F_p)/2E(F_p)` matrices of full column rank 17, exact
rational-2-torsion checks, point membership, and direct minimal-model,
conductor, and root-number replays for `u=135/2`, `471/11`, `42`, and `74`.
It proves `rank >= 17` unconditionally for all four.  It also proves the
strict conductor comparison without floating logarithms: every exact
conductor is below `10^66`, and the recorded positive Taylor partial sum proves
`ln(10)<231/100`, so `ln(N)<7623/50=152.46<182.72`.

The smaller, independent `u=42` certificate replay is:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_u42_rank17.py
```

[`elliptic_nagao_u42_rank17_certificate.json`](../artifacts/generated-results/elliptic_nagao_u42_rank17_certificate.json)
stores an `18 x 17` exact binary matrix of rank 17 and a trivial rational
2-torsion certificate at `p=31`.  The proof is elementary finite reduction;
it does not rely on numerical heights, BSD, parity, full 2-descent, or the
finite-index premise of `ellsaturation`.

Run the leakage-free integer scans through 200 and 2,000:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank13_integer_u.py
```

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank13_integer_u.py \
  --u-bound 2000 \
  --stages 200,2000,20000 \
  --keep-counts 100,25,15 \
  --score-timeout 600 \
  --conductor-timeout 35 \
  --stack-bytes 512000000 \
  --output artifacts/generated-results/elliptic_nagao_rank13_integer_u2000.json
```

These write
[`elliptic_nagao_rank13_integer_u.json`](../artifacts/generated-results/elliptic_nagao_rank13_integer_u.json)
and
[`elliptic_nagao_rank13_integer_u2000.json`](../artifacts/generated-results/elliptic_nagao_rank13_integer_u2000.json).
The smaller scan supplies `u=42`, `u=84`, and the eliminated `u=50`; the
larger scan retains only `u=1256` and `u=42` below the strict conductor target
among completed finalists.  Scores and generic rank never populate
`target_hits`.

Discover clean local balls directly in the base-change variable, combine all
declared symbols by CRT/Gauss reduction, and triage the two integer survivors:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank13_local_crt.py
```

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/triage_nagao_rank13_local_candidates.py
```

[`nagao_rank13_local.py`](cas/nagao_rank13_local.py) is the exact local library.
The first command writes
[`elliptic_nagao_rank13_local_crt.json`](../artifacts/generated-results/elliptic_nagao_rank13_local_crt.json),
covering 2,048 symbol choices at 7, 11, 13, 19 and 31.  The second writes
[`elliptic_nagao_rank13_local_candidate_triage.json`](../artifacts/generated-results/elliptic_nagao_rank13_local_candidate_triage.json):
at height 50,000, `u=118` has 10 nonvisible images and stable numerical rank
15, while `u=316` has one nonvisible image and numerical rank 13.  Their
respective conductors have `ln N=128.027255994266...` and
`177.107241730402...`.

Extend `u=118` once to height `10^6`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/extend_nagao_u118_height.py
```

[`elliptic_nagao_u118_height_1000000.json`](../artifacts/generated-results/elliptic_nagao_u118_height_1000000.json)
records 43 nonvisible exact images but stable numerical rank 15 at 72 and 120
digits.

Run exact-section and bounded point-pool triage on the integer-scan finalists:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/triage_nagao_rank13_finalists.py
```

The artifact
[`elliptic_nagao_rank13_finalist_triage.json`](../artifacts/generated-results/elliptic_nagao_rank13_finalist_triage.json)
stores exact membership for every retained point and two-precision height
matrices.  `u=42` reaches stable numerical rank 17 and `u=84` rank 16 through
height `10^6`; `u=50` returns PARI effort-zero computational bounds `[13,13]`.

The final capped `u=42` checkpoint is reproduced separately by

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/extend_nagao_u42_frontier.py
```

and writes
[`elliptic_nagao_u42_height_10000000.json`](../artifacts/generated-results/elliptic_nagao_u42_height_10000000.json).
It records 17 exact checked points returned by small-prime `ellsaturation`,
height-determinant ratio `2^32` (consistent with within-span index `2^16`),
stable numerical rank 17, an effort-zero rank timeout after 60 seconds without
bounds, and the sole height-`10^7` point-search timeout after 120 seconds.
These bounded failures and the numerical height calculation are not an upper
bound.  The exact rank lower bound is supplied separately by the
finite-reduction certificates above, not by `ellsaturation`'s finite-index
premise.

Search nonuniform height regions for new `u=42` points:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_u42_skew_height.py
```

[`elliptic_nagao_u42_skew_height.json`](../artifacts/generated-results/elliptic_nagao_u42_skew_height.json)
records ten disjoint reduced numerator/denominator boxes through denominator
128,000 and 76 determinant-one Möbius charts at transformed height 50,000.
The run finds 40 abscissas outside the old uniform height-`10^6` checkpoint.
Every returned point is checked exactly on the quartic and Jacobian, and its
relation in the pinned certified rank-17 subgroup is replayed exactly.  The
augmented pool retains numerical height rank 17 at 72 and 120 digits.  This is
a bounded negative search result, not a rank upper bound.

Run the matching search on the smallest-conductor certified frontier curve,
`u=135/2`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_u135_skew_height.py
```

[`elliptic_nagao_u135_skew_height.json`](../artifacts/generated-results/elliptic_nagao_u135_skew_height.json)
records the same ten skew boxes and 76 determinant-one Möbius charts.  The
boxes find 74 exact points outside height `10^6`, with maximum absolute
numerator 43,277,563 and maximum denominator 107,700; every point has an
exactly replayed relation in the certified rank-17 subgroup, with coefficients
of absolute value at most 8.  The charts add no outside point, and the stable
numerical height rank stays 17.  This bounded negative result motivates trying
alternate 2-covers before extending the same quartic boxes further; it is not
a rank upper bound.

Generate the pinned Magma programs locally for inspection:

```sh
python3 elliptic-curves/tools/generate_u42_magma.py --mode verify
python3 elliptic-curves/tools/generate_u42_magma.py --mode rankbounds
python3 elliptic-curves/tools/generate_u42_magma.py --mode twoselmer
python3 elliptic-curves/tools/generate_u42_magma.py --mode twodescent
```

The generator writes the selected program to standard output and changes no
files.  The historical results are recorded in
[`elliptic_nagao_u42_magma_probe.json`](../artifacts/generated-results/elliptic_nagao_u42_magma_probe.json):
Magma V2.29-9 exactly verified the 17 points and their independence, but the
anonymous public calculator reached its observed 311.34 MB memory limit on
the upper-bound paths.  No Selmer group, covers, or rank interval was returned.
The commands above only generate auditable input; they are not instructions to
submit work to an external service.

Replay the bounded local PARI diagnostic under explicit time and memory caps:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/tools/run_u42_pari_rank.py \
  --timeout 600 \
  --stack-bytes 8000000000 \
  --rss-limit-bytes 8000000000 \
  --heartbeat 45
```

The expected recorded outcome is a strict timeout (exit status 124) without
an `ellrank` interval; the runner writes no files.  PARI 2.17.4 reached about
1,895,056 KiB recorded RSS, while its largest reported allocated stack was
1.024 GB.  The combined bounded PARI, eclib/mwrank, Sage-launcher, and Magma
diagnostics are preserved in
[`elliptic_nagao_u42_descent_toolchain.json`](../artifacts/generated-results/elliptic_nagao_u42_descent_toolchain.json).
Every upper-bound route ended at a recorded resource or implementation limit,
so this artifact supplies no rank upper bound and no mathematical evidence for
one.

## Nagao rank-21 record neighborhood

Search the CRT/Gauss neighborhood that preserves five local features of
Nagao's published rank-21 curve, then triage its best low-conductor neighbor:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank21_neighborhood.py
```

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/triage_nagao_rank21_neighbor.py
```

The outputs are
[`elliptic_nagao_rank21_neighborhood.json`](../artifacts/generated-results/elliptic_nagao_rank21_neighborhood.json)
and
[`elliptic_nagao_rank21_neighbor_triage.json`](../artifacts/generated-results/elliptic_nagao_rank21_neighbor_triage.json).
Among 110 bounded parameters, constructor `T=6041/198` has
`ln N=170.765123121845...`, but it does not beat the published record's final
score.  Its height-50,000 search finds no new Jacobian sign-pair, the visible
and augmented pools both have numerical rank 11, and effort-zero `ellrank`
times out.  For calibration, the published specialization's visible pool also
has numerical rank 11 while its full printed 21-point set has numerical rank
21 and cited published independence.

## Bounded Hensel--CRT--Gauss pilot

The default pilot forces two split multiplicative discriminant roots to square
moduli, uses separate good-reduction residues at 7 and 11, and performs PARI
conductor work on a bounded candidate subset:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_crt_lattice.py \
  --split-only \
  --power-primes 89:2,131:2 \
  --rank-primes 7,11 \
  --beam-width 128 \
  --coefficient-radius 12 \
  --score-bound 200 \
  --score fermigier-good \
  --keep 30 \
  --pari-count 1 \
  --pari-timeout 60 \
  --output artifacts/generated-results/elliptic_fermigier_crt_lattice_pilot.json
```

The output must say `target_hits=0` unless a candidate has both a recorded
PARI rank lower bound at least 21 and `ln N < 182.72`.  The default command does
not request rank work, so a low-conductor candidate alone cannot become a hit.

To request PARI 2-descent for the first selected candidate, add

```text
--pari-rank-count 1 --pari-rank-effort 0
```

This can be substantially slower.  `ellrank` effort greater than zero includes
a randomized point search, so exact reproduction then also requires retaining
the output, not just the command.

Useful matched controls are:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_crt_lattice.py \
  --power-primes '' \
  --rank-primes 7,11 \
  --pari-count 1 \
  --output artifacts/generated-results/elliptic_fermigier_rank_only_control.json
```

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_crt_lattice.py \
  --split-only \
  --power-primes 89:2,131:2 \
  --rank-primes '' \
  --pari-count 1 \
  --output artifacts/generated-results/elliptic_fermigier_power_only_control.json
```

These runs are **bounded experiments**.  Changing the beam width, coefficient
radius, prime sets, scoring cutoff, or PARI subset changes the enumerated
search and must produce a separately named artifact or an explicitly recorded
replacement.

## Simple-root power-pair control

Enumerate every distinct-prime pair of clean simple roots modulo `p^2` for
`5<=p<=500`, reconstruct bounded lattice neighborhoods, score the lowest 256,
and compute two conductors:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_power_pairs.py \
  --pari-count 2 \
  --pari-timeout 100 \
  --score fermigier-good
```

This writes
[`elliptic_fermigier_power_pairs.json`](../artifacts/generated-results/elliptic_fermigier_power_pairs.json).
The pinned run covers 2,084 root-choice pairs and 4,165 unique nonsingular
representatives.  The two conductor calls give `ln N=239.682887...` at
`T=147/2210` and `269.592630...` at `T=4947/877`; neither meets the conductor
bound, and no rank was computed.

## Multiple-root CRT experiment

Run all 144 combinations of the five cheap local-symbol groups, enumerate the
declared reduced-basis boxes, score the lowest-height pool, and make four PARI
conductor calls:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_multiple_root_crt.py \
  --score fermigier-good
```

This writes
[`elliptic_fermigier_multiple_root_crt.json`](../artifacts/generated-results/elliptic_fermigier_multiple_root_crt.json).
The expected conductor-only success is

```text
T=70/223
v_7, v_11, v_17, v_19, v_37 of H(T) = 18, 5, 5, 4, 3
ln N = 165.979271710943...
target hits = 0
```

The run checks the five local conditions against PARI's minimal model.  It
does not call `ellrank`; the twelve supplied Jacobian seeds and their numerical
height determinant do not turn this into a rank claim.

## Record-residue-class scan

Exhaust every primitive `a/b` of projective height at most 5,000 in
`a=2142*b mod 4403`, the cheap multiple-root class containing Fermigier's
record, and compute conductors for the top two scores:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_record_residue_class.py \
  --pari-count 2 \
  --pari-timeout 30 \
  --score fermigier-good
```

This writes
[`elliptic_fermigier_record_residue_class.json`](../artifacts/generated-results/elliptic_fermigier_record_residue_class.json).
The expected first record is

```text
T=1666/9
score through p<=500 = 40.048807038898...
v_7, v_17, v_37 of H(T) = 18, 4, 3
ln N = 128.959882907388...
target hits = 0
```

The scan is exhaustive only in its declared height box.  No rank call was made;
the conductor inequality and high heuristic score are insufficient.

Verify `T=1666/9` independently, including its exact specialization, twelve
point equations, conductor factorization, and local data at every bad prime:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_fermigier_1666_9.py
```

This writes
[`elliptic_fermigier_1666_9.json`](../artifacts/generated-results/elliptic_fermigier_1666_9.json).
It reproduces the exact conductor
`101523255017246417712694892860237179024105368632978033830`,
`ln N=128.959882907388...`, and `v_7(H),v_17(H),v_37(H)=18,4,3`.
It does not invoke `ellrank` and certifies no target hit.

Search the associated quartic for extra rational points:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_extra_points.py
```

The bounded artifact
[`elliptic_fermigier_extra_points.json`](../artifacts/generated-results/elliptic_fermigier_extra_points.json)
records height bound `10^6`, 114 signed points, 57 distinct `x`-values, and 44
new `x`-values beyond the thirteen visible sections.  Exact mapping checks pass.
Both 96- and 192-digit height calculations give numerical rank 16 on the same
subset with one-based indices
`[1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17]`; at 192 digits its determinant is
`3.4790082779825878e23` and smallest eigenvalue is `0.999542...`.  This is not
an exact rank certificate.  The artifact also records a prior twelve-seed
`ellrank` attempt that exhausted a 1 GB PARI stack after about 300 seconds
before returning bounds; that failed call is not rerun by this command.

Compare the same good-prime statistic through growing numerical cutoffs:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/compare_score_cutoffs.py
```

This writes
[`elliptic_fermigier_score_cutoffs.json`](../artifacts/generated-results/elliptic_fermigier_score_cutoffs.json).
At `B=500,2000,10000,100000`, `T=1666/9` scores respectively
`40.048807...`, `60.519323...`, `81.634378...`, `122.707444...`; `E22` scores
`40.913185...`, `69.525178...`, `106.746181...`, `163.165765...`.  The ten
record candidates were selected at `B=500`, so comparisons at that cutoff are
selection-biased.

Run the leakage-free staged height-5,000 rescore:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/staged_record_rescore.py \
  --conductor-count 0
```

The artifact
[`elliptic_fermigier_record_rescore_h5000.json`](../artifacts/generated-results/elliptic_fermigier_record_rescore_h5000.json)
records `5520 -> 50` at `B=2000`, `50 -> 10` at `B=10000`, and `10 -> 10`
at `B=100000`.  `T=1666/9` is not a finalist.  Final leader `T=1547/492`
scores `124.536543...` at `B=100000`, versus `163.165765...` for separately
evaluated `E22`.  With `--conductor-count 0`, the canonical run computes no
conductors or ranks and reports no hit.

## Nagao section-7 generic theorem and rank-20 fiber

Classify the generic polynomial sections through abscissa degree five, prove
the twelve-section mod-3 independence certificate, and verify the K3 fiber
geometry:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_section7_linear_sections.py
```

Then count the good reduction over `F_29` and `F_(29^2)` and reconstruct the
residual Frobenius factor:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_section7_picard_bound.py
```

Expected output includes `#S(F_29)=1212`, `#S(F_29^2)=723600`, arithmetic
generic rank exactly 12 over `Q(T)`, and geometric generic rank in `[12,13]`.
The artifacts are
[`elliptic_nagao_section7_linear_sections.json`](../artifacts/generated-results/elliptic_nagao_section7_linear_sections.json)
and
[`elliptic_nagao_section7_picard_bound.json`](../artifacts/generated-results/elliptic_nagao_section7_picard_bound.json).

Replay the exact rank-20 specialization at constructor `T=5081/47`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank20_t5081.py
```

This verifies twenty exact points, their finite-reduction independence,
`ln N=174.249816228548...`, and root number `+1`.  The result is rank at least
20, not a target hit.  The complete `2^20-1` cover-class search and the
follow-up extreme skew pass are substantially more expensive bounded runs:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank20_t5081_direction.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank20_t5081_cover_skew.py
```

They find 224 decontaminated points in total, all with exactly replayed
relations in the certified rank-20 subgroup.  Finally, reproduce the
conditional explicit-formula diagnostic:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/explicit_formula_rank20_t5081_delta22.py
```

Its conservative bound is below 22.  With root number `+1`, this gives
analytic rank at most 20 under GRH; the algebraic exact-rank conclusion also
uses BSD and is explicitly conditional.

Construct auxiliary curves through the accidental points on that fiber:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank21_accidental_slices.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_section7_remaining_auxiliary_slices.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_section7_accidental_genus2_slices.py
```

The two higher-yield genus-one slices `a04` and `a09` have a separate exact
group-orbit replay:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_section7_auxiliary_group_orbits.py \
  --gp /private/tmp/pari-map-src.33iJSU/pari/Odarwin-aarch64/gp-dyn
```

This command uses a locally pinned PARI 2.19 development binary because the
required `ellfromeqn` map is absent from the system PARI 2.17.  The artifact
records version commit `776ac5d`, binary SHA-256
`3bff0db14041b12b1af88ecd13b73ba09829c3abdde5ed9bd2b3112b368d7f88`,
the exact maps, and all pullbacks.  The temporary binary path is not portable;
if it is absent, the generated artifact remains the replay checkpoint and the
development build must be reconstructed before rerunning the command.

The first command classifies the `x=m*T+n` slices.  The second exhausts its
declared doubled-point ternary boxes on the remaining genus-one Jacobians.  The
third completes all 240 genus-two slices at height 5,000 and finds the exact
subtarget-conductor parameters `T=163` and `T=1049/10`; their specialized
height-50,000 screens remain numerical rank 12.  The optional deeper genus-two
extension is intentionally capped and reports seven timeouts, not a negative
height-one-million enumeration:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/extend_nagao_section7_a10_genus2.py
```

## New rank-19 and rank-18 specialization certificates

Replay the exact rank-19 certificate at constructor `T=6793/64` and the
historical-finalist rank-18 certificate at `T=6629/174`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank21_t6793.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank21_historical_finalists.py
```

The first gives exact rank at least 19, `ln N=158.572648489303...`, root
number `-1`.  The second gives exact rank at least 18 at `T=6629/174`,
`ln N=154.795114152374...`, root number `+1`.  The exhaustive primitive box
that found the additional exact rank-18 fibers `3137/72` and `5783/16` is:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_rank21_unbiased.py
```

The section-7 global scan is much larger: it exhausts 18,244,819 primitive
parameters in its declared `30000 x 1000` box before its staged point work.
Run it only when that cost is intended:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_nagao_section7_global.py
```

Its exact new rank-17 leaders are `T=599/2`, `ln N=124.061012256948...`, and
`T=426`, `ln N=138.825822315291...`; both have root number `-1`.  None of
these commands reports a rank-21 target hit.

## Expanded family and terminal bounded searches

The root-tuple census commands compile their pinned C++ enumerators and then
run the declared conductor/point tranches.  They are bounded searches, not
quick unit tests:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_mestre_root_tuple_scale.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_mestre_root_tuple_scale_max100.py
```

The second command exhausts the normalized maximum-root-100 census and its
pinned 64-fiber point tranche.  It finds roots `(0,6,47,55,70,80)`, `T=8`.
Replay the resulting exact rank-at-least-13 certificate independently with:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_mestre_0647557080_t8_rank13.py
```

This checks all thirteen exact points, full finite-reduction rank 13, absence
of rational 2-torsion, `ln N=82.351544058010...`, and root number `-1`.

The other new family has separate global and prime-power lanes:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_mestre_0430313946_frontier.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_mestre_0430313946_power_crt.py
```

The first exhausts 18,244,819 primitive positive rational parameters in its
declared box and contains an exact rank-at-least-12 certificate at `T=5`.
The second builds 144 exact CRT classes from complete `p^4` profiles through
199.  Its sole completed subtarget finalist, `T=209001/3868`, stays at
numerical rank 10 through height one million.

The following searches are materially longer and should only be rerun when
their complete bounded populations are intended:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_global.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_fermigier_rank22_record_group_triples_remainder.py
```

They respectively enumerate the 60,815,684-member Fermigier global box and
finish the 5,761 weight-three record-group directions left after the pilot.
Neither produces a target hit.  The second command writes an append-only JSONL
stream as well as its summary artifact.

Finally, reproduce the two disjoint rank-29 denominator sieves with:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30_denominator_sieve.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/search_elkies_klagsbrun_rank30_companion_center_sieve.py
```

Together they cover exactly 55,267,250,510 primitive rational abscissas and
find no thirtieth point.  Their conclusion is confined to those two declared
regions.

## Fermigier `E22` conductor benchmark

This PARI/GP replay minimizes the published model and prints its conductor and
natural logarithm:

```sh
gp -q <<'GP'
default(realprecision, 60);
E=ellinit([1,0,1,-940299517776391362903023121165864,10707363070719743033425295515449274534651125011362]);
Em = ellminimalmodel(E);
N = ellglobalred(Em)[1];
print(N);
print(log(N));
quit
GP
```

Expected output begins

```text
22720638514787473197194583889675055980109503436060704437972911338086049759883790
182.724910950637428796
```

This GP snippet verifies the conductor benchmark, not 22-point independence.
Fermigier's paper proves that theorem, and the earlier repository verifier now
also supplies an independent exact finite-reduction certificate.

## Historical score-table compatibility

The following replay uses the published `E22` model.  Here `M` is a prime
ordinal: the loop omits 2 and reaches the `M`-th prime.

```sh
gp -q <<'GP'
default(realprecision, 30);
E=ellinit([1,0,1,-940299517776391362903023121165864,10707363070719743033425295515449274534651125011362]);
S=0.;p=2;
for(n=2,10000,p=nextprime(p+1);ap=ellap(E,p);S+=(2-ap)/(p+1-ap)*log(p);if(n==50||n==100||n==200||n==400||n==1000||n==2000||n==4000||n==10000,print(n," ",p," ",S)));
quit
GP
```

Rounded score values should be

```text
29.49 44.12 57.54 81.51 105.17 122.76 143.84 166.47
```

matching Fermigier's table.  The second printed column records the actual last
prime (`229, 541, 1223, 2741, 7919, 17389, 37813, 104729`) and exposes the
difference from a literal cutoff `p <= M`.

## Reading a JSON run

Each search artifact records:

- the explicit target and its strict inequality;
- status as a bounded experiment;
- the family and external source;
- every search parameter and the reproducing command;
- Python and PARI/GP versions;
- chosen CRT constraints and actual forced valuations;
- rational parameter, height, and local traces; and
- PARI output or a timeout/error for the bounded subset.

Before citing a candidate, verify that `target_hits` contains it and inspect the
underlying `pari_ellrank.lower_bound`; a high score or
`below_log_conductor_target: true` is insufficient.

## Additional reproduction commands

Run commands from the repository root.  The exact lattice code uses only the
Python standard library.  The curve replay additionally requires PARI/GP;
version 2.15.4 produced the pinned manifest.

### Unit tests

```bash
python3 -m unittest discover -s elliptic-curves/tests -v
```

These include simple and singular Hensel lifts, generalized CRT, skew Gauss
reduction, the primitive-vector cancellation trap, complete asymmetric height
boxes, exact Fermigier quartic-to-Weierstrass identities, and finite-reduction
independence certificates.

### Benchmark arithmetic

```bash
python3 elliptic-curves/scripts/verify_benchmarks.py
```

This recomputes the exact conductor of the Fermigier E22 benchmark and the
literal integer cutoff.  It also pins the unresolved source normalization:
the printed shift `19754/39` gives a different model and conductor, while the
doubled shift and canonical adapter reproduce E22 exactly.  It does not
itself reproduce point independence; the exact lower-bound checker below does.
No command here supplies an unconditional rank upper bound.

Cross-check both family metadata files against the executable equations and
the calibration family's rank-two specialization with:

```bash
python3 elliptic-curves/scripts/verify_family_data.py
```

### CRT--lattice calibration

Replay the checked-in manifest:

```bash
python3 elliptic-curves/scripts/verify_crt_lattice_calibration.py
```

Generate an unpinned fresh copy in the ignored local cache:

```bash
python3 elliptic-curves/scripts/run_crt_lattice_calibration.py \
  --output artifacts/local/elliptic-curves/crt_lattice_calibration.json
```

The generator refuses to overwrite an existing file.  The pinned artifact was
created with:

```bash
python3 elliptic-curves/scripts/run_crt_lattice_calibration.py \
  --output artifacts/generated-results/elliptic-curves/crt_lattice_calibration_v1.json
```

Its SHA-256 is
`eb1543031e68026042c921ee2b93e765070b65340b8129b74f0629a9b3d5c8fa`.

### Fermigier high-family CRT seed

Replay the checked-in three-prime seed:

```bash
python3 elliptic-curves/scripts/verify_fermigier_crt_seed.py
```

Generate a fresh copy in the ignored cache:

```bash
python3 elliptic-curves/scripts/run_fermigier_crt_seed.py \
  --output artifacts/local/elliptic-curves/fermigier_crt_seed.json
```

The pinned command uses output
`artifacts/generated-results/elliptic-curves/fermigier_crt_seed_v1.json`.
Its SHA-256 is
`a4f2e27d63bbf2160cb8afaed1b171295bf941e99ac8db8f3d2bb85424edaf0c`.
The replay certifies the exhaustive CRT/height result and local reduction at
89, 131, and 137.  It intentionally does not compute a global conductor or a
rank certificate.

### Fermigier rank evaluator and certificates

Replay the reconstructed thirteenth point, the twelve independent generic
section differences, and all 22 independent published E22 points:

```bash
python3 elliptic-curves/scripts/verify_fermigier_rank_certificates.py
```

The replay uses exact finite-field arithmetic and independently checks every
stored group order, generator order, and discrete-log equality with PARI/GP.
The pinned artifact SHA-256 is
`94fc64d7f1744f6a20a0396d32914cd36330107db2538e03ee95cc3e32927051`.
Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_fermigier_rank_certificates.py \
  --output artifacts/local/elliptic-curves/fermigier_rank_certificates.json
```

Evaluate another adapter parameter, optionally run PARI's bounded quartic
point search, and certify a modularly independent subset of the point cloud:

```bash
python3 elliptic-curves/scripts/evaluate_fermigier_specialization.py 19754/39 \
  --quartic-height 100000 --certify-searched-subset \
  --output artifacts/local/elliptic-curves/e22_evaluation.json
```

The `hyperellratpoints` height limit is a bounded search, not a completeness
claim.  The exact certificate proves only the rank of the selected subset.
For larger batches, an installed `ratpoints` executable can replace PARI and
apply a denominator cap:

```bash
python3 elliptic-curves/scripts/evaluate_fermigier_specialization.py 3251/16 \
  --search-engine ratpoints --quartic-height 2000000 \
  --denominator-bound 13000 --certify-searched-subset \
  --output artifacts/local/elliptic-curves/3251_16_evaluation.json
```

`ratpoints` is optional and is not vendored or installed by this repository.
An existing quiet abscissa-only output can be replayed without rerunning the
search by replacing the search options with
`--ratpoints-output artifacts/local/elliptic-curves/POINTS.out`.

### Rank-20 low-conductor near miss

Replay the pinned 58-abscissa search output, exact 20-point independence
certificate, global minimal model, and conductor:

```bash
python3 elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py
```

Regenerate with an installed `ratpoints` 2.1.3:

```bash
python3 elliptic-curves/scripts/run_fermigier_rank20_near_miss.py \
  --output artifacts/local/elliptic-curves/fermigier_rank20_near_miss.json
```

The generator runs height `2000000` with denominator bound `13000`.  To replay
an already captured quiet output instead, add `--ratpoints-output PATH`.  The
pinned artifact SHA-256 is
`8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1`.
It is a rank-at-least-20 near miss, not a target solution.

### Staged Fermigier score sweep

Build and run the dependency-light C++ ranking pass with:

```bash
g++ -O3 -march=native -fopenmp -std=c++20 \
  -o /tmp/fermigier-score-sweep \
  elliptic-curves/ecsearch/fermigier_score_sweep.cpp
OMP_NUM_THREADS=32 /tmp/fermigier-score-sweep 100000 500 20000 \
  > artifacts/local/elliptic-curves/fermigier_score_sweep.tsv
```

The three integer arguments are the maximum numerator, maximum denominator,
and retained output count.  The score is a search heuristic, not a rank or
conductor computation, and the command intentionally writes only to the
ignored local cache.

### Kihara rank-14 baseline

Replay the exact specialization and independence certificate with:

```bash
python3 elliptic-curves/scripts/verify_kihara_rank14.py
```

Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_kihara_rank14.py \
  --output artifacts/local/elliptic-curves/kihara_rank14_t2.json
```

The pinned artifact SHA-256 is
`851ff6da6ccf4f4dca947048edd43846ff7da41161e83fde419747e715a0df46`.
This is a rank-at-least-14 family baseline, not a rank-30 candidate.

### Public rank-29 baseline

Replay all 29 published points and their exact finite-reduction independence
certificate with:

```bash
python3 elliptic-curves/scripts/verify_e29_independence.py
```

Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_e29_independence.py \
  --output artifacts/local/elliptic-curves/elkies_klagsbrun_e29.json
```

The pinned artifact SHA-256 is
`a585a8bc081c67fc6314b7be8ea29721b465fcd8f147d170b534ecb52395891e`.
It proves the public lower bound 29 locally and exactly; it neither supplies a
thirtieth point nor replays the conditional upper bound.

### Combined gate

```bash
make verify-elliptic-curves PYTHON=python3
```

The repository's default `.venv/bin/python` is not present in every checkout;
the explicit override is sufficient for the dependency-free Python layer.
The curve and local-reduction replays additionally require `gp` on `PATH`.

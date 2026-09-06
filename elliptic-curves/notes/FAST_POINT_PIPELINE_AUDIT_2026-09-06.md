# Faster pointed searches and complete fixed-curve follow-ups

The coordinate control exposed two avoidable costs: the point-search engine
and repeated finite-prime admission work. Both now have separately versioned
replacements with exact retained witnesses. The new rank-at-least-26 curve
has completed all 301 fixed adaptive PARI boxes; its full 3,756-point cloud
still certifies 26 modulo 2, 3 and 5. The small-conductor curve's 301 PARI
boxes and all archived histories also replay; its 25,068-point cloud still
certifies 22 modulo 2, 3 and 5. All standalone replays pass. No new rank-at-least-28/32 curve is claimed.

## Backend calibration and the infinity correction

The [previous native control](EXACT_PARITY_AND_COORDINATE_AUDIT_2026-09-06.md)
verified identical affine square-coordinate sets on 49 fixed boxes. On the
new rank-26 curve, the comparison now passes on another **43 boxes and 334
affine square coordinates**. PARI reports 14.147 seconds of search CPU time.
These are finite comparisons, not universal completeness or speed theorems.

The first new-curve comparator failed on a projective infinity entry because
it attempted `Fraction(1,0)`. Its source, partial artifact and failure log
remain intact. Version 2 excludes denominator zero from the affine comparison;
infinity remains separately checked by exact evaluation of the leading
quartic coefficient. Four tiny independently enumerated quartic boxes also
verify inclusive numerator and denominator height boundaries in the pinned
PARI executable, despite the loose wording of its short help.

A stricter calibration replayer now also checks the complete ordered chart
rosters, raw programs, return codes and hit-count fields. It rejects empty,
duplicated or reordered records. The old square-set checks did not enforce
all of those fields; their original artifacts remain unchanged. Both
calibration rosters and all four boundary boxes pass the stricter replay.

The optional [PARI backend](../cas/pari_pointed_backend.py) consumes an exact
shared pointed chart and a frozen reduced model. It checks the two horizontal
maps, their composition, the raw quartic identity and rational square scale.
PARI performs only `hyperellratpoints` on the fixed `P,Q` model. Raw program,
stdout, stderr, executable hash, limits and return code are retained. A
completion marker is required. Every returned coordinate is checked against
the final integer quartic; both square-root signs are mapped exactly back to
the curve. Infinity is checked independently in Python.

These records identify a **separate backend**. They do not pretend to be GMP
transcripts. A timeout has no asserted denominator prefix; only its retained
raw output and independently computed infinity witnesses survive. A replay
checks maps, squares and points without calling either search engine. Finite
coverage still trusts the pinned PARI execution. Focused tests reject truncated
output and out-of-box coordinates and cover a nontrivial point at infinity.

## Prime banks and archived observation history

The old admission method starts an ambiguous point from a small certificate
prime set, then repeatedly extends the whole old point basis through the
remaining primes. That work is repeated for the next ambiguous point.

An immutable first-chart benchmark on the small-conductor curve compares the
same ordered **81 point admissions**. Preloading all available primes through
997 once preserves every basis, rank, observation status and finite relation
mask. Existing admission took 19.472 seconds; preload took 0.819 seconds and
subsequent admissions took 2.782 seconds. A separate exact final-rank check
passes. This is a measured comparison on one retained chart, not a general
speed guarantee. Tests additionally include a new independent point,
dependence and a known point with the opposite sign.

The [preloader](../cas/research_runtime/preloaded_prime_state.py) preserves the
rational basis and records every added or excluded prime. State keys and
finite-certificate presentations change explicitly. Separately, the
[observation rotation](../cas/research_runtime/rotated_observation_state.py)
returns the **entire old state for archival** before producing the next live
state with an empty observation list. The continuation writes that archive
before searching. Each new admission remains in its chart's observations;
the next archive, or the final state, retains it.

This changes bookkeeping only. In the hash-pinned `MWState.adjoin` source,
observations enter appended history and parent keys, while the admission
choices use the preserved model, basis, finite reductions, torsion witness
and geometry fields. The raw points, all archived observations and exact
parent-state keys remain available. Tests check preserved decisions and
archive restoration; the continuation replayer reconstructs and compares
every archive and admission. No frozen worker or earlier artifact was edited.

## Rank-26 follow-up

The 301 centres are the original fixed adaptive words in its certified
26-point subgroup. PARI maps were constructed before the new search. The
first backend run reached its 600-second outer limit with **285 complete
boxes**; every retained chart and the complete point cloud replayed.
A 45-chart immutable prefix had also replayed while that worker was running.
That snapshot was not treated as a terminal result or a restart signal.

A separately frozen prime-bank continuation searched only indices 285–300.
All **16 remaining boxes** completed in 13.454 seconds and their map,
archive and admission replay passed in 2.271 seconds. The combined roster
contains each of the 301 fixed maps exactly once, all at height 100,000.
The original metric attempts, their tails, the 43 initial PARI boxes and
this adaptive PARI follow-up give **989 retained source charts**. Their
point-only concatenation contains 3,756 unique points up to sign and still
certifies rank at least 26 through prime 997.

A separate odd-modulus audit checks the same complete cloud, rather than
adding a quotient-escape count to a baseline. Both moduli 3 and 5 give full
finite column rank 26 with exact no-rational-ell-torsion witnesses. All
selected finite signatures and independent columns replay. Streaming row
reduction is tested against independent dense elimination. These results
supply lower bounds only; none proves exact rank or saturation.

## Small-conductor follow-up

The curve at `3/17` retains its previously certified rank-at-least-22 bound
and exact 76-digit conductor. The experiment uses all 301 previously fixed
adaptive centres in its 22-point subgroup, with a stop at the first certified
increment pending independent replay. There are no new parameters or public
exceptional points in the search.

After the admission benchmark, the slower run was intentionally terminated
at 822.832 seconds. The stop decision was recorded before signalling the
verified owned process group. Its terminal checkpoint contains **39 complete
charts**, and all retained maps, admissions and points replay. This is an
explicitly censored experiment, not a completed 301-chart run. Its original
protocol and the in-flight uncheckpointed chart were not silently rewritten.

The new continuation begins at index 39, retaining the exact same remaining
maps. Its **262 boxes** completed in 551.671 seconds under a separately frozen
600-second / 1.5-GiB cap, with worker lower bound 22. All previous observations
and subsequent chart histories are archived. All 262 exact map/admission and archive histories replayed in 313.005
seconds. Together with the 39-chart prefix, every fixed map appears once.
The initial 43 charts, older 127-chart metric attempt and new 301 PARI boxes
give 471 source charts and a complete 25,068-point cloud up to sign. Its
mod-2 certificate and standalone replay pass, still at 22. The mod-3 and
mod-5 audits and standalone replays also give 22.

## Reproduction

```sh
python3 elliptic-curves/cas/compare_new_rank26_search_engines_v2.py --check
python3 elliptic-curves/cas/check_pari_height_boundaries.py --check
python3 elliptic-curves/cas/continue_fixed_pari_search.py replay --case rank26
python3 elliptic-curves/cas/audit_retained_cloud_modl.py \
  --check artifacts/generated-results/elliptic-curves/rank26_all_retained_modl_v1.json
```

The source point-cloud concatenations are explicitly point-only inputs,
not synthetic chronological transcripts. Earlier evidence remains in its
original bundles. Mathematical status remains in `MATH_STATUS.json`.

The [evidence supplement](../../artifacts/generated-results/elliptic-curves/fast_point_pipeline_evidence_v1.json) names four pinned base archives. Its [isolated verifier](../cas/verify_fast_point_pipeline_bundle.py) checks all 602 new adaptive geometries and raw point sets, all 278 continuation admission/archive histories, both complete clouds modulo 2, 3 and 5, calibration rosters and targeted tests. The original slow 285+39 admission histories passed separately and are retained, without claiming they were rerun in the isolated check.

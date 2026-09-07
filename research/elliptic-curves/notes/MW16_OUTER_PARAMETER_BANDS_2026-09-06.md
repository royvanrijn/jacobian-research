# New MW16 fibres beyond compact parameter height

**Completed: two new curves with rank at least23. All2580 point boxes,
local proofs and182 isolated replay stages pass.**

| Local ID | MW16 family | Parameter | Certified lower bound |
|---|---|---|---:|
| new-20260906-190 | a1-fibration-02 | -5383/2146 | 23 |
| new-20260906-191 | a1-fibration-03 | 9053/2 | 23 |

Both equations are unmatched among593 pinned catalogue and917 prior equations.
The [point certificates](../../artifacts/generated-results/elliptic-curves/outer60_mw16_results_v1.json)
and [complete V16 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v16.json)
are available. All191 inventory point proofs and the equation CSV replay.
Exact ranks, records and universal novelty remain open.

| Band | Completed curves | Certified lower bounds |
|---|---:|---|
| 4096 < H <=16384 | 20 | 11 at16,5 at17,1 at18,1 at20,2 at23 |
| 16384 < H <=65536 | 20 | 17 at16,2 at17,1 at18 |
| 65536 < H <=262144 | 20 | 19 at16,1 at17 |

Every retained-cloud bound agrees modulo2,3,5. The higher-band results measure
this fixed point exposure, not exact ranks or the absence of higher-rank
fibres. The corrected57,439,050-byte standalone point bundle passes all182
isolated stages. Its v1 verifier used a leftover directory name and failed
before its first logical proof stage; v1 evidence is preserved, and v2 fixes
the path. The bundle does not rerun the broad scan or15360 scalar scores.

One new23-point curve is504th of512 in its signed short-score list and1016th
of1024 after the two signs merge. This motivates the
[fresh deeper-retention trial](MW16_FRESH_OUTER_RETENTION_2026-09-06.md),
which stays beyond4096 and uses disjoint denominator slices.

The five compact MW16 families have broad completed parameter scans through
height4096. Their million-height follow-ups increase point-search height on
already selected fibres; they do not enlarge this parameter population.
The user requested priority for new parameter territory and new fibres.

For primitive nonzero n/d with d positive, put H=max(abs(n),d).
This campaign uses three disjoint bands:

| Parameter band | Denominator modulus | Signed slices | Retained candidates |
|---|---:|---:|---:|
| 4096 < H <=16384 | 16 | 10 | 5120 |
| 16384 < H <=65536 | 256 | 10 | 5120 |
| 65536 < H <=262144 | 4096 | 10 | 5120 |

Each band gives one slice to each family and sign. SHA256 of a frozen literal
salt, band, family and sign chooses the denominator residue, with five odd
and five even residue classes per band. The inner square is excluded before
heap admission. These are deterministic stratified populations, not complete
outer squares or uniform samples. Parameter height depends on this fixed
family coordinate; new parameters are not by themselves new Q-isomorphism
classes.

The [independent model audit](../../artifacts/generated-results/elliptic-curves/mw16_outer_parameter_models_v1.json)
recomputes all15360 equations by homogeneous integer arithmetic, checks every
primitive parameter and band assignment against the scanner, and finds
**15360 different j-invariants**. Thus these retained models are mutually
nonisomorphic even over an algebraic closure. This compares the new pool with
itself; comparison with prior work remains post-terminal. Displayed integral
coefficient sizes have medians301,325 and349 bits across the three bands;
these are not claims about globally minimal heights. The audit build and
independent command replay pass.

All30 slices complete, covering **286812899 primitive addresses**. Every one
of the15360 retained562-prime scores replays against the canonical projective
trace tables. Primitive counts independently replay by inclusion-exclusion.
Forty small signed/frame tests compare the complete reference population,
annulus filtering and top-seven ordering; all pass. One predetermined slice
per band passes a20-second runtime gate before the remaining27 slices.
The scanner has two workers,60 seconds per slice and a900-second outer cap.

The generic mathematical input is the existing exact16-section atlas and its
successful target-free prospective point detector. The experiment tests
candidate incidence at new heights; it does not assume that larger height
improves rank or that existing point visibility extrapolates without loss.
All five families have equal retention and finalist budgets.

The frozen extended-score protocol evaluates all15360 survivors at primes
4099..65521, with two independent direct character sums per curve. The first15
rows cover every band/family once and pass the fixed runtime gate, projecting
1921.070 seconds serial against4800 allowed. Two workers,20 seconds per curve,
a3600-second outer cap and30-address checkpoints apply. Raw calls are
immutable and failures/censoring are retained without automatic retry.

After exact score replay, select four Q-isomorphism-distinct curves per
family per band,60 total, by combined quantized S1 through65521, good count,
denominator and signed numerator. Cross-group deduplication uses fixed band
then family order and refills only from the same frozen1024-address group.
Wholly disjoint65537..131071 validation follows the frozen selection and never
changes it. Known-record equations, parameters, points, ranks, target j-values
and jump labels stay outside selection and prospective execution.

All15360 extended scores pass exact replay, including30720 direct character
sums. The60 finalists have distinct equations with no deduplication skips;
all180 additional direct sums in their wholly disjoint validation pass.
Their parameter heights range5383..15634 in band1,23175..62635 in band2,
and67511..260076 in band3. None of these scores is a rank bound.

The separately frozen point controller passed all selection, validation and
map gates before completing60 generic16-only
attempts,43 recorded generic parity labels per curve, point height125000 and
ten seconds per chart. Every map file must precede any points. This does not
claim exhaustive specialized parity coverage. Catalogue and prior-equation
comparison follows the terminal point batch and independent exact proofs.

Sources: `../cas/scan_mw16_outer_bands.py`,
`../cas/newfamily/scan_rational_nagao_annulus.cpp`,
`../cas/finish_mw16_outer_scan.py`, `../cas/score_mw16_outer_bands.py` and
`../cas/finish_mw16_outer_scores.py`, `../cas/finish_outer60_mw16_points.py`
`../cas/finalize_outer60_mw16_results.py` and
`../cas/audit_mw16_outer_parameter_models.py`. Protocols, raw calls and ledgers live
under `artifacts/local/elliptic-curves/mw16-outer-bands-v1` and
`mw16-outer-band-scores-v1`. Heuristic scores are not rank lower bounds;
no bounded miss establishes point absence, exact rank or an upper bound.

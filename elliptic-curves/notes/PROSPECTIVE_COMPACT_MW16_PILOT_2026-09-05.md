# Prospective searches on five compact MW16 families

A later [fixed twelve-address follow-on](NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md)
adds four further curves, including a rank-at-least-22 curve with an exact
76-digit conductor. The completed pilot and wider experiment below retain
their original cohort counts and protocols.

The [five compact inputs](COMPACT_FIVE_MW16_ATLAS_2026-09-05.md) now have a
separate, frozen prospective experiment. Its selector and point workers read
no catalogue equations, known record parameters, target invariants, public
points, rank labels or jump labels. Catalogue comparison is a post-batch
operation. The completed experiment adds ten curves of independently
certified rank at least 22–25 absent from the pinned 584-equation catalogue
and twenty-one earlier certified discoveries. Its best new curve has rank
at least 25, a gain of nine beyond its sixteen generic sections. Together
the first six certificate exports contain thirty-one distinct curves with
lower bounds at least 22, including two at least 25. The completed wider
extension below adds a third rank-at-least-25 curve, giving thirty-two
distinct curves across seven exports. No new curve reaches
the rank-at-least-28 near-record or rank-at-least-32 record target.

## New certified curves

| Anonymous family | New compact parameter | Rank at least |
|---|---|---:|
| `a1-fibration-05` | `307/206` | 25 |
| `a1-fibration-01` | `787/103` | 24 |
| `a1-fibration-05` | `509/149` | 24 |
| `a1-fibration-05` | `-227/647` | 24 |
| `a1-fibration-03` | `478/49` | 23 |
| `a1-fibration-01` | `32/253` | 22 |
| `a1-fibration-01` | `699/19` | 22 |
| `a1-fibration-04` | `-338/47` | 22 |
| `a1-fibration-04` | `-209/118` | 22 |
| `a1-fibration-04` | `-20/23` | 22 |

A convenient integral equation for the rank-at-least-25 example is

```text
y^2 = x^3 + x^2
 - 178583656613609913632588264590135324867787225*x
 + 910711367539034324827860968784624386747127915714766155601502590839.
```

From its short certificate coordinates `(X,Y)`, set `x=X-1/3`, `y=Y`.
All twenty-five certified points transport to this equation. No global
minimality, exact rank or conductor assertion is made.

The [standalone certificate](../../artifacts/generated-results/elliptic-curves/prospective_mw16_results_v1.json)
retains all twenty measurements, including the catalogue match and smaller
lower bounds. Its Sage-free checker recomputes point membership, full column
rank in a product of finite mod-2 quotients, absence of rational 2-torsion,
all generic section transports and exact rational-isomorphism comparisons.
Every integral relation must then have even coefficients, giving independence
by infinite descent. All twenty certificate replays passed. The checker
shares group-law primitives with discovery; external or formal verification
is not claimed.

```sh
python3 elliptic-curves/cas/certify_prospective_mw16_results.py --check \
  artifacts/generated-results/elliptic-curves/prospective_mw16_results_v1.json
```

Post-batch comparison finds one known equation: `a1-fibration-05,-34/87`
matches ICARM 548 and recovered lower bound 24. All nineteen other equations
are absent from the pinned catalogue and earlier twenty-one certificates;
none of the twenty are rationally isomorphic to each other. The ten new
rank-at-least-22 examples have distinct exact `j`-invariants, also distinct
from those of the earlier twenty-one. These finite comparisons do not prove
absence from every publication or unpublished computation. The snapshot is
dated 2026-09-05, SHA-256
`7e80549befa11a07422a3960967f4cd80264d8675cb3e0a99f0c9c5afb340f72`.
This work used ICARM, supported by NSF Grant DMS 2425401.

## Completed batch and limits

All twenty selected fibres certified the full sixteen specialized generic
points. Nineteen workers completed 43 charts. `a1-fibration-02,39/10` hit
its 300-second limit after 36 charts and retains certified lower bound 21.
Across all twenty retained measurements the lower-bound counts are
`16:3, 20:3, 21:3, 22:5, 23:1, 24:4, 25:1`; one of the four rank-24
measurements is the known catalogue curve.

All 853 retained chart/admission records passed exact replay without
resieving. The 36-chart transcript needed a separately declared 300-second
retry after its initial 120-second replay allowance expired. None of the
853 chart boxes reached denominator 100000 within the four-second allowance.
Completing nineteen chart plans therefore does not mean their boxes were
exhaustive. The [evidence manifest](../../artifacts/generated-results/elliptic-curves/prospective_mw16_pilot_evidence_v1.json)
and [portable bundle](../../artifacts/generated-results/elliptic-curves/prospective_mw16_pilot_evidence_v1.zip)
preserve both successful checks and censored attempts.

A separately frozen adaptive follow-up on the new `307/206` rank-25 curve
uses only its discovered and independently certified points. It lifts the
first 301 generic parity classes from the fresh census through nonzero
nine-bit quotient words, then orders representatives by specialized numerical
height. Its allowance is one 1500-second/1.5-GiB worker, 301 charts at the
same height and four-second chart limit, stopping at lower bound 28.
The original worker stopped at 281 charts. A separately declared 300-second
continuation completed the remaining twenty frozen charts, with no further
gain. All 301 chart/admission records passed exact replay. Protocols are in
`artifacts/local/elliptic-curves/prospective-mw16-rank25-followup-v1/`.
Separate mod-3/mod-5 audits of the four fresh rank-24/25 initial point
clouds retained ranks 24/24, 24/24, 24/24 and 25/25. They checked 615, 455,
450 and 375 points respectively, up to sign, for `787/103`, `509/149`,
`-227/647` and `307/206`. All four certificate replays passed; none proves
dependence or a rank upper bound.

The [diagnostic manifest](../../artifacts/generated-results/elliptic-curves/prospective_mw16_followup_diagnostics_v1.json)
and [portable bundle](../../artifacts/generated-results/elliptic-curves/prospective_mw16_followup_diagnostics_v1.zip)
retain the final adaptive transcript, allowances, failed gates and checks.
They also preserve the bounded conductor work: generic discriminant-factor
splitting leaves composite cofactors of 318 and 309 bits for the first two
new rank-25 curves, beyond the declared 192-bit factorization gate. An exact
local screen of the then-current 31 new curves computes conductor upper
bounds from primes through 10000 and the remaining discriminant cofactor.
None proves a conductor below its listed ICARM rank-threshold benchmark.
Testing known repeated polynomial factors does not sharpen those bounds on
the ten new MW16 inputs. These are inconclusive upper bounds; exact
conductors remain unknown and no conductor record is established.

## Completed balanced height-4096 extension

The positive rank-25 yield opened a separately frozen wider population gate.
All five families scored 20400078 signed primitive parameters each, totalling
102000390 addresses, with the identical 562-prime tables. No generic census
or trace computation was repeated. The complete input bytes were reused and
hashed. The wider box contains the earlier box and is not an independent
population. Scanners were capped at 180 seconds each, with two family
selection workers capped at 400 seconds/1 GiB.

Four finalists per family were frozen. Two exactly repeated addresses reused
their immutable initial measurements without refilling. All eighteen new
workers completed their 43-chart plans under the original point limits;
all 774 new chart/admission records passed replay. All twenty consolidated
point certificates passed independent replay. The complete roster has lower
bound counts `16:10, 17:4, 18:1, 20:1, 24:2, 25:2`, including reused results.

One additional curve is new to the pinned catalogue and all 41 equations
retained in the earlier certificates: `a1-fibration-04` at `-1647/91` has
unconditional rank at least 25. A convenient integral equation is

```text
y^2 = x^3 + x^2
 - 229492257703712367397273830141349998161108*x
 + 43151908014471055892737153120239384952325513367104943738588288.
```

Again `x=X-1/3`, `y=Y` transports all twenty-five short-model points to this
equation. Its exact `j`-invariant is distinct from those of the earlier
thirty-one rank-at-least-22 curves. The other new rank-25 measurement,
`a1-fibration-04,-1905/52`, matches ICARM 542. The repeated `-34/87` matches
ICARM 548; repeated `787/103` is our earlier rank-24 discovery. Among the
seventeen noncatalogued, nonrepeated new workers, only `-1647/91` has a
retained lower bound at least 22. No failed gain is a rank upper bound.

The [wider certificates](../../artifacts/generated-results/elliptic-curves/prospective_mw16_wide_results_v1.json)
replay with `certify_prospective_mw16_wide_results.py --check FILE`.
The [manifest](../../artifacts/generated-results/elliptic-curves/prospective_mw16_wide_evidence_v1.json)
and [portable evidence](../../artifacts/generated-results/elliptic-curves/prospective_mw16_wide_evidence_v1.zip)
preserve the wider protocols, source, populations, reused measurements and
fresh transcripts. The ignored source directory is
`artifacts/local/elliptic-curves/prospective-mw16-h4096-v1/`.

## What remains missing

The wider all-prime score still drops a strong already measured fibre:
`a1-fibration-05,307/206`, rank at least 25 with score `93.223654229337`,
falls below that family's new fourth-place cutoff `93.662759529635`.
Its three new higher-scoring finalists returned lower bounds 16, 17 and 17.
This demonstrates a retention/visibility limitation of the finite experiment,
not an ordering of the curves' true ranks. The subsequent
[retained-point audit](RECORDED_POINT_ADMISSION_AUDIT_2026-09-05.md) resolves
the apparent ICARM 542 loss: its original stored points already certify 26
when the finite admission budget reaches prime 257. The worker's prime-251
cutoff had certified only 25. This diagnosis uses no public oracle points
in the new rank proof.

The new `-1647/91` curve has now completed its own frozen 301-chart adaptive
follow-up. The initial 1,500-second allowance stopped at 217 charts; a
600-second continuation reached 290, and a final 180-second allowance
completed the remaining eleven. Earlier checkpoints and all supervision
logs remain preserved. All 301 chart/admission records passed independent
replay. Re-certifying its 3,739 retained points up to sign through prime 997
still gives lower bound 25, without an upper-rank conclusion.

A separate fixed experiment sampled 2,048 new full specialized rank-25
parity classes on each of the three new rank-25 curves. All 6,144 integral
representatives, their rounded-form norms and deterministic masks passed
exact replay. None met the predeclared gate requiring a 5% improvement in
the top-43 median proposed norm over the earlier 301-class pool. No extra
point searches were launched from this sample, and no CVP optimality or
rank exclusion follows. The
[follow-up diagnostic manifest](../../artifacts/generated-results/elliptic-curves/rank25_followup_diagnostics_v2.json)
and accompanying ZIP retain the completed adaptive transcript, all replay
logs and the full parity samples.

All 32 current rank-at-least-22 discoveries are available in one
[JSON index with points and exact rank certificates](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v1.json)
and [equation CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v1.csv).
`export_new_high_rank_curve_index.py --check FILE` replays every rank proof,
all 32 distinct exact `j`-invariants and the pinned catalogue comparison.

## Frozen scope and mathematical gate

The exact generic rank-16 markings and positive definite height Grams have
determinant 474. Each Gram is scaled by its exact common denominator before
enumerating all 65536 parity classes. The 43 nonzero classes with greatest
computed representative norm, breaking ties by mask, supply that family's
pointed-chart classes. Floating CVP chooses representatives; exact norm and
parity checks do not prove numerical optimality or a covering-radius theorem.
The [exact census audit](../../artifacts/generated-results/elliptic-curves/prospective_mw16_census_audit_v1.json)
checks all 327680 representatives and the retained class sets.

Each family scores 1275854 signed primitive rational parameters `n/d`, with
`1<=abs(n),d<=1024`, using all 562 primes 5 through 4093 before selection.
Zero and infinity are outside this finite population. Every good local trace
contributes `round(10^12*(2-a_p)*log(p)/(p+1-a_p))`. Integer total score, good
prime count, denominator and signed numerator determine retention. The final
four scores are strictly above either scanner shard's retention boundary;
the negative shard's opposite numerator tie convention cannot affect these
four. The duplicate prime-5 H band required by the scanner is unused and is
not validation. Sage cardinalities independently checked all 200 projective
rows at primes 5, 7, 11 and 13 across the five families.

Every family has its own fresh parity census. The geometry helper is
dimension-independent and loads no record fixture. Exact diagonal-Gram
regressions at dimensions 16, 25 and 28 passed. All eighty generic point
transports through the normalization used by the worker passed a separate
fixed-input check before the point protocol was frozen.

The point batch runs all twenty fixed addresses in round-robin family order,
without catalogue filtering, replacement or a rank-based batch stop. Each
worker first specializes all sixteen sections and certifies their independence
by finite mod-2 quotients. A deficient certificate retains all input points
and the certified subset and censors this chart plan; it is neither a rank
upper bound nor a permanent exclusion. Successful inputs receive 43 charts
at height 100000, four seconds per chart. At most four point workers run
concurrently, each capped at 300 seconds and 1.5 GiB. A worker stops early
only at certified lower bound 32. The fixed batch still finishes its other
addresses. Complete chart plans do not imply exhaustive boxes.

This experiment measures candidate incidence under one score and conditional
point visibility under a finite chart policy. It does not settle global
solubility of unobserved covers or an upper bound on any fibre.

## Reproduction

All protocols, software bindings, tables, full generic censuses, populations,
supervision logs and chart checkpoints live in
`artifacts/local/elliptic-curves/prospective-mw16-h1024-v1/`.

- `select_prospective_mw16_atlas.sage`: frozen selection protocol, per-family
  census and per-family score scan; two workers, 300 seconds and 1 GiB each.
- `search_prospective_mw16_atlas.sage`: frozen point protocol and one worker.
- `run_prospective_mw16_batch.py`: the fixed twenty-address supervisor.
- `replay_prospective_mw16_search.py RESULT.json`: Sage-free chart and rank
  admission replay, without resieving.
- `certify_prospective_mw16_results.py`: post-batch standalone finite rank
  certificates and exact rational-isomorphism catalogue comparison. It refuses
  to read the catalogue until all twenty batch attempts are terminal.

Never overwrite a frozen protocol or change a running worker's source. A
continuation requires a new explicit finite allowance and preserved earlier
logs. All timeouts and memory failures remain visible in the point ledger.

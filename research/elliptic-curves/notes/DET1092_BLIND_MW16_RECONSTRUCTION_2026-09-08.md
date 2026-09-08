# Determinant1092: blinded MW16-core reconstruction

Authority: `EC-DET1092-BLIND-MW16-RECONSTRUCTION-20260908`.

**CORE_DEPENDENT_16_TO_17_RECOVERY: 9 of 17 frozen leave-one-out cores recover a rational 17th generic direction.** Every promoted claim passes two serializations of an independent algebraic replay. This answers the positive-control question affirmatively for the successful cores. It does not discover a new MW17 surface or increase the known generic rank.

Starting repository commit: `5257d634a390b3ce6d65790d93a7817cfea74b3b`. The checkout already had uncommitted research changes. Their exact bytes, the parent, governing instructions and background notes are bound in the [frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/protocol.json) and copied without alteration under `frozen_sources/`. The initial worktree record is retained. No historical certificate was overwritten.

## Fixed population and blinding

The [immutable roster](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/roster.json) contains exactly the 17 coordinate omissions, numbered **one-based**. Each `fixtures/arm-NN.json` contains the same original a-invariants, precisely 16 retained rational sections, their principal16-by16 Gram, and the same declared policy. Numerators and denominators use ascending rational coefficient lists. No omitted row/column pairings, full-basis words, specializations or exceptional points are supplied.

Fixture preparation, isolated reconstruction, the all 17 terminal barrier, full-basis evaluation and the later specialization diagnostic are separate stages. Each fresh worker runs behind a Linux Landlock filesystem allowlist and in separate network/PID namespaces. It can read its own fixture, its frozen worker source and the CAS/system runtime; the full parent, other fixtures and outputs, repository research data, personal caches and network are inaccessible. A read-denial probe is part of every arm. The designer read the requested historical background; the claim is audited algorithm/process input blindness, not human ignorance of the known parent.

Historical orbit8044, historical winning centres, saved RR nets, cover parameters and exceptional witnesses are **not** reconstruction inputs. The equation itself is expressly allowed even though it was historically recovered using the full parent. Every worker verifies the 24-I1 K3 height hypotheses directly from that equation, all section identities, the exact retained Gram and positive-definite rank16.

## One common reconstruction policy

1. Reduce the retained Gram using a verified unimodular LLL transformation. Generate 512 deterministic SHA256 parity masks, with at most 128 strictly norm-decreasing coordinate/signed-pair moves per mask. These are bounded heuristic representatives, **not certified coset minima**.
2. Retain at most 32 distinct parities of heights 8 through 14, ordered by decreasing height, original-word L1 size and lexicographic word. The complete sampled ledger and selected words are saved before any RR factorization in that arm.
3. For each retained centre `T`, of height `h`, solve the exact system `f0+f1*x(T)+f2*y(T)=0` with degree bounds `(n,n-4,n-6)`, first at `n=ceil(3h/4)+1`, then at `n+1`. Inspect the rational RREF kernel basis rows in their returned order. There are no target-fitted linear combinations.
4. Remove `T` from each cubic-line intersection and test the residual quadratic for splitting over **Q(t)**. Retain both distinct rational roots of each split member. Test independence by the exact 17-section height Gram, stopping after the first successful member and retaining both its roots.

Every arm has the same 180-second CPU cap (185 hard), 240-second wall cap, 8 GiB address-space cap and 128MiB per-file cap. At most 64 RR systems are attempted. The stop after success censors the remaining selected centres, degree bounds and kernel rows; it is not a completed negative exposure. No failed arm was replaced, refilled or rerun, and no policy was changed after the roster freeze.

The mathematical motivation is the effective residual divisor

```text
D = 3O+nF-T = 2O+(n-h/2)F-phi(T),
D^2 = 4n-3h-8.
```

A rational splitting can expose section classes absent from the retained lattice, even though the divisor and its RR equations were constructed from that lattice. The method tests this actual splitting, rather than inferring a point from an abstract lattice class. A nonsplit residual curve, including a soluble nonconstant base cover, gives no primary success.

## Seventeen terminal outcomes

| Omitted | Outcome | Recovered subgroup rank | Quotient | Wall seconds | Observed geometry |
|---|---|---:|---|---:|---|
| P1 | bounded miss | 16 | no gain | 90.991 | 64 RR stages; det(core)=1965 |
| P2 | bounded miss | 16 | no gain | 108.794 | 64 RR stages; det(core)=1853 |
| P3 | success | 17 | recovered, coefficient 1 | 10.442 | centre 4, h=14, n=12; det(core)=4392 |
| P4 | success | 17 | recovered, coefficient -1 | 40.670 | centre 14, h=14, n=12; det(core)=2400 |
| P5 | bounded miss | 16 | no gain | 96.821 | 64 RR stages; det(core)=1916 |
| P6 | success | 17 | recovered, coefficient 1 | 13.048 | centre 5, h=14, n=12; det(core)=3876 |
| P7 | success | 17 | recovered, coefficient 1 | 10.092 | centre 4, h=14, n=12; det(core)=2420 |
| P8 | bounded miss | 16 | no gain | 106.895 | 64 RR stages; det(core)=1805 |
| P9 | bounded miss | 16 | no gain | 92.304 | 64 RR stages; det(core)=1776 |
| P10 | success | 17 | recovered, coefficient 1 | 60.168 | centre 22, h=14, n=12; det(core)=1853 |
| P11 | success | 17 | recovered, coefficient 1 | 4.076 | centre 2, h=14, n=12; det(core)=3356 |
| P12 | success | 17 | recovered, coefficient -1 | 9.891 | centre 4, h=14, n=12; det(core)=2933 |
| P13 | bounded miss | 16 | no gain | 99.073 | 64 RR stages; det(core)=1461 |
| P14 | success | 17 | recovered, coefficient 1 | 4.276 | centre 2, h=14, n=12; det(core)=8237 |
| P15 | bounded miss | 16 | no gain | 95.454 | 64 RR stages; det(core)=2621 |
| P16 | success | 17 | recovered, coefficient -1 | 12.346 | centre 5, h=14, n=12; det(core)=1944 |
| P17 | bounded miss | 16 | no gain | 94.053 | 64 RR stages; det(core)=1445 |

Totals: **627 RR stages**, **224 returned-section records, deduplicated within each arm**, **0 failed RR stages**, **0 wall timeouts**; worker wall time summed over the sequential dispatch is **949.393s**, and CPU time is **940.055s**. These include startup and input certification. They are observed costs on a shared host, not dedicated-machine benchmarks. Every candidate, including every dependent one, is recorded in the per-arm certificates; all completed nonsplitting members are retained. Exact resource samples and stop/censor locations are in the stage ledgers. Administrative synthetic/runtime probes and their failures are retained separately under `preflight/`.

The [full-basis evaluation](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/evaluation.json) was run only after `primary_complete.json` sealed all 17 terminal outputs. For every returned section it solves the exact full height system, verifies the integral group-law identity, records all 17 integer coefficients and checks the omitted coefficient. Any nonzero residual height or nonintegral word would fail closed as an inconsistency with the saturated full MW17 theorem.

For example, omitting P3 produces a certified quotient representative with ordered full-basis word

```text
(-1,0,1,1,-1,1,0,0,0,1,-1,-1,0,0,0,0,0)
```

Its rational-function coordinates, RR line, retained-centre word and exact Gram certificate are in the corresponding arm certificate. Both signs of a quotient representative and translations by retained sections count as recovery; equality to the literal omitted basis point is unnecessary.

Recovered candidate relations: `{"DIFFERS_BY_RETAINED_SUBGROUP": 9, "NEGATIVE_QUOTIENT_PLUS_RETAINED_SUBGROUP": 9}`.

The same recovered pair appears independently in the P12-omission and P16-omission arms, with nonzero coefficients in both missing directions. This is a repeated alternative quotient representative across different cores; it was not injected from one worker into the other.

The audit trail for each candidate is: **frozen equation + retained 16 + recomputed Gram → retained-lattice centre → exact RR kernel row → rational residual factor → exact section → core height independence**. The later full-basis word is an evaluation endpoint, never an input. This trail is witnessed numerically and symbolically by the per-arm coefficients, source/fixture hashes, restricted workers and the frozen stage barrier.

A structural discriminator is visible after evaluation: every successful first member has centre height 14 and degree n=12, hence residual arithmetic genus 0. Its two rational components have a sum in the retained core but individually leave its rational span. Core determinant alone does not completely separate outcomes: P2 and P10 both have retained determinant 1853, but only the latter recovers within the fixed budget.

## Independent replay and evidence layout

The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/replay.json) uses separately written affine chord/tangent formulas and pole-intersection heights, without importing the producer/evaluator or using Sage elliptic group operations. It verifies every returned section and full-basis identity, the rank16 inputs, every recorded RR-line identity, every residual factor/square decision and all promoted ranks. For completed RR systems, rational independent kernel rows give an upper rank bound and exact finite-field matrix ranks give the matching lower bound. The latter certify rational **matrix rank**, not elliptic independence.

The [compact replay](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/compact_replay.json) independently regenerates omitted expanded polynomials and checks their canonical digests before repeating the exact claims. The full raw output bytes are preserved without editing under `research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1/runs/`, as required by the repository raw-evidence policy. `raw_archive.json` binds those bytes to the original all 17 barrier. Compact `arms/` records retain all RR lines, candidates, stage outcomes, source dependencies and resource ledgers; no unsuccessful member is silently discarded.

A first independent-replay attempt was deliberately censored after seven passing arms because it redundantly checked every internal group-operation equation and formed unused ordinates. Its source, log, timing and censor reason are preserved. The corrected checker retains every final equation, group identity, rank and RR proof gate; no primary arm was rerun.

The pre-run infrastructure failures were a refused WSL bind mount, a sign typo in a synthetic toy chord and a duplicate basename during source copying. All occurred before any primary reconstruction, are preserved, and did not change any arm outcome. Linux Landlock replaced the unavailable mount isolation before the protocol freeze.

Replay from the repository root:

```sh
sage -python research/elliptic-curves/cas/det1092_blind_mw16/replay.py --compact
```

This verifies the portable compact package without local raw files. The original full replay command omits `--compact` and uses the preserved local raw records. Existing certificates are compared without overwriting them. The original `prepare.py` and `run.py` refuse to overwrite the frozen primary experiment; reproducing a new search execution requires a separately named package and must not be merged into these17 outcomes.

## Secondary specialization and the selector question

Only after the generic terminal freeze and independent replay, the [secondary diagnostic](../../artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1/diagnostic.json) evaluates recovered sections at t=0. Each is checked against its exact full-generic word on curve302 and mapped to the established 31-point embedding by an exact integer vector. No exceptional point coordinates or point-cloud searches are used. The canonical injective specialization of MW17 supplies the rank/quotient boundary.

- arm-03: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 415 bits (original 381).
- arm-04: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 465 bits (original 381).
- arm-06: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 423 bits (original 381).
- arm-07: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 433 bits (original 381).
- arm-10: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 450 bits (original 381).
- arm-11: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 396 bits (original 381).
- arm-12: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 431 bits (original 381).
- arm-14: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 138 original classes and changes 15; median raw quartic coefficient size is 386 bits (original 381).
- arm-16: the replacement is unimodular and preserves the entire integral MW17 specialization. Its fixed 153-class singleton/pair footprint keeps 136 original classes and changes 17; median raw quartic coefficient size is 440 bits (original 381).

For unimodular replacements the intrinsic half-lattice is unchanged. The finite chart footprint and raw coefficient sizes can change; whether that improves actual exceptional-point visibility is **UNMEASURED**. This 289-centre diagnostic uses no quartic reduction or point search and cannot retroactively affect the primary classification.

A post-run generic-data-only score, maximum retained Gram determinant, chooses `arm-14` with determinant 8237. It succeeds in 4.276s, cost rank 2 among successful arms. In this core the first recovered representative is exactly `P14-P16`; the missing P14 was not supplied to the producer. The determinant score can be computed from the retained data alone, and the orthogonal missing height is 1092/det(core), but this is a post-hoc observation, not prospective validation.

The fastest successful frozen core was `arm-11`. Core determinants, selected centre heights and construction cost are available from generic data, but the observed success/cost comparison is post-run. No generic-only chooser has been validated, no selector was optimized on these outcomes and no additional reconstruction arms were run. Transferring a prespecified rule to a genuine MW16 parent is the next experiment, not a result of this calibration.

## Claim boundary

The conclusion is **CORE_DEPENDENT_16_TO_17_RECOVERY**. The fresh RR split-component procedure can manufacture a missing generic quotient representative from some lower-rank cores without receiving the missing section. The successful quotient recoveries are real sections over the original Q(t), not base-change points or specialization-only points. The bounded misses are properties of this finite selector/degree/kernel-row policy; they do not prove geometric inaccessibility or implicate every other determinant1092 construction.

The parent remains the already proved full geometric/arithmetic MW17 surface of determinant1092 and Picard19. Curve302 retains its established rank-at-least31 boundary. This experiment supplies no eighteenth original generic direction, no explanation of the fourteen exceptional directions, no new surface, and no universal MW16→MW17 theorem.

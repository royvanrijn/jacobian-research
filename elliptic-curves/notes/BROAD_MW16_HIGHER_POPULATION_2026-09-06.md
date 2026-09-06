# Broader initial MW16 territory with the selector held fixed

**All 3,058,897,488 initial addresses and 1,310,720 retained scores replay.
The 10,240 distinct scalar models and strict initial-score retention coverage
pass independent reconstruction and replay. All 2,580 point boxes and 182
isolated point-proof stages pass. Three new inventory curves have certified
rank lower bounds 24, 24 and 22. Exact ranks and record claims remain open.**

The user requests widening the earlier search stage because the later
selection appears useful. This campaign changes initial parameter coverage;
it preserves the score expressions, prime bands, retention per slice,
scalar candidate count and final point-search budget. The previously frozen
[narrow higher-band trial](MW16_JOINT_HIGHER_ANNULI_2026-09-06.md) remains a
separate experiment with its original evidence and outcomes.

| Stage | Preserved narrow trial | Broader trial |
|---|---:|---:|
| Untouched denominator slices | 20 | 320 |
| Primitive parameter addresses | 191,215,782 | 3,058,897,488 |
| Retained candidates per signed slice | 4,096 | 4,096 |
| Retained full-score candidates | 81,920 | 1,310,720 |
| Scalar-score candidates | 10,240 | 10,240 |
| New-fibre point finalists | 60 | 60 |
| Maximum initial point boxes | 2,580 | 2,580 |
| Strongest certified lower bound | 19 | 24 |
| Fibres gaining beyond the generic sixteen | 3 | 10 |
| New inventory curves with bound at least 22 | 0 | 3 |

## Certified curves and finite comparison

| Inventory ID | Family | Parameter | Certified rank lower bound |
|---|---|---:|---:|
| `new-20260906-196` | `a1-fibration-02` | -22059/7204 | 24 |
| `new-20260906-197` | `a1-fibration-03` | 24405/86 | 24 |
| `new-20260906-198` | `a1-fibration-05` | -37260/16691 | 22 |

All three lie in 16384 < H <= 65536. The sixty final bounds are fifty 16s,
three 17s, three 18s, one 21, one 22 and two 24s. All sixty curves are
mutually nonisomorphic over Q and unmatched among the 593 pinned catalogue
equations and 1,225 prior equations. The complete
[point results](../../artifacts/generated-results/elliptic-curves/broad60_mw16_results_v1.json)
and [experiment aggregate](../../artifacts/generated-results/elliptic-curves/broad60_mw16_experiment_v1.json)
include exact admission histories, rational maps, generic transports and
full-cloud independence checks modulo 2, 3 and 5. All
[182 isolated replay stages](../../artifacts/generated-results/elliptic-curves/broad60_mw16_point_portable_replay_v1.json)
pass from a 56,982,329-byte standalone point-evidence bundle.

The broad cohort has 36 summed directions beyond the respective sixteen-point
generic subgroups, versus five in the narrow cohort. This sum concerns sixty
different curves; it is not the rank of one curve. With the selector and point
budgets held fixed, this finite outcome supports exploring a broader initial
population. The two disjoint populations are not replicated trials and do
not establish a discovery density or selector optimality. The upper band
produces no bound at least 22 in this cohort.

All three global minimal models and their transported independent points
are [proved](../../artifacts/generated-results/elliptic-curves/broad60_high_rank_minimal_proof_v1.json)
and available in the executable
[Sage export](../../artifacts/generated-results/elliptic-curves/new_broad60_high_rank_curves.sage).
The [V18 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v18.json)
preserves the earlier IDs and adds these three. Its 198 finite independence
certificates, rational-isomorphism exclusions and equation CSV are separately
replayed. Catalogue absence does not prove universal novelty.

A fixed [conductor audit](../../artifacts/generated-results/elliptic-curves/broad60_conductor_bounds_v1.json)
checks all discriminant primes through 10000. Sage Tate data and independent
PARI local reductions agree on all 35 local pairs. The upper bounds have
424, 380 and 439 bits respectively, exceeding the pinned comparison thresholds
of 288, 288 and 247 bits. No conductor record follows. Residual cofactors of
114, 101 and 111 decimal digits remain unfactored; upper bounds above the
benchmarks do not exclude smaller true conductors. The
[independent replay](../../artifacts/generated-results/elliptic-curves/broad60_conductor_bounds_pari_replay_v1.json)
passes without a factorization escalation.

## Fresh territory and fixed ordering

Both experiments use 16384 < H <= 65536 and 65536 < H <= 262144, where
H=max(|n|,d), gcd(n,d)=1 and d>0. Each of the five families receives both
signs in each band. The broader trial takes sixteen new denominator residues
per band/family/sign, modulo 256 and 4096 respectively. A frozen SHA256 seed
orders the same-parity choices. Both the earlier outer-band residue and the
preserved narrow-trial residue are excluded before taking sixteen.

The resulting 320 slices are mutually disjoint within each family and
exclude all previous slices in those bands. They are all outside the compact
region. These are primitive parameter addresses, not a claim of three billion
nonisomorphic or novel curves.

All 3,510 selection primes through 32749 contribute before retaining 4,096
per slice. The combined pool has 131,072 survivors per band/family. As before,
1,024 distinct prospective equations per band/family receive scalar traces
through 65521, and six per band/family become point finalists. Exact
within-prospective-roster rational-isomorphism deduplication, tie rules and
wholly disjoint validation at 65537 through 131071 are unchanged. No public
record equation, parameter, point, rank, j-invariant or jump label enters
selection or execution.

The [model and retention audit](../cas/audit_higher_mw16_score_coverage.py)
checks all 10,240 homogeneous models independently and finds 10,240 distinct
j-invariants. In every band/family, the worst admitted scalar candidate is
strictly ahead of all 32 signed-slice retention boundaries in the first
three ordering fields. Thus the 4,096-per-slice truncation removes no address
ahead of the selected initial-score prefix within this frozen population.
The complete certificate and its read-only replay pass in
[`broad_higher_mw16_score_coverage_v1.json`](../../artifacts/generated-results/elliptic-curves/broad_higher_mw16_score_coverage_v1.json).
This does not certify the globally strongest scores through 65521 or predict
point-search outcomes.

## Finite budgets and evidence

The expanded scanner checks complete signed actual-modulus frames and
first-seven ordering for all 320 slices against the independent cached
reader. The first full slice per band has a 45-second gate and is reused
once. The main scan has four workers, 120 seconds per call and a 7,200-second
cap; raw output and scores are checkpointed per slice. The aggregate stores
file hashes and metadata, avoiding a growing in-memory copy of all 1.31
million survivors. Full replay has a 1,800-second cap. No retry or replacement
slice is included.

Later scalar and point budgets remain those of the preserved trial: four
scalar workers, 10,240 candidates, first-twenty cost gate, independent cached
extension agreement and two direct character sums per curve; then sixty
point attempts using only the sixteen generic sections, 43 fixed generic
parity labels per curve, height 125000 and ten seconds per chart. Every map
must pass before any point search. No adaptive follow-up is included.

Terminal exact histories, rational geometry and point-cloud independence
checks precede any catalogue comparison. Final certification also waits for
the preserved narrow cohort, so the prior-equation snapshot includes its
results. Standalone point-only replay is a separate final stage. Wider
coverage tests the initial-population hypothesis; it does not guarantee a
higher rank or establish that the selector is optimal.

The prepared [population comparison](../cas/compare_higher_mw16_populations.py)
requires both complete scalar-selection and disjoint-validation gates before
writing a result. It compares the ten band/family groups at the initial
scalar cutoff, final selection and withheld prime band, with exact score
means and medians. The scalar scoring, finalist selection and validation
program functions have identical syntax trees in the two implementations;
both trials' 10,240 scalar and sixty validation inputs pass the
comparison's correspondence checks. The complete
[comparison certificate](../../artifacts/generated-results/elliptic-curves/higher_mw16_population_comparison_v1.json)
and read-only replay pass.
This is a comparison of disjoint finite populations, not a nested-population
experiment or a calibrated estimate of discovery yield. It reads no point
outcomes or catalogue labels and cannot change either frozen campaign.

Across the ten equal-size band/family groups, broad minus narrow mean scores
are +3.73625 at the initial scalar stage and +3.40007 among the final six,
with increases in all ten groups. On the wholly withheld prime band, only
three groups increase and the overall difference is -0.21276. These numbers
use score units divided by 10^12; exact rational means are in the certificate.
The stronger selection scores therefore do not show an improvement on the
withheld primes. This is not a rank comparison or a significance test.

A separate [local-reduction audit](MW16_OMITTED_GOOD_PRIMES_2026-09-06.md)
finds good-prime terms omitted by displayed singular reductions at 5 and 13.
This is a defined score limitation, distinct from initial population size.
Its diagnostic changes five finalists on the preserved scalar roster; both
completed experiments retain their frozen scores and exposures. A separate
[corrected campaign](CORRECTED_MW16_HIGHER_POPULATION_2026-09-06.md) applies
the missing local terms before retention on another fresh population.

Sources: `../cas/scan_broad_mw16_higher_annuli.py`,
`../cas/score_broad_mw16_higher.py`, the `finish_broad_mw16_higher` controllers,
and the `broad60_mw16` point/proof sources. Protocols and checkpoints are under
`broad-mw16-higher-annuli-v1`, `broad-mw16-higher-scores-v1` and
`broad60-mw16-pari-v1` in local artifacts.

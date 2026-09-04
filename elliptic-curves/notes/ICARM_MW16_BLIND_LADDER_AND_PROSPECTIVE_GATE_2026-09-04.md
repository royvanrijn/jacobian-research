# MW16 blind ladder and prospective parent-pencil gate

Date: 2026-09-04

Status: exact five-curve calibration complete; bounded prospective Nagao stage
complete on nine coordinate presentations; first prospective half-lattice
attempt complete but wholly timeout-censored.

## Exact presentation count

The nine atlas hit labels are not nine fibrations.  Exact rational base changes,
constant square Weierstrass scalings, target-parameter transport, and integral
specialized-MW16 comparisons give:

| target curve | atlas priorities | exact fibration classes |
|---:|---|---:|
| 398 | 16,875; 63,669 | 1 |
| 400 | 53,042; 62,992 | 1 |
| 401 | 57,487 | 1 |
| 542 | 30,486 | 1 |
| 548 | 31,627; 54,835; 63,647 | 1 |

Thus the correct statistical unit is five target curves.  The nine labels are
still useful as nine search coordinate systems: affine changes of the base do
not preserve a bounded projective-height box.  Results from those coordinate
systems are nested trials and are never counted as independent observations.

The five representative fibrations are mutually separated.  Nine of ten
pairs have frozen modular `j`-map landmark witnesses; the remaining
curve-398/curve-548 pair has different exact generic MW16 half-lattice coset
histograms, which is incompatible with a fibration isomorphism.

## Complement-blind five-curve ladder

The search fixture contains each target equation, its sixteen specialized
generic sections, and the exact generic MW16 Gram.  It contains no public
exceptional point, complement coordinate, target rank lower bound, or jump
label.  The first stage completely enumerates `M/2M` and searches every class
in the exact maximum-depth stratum.

| target | demonstrated atlas jump | deepest initial charts | blind initial gain | blind adaptive gain | best blind gain |
|---:|---:|---:|---:|---:|---:|
| 398 | +14 | 12 | +5 | +9 | +14 |
| 400 | +12 | 4 | +5 | +7 | +12 |
| 401 | +11 | 8 | +10 | not run | +10 |
| 542 | +10 | 10 | +10 | not needed | +10 |
| 548 | +8 | 8 | +8 | not needed | +8 |

Repeated coordinate presentations give the same specialized integral MW16
subgroup and the same initial exact recovered points.  They therefore do not
inflate the table.

At the five-curve level, the initial wave recovers 38 of the 55 demonstrated
quotient directions, completely recovers two curves, and comes within one
direction on three.  The previously completed curve-398 five-bit adaptive
wave and the new curve-400 five-bit adaptive wave raise this to 54 of 55.
The curve-400 run searches all `4(2^5-1)=124` adaptive charts at height
100,000; all 124 complete, and exact group-law plus finite-reduction
certificates give `M16 -> M21 -> M28`.

Only curve 401 remains one direction short of its demonstrated atlas lower
bound.  Its literal complete next wave has `8(2^10-1)=8,184` charts and has
not been run.  This purposive five-curve panel is detector calibration, not a
population success-rate estimate.  Matching an atlas lower bound also does
not prove that a target has no further rational directions.

## Prospective ladder

The enforced pipeline is:

```text
cheap Nagao/local ordering
    -> exact bounded half-lattice jump recovery
    -> complete same-minimal-curve residual 2-Selmer gate
    -> only then any expensive continuation
```

There is no Nagao-to-unrestricted-point-search edge.

The first prospective pass evaluates every primitive projective parameter of
height at most 300 in each of the nine coordinate presentations: 109,592
parameters per presentation.  Three disjoint prime blocks and frozen
within-bucket caps leave 104 finalists across the nine presentations.  Nagao
is used only for ordering; it contributes no rank evidence.  Exact
specialization transports all sixteen generic sections, rechecks the
candidate ledger, and finds 104 distinct `QQ`-isomorphism classes inside this
bounded selection.

The first half-lattice attempt covers all 104 finalists and all 856 exact
maximum-depth charts at rational height 10,000 with a two-second per-chart
budget.  Every specialized generic MW16 remains independent, but every one of
the 856 quartic processes times out.  Consequently:

- exact new quotient directions recovered: 0;
- fail-closed structural rejections: 0;
- completed height boxes: 0;
- timeout-censored chart attempts: 856;
- residual-Selmer calls authorized: 0.

This is a complete execution of the declared attempt, not a complete
height-10,000 search.  It supplies no negative rank evidence and does not
evaluate the Nagao ordering.  The raw specialized short models have maximum
rational coefficient size 5,871--6,302 bits, and the generated quartic
coefficients reach 10,140--11,302 bits.  Larger search budgets are therefore
not the next gate.  The next engineering requirement is an exact
`QQ`-isomorphic model and section transport that materially reduces this
arithmetic size, or a half-lattice backend designed for these large
coefficients.  Any transformed candidate must re-prove all sixteen section
identities and independence before its bounded chart search is comparable to
the controls.

## Replay

Build the complement-blind nine-presentation fixture and replay the exact
initial calibration:

```bash
sage -python elliptic-curves/cas/prepare_icarm_mw16_parent_ladder_inputs.sage --check
sage -python elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage
sage -python elliptic-curves/cas/audit_icarm_mw16_parent_presentations.sage --check
```

Replay the curve-400 adaptive calibration and the held-out five-curve
comparison:

```bash
sage -python elliptic-curves/cas/run_icarm_mw16_curve400_adaptive_calibration.sage
python3 elliptic-curves/cas/verify_icarm_mw16_blind_ladder_calibration.py --check
```

Replay the prospective local ordering and exact specializations:

```bash
python3 elliptic-curves/cas/sieve_icarm_mw16_parent_presentations_nagao.py --check
sage -python elliptic-curves/cas/specialize_icarm_mw16_nagao_finalists.sage --check
```

The 104 half-lattice candidates are independently checkpointable with
`--candidate-start` and `--maximum-candidates`.  The completed campaign used
eight contiguous 13-candidate shards:

```bash
for shard in $(seq 0 7); do
  sage -python elliptic-curves/cas/run_icarm_mw16_nagao_finalist_half_lattice.sage \
    --candidate-start $((13 * shard)) --maximum-candidates 13 \
    --height-bound 10000 --timeout-seconds 2 \
    --output "artifacts/local/elliptic-curves/mw16-finalist-half-lattice/shard-${shard}.json"
done

python3 elliptic-curves/cas/merge_icarm_mw16_nagao_finalist_half_lattice_shards.py --check
```

The raw checkpoints stay under ignored `artifacts/local/`; the compact merged
certificate retains candidate order, exact zero-discovery group
classifications, chart timeout counts, source hashes, and full records for any
future positive or fail-closed result.

## Claim boundary

The five atlas jumps are rank lower-bound material, not exact target ranks.
The blind recovery certificates prove only the independent subgroup they
construct.  The prospective zero is entirely timeout-censored.  No candidate
has passed the residual-Selmer gate, no expensive continuation is authorized,
and rank 32 remains open.

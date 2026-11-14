# Arithmetic realization of geometric MW17 lattices

The foundry is now a classification experiment: **which arithmetic obstruction
prevents a geometrically admissible rank-17 Mordell–Weil lattice from being
realized over Q?** The first experiment joins retained evidence; it does not
start a new lattice enumeration or compile another K3 equation.

Mathematical authority remains [MATH_STATUS.json](../MATH_STATUS.json).
The census claim certifies the evidence accounting, not an independent replay
of the underlying geometric, modular-curve or descent proofs.

## Population and identity

The frozen [catalogue](../artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json)
has 827 `(T,NS)` classes and 1,840 positive rank-17 frame classes. Its 171
rootless frames belong to **66** NS classes. These are the geometric MW17
cohort; the other frames have geometric MW ranks 12 through 16. Every
retained row is included in the experiment.

This is an imported population assembled from several bounded construction
languages. All 96 determinant-band/backend enumeration shards remain open.
It is not every admissible lattice, a uniform sample, or a complete census
below a determinant bound. No population-wide frequency follows from it.

Use the content-derived surface ID, actual integral T, and full NS
discriminant form as identity. A determinant is a display field, not a key.
Frames remain separate child records with their original Gram hashes and
embedding provenance. Count NS classes and frame classes separately: a
surface with many retained frames must not dominate the arithmetic statistic.

## Gates and certificate endpoints

All statements concern geometric Picard rank **19**. The experiment asks for
a rootless elliptic fibration with section and full MW group Z^17 over Q(t).

| Gate | Exact positive endpoint | Exact negative endpoint |
|---|---|---|
| Geometric admissibility | Primitive K3-compatible NS and a rootless primitive U complement | A root in this particular frame; excluding its NS requires an all-frames theorem |
| Galois-compatible marking | A rational non-CM period on the full canonical marked curve, with the literal discriminant action and orientation accounted for | A necessary real/local marking condition fails, or the relevant full curve/necessary quotient has no admissible rational point |
| Rational NS | Effective descent of the K3 and of every class of its full saturated geometric NS as an actual Q-line bundle | A proved obstruction to that effective descent; invariance alone does not supply line bundles |
| Rational elliptic fibration | A primitive nef isotropic Q-line bundle and a degree-one rational section class, giving a map to P1 over Q | A proved failure for that marking/pencil |
| Arithmetic MW17 | Rootlessness, Picard19, trivial torsion, and the full rational saturated rank-17 height lattice | An exact rank or necessary-gate obstruction |

The gates describe cumulative requirements, but evidence may arrive out of
order. An arithmetic obstruction can exclude a row whose geometric rootless
gate is still unknown. The ledger retains that unknown and the later proof;
`first_certified_failure` means the earliest **proved** failure in the displayed
order, not that every earlier gate passed. For the 66-class rootless cohort,
the geometric gate has already passed, so its arithmetic attrition is directly
comparable.

`PASS`, `FAIL`, `UNKNOWN`, and `NOT_REACHED` are distinct. A later gate is
`NOT_REACHED` when an exact necessary gate already excludes the target. A
bounded miss, missing model, finite-field point, formal lift, surviving
congruence, or unfinished rational-point calculation cannot become `FAIL` or
`PASS`. Independent later evidence remains visible: 388 and 852 pass rational
NS even though their geometric MW17 gate fails.

### Why the last gates are often consequences

For a Q-defined elliptic fibration with section, Galois-equivariant
Shioda–Tate gives

\[
\operatorname{rank}MW(\mathbf Q(t))
=\dim(NS(X_{\overline{\mathbf Q}})\otimes\mathbf Q)^{G_{\mathbf Q}}
-2-\dim(R\otimes\mathbf Q)^{G_{\mathbf Q}},
\]

where R is the geometric fibre-root lattice. See the repository's
[Theorem A2](RANK_MUTATION_AND_LIFT_THEOREMS.md)
and the underlying [Schütt–Shioda theory](https://arxiv.org/html/0907.0298v3).

If the full integral NS descends as Q-line bundles and contains a rootless U,
choose its effective isotropic sign and reflect into the nef chamber. The
descended line bundle has a two-dimensional space of sections over Q, hence
a pencil to P1 over Q. Divisibility one supplies a section. The primitive
unimodular U and rootlessness give torsion-free MW of rank 19−2=17. Every
geometric section is determined by its divisor class, so all descend. The
[determinant1020 proof](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md)
states these steps with the divisor-descent hypotheses.

Conversely, in Picard rank19 an arithmetic rank17 fibration forces geometric
root rank zero and trivial Galois action on NS tensor Q. Its sections and U
then supply the rational divisor marking. Thus the central arithmetic
question is realization of the **full marked NS**, with effective descent;
finding a small equation is a separate representation problem. Computing
section coordinates is not an extra arithmetic-existence obstruction.

## Frozen baseline

<!-- BEGIN CLASSIFICATION COUNTS -->

| Cohort | Total | Pass | Excluded | Unknown |
|---|---:|---:|---:|---:|
| Retained NS: full rational marking | 827 | 3 | 224 | 600 |
| Geometric MW17 NS classes: arithmetic MW17 | 66 | 1 | 19 | 46 |
| Geometric MW17 frames: arithmetic MW17 | 171 | 2 | 34 | 135 |

Primary arithmetic causes below are disjoint; all overlapping proofs remain in the JSON.

| First arithmetic obstruction | All retained NS | Geometric MW17 NS | Geometric MW17 frames |
|---|---:|---:|---:|
| REAL_CONTENT | 70 | 11 | 11 |
| REAL_PERIOD | 31 | 0 | 0 |
| LOCAL_SQUARE_TWO | 36 | 5 | 6 |
| SCALED_SPLIT_MARKED_CURVE | 84 | 0 | 0 |
| SPECIFIC_MARKED_CURVE | 3 | 3 | 17 |

<!-- END CLASSIFICATION COUNTS -->

The 19 excluded geometric MW17 NS classes have three primary causes:
content at the real place (11), the required square-two vector obstructed
locally (5), and a global marked-curve rational-point obstruction (3).
These proofs already existed; the new result is their joined accounting.
The remaining 46 are unresolved, so **“almost all fail” is not established**.

The full 827-row admission partition is 224 excluded, 3 admitted, 600 unknown.
The admitted rows are 948, 388 and 852. Only 948 has geometric MW17; both other
NS genera are globally rootful. The 759 further rows without a retained
rootless frame or an all-frames exclusion retain geometric status UNKNOWN.
In particular, a source frame of rank12 does not establish a maximum of12
for its NS.

The obstruction ordering is frozen as content, real-period dichotomy,
local square-two representation, scaled-split marked curve, and other
marked-curve proofs. Changing that order can change the primary-cause
histogram; it cannot change the union of excluded rows. For example, the
500 and 750 full-curve exclusions overlap the earlier real content condition.
All applicable retained proofs appear in `obstructions`.

### Controls outside the frozen population

| Actual NS | Full rational NS | Geometric MW17 | Arithmetic MW17 |
|---|---|---|---|
| [622 Inose](DET622_INOSE_RATIONAL_MARKING_SOURCE_2026-09-15.md) | Proved existentially after a constant twist | [Globally excluded](DET622_GLOBAL_ROOT_OBSTRUCTION_2026-09-15.md) | Excluded |
| [1020](DET1020_RATIONAL_MARKING_SOURCE_2026-09-15.md) | Proved, with a later explicit source | [Rootless frame proved](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md) | Proved; [physical pencil and section classes](DET1020_PHYSICAL_MW17_FIBRATION_2026-09-15.md) certified |

Neither changes the 827-row denominators. The determinant1092 recovered
parent and its distinct class1 fibration remain additional known constructions
in the [K3 map](README.md); no unproved identity with a catalogue row is used.

## Reproduction and inspection

- [Pinned input and claim manifest](data/arithmetic/mw17-classification-inputs-v1.json)
- [All NS rows](../artifacts/generated-results/elkies-k3-mw17-classification-v1/surfaces.tsv)
- [All frame rows](../artifacts/generated-results/elkies-k3-mw17-classification-v1/frames.tsv)
- [Full evidence, controls and summaries](../artifacts/generated-results/elkies-k3-mw17-classification-v1/classification.json)
- [Dependency-free builder/checker](scripts/classify_different_ns_mw17.py)

The TSV exports spell JSON nulls as `NONE`; `UNKNOWN` remains a distinct gate state.

Run from the repository root:

```sh
python3 research/elkies-k3/scripts/classify_different_ns_mw17.py --check
```

After an intentional evidence review, `--write` regenerates the three data
files and the count tables above. It does not refresh the input manifest.
Missing or changed input bytes, changed supporting claim scopes, duplicate
identities, contradictory positive/negative evidence, or incomplete coverage
stop the build. There is no CAS import or automatic reconstruction. The
manifest pins both historical packets and their later successor certificates;
old producer scripts and certificates are preserved.

## Next experiment and interpretation

The next bounded theorem target is the retained
[determinant1236 branch identification](DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md).
Its exact quotient points and conditional CM lifts are already computed.
Identifying the candidate's particular cubic branch orbit with the actual
marked Shimura branch would decide one of the 46 unresolved geometric rows.
Until that identification is certified, it remains UNKNOWN.

For each later case, freeze the actual T/NS identity, gate, theorem hypotheses,
finite input/coverage boundary and certificate endpoint before computation.
Record failure of a proposed model separately from failure of the moduli
problem. A complete determination of a marked curve's rational points is
different from a height-bounded scan. Added lattices enter a new cohort version;
they do not retroactively change this denominator.

The current evidence suggests studying **real involutions, local integral
representation constraints, and rational non-CM points on full marking curves**.
It does not identify a single complete obstruction. These mechanisms, effective
surface/divisor descent, and geometric root avoidance are distinct questions.
Rank-19 periods and CM rank20 specializations are treated in
[Elkies's K3/Shimura construction](https://arxiv.org/abs/0802.1301);
the scaled-split exclusions use the rational cyclic-isogeny classification
recorded in [Banwait–Najman–Padurariu, Theorem1.1](https://arxiv.org/html/2206.08891v3).
This experiment makes no literature-priority or statistical prevalence claim.

<!-- status-consumer: EC-K3-DIFFERENT-NS-CLASSIFICATION-BASELINE-20260915 7de9b16860a4bb0f -->

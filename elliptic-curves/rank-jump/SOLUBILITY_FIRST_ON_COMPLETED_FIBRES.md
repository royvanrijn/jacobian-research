# A rational-solubility condition explains part of four new jumps

The organizing question is now **what makes several pre-existing covers
rationally soluble at this parameter?** Describing recovered points, their
halving fields, or their lattices is secondary unless that description supplies
an obstruction or a construction available before those points.

A completed bounded test supplies a partial positive chain on the original
published R17 family:

\[
q_1(t_0),\ldots,q_k(t_0)\in\mathbf Q^{\times2}
\Longrightarrow P_1(t_0),\ldots,P_k(t_0)\in E_{t_0}(\mathbf Q)
\overset{\text{separate certificate}}{\Longrightarrow}
\operatorname{rank}\langle G,P_1,\ldots,P_k\rangle/G\ge d.
\]

Here \(G\) is the marked rank-17 generic subgroup. On four completed fibres,
the certified \(d\) is 3, 3, 3, and 2. This explains a **subblock**, not the
whole observed jump. The missing implication is a simpler arithmetic event
forcing several square conditions together, followed by an explanation of
the uncovered directions. No common low-degree rational construction for the
whole +8, +10, or +14 quotient has been found.

## Frozen test and results

The [protocol](SOLUBILITY_FIRST_PROTOCOL.json) fixes all **39,119 finite-chart
bisections** in the existing generic lattice atlas, excluding its sole inverted
chart explicitly. It uses all 32 completed initial `08234` fibres of Compact192,
with the proved change of parameter \(t=-(s+50)/26\). This is the same generic
family and integral subgroup as published R17; see the existing
[exact transport](../notes/COMPACT_CROSS_FAMILY_INCIDENCE_2026-09-05.md).

The dictionary does **not** select the 37 historically successful covers.
Its construction uses generic sections and lattice data. All 1,251,808 tests
are exact rational square evaluations. Explicit lift maps are retained based
only on those square values. The later quotient check uses the generic prefix
and constructed points; exceptional witness points serve only to report the
previously measured quotient and compare finite-character spans. No new
parameters, point-search boxes, scoring policies, or workers were introduced.

| Fibre | Compact parameter s | Observed quotient rank in retained subgroup | Soluble covers | Certified constructed quotient rank |
|---|---:|---:|---:|---:|
| 08234-003 | 326/5 | 7 | 4 | at least 3, at most 4 |
| 08234-009 | 774/149 | 8 | 4 | at least 3, at most 4 |
| 08234-020 | -2185/106 | 5 | 3 | exactly 3 |
| 08234-027 | 2627/65 | 5 | 2 | exactly 2 |
| 08234-002 | -20/2827 | 2 | 1 | exactly 1 |
| 08234-026 | -11/235 | 6 | 0 | exactly 0 for this dictionary |
| All 18 observed zero-gain fibres | full cohort in certificate | 0 | 0 | exactly 0 for this dictionary |

The other eight observations also have no split covers. There are 14 distinct
successful labels and 14 total hits; these fibres do not reuse a successful
label. No tested value is a branch degeneracy. The constructed points do not
increase the finite-character rank of the retained witness bases; that fact
alone does not prove containment in their rational spans. In the four-cover
cases the character calculation proves three independent images, **not** a
rational relation among all four.

All full fibre ranks and full quotient ranks remain **UNKNOWN**. A zero-gain
observation means that the completed bounded search certified only the generic
rank, not that the curve has no exceptional points.

## An explicit condition on the +8 fibre

For `08234-009`, set \(t_0=-4112/1937\). Four polynomials from the generic
dictionary, each obtained by removing a rational square factor from its cover
equation, are

\[
\begin{aligned}
f_1(t)&=47184481444+25665872604t+3559010841t^2,\\
f_2(t)&=531230427916489-49864417795362t-1905389166951t^2,\\
f_3(t)&=409689-1439214t+328441t^2,\\
f_4(t)&=181288768144+124961177432t+19675777849t^2.
\end{aligned}
\]

Their positive square roots at \(t_0\), in order, are
\(181066742/1937\), \(48560375165/1937\),
\(4307419/1937\), and \(10196244/149\). The atlas labels are respectively
`orbit-0911e`, `orbit-0a037`, `orbit-1795d`, and `orbit-18f5d`.
Their explicit maps have the form

\[
x=x_0(t)+u x_1(t),\quad y=y_0(t)+u y_1(t),\quad u^2=q(t).
\]

Both coefficients in the identity \(y^2-x^3-Ax-B=0\) are verified as
polynomials in \(t\), before substituting the parameter. Thus simultaneous
square values really are a sufficient global rational-solubility event.
They yield a certified three-dimensional subblock of the observed
eight-dimensional quotient.

The **whole dictionary** predates these exceptional points. This particular
successful quartet was extracted retrospectively and is not a pre-registered
four-condition predictor. Its usefulness is that it gives a concrete system
whose simultaneous solubility can now be investigated. Merely renaming the
conjunction of four square tests a shared mechanism would not explain it.

## Comparisons and their limits

The [new comparison panel](../../artifacts/generated-results/elliptic-curves/rank_jump_completed_cohort_panel_v1.csv)
contains 215 exposure observations of 212 distinct initial curves: 192 R17,
20 MW16, and three later MW16 observations of already included equations.
Every row rechecks curve membership and independent finite Kummer fingerprints
for the ordered generic prefix and witness subgroup. Ranks 17/16 refer to
those generic subgroups; the quotient count is the exact difference of
independent subgroup ranks, not an exact curve rank.

The following three pairs were fixed by the
[panel protocol](COMPLETED_COHORT_PANEL_PROTOCOL.json), using each qualifying
frame's rank26/observed17 pair with closest multiplicative parameter height,
then parameter distance and identifier. They preserve family, point-search
height 125000, and completed initial box count. They are not score matched.

| Frame | High parameter (quotient +9) | Low parameter (observed quotient 0) | Height ratio | Boxes each |
|---|---:|---:|---:|---:|
| 07ca9 | 3307/1128 | 2935/1939 | 3307/2935 | 43 |
| 103b2 | 726/761 | -1049/2296 | 2296/761 | 43 |
| 11952 | -1826/2583 | -3891/1396 | 1297/861 | 49 |

These are three new controlled comparison cases, **not** three solved
solubility comparisons. All six cubics have S3 Galois group, as do all 215
observations. That supplies no discriminating solubility condition. The
published-R17 bisection maps have not been transported to these different
frames; applying their equations there would be invalid. Their simultaneous
cover, Selmer, and CT computations remain open. The 103b2 height match is
particularly weak, about a factor of three.

The direct solubility test instead covers the entire equivalent `08234`
frame. Its high/zero comparison has a substantial parameter-scale imbalance.
For its +7, +8, and +6 fibres, even the closest observed zero by multiplicative
compact-parameter height is `08234-013`, \(s=-1569/1505\), with ratios about
4.81, 2.03, and 6.68 respectively. The +7 parameter is also far away on the
real parameter line. A change to published coordinates changes these heights
again. The cohort was selected by Agent 1's search, not by random sampling.
The absence of splits in its observed-zero controls is therefore a descriptive
result, not a validated rank-prediction advantage.

The newer MW16 data also corrects any blanket statement that adaptive waves
always find nothing: `a1-fibration-01-052`, parameter \(-1867/270\), progresses
from rank-at-least 26 to 27, so its observed quotient goes from +10 to +11.
The panel keeps those exposures separate. This is evidence about recovery,
not by itself about why the additional class is soluble.

## Research priorities after this test

1. **Most concrete mechanism: simultaneous rational splitting of explicit
   generic covers — solubility.** It now constructs certified subblocks on
   four completed fibres without their exceptional points. Seek a common
   arithmetic condition explaining several splits; require independence
   modulo the original generic subgroup as a separate gate.
2. **Block obstruction changing from Sha to rational — solubility.** The
   existing norm-defect and twist controls show why matching cubic fields,
   local supports, or Selmer dimensions cannot settle this. A parameter-level
   trivialization of the shared torsors is still missing on production fibres.
3. **Extra rational sections after a specified base change — incidence.**
   This can provide candidate directions; it becomes a jump mechanism only
   after a rational lift of the specialization and surviving quotient
   independence have been proved. Geometric rank capacity alone does neither.

Weak explanations include a smaller torsion Galois group (uniform S3 in this
panel), a shared halving field inferred from independent points, and cover
counts without quotient accounting. The initial half-lattice chart success
and adaptive search yields remain **visibility** evidence. They must not be
promoted to incidence or solubility certificates.

The next theorem-directed computation should target one of the two
four-cover systems above: determine whether a point-independent shared norm
or descent obstruction explains its simultaneous square values, and whether
that explanation is absent on fixed same-family controls. It must not fit a
new auxiliary curve through the supplied points and then call the resulting
anchor soluble by prediction. A positive result must state the common
condition explicitly and prove that it produces at least two independent
quotient images. A negative result excludes only that proposed construction.

Agent 1 can eventually use an exact pre-point simultaneous-splitting feature
as a **sufficient construction gate**, once its coverage, cost, and independent
validation are understood. The present cohort is insufficient to recommend
changing selection. No active search policy has been changed.

## Reproduction and evidence

The gzip geometry input contains all tested equations and is deterministic;
replay does not need the large untracked atlas or current Agent 1 outputs.
Capture provenance pins the atlas bytes. The compact input records all selected
maps and model coefficients. The independent verifier evaluates the square
conditions using rational fractions, separately from the producer's homogeneous
integer test, and checks the 14 generic polynomial lift identities.

```sh
python3 elliptic-curves/rank-jump/completed_cohort_panel.py check
python3 elliptic-curves/rank-jump/solubility_first.py check
python3 elliptic-curves/rank-jump/verify_solubility_first.py check
```

- [Completed panel certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_completed_cohort_panel_v1.json)
- [Solubility and quotient certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_solubility_first_v1.json)
- [Independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_solubility_first_verification_v1.json)
- [Earlier historic panel, paired cases, and mechanism audit](ANALYSIS.md)
- [Historical split-cover collapse](BRANCH_BLOCKS_AND_SPECIALIZATION.md)
- [Norm defects and production Sha blocks](SCALAR_TWIST_BLOCKS_ARE_ELLIPTIC_NORM_DEFECTS.md)

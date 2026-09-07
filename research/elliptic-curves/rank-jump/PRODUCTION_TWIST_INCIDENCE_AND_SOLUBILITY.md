# Production controls: exact incidence contraction and persistent Sha blocks

Follow-up: [generic classes alone force all three relative contractions](GENERIC_SUBGROUP_FORCES_TWIST_CONTRACTION.md).
The exceptional-point dependency below concerns the individual full boundary
dimensions; the relative R17-high comparison no longer needs it.

The three production fibres with a certified complete relaxed boundary now
have **exact full-Selmer dimension comparisons with their scalar −1 twists**.
Their strict incidence space stays fixed, while part of the boundary disappears
and a retained rational block acquires a nonzero CT pairing.

| Original fibre | Marked generic / independent witness / observed quotient | Full Selmer dimension, original → −1 twist | Twist Sha[2] dimension, at least | Whole twist rank, at most |
|---|---:|---:|---:|---:|
| A1/MW16-05, 307/206 | 16 / 25 / +9 | 25+ε → 17+ε | 8 | 9+ε |
| Published R17, −2300/843 | 17 / 24 / +7 | 24+ε → 18+ε | 6 | 12+ε |
| Published R17, −1561/3133 | 17 / 17 / observed 0 | 17+ε → 14+ε | 6 | 8+ε |

There is a **separate unknown ε for each row**, unchanged within its
original/twist pair. Respectively,
ε = dim Cl(O_K,S)/2 − 10, − 8, or − 6, and ε≥0.
No numerical upper bound on ε is supplied. The observed-zero control is
censored, and none of the original ranks becomes exact from this calculation.
The original family's generic subgroup does not survive this twist as a
rational subgroup, so a twist rank is not a jump count in that family.

The comparison separates two effects that a rank proxy would otherwise mix:

1. **Incidence:** the full Selmer dimension falls by exactly 8, 6 and 3.
2. **Solubility:** inside the strict classes that remain Selmer on both curves,
   the twist has obstruction ranks at least 8, 6 and 6.

No visibility endpoint is measured. In particular, the lost local boundary
directions and the persistent but obstructed strict directions are different
spaces; their counts must not be interpreted as two counts of the same points.

## The complete boundary makes the comparison exact

Use the retained integral short cubic f(x)=x³+Ax+B, K=Q(θ), and the complete
set S containing 2, infinity and every prime dividing its discriminant.
The twist is y²=x³+Ax−B, with labelled root −θ. Its Kummer class at
(x,y) is x+θ. Thus both curves have the same labelled module E[2], the
same S-unramified relaxed group R_S, and the same strict kernel

\[
 U=\ker(\operatorname{loc}_S:R_S\to\Omega_S)
   =\operatorname{Hom}(\operatorname{Cl}(\mathcal O_{K,S}),\mathbf F_2).
\]

The [derivative certificate](DERIVATIVE_RECIPROCITY_AND_COMPLETE_BOUNDARY.md)
already proves that the **entire** relaxed boundary is

\[
 \Lambda=\operatorname{loc}_S(R_S)
 =\operatorname{loc}_S(W)+
   \langle\operatorname{loc}_S(-\operatorname{disc}(f)f'(\theta))\rangle.
\]

Its dimension is ℓ=16,17,12. Here W is the retained rational subgroup,
and its image has dimension ℓ−1. Let L₀ and L₋ be the products of the
original and twist local point images. By the definition of Selmer conditions,

\[
 0\longrightarrow U\longrightarrow\operatorname{Sel}_2(E^{(d)})
 \longrightarrow\Lambda\cap L_d\longrightarrow0,
 \qquad d=1,-1.
\]

Outside S both local conditions are the same unramified condition: the
twist adds no bad prime there. Local self-duality and the usual twist
identification are reviewed in
[Morgan–Paterson, §3.1 and Lemma 4.5](https://arxiv.org/pdf/2011.04374).
Here the actual boundary intersections are computed in the retained cubic
squareclass coordinates, rather than estimated by a parity formula.

| Fibre | dim(Λ∩L₀) | dim(Λ∩L₋) | dim(Λ∩L₀∩L₋) | dim L₀ − dim(L₀∩L₋) |
|---|---:|---:|---:|---:|
| A1 +9 | 15 | 7 | 7 | 8 |
| R17 +7 | 16 | 10 | 10 | 6 |
| R17 observed 0 | 11 | 8 | 8 | 3 |

Consequently Sel₂(E^(−1)) is an actual **subspace** of Sel₂(E), under the
fixed labels, in all three cases. The dimension drop attains the maximum
allowed by the change of local conditions. For any two such local products,
intersection with L₀∩L₋ can remove at most
dim L₀−dim(L₀∩L₋) from Λ∩L₀. Our equality is a property of these
certified global boundaries, not a general assertion about negative twists.

## The quotient directions in the persistent obstruction block

The [scalar cup certificate](INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md)
gives the twist's CT form on the retained strict rational space W∩U:

\[
 M_{-1}=A+A^{\mathsf T},\qquad
 A_{ij}=\chi_{\beta_i}([\mathfrak J_{\beta_j}]),
 \quad (\beta_j)=\mathfrak J_{\beta_j}^{2}\text{ outside }S.
\]

Its ranks are 8,6,6. These are unchanged previously certified values;
the new calculation supplies the **full Selmer dimensions** to which they
can now be applied. Since rational Kummer classes lie in the full CT
radical, rank E^(−1) ≤ dim Sel₂(E^(−1)) − rank M₋. This proves the
whole-twist formulas in the first table, without assuming that the retained
strict space is all of U.

The new [generic-adapted block certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_production_twist_blocks_v1.json)
retains every basis transformation and its exceptional quotient mask.
Writing H for the binary alternating plane with matrix [[0,1],[1,0]],
the decompositions are:

* **A1 +9:** W∩U has dimension 10; its generic strict subspace has
  dimension 1. That generic class pairs nontrivially with an exceptional
  class. Split off their H. The orthogonal complement has dimension 8,
  rank 6 and radical dimension 2, and maps injectively into the exceptional
  quotient. Thus the nine strict quotient directions consist, in this
  chosen basis, of one partner of the generic class, three further H pairs,
  and two retained radical directions.
* **R17 +7:** W∩U has dimension 8; its two-dimensional generic strict
  subspace is already a nonsingular H. Its orthogonal complement has
  dimension 6, rank 4 and radical dimension 2. This accounts for the six
  strict quotient directions as two H pairs and two retained radical
  directions. The seventh observed quotient direction is the additional
  global boundary direction, as in the earlier support audit.
* **R17 observed 0:** all six retained strict classes belong to the generic
  subgroup, and form three H pairs. There is no observed exceptional
  quotient to decompose. This control rules out treating a large scalar-cup
  block by itself as a high-jump signature.

These are decompositions of a solubility obstruction, not canonical
partitions of the original point lists. Changing basis changes the displayed
pairs. A retained radical direction may pair with an uncomputed class or
remain Sha for another reason; it is not a soluble-cover certificate.

Equivalently, with s₀=dim Sha(E)[2], the rank differences satisfy

\[
 \operatorname{rank}E-\operatorname{rank}E^{(-1)}
 \ge 16-s_0,\quad12-s_0,\quad9-s_0
\]

respectively. This follows from the exact Selmer differences and the twist
Sha lower bounds. It is conditional on the displayed unknown s₀, not an
unconditional numerical rank-drop claim. If ε=0 were independently proved,
the original ranks would be exactly 25,24,17, and the twist bounds would be
9,12,8.

## What can be known before exceptional points are supplied

The local intersection cost d=8,6,3 uses only the equations and full local
point images. Those images can be generated with local arithmetic; they do
not require global exceptional points. This is an **incidence comparison
feature** for a specified twist, not a rank predictor. In general it only
bounds a possible Selmer dimension difference. It neither supplies a large
strict class space nor establishes rational solubility.

For A1 +9 and R17 observed 0, the marked generic points together with the
point-blind derivative class already span the full relaxed boundary:
their dimensions are 16/16 and 12/12. Their exact formulas can therefore
be stated without exceptional point input as

\[
\begin{array}{c|cc}
 &\dim\operatorname{Sel}_2(E)&\dim\operatorname{Sel}_2(E^{(-1)})\\
 \text{A1 }307/206&c_S+15&c_S+7\\
 \text{R17 }-1561/3133&c_S+11&c_S+8.
\end{array}
\]

Here c_S=dim Cl(O_K,S)/2 is still unknown. For the R17 +7 fibre, the
generic-plus-derivative image has dimension 16/17. One additional boundary
direction from the exceptional witnesses is used in the full certificate.
Its exact comparison must remain labelled retrospective.

Computing the full strict characters and their cup pairing directly from
the field would also be point-blind. The present production character bases
were extracted from known points, so neither their dimensions nor their
cup ranks may enter a prospective selector as if obtained independently.

## Frozen experiment and verification

The [first protocol](PRODUCTION_MINUS_TWIST_PROTOCOL.json) allowed three
60-second workers and at most 2048 local x candidates per place. It completed
the observed-zero control, but left one missing local direction at A1 prime
3 and R17-high prime 2. Those first-pass records remain incomplete.

The [two-place completion](PRODUCTION_MINUS_TWIST_COMPLETION_PROTOCOL.json)
allowed 30 seconds per target, depth 48 and 16384 candidate evaluations.
It closed the gaps after only 110 and 100 evaluations, at root-tree depths
7 and 11. The new minimal-chart x witnesses are 1127 and 1405. Each retained
cubic value is an exact local square; it need not be a rational square.
These are local image witnesses, not new rational points on the curves.

Independent verification reevaluates every global localization, verifies
each local point-image dimension, checks every nonempty subset of each
finite local basis using PARI `nfislocalpower`, and checks the point witnesses
with PARI p-adic square testing. Sage vector spaces independently reproduce
the boundary intersections and subgroup inclusions. No class group,
norm-equation campaign, global point search, or new parameter sweep ran.
Two verifier harness failures (a dependent input basis and Sage integer
serialization) were corrected before sealing the verification certificate;
their local logs and partial output are retained.

Certificates:
[first pass](../../artifacts/generated-results/elliptic-curves/rank_jump_production_minus_twist_inputs_v1.json),
[local completion](../../artifacts/generated-results/elliptic-curves/rank_jump_production_minus_twist_completion_inputs_v1.json),
[full comparison](../../artifacts/generated-results/elliptic-curves/rank_jump_production_minus_twist_completion_v1.json),
[independent verification and rank bounds](../../artifacts/generated-results/elliptic-curves/rank_jump_production_minus_twist_verification_v1.json).

Replay, from the repository root:

```sh
python3 elliptic-curves/rank-jump/production_minus_twist.py check
python3 elliptic-curves/rank-jump/production_minus_twist_completion.py check
sage -python elliptic-curves/rank-jump/verify_production_minus_twist.py check
python3 elliptic-curves/rank-jump/production_twist_blocks.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_production_twist_blocks.py
```

## Ranked conclusions and next missing implication

1. **Strongest supported mechanism:** a large strict arithmetic incidence
   space can survive a change of curve while its simultaneous rational
   solubility changes in a high-rank CT block. The production comparisons
   now quantify this together with an exact contraction of the full Selmer
   boundary. The [small exact-rank control](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md)
   supplies the closed rational-versus-Sha example.
2. **Useful but insufficient incidence event:** local-condition changes
   can remove several global boundary directions at once. The equation-only
   cost d is interpretable, but the maximal contraction established here
   requires the global boundary certificate. These three selected controls
   do not establish a statistical rank predictor.
3. **Weak explanations:** a shared cubic field, a large local product, or
   an alternating obstruction block alone does not distinguish high jumps.
   The observed-zero control has a nondegenerate six-dimensional block.
   The earlier chart concentration remains evidence about visibility.
4. **Missing computations and theorem:** certify or bound c_S and the full
   CT structure for a production field without exceptional-point input;
   then explain why many surviving classes are simultaneously rational.
   Even unramified norm lifts and Jacobian Selmer lifts do not finish that
   last implication, as the [closed Jacobian control](JACOBIAN_SELMER_LIFTS_CAN_BE_SHA.md)
   demonstrates. No full Selmer radical may be declared rational without
   a further argument.
5. **Information Agent 1 could eventually use:** the known generic local
   boundary and the derivative class can sometimes close the entire boundary
   before a point search. This isolates the missing incidence computation
   to c_S and the missing solubility computation to global classes. Local
   twist cost can be recorded as an incidence comparison; an independently
   computed cup obstruction can exclude rational subspaces. Neither should
   change current candidate scores on this retrospective evidence alone.

The remaining chain is precise: an equation-defined condition must force a
large strict class space **and** a large rational image inside it. These
certificates identify and separate both gates; they do not yet prove a
specialization condition forcing the soluble image.

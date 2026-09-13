# MW18 centre selection in the generic height geometry

The controlled experiment compares nearest-first, deepest, and geometrically
diverse deep classes on the five retained MW18 anchor presentations. They
represent four distinct curves: 531, 534 (two covers), 545, and the historical
rank-28 curve. Each has **ten demonstrated directions remaining beyond the
displayed MW18 subgroup**. Fourteen is the increment needed to reach 32 from
18; it is not an observed jump in these families.

**Completed result:** all 600 declared charts finish without timeout. Deepest
selection certifies 22 additional-direction recoveries across the five
presentations; diverse deep selection certifies 21; nearest-first certifies
none. The [retained comparison](../../artifacts/generated-results/elliptic-curves/mw18_deep_centre_comparison_v1.json)
and [witness bundle](../../artifacts/generated-results/elliptic-curves/mw18_deep_centre_comparison_v1.zip)
preserve the frozen protocol, exact geometry, full parity witnesses, chart
transcripts, and independent point certificates.
The [fresh-cache replay](../../artifacts/generated-results/elliptic-curves/mw18_deep_centre_replay_v1.json)
passes for the compact bundle, including every geometry, centre-selection,
chart, and subgroup check.

| Anchor presentation | Nearest-first | Deepest | Diverse deep | Demonstrated remaining |
|---|---:|---:|---:|---:|
| 531 / `08234-12f61` | 0 | 6 | 6 | 10 |
| 534 / `08234-13d7a` | 0 | 2 | 2 | 10 |
| 534 / `08234-1a371` | 0 | 3 | 3 | 10 |
| 545 / `07ca9-08c1e` | 0 | 4 | 3 | 10 |
| historical / `15a68` | 0 | 7 | 7 | 10 |
| Presentation-weighted total | 0 | 22 | 21 | 50 |

The two curve-534 rows are correlated presentations of one curve. These totals
are not fifty independent curve observations. Deterministic deepest selection
has the highest recovery in this comparison, but **the prospective gate fails**:
22 is below 35 and the weakest presentation recovers two, below five. No
winning policy is promoted and none of the 27 frozen prospective candidates
is searched. The result establishes a strong bounded contrast with shallow
selection, while leaving further recovery unresolved.

An additional finite-reduction audit of the unions of the retained independent
policy bases gives joint certified gains `6,2,3,4,7`, again totaling 22. Combining
the two deep arms therefore produces no further certified direction in this
diagnostic. Equal finite ranks alone do not establish equality of their
rational spans or integral subgroups.

Mathematical rank status remains in [MATH_STATUS.json](../../MATH_STATUS.json).
This experiment changes search selection and tests bounded recovery. It proves
neither an exact generic rank nor a specialized rank upper bound.

The subsequent [visibility diagnostics](RANK_JUMP_DIAGNOSTICS_2026-09-05.md)
verify exact recovered/public rational-span relations, locate translated
missing representatives, and find no omission inside completed coverage.
They do not rerun larger boxes or change this frozen gate.

## Exact generic Gram

The [geometry certificate](../../artifacts/generated-results/elliptic-curves/mw18_generic_height_geometry_v1.json)
checks all nine covers, in precisely the ordered basis used by the
specialization certificate: the seventeen inherited sections followed by the
cover section `T`.

Let `G` be the inherited Gram over the old base and `Q = T + sigma(T) = sum q_i P_i`.
Galois invariance and degree-two base change give

\[
2\langle T,P_i\rangle_L=\langle Q,P_i\rangle_L=2(Gq)_i.
\]

For every cover, exact polynomial checks show:

- the old discriminant has degree 24 and is squarefree;
- the smooth quadratic branch polynomial is coprime to that discriminant;
- infinity is smooth and unramified;
- `x=x0+x1*u`, `y=y0+y1*u`, `u^2=q(t)` obey the curve equation;
- the degrees of `x0,x1,y0,y1` are `4,3,6,5`;
- the rational chord trace equals the stated inherited basis word exactly,
  and that word has old height 10.

Thus the pullback has 48 `I1` fibres and `chi=4`. The section is disjoint from
the zero section, including at infinity, and there are no reducible-fibre
corrections. The intersection formula gives `h_L(T)=2*4=8`. Consequently

\[
G_{18}=\begin{pmatrix}2G&Gq\\q^TG&8\end{pmatrix},\qquad
8-\tfrac12q^TGq=3,\qquad
\det G_{18}=2^{17}\cdot948\cdot3=372768768.
\]

These are deductions from the standard height properties, described in
[Schütt–Shioda, *Elliptic Surfaces*, §§11.8, 11.19–11.20](https://arxiv.org/abs/0907.0298).
The certificate establishes the Gram of the displayed generic subgroup.
It makes no claim that this rank-18 subgroup is saturated in the full generic
Mordell–Weil group.

The [builder](../cas/build_mw18_height_geometry.sage) evaluates trace relations
over the function field. It adds signed basis sections in increasing
intermediate height; this keeps the intermediate heights at most ten and
avoids enormous rational functions which would cancel later. No specialized
numerical height computation is needed.

## What is held fixed

All three policies use the same unimodularly reduced generic Gram, the same
retained closest-representative table, `PointedQuarticSearch`, raw curve model,
`metric:16` coordinates, height **100,000**, **40 charts per presentation**, a
20-second worker cap per chart, and finite-independence primes through 1,000.
There are 600 declared charts. At most three independent cells run concurrently;
each has a 30-minute wall limit and a 2 GiB memory limit. Identical chart boxes
can reuse exactly checked search records.
Centres remain fixed while discoveries accumulate; this comparison includes
no adaptive quotient wave and does not establish general backend insensitivity.

The policies are:

1. **Nearest-first:** the existing `VoronoiIterator` selects masks with its
   current diversity window of two, using the correct generic Gram.
2. **Deepest:** the first forty maximum-depth classes, ordered by mask in the
   frozen reduced basis.
3. **Diverse deep:** forty classes from that same complete maximum-depth
   stratum, selected greedily by greatest minimum distance to earlier classes.
   Distances are flat-torus distances: the squared distance between classes
   `a,b` is one quarter of the minimum norm in parity `a xor b`. This removes
   the arbitrary signs of shortest representatives from the diversity score.

All selected representatives are taken from one common table, so changing the
CVP tie solver is not a second experimental variable. A fixed deterministic
tie rule resolves equal diversity scores. The initial proposal of 43 charts
was reduced to 40 **before any point searches**, because one lattice has only
41 deepest classes. The gate and all centre lists were frozen before discovery.

| Anchor presentation | Complete deepest stratum | Squared half-lattice depth |
|---|---:|---:|
| 531 / `08234-12f61` | 47 | 6 |
| 534 / `08234-1a371` | 41 | 6 |
| 534 / `08234-13d7a` | 47 | 6 |
| historical / `15a68` | 51 | 6 |
| 545 / `07ca9-08c1e` | 67 | 6 |

The nearest-first batches remain at depth 2. Lazy nearest-representative
enumeration is therefore not a substitute for maximum-depth selection.
The earlier identity-Gram, twelve-centre, height-10,000 trial remains a valid
bounded experiment, but does not answer this controlled transfer question.

## Exactness of the finite selection

The [preparer](../cas/prepare_mw18_deep_centres.sage) retains one representative
of every one of the `2^18 = 262144` parity classes. Every parity and integral
norm is checked. These norms initially provide **upper bounds** on the true
coset minima; the floating CVP engine is not trusted to prove minimality.

The [exact auditor](../cas/research_runtime/deep_centres.py) uses rational LDL
costs, integer square roots, and exhaustive parity-constrained branching to
exclude all shorter representatives for every proposed maximum. Every other
class already has a representative strictly below that maximum. Together
these checks prove both the maximum and the **complete** deepest stratum.
All pairwise coset minima queried by the diversity policy, and all selected
nearest-first minima, receive the same exact audit. Other table entries remain
labelled upper bounds. A node-budget failure stops certification.

## Recovery and prospective gate

The frozen gate requires at least **35 of 50** additional-direction recoveries,
at least **five on every presentation**, and completion of every declared box.
Policies are ranked by total certified gain, minimum gain on an anchor, then
policy name. Timing does not select the winner. Two presentations on curve 534
are reported together as well as individually; they are not independent curve
observations.

Each result preserves exact chart maps, square hits, curve points, and
incremental finite-reduction independence witnesses for the enlarged subgroup.
Search and policy selection do not read the public complements. A reported gain
means certified independent directions relative to the eighteen input sections.

The candidate roster was also frozen before discovery: three equal-count
parameter-height strata in each of nine covers, with a seeded hash selecting
one candidate in each stratum. This gives 27 candidates drawn from the 178-row
certificate. Both points over every known anchor `t` are excluded, including
parameter infinity. Only a completed, independently replayed anchor gate can
freeze the winning policy and start this roster. A failed gate leaves that
trial unrun; a bounded miss is not an arithmetic exclusion.

## Reproduction

The frozen local work is under `artifacts/local/mw18-deep-centres/`.
The [experiment runner](../cas/run_mw18_centre_experiment.py) freezes inputs,
searches cells, replays them, and enforces the prospective gate. The
[supervisor](../cas/supervise_mw18_centre_cells.py) runs at most three independent
cells, retaining per-chart checkpoints and strict worker logs.

```sh
sage -python elliptic-curves/cas/build_mw18_height_geometry.sage \
  --output artifacts/local/mw18-deep-centres/geometry.json
sage -python elliptic-curves/cas/prepare_mw18_deep_centres.sage \
  --geometry artifacts/local/mw18-deep-centres/geometry.json \
  --directory artifacts/local/mw18-deep-centres/centres \
  --cover 07ca9-orbit-08c1e --cover 08234-orbit-12f61 \
  --cover 08234-orbit-1a371 --cover 08234-orbit-13d7a --cover orbit-15a68
.venv/bin/python elliptic-curves/cas/run_mw18_centre_experiment.py \
  --freeze artifacts/local/mw18-deep-centres/protocol.json \
  --geometry artifacts/local/mw18-deep-centres/geometry.json \
  --centres artifacts/local/mw18-deep-centres/centres
.venv/bin/python elliptic-curves/cas/supervise_mw18_centre_cells.py \
  --protocol artifacts/local/mw18-deep-centres/protocol.json \
  --directory artifacts/local/mw18-deep-centres/anchor-trial --workers 3
```

`--verify` on the preparer checks retained witnesses without floating CVP or
nearest enumeration. `--verify` on the supervisor checks retained charts and
independence without repeating point enumeration. The portable bundle tool is
[`retain_mw18_centre_experiment.py`](../cas/retain_mw18_centre_experiment.py).

The complete portable replay starts from an empty arithmetic cache, rebuilds
the generic height proof from pinned canonical inputs, audits the retained
deep strata and torus distances, and rechecks every chart and subgroup
transition without point enumeration or floating CVP:

```sh
.venv/bin/python elliptic-curves/cas/retain_mw18_centre_experiment.py \
  --verify-bundle artifacts/generated-results/elliptic-curves/mw18_deep_centre_comparison_v1.zip \
  --output artifacts/local/mw18-deep-centres-fresh-replay
```

A fresh `/dev/shm/` output directory also works on Linux and avoids journal
flushes while materializing the disposable arithmetic cache. The successful
audit used that option; its logs are retained under
`artifacts/local/mw18-deep-centres/final-portable-audit/`.

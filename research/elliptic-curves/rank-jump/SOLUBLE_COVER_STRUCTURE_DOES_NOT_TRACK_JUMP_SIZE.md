# The matched panel does not support compressed native defects as the jump mechanism

The requested discrimination test changes the conclusion. Soluble native
covers detect some retained exceptional directions, but neither their
number nor a compressed collision-defect space explains the largest gains
in the tested panel. **All 17 nontrivial compatible blocks have the maximum
realizable defect span n−1.** None has a smaller span than its number of
covers allows.

The [joined comparison panel](../../artifacts/generated-results/elliptic-curves/rank_jump_fibre_discrimination_metrics_v1.csv)
contains every requested metric, with unavailable entries explicitly marked.
It keeps 542 exposure observations, tests 336 of them at 326 distinct
family/parameter addresses, and retains 206 untested observations for which
no valid family dictionary or parameter transport was available. Addresses
are not asserted to be pairwise nonisomorphic curves. The computations use
existing parameters and generic cover equations; no elliptic point search,
new parameter selection, or active-search modification occurs.

The inherited rank labels describe certified **retained subgroups** relative
to the specified generic subgroup. A reported +0 is a censored observation,
not a proved rank-16/17 curve. Full ranks remain UNKNOWN throughout this
panel. The existing point certificates are inherited; 239 stored finite-rank
matrices are independently replayed here, including their rank-17 prefixes.

## Fixed dictionaries and coverage

- Published R17: all 39,120 native bisections. This restores the one previously
  omitted construction chart. Its equation was already transported back to
  the published coordinate by the original constructor; it is evaluated in
  that coordinate. The five historic full-atlas counts agree exactly with
  the existing specialization certificate.
- The equivalent compact 08234 frame uses the proved affine transport
  t=−(s+50)/26 and the same dictionary.
- Native 11952: the existing 1,024 cheapest constructed bisections. This is a
  partial dictionary, not a census of its native geometry. The exact compact
  parameter transport is frozen, and every coefficient of the A and B model
  identities is independently checked. Counts are never pooled with R17's
  much larger dictionary.
- Other R17 frames and A1/MW16: untested where no corresponding dictionary is
  available in this audit. In particular the historic +12…+14 examples are
  retained as coverage gaps, not filled using unrelated R17 equations.

The parameter-only pass performs **6,619,664 exact rational-square tests**.
It receives no exceptional coordinates or gain labels. Here “compatible”
means a nonzero square value of a stored native cover equation over Q, not
merely local compatibility. No branch degeneracy occurs in the tested rows.
Rank labels are joined only after the square census.

## High gains do not consistently mean many covers

All parameters in this table are in published R17 coordinates. The gains
are relative to the same marked rank-17 subgroup.

| Parameter | Retained gain | Soluble covers n | Exact realizable defect span | Collision primes | Full native carrier (genus, degree over t) |
|---|---:|---:|---:|---:|---|
| 3/8, historic rank≥21 | +4 | 25 | 24 | 840 | (385,875,969; 33,554,432) |
| −70/61 | +5 | 13 | 12 | 232 | (45,057; 8,192) |
| 28/117 | +6 | 11 | 10 | 179 | (9,217; 2,048) |
| −288/65, 08234-003 | +7 | 4 | 3 | 18 | (17; 16) |
| −4112/1937, 08234-009 | +8 | 4 | 3 | 23 | (17; 16) |
| −2/377, historic rank≥25 | +8 | 6 | 5 | 62 | (129; 64) |
| −308/251, historic rank≥26 | +9 | 3 | 2 | 17 | (5; 8) |
| 2456/135, historic rank≥27 | +10 | 2 | 1 | 4 | (1; 4) |
| −9529/5471, historic rank≥28 | +11 | 1 | 0 | 0 | (0; 2) |

The large genera are exact consequences of the verified disjoint branch
sets, not the result of constructing or searching enormous curves. For n
independent quadratics with disjoint branch pairs,

\[
\deg(C_n\to\mathbf P^1)=2^n,\qquad
 g(C_n)=1+2^{n-1}(n-2).
\]

These are carriers of the selected native lifts, not necessarily minimal
carriers of the full exceptional quotient. Their n native sections can
specialize to fewer than n independent quotient directions. In particular,
25 split covers on the +4 control are compatible with only four known
quotient directions. Genus and degree here mostly restate n; they are not
independent statistical features.

Within the complete-atlas panel, the three distinct addresses whose best
retained gains are at least nine have only one, two and three split covers.
Among 96 addresses with best retained gain zero, 94 have no split cover and
two have one. This suggests some association with the presence of retained
exceptional points, but does not establish a monotone relationship with jump
size or true-rank specificity. Censoring and the original score selection
remain material.

## Every nontrivial block has full local defect span

The arithmetic pass includes **all 17** compatible blocks with at least two
covers, regardless of rank label. It factors all **519 pair resultants**
within the declared five seconds per block, retaining a checkpoint for each
completed pair. Every factor is checked for primality with proof enabled.
All primitive normalizations remove a rational square, so the normalized
radicands preserve the native squareclasses.

For each block, exact independent local product-curve witnesses certify

\[
\operatorname{span}_{\mathbf F_2}
 \{(v_p(f_1(t)),\ldots,v_p(f_n(t)))\bmod2:
     (t,y)\in H_n(\mathbf Q_p)\}
 =\{\epsilon:\sum_i\epsilon_i=0\}.
\]

The span is taken across places, with the native coordinates retained.
Its dimension is n−1. This is a span statement, not an assertion that every
mask occurs at one fixed prime. Nor does it prove realization by H_n(Q).
The verifier records **82 independent local witnesses**, one per basis
vector across the 17 blocks, and checks their odd valuations and square
product units exactly. These are local residues, not new global points.

The first small-prime audit alone provided lower bounds on realizable span;
its potential collision masks were only upper possibilities. The final
verification uses the completely factored support and attains the universal
upper bound in every block. Thus the table reports actual certified local
span, not just a graph-theoretic potential span.

There is no observed event
“many covers → smaller-than-allowed defect rank” in this census. A global
image or rational relation could still synchronize local classes, but that
would be additional structure not supplied by these collision spans.

## Matched cases, including newer best gains and score false positives

| Comparison | High retained gain | Soluble covers, high vs observed-zero control | Matching and limitation |
|---|---:|---|---|
| R17 −2300/843 vs 1576/2331 | +7 | 3 vs 0 | Same recorded cohort and 43 charts; parameter-height ratio 2331/2300. Full-box completion is not inferred from chart count. |
| R17 −2300/843 vs −1929/3242 | +7 | 3 vs 0 | Same scored cohort and chart count; score gap about 0.0612%; not height matched. |
| R17 −1264/1047 vs 444/2623 | +6 | 0 vs 0 | Same scored cohort and chart count; score gap about 0.0133%; the cover test misses the high member too. |
| Compact11952 4286/1881 vs −8897/3736 | +10 | 0 vs 0 | Same completed full11952 cohort, 49 boxes each; height ratio about 2.08. Partial atlas only. |
| Compact11952 110314/102227 vs 98624/110077 | +10 | 0 vs 0 | Same completed late11952 cohort, 49 boxes each; height ratio about 1.0022. Partial atlas only. |
| Compact11952 2618/26913 vs 9201/29507 | +9 | 0 vs 0 | Same completed late11952 cohort, 49 boxes each; height ratio about 1.0964. Partial atlas only. |

The full11952 observation at 89074/31895 recovers +11, has zero dictionary
hits, and is catalogue curve 12, already known at rank≥29. It is not a new
rank-28 curve. Its same-cohort observed-zero match 89032/35897 has almost
identical parameter height and also zero hits. The portable verification
keeps this catalogue context separate from the recorded subgroup bound.

Across all 161 tested 11952 addresses, the partial dictionary has exactly
one hit, on a +4 observation. All five addresses with retained gain at
least nine have zero hits. This rejects usefulness of **this 1,024-cover
subset** for these cases, not the existence of other native covers or a
soluble block on 11952.

The verification retains eight score-matched R17 comparisons, not just the
illustrative two. Four of those eight comparisons have zero covers on both
members. No score-matched comparison is labelled parameter-height matched.
Candidate pairs in the initial table whose cohort or box metadata is null
are explicitly excluded from the accepted matched-case list. Equal missing
values do not prove comparable exposure.

## Which requested measurements are informative?

| Requested measurement | Outcome in this panel | Layer and interpretation |
|---|---|---|
| Compatible native covers | Exact for the declared dictionaries | Solubility; useful for some subblocks, not a consistent jump-size ordering |
| Pair-mask/collision-defect span | Exact n−1 for all 17 nontrivial blocks | Solubility obstruction structure; no observed compression |
| Collision-prime support | Exact for all nontrivial blocks; complete factor lists retained | Solubility; depends strongly on the selected cover set and its size |
| Rational low-degree relation components | UNKNOWN: no common fixed relation dictionary yet | Potential solubility mechanism; not measured by counting split covers |
| Three-cover blocks extending to a fourth | Exact combinatorial counts at each tested t | Solubility incidences; determined by n for the compatible subset |
| Remaining full simultaneous carrier genus/degree | Exact for every tested compatible set | Generic geometry of these lifts; not full-jump incidence or minimality |

For clarity, there are binomial(n,3) split triples and binomial(n,4) split
quartets. If n≥4 every such triple extends, with n−3 possible fourth labels.
Those counts are reported in the CSV, but do not provide another predictor
beyond n. Likewise the specialized squareclass defect of a subset selected
because all its radicands are squares is automatically zero. The nontrivial
span above belongs to the fixed carrier across its local points; it is not
that tautological specialized zero.

A rational point at a fixed t is also not a rational component of a relation
scheme over the parameter line. The earlier two successful relation examples
were selected retrospectively. Their component counts cannot be copied into
unrelated fibres or counted as a common pre-point feature.

## Revised mechanism ranking and next experiment

1. **Global rational relation or auxiliary-point synchronization remains the
   strongest open mechanism.** It must predict that a rational product or
   relation point reaches the zero affine class. Neither the present complete
   local span nor a large Jacobian/Selmer group supplies this implication.
2. **Native cover splitting explains subblocks**, sometimes substantial ones,
   but the full published atlas misses most known directions on the +10/+11
   controls. Its number of split covers does not explain the ranking of gains.
3. **A small local collision-defect span is unsupported here.** All tested
   nontrivial blocks attain n−1, including moderate-gain examples with many
   more covers than the largest published-R17 controls. Carrier genus and
   simple triple/fourth counts add no independent evidence beyond n.

The next discriminating experiment should fix a **family-wide finite
relation dictionary using generic lattice data**, then count rational
components and surviving fourth lifts across the same high/low panel.
Selection must precede access to which components contain the high points.
The degree bound and trace-word roster must be explicit; an unbounded count
of “all low-degree relations” is not a defined feature. Coverage expansion
for 11952 and a valid A1/MW16 dictionary are also needed before this negative
result can be extended to the new MW16 +10/+11 and historic +12…+14 fibres.
No candidate selector, worker limit, point budget, or rank status is changed.

## Replay

The inputs, generic equations, two arithmetic passes, source hashes and
joined metrics are immutable. Protocols are
[FIBRE_DISCRIMINATION_PROTOCOL.json](FIBRE_DISCRIMINATION_PROTOCOL.json),
[FIBRE_COLLISION_PANEL_PROTOCOL.json](FIBRE_COLLISION_PANEL_PROTOCOL.json), and
[FIBRE_DISCRIMINATION_VERIFICATION_PROTOCOL.json](FIBRE_DISCRIMINATION_VERIFICATION_PROTOCOL.json).

```sh
sage -python elliptic-curves/rank-jump/verify_fibre_discrimination.py check
python3 elliptic-curves/rank-jump/report_fibre_discrimination.py check
```

The first command replays the complete square census, exact parameter/model
transport, independent Sylvester determinants, prime products/primality,
local valuation witnesses, and stored rank matrices. It does not rerun the
underlying elliptic point searches or replace their point-membership proofs.
The second checks the joined CSV and comparison summary.

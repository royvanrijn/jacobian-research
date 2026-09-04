# Residual-Selmer fingerprints and the monotone search sieve

## Result and boundary

The complete groups

\[
  \operatorname{Sel}_2(E_t/\mathbf Q)/
  \operatorname{im}(MW17/2MW17)
\]

for ICARM curves 356 and 385 remain `UNKNOWN`.  The exact cubic BNF and full
local-solubility calculation have not completed, so neither a Selmer upper
bound nor an exact-rank-29 claim is made.

What is now complete is the richer control layer needed around that descent:

- the two rank-29 fibres are exact inputs to the checkpointed complete
  2-descent;
- the twelve certified point classes modulo generic MW17 are embedded as the
  displayed residual lower-bound space;
- the two rigid-character directions are recorded inside that twelve-space;
- every exact control retains its known residual local Kummer subspaces,
  cumulative intersections, delete-one-place ranks, component data, available
  Hilbert information, and sampled `t mod p^k` strata;
- the complete descent worker now emits the analogous delete-one-place ranks
  on the *full* global norm-square ambient space when the BNF is available;
- an incomplete monotone residual sieve may authorize only an explicitly
  bounded point search.  It cannot authorize a Selmer or rank claim.

The machine-readable record is
[`../artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json`](../artifacts/generated-results/elkies-k3-r17-residual-selmer-fingerprints-v1.json).

The exact place-block support code sharpens the otherwise flat delete-one
table. On controls `351,356,376,377,385,12`, the minimum numbers of places
supporting a nonzero known residual class are `4,5,7,5,6,4`. Consequently
every deletion of at most four places for 356 and at most five for 385 still
leaves the known residual localization injective. Prime attribution for the
record blocks is therefore intrinsically a multi-place support problem; a
zero delete-one rank drop does not mean that the place is arithmetically
irrelevant. After quotienting by the two rigid directions, the selected
ten-dimensional blocks have minimum support weights `6` for 356 and `7` for
385, so the non-rigid phenomenon becomes more, not less, distributed across
places. Their block-localization dual matroids are both connected of dimension
ten. Thus neither known ten-block admits a direct-sum partition separated by
audited places. This is an exact support-theoretic indecomposability result on
the known point quotient, not yet a Cassels-pairing or complete-Selmer
indecomposability theorem.

## Complete-descent execution status

The record-pair jobs were run separately. The open-source backends did not
return a Selmer dimension:

- eclib/mwrank on curve 356 reached a strict 3,600-second limit without
  output;
- eclib/mwrank on curve 385 reached a strict 180-second retained run after
  the 2-adic initialization, without completing the descent;
- factor-hinted PARI 2.19 BNF runs used eight relation threads,
  `c1=0.3,c2=4,nrpid=20,max_fact=10000,idex=4`, and strict 1,200-second
  limits. Curve 356 formed a 3,026-ideal factor base and stopped with 2,439
  unresolved ideals; curve 385 formed a 3,092-ideal factor base and stopped
  with 2,511 unresolved ideals.

Both PARI runs retained logs and telemetry but deleted their uncertified BNF
checkpoints. Backend timeout and relation deficit are not mathematical upper
bounds. This is why the complete residual dimensions, the full delete-one
matrix, and the exact-rank question stay `UNKNOWN`.

## Why the full BNF is the wrong front door

The known points themselves now give an exact localization of the class-group
wall.  The comparison covers curves 351, 356, 376, 377, and 385 and native
alternate-Q80 curve 12.  Form the matrix of prime-ideal valuation parities of
all displayed Kummer classes `4*x(P)-zeta` at every prime above the bad
rational primes.  The result is:

| curve | gain over MW17 | full valuation rank | everywhere-even kernel | proved `dim Cl(K)[2]` lower bound | adjusted residual image lower bound |
|---:|---:|---:|---:|---:|---:|
| 351 | 8 | 5 | 20 | 18 | 6 |
| 356 | 12 | 7 | 22 | 21 | 11 |
| 376 | 5 | 3 | 19 | 17 | 3 |
| 377 | 6 | 7 | 16 | 15 | 5 |
| 385 | 12 | 12 | 17 | 15 | 10 |
| 12 (alternate Q80) | 12 | 6 | 23 | 22 | 11 |

These are unconditional lower bounds for the full ideal class groups of the
completed-square cubic fields, not S-class groups and not Selmer upper bounds.
The proof uses only the exact point classes.  Outside 2 and the discriminant a
Kummer class is unramified.  A combination in the displayed valuation kernel
therefore maps to `Cl(K)[2]`.  Because its norm is a positive square, the
kernel of this class-group map is contained in the norm-positive unit
squareclasses, of dimension `r1+r2-1` (the missing unit direction is `-1`,
whose norm is `-1` in degree three).

More significantly, the exceptional rows add **zero** valuation rank modulo
the generic MW17 rows on all six fibres.  Thus every exceptional direction
can be adjusted by a certified generic Kummer class to become everywhere
even.  After the norm-positive unit ambiguity, the three `+12` residual blocks
force at least 11, 10, and 11 cubic 2-class directions, while the `+5` control
forces only three.  Across the observed jump strata the bounds are strictly
ordered as `+5 -> 3`, `+6 -> 5`, `+8 -> 6`, and `+12 -> 10..11`.  The exact
quotient statement and complete data set are canonical in
[`R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md`](R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md).
Because the exceptional points are inputs, this is currently an explanation
of the computational class-group wall, not an explanation of why the rational
points occur or an out-of-sample rank predictor.  The zero incremental
valuation rank makes the observed lower-bound separation largely formal.  In
the prospective direction, a large residual Selmer space still combines
Mordell--Weil and `Sha[2]` contributions; the missing constructive step is to
put rational points on enough compatible explicit coverings.

The replay additionally constructs for every point the exact integral
half-ideal

```text
A_P = (d^2*(4*x(P)-zeta), d^3*4*(2*y(P)+a1*x(P)+a3)),
d^2 = denominator(4*x(P)),
```

and verifies that `A_P^2/(d^2*(4*x(P)-zeta))` is supported only above the bad
primes.  These are ready-made generators for the quotient that the descent
actually needs.  The original five-fibre certificate is
[`../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-v1.json`](../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-v1.json).
The six-fibre comparative certificate is
[`../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json`](../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json).
<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->

The resulting algorithmic target is not `Cl(K)`: it is

```text
Cl(K) / (2*Cl(K) + <bad-prime ideals> + <known Kummer half-ideals>).
```

`r17_kummer_quotient_search.py` makes the companion-selection policy and the
two quotient projections deterministic and independently testable.  The Sage
driver `run_r17_kummer_quotient_sclass_collector.sage` now cycles through
`alpha*I_i`, `alpha*I_i*I_j`, and short `alpha*prod(I_k^e_k)` candidate
lattices, with a declared bias toward the exceptional point block.  It keeps
simultaneous presentations modulo generic MW17 (retaining exceptional
half-ideals as formal exact class columns) and modulo all known points.  The
active target rotates through nonzero products of unresolved columns of the
chosen objective, and each closed relation records actual rank gain in both
projections.  An optional `idealredmodpower(...,2)` engine preserves the
mod-two ideal class by construction while asking PARI for a smaller
representative; on curve 356 it produced no closed row in a 3,000-trial,
factor-base-5000 pilot, so it remains a comparison lane rather than evidence
of progress.

`run_r17_kummer_quotient_sclass_suite.py` applies identical bounded settings
to curves 351, 356, 376, 377, and 385 in both objectives.  Its relation-yield,
row-weight, and descriptive correlation output is the controlled experiment
for the proposed middle layer

```text
t -> 2-division cubic -> residual S-class structure
  -> residual Selmer -> rank jump.
```

The factor-base-240 smoke run completed 1,500 requested reductions per curve
and objective.  It produced seventeen smooth rows across the five
generic-objective lanes and none in the full-known lanes, but every row was
dependent, hence there was no quotient-rank gain.  The raw starting dimensions
are dominated by different factor-base widths and, in the generic
presentation, the deliberately adjoined exceptional columns; they are not
evidence for a correlation.  A separate curve-356 diagnostic at
factor-base bound 5000 found that allowing exponent two produced 33 smooth
relations in 10,000 rotating-target trials, but every row was dependent in
both quotients; odd-exponent-only trials found no smooth row in their stated
budgets.  This directly confirms that norm smoothness is not the optimization
target.  The earlier 100,000-trial single-target miss is a further design
constraint in favour of rotating or batched targets and a larger-prime sieve.

A common factor-base-5000 comparison with exponents one and two and target
products of widths one through three gave the following sharper bounded
fingerprint.  Counts are normalized because two generic and one full-known
lane reached their declared ten-second caps before 3,000 trials:

| curve | displayed MW gain | generic smooth rows / 1000 trials | full-known smooth rows / 1000 trials | rank gain in either objective |
|---:|---:|---:|---:|---:|
| 351 | 8 | 1.00 | 0.67 | 0 |
| 356 | 12 | 2.54 | 2.19 | 0 |
| 376 | 5 | 0 | 0 | 0 |
| 377 | 6 | 0 | 0 | 0 |
| 385 | 12 | 0.40 | 0.33 | 0 |

The descriptive smooth-yield ordering is suggestive, especially for 356, but
nineteen of twenty rows come from the sparse lane in the exponent-two
configuration, and every row dies in the objective quotient.  Thus it
currently measures easy square-derived norm
relations, not the proposed residual S-class layer.  Five samples, unequal
factor-base widths, and the deliberately different exceptional-block sizes
also rule out a predictive interpretation.  The next useful escalation is a
larger-prime collision or batched short-vector sieve scored only by quotient
rank gain, not more exponent-two smoothness.

PARI's square-class-specific `idealredmodpower(...,2)` was also tested on the
curve-356, bound-5000, 10,000-trial, width-at-most-three lane. It removes the
square-ideal contribution before reduction and produced no smooth row. In
combination with the dependent ordinary-reduction rows, this confirms that
the apparent exponent-two yield is representative-shaping noise rather than
new mod-two class information. The option remains useful as a negative
control, not as the main collector.

The full-ideal special-`q` route now has a stricter large-prime pipeline. Full
integer factorization was observed to stall on the first large norm cofactor,
so hybrid collection instead uses independent bounded trial division,
prime-first proof, retained unresolved cofactors, and an optional exact
batch-GCD pass. Primitive projective power-basis normalization removed 184
apparent cycles from an early diagnostic: 180 had the canonical rational row
`0x3` and four had zero row, all caused by rational multiples of the same
sampled element. After that correction, the larger curve-356 run has 1,079
exact partial edges on 1,355 vertices and nullity zero. A disjoint-ideal run
has 666 edges on 849 vertices, also nullity zero. Batch GCD proved complete
factorizations for 26 additional cofactors in the current smaller lane,
giving 385 edges on 525 vertices with no dependency. Its product/remainder
tree produced the same exact factor splitting as the all-pairs algorithm on
all 3,017 unresolved cofactors, while reducing pairwise fallback comparisons
from 4,549,636 to 18,096; only six composite shared aggregates needed that
fallback. This removes the quadratic GCD front door without weakening the
subsequent proved-primality gate.

Every exact partial edge is now retained, so independent ledgers can be
replayed in one global hypergraph. A fresh two-run curve-356 replay combined
361 projectively distinct edges on 473 vertices; adding a feedback run gave
539 edges on 685 vertices. Both combined matrices still had full row rank and
zero quotient-rank gain. The feedback selector exactly identified the five
repeated small residual ideals at rational primes
`1777,1997,2129,3347,7901` and recycled them in a paired special-ideal
cycle, but its remaining residual vertices were new leaves. Thus the present
bottleneck is measured precisely: increase residual-vertex reuse or true
factor-base smoothness, rather than widen the special-`q` interval or count
raw partials. The natural next engine is a bucketed norm sieve/product-tree
batch pipeline aimed at a controlled residual-prime bound. These are bounded
local checkpoints, not a class-group or Selmer upper bound.

The generic ledger audit and canonical-row augmenter now re-prove and supply
the declared Selmer primes as PARI factor hints. On the merged curve-356
factor-base-1000 checkpoint this avoids repeating the discriminant
factorization wall: all 171 canonical principal rows verify, leaving displayed
dimension 147 after the 30 `S` columns. Since 1,000 is far below the
Bach/ERH generation bound 1,056,719, the classification remains
`UNCERTIFIED_FACTOR_BASE` and the number 147 is not a global upper bound.

Until factor-base generation is separately proved, every materialized
dimension remains a checkpointed relation fingerprint and must not enter the
theorem gate as an upper bound.

The eclib route has a separate structural obstruction.  Its general descent
enumerates binary-quartic coefficients in machine `long` intervals;
`selmer_only` still performs that enumeration, and comparable high-height
fibres require coefficient bounds around `10^29`--`10^32`.  Replacing the
counter type would remove the overflow but leave an impossible exhaustive
loop, so that backend is retained as a small-height control rather than the
record-fibre engine.  A 300-second four-way PARI 2.19 comparison of ideal
powers `1,6,8,12` made the same initial progress in every lane (six trials and
three resolved ideals), so changing `idex` alone does not repair the BNF
architecture.

## Exact parity squeeze

The same replay proves trivial rational 2-torsion and exact root number `-1`
for both 356 and 385.  The 2-parity theorem over `Q` therefore makes the total
2-Selmer dimension odd.  This improves both gates without assuming BSD or the
finiteness of Sha:

- a proved total Selmer upper bound `30` already forces the exact value `29`;
- the dimension beyond the 29 known point classes is even;
- for a rank-32 search the residual quotient modulo MW17 must have dimension
  at least `16`, since its parity is even, rather than merely meeting the raw
  threshold `15`.

The monotone gate accepts this refinement only with explicit parity evidence,
rounds proved upper bounds down to the allowed parity, and rounds target
thresholds up.  Without such evidence its previous behavior is unchanged.

## Exact known residual and rigid quotients

For both record fibres, the first seventeen displayed points are the exact
specialized generic MW17 subgroup and `P18,...,P29` are independent modulo it.
Thus the certified realized residual subgroup is `F_2^12`.

For curve 356 the rigid plane has rows

```text
000000010000
011010010100
```

in `P18,...,P29` coordinates.  Its leftmost-pivot complement is

```text
P18,P20,P21,P22,P23,P24,P26,P27,P28,P29.
```

For curve 385 the rigid rows are

```text
000100000000
000100110000
```

and the complement is

```text
P18,P19,P20,P22,P23,P25,P26,P27,P28,P29.
```

Consequently the known subgroup modulo the rigid plane is exactly
ten-dimensional.  This is a lower-bound statement inside the residual Selmer
quotient.  Once a complete descent returns total 2-Selmer dimension `d`, the
three decisive dimensions are

```text
residual modulo MW17              d - 17
after quotienting the rigid plane d - 19
additional beyond all 29 points   d - 29.
```

If `d=29`, the known twelve classes exhaust the residual Selmer group and the
rigid quotient is exactly the displayed ten-space.  Together with the already
certified trivial rational 2-torsion and saturation/independence gates, this
would supply the desired rank upper bound 29.  No such value of `d` is asserted
here.

## Control fingerprint: `+12` versus `+5`

The exact control dimensions are

```text
curve     351  356  376  377  385   12
gain        8   12    5    6   12   12
```

At the eleven places audited for both record fibres and the `+5` control 376,
the following scalar features separate *both* `+12` samples from 376:

| place | separating exact features |
|---:|---|
| 2 | Kodaira symbol, minimal-discriminant valuation, Tamagawa 2-part, known-residual localization-kernel dimension |
| 13 | known-residual localization-kernel dimension |
| 37 | localization-kernel dimension, nontrivial component-Hilbert pair count |
| 53, 67, 71, 83, 113 | localization-kernel dimension |
| 79, 97, 101 | ambient local Kummer dimension, known image dimension, kernel dimension, selected-block image dimension |

This confirms that 2 is structurally discriminating in the current 356
comparison, but not that it alone creates the jump.  The odd-place pattern is
substantial, including for 385.  With two high samples and one low sample this
is a theorem about the controls, not a sufficient local criterion or a
statistical classifier.

The stacked known-point localization matrix has full source rank on all six
controls.  Deleting any one audited place leaves that rank unchanged.  This
means every individual place is redundant *for separating the certified known
classes in the full stored place set*.  It does not answer which place cuts
down the unknown global Selmer ambient space.  That requested matrix is now
computed by `run_elkies_2026_relative_2selmer_checkpointed.py` after a
certified BNF; until then its ranks remain `UNKNOWN`.

## Pairing and decomposition boundary

The available componentwise Hilbert forms do not descend through the rigid
plane: the exact obstruction primes are `13,23,37,139` for curve 356 and
`5,29,37,41,73,109,127` for curve 385.  Their corestricted local Tate controls
are zero.  Therefore there is no canonical pairing on the ten-dimensional
quotient from these data and no invariant meaning to an “indecomposable ten”
or a split into smaller blocks.  The certificate records this as `NOT_DEFINED`
instead of diagonalizing an arbitrary coordinate complement.

A complete Cassels pairing on the complete residual Selmer group would make
the decomposition question meaningful and is retained as a separate field in
the fingerprint schema when available.

## Parameter strata and CRT prototypes

For every common audited prime the artifact records the exact affine or
infinity-chart residue of each control modulo `p`, `p^2`, and `p^3`, together
with the complete local-feature hash.  These are finite sampled strata; the
data do not prove that a fingerprint is constant on an entire residue
cylinder.

Two exact CRT classes preserve five discriminating `p^3` residues:

```text
356 prototype:
t = 14503794234702288112 + 47438163879590960216*n
places 2,13,37,53,71

385 prototype:
t = 13329277794157146704 + 39863665779809550328*n
places 2,13,37,53,67
```

The congruences are exact; preservation of the local Kummer/Selmer fingerprint
is a search hypothesis.  Every manufactured parameter must be re-audited
before it inherits any arithmetic label.

## Monotone residual sieve

`elliptic-curves/cas/elkies_residual_selmer_gate.py` now separates two kinds of
authorization.

1. A complete unconditional descent may issue the existing exact residual
   gate and support theorem-directed cover/search work.
2. An incomplete sieve stores a sequence of proved residual upper bounds.
   Bounds must be nonincreasing and carry evidence provenance.  A missing BNF
   is stored as “no finite upper bound yet”, not as a numerical estimate.  The
   fibre is rejected as soon as the proved bound is below 15.  Otherwise only
   a point search with explicit height/time/resource limits is authorized.

The point-search entrypoints report their requested limits back to the gate;
an open authorization is rejected if a required limit is absent or the run
would exceed its recorded allowance.

An open sieve always stores `theorem_claim_authorized=false` and
`expensive_search_authorized=false`.  A bounded search can improve the
Mordell--Weil lower bound; it cannot turn missing global descent data into an
upper bound.

## Replay

```bash
sage -python \
  elkies-k3/scripts/build_r17_residual_selmer_fingerprints.sage --check

python3 -m unittest -v \
  elliptic-curves/tests/test_elkies_residual_selmer_gate.py \
  elliptic-curves/tests/test_elkies_relative_2selmer_checkpointed.py
```

The record-pair complete-descent input suite can be regenerated with

```bash
python3 elliptic-curves/cas/build_elkies_2026_relative_2selmer_suite.py \
  --record-pair-only \
  --output-dir \
  artifacts/local/elliptic-curves/r17-074d9-record-residual-2selmer-v1/programs \
  --manifest \
  artifacts/local/elliptic-curves/r17-074d9-record-residual-2selmer-v1/manifest/input.json \
  --overwrite
```

If the PARI 2.19 relation collector produces a `bnfcertify`-validated binary
checkpoint, continue in that same GP build with
`elliptic-curves/cas/run_elkies_2026_pari219_selmer_from_bnf.py`. This avoids
binary-checkpoint compatibility assumptions and emits the requested full
leave-one-place-out matrix before any known-point alignment.

Long BNF/descent outputs remain under `artifacts/local/elliptic-curves/`; only
a completed certified result should be promoted to the generated-results
tree.

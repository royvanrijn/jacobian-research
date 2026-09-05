# What distinguishes a large specialization jump?

The evidence identifies **a large simultaneously rational Kummer subspace**, but
does not yet identify a point-blind specialization condition that forces it.
There is substantial block behaviour in *visibility*, and an exact block
structure in *obstructions*. Neither establishes a common construction of the
new rational directions. “Independent accidents” and “a common arithmetic
source” remain unresolved alternatives: linear independence does not imply
probabilistic independence.

This analysis uses commit `7471645748ae07b8d942f6957e314378ca8aed1f`, including
the terminal R17 `07ca9,-2507/3068` rank-at-least-26 certificate. It does not
consume later changing search exports. It adds exact, cheap retrospective
computations and one frozen 96-evaluation experiment. Mathematical status is
inherited from the repository certificates; no new rank or upper bound is
claimed and no search selector is changed.

## Panel and accounting

The [CSV](../../artifacts/generated-results/elliptic-curves/rank_jump_comparison_panel_v1.csv)
has 91 main observations, three earlier MW12 controls, and two explicitly
unresolved-parent records. The 91 main rows represent **89 distinct curves**:
MW16 `-34/87` repeats ICARM 548 and MW16 `-1905/52` repeats ICARM 542. The latter
worker observation certified 25, while the historical curve certificate
certifies 26. These are observation rows, not independent statistical trials.
The CSV preserves those separate lower bounds and observations.

The main panel includes all 36 discoveries in the pinned v2 index, the new
rank-26 initial certificate, the completed MW16 pilot/wider control cohort,
matched R17 controls, eleven historical soluble-subspace records and fifteen
historical R17 Kummer records. All high-gain rows receive a
[per-fibre report](FIBRE_REPORTS.md); full coordinates, signatures, exact
independent-subset indices and subgroup labels are retained in the
[portable inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_retrospective_inputs_v1.json)
and [report](../../artifacts/generated-results/elliptic-curves/rank_jump_retrospective_report_v1.json).

For a specified reference group (M\subset A=E(\mathbf Q)), put
(G=\delta(M)), (W=\delta(A)/G), and (\epsilon_M=\dim(A/M)[2]).
The [canonical residual exact sequence](../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md)
gives

\[
 \dim W=\operatorname{rank}(A/M)+\epsilon_M,
 \qquad \dim(\mathrm{Sel}_2/G)=\dim W+\dim\Sha[2].
\]

On all 91 main rows, absence of rational 2-torsion and full generic Kummer
rank certify **epsilon_M=0**. The new replay checks these facts using
irreducibility witnesses and split-root characters; the fifteen earlier
signature-only records replay their supplied exact matrix algebra rather than
recomputing their original points. For the other 76 rows, all supplied point
coordinates are re-evaluated on the short model. Exact full-rank signatures
select an independent subgroup containing the marked generic basis. Its
quotient rank is exact; the rank of the *whole curve* remains a lower bound.
A zero recovered quotient is a censored control, never a proved ordinary rank.

| Historical fibre | Specified reference rank | Certified free quotient |
|---|---:|---:|
| Published R17, ranks at least 25/26/27/28 | 17 | 8 / 9 / 10 / 11 |
| ICARM 356, 385, 543 | 17 | 12 each |
| ICARM 398 | 16 | 14 |
| ICARM 400 / 401 / 542 / 548 | 16 | 12 / 11 / 10 / 8 |
| Fermigier E22 / Fermigier rank-20 control | 12 | 10 / 8 |
| Nagao section 7, T=5081/47 | 12 | 8 |

The earlier MW12 examples expose an important different kind of “block”.
Recomputed integer Smith forms give:

- Fermigier E22: (L/M\cong\mathbf Z^{10}\oplus(\mathbf Z/2)^{11}\oplus\mathbf Z/12).
  Thus ((L/M)/2) has dimension **22**, not a free jump of 22.
- Fermigier rank 20: the displayed quotient is free of rank 8.
- Nagao section 7: (L/M\cong\mathbf Z^8\oplus(\mathbf Z/2)^{11}).
  Its mod-2 tensor quotient has dimension **19**, not a free jump of 19.

Here (L) is the certified displayed lattice, not the full MW group. Full
saturation is not supplied by these Smith forms. The conspicuous twelve- and
eleven-dimensional divisibility blocks must not be relabelled new rank.
Fermigier's K11 graphic matroid comes from taking weight-two combinations of
an already supplied ten-direction basis; that combinatorial structure is
not a pre-point explanation of its existence.

ICARM 273/302 have lower bounds 30/31. Their arbitrary rank-17 masked coordinate
subgroups leave displayed quotients 13/14, but are **not certified generic
subgroups of a recovered parent family**. No “31=17+14 specialization” claim
is made for 302. Prefix-d17 recovery of 10/13 versus 3/14 is evidence about
chosen coordinates; changing the mask changes recovery. Promoting anchor
points by base change likewise changes (r+j=(r+k)+(j-k)), not the fibre rank.

## What the quotient blocks actually show

**Kummer and local support.** All **94 equation-bearing rows**, including the
three MW12 rows, have irreducible 2-division cubic and nonsquare discriminant,
hence Galois group S3. This exactly rejects a cubic splitting, cyclic-cubic
transition or rational-2-torsion event as the explanation of these hits.
It does not identify their number fields with each other, nor exclude subtler
Galois structure of higher covers.

In all **76** main rows with re-evaluated coordinates, **every selected good
auxiliary prime** has

\[
 \operatorname{im}(M\to E(\mathbf F_p)/2)
 =\operatorname{im}(\langle M,\text{supplied points}\rangle\to E(\mathbf F_p)/2).
\]

Nevertheless their combined prime fingerprints distinguish up to fourteen
extra global directions. There is no contradiction: expressing a point in
the generic image separately at each prime permits different coefficient
words; membership in a global generic span requires one word working at
all primes. This is a reason not to assign each extra direction a unique
“support prime”. These good-prime characters are **not** a bad-place local
Selmer calculation, ideal-factorization support, or a completed descent.

**Cubic ideals and descent fields.** Every direction on a fibre naturally
uses its one cubic algebra. Norm-square identities and common field membership
are therefore baseline properties of Kummer descent. The recent high/low
pairs have no completed ideal-class/S-unit presentation or full Selmer basis
in this audit. Those fields are explicitly UNKNOWN, not empty blocks. The
reassessment also explains that generic sections themselves contribute cubic
class information; subtract the marked generic image and component/valuation
corrections before interpreting class pressure as exceptional incidence.

**Cassels–Tate blocks.** Exact symplectic reduction of the seven existing
fixed-field matrices gives the following noncanonical decompositions; H is a
nondegenerate alternating plane over F2 and 0 denotes the restricted radical.

| u | Locally admissible inherited dimension | Restricted CT decomposition | Witnessed soluble dimension |
|---:|---:|---|---:|
| -3 | 17 | H^8 + 0^1 | 0 |
| -2 | 13 | H^6 + 0^1 | 0 |
| -1 | 18 | H^8 + 0^2 | 0 |
| 0 | 20 | 0^20 | 20 |
| 1 / 2 | 13 | H^6 + 0^1 | 0 |
| 3 | 15 | H^7 + 0^1 | 0 |

The zero witnessed dimensions are lower bounds. These are *transported class
subspaces*, not full Selmer groups. The independently known point at u=-1
lies outside that inherited span. The nonsquare-B pencil has arithmetic generic
rank zero: u=0 is itself an exceptional fibre of that pencil, not a transported
generic rank-20 group. This comparison controls the cubic field exactly while
showing that simultaneous solubility can disappear in many dimensions.

At u=-1, rank 16 excludes **262140** of 262144 inherited classes by pairing;
only a two-dimensional restricted annihilator survives. Symplectic pairs are
obstruction blocks, not independent geometric causes of those obstructions.
On the historical point-built high-gain subspaces the CT rows vanish against
every Selmer class because rational witnesses exist. That vanishing is
mathematically forced and not a point-blind classifier trained on independent
class samples. Nonzero pairing detects Sha; a zero radical can still hide
4-divisible Sha. See [Fisher, Introduction and Theorem 3.1](https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf)
for the pairing and its relation to higher descent.

**Height Gram structure.** For the three recent rank-25 curves, project out
the generic span with (S=H_{QQ}-H_{QM}H_{MM}^{-1}H_{MQ}), using the retained
25-by-25 decimal height matrices. Exact rational arithmetic on those decimal
inputs gives no zero off-diagonal entries. Largest absolute correlations are
approximately 0.536, 0.588 and 0.670 for R17 505/794, MW16-05 307/206 and
MW16-04 -1647/91. At threshold 0.5 there are several small graph components;
at 0.75 every component is a singleton. Those thresholds and components are
basis-dependent diagnostics, not a certified orthogonal lattice factorization.
The projected-to-raw point-height ratios range approximately 0.14–0.38.
Much height lies in generic projections. This supports compression/visibility,
not a causal rank predictor. These are numerical height inputs, not rigorous
intervals, and do not rule out an orthogonal decomposition after changing basis.

## Initial half-lattice incidence maps

The new audit processes all returned points in **344 retained initial-wave
charts** (eight fibres, 43 each), rechecks membership, and computes their finite
Kummer signatures modulo the generic image. It retains every chart's parity
mask, observed signature set, span, prefix contribution, coverage, and every
pairwise span intersection. Masks of quotient classes use a deterministic
independent reference basis. A matching finite signature does not prove a
rational relation for an arbitrary returned point; the map is injective on
the certified reference span only. Independent signature images do certify
independent directions. “31 signatures” is not an assertion that the returned
cloud has exactly 31 global classes.

| Recent fibre | Certified extra rank | Maximum single-chart image rank | Greedy charts spanning observed image | Initial prefix needed |
|---|---:|---:|---|---:|
| R17 07ca9, 505/794 | 8 | 5 | 30, 21 | 8 |
| R17 07ca9, -2507/3068 | 9 | 5 | 1, 12 | 6 |
| MW16-05, 307/206 | 9 | 4 | 43, 3, 13 | 43 |
| MW16-04, -1647/91 | 9 | **9** | **36** | 12 |
| Published R17, -2300/843 | 7 | 2 | 25, 33, 39, 2 | 33 |

Chart indices here are one-based. The greedy sets are sufficient, not asserted
minimal. A chart late in the wave may see directions already admitted earlier.
The worker's reference state grows during this initial wave, but the centres
remain combinations of the original generic basis: newly adjoined coordinates
are padded with zero in centre words. This is not an experiment running every
chart against an independently frozen initial state.

One chart seeing nine directions is strong **visibility clustering**, but it
is not nine soluble exceptional covers sharing one Selmer label. For Q in M,
the pointed quartic's birational map \(\phi_Q:C_Q\simeq E\) can reach all rational
points. Its **2-cover map** is \(\pi_Q=2\phi_Q-Q\), whose residual Kummer class
is zero. Under the chart involution, \(R\mapsto Q-R\), so \(\bar R\mapsto-\bar R\)
in A/M and the mod-2 quotient signature is unchanged. The two roots above a
slope therefore supply the same direction, up to sign. This is also proved in
[the canonical chart/cover distinction](../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md).
A deep generic parity mask labels a midpoint/coordinate system, **not an
exceptional residual torsor**.

Hamming-distance bins versus span intersections are retained for every chart
pair, but no arithmetic clustering is inferred from them: Hamming distance
changes with the marked generic basis. Likewise ordinary point-cloud density
can reflect combinations and the involution. Large adaptive waves may revisit
the same finite quotient span after the initial charts have exposed its visible
part. The retained 301-chart zero-gain follow-ups on both new MW16 rank-25 curves
are consistent with this explanation, but cannot exclude further invisible
rational directions or prove that all their spans are exhausted globally.
The existing historic rank-28 oracle exactly reconstructs visibility of its
forty targets and eleven directions; it remains an oracle audit.

## Three high-versus-low paired case studies

Coverage below is the mean completed denominator fraction, not an area fraction,
a point-detection probability, or an exposure correction. All six rows have
certified full generic subgroups and no observed new local image at any tested
good prime. Every lower endpoint is a search observation without a rank upper
bound.

| Pair | High / low parameter | Extra rank observed | Short coefficient bits | Score | Mean coverage |
|---|---|---|---|---|---|
| MW16-05 | 307/206 / -3158/1291 | 9 / 0 | 224 / 265 | 93.224 / 94.196 | .684 / .670 |
| MW16-04 | -1647/91 / -2177/2397 | 9 / 0 | 210 / 268 | 96.691 / 93.614 | .754 / .761 |
| Published R17 | -2300/843 / -1561/3133 | 7 / 0 | 250 / 248 | 94.548 / 93.955 | .704 / .726 |

**MW16-05:** the lower-gain fibre actually has the higher retained score. Family,
prime table and 43-chart policy agree, but the parameter boxes differ (1024
versus 4096) and the low-gain coefficients are 41 bits larger. The high fibre's
422 returned point records expose 27 nonzero finite quotient signatures of
rank nine; the low fibre's 172 records expose none. No CT or ideal-class
comparison distinguishes incidence from global solubility on this pair yet.
Thus this is a score-retention counterexample for the measured endpoint,
not proof that the score inversely predicts true rank.

**MW16-04:** both fibres are in the same 4096 cohort. Denominator coverage is
nearly equal, but coefficient sizes differ by 58 bits. The high fibre gives
57 nonzero finite signatures of rank nine, with all nine visible in chart 36;
the low fibre's 164 returned records add zero detected quotient directions.
Neither a different E[2] representation nor a good-prime generic-image deficit
separates them. The coefficient mismatch still permits a visibility explanation
for the null control. Initial masks and a high score do not supply a common
soluble auxiliary curve.

**Published R17:** this is the strongest scale control. Within the same cohort,
choose the zero-observed-gain row nearest in coefficient bits to the high row,
then nearest score, then lexical parameter. It has 248 versus 250 bits, similar
parameter height (3133 versus 2300), slightly greater mean coverage, and the
same chart count. Seven detected directions versus zero therefore cannot be
explained solely by gross coefficient size or fewer completed chart plans.
Specific phases, coordinate distortion and unseen solubility still differ.
This pair is a suitable future arithmetic-descent comparison; it does not
supply a proved zero-jump fibre.

A fourth structural comparison is u=0 versus u=-2 in the fixed-field pencil:
the cubic field is identical, but inherited local dimensions are 20 versus 13
and restricted CT ranks are 0 versus 12. This separates local survival from
simultaneous global solubility far more directly than the search-score pairs.
It is not a small change inside the original generic MW17/MW16 family.

## A sufficient construction and a falsifiable test

A precise possible mechanism is a **common rational carrier**. Suppose a finite
base map (C\to\mathbf P^1\) and k sections on the pullback are defined before
exceptional fibre points are supplied. If those sections are independent modulo
the pulled-back generic subgroup, a rational point of C specializes all of them
together. At a fibre where the combined specialization remains independent,

\[
 t\in\operatorname{image}C(\mathbf Q)
 \quad+\quad\text{certified independent specialized sections}
 \quad\Longrightarrow\quad\text{a soluble k-dimensional block}
 \quad\Longrightarrow\quad\operatorname{rank}E_t\ge r+k.
\]

The independence clause is essential. Several bisections can have dependent
images or differ by generic sections; a single rational lift does not certify
k new directions. A positive-rank genus-one carrier or a rational carrier with
many rational points would make the condition usable; a high-genus fibre
product with no known rational points merely relocates the solubility problem.
A common squareclass (d(t)) in equations (u_i^2=d(t)f_i(t)^2) is a concrete
way multiple covers could share a lift. Degree and genus of the *parameter
carrier* matter. Smooth genus-one 2-coverings of a nonsingular elliptic curve
themselves do not become unions of rational components or drop genus; any
such event must concern a different auxiliary incidence model or a degenerate
fibre, with its maps checked.

The [frozen experiment](EXPERIMENT.json) selects sixteen bisections by lexical
label from the generic equation atlas, without scores or exceptional points.
It compares all 120 pairs of discriminant squareclasses in Q(t), then makes
96 exact evaluations at three already studied high/low R17 pairs. Products
of two nonzero discriminant polynomials are tested for being exact rational
polynomial squares, which is equivalent to equality of their Q(t) squareclasses.
The primary gate requires at least three bisections on a common nontrivial
quadratic carrier splitting at a high fibre; any positive would then need
independent section images.

**Result:** all sixteen generic squareclasses are distinct, and none of the
96 specializations is a rational square or degenerate. The precise
sixteen-bisection shared-carrier hypothesis is **refuted**. This is an exact
finite statement, not a refutation of other carriers or of coordinated global
solubility. The complete earlier degree-two atlas and sampled degree-three/four
extensions already explain only **5/8, 3/9, 2/10, 1/11** of the four published
R17 quotient ranks; the higher-degree samples add none. The new test does not
expand that atlas or launch a prospective search. The ignored large atlas was
loaded by the hash bound in the committed deep-cover certificate; an initial
capture adapter incorrectly expected it to be tracked directly and was fixed
before any evaluations. The compact captured input makes replay independent
of that large local file.

The **strongest observed arithmetic mechanism** remains coordinated survival
of global obstructions, exhibited by the fixed-field pairing contrast. The
common-carrier construction is a sufficient, testable route to such survival,
but currently has no witness large enough for the observed +8 to +14 fibres.
The missing implication is not “many points are visible”: it is

\[
 \text{point-blind specialization condition}
 \ \not\Rightarrow_{\text{currently proved}}\
 \text{large independent subspace of globally soluble residual classes}.
\]

Even a large point-blind CT annihilator would leave the further implication
“surviving classes are rational rather than higher-divisible Sha” unresolved.

## Ranked conclusions and information for Agent 1

1. **Most plausible structural mechanism: coordinated global solubility.**
   The fixed-field experiment proves that many local classes can become
   obstructed together while E[2] stays fixed. A common carrier or other
   low-degree identity supplying several independent sections is the strongest
   sufficient construction to seek. Neither is yet a prospective predictor
   for the recent hits.
2. **Established secondary mechanism: correlated visibility.** A few charts
   span the observed quotient, sometimes a single chart spans all nine.
   Height projection and midpoint phases help explain recovery waves; they
   do not explain existence or imply any additional rank after a null follow-up.
3. **Established accounting mechanism: specialization divisibility.** The
   MW12 Smith examples show large 2-primary quotient blocks without equally
   large free-rank jumps. This must be separated from soluble exceptional
   directions before interpreting a descent dimension.

The following explanations are disproved in their stated scope or weak:

- **Cubic splitting/Galois drop:** absent on all 94 equation-bearing rows.
- **One new good-prime support per direction:** absent on every selected good
  prime in all 76 re-evaluated main rows; bad-place behaviour remains open.
- **Preserving the descent field preserves rank/solubility:** directly
  contradicted for the transported class space by the fixed-field CT data;
  the pencil's generic rank is zero.
- **The retained degree-2–4 atlas explains the whole jump:** false on the four
  published R17 controls. The new sixteen-carrier subhypothesis also fails.
- **Different retained fibrations contribute overlapping generic blocks:** the
  prior 384-pair cross-family incidence certificate found no additional generic
  directions on the first 32 discoveries; duplicate R17 presentations generate
  the same integral subgroup. Unknown families remain outside that exclusion.
- **Small quartics, many points, a single-chart burst, or a large odd CT radical
  are rank predictors:** none establishes incidence; the last can be parity
  forced, while the others depend on presentation or already known points.

Missing computations/theorems, in priority order:

1. An equation-only residual Selmer/S-unit space on a matched high/low pair,
   with exact generic-image removal and full bad-place local conditions.
   Until then there is no pre-point incidence comparison for the recent hits.
2. CT pairings against independently constructed complement classes, or a
   certified Selmer upper envelope plus a nondegenerate minor. Point-built
   zero rows cannot substitute. Complete descents are currently costly and
   were deliberately not relaunched here.
3. Higher-descent control of surviving Sha or explicit rational witnesses on
   several independently identified covers. Restricted radicals alone fail.
4. A common auxiliary carrier with verified maps and at least three independent
   pullback directions, then an explanation of how its block could grow to
   eight or fourteen. The fixed sixteen-cover dictionary is now a negative
   regression, not a queue to enlarge automatically.
5. Full cubic ideal/S-unit labels on the recent pairs, and rigorous height
   reduction if a basis-invariant lattice splitting is proposed. Common field
   membership and decimal correlation graphs do not provide these.

Potential information for Agent 1 is limited to correctly typed quantities:

| Proposed information | Layer | Admissible interpretation |
|---|---|---|
| Residual locally admissible dimension after generic/saturation corrections | incidence | Selmer capacity; no point claim |
| Equation-only common-carrier condition with independent pullback sections | incidence | A construction capable of producing several directions at rational lifts |
| Exact rational lift / cover witness, followed by independence | solubility | Certifies a soluble block and its rank lower bound |
| Rectangular CT obstruction rank and annihilator dimension | solubility | Excludes classes; a surviving kernel is only a necessary condition |
| Certified Selmer upper envelope minus obstruction rank | incidence | An unconditional upper-rank exclusion when all hypotheses hold |
| Coefficient size, chart cost, midpoint phase, height projection, mask overlap | visibility | Coordinate/exposure diagnostics only; phase involving exceptional points is oracle-only |
| Nagao/Mestre score | incidence (heuristic) | Curve-level heuristic, neither a rank nor solubility certificate; this panel does not validate it |

No item is promoted into Agent 1's selector. The most useful eventual handoff
would be a compact, point-blind certificate of a shared soluble construction,
or an unconditional obstruction upper bound. The present handoff is the paired
input panel, exact chart incidence maps, negative regressions, and the precise
remaining global-solubility gap.

## Reproducibility and limits

The [README](README.md) contains the five narrow replay commands. Six small
arithmetic regression tests independently check finite-character homomorphisms
(including rational 2-torsion reductions and reduction to O), Smith forms via
2-by-2 determinantal divisors, alternating kernels by exhaustive enumeration
through dimension four, finite-signature aliasing, and polynomial-square gates.
Deterministic report replays passed. These are not formal proofs or reruns of
the historic expensive descents/searches. Prior group-law embeddings in the
MW12 source are inherited; their Smith algebra is recomputed.

All created files are rank-jump-specific. No existing active proof, search
script, protocol, worker setting, candidate population, mathematical-status
entry or Agent 1 output was modified. The full chart transcripts are bound by
hash, while compact projections retain exactly the data needed for this audit.
Newly changing search results can later enter a separately versioned panel.

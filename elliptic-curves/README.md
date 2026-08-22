# Elliptic-curve rank/conductor search

This directory is an independent research track for finding an elliptic curve
over `Q` satisfying either

- `rank(E(Q)) >= 21` and `ln N(E) < 182.72`, where `N(E)` is the arithmetic
  conductor and `ln` is the natural logarithm; or
- `rank(E(Q)) >= 30`.

The inequality is strict.  A heuristic rank score, an analytic-rank estimate,
or a list of points that has not been checked for independence does not satisfy
either target.

## Current status

**Second target met.**  The 2026 public ICARM curve 273 has 30 displayed
rational points whose exact finite-reduction images have full binary rank 30.
The repository's primary checker and an independent Sage implementation prove
`rank(E(Q)) >= 30` unconditionally.  No unconditional exact-rank-30 statement
is claimed.  The first, low-conductor rank-at-least-21 target remains open.

The research program now has exact-arithmetic implementations of the proposed
pipeline

```text
local discriminant roots -> Hensel lifting -> CRT -> 2D lattice reduction
                           -> exact conductor calculation -> rank work
```

in both the Mestre--Fermigier generic-rank-12 family and Nagao's published
generic-rank-13 base change.  It also contains the rank-nine
Elkies--Klagsbrun K3 family as a smaller local-arithmetic test bed and an exact
replay of Nagao's published rank-21 curve.

The family-design branch now also contains a new six-root Mestre family with
centers `(0,23,93,128,133,175)`.  Six nonvisible linear-abscissa square
identities are verified symbolically.  After
`T=(14406-u^2)/(2u)`, one chosen companion and split infinity raise the exact
generic lower bound to 13.  Its primitive base-changed discriminant frontier
has degree 40, versus degree 398 for the reproduced Kihara rank-14 template.
This is a proved generic construction and materially better conductor
geometry, although its own bounded specialization scans have not produced the
strongest six-root fiber.  A second split-infinity family, with centers
`(0,25,95,143,168,205)`, also has generic rank at least 13.  Its specialization
`u=197`, `T=337/394` has an exact rank lower bound of 17 and exact
`ln(N)=173.594891144976...<182.72`.
Bounded integer specialization searches in this and a companion family found
two stronger exact frontiers: family 2 at `u=483` has rank at least 19 and
`ln(N)=157.759935999987...`, while family 3 with centers
`(0,43,128,197,231,289)` at `u=660` has rank at least 19 and
`ln(N)=164.053646218834...`.  Their conservative `Delta=11/5`
explicit-formula values are respectively `20.3010...` and `20.4139...`;
since both root numbers are `-1`, GRH would bound each analytic rank by 19.

Evidence is classified as follows throughout this directory:

| Label | Meaning |
| --- | --- |
| **Theorem** | A proved mathematical statement, with a proof or a cited source. |
| **Verified computation** | An exact or independently replayed calculation, with its command and software assumptions recorded. |
| **Bounded experiment** | A finite search with stated parameters; it proves nothing outside the enumerated set and a score is not a rank certificate. |
| **Open problem** | A conjecture, proposed method, or unfinished certification step. |

### Benchmarks

- **Theorem (published) and new exact replay:** Fermigier's curve `E22` has
  rank at least 22.  In addition to checking all printed coordinates, this
  repository now transports them to the normalized short model and gives a
  finite-reduction independence certificate; the lower bound no longer rests
  only on the published regulator computation.
- **Verified computation:** its conductor is

  ```text
  22720638514787473197194583889675055980109503436060704437972911338086049759883790
  ```

  and `ln N = 182.724910950637...`.  It therefore misses the strict
  `182.72` threshold by about `0.00491`; it is a calibration benchmark, not a
  hit.
- **Published record data and independent exact replay:** the 2026 ICARM curve
  273 has 30 independent rational points.  The primary checker exhaustively
  constructs the finite quotients `E(F_p)/2E(F_p)`; a separate Sage replay uses
  invariant factors and discrete logarithms.  Both obtain binary rank 30, and
  a modulo-23 witness excludes rational 2-torsion.  This proves the second
  target unconditionally as a rank lower bound.
- **Published low-conductor data and independent exact replay:** ICARM curve
  245 has twenty rational points whose finite-reduction images have binary
  rank 20.  Its exact conductor has `ln N=150.668907152237...<182.72`, making
  it the smallest-conductor exact rank-at-least-20 curve currently recorded
  here.  This improves the one-point-short frontier but is not a target hit.
- **Historical record replay:** the 2024 Elkies--Klagsbrun curve has 29
  independent rational points.  Its exact rank 29 is conditional on GRH; the
  unconditional lower bound 29 does not use GRH.
- **Theorem (published) and verified replay:** Nagao's 1994 curve has 21
  independent rational points.  The repository checks the printed model and
  all 21 points exactly, and PARI gives `ln N=196.679545735892...`; it is a
  family and score calibration point, not a target hit.

The rank-29 curve is also strong evidence for the conductor-engineering idea:
its discriminant contains many repeated powers of small primes, while its
conductor contains each of its 17 bad primes only once.  This observation is
motivation, not evidence that forcing prime powers will increase rank.

The record-fiber replay uses the rational isomorphism
`X=36*x+3`, `Y=108*(2*y+x)` to an integral short model.  Exact images in
`E(F_p)/2E(F_p)` at 22 good primes have full binary column rank 29, and the
reduced 2-division cubic at `p=67` has no root.  Infinite descent therefore
proves the 29 public points independent without numerical heights or GRH.
Two complete direct-search tiers then enumerated 1,647 affine x/slope charts
each, the deeper at height 50,000 (and ten million in the 29 x-offset charts),
without finding a nonpublic image.  A separate search constructed all 4,060
weight-two/three subgroup-sum covers and completed 448 charts on the best 64,
again finding no nonpublic image.  These are bounded negative searches, not a
rank upper bound.  A disjoint higher-weight pass then scored 2,000 signed
representatives across five weight bands, searched 50 charts on ten retained
covers, and found 57 exact images, all public points or exact companions in the
rank-29 subgroup.  Those negative searches remain useful historical
calibration, but curve 273 now supplies the thirtieth independent point on a
different curve.  Reconstruction of the unpublished rank-17 K3 fibration
remains relevant for understanding the record-search mechanism and for
producing new high-rank families.

Two denominator-normalized modular sieves subsequently searched disjoint
neighborhoods of the 29 public abscissas and 32 exact nonpublic subgroup
companions.  Together they covered exactly `55,267,250,510` primitive rational
abscissas.  Modular quadratic-residue filtering left 217,309 candidates for
exact integer-square tests, and none was a point.  This is a very large bounded
negative result, not a rank upper bound.  A pinned public-source audit through
2026-08-14 likewise found no reproducible public rank-30 curve or release of
the missing K3 construction data; the audit predates the 2026 curve and is now
retained only as a time-bounded historical snapshot.

### Current certified and bounded-search frontier

The first row is the rank-at-least-30 target hit. Every row labelled “exact” has an unconditional,
portable rank lower bound: exact reduction maps into products of
`E(F_p)/2E(F_p)` give full-column-rank binary matrices, and a separate good
prime proves that the curve has no rational 2-torsion.  “Numerical rank” in
the remaining rows still means only the stable rank of a floating-point
canonical-height matrix.

| Family/search | Parameter | `ln N` | Current rank evidence |
| --- | ---: | ---: | --- |
| ICARM public record | curve 273 | `339.347931713664...` | **exact unconditional rank at least 30**; second operational target met, independently replayed; no unconditional exact-rank claim |
| Fermigier benchmark | normalized `T=39508/39` | `182.724910950637...` | **exact unconditional rank at least 22**; misses the strict conductor bound by `0.004910950637...` |
| ICARM low-conductor record | curve 245 | `150.668907152237...` | **exact unconditional rank at least 20**; independently replayed finite-reduction certificate; one point short of the target |
| Fermigier--Mestre adapter | `u=28917/20` | `159.934825225525...` | **exact unconditional rank at least 20**; imported 58-abscissa search and finite-reduction certificate, one point short of the target |
| Nagao section-7 family | constructor `T=5081/47` | `174.249816228548...` | **exact unconditional rank at least 20**; full mod-2-cover and skew searches found no 21st direction |
| Split-infinity Mestre family 2 | `u=483`, `T=-8441/42` | `157.759935999987...` | **exact unconditional rank at least 19**; GRH-conditional analytic closure at 19 |
| Split-infinity Mestre family 3 | `u=660`, `T=-12655/44` | `164.053646218834...` | **exact unconditional rank at least 19**; GRH-conditional analytic closure at 19 |
| Nagao rank-21 family | constructor `T=6793/64` | `158.572648489303...` | **exact unconditional rank at least 19**; alternate-cover and skew searches remained in that subgroup |
| Nagao rank-21 family | constructor `T=6629/174` | `154.795114152374...` | **exact unconditional rank at least 18**; historical-specialization replay |
| Nagao rank-21 family | constructor `T=3137/72` | `149.535359251913...` | **exact unconditional rank at least 18**; exhaustive small-denominator search certificate |
| Nagao rank-21 family | constructor `T=5783/16` | `177.228672913042...` | **exact unconditional rank at least 18**; exhaustive small-denominator search certificate |
| Nagao section-7 global scan | constructor `T=599/2` | `124.061012256948...` | **exact unconditional rank at least 17**; lowest-conductor certified rank-17 frontier |
| Nagao rank-gain search | `u=135/2`, `T=5065/36` | `144.927455914577...` | **exact unconditional rank at least 17**; 74 deeper skew-box points all replay in that subgroup |
| Nagao rank-gain search | `u=471/11`, `T=5579/22` | `146.678928806750...` | **exact unconditional rank at least 17**; 17-point reduction certificate |
| Nagao rank-13, integer scan | `u=42`, `T=3631/14` | `148.621053634068...` | **exact unconditional rank at least 17**; broad skew/chart search found no 18th direction |
| Nagao mutation search | `u=74`, `T=9037/74` | `151.423206831026...` | **exact unconditional rank at least 17**; 17-point reduction certificate |
| Mestre roots `(0,25,95,143,168,205)` | `u=197`, `T=337/394` | `173.594891144976...` | **exact unconditional rank at least 17**; the split-infinity family has generic rank at least 13 and the pinned certificate is four points short of the target |
| Nagao rank-13, integer scan | `u=84`, `T=2749/28` | `139.773456475157...` | stable numerical rank 16 through height `10^6`; effort-zero `ellrank` timed out |
| Nagao rank-13, local CRT | `u=118`, `T=4813/118` | `128.027255994266...` | 43 nonvisible images at height `10^6`; stable numerical rank remains 15 |
| Nagao rank-13, integer scan | `u=50`, `T=421/2` | `89.115263351204...` | PARI effort zero returned computational bounds `[13,13]`; eliminated as a target |
| Mestre roots `(0,6,47,55,70,80)` | `T=8` | `82.351544058010...` | **exact unconditional rank at least 13**; three accidental quartic points raise the visible span from 10 to 13 |
| Mestre roots `(0,4,30,31,39,46)` | `T=5` | `79.729318123910...` | **exact unconditional rank at least 12**; separate global and CRT searches found no target signal |
| Fermigier record class | `T=1666/9` | `128.959882907388...` | stable numerical rank 16 through height `10^6`; long-cutoff score plateaued |
| Fermigier expanded CRT | `T=644/87` | `153.964400023010...` | stable numerical rank 12 through height `50000` |
| Nagao rank-21 neighborhood | constructor `T=6041/198` | `170.765123121845...` | its six visible sign-pairs have numerical rank 11; height-`50000` search found no new sign-pair |

The smallest-conductor curve currently certified to have rank at least 17 is
the section-7 global-scan specialization `T=599/2`, with global minimal model

```text
y^2 + x*y + y = x^3 + x^2 - 399081217162040000565*x
                         + 2932271975048532129776925923547,
```

with

```text
N = 756855624206125019617445466192515119979543579877522270
rank >= 17.
```

Its 17 exact points and finite-reduction matrices are stored in the section-7
global-search artifact.

The strongest conductor-qualified rank frontier now has three exact rank-at-
least-20 fibers.  ICARM curve 245 has the smallest conductor, with
`ln N=150.668907152237...`; its twenty exact rational points have full rank in
finite-reduction quotients.  The independent certificate and bounded
next-point search are documented in
[`ICARM_CURVE245_RANK20.md`](notes/ICARM_CURVE245_RANK20.md).

The imported Fermigier--Mestre adapter specialization `u=28917/20` remains a
separate exact rank-at-least-20 anchor at `ln N=159.934825225525...`.
The canonical record
[`elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json`](../artifacts/generated-results/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json)
uses `fermigier-mestre-v1:u=28917/20` as the stable identity and records
`T=28917/10` only as an alias.  It unifies all exact model transformations,
the full discovered pool and selected basis, certificate provenance, bounded
saturation status, and the experiment-level promotion/rejection ledger.
A `Delta=11/5` explicit-formula replay has conservative value
`21.033532822984...<22`; its root number is `+1`, so under GRH parity bounds
the analytic rank by 20 and BSD+GRH predicts exact rank 20.  Unconditionally,
this remains only a rank-at-least-20 statement.  The fixed fiber is therefore
deprioritized; the structural priority is to preserve one or more exceptional
quotient directions under a base change.

The independent section-7 specialization `T=5081/47` has
`ln N=174.249816228548...` and the same exact rank lower bound.  All
`2^20-1` nonzero classes of its certified subgroup were scored for
alternate-cover search; 220 completed bounded chart/skew calls found 224
decontaminated points, every one exactly dependent in the same rank-20
subgroup.  Its sharpened explicit-formula calculation also conditionally
closes the fixed fiber at rank 20.  Neither conditional statement is an
unconditional algebraic-rank upper bound.

At the generic level, exact counts of the section-7 elliptic K3 surface over
`F_29` and `F_(29^2)` prove that its Mordell--Weil rank over `Q(T)` is exactly
12 and its geometric rank over `Qbar(T)` is either 12 or 13.  The same symbolic
audit classifies every polynomial quartic abscissa through degree five: the
extra six linear and three quadratic formulas are all dependent on the
rank-12 basis.  Exceptional specialization rank, rather than a missing
`Q(T)`-rational section, is therefore the relevant mechanism for this family.

The same question is now settled for the Fermigier--Mestre adapter.  Exact
counts over `F_41` and `F_(41^2)` prove unconditionally that its arithmetic
generic rank over `Q(u)` is exactly 12, while its geometric rank is 12 or 13.
Thus E22 and the low-conductor rank-20 anchor have exceptional quotient ranks
10 and 8.  Exact comparison of those quotients found no low-genus shortcut:
all 88 affine transports are genus 2; the complete quadratic and projective
Möbius pencils contain no genus-0/1 member; and all 3,160 genuine two-direction
fiber products are connected genus-9 covers with genus-5 third quotient.
An independent quotient-ball audit then enumerated every signed support-at-
most-two representative at both anchors and found all 25,600 additional
affine interpolants irreducible of genus 2.  This closes those finite ansatzes
and that representative ball, not arbitrary higher-degree base changes.  A
complete bidegree-(2,1) classification of all 80 independent pairs likewise
found only genus-2/3/4 rational components; all 80 residual degree-32 factors
are irreducible over `Q`.  A genuine simultaneous-square search on all 3,160
pair covers through projective height 200,000 returned only the two prescribed
anchors.  A separate projective-height-1,024 rational-point sieve on the
pilot's irreducible degree-32 discriminant curve found no point or known-line
intersection.  Rational points beyond those boxes on the discriminant and
genus-9 cover curves remain open.

The smallest-conductor member is `u=135/2`.  It produced 66 unexpected
quartic abscissas by height `10^6`; a subsequent ten-box search through
denominator 128,000 found 74 more exact images, but exact group-law relations
put all of them in the certified rank-17 subgroup.  Its 76 determinant-one
Möbius charts found nothing beyond those boxes.  The nearby `u=471/11`, the
original `u=42`, and the second-generation mutation `u=74` show that the
rank-gain objective is not a one-candidate accident.  For all four curves,
small-prime saturation reduced the height determinant by `2^32`, the square
of the recorded within-span index `2^16`.

The saturation routine's finite-index premise is no longer part of the rank
claim.  Every returned coordinate is checked exactly, and the independent
finite-reduction certificate proves rank at least 17 directly.  For `u=42`,
ten disjoint skew boxes through denominator 128,000 and 76 determinant-one
Möbius charts found 40 points beyond the old uniform box; exact group-law
relations put all 40 in the certified rank-17 subgroup.  This is substantial
negative search evidence, not a rank upper bound.

The same certificate run directly replays all four exact conductors.  The
strict logarithmic inequality is also exact rather than a floating-point
decision: every conductor is less than `10^66`, while the degree-seven
positive Taylor partial sum proves `exp(2.31)>10`; hence
`ln N < 66*2.31 = 152.46 < 182.72`.  An independent exact Magma replay for
`u=42` also verifies the 17 points and their independence after transport to
the minimal model.  Its full 2-descent and rank-bound attempts exceeded the
public calculator's memory limit, and a capped 600-second PARI run returned no
upper bound.

**Verified computation and bounded experiments (Fermigier).**  Automatic
local-condition discovery scanned all 44 primes from 5 through 199, classified
188 compressed `p`-adic balls, and selected five clean split-multiplicative
groups.  It rediscovered the 7, 11, 17 and 19 conditions and replaced the
hand-selected 37-condition by the cheaper new condition
`T=+2,-2 mod 13`, which forces `v_13(H)>=3`.  The resulting 144-class search
found `T=154/103` with `ln N=162.234032455648...`; its height-`50000` point
search remained at numerical rank 12.

A newer construction starts from the benchmark's accidental quartic points
rather than from local scores.  Its height-one-million record replay exposes
fourteen search-relative accidental abscissas and gives an exact rank-22
finite-reduction certificate after adjoining the missing published direction.
All 28 lines `x=+/-T+n` through those abscissas are genus-one quartics.  Their
height-200,000 searches produced fourteen new parameters and five completed
conductors below the strict target; the best rank signal was `T=3115/3`, with
`ln N=133.171856293608...` and stable numerical rank 15 through height one
million.  Thus the construction succeeds at moving far inside the conductor
boundary, but not yet at retaining rank 21.

The exact inverse map also recovered five published-point preimages absent
from that height-one-million source set (`P14`, `P15`, `P20`, `P21`, `P22`).
All ten additional `x=+/-T+n` slices were searched through height 200,000.
After 78-fiber decontamination they produced four new parameters; the two
completed conductor calls had `ln N=191.834...` and `185.836...`, while the
other two hit their declared one-shot caps.  None reached the specialized
point-search gate.

To force two published directions simultaneously, a pairwise quotient screen
then formed all 220 cross-label products of the 22 signed preimage slices.
Every height-50,000 search recovered the record fiber exactly.  Three other
parameters made a product square, but in each case both individual factors
were nonsquares, so none actually forced two quartic points.  This closes the
declared pairwise quotient tranche with zero new double-forced fiber; it does
not rule out rational points beyond that height.

A disjoint deep tranche then replayed all 23,769 members of the original
height-50,000 multiple-root population and selected 48 fibers using a
discovery prime band and two held-forward bands, without using the cumulative
`B=500` score or any point result.  Twenty-eight conductor calls completed.
Two new fibers met the strict conductor bound: `T=3206/265`, with
`ln N=168.031754726474...`, root number `-1`, and stable numerical rank 14
through quartic height one million; and `T=1925/157`, with
`ln N=174.420464807504...`, root number `+1`, and stable numerical rank 12 at
height 50,000.  This leakage-controlled tranche is bounded negative evidence,
not a rank bound.

An independent exhaustive height-box run enumerated every primitive member of
the original 144-class union through projective height 50,000, modulo the
`T -> -T` curve symmetry: 23,769 nonsingular specializations.  Leakage-free
scoring reduced these to `256 -> 32 -> 12` at numerical prime cutoffs
`200 -> 2000 -> 10000`.  All three leaders scored below `E22`; five of six
requested conductor calls completed, with best `ln N=192.051614237934...`.
Thus this finite box contains no completed conductor success and no rank claim.

The expanded Fermigier conductor pass completed 21 of the first 24 pinned CRT
candidates and found four below the strict threshold: `644/87`, `847/184`,
`70/223`, and `1057/218`.  Uniform height-`50000` triage of these and the score
leaders produced no gain beyond numerical rank 12; only `1666/9` retained its
previous numerical rank 16 after escalation.  This is useful negative evidence
against spending more descent time on that Fermigier pool.

Two later searches close much larger, disjoint pieces of the Fermigier space.
An exact global enumeration covered 60,815,684 primitive sign-quotiented
parameters with `0<=a<=100000` and `1<=b<=1000`.  Fresh discovery and held
prime bands selected a conductor-first tranche; 22 completed conductors were
below the target, but the complete point stages topped out at numerical rank
13.  Independently, all 6,160 signed weight-three combinations of Fermigier's
22 certified record points were converted into exact auxiliary directions.
The initial 399-direction pilot produced one second fiber, `T=29771/78`, with
`ln N=158.467530623289...` and numerical rank 12.  All remaining 5,761
directions then completed both height-50,000 slope searches—11,522 calls in
total—and returned only the record-fiber calibration.  This exhausts the
declared weight-three direction population exactly; it does not exhaust
higher-weight combinations or arbitrary auxiliary curves.

**Verified computation and bounded experiment (family generalization).**  The
arbitrary six-root Mestre implementation has now moved far beyond the original
maximum-root-14 survey.  Exact affine normalization through maximum root 50
enumerated 1,032,506 tuples and left 44 nonsingular nonreflection families.
The maximum-root-100 extension enumerated 36,475,792 normalized tuples, with
777 nonreflection obstruction solutions and 235 generically nonsingular
families; 191 were genuinely new beyond the maximum-root-50 prefix.

A rank-blind conductor-first screen of the new families produced an exact
rank-at-least-13 specialization for roots `(0,6,47,55,70,80)` at `T=8`:

```text
N = 581863561133867566935518764040599206
ln N = 82.351544058010...
```

The twelve displayed quartic points span only numerical rank 10.  The exact
points with quartic abscissas `75/2`, `175/37`, and `243/4` add three
directions, and finite reductions certify the resulting thirteen points
independent.  A second new family, roots `(0,4,30,31,39,46)`, has an exact
rank-at-least-12 specialization at `T=5`, `ln N=79.729318123910...`.  Its
exhaustive 30,000-by-1,000 rational box contained 18,244,819 primitive
parameters; the strongest new signal was only numerical rank 13 at
`T=151/40`.  A separate 144-class prime-power CRT search produced the
subthreshold fiber `T=209001/3868`, `ln N=179.289493580807...`, but its point
span stayed at numerical rank 10 through height one million.  These results
show that new low-conductor families can be manufactured and certified, while
also confirming that visible point count and discriminant smoothness alone do
not predict the exceptional rank jumps needed for the target.

**Source-recovery boundary.**  The published Elkies rank-18 construction
cannot yet be instantiated from its official arXiv bundles because the model,
base change, and section coordinates are omitted.  The separate
[`ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md`](ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md)
records that archive search and the strongest explicit alternative found.
For Kumar--Kuwata's geometric rank-18 `F^(6)` example, an exact finite-fibre
Galois/lattice audit recovers fixed rank 5 over `Q(t)` and maximum quadratic
twist rank 3.  It is therefore not an arithmetic rank-18 specialization
family; the audit is computational evidence at four exact fibres, not a
formal symbolic function-field proof.

**Verified computation (Nagao calibration).**  The exact replay implements
Nagao's root tuples, the quadratic rank-13 base change, and the printed
rank-13 and rank-21 points.  It checks 13 and 21 printed points respectively,
matches both printed minimal models, and reproduces
`ln N=165.406045732331...` for the rank-13 `u=1` specialization and
`ln N=196.679545735892...` for the rank-21 record.  Published independence is
cited; the repository does not replace it with a new exact proof.

Leakage-free integer scans covered every positive integer `u<=200` and then
`u<=2000`.  The larger run used cutoffs `200 -> 2000 -> 20000`; only `u=1256`
and the already-known `u=42` among its completed finalists fell below the
conductor threshold.  A separate local-CRT search over 2,048 symbol choices at
primes 7, 11, 13, 19 and 31 found `u=118` above and also `u=316`, with
`ln N=177.107241730402...` but no point gain beyond the 13-section baseline.
Escalating `u=118` to height `10^6` found 43 nonvisible Jacobian images, all
checked exactly, but the height matrix remained numerical rank 15 at both 72
and 120 digits.  Point count is not rank gain.

Finally, a 110-member neighborhood of Nagao's rank-21 record produced the
conductor-only candidate `T=6041/198`.  It did not beat the published record's
long-cutoff score; a bounded quartic search found only the six visible
sign-pairs, whose numerical span has rank 11.  This is a negative local probe,
not evidence that the candidate's Mordell--Weil rank is 11.

### First end-to-end fixture

**Verified computation, not a target hit.**  In the Elkies--Klagsbrun family
at `u=2/5`, the constraints

```text
t = 5 mod 11^2,   t = 109 mod 17^2,   t = 102 mod 19^2
```

combine to `r=7814669 mod 12623809`; exact Gauss reduction recovers
`t=-1468/21`.  At each forced prime the specialization has
`v_p(Delta)=4`, type `I_4` nonsplit multiplicative reduction, and conductor
exponent one.  PARI/GP reports

```text
ln N = 148.626486493304...
ellrank bounds = [10,10]
```

Nine published sections were supplied and checked on the curve.  This fixture
validates the local-to-global mechanism and lands below the conductor
threshold, but its rank is far too small.  The full machine record is
[`elliptic_ek_k3_crt_fixture.json`](../artifacts/generated-results/elliptic_ek_k3_crt_fixture.json).

### Algebraic six-root family design

The Kihara replay has now been separated into its essential square identities
and its high-height parameterization.  For a mean-centered center polynomial

```text
B(z)=z^6+c4*z^4+c3*z^3+c2*z^2+c1*z+c0,
```

the degree-five obstruction in the paired-root Mestre remainder is exactly
`2*c1=c3*c4`.  On the general affine-normalized Mestre moduli, the smallest
nonvisible section ansatz is

```text
x=x0+x1*T,  y=z0+z1*T+z2*T^2+z3*T^3.
```

Exact coefficient elimination gives the Mestre equation, a leading-square
condition `D=w^2`, and three residual equations after eliminating all four
ordinate coefficients.  The executable equations and a verified
two-dimensional Fermigier component are documented in
[`notes/MESTRE_AFFINE_SECTION_MODULI.md`](notes/MESTRE_AFFINE_SECTION_MODULI.md).

A targeted search on this section locus produced the explicit centers
`(0,23,93,128,133,175)`.  Their primitive quartic has fixed square content
`2400^2`, twelve paired-root sections, and six verified nonvisible companion
identities.  Its leading coefficient is `9*(T^2+14406)`, so
`T=(14406-u^2)/(2u)` splits quartic infinity.  At `u=1`, exact good-reduction
quotients modulo 3 give dimensions 11 for the twelve visible images, 12 after
one chosen companion, and 13 after split infinity.  This proves generic rank
at least 13 by specialization.

The base-changed discriminant frontier is primitive, irreducible, squarefree,
and degree 40, versus degree 398 for Kihara's rank-at-least-14 template.  See
[`notes/MESTRE_RANK13_02393128133175.md`](notes/MESTRE_RANK13_02393128133175.md).
This is a new proved family with materially leaner conductor geometry.  It
does not prove generic rank equality, independence of all six companions, or
a high-rank specialization in that particular family.

A second exact D-square family has centers `(0,25,95,143,168,205)`, leading
coefficient `T^2+39146`, and the split-infinity base change
`T=(39146-u^2)/(2u)`.  Six nonvisible companions are classified exactly; at
`u=197` one raises the visible-plus-infinity mod-3 dimension from 12 to 13,
proving the generic lower bound.  The same specialization has 17 rigorously
independent rational points and exact conductor
`2462086522751621334987931952469307556796057284118717977320345864383117775914`.
The tracked certificate is
[`elliptic_mestre_dsquare_four_u197_rank17.json`](../artifacts/generated-results/elliptic_mestre_dsquare_four_u197_rank17.json),
and the complete bounded four-family screen is documented in
[`notes/MESTRE_DSQUARE_FOUR_SCREEN.md`](notes/MESTRE_DSQUARE_FOUR_SCREEN.md).

## Pilot family

[`cas/fermigier_mestre.py`](cas/fermigier_mestre.py) implements the normalized
root tuple

```text
(0, 55, 314, 378, 1007, 1036).
```

If `q(X)` has these six roots, the Mestre identity

```text
q(X-T) q(X+T) = g(X,T)^2 - r(X,T)
```

produces a genus-one quartic with twelve visible points and an additional
section in this member.  The code uses the binary-quartic Jacobian

```text
y^2 = x^3 - 27 I(T) x - 27 J(T).
```

Its internal normalization has `T = 2t`, so Fermigier's printed specialization
`t = 19754/39` is `T = 39508/39` here.  This is pinned by reducing the
specialization to the published `E22` model, rather than by comparing parameter
names alone.

The classical binary-quartic covariants also map all thirteen visible quartic
points exactly onto this Jacobian.  This gives the search a concrete seed set
for descent and point work instead of relying on generic-rank metadata alone.
The implementation checks the source equation, target equation, and covariant
syzygy with exact rational arithmetic for every mapped point.  The formulas are
placed in the context of the primary invariant-theory references in
[SOURCES.md](SOURCES.md).

[`cas/search_crt_lattice.py`](cas/search_crt_lattice.py) then:

1. finds simple roots of the degree-20 discriminant factor modulo selected
   primes and lifts them to prime powers;
2. keeps good-reduction rank-score primes disjoint from those forced bad
   primes;
3. combines one local choice per prime with CRT;
4. applies exact two-dimensional Gauss reduction to
   `L(r,M) = {(a,b): a-rb = 0 mod M}`;
5. scores the resulting `T=a/b`; and
6. asks PARI/GP for a minimal model, conductor, discriminant and root number on
   a bounded subset.

The stronger [`cas/multiple_root_lifting.py`](cas/multiple_root_lifting.py)
handles derivative-zero roots by enumerating every valid next `p`-adic digit.
Complete sibling sets are compressed into maximal residue balls, exposing the
actual congruence cost.  For example, the clean 7-adic condition used above
forces `v_7(H)>=18` using only six classes modulo 49, rather than one class
modulo `7^18`.

See [THEORY.md](THEORY.md) for the exact local lemmas and failure modes, and
[REPRODUCE.md](REPRODUCE.md) for commands.

## What is and is not certified

The exact Python layer verifies rational identities, modular roots, Hensel
lifts, CRT congruences, and forced valuations.  PARI/GP computes minimal models
and conductors.  When requested, PARI's `ellrank` returns unconditional lower
and upper bounds from 2-descent, but those outputs remain recorded here as
software computations.

The following do **not** certify a target curve:

- generic rank of a family;
- a Mestre--Nagao or local likelihood score;
- a root number or parity prediction;
- a numerical analytic-rank estimate;
- a bounded point search that fails to find another point; or
- a nonzero floating-point height determinant without precision and
  independence checks appropriate to the claim.

A target candidate must retain the explicit model, the minimal-model conductor
calculation, all claimed rational points, curve-membership checks, and a
reproducible independence certificate.  An exact upper bound is unnecessary
for either target: the required claim is a rank lower bound.

## Layout

- [`cas/fermigier_mestre.py`](cas/fermigier_mestre.py): exact normalized
  Mestre--Fermigier family and degree-20 discriminant factor.
- [`cas/ek_k3.py`](cas/ek_k3.py): exact Elkies--Klagsbrun rank-nine K3 family,
  local point counts, and linear-factor lifts.
- [`cas/elkies_klagsbrun_rank29.py`](cas/elkies_klagsbrun_rank29.py): public
  rank-29 model and points, exact short-model transport, and invariants.
- [`cas/verify_elkies_klagsbrun_rank29.py`](cas/verify_elkies_klagsbrun_rank29.py):
  portable exact finite-reduction rank-at-least-29 certificate.
- [`cas/search_elkies_klagsbrun_rank30.py`](cas/search_elkies_klagsbrun_rank30.py):
  bounded direct x-, secant-, and cross-ratio-chart point search.
- [`cas/search_elkies_klagsbrun_rank30_alternate_covers.py`](cas/search_elkies_klagsbrun_rank30_alternate_covers.py):
  weight-two/three subgroup-sum alternate-cover search.
- [`cas/search_elkies_klagsbrun_rank30_denominator_sieve.py`](cas/search_elkies_klagsbrun_rank30_denominator_sieve.py)
  and [`cas/search_elkies_klagsbrun_rank30_companion_center_sieve.py`](cas/search_elkies_klagsbrun_rank30_companion_center_sieve.py):
  disjoint exact modular-square sieves around public and subgroup-companion
  abscissas.
- [`cas/crt_lattice.py`](cas/crt_lattice.py): exact CRT and rank-two Gauss
  reduction.
- [`cas/search_crt_lattice.py`](cas/search_crt_lattice.py): bounded pilot
  search and JSON artifact writer.
- [`cas/multiple_root_lifting.py`](cas/multiple_root_lifting.py): complete
  digit lifting for simple and multiple roots, maximal-ball compression, and
  exact fixed-divisor valuations.
- [`cas/local_condition_discovery.py`](cas/local_condition_discovery.py):
  automatic discovery, compression, local classification, and efficiency
  ranking of Fermigier `p`-adic balls.
- [`cas/search_discovered_local_conditions.py`](cas/search_discovered_local_conditions.py):
  bounded CRT search driven entirely by the automatically selected local
  groups.
- [`cas/search_multiple_root_crt.py`](cas/search_multiple_root_crt.py): all 144
  declared multiple-root CRT symbol combinations and bounded lattice
  neighborhoods.
- [`cas/exhaustive_multiple_root_height.py`](cas/exhaustive_multiple_root_height.py):
  exact projective-height enumeration of the full 144-class union and
  leakage-free staged scoring.
- [`cas/screen_multiple_root_frontier.py`](cas/screen_multiple_root_frontier.py):
  expanded conductor and root-number pass over the pinned Fermigier CRT pool.
- [`cas/search_record_residue_class.py`](cas/search_record_residue_class.py):
  exhaustive projective-height scan inside the residue class containing
  Fermigier's record specialization.
- [`cas/compare_score_cutoffs.py`](cas/compare_score_cutoffs.py): pinned
  comparison of `T=1666/9`, `E22`, and the record scan's top ten through
  numerical cutoffs up to 100,000.
- [`cas/staged_record_rescore.py`](cas/staged_record_rescore.py): leakage-free
  staged rescoring of the full declared record-residue height box.
- [`cas/search_fermigier_global.py`](cas/search_fermigier_global.py): exact
  60.8-million-parameter global-box enumeration and held-band point triage.
- [`cas/search_fermigier_rank22_record_group_triples_remainder.py`](cas/search_fermigier_rank22_record_group_triples_remainder.py):
  terminal exact search of the record subgroup's weight-three auxiliary
  directions not covered by the pilot.
- [`cas/search_power_pairs.py`](cas/search_power_pairs.py): bounded control
  search over pairs of ordinary simple roots modulo `p^2`.
- [`cas/search_extra_points.py`](cas/search_extra_points.py): bounded quartic
  point enumeration, exact covariant mapping, and two-precision numerical
  height-matrix experiment.
- [`cas/batch_rank_triage.py`](cas/batch_rank_triage.py): uniform staged point
  search and numerical-height triage of pinned Fermigier candidates.
- [`cas/mestre_root_tuples.py`](cas/mestre_root_tuples.py): exact arbitrary
  six-root Mestre construction, quartic obstruction, visible points, and
  discriminant proxies.
- [`cas/survey_mestre_root_tuples.py`](cas/survey_mestre_root_tuples.py):
  bounded affine-normalized root-tuple survey and PARI specialization probes.
- [`cas/search_mestre_root_tuple_scale.py`](cas/search_mestre_root_tuple_scale.py)
  and [`cas/search_mestre_root_tuple_scale_max100.py`](cas/search_mestre_root_tuple_scale_max100.py):
  exact root-tuple censuses through maximum roots 50 and 100.
- [`cas/certify_mestre_0647557080_t8_rank13.py`](cas/certify_mestre_0647557080_t8_rank13.py):
  exact rank-at-least-13 certificate for the strongest new tuple-family
  specialization.
- [`cas/verify_mestre_02595143168205_rank13_section.py`](cas/verify_mestre_02595143168205_rank13_section.py)
  and [`cas/verify_mestre_dsquare_four_u197.py`](cas/verify_mestre_dsquare_four_u197.py):
  the generic-rank-at-least-13 companion replay and the standalone pinned
  rank-at-least-17 specialization certificate for the second D-square family.
- [`cas/search_mestre_0430313946_frontier.py`](cas/search_mestre_0430313946_frontier.py)
  and [`cas/search_mestre_0430313946_power_crt.py`](cas/search_mestre_0430313946_power_crt.py):
  global-box and prime-power searches in the new exact rank-12 family.
- [`cas/nagao_1994.py`](cas/nagao_1994.py): exact Nagao rank-13 base change,
  rank-21 family, published models, and printed point data.
- [`cas/nagao_linear_sections.py`](cas/nagao_linear_sections.py): all six
  linear companion sections and their exact generic Mordell--Weil relations.
- [`cas/verify_nagao_linear_sections.py`](cas/verify_nagao_linear_sections.py):
  exact elimination certificate classifying the 18 linear-abscissa sections.
- [`cas/mod2_reduction_independence.py`](cas/mod2_reduction_independence.py):
  reusable exact independence certificates from finite reduction quotients.
- [`cas/certify_nagao_rank17_frontier.py`](cas/certify_nagao_rank17_frontier.py):
  unconditional 17-point certificates for the four current frontier curves.
- [`cas/verify_nagao_1994.py`](cas/verify_nagao_1994.py): model, point,
  conductor, local-reduction, and discriminant-geometry replay.
- [`cas/search_nagao_rank13_integer_u.py`](cas/search_nagao_rank13_integer_u.py):
  leakage-free integer-`u` scoring and conductor search.
- [`cas/nagao_rank13_local.py`](cas/nagao_rank13_local.py): exact local
  discriminant balls and clean multiplicative classification in the Nagao
  base-change parameter.
- [`cas/search_nagao_rank13_local_crt.py`](cas/search_nagao_rank13_local_crt.py):
  Nagao local-condition/CRT/Gauss search.
- [`cas/triage_nagao_rank13_finalists.py`](cas/triage_nagao_rank13_finalists.py):
  exact section checks, bounded point searches, numerical-height replay,
  bounded `ellrank`, and small-prime subgroup saturation probes.
- [`cas/triage_nagao_rank13_local_candidates.py`](cas/triage_nagao_rank13_local_candidates.py):
  the same exact/numerical point triage for local-CRT survivors.
- [`cas/search_nagao_rank13_rank_gain.py`](cas/search_nagao_rank13_rank_gain.py):
  staged rare-event search whose selection objective is bounded extra-point
  yield and stable rank gain rather than a Nagao prime score.
- [`cas/extend_nagao_u42_frontier.py`](cas/extend_nagao_u42_frontier.py):
  capped height-`10^7`, exact saturated-basis, and effort-zero descent
  checkpoint for `u=42`.
- [`cas/extend_nagao_u118_height.py`](cas/extend_nagao_u118_height.py): capped
  height-`10^6` point and numerical-height checkpoint for `u=118`.
- [`cas/search_nagao_u42_skew_height.py`](cas/search_nagao_u42_skew_height.py):
  disjoint skew-height boxes, Möbius charts, and exact relation replay beyond
  the uniform `u=42` search.
- [`cas/search_nagao_u135_skew_height.py`](cas/search_nagao_u135_skew_height.py):
  the same bounded geometry at the smallest-conductor certified rank-17
  specialization.
- [`cas/nagao_skew_height.py`](cas/nagao_skew_height.py): parameter-independent
  target/checkpoint API for applying that search geometry to other certified
  candidates.
- [`tools/generate_u42_magma.py`](tools/generate_u42_magma.py) and
  [`tools/run_u42_pari_rank.py`](tools/run_u42_pari_rank.py): exact-input
  generators and a process-capped descent/rank-toolchain checkpoint.
- [`cas/search_nagao_rank21_neighborhood.py`](cas/search_nagao_rank21_neighborhood.py):
  bounded search in the record's forced local class.
- [`cas/triage_nagao_rank21_neighbor.py`](cas/triage_nagao_rank21_neighbor.py):
  exact visible-seed and bounded point-pool replay for the best low-conductor
  rank-21-family neighbor.
- [`cas/pari_bridge.py`](cas/pari_bridge.py): narrow PARI/GP subprocess bridge.
- [`cas/verify_fermigier_benchmark.py`](cas/verify_fermigier_benchmark.py):
  family normalization, conductor, and historical score replay.
- [`cas/verify_ek_k3_fixture.py`](cas/verify_ek_k3_fixture.py): end-to-end local
  power/CRT/Gauss/conductor/rank regression fixture.
- [`cas/verify_fermigier_rank22_points.py`](cas/verify_fermigier_rank22_points.py):
  exact replay of the 22 published points, exact transport to the normalized
  short model, and a finite-reduction independence certificate; the
  two-precision height matrix is retained only as a secondary replay.
- [`cas/verify_fermigier_1666_9.py`](cas/verify_fermigier_1666_9.py): exact
  specialization, point, conductor-factorization, and bad-prime replay for the
  best conductor-only candidate.

The corresponding new machine records are:

- [`elliptic_fermigier_discovered_local_conditions.json`](../artifacts/generated-results/elliptic_fermigier_discovered_local_conditions.json),
  [`elliptic_fermigier_multiple_root_height_h50000.json`](../artifacts/generated-results/elliptic_fermigier_multiple_root_height_h50000.json),
  [`elliptic_fermigier_multiple_root_frontier.json`](../artifacts/generated-results/elliptic_fermigier_multiple_root_frontier.json), and
  [`elliptic_fermigier_batch_rank_triage.json`](../artifacts/generated-results/elliptic_fermigier_batch_rank_triage.json);
- [`elliptic_mestre_root_tuple_survey.json`](../artifacts/generated-results/elliptic_mestre_root_tuple_survey.json);
- [`elliptic_nagao_1994.json`](../artifacts/generated-results/elliptic_nagao_1994.json),
  [`elliptic_nagao_rank13_integer_u.json`](../artifacts/generated-results/elliptic_nagao_rank13_integer_u.json), and
  [`elliptic_nagao_rank13_integer_u2000.json`](../artifacts/generated-results/elliptic_nagao_rank13_integer_u2000.json);
- [`elliptic_nagao_rank13_local_crt.json`](../artifacts/generated-results/elliptic_nagao_rank13_local_crt.json),
  [`elliptic_nagao_rank13_local_candidate_triage.json`](../artifacts/generated-results/elliptic_nagao_rank13_local_candidate_triage.json), and
  [`elliptic_nagao_rank13_finalist_triage.json`](../artifacts/generated-results/elliptic_nagao_rank13_finalist_triage.json);
- [`elliptic_nagao_u42_height_10000000.json`](../artifacts/generated-results/elliptic_nagao_u42_height_10000000.json) and
  [`elliptic_nagao_u118_height_1000000.json`](../artifacts/generated-results/elliptic_nagao_u118_height_1000000.json); and
- [`elliptic_nagao_linear_sections.json`](../artifacts/generated-results/elliptic_nagao_linear_sections.json),
  [`elliptic_nagao_rank13_rank_gain_search.json`](../artifacts/generated-results/elliptic_nagao_rank13_rank_gain_search.json),
  [`elliptic_nagao_rank13_rank_gain_mutations.json`](../artifacts/generated-results/elliptic_nagao_rank13_rank_gain_mutations.json),
  [`elliptic_nagao_rank17_frontier_certificate.json`](../artifacts/generated-results/elliptic_nagao_rank17_frontier_certificate.json), and
  [`elliptic_nagao_u42_skew_height.json`](../artifacts/generated-results/elliptic_nagao_u42_skew_height.json), and
  [`elliptic_nagao_u135_skew_height.json`](../artifacts/generated-results/elliptic_nagao_u135_skew_height.json); and
- [`elliptic_nagao_u42_magma_probe.json`](../artifacts/generated-results/elliptic_nagao_u42_magma_probe.json) and
  [`elliptic_nagao_u42_descent_toolchain.json`](../artifacts/generated-results/elliptic_nagao_u42_descent_toolchain.json); and
- [`elliptic_nagao_rank21_neighborhood.json`](../artifacts/generated-results/elliptic_nagao_rank21_neighborhood.json) and
  [`elliptic_nagao_rank21_neighbor_triage.json`](../artifacts/generated-results/elliptic_nagao_rank21_neighbor_triage.json).

Documentation and test support:

- [`tests/`](tests/): exact arithmetic and search-hazard regressions.
- [`THEORY.md`](THEORY.md): mathematical design and limitations.
- [`REPRODUCE.md`](REPRODUCE.md): executable commands and expected outputs.
- [`SOURCES.md`](SOURCES.md): primary sources and the status of each external
  claim.

Generated runs belong in `../artifacts/generated-results/`.  Do not silently
replace a pinned run: record the full command, parameters, software versions,
and changed hash.

## Structural-search groundwork

The next search layer is now implemented as reusable exact infrastructure, not
just a list of suggestions. See
[`notes/STRUCTURAL_SEARCH_GROUNDWORK.md`](notes/STRUCTURAL_SEARCH_GROUNDWORK.md)
and the pinned
[`elliptic_structural_search_groundwork.json`](../artifacts/generated-results/elliptic_structural_search_groundwork.json).

The new layer provides:

- finite-quotient **escape scoring** over arbitrary prime moduli, so a point
  that adds a certified direction is promoted before another dependent point;
- the exact generalized-Weierstrass 2-division cubic and residual-Selmer
  dimension bookkeeping for a relative rank-20 descent;
- bounded integral Néron--Severi lattice enumeration for low-degree K3
  multisections and alternate isotropic fibrations;
- exact branch/genus bookkeeping for the genus-9 Fermigier `V4` pair covers,
  exposing their genus-`2,2,5` quotient-Jacobian work queue;
- quadratic-twist character and projective `p`-adic chart primitives, including
  exact mixed affine/infinity congruence lattices; and
- strict CAS task schemas for Selmer/class groups, moduli components, isogeny
  hopping and Frobenius/Picard filters.

The real rank-20 mod-5 certificate is the first quotient calibration: its
20-column matrix has rank 12 on the first twelve columns and all eight
remaining columns individually escape that prefix image. This is an exact
filter calibration, not a new rank result.

## Near-term research queue

These are **open problems**, in descending order of leverage:

1. Run a **relative 2-Selmer/cubic-class-group computation** on the canonical
   rank-20 Fermigier--Mestre fiber. Quotient the Selmer data by the twenty
   pinned Kummer images and minimize only the residual covers. This can either
   close the fixed fiber at rank 20 or expose the precise classes in which a
   twenty-first generator could live.
2. Export and saturate the known rank-17 Néron--Severi sublattice of the
   Fermigier K3, then enumerate and geometrically realize low-degree, low-genus
   multisections and alternate elliptic fibrations. Do not restrict the search
   to simple rational formulas for `x(T)`.
3. Choose one of the 3,160 exact Fermigier `V4` pair covers and replace another
   projective-height extension by quotient-Jacobian arithmetic: the two
   genus-2 factors and genus-5 product quotient, followed by Selmer, Chabauty
   or a Mordell--Weil sieve.
4. Turn quadratic base changes into a **twist-section engine**. Search
   low-branch squareclasses `d(T)`, reject twists by an early Picard/Frobenius
   bound, and require a specialization that escapes the old subgroup. Making
   an old section divisible is saturation, not rank gain.
5. Replace larger integer root-tuple boxes by algebraic component searches in
   Mestre moduli: combine the quartic obstruction, companion-section and
   split-infinity equations, decompose over several finite fields, then lift
   and parameterize surviving components over `Q`.
6. Run complete rational isogeny graphs for every rank-18/19/20 finalist and
   move point/descent work to the easiest model. Rank and conductor remain
   invariant, but cover and coordinate complexity need not.
7. Extend local-condition discovery from affine `T` to all of
   `P^1(Q_p)`, especially denominator-divisible neighborhoods of multiplicative
   infinity fibers. Re-minimize every resulting specialization.
8. Keep quotient escape and finite-reduction certification in the inner loop.
   Point count, score and numerical height rank remain secondary filters;
   statistical rare-event models stay gated until the certified positive label
   set grows.

## Complementary certified pipeline

This is a research track separate from the Keller-map programme.  Its two
operational targets are:

1. exhibit an elliptic curve over \(\mathbb Q\) with **at least 21 independent
   rational points** and
   \(\log N<182.72\), where `log` is the natural logarithm; or
2. exhibit one with **at least 30 independent rational points**, improving the
   public rank-at-least-29 record.

The literal first cutoff is

```text
N <= 22609332114411420624526008180120083289443726642551045721326266745468787166869234.
```

The number `182.72` comes from the rounded value for Fermigier's rank-at-least-22
curve in Bober's table.  Its exact conductor has
`log(N) = 182.7249109506374287961...`, so that benchmark narrowly misses the
literal cutoff.  See [baselines and literature](notes/BASELINES_AND_LITERATURE.md).

### What is implemented

The first end-to-end calibration implements the proposed

```text
p-adic roots -> CRT -> exact 2D Gauss reduction -> small rational t
             -> integral curve -> minimal model/conductor/rank checks
```

for

\[
E_t:y^2=x^3-t^2x+t^2,
\qquad
\Delta(a,b)=-16a^4b^6(27b^2-4a^2),
\]

where \(t=a/b\).  Lifting roots at \(23^3,47^2,73^2\), considering all eight
root-sign combinations, and exhaustively searching the relevant rational-height
boxes gives

\[
t=-\frac{110627}{84367},\qquad
27b^2-4a^2=23^3\,47^2\,73^2.
\]

The resulting integral model is

\[
y^2=x^3-87109893594476435881x
 +620029989546545117143687312009.
\]

PARI/GP 2.15.4 verifies that the displayed model is already global minimal,
that the three shaped primes have conductor exponent one, that torsion is
trivial, and that `ellrank` returns bounds `[3,3]`.  The pinned manifest is
[here](../artifacts/generated-results/elliptic-curves/crt_lattice_calibration_v1.json).
This is an exact mechanism calibration, not a candidate for either target and
not an independent rank computation.

### First high-family seed

The programme now has an exact adapter for Fermigier's fixed-root Mestre
family.  Exact binary-quartic invariants connect its quartic construction to a
canonical model

```text
[1, a2(u), 1, a4(u), a6(u)],     u = s/2,
```

whose discriminant is one primitive even polynomial \(\Phi(u)\) of degree 20.
The apparent \(u^{12}\) and twelfth-power constant in the raw quartic
conversion are removable coordinate artifacts and are not used as shaping
factors.  Roots of \(\Phi\) at 89, 131, and 137 are simple and give split
multiplicative reduction.

Forcing all three primes to the second power, checking all eight root choices,
and exhaustively reconstructing rational parameters through the first occupied
height box gives

\[
u=\frac{673709}{29965},\qquad M=89^2 131^2 137^2.
\]

The homogeneous degree-20 discriminant factor has valuation exactly two at
each shaped prime.  PARI verifies split type \(I_2\), Tamagawa number 2, and
conductor exponent one at all three.  The
[pinned seed](../artifacts/generated-results/elliptic-curves/fermigier_crt_seed_v1.json)
does **not** claim a global conductor or a rank: its 469-bit remaining cofactor
is deliberately left unfactored.  Separately, the reconstructed thirteenth
quartic point and an exact finite-reduction certificate prove the twelve
generic section differences independent.  A second certificate proves all 22
published E22 points independent.  The
[pinned rank replay](../artifacts/generated-results/elliptic-curves/fermigier_rank_certificates_v1.json)
is a lower-bound certificate only: it supplies neither saturation nor an
unconditional rank upper bound.

There is also a reproducibility discrepancy in the source.  Fermigier prints
the symmetric shift `s=19754/39`, but that literal substitution gives a
different curve.  The displayed E22 curve is recovered exactly at
`s=39508/39`, or `u=19754/39`; the two literal specializations have different
\(j\)-invariants.  No erratum or derivation of this factor two was found, so the
coordinates remain explicitly separated rather than silently identified.

### Best low-conductor near miss

A bounded `ratpoints` search at

\[
u=\frac{28917}{20}
\]

found a point cloud from which exact finite-reduction arithmetic selects 20
independent points.  PARI gives the global minimal model

\[
y^2+xy+y=x^3+x^2
-4437412060110743641525245114305x
+3586842216822165612930264910099076801587288127
\]

with exact conductor

```text
2876153493562761211278364526603564191699143885403233935132057708367930
```

and `log(N) = 159.9348252255254533984...`, comfortably below the target
cutoff.  The [pinned near-miss certificate](../artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json)
records the 58 searched abscissas, exact search bounds, 20 reduction-log rows,
and minimal-model data.  This is one independent point short of the target;
it is not a solution and no rank upper bound is claimed.

### Exact high-rank baselines

Two further replays keep the search calibrated against stronger public
families and records:

- Kihara's fully published family is reconstructed at (t=2).  Fourteen
  rational-function sections specialize to fourteen points with an exact
  full-rank certificate modulo 5.  This proves generic rank at least 14, but
  it supplies no rank-30 candidate and its simple specialization is not
  conductor competitive.
- The 29 displayed points on the Elkies--Klagsbrun record curve are checked
  exactly.  A full-rank matrix modulo 2 and an independent torsion witness
  prove rank at least 29.  This local replay does not claim exact rank 29 and
  does not produce a thirtieth point.

Their equations, sources, and certificate commands are recorded in
[baselines and literature](notes/BASELINES_AND_LITERATURE.md).

### Search drivers

`ecsearch/fermigier_score_sweep.cpp` is a deterministic staged rank-jump
ranking pass.  At every odd scoring prime it computes the exact mean trace on
the good projective fibres of the Fermigier family and scores a specialization
by $A_p-a_p(E_t)$, not by its raw trace.  It retains cumulative and disjoint
windowed residual $S_0$ and $S_5$ features through 2000; a bad
specialization prime is excluded from ordinary-trace scoring.  Its output is
only a heuristic ordering: every survivor still requires a bounded point
search, an exact independence certificate, global minimization, and an exact
conductor.  The current bounded sweep produced the rank-20 near miss above;
it did not produce a target curve.

### Important separation of local roles

If a prime \(p\) is forced into \(\Delta\), the specialized curve has bad
reduction at \(p\).  It therefore cannot simultaneously contribute an ordinary
good-reduction trace \(a_p\) to a Nagao score.  The search uses two roles:

- **shaping primes:** force prime powers in discriminant factors and score
  reduction type, local Euler factor, Tamagawa data, or root number;
- **scoring primes:** require good reduction and score finite-field traces.

The local constraints can still be assembled globally by CRT.  The final
specialization is always minimized before its conductor is accepted.

### Programme map

- [Mathematical pipeline](notes/CRT_LATTICE_PIPELINE.md) gives the exact
  congruence, lattice, primitive-pair, and conductor checks.
- [Baselines and literature](notes/BASELINES_AND_LITERATURE.md) records what is
  external and conditional, together with the exact Kihara rank-14 and public
  record replays.
- [ICARM curve 273 rank-30 certificate](notes/ICARM_CURVE273_RANK30.md)
  records the exact lower-bound proof, independent Sage replay, and claim
  boundary for the 2026 public record.
- [Fermigier reproduction](notes/FERMIGIER_REPRODUCTION.md) records the exact
  canonical-family bridge and unresolved factor-two source discrepancy.
- [`families/`](families/) contains the rank-at-least-two calibration family
  and the exact Fermigier--Mestre rank-at-least-twelve adapter.
- [`ecsearch/`](ecsearch/) contains dependency-free exact arithmetic.
- [`scripts/`](scripts/) contains generation and replay commands.
- [Reproduction commands](REPRODUCE.md) are the public entry points.

The rank-at-least-30 target is complete. The immediate arithmetic priorities
are a certified residual 2-Selmer calculation for curve 273, which may prove
exact rank 30 or expose a 31st direction, and a twenty-first point on the
rank-20 low-conductor near miss. In parallel, the 17-by-17 lattice programme
is reconstructing a lower-rank neighbor fibration from which the missing
rank-17 K3 model may be recovered.

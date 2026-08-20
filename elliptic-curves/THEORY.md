# Mathematical design and limits

This note records the mathematical core of the pilot.  It deliberately
separates exact implications from search heuristics.

## 1. Why prime powers can be conductor-friendly

For a minimal elliptic curve over `Q`, a prime of multiplicative reduction has
conductor exponent one, while its minimal discriminant valuation can be much
larger.  Thus forcing

```text
v_p(Delta) >> 1
```

may make the discriminant large without paying the same exponent in the
conductor.  The useful event is not merely `p^k | Delta`: one must also verify
minimality and the local reduction type.  Additive reduction, denominators, and
the primes 2 and 3 require separate treatment.

This mechanism occurs conspicuously in the public rank-29 curve: its
discriminant has many repeated small-prime factors, whereas every bad prime
appears in the conductor with exponent one.  That is empirical motivation for
the strategy, not a theorem that powerful discriminants cause high rank.

## 2. A clean local theorem in the K3 test family

The Elkies--Klagsbrun test family is

```text
E: y^2 = x^3 + 2 A x^2 + B x,
```

with `B` a product of eight linear factors in the specialization parameter.
Its invariants are

```text
c4    = 16(4 A^2 - 3 B),
c6    = 64 A(9 B - 8 A^2),
Delta = 64 B^2(A^2 - B).
```

**Theorem (elementary local calculation).**  Let `p` be odd and suppose a
specialization is `p`-integral, `v_p(A)=0`, and `v_p(B)=n>0`.  Then the displayed
model is minimal at `p`,

```text
v_p(Delta) = 2n,    v_p(c4) = 0,
```

so the reduction is multiplicative of Kodaira type `I_(2n)` and the conductor
exponent is one.  It is split multiplicative exactly when `2A` is a square
modulo `p`.

The implementation only labels a root “clean” when exactly one `B` factor
vanishes modulo `p` and `A` does not.  This avoids collisions of factors and
additive fibers.  See [`cas/ek_k3.py`](cas/ek_k3.py).

For the Fermigier pilot, a simple root of its degree-20 discriminant factor is
retained only when the binary-quartic invariant `I` is nonzero modulo `p`.
For `p >= 5`, this is the corresponding `c4 != 0` multiplicative case.  The
actual minimal model and conductor are still recomputed after specialization.

## 3. Simple and multiple-root lifting

Let `F(T)` be an integral discriminant factor.  If

```text
F(r_1) = 0 mod p,       F'(r_1) != 0 mod p,
```

Hensel's lemma gives a unique compatible root `r_k mod p^k` for every `k`.
The original pilot implements exact Newton lifting and verifies the final
congruence.

The discriminant polynomial has many more useful derivative-zero roots.
[`cas/multiple_root_lifting.py`](cas/multiple_root_lifting.py) lifts these by
digits.  If `r` is a root modulo `p^j`, every child has the form

```text
r + d*p^j,    0 <= d < p,
```

and the exact first-order congruence is

```text
F(r+d*p^j)/p^j = F(r)/p^j + d*F'(r) mod p.
```

When `F'(r)=0 mod p`, either all `p` children survive or none do.  The routine
enumerates every surviving child and raises an exception if its safety cap is
exceeded; it never returns a silently truncated root set.  Complete sibling
sets are then compressed into maximal `p`-adic balls.

There is a second source of cheap divisibility.  For an integral polynomial of
degree `d`, its fixed divisor is the gcd of `F(0),...,F(d)`, by the Newton
expansion in integer-valued binomial polynomials.  Consequently the exact
minimum of these valuations proves divisibility on an entire affine residue
class, even when coefficient content does not expose it.

**Verified computation for the Fermigier factor `H`.**

- `H(T)` is always divisible by `5^2` for integral `T`.
- Modulo 7, `T=0` is a root of multiplicity 16; if `T=7s`, then
  `v_7(H(T))>=18`.  Clean multiplicative reduction after one 7-adic
  minimalizing scaling occurs for `s!=0 mod 7`, giving six usable classes
  modulo 49 and minimal discriminant valuation at least 6.
- `T=0,5,6 mod 11` forces `v_11(H)>=4` with clean split multiplicative
  reduction.
- `T=+1,-1 mod 17` and `T=+5,-5 mod 19` each force valuation at least 4.
- `T=+4,-4 mod 37` forces valuation at least 3.

Thus the five default local groups have `6*3*2*2*2=144` combinations modulo
`49*11*17*19*37=6441589`.  This is the central multiple-root insight: large
discriminant powers can cost only a low-power residue condition.

The record-containing class makes the same structure visible with fewer
primes: its residues at `7,17,37` have multiplicities `16,4,3` modulo those
primes and force valuations at least `18,4,3`, respectively.

The first conditions were found by hand.  The generalized discovery layer now
enumerates every root ball through a declared depth, compresses full sibling
sets, evaluates the minimal-model scaling and split symbol on the coarse ball,
and ranks clean groups by radical saving divided by congruence cost.  In the
bounded scan of all 44 primes from 5 through 199 it classified 188 balls and
found 16 eligible clean split groups.  Its automatic five-prime objective kept
7, 11, 17 and 19 and selected

```text
T = +2,-2 mod 13  =>  v_13(H) >= 3
```

instead of the hand-chosen 37-condition.  The classification and the final
PARI local-reduction replay are separate checks: automatic discovery does not
infer the global conductor from a polynomial valuation.

**Bounded exhaustion.**  For the original five-group union, exact congruence
enumeration through projective height 50,000 contains 23,769 primitive,
nonsingular curves after quotienting the even-family sign symmetry.  A
leakage-free `B=200 -> 2000 -> 10000` score cascade retained
`23769 -> 256 -> 32 -> 12`.  Every stage leader scored below the separately
evaluated `E22` benchmark.  Five of six requested conductors completed; the
best had `ln N=192.051614...`.  This proves only a negative statement about
that explicit finite population and those completed conductor calls.

### Rational specialization

For a rational parameter `T=a/b` with `p` not dividing `b`, the homogeneous
condition is

```text
a - r_k b = 0 mod p^k.
```

It follows that the homogenization of `F` at `(a,b)` is divisible by `p^k`.
The code then evaluates `F(a/b)` exactly and records the actual valuation; it
does not rely on this implication alone.

For pairwise distinct primes, CRT combines selected roots into

```text
a - r b = 0 mod M,      M = product(p_i^k_i).
```

For a compressed ball, the corresponding lower ball modulus replaces the
naïve `p^k`.  Exact evaluation of `H(a/b)` and PARI local reduction still occur
after reconstruction; fixed-divisor guarantees do not replace minimal-model
checks.

## 4. Exact two-dimensional lattice reduction

The rational representatives satisfying the global congruence form the
rank-two lattice

```text
L(r,M) = {(a,b) in Z^2 : a-rb = 0 mod M}
```

with basis `(M,0)`, `(r,1)` and determinant `M`.  In dimension two, exact
Gauss reduction is sufficient; floating-point LLL is unnecessary.  The pilot
reduces the basis using integer dot products and enumerates bounded linear
combinations, retaining primitive vectors with `gcd(b,M)=1`.

The reduction step is exact, but the bounded coefficient-radius enumeration is
**not** a proof that every useful representative has been found.  A search run
must record that radius.

## 5. Bad-prime and good-prime information are different

A prime deliberately forced into the discriminant is a bad-reduction prime.
It cannot simultaneously supply a good-reduction trace observation for a
Nagao-style score.  At a multiplicative prime the Euler-factor coefficient is
`+1` or `-1`; that is useful local information, but it is not the same random
variable as the good-reduction trace `a_p`.

The pilot therefore uses disjoint sets:

- **power primes** impose lifted discriminant roots and record split versus
  nonsplit multiplicative reduction;
- **rank primes** admit only good residues and contribute a local score; and
- the full candidate score skips bad primes and primes dividing the rational
  denominator.

This separation matters for any later Bayesian model.  Mixing the two types as
though they had one likelihood would double-count deliberately imposed data.

## 6. The current score is only a ranking statistic

The `fermigier-good` option in
[`cas/search_crt_lattice.py`](cas/search_crt_lattice.py) adds, over good primes
in the requested numerical range,

```text
(2-a_p)/(p+1-a_p) * log(p).
```

The search objective during CRT assembly is

```text
sum(local scores) - lambda * log(max(|a|,b)).
```

Neither expression is a probability or a rank bound.

The name `fermigier-good` is intentional.  It denotes this repository's
good-reduction-only, numerical-prime-bound variant and is not Fermigier's
historical table statistic.

**Verified computation (historical compatibility).**  Fermigier's paper
prints a specialization score with `p <= M`, but its `E22` table

```text
M       50     100    200    400    1000   2000   4000   10000
S(E22)  29.49  44.12  57.54  81.51  105.17 122.76 143.84 166.47
```

is reproduced by treating `M` as a prime ordinal: omit `p=2`, then sum from
the second prime through the `M`-th prime.  It is not reproduced by taking a
numerical cutoff `p <= M`.  The command is in [REPRODUCE.md](REPRODUCE.md).
Accordingly, the pilot's `--score-bound 200` means the literal prime bound 200
and must not be compared with the historical table's column `M=200`.

**Bounded experiment (cutoff stability).**  The height-5,000 scan's top ten
were selected at numerical cutoff 500, so their values there are
selection-biased.  When the selected leader `T=1666/9` is evaluated without
further selection at cutoffs 500, 2,000, 10,000 and 100,000, its scores are
`40.048807...`, `60.519323...`, `81.634378...` and `122.707444...`; published
`E22` scores `40.913185...`, `69.525178...`, `106.746181...` and
`163.165765...`.  The early near-tie does not persist.

A separate leakage-free staged run starts with all 5,520 declared parameters
at `B=2000`, retains 50, then retains 10 at `B=10000` and the same 10 at
`B=100000`.  Later cutoffs see only the immediately preceding survivors.
`T=1666/9` is not a finalist; the final leader `T=1547/492` scores
`124.536543...`, versus `163.165765...` for the separately evaluated `E22`
benchmark.  No conductor or rank was computed in this experiment.  Staging
prevents retrospective cutoff leakage, but pruning is still a bounded search
choice and the scores still imply no rank.

The more principled future statistic is a calibrated tail likelihood such as

```text
log P(local data | rank >= 21) - log P(local data | baseline rank),
```

evaluated out of sample.  The Elkies--Klagsbrun search already discusses a
Bayesian version of its score and staged cutoff searches, so the novelty here
would be calibration against certified tail labels and joint optimization with
conductor constraints, not merely the words “Bayesian” or “staged”.

## 7. Exact binary-quartic map and rank seeds

For the homogenized binary quartic `U`, the implementation forms normalized
covariants

```text
g = (U_XY^2 - U_XX*U_YY)/144,
h = (U_X*g_Y - U_Y*g_X)/8
```

and maps an affine point `(x,z)` on `z^2=U(x,1)`, with `z!=0`, to

```text
(X,Y) = (36*g(x,1)/z^2, 108*h(x,1)/z^3)
```

on `Y^2=X^3-27IX-27J`.  The normalization is checked directly: Python exact
rational arithmetic verifies both curve equations and the resulting covariant
syzygy for all thirteen visible points at test specializations.  Primary
references for the invariant-theory construction and covering map are listed
in [SOURCES.md](SOURCES.md).

Twelve images are normally passed to PARI as a rank-search seed because the
generic construction supplies rank at least 12.  Exact membership is a point
certificate; it does not prove their independence after specialization.

**Bounded experiment.**  At `T=1666/9`, a quartic search through naive height
`10^6` found 57 distinct rational `x`-values, including 44 beyond the visible
sections.  Every image was checked exactly on the Jacobian.  At both 96 and
192 digits the height matrix has numerical rank 16; the recorded 16-point
subset has determinant about `3.479008e23` and minimum eigenvalue
`0.999542...`.  These stable floating-point data remain below the target and
are not an exact independence certificate.  A recorded `ellrank` attempt with
a 1 GB stack ran about 300 seconds and overflowed before producing bounds.

## 8. Two search hazards

### Singular exact roots

If a discriminant factor is linear, its lifted residue can rationally
reconstruct to the exact root.  That vector is extremely short but produces a
singular specialization.  For example, at `u=2/5` the fifth primitive `B`
factor is

```text
67t + 28.
```

Lifting its root modulo `11^4` makes `t=-28/67` an exact lattice vector.  A
search that accepts only the shortest vector can therefore select precisely
the value it must reject.  Every candidate is checked for nonzero
discriminant, and several short vectors should be enumerated.

### Height is not monotone under incremental CRT

The height of the shortest rational representative can fall after another
congruence is added.  Therefore pruning a partial CRT state solely because its
current representative is large is unsafe.

**Verified counterexample.**  For residues at primes `19,23,29,37`, a greedy
width-one path encountered height `1409`, while the complete choice

```text
r = 238875 mod 468901
```

has the representative `48/53`, of height 53.  Indeed

```text
48 - 238875*53 = 0 mod 468901,
```

and its residues are respectively `7,20,2,3`.  Beam width is therefore a
heuristic resource bound, not a correctness parameter.  Meet-in-the-middle,
larger retained frontiers, or delayed height pruning are the next remedies.

## 9. Six-root constructions are a search space, not a rank proxy

[`cas/mestre_root_tuples.py`](cas/mestre_root_tuples.py) makes the six-root
input arbitrary.  For a tuple `A=(a_1,...,a_6)`, it verifies the polynomial
remainder condition that makes

```text
q(X-T)q(X+T) = g(X,T)^2 - r(X,T)
```

have `deg_X(r)<=4`, removes only exact square content, and constructs the
binary-quartic Jacobian and twelve displayed points.  These identities are
exact.  They do not say that the displayed points are distinct modulo sign or
independent.

Reflection symmetry is a particularly important failure mode.  If the six
roots pair around a common center, six displayed points can be the negatives
of the other six on the Jacobian.  Such a family may have an attractive
discriminant radical but fails the intended visible-rank gate.

**Bounded survey.**  Among the 1,023 affine-normalized integer tuples with
largest root at most 14, 68 pass the quartic condition and 59 are generically
nonsingular.  Of those, 57 are reflection-symmetric.  The two remaining root
tuples are

```text
(0,2,8,9,11,14),    (0,1,7,8,9,11).
```

One effort-zero PARI probe for each returned computational rank bounds `[9,9]`
and `[4,4]`; the other 57 probes also had lower bound at most 9.  This is not a
classification beyond the declared height box, but it decisively rejects
“twelve displayed abscissas plus low discriminant degree” as a sufficient
family score.  Future family comparison should measure verified independent
sections and extra-point yield per logarithmic conductor growth.

## 10. Nagao's rank-13 base change and rank-21 neighborhood

Nagao starts with the root tuple

```text
(148,116,104,57,25,0)
```

and applies the quadratic base change

```text
T = (23550-u^2)/(2u).
```

The twelve Mestre sections and the additional base-change section are
implemented with exact rational arithmetic.  Nagao's generic-rank-13 result is
a cited theorem; each specialization still gets exact point-membership and
numerical-height checks because a generic theorem does not by itself store a
specialized independence certificate.  The primitive discriminant factor has
degree 20 in `T` and the cleared base-changed factor has degree 40 in `u`; PARI
finds both relevant factors irreducible over `Q` in the recorded replay.

The local search works directly in homogeneous `u=a/b`.  A retained ball must
force a discriminant valuation, keep the short-model `c4` a unit, keep `u` a
unit, and have a constant split-multiplicative tangent symbol.  Under those
hypotheses the same elementary local implication applies: the model is
multiplicative at the odd prime and the conductor exponent is one.  The final
specialization is nevertheless minimized and replayed prime by prime.

**Exact finite-reduction independence lemma.**  Let
`P_1,...,P_r in E(Q)`, and for finitely many good primes reduce their images
in the exact binary quotients `E(F_p)/2E(F_p)`.  If the matrix obtained by
stacking those quotient coordinates has column rank `r` over `F_2`, every
integral relation among the `P_i` has all coefficients even.  After division
by two, the resulting half-relation is a rational 2-torsion point.  If
`E(Q)[2]=0`, it is zero; iteration makes every original coefficient divisible
by every power of two, hence zero.  This proves independence using only exact
finite group arithmetic and a trivial-rational-2-torsion check.  It does not
use heights, BSD, parity, a full 2-descent, or the finite-index premise of
`ellsaturation`.

**Verified rank-17 frontier.**  Applying that lemma gives four unconditional
rank lower bounds below the strict conductor threshold:

```text
u=135/2, T=5065/36: ln N=144.927455914576..., rank >= 17;
u=471/11, T=5579/22: ln N=146.678928806750..., rank >= 17;
u=42,    T=3631/14: ln N=148.621053634068..., rank >= 17;
u=74,    T=9037/74: ln N=151.423206831026..., rank >= 17.
```

For each curve, the certificate checks all 17 rational points, exact reduction
maps and full `F_2` column rank, and trivial rational 2-torsion.  It also
reconstructs the rational short model and directly replays the minimal model,
conductor and root number instead of trusting the discovery artifact.  The
strict logarithmic comparison is exact rather than a floating-point decision:
every replayed `N` is below `10^66`, while the positive degree-seven Taylor
partial sum gives
`exp(231/100) > 80381233705038797/8000000000000000 > 10`.  Hence
`ln(N) < 66(231/100)=7623/50=152.46 < 182.72`.  For `u=42` the combined matrix
is `18 x 17` of rank 17, and the reduced 2-division cubic at `p=31` has no
root.  These are exact, portable rank-lower-bound certificates, but none
reaches rank 21: four further independent points are still required.

**Exact generic companion audit.**  Solving symbolically for all polynomial
sections with `x=mT+n` and `deg(y)<=3`, up to changing the sign of `y`, finds
exactly

```text
(m,n)=(+/-1/15,703/15), (+/-7/15,928/15), (+/-5/3,3628/15).
```

The `+1/15` section is Nagao's thirteenth section.  Exact group-law identities
over `Q(T)` express each of the other five sections in the pinned generic
basis; the certificate records the relation vectors.  They must therefore be
excluded from specialization-only point yield.  This changes point labels,
not the previous numerical height ranks.  The classification is deliberately
narrow: it does not classify arbitrary rational sections.

**Exact rank-29 record-fiber replay and bounded rank-30 attack.**  For the
public Elkies--Klagsbrun curve

```text
y^2+x*y=x^3-27006183241630922218434652145297453784768054621836357954737385*x
 +55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497,
```

all 29 published point pairs are checked exactly.  The change of variables

```text
X=36*x+3,  Y=108*(2*y+x)
```

is an isomorphism over `Q` to an integral short model.  Reductions at

```text
19,23,29,37,47,53,59,73,79,83,97,101,103,107,109,
127,131,151,157,163,173,179
```

give a stacked `F_2` matrix of column rank 29.  The short 2-division cubic has
no root modulo 67, so `E(Q)[2]=0`; the finite-reduction lemma above proves the
29 points independent unconditionally.  The published claim that the exact
rank is 29 assumes GRH and is not used here.

A bounded search on the fiber completed two 1,647-chart manifests.  The
deeper manifest searched 406 affine abscissa-pair charts, 812 oriented
secant-offset charts and 400 selected slope cross-ratio charts to naive chart
height 50,000, plus 29 abscissa-offset charts to height ten million.  It found
no affine point outside the chart seeds.  A materially different pass formed
every 4,060 positive weight-two/three subset sum `Q` in the certified
subgroup and used the degree-two coordinate
`t_Q(P)=(Y(P)+Y(Q))/(X(P)-X(Q))`.  It retained the 64 covers with smallest
public `t_Q` heights and completed 448 offset/cross-ratio charts, again with
no nonpublic image.  A third, disjoint tranche used signed representatives of
higher-weight subset sums.  It performed 2,000 signed-representative
evaluations across five weight bands, retained ten covers of weights 4
through 27, and completed 50 affine,
cross-ratio, three-point Möbius and skew charts.  The 57 exact returned images
all replay as public seeds or known companions in the rank-29 subgroup.  A
strictly capped 120-second PARI effort-zero 2-descent
accepted all 29 points but returned no rank interval.  None of these negative
computations is a rank upper bound, and no 30th point was certified.

**Exact rank-at-least-30 public-record replay.** The later 2026 public ICARM
curve 273 is a different curve. All 30 displayed rational points are checked
exactly and transported by

```text
X=36*x+3,  Y=108*(2*y+x)
```

to an integral short model. Exact exhaustive reductions at 25 good primes
give a stacked `31 x 30` matrix in products of `E(F_p)/2E(F_p)` with binary
rank 30. The short 2-division cubic has no root modulo 23, hence the curve has
no rational 2-torsion. The same infinite-descent lemma therefore proves the
30 points independent and

```text
rank E(Q) >= 30
```

unconditionally. A second Sage implementation reconstructs the finite groups
through invariant factors and discrete logarithms and obtains the same rank.
PARI separately checks the global minimal model, exact conductor, trivial
torsion, and root number `+1`. The complete proof record and claim boundary are
in [`notes/ICARM_CURVE273_RANK30.md`](notes/ICARM_CURVE273_RANK30.md). No
unconditional upper bound or rank-at-least-31 claim follows.

The source-recovery audit also prevents a misleading generic-rank shortcut.
Kuwata's explicit rank-18 K3 examples have geometric Mordell--Weil rank over
an algebraic closure, not eighteen sections over `Q(t)`.  For the exact
Kumar--Kuwata `F^(6)` basis, finite-fibre Galois matching at `t=2,3,5,7`,
followed by exact group-relation and height-isometry checks, gives fixed rank
5 and quadratic-character ranks `3,2,2,1,1,0,0`.  Its maximum quadratic twist
rank is therefore 3 in this computational audit.  Because the section
identities were not proved symbolically over the full function field, this is
recorded as an exact finite-fibre audit rather than a theorem.

**Bounded rare-event searches.**  Leakage-free scans of every positive
integer `u<=200` and `u<=2000` produced many conductor-only successes but no
target.  Earlier uniform exact/numerical point triage found

```text
u=42:  ln N=148.621053..., stable numerical rank 17 at height 10^6;
u=84:  ln N=139.773456..., stable numerical rank 16 at height 10^6;
u=50:  ln N= 89.115263..., PARI effort-zero bounds [13,13].
```

The last line is a stored PARI software computation and eliminates `u=50` as
a target; it is not being promoted to a portable descent proof.

A subsequent rank-gain search optimized unexpected exact quartic point yield
before stable numerical height rank and only then computed conductors.  The
primary declared population had 9,196 rational parameters; a second declared
mutation population had 5,133.  After excluding the five generic companions,
the maximum stable numerical rank was 17 and no rank-18 candidate appeared.
At height `10^6`, the primary run found respectively 66, 33, and 32 unexpected
abscissas for `u=135/2`, `42`, and `471/11`; the mutation run added `u=74`
with 20.  The search artifacts themselves make no exact rank claim; the four
rank lower bounds above come from the separate finite-reduction certificates.

For `u=42`, ten disjoint skew boxes with denominators through 128,000 and 76
determinant-one Möbius charts found 40 abscissas outside the previous uniform
height-`10^6` box.  Every returned point lies exactly on the quartic and
Jacobian, and every one has an exactly replayed relation in the certified
rank-17 subgroup.  The augmented numerical height rank remains 17 at 72 and
120 digits.  This is useful bounded negative evidence, not an upper bound.

The same search geometry on the smallest-conductor frontier curve, `u=135/2`,
found 74 exact points outside height `10^6` in the ten skew boxes, reaching
absolute numerator 43,277,563 and denominator 107,700.  All have exact
relations in its certified rank-17 subgroup, with relation coefficients of
absolute value at most 8.  The 76 Möbius charts added no outside point, and the
augmented numerical height rank again stayed 17.  This remains bounded
negative evidence, but it changes the immediate search priority: explore
alternate 2-covers before merely deepening this same quartic cover.

A later rare-event model deliberately used only four finite-reduction-certified
positive fibers, three-fold cross-fitted controls, and a published rank-21
fiber held out from training.  Both its Dirichlet-smoothed local Naive Bayes
score and a separately declared paper-informed comparator failed the
predeclared leave-one-positive-out recovery thresholds.  The gate therefore
authorized no broad population or conductor search.  With so few certified
tail labels, a tuned classifier would currently encode selection leakage more
readily than exceptional-rank information.

Small-prime `ellsaturation` had separately returned 17 checked points with
height-determinant ratio `2^32`, consistent with index `2^16` between the
recorded subgroups under its finite-index premise.  That premise is no longer
used for the rank-17 claim.  A basis-assisted effort-zero PARI `ellrank` run
timed out after 60 seconds, and the quartic height-`10^7` search timed out after
120 seconds.

**Bounded descent diagnostics.**  No rigorous upper bound for `u=42` was
obtained.  PARI 2.17.4 emitted no `ellrank` interval before a strict 600-second
timeout under an 8 GB stack/RSS cap.  The pinned eclib/mwrank source failed on
a machine-integer range in 2-descent; a temporary arbitrary-integer diagnostic
removed that cast failure but exposed an infeasible enumeration and was not
retained.  The Sage launcher present on the machine was unusable, and no
installation or system mutation was attempted.  Magma V2.29-9 independently
verified the 17 points and their independence, but its anonymous public
calculator reached its observed 311.34 MB memory limit on `RankBounds`,
`TwoSelmerGroup`, and `TwoDescent` paths.  Those resource failures returned no
Selmer group, covers, or rank upper bound and are not mathematical evidence for
one.

**Bounded local CRT search.**  Clean groups at 7, 11, 13, 19 and 31 give 2,048
symbol choices, identified in 1,024 sign-paired curves.  The strongest
conductor/rank tradeoffs after the declared score stages are

```text
u=118: ln N=128.027255..., stable numerical rank 15 at height 10^6;
u=316: ln N=177.107241..., stable numerical rank 13 at height 50000.
```

The local valuations and reduction types are exact verified computations; the
finite candidate enumeration and height ranks remain bounded experiments.  At
`u=118`, the height-`10^6` run found 43 nonvisible images but no numerical-rank
gain at either 72 or 120 digits.

Nagao's second tuple `(399,380,352,47,4,0)` specializes to his published
rank-21 curve.  The exact replay checks all 21 printed points and gives
`ln N=196.679545...`.  Preserving its five conspicuous local conditions and
enumerating 110 nearby rational parameters found constructor
`T=6041/198` with `ln N=170.765123...`.  At height 50,000, however, only the
six visible Jacobian sign-pairs were found and their numerical span has rank
11.  This says nothing about points outside the bounded search; it does show
that inheriting the record's local factorization pattern does not inherit its
21 known points.

## 11. Nagao's section-7 K3 and the rank-20 frontier

Nagao's section-7 tuple

```text
(346,260,255,146,55,0)
```

defines an elliptic K3 surface in the repository's constructor parameter
`T=2t`.  The short Weierstrass coefficients have degrees 8 and 12.  Its affine
discriminant has degree 20 and is squarefree, giving twenty `I1` fibers; the
fiber at infinity is split `I4`.  Thus the trivial lattice has rank
`2+(4-1)=5`.

Exact symbolic elimination finds twelve independent generic sections and
classifies all polynomial quartic abscissas through degree five.  In addition
to the twelve visible Mestre formulas, there are six linear and three
quadratic companion formulas, but eight exact group-law identities put every
new companion in the same rank-12 group.  The degree-three, degree-four and
degree-five leading branches are eliminated by exact Groebner bases.  This is
a classification of that polynomial ansatz, not of arbitrary rational
sections.

There is also an unconditional generic-rank theorem.  At the good prime 29,
direct exact point counting gives

```text
#S(F_29)   = 1212
#S(F_29^2) = 723600,
```

so the traces on `H^2` are 370 and 16318.  The split `I4`, fiber, zero section,
and twelve independent sections give seventeen divisor classes defined over
`Q`, hence seventeen Frobenius eigenvalues `+29`.  Removing them leaves traces
`-123` and `2021`.  Reciprocity and the Weil bounds reconstruct the residual
factor uniquely as

```text
(X+29) * (X^4 + 94*X^3 + 3828*X^2 + 79054*X + 707281).
```

The quartic has no normalized cyclotomic factor, and the full residual factor
evaluates to `534704436` at `X=29`.  A divisor defined over `Q` reduces to a
class in the `+29` eigenspace, while smooth proper specialization injects the
characteristic-zero Neron--Severi group.  Consequently the seventeen explicit
rational classes exhaust `NS(S/Q)`.  Shioda--Tate therefore proves

```text
rank E(Q(T)) = 17 - 5 = 12.
```

Geometrically, the single residual eigenvalue `-29` is also `29` times a root
of unity, while the other four are not.  Hence `rho(S/Qbar)<=18` and

```text
12 <= rank E(Qbar(T)) <= 13.
```

This last interval does not assert that a thirteenth geometric section exists.
In particular, there is no missing thirteenth section over `Q(T)`; the search
must exploit specialization rank jumps.

The strongest such jump currently certified is Nagao's paper parameter
`t=5081/94`, constructor parameter `T=5081/47`.  Exact finite-reduction
matrices certify twenty independent rational points, and PARI gives

```text
ln N = 174.2498162285480383539...
root number = +1.
```

This misses the target by one independent point.  A complete scan of all
`2^20-1` nonzero classes of the certified subgroup followed by bounded
alternate-cover, uniform and skew-chart searches found 224 decontaminated
points, all exactly dependent in the same rank-20 subgroup.  Separately, a
sharpened explicit-formula computation gives a conservative bound below 22;
with root number `+1`, GRH implies analytic rank at most 20, and BSD+GRH would
make the algebraic rank exactly 20.  That conditional closure is useful for
search allocation but is not an unconditional upper bound.

The fixed fiber also suggests a more general search-space construction.  After
removing the twenty-one known generic quartic abscissas, its bounded point pool
contains sixteen accidental abscissas.  Through an accidental point
`(T0,x0)` put

```text
x = m*T + (x0-m*T0).
```

For the section-7 quartic, the leading coefficient after substitution is a
square multiple of `(m^2-1)^2`.  Hence `m=+1` and `m=-1` cancel the degree-six
and degree-five terms and produce genus-one quartics in `T`; the other fifteen
integer slopes `-8<=m<=8` produce square-free genus-two sextics.  This turns an
exceptional specialization point into a new auxiliary Diophantine curve whose
rational points deliberately force an extra point on later specializations.

The first exact bounded tranche classified all 272 slices.  Pinned
doubled-point ternary boxes on thirteen of the genus-one Jacobians produced no
new conductor-plausible parameter.  A separate search completed all 240
genus-two slices through height 5,000.  Only one slice produced non-generic
parameters: `T=163` and `T=1049/10`, with natural logs of their exact
conductors `104.922002135807...` and `115.105316156014...`.  Both specialized
height-50,000 screens stayed at
numerical rank 12.  The deeper height-one-million calls on that genus-two
slice timed out, so they are incomplete, not negative enumerations.  These
experiments establish the slice mechanism and two unusually small conductors,
but no rank-21 certificate.

## 12. Exact finite-reduction certificate for Fermigier's rank-22 curve

Fermigier printed twenty-two rational points on the minimal model

```text
y^2 + x*y + y = x^3
  - 940299517776391362903023121165864*x
  + 10707363070719743033425295515449274534651125011362.
```

The published regulator calculation proves that they are independent.  The
repository now has a second, wholly exact certificate.  The inverse of the
normalized-short-to-minimal point change `v=(u,r,s,t)` is

```text
x_short = u^2*x_min + r,
y_short = u^3*y_min + s*u^2*x_min + t,

v = (14/507, 49/771147, 7/507, 1372/130323843).
```

All twenty-two transported points satisfy the short equation exactly.  Their
images in the product of `E(F_p)/2E(F_p)` for

```text
p = 29,43,67,73,79,83,89,101,103,107,109,127,131,137,149,191,223
```

give a binary matrix of column rank twenty-two.  At the separate good prime
31 the reduced 2-division cubic has no root, so `E(Q)[2]=0`.  Any integral
relation among the points therefore has all coefficients even.  Dividing the
relation and iterating would otherwise produce rational 2-torsion, so infinite
descent forces every coefficient to vanish.  Thus

```text
rank E(Q) >= 22
```

unconditionally, without using numerical heights or GRH.  The exact conductor
is the benchmark value recorded above and
`ln N=182.724910950637...`; consequently this curve misses the strict
`ln N<182.72` requirement by about `0.00491095`.  The exact independence
certificate is a stronger replay of the benchmark, not a target hit.  The
strict comparison is exact as well: a rational exponential-series bound gives
`e<1359141/500000`, and direct integer arithmetic verifies
`1359141^4568 < 500000^4568*N^25`, hence `ln N>4568/25=182.72`.

The same accidental-slice construction is especially natural here because
the benchmark is so close to the conductor boundary.  An exact
`hyperellratpoints` replay through quartic height one million finds 27
signless abscissas and a stable numerical image rank 21: thirteen selected
generic abscissas and fourteen further search-relative abscissas.  One of the
latter is the `T -> -T` conjugate of Mestre's extra generic section, leaving
thirteen genuinely specialization-accidental sources.  All fourteen sources
nevertheless give nonsingular genus-one quartics from `x=+/-T+n`, because the
Fermigier leading term again contains `(m^2-1)^2`.

The first bounded height-200,000 slice pass completed all 28 curves.  It found
fourteen new parameters, including three lying on two slices.  Five completed
conductor calls were below the strict target:

```text
T=3115/3       ln N=133.171856293608...
T=11305/6      ln N=132.718858978907...
T=9191/30      ln N=145.544405660627...
T=121919/260   ln N=165.609479503410...
T=88893/26     ln N=176.118486741750...
```

Their exact forced pools contain at most fifteen points.  Direct specialized
quartic searches put the strongest member, `T=3115/3`, at stable numerical
rank 15 through height one million (134 signed points, 54 non-generic images);
the other completed screens were no stronger.  Thus the construction crosses
the conductor boundary cleanly but has not yet preserved enough of the record
fiber's exceptional rank.  These are bounded numerical triage results, not
rank upper bounds.

Five of the eleven accidental preimages of Fermigier's published points were
not among those height-one-million sources.  Their ten additional
`x=+/-T+n` quartics were therefore searched separately through height 200,000.
Exact decontamination left four new parameters.  The completed conductors at
`T=56441/810` and `T=56441/240` have natural logarithms about 191.834 and
185.836, both above the target; the other two one-shot conductor calls timed
out.  No parameter passed the specialized point-search gate.

For simultaneous persistence, all 220 unordered cross-label products of the
22 signed published-preimage slice polynomials were searched through height
50,000.  The product-square condition is necessary for both factors to be
squares.  Every pair recovered the record parameter as an exact positive
calibration.  Three further product-square parameters appeared, but direct
exact tests showed both individual factors nonsquare in every case.  Hence
this tranche produced no new fiber carrying two forced published directions.

A separate leakage-controlled tranche regenerated all 23,769 nonsingular
members of the height-50,000 multiple-root population.  It excluded 111
previously examined parameters and selected 48 fibers from a discovery band
`p<=199` and held-forward bands `211<=p<=349` and `353<=p<=499`; neither the
cumulative `B=500` score nor point data entered selection.  Of 28 completed
conductor calls, only `T=3206/265` and `T=1925/157` were below the strict
target.  Their stable numerical ranks were respectively 14 through height
one million and 12 at height 50,000.  The other twenty conductor attempts hit
their declared one-shot caps.  This is a bounded out-of-sample negative result
and does not exclude higher rank elsewhere in the population.

## 13. Expanded root-tuple space and exact bounded frontiers

The six-root construction should be treated as a moduli search, not as a
single named family.  After translating and scaling a primitive integer tuple,
fixing its first two roots, quotienting reflection, imposing the Mestre
quartic obstruction, and rejecting a generically singular discriminant gives
a finite exact census at any maximum root.  The maximum-root-100 census has

```text
36,475,792 affine-normalized tuples,
33,945 obstruction solutions,
777 nonreflection solutions,
235 generically nonsingular nonreflection families.
```

Of these, 191 lie beyond the maximum-root-50 prefix.  These counts concern
families, not ranks.  In particular, the twelve displayed quartic points can
be dependent after specialization.

For roots

```text
(0,6,47,55,70,80)
```

at `T=8`, the displayed points have stable numerical height rank 10.  Exact
quartic points with abscissas

```text
75/2, 175/37, 243/4
```

raise the selected point span successively to 11, 12, and 13.  After exact
curve-membership checks, their reductions in products of
`E(F_p)/2E(F_p)` at

```text
p = 5,13,19,23,29,31,41,47,67,71,73
```

have binary column rank 13; reduction at `p=11` rules out rational
2-torsion.  The infinite-descent lemma of Section 10 therefore proves

```text
rank E(Q) >= 13
```

unconditionally.  PARI saturation supplied smaller exact coordinates but is
not a premise of the independence proof.  The conductor is

```text
581863561133867566935518764040599206,
```

with `ln N=82.351544058010...` and root number `-1`.

This expands the certified family frontier without approaching the requested
rank.  A second new family, roots `(0,4,30,31,39,46)`, gives an exact
rank-at-least-12 fiber at `T=5`, while its 18,244,819-parameter global box and
an independent 144-class prime-power CRT tranche never exceeded numerical
rank 13 and 10 respectively.  The lesson is structural: repeated prime powers
can make conductors exceptionally small, but they do not by themselves retain
exceptional Mordell--Weil directions.

The same distinction appears in the two record-fiber searches.  For the
Fermigier rank-22 benchmark, every one of the 6,160 signed weight-three
subgroup directions was converted into two exact auxiliary slices.  Across
the pilot and remainder, 12,320 height-50,000 calls were completed; only the
pilot fiber `T=29771/78` was new, and its stable numerical rank was 12.  For
the rank-29 record, two exact denominator-normalized sieves covered
55,267,250,510 primitive rational abscissas around public and subgroup-
companion centers.  All 217,309 exact square-test survivors were nonsquares.
Both are exact bounded exclusions of their declared search regions, not rank
upper bounds.

## 14. Certification ladder

Candidates should advance through these stages without promoting evidence:

1. exact family identity and nonsingularity;
2. forced valuation replay and local-reduction checks;
3. minimal integral model and exact conductor;
4. rational point discovery and exact curve-membership checks;
5. exact independence of at least 21 or 30 points, for example by a portable
   finite-reduction certificate or a rigorously documented descent; and
6. optional rank upper bound, clearly marked unconditional or conditional.

PARI's `ellrank` returns `[r1,r2,s,L]` with unconditional
`r1 <= rank(E(Q)) <= r2`; `L` is a list of independent non-torsion points.  It
can determine `r1=r2` without finding every point, so the number of returned
points must not be substituted for either bound.  Expensive or randomized
effort settings must be recorded.

# ICARM curve 245: exact rank-at-least-20 replay

## Status

This note records a **verified computation**, not an exact-rank theorem.  The
public ICARM curve 245 has twenty rational points that are proved
`Z`-independent here, so its Mordell--Weil rank over `Q` is at least 20.  No
twenty-first point and no rank upper bound is claimed.

The global minimal model is

```text
y^2 + x*y + y = x^3 - x^2
  - 25880411472355347134118026792*x
  + 1606663697747901005185875883284420820193259.
```

Its exact conductor is

```text
272066437942638823321634957004224153562929337633497250319389959310
```

and PARI records

```text
log(N) = 150.668907152237...
```

This is below the strict `182.72` target and improves the repository's
smallest-conductor exact rank-at-least-20 frontier.  It remains one independent
point short of the low-conductor target.

## Exact certificate

The checker first verifies all twenty coordinates in rational arithmetic.  It
then uses the integral short model obtained from

```text
X = 36*x - 9,
Y = 108*(2*y + x + 1).
```

The point images in the exact quotients `E(F_p)/2E(F_p)` at

```text
11, 23, 29, 41, 47, 59, 61, 67,
73, 83, 97, 101, 113, 127, 139, 149
```

give a stacked binary matrix of column rank 20.  The short 2-division cubic
has no root modulo 7, so `E(Q)[2]=0`.  Infinite descent therefore proves that
the twenty points are `Z`-independent.

The discriminant factorization is pinned completely:

```text
-Delta = 2^17 * 3^7 * 5^4 * 13^4 * 19^5 * 37^2
         * 7770053
         * 763973980372286963203
         * 55722582408764114465841769948159.
```

PARI local reduction at exactly these primes gives conductor exponents

```text
1, 2, 1, 1, 1, 1, 1, 1, 1,
```

whose product replays the displayed conductor.  Every local minimal change is
`[1,0,0,0]`, and PARI independently returns the displayed global minimal
model and root number `+1`.  The conductor inequality does not depend on a
floating logarithm: `N<10^66` and the rational exponential-series proof
`log(10)<231/100` give `log(N)<152.46<182.72`.

Run

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve245_rank20.py --check
```

against
[`icarm_curve245_rank20_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve245_rank20_v1.json).
The pinned artifact has SHA-256
`487d6e072ed7a2508d7ab12663910b3028c8b23362039c3e8b93a278809a2cbd`.

## Conditional fixed-fiber closure

A separate sinc-squared explicit-formula diagnostic uses `Delta=11/5` and
every prime through `1007525`.  Its conservative value is

```text
21.018943490740643741500498320491540766844285216155...
```

which is strictly below 22.  The exact root number is `+1`, so GRH would force
the analytic rank to be even and at most 20.  Together with the unconditional
twenty-point lower bound, GRH+BSD would make the algebraic rank exactly 20.
This is explicitly a conditional diagnostic, not an unconditional rank upper
bound.  It redirects the constructive search from this fixed fiber toward its
Mestre neighborhood or another surface.

Replay it with

```sh
PYTHONPATH=elliptic-curves/cas .venv/bin/python \
  elliptic-curves/cas/explicit_formula_icarm_curve245_delta22.py --check
```

against
[`icarm_curve245_explicit_formula_delta22_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve245_explicit_formula_delta22_v1.json).

## Bounded next-point searches

The initial exact search completed 990 x-coordinate and secant-slope charts
with no timeout and no nonbasis image.  A separate alternate degree-two-cover
pass completed 720 charts based on 120 short Mordell--Weil combinations.  It
recovered 199 distinct non-witness points, all of which replay exactly in the
certified rank-20 subgroup.  These are finite negative experiments, not a rank
upper bound or a saturation proof.

Public source: [ICARM curve 245](https://elliptic-rank.icarm.cloud/curve/245).

## Recovered Mestre parent

The terse public commentary can now be completed exactly.  Fermigier's
six-root formulas at

```text
(u,v) = (3/2,2)
```

give the labelled root set

```text
{-375/16,-269/16,-31/16,25/4,219/16,89/4}.
```

At the public parameter `T=5801/160`, its primitive quartic Jacobian has the
same exact `j`-invariant as curve 245.  Scaling roots and parameter by 16 and
translating the least root to zero gives the integral presentation

```text
roots = (0,106,344,475,594,731),   T = 5801/10.
```

This also explains why the earlier affine-normalized census through diameter
300 could not identify the parent: the primitive diameter is 731.  In the
integral presentation the primitive short Jacobian at the anchor is

```text
[0,0,0,
 -33541013268172529885816962722027/6250000,
  37480250741062883714416405828874482717136103027/7812500000].
```

The exact Weierstrass change

```text
(u,r,s,t) = (3/25,-9/2500,3/50,27/31250)
```

sends this model to the public minimal model.  Fermigier's additional generic
quartic section becomes

```text
x(T) = (4558+7*T)/29,
y(T) = 4801853 - 123478023*T/841
       + 88438*T^2/841 + 792*T^3/841.
```

All identities, the degree-eight/twelve short-Jacobian coefficient
polynomials, and the public-model change are replayed by
`test_icarm_curve245_mestre.py`.  This is an exact reconstruction of the
parent family, not a claim that nearby fibers retain rank 20.

## Conductor-first parent-family experiments

Three bounded searches use the recovered normalization.

1. The local interval `500 <= T <= 660`, denominators at most 320, contains
   `4,997,120` reduced non-anchor parameters.  A compiled exact-local sieve
   retained 4,096 before conductor or point data were read.  Of 118 completed
   conductor calls in the fixed 144-fiber panel, 50 were below `182.72`.
   Fixed quartic tiers at heights 5,000, 50,000 and 250,000 reached maximum
   numerical rank 17, at `T=7611/13` with
   `log(N)=154.482708820694...`.  No finite-reduction rank-18 trigger fired.

2. The complete canonical rectangle `a<=30000`, `b<=1000` evaluates
   `18,244,818` reduced non-anchor parameters.  Its corresponding fixed panel
   had 124 completed conductors, 60 below the target, and maximum numerical
   rank 16 in the same three point tiers.

3. An exact mod-3 replay of the anchor's height-50,000 quartic pool certifies
   rank 20 and separates eight accidental pivot points from the twelve
   visible sections and Fermigier's extra generic section.  The 16 lines
   `x=+/-T+n` through those pivots give genus-one slices.  Their bounded
   height-50,000 searches produced 236 non-anchor fibers.  In the fixed
   conductor panel, 113 calls completed and 83 were below the target, but the
   promoted point tiers reached only numerical rank 14.

The three artifacts are local experimental output:

```text
artifacts/local/elliptic-curves/icarm245-mestre-neighborhood-v1.json
artifacts/local/elliptic-curves/icarm245-mestre-global-a30000-b1000-v1.json
artifacts/local/elliptic-curves/icarm245-accidental-slices-v1.json
```

They contain no new rank-20 curve and no target hit.  In particular, the
numerical ranks are not algebraic-rank certificates.  The exact parent
reconstruction is the durable result; the bounded experiments justify moving
away from unconditioned one-parameter widening.

## Cross-shape comparison with curve 275

Public ICARM curve 275 is a second exact rank-at-least-20 anchor in the same
two-parameter Fermigier construction, but with a different six-root shape:

```text
(u,v)=(-3,-1/2),
native T=10239/176,
canonical roots=(0,113,550,753,868,1058),
canonical T=3413/11.
```

The exact change `(u,r,s,t)=(18/121,27/14641,9/121,2916/1771561)` sends its
primitive short Jacobian to the public minimal model.  Reconstructing 54
quartic images from the pinned height-2,000,000, denominator-13,000
`ratpoints` output and injecting Fermigier's extra generic point gives a
mod-3 finite-reduction certificate of rank 20.  The selected basis separates
as twelve generic and eight exceptional directions.  Replay it with
`verify_icarm_curve275_mestre_rank20.py`; its exact score is

```text
log(N)=162.6382563237891921510260149605792965012...
```

This second shape permits a transport test unavailable inside the fixed
curve-245 family.  For every one of the `8*8=64` pairs of exceptional
abscissas, interpolate `u,v,T,x` affinely between curves 245 and 275 and
evaluate the complete Fermigier quartic.  After removing all rational-square
factors, every residual squareclass is an irreducible degree-30 polynomial,
hence a genus-14 hyperelliptic condition.  The exact finite classification is
`analyze_icarm245_275_cross_shape_transport.sage`.  It rules out only this
smallest affine cross-shape ansatz, not nonlinear transports.

A complementary alternate-cover search then completed 790 charts on curve
275: all x-pair and slope-offset charts at height 10,000, the integer-x charts
at height 100,000, and the best 200 slope-pair charts at height 10,000.  It
found no nonbasis image.  The ignored local artifact is
`artifacts/local/elliptic-curves/icarm275-alternate-covers-h10000-v1.json`.
This is a bounded negative search, not a rank upper bound.

The same `Delta=11/5` conditional explicit-formula diagnostic used for curve
245 gives

```text
curve 262   21.0378455960178433187028971023...
curve 275   21.3153055869002603980901987375...
```

Both exact root numbers are `+1`.  Thus GRH forces analytic rank 20 for each
fixed curve, and GRH+BSD would force algebraic rank 20.  This is conditional
fixed-fibre closure, not an unconditional upper bound.  Replay the two exact
prime sums with `explicit_formula_icarm_curve262_275_delta22.py`.

The remaining three public rank-at-least-20 curves below the cutoff have the
same exact parity and conditional diagnostic:

```text
curve 243   21.3580557352708097572321413297...
curve 226   21.1934326698192049111537092186...
curve 7     21.3153720855424697102270624343...
```

Their exact root numbers are also `+1`.  Hence, under GRH, all seven public
sub-cutoff rank-20 fibres in the comparative screen have analytic rank 20;
under GRH+BSD they have algebraic rank 20.  Unconditionally each statement
remains only rank at least 20.  The three new prime sums replay with
`explicit_formula_icarm_curve7_226_243_delta22.py`.

## Unused-slope deformation experiment

The record E22 auxiliary searches use the two slopes `x=+/-T+b`; exactly
those slopes cancel the degree-six term and produce genus-one slices.  The
five smallest unused integral slopes `0,+/-2,+/-3` instead give 55 genus-two
slices through the eleven E22 accidental points.  Exact bounded
`hyperellratpoints` searches at height 200,000 found six genuinely
non-generic fibres, all on the `P13`, slope-3 slice, and no parameter collision
between distinct E22 accidental sources.

Two of the six fibres lie below the conductor target.  Their numerical height
triage is

```text
T=2429/6    log(N)=115.1074986574338...   stable numerical rank 14
T=25753/60  log(N)=162.0039130466744...   stable numerical rank 12
```

These ranks are numerical diagnostics, not algebraic-rank certificates.  At
`T=2429/6` the forced direction and the height-50,000 search contribute two
independent numerical directions beyond the generic rank 12.  Taking those
two exact quartic points as new sources, all four `+/-1` auxiliary slices were
enumerated through signed support four.  The resulting 16,824 exact parameter
images were pairwise disjoint across the two sources, so that iterative lane
was stopped rather than widened.  The reproducing scripts are
`search_fermigier_rank22_multislope_collisions.py` and
`search_fermigier_t2429_two_source_collisions.py`; their outputs remain local
bounded experiments.

## Comparative low-conductor screen

The same exact chart engine was applied to the other public rank-at-least-20
curves below the `182.72` cutoff.  The smallest exact scores encountered were

```text
curve 245   log(N) = 150.668907152237...   rank >= 20
curve 262   log(N) = 154.605048...          rank >= 20
curve 92    log(N) = 159.934825225525...   rank >= 20
curve 275   log(N) = 162.638256323789...   rank >= 20
curve 243   log(N) = 163.698567860195...   rank >= 20
curve 226   log(N) = 164.520926899998...   rank >= 20
curve 7     log(N) = 170.087664842249...   rank >= 20
```

No twenty-first independent point was found.  Curve 262's public basis was
not 2-saturated: an exact recovered point `Q` satisfies

```text
2*Q = -2*P1-P2-P3-P4+P6-2*P8-P9-P10+P12-P13
      +3*P14-P15+P16+3*P18+P20,
Q = (-1493957509135121/9, -9186675491109878322692/27).
```

Replacing `P2` by `Q` gives an exact rank-20 mod-2 certificate.  A reduced
basis then completed 720 alternate-cover charts, producing 668 distinct
nonbasis images; all replay in the same subgroup.  This diagnoses a basis
index, not a new rank direction.

Curve 243 has a separate exact six-root parent, with roots
`{-1851/4,-1455/4,-687/2,-1149/4,1437/2,2955/4}` and anchor
`T=3895/6`.  A bounded height-200,000 replay finds all twelve visible
abscissas and sixteen non-generic abscissas.  Searching all 32 associated
`x=+/-T+n` genus-one slices in the same box finds exactly two non-anchor
cross-source collisions:

```text
T=27265/144   log(N)=205.134568539366...   root number -1
T=15580/7     log(N)=206.134128500343...   root number +1
```

Each collision forces two distinct accidental quartic directions, but both
conductors exceed the target and no rank or independence claim is made.  The
bounded experiment is reproduced by
`search_icarm_curve243_accidental_slices.py` and stored under
`artifacts/local/elliptic-curves/`.

Curve 226 similarly decodes to the exact six-root parent
`{-138,-90,-60,-12,138,162}` at `T=10167/350`.  The height-200,000 anchor
search has twelve visible and nineteen non-generic abscissas.  Its 38
`x=+/-T+n` slices have eight cross-source collisions: the all-source values
`T=+/-69` are singular, while the six regular values are

```text
T=10167/2800   log(N)=201.906197622792...   root number -1
T=10167/1225   log(N)=169.164465132563...   root number -1
T=10167/550    log(N)=157.419681221708...   root number -1
T=111837/2450  log(N)=182.289068403305...   root number +1
T=10167/100    log(N)=143.707349732487...   root number +1
T=40668/175    log(N)=159.810097946694...   root number -1
```

After carrying the exact slice ordinates into each collision fibre, all six
height matrices have stable numerical rank 11 at 96 and 192 digits.  Thus
even the four sub-cutoff fibres do not give a rank-jump candidate in this
bounded lane.  The reproducing script is
`search_icarm_curve226_accidental_slices.py`; its output remains under
`artifacts/local/elliptic-curves/` because this is an experiment, not a
Mordell--Weil certificate.

For comparison, the six screened rank-at-least-19 curves all have exact root
number `-1`.  Curve 90's initial 24 unresolved images were likewise halves of
the published subgroup.  The exact identity stored in `icarm_curve90.py`
replaces one public generator by a half-point, and the resulting reduced basis
has an exact rank-19 mod-2 certificate.  An expanded 1,440-chart run completed
without timeout, found 840 distinct nonbasis images, and resolved all 840 in
that subgroup.  These finite searches justify changing routes; they do not
prove exact rank 19 or 20.

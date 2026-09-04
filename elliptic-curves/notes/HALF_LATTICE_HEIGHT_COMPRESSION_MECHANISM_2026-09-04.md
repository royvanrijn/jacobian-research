# Why the deep half-lattice charts work (2026-09-04)

## Result

The half-lattice search has a precise geometry-of-numbers interpretation.  It
does not depend on local solubility and it is not a disguised Selmer search.
For

\[
 E:y^2=x^3+Ax+B,
 \qquad Q=(x_Q,y_Q),
\]

the pointed chart

\[
 C_Q:\quad
 w^2=t^4-6x_Qt^2-8y_Qt-3x_Q^2-4A
\]

is a degree-two coordinate on `E` whose fiber involution is

\[
 R\longmapsto Q-R.
\]

Consequently the chart is canonically centered at the half-lattice point
`Q/2`.  A deep class creates a large empty ball around that center containing
no point of the subgroup already known.  A new point is compressed when it
falls unusually close to the center.  PARI minimization and reduction choose a
small rational coordinate for the same degree-two linear system; they change
the bounded height distortion but not the midpoint.

This explains the search mechanism.  It does **not** explain why a particular
specialization has new Mordell--Weil directions near those holes.  That is the
separate rank-jump arithmetic question.

## Exact identities

Put

\[
 t_Q(R)=\frac{y_R+y_Q}{x_R-x_Q}.
\]

The line of slope `t_Q(R)` passes through `-Q` and `R`.  Its third
intersection with `E` is `Q-R`.  Comparing the `x^2` coefficient after
substituting the line into the Weierstrass equation gives the exact identity

\[
 \boxed{t_Q(R)^2=x(R)+x(Q-R)+x(Q).}
\]

The sign change `w -> -w` exchanges the two points `R` and `Q-R` in the
fiber.  The parallelogram law for the Neron--Tate pairing then gives

\[
 \boxed{
 \widehat h(R)+\widehat h(Q-R)
 =\frac{\widehat h(Q)+\widehat h(2R-Q)}2.}
\]

In `E(Q) tensor R`, the target-relative part is

\[
 \boxed{
 2\widehat h(R-Q/2)=\frac12\widehat h(2R-Q).}
\]

Now let `M` be the known Mordell--Weil lattice, let `c` be a class in
`M/2M`, and choose a shortest representative `Q_c` of `c`.  Since
`2P-Q_c` runs through the same parity class as `P` runs through `M`, one has

\[
 \boxed{
 \min_{P\in M}\widehat h(P-Q_c/2)
 =\frac14\min_{Q\in c}\widehat h(Q)
 =\frac14\widehat h(Q_c).}
\]

Thus the legacy half-lattice depth is exactly the squared radius of the ball
about `Q_c/2` that is empty of known lattice points.  “Deep” means good
old-point exclusion, not Selmer depth and not a larger intrinsic probability
of rank gain.

Every horizontal Mobius change used by `hyperellminimalmodel` and
`hyperellred` preserves the degree-two divisor class.  Its pole fiber still
has two points whose sum is `Q`, so its canonical midpoint is still `Q/2`.
For the reduced coordinate `s` there is a chart-dependent constant `C_s` with

\[
 \left|h(s(R))-\frac12\widehat h(2R-Q)\right|\le C_s.
\]

This is the height-machine form of the compression statement.  The useful
next theorem is a *tight explicit* local decomposition of `C_s`, not another
mask census.

The rational cancellation is also completely explicit.  If the reduced-to-
raw horizontal map is

\[
 t=\frac{as+b}{cs+d}
\]

and the primitive raw coordinate is `t=[T:U]`, then

\[
 s=[dT-bU:-cT+aU].
\]

Therefore

\[
 H(s)=
 \frac{\max(|dT-bU|,|-cT+aU|)}
 {\gcd(dT-bU,-cT+aU)}.
\]

This gcd is the exact denominator/pole cancellation that can turn a huge raw
slope into a tiny reduced coordinate.

Finally, for a binary quartic
`a t^4+b t^3+c t^2+d t+e`, direct substitution in the classical formulas
gives, for every `Q`,

\[
 \boxed{I=-48A,\qquad J=-1728B.}
\]

The invariants are constant across the chart family and cannot rank parity
classes.  Integralization scales them predictably; the reduced generalized
quartics again have one invariant pair on each curve.

## Ledger audit

The reproducible audit covers 3,865 completed detailed chart records:

| ledger slice | charts | first-independent event charts | charts returning finite points |
|---|---:|---:|---:|
| rank-28 selected union | 64 | 11 | 43 |
| curve 385, first quotient iteration | 301 | 9 | 118 |
| curve 398, deepest MW16 classes | 12 | 5 | 8 |
| curve 398, first quotient iteration | 372 | 8 | 253 |
| curve 385 `M29`, natural weight one | 394 | 0 | 111 |
| curve 385 `M29`, natural weight two | 2,722 | 0 | 554 |

For every chart it recomputes before-search quantities: canonical depth, raw,
integral and reduced coefficient heights, the exact invariant pair, the
horizontal reduced-to-raw matrix and determinant, stage-map sizes, modular
density, and the reduced-coordinate/canonical-height distortion on the signed
known basis.  When `hyperellred` lands on the other LLL tie, a small exact
unimodular equivalence recovers the horizontal coordinate used in the sealed
ledger.

The audit also imports the four earlier published-MW17 positive controls and
the three rank-29 control searches.  Their 394 compact chart records retain
the common presearch fields and their separately sealed verification ledgers
supply exact prefix-gain labels:

| control order | displayed jump | blindly recovered quotient rank |
|---|---:|---:|
| R17 control A | 8 | 8 |
| R17 control B | 9 | 9 |
| R17 control C | 10 | 10 |
| R17 control D | 11 | 9 |
| curve 12 from rank 29 | 12 | 10 |
| curve 356 from rank 29 | 12 | 12 |
| curve 385 from rank 29, before iteration | 12 | 4 |

No common target-free scalar is stable across these seven controls and the
four detailed positive chart orders at even the modest requirement “same
direction and worst-case AUC at least 0.7.”  Depth, integral coefficient
height, reduced coefficient height, and modular density all reverse direction
somewhere.  Within the four detailed orders, the best individual effects are
weak and change identity or sign:

- rank 28: depth, AUC `0.626`;
- curve 385 iteration: modular density, AUC `0.654`;
- curve 398 initial twelve: smaller depth, AUC `0.743` on only twelve charts;
- curve 398 iteration: smaller reduced coefficient height, AUC `0.605`.

So coefficient height, reduction gain, map size, sieve density, and depth do
not individually predict first independent appearance within the selected
strata.  The expanded eleven-order comparison includes the full `+8`, `+9`,
and `+10` recoveries, the `9/11` and `10/12` cases, curve 356's `12/12`, and
curve 385's pre-iteration `4/12`; it is not an inference from the spectacular
later cases alone.

The clean rank-28 ledger supplies the decisive target-relative check.  Across
64 charts and all 40 blindly found non-generic points, reduced-coordinate
visibility at height `100000` reproduces the sealed source relation exactly:

```text
true source-visible pairs   40
false positives              0
false negatives              0
true source-invisible pairs 2520
```

Those coordinate-visible quotient masks reproduce every prefix rank gain and
the final quotient rank eleven.  More importantly, the minimum
`hhat(2R-Q)/4` over the forty targets perfectly separates productive and
nonproductive charts:

```text
largest productive minimum       37.9881458706605
smallest nonproductive minimum   42.7115771810219
strict gap                         4.7234313103614
AUC for smaller centered height    1.0
```

This is a retrospective fixed-ledger statement, not a prospective success
probability.  It identifies the operative variable: closeness of the new
point to `Q/2`, followed by exact rational compression through the reduced
coordinate.

The curve-385 `M29` failure is compatible with this mechanism.  In its natural
weight-one and weight-two stages, larger depth predicts *absence of already
known finite hits* with AUC `0.739` and `0.690`, respectively.  That is exactly
the empty-ball effect.  It does not manufacture an unknown point inside the
ball.  The 1,445 returned points show that thousands of extra coordinate
charts can expose old `M29` while adding no arithmetic direction.

## Builder pilot on curve 385

The mechanism was also tried as an actual small builder, rather than only as a
retrospective ranking.  The protocol starts from the frozen `M29` lattice but
does not reuse the old 43 classes or enumerate the next quotient-weight shell.
It takes 256 deterministic seeds in the full `2^29` parity space and performs
strict one-bit ascent using the current canonical-height form.  The 256 seeds
land at 256 distinct local maxima; the mean ascent length is `2.09` and the
maximum is six steps.

For the 32 deepest unsearched maxima whose shortest representatives agree at
rounding scales `10^5` and `10^6`, the builder constructs and reduces the
quartic, records its exact horizontal map, and calibrates

\[
 2D(c)+\min_{P\in\{\pm P_1,\ldots,\pm P_{29}\}}
 \left(h(s(P))-\tfrac12\widehat h(2P-Q_c)\right)-\log(100000).
\]

This is deliberately marked as an empirical signed-basis margin, not a global
height bound.  The 24 best compression candidates are then thinned by maximin
distance in the current lattice torus, and a 16-chart order is frozen before
search.  The selected depths range from `42.4865` to `44.3395`, compared with
maxima `40.4667` and `42.7371` in the completed natural weight-one and
weight-two stages.  None of the 32 pool charts repeats any of the 3,116 old
natural-weight charts or any of the 3,460 previously searched base points in
the full source ledger.

The committed 16 charts were searched at the unchanged reduced-coordinate
height `100000`, 15-second wall limit, 1 GB GP stack, and zero retries.  All
sixteen completed, all returned zero finite points, and exact group
classification left `M29` unchanged.  Descriptively, this is `0/16` charts
with an old finite hit versus `665/3116` in the preceding natural-weight
stages; the fixed-table lower-tail Fisher value is `0.02173`.  The arms were
selected by different deterministic policies, so that number is not a
randomized inferential p-value.

The pilot therefore does one useful thing cleanly: it finds genuinely fresh,
deeper regions and suppresses the old-point rediscovery that dominated the
3,116-chart extrapolation.  It did not find a new arithmetic direction.  This
supports the empty-ball half of the mechanism, while again showing that the
existence of a new point in the empty region is a separate gate.  Because all
signed-basis margin estimates remain negative, the result does not yet justify
a larger campaign; the next builder step is a tight global or proved-shell
bound for `C_s`, followed by a precommitted multi-height test.

## Rank-32 scheduler consequence

The next scheduler should optimize coverage, not mask count.

1. Recompute large Voronoi/Delaunay holes for the **current** lattice.  Do not
   use low Hamming weight in an old quotient basis as a proxy for depth.
2. Reduce every candidate quartic before search and compute a tight local
   height-comparison constant for its reduced coordinate.
3. Rank by an old-point exclusion margin such as

   \[
   2D(c)-C_s-\log B,
   \]

   where `B` is the reduced-coordinate bound.  A certified positive margin
   excludes non-pole points of `M` from the search box.  Until `C_s` is made
   rigorous, enumerate a proved short-vector shell and label the margin as a
   numerical scheduling score only.
4. Enforce geometric diversity between chosen midpoints so that many charts
   do not probe the same part of `M tensor R / M`.
5. After every rank gain or finite-index enlargement, recompute the lattice,
   holes, maps, and distortion bounds.  The old ordering is invalid.

Residual Selmer/class-group information remains the arithmetic gate deciding
which fibers deserve this search.  The midpoint/exclusion calculation decides
how to search a chosen fiber.  The two roles should remain separate.

## Evidence and replay

The compressed chart census and posthoc oracle audit are
[`half_lattice_height_compression_analysis_v1.json.gz`](../../artifacts/generated-results/elliptic-curves/half_lattice_height_compression_analysis_v1.json.gz),
with SHA-256

```text
b463e3af5f9262a774f06c189eab66533c185d7a6fc7dd427aa3d49f356fa888
```

Rebuild or check it with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/analyze_half_lattice_height_compression.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/analyze_half_lattice_height_compression.sage --check

python3 -m unittest -v \
  elliptic-curves/tests/test_half_lattice_height_compression.py
```

The curve-385 builder protocol and bounded pilot are
[`curve385_height_compression_pilot_protocol_v1.json`](../../artifacts/generated-results/elliptic-curves/curve385_height_compression_pilot_protocol_v1.json)
and
[`curve385_height_compression_pilot_blind_v1.json`](../../artifacts/generated-results/elliptic-curves/curve385_height_compression_pilot_blind_v1.json),
with SHA-256 values

```text
93a1e736154e05d21da809b21fff8e669c6039e0b74c8411fe3bc62e9795fef9  protocol
94f3ad8fb7f02961ccfab049f1ad01c3846706eb066a2aea9995bf783a52f3f7  result
```

Replay its phase boundary and verification with:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_height_compression_pilot.sage --phase build

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_height_compression_pilot.sage --phase search

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_height_compression_pilot.sage

python3 -m unittest -v \
  elliptic-curves/tests/test_curve385_height_compression_pilot.py
```

The quartic formulas, rational maps, invariant identities, and coordinate-
visibility replay are exact.  Canonical-height decimals and AUC values are
high-precision numerical diagnostics.  No bounded miss is a point-absence,
saturation, Selmer, exact-rank, or rank-upper-bound result.

<!-- status-consumer: EC-HALF-LATTICE-HEIGHT-COMPRESSION 3baeaf370aec751c -->
<!-- status-consumer: EC-CURVE385-HEIGHT-COMPRESSION-BUILDER-PILOT c0d6f2d67018def4 -->

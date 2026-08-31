# Elkies 2026 rank-17 paper: certified integration and search handoff

Date: 2026-08-31

Source: Noam D. Elkies, *An elliptic K3 surface X/Q(t) with Mordell-Weil
rank 17, I: Formulas for X and base changes of ranks 18 and 19*,
[arXiv:2608.25406v1](https://arxiv.org/abs/2608.25406), submitted 2026-08-26.

## Outcome

The paper is now integrated into the exact H3 theorem package and the active
search framework. The repository has:

1. the compact published `A(t),B(t)` model;
2. all seventeen published sections, stored as the displayed x-coordinates,
   `P1`'s y-coordinate, and the sixteen quadratic chord descriptions;
3. an exact Mobius and Weierstrass identification with q12/orbit5867;
4. exact positive controls at all four disclosed high-rank parameters, with
   quotient gains `8,9,10,11` beyond the generic rank 17;
5. a complete compact-`t`, height-10000, three-ensemble Nagao calibration;
6. a fail-closed residual 2-Selmer gate before expensive point search;
7. an explicit eighteenth section and rational parameter on the first conic;
8. both new sections on the paired cover, explicit `E0 -> t` maps, four
   independent `E0` generators, and a modular Mordell--Weil sieve.

The independent q8/q12 construction remains the canonical marked route proof.
The paper supplies the preferred arithmetic coordinate and the base changes.

## Published rank-17 data

The compact equation is

```text
y^2 = x^3 + A(t)*x + B(t),  deg(A,B,Delta)=(8,12,24).
```

The exact verifier reconstructs every ordinate from the published chords and
checks all seventeen Weierstrass identities. It also proves:

```text
height Gram determinant                    948
height-2 vector pairs                        0
unoriented height-4 pairs                 1311
published-basis -> pinned-R17 determinant    1
```

The coordinate matcher proves exact identities

```text
t = (a*u+b)/(c*u+d),
A_q12(u) = s^4 (c*u+d)^8  A(t),
B_q12(u) = s^6 (c*u+d)^12 B(t)
```

over `QQ`, with trivial twist. The compact published `t` chart is therefore
the default specialization and Nagao chart. The raw q12 coordinate remains a
construction regression; its images of the four disclosed fibres have
221--234-bit numerators.

## Calibration anchors

Use these before assessing a search change:

```text
certified rank >= 25: t = -2/377      quotient gain >= 8
certified rank >= 26: t = -308/251   quotient gain >= 9
certified rank >= 27: t = 2456/135   quotient gain >= 10
certified rank >= 28: t = -9529/5471 quotient gain >= 11
```

All four exact fibres equal the corresponding public minimal models in
Dujella's rank-record table. The imported exact public point sets have lengths
25, 26, 27 and 28. Finite-reduction pivoting keeps the generic seventeen
first, selects public complements of dimensions `8,9,10,11`, and certifies
each combined list in one matrix. These are unconditional rank lower bounds,
not exact-rank results. The rank-28 fibre is a positive control, not permission
for an indefinite four-point search. ICARM curves 302/351/356 are no longer
needed to infer the family.

## Compact-t scoring calibration

The accepted scorer partitions the 102 primes from 19 through 599
round-robin into three disjoint ensembles of 34. Each local Nagao contribution
is centered and population-standardized over the good fibres of
`P^1(F_p)`; singular fibres receive the prime mean and remain counted. The
primary ranking key is the weakest normalized block.

The complete height-10000 scan contains 121,589,944 primitive projective
parameters. The four positive-control ranks are:

```text
t=-2/377      rank  54,624   fraction 0.0449%
t=-308/251    rank 593,936   fraction 0.4885%
t=2456/135    rank 422,873   fraction 0.3478%
t=-9529/5471  rank  55,387   fraction 0.0456%
```

All pass the declared top-one-percent gate. The ranking remains heuristic; it
only decides which candidates deserve descent.

## Residual 2-Selmer gate

For a candidate fibre the required object is

```text
Sel_2(E_t) / <P1,...,P17>.
```

Rank 32 requires residual dimension at least 15. A completed unconditional
2-descent with residual dimension below 15 is therefore an exact rejection.
Only a result at least 15 authorizes two-cover solving or expensive point
search on the same minimal model. The BNF-free Kummer signature, norm-one
cubic elements, incomplete relation ledgers, `K(S,2)` envelopes and candidate
local classes do not pass this gate.

The first exact-backend attempt used PARI `ellrank` through Sage on the public
rank-28 fibre, with all 28 certified points supplied. It reached the strict
300-second limit at 230,338,560 bytes peak observed RSS without returning a
Selmer dimension. The pinned result is therefore
`INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN`: it is not an upper bound and it
does not authorize a point search. The rank-28 control already supplies eleven
quotient directions and would need four more for rank 32.

## First rank-18 cover

The published conic is

```text
u^2 = 4225*t^2 + 38636*t + 289444.
```

Branch reductions isolate a single nonzero mod-2 trace coset. A shortest trace
has height 10. Exact chord recovery factors the residual discriminant as

```text
(127170526080*h(t))^2 * (4225*t^2+38636*t+289444).
```

The stored section has the form

```text
x = x0(t) + x1(t)*u,
y = y0(t) + y1(t)*u,
```

and passes both coefficient identities modulo the conic equation and the
Galois trace identity. The conic parameterization is

```text
t = (289444-r^2)/(130*r-38636),
u = 65*t+r.
```

The paper's anti-invariance lemma then proves generic rank at least 18 over
`QQ(r)`. `search_elkies_2026_rank18_conic_nagao.py` pulls the existing local
R17 tables back along this rational map and sieves primitive rational `r`.

## Paired rank-19 cover

The second cover is

```text
u2^2 = 54756*t^2 - 3269604*t + 22473889.
```

Its second height-10 trace and polynomial section are also recovered exactly.
After parameterizing the first conic, the fibre product becomes

```text
v^2 = 54756*r^4 + 425048520*r^3 + 221786712628*r^2
      - 348786049427920*r + 74698868489239696.
```

A rational plane-cubic transformation gives the paper's curve

```text
E0: y^2 = x^3 + 1029367969*x^2 - 42900734074705920*x.
```

The exact data file records rational functions `r(x,y)` and `v(x,y)`, followed
by

```text
t  = (289444-r^2)/(130*r-38636),
u1 = 65*t+r,
u2 = v/(130*r-38636).
```

All identities reduce to zero modulo the `E0` equation. Four displayed points
are unconditionally independent, proving `rank(E0(Q)) >= 4`; the paper states
exact rank 4. The two anti-invariant characters give generic rank at least 19.

`search_elkies_2026_E0_mw_nagao.py` reduces the four-generator lattice once
per prime, enumerates its finite image in at most `#E0(F_p)` steps, caches the
`E0 -> t` Nagao symbol, and scans a bounded coefficient box. Finalists include
the exact rational `E0` point and exact `t` value.

## Active search order

1. Keep the four exact fibres in every score calibration report.
2. Search compact `t=a/b` with three or more disjoint prime ensembles and rank
   by weakest-block performance.
3. Compute the actual residual 2-Selmer quotient for a survivor. Reject it
   exactly when the residual dimension is below 15.
4. Solve two-covers or run expensive point search only after the same minimal
   curve passes that gate.
5. Promote rank 32 only after fifteen certified quotient directions and one
   finite-reduction independence matrix of rank 32.

The first-conic and `E0(Q)` lattice sieves remain supporting base-change
routes. They do not bypass the residual gate.

Nagao scores, bounded coefficient boxes, and absence of points are experiments,
not rank bounds.

## Parked paths

The following are provenance or regression paths, not current priorities:

- q8/q12 coefficient discovery after the exact endpoint;
- q12/orbit4484 equation lifting;
- fixed-corridor reverse lifts, q323, changed-zero reranking, and compiler
  optimization without a direct specialization use;
- ICARM fingerprint fitting as a way to infer the now-published family;
- ungated raw `ratpoints`, slope-box, and two-cover point searches;
- further searching of the rank-28 calibration fibre merely to rediscover its
  known 28 independent points.

## Reproduction

```bash
SAGE=/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python

$SAGE elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
$SAGE elkies-k3/scripts/match_h92_q12o5867_to_elkies_2026_qq.sage
python3 elliptic-curves/scripts/verify_elkies_2026_high_rank_calibrations.py
python3 elkies-k3/scripts/calibrate_elkies_2026_positive_controls_nagao.py
python3 elliptic-curves/cas/run_elkies_2026_rank28_residual_selmer.py \
  --timeout 300 --overwrite
$SAGE elkies-k3/scripts/verify_elkies_2026_rank18_first_cover.sage
$SAGE elkies-k3/scripts/verify_elkies_2026_rank19_paired_cover.sage

python3 elkies-k3/scripts/search_elkies_2026_rank18_conic_nagao.py \
  --numerator-bound 1000 --denominator-bound 1000 \
  --output artifacts/local/elkies-k3/elkies-r18-conic-nagao-h1000.json

python3 elkies-k3/scripts/search_elkies_2026_E0_mw_nagao.py \
  --coefficient-bound 6 --primes 11-97 \
  --output artifacts/local/elkies-k3/elkies-E0-mw-nagao-b6.json
```

The model, coordinate, control, and cover verifiers are exact replays. The
Nagao command is a complete bounded heuristic ranking. The residual-descent
command is an exact backend attempt, but its pinned timeout is incomplete and
authorizes no search. The last two commands are bounded supporting heuristic
searches.

<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-ELKIES-2026-HIGH-RANK-CALIBRATIONS 345b9fb977057133 -->
<!-- status-consumer: EC-K3-ELKIES-2026-NAGAO-POSITIVE-CONTROL f99c98cdb6b8cd7d -->
<!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE f5600026fe1e9656 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R18-COVER 6b4ee5bbc1afc01e -->
<!-- status-consumer: EC-K3-ELKIES-2026-R19-PAIRED f1e135d2ba803e80 -->

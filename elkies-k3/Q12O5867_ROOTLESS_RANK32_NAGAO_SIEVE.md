# q12/orbit5867 rootless rank-32 Nagao sieve

<!-- status-consumer: EC-CRT-BEAM-NONMONOTONE 5ae7e135da8cc80f -->

This bounded search utility reads the exact short-model coefficients `A(u)`
and `B(u)` from
`../artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json`. It does not read or
evaluate the seventeen large section coordinates.

For every usable prime `p`, the output contains all `p+1` entries of
`P^1(F_p)`: `#E_u(F_p)`, `a_p(u)`, the singular/bad-reduction flag, and the
repository Nagao contribution

```text
((2-a_p)/(p+1-a_p))*log(p).
```

The contribution is stored in deterministic units of `10^-12`. Primes that
divide a coefficient denominator of this particular rational short model are
reported and skipped. The rational scan includes primitive `(a:b)` with
negative, zero, positive, and infinite parameters. After each prime block it
keeps a capped Pareto frontier in every projective-height bucket, using the
cumulative score, number of good primes, number of bad primes, and height as
objectives. Individual block scores remain in the output for held-block
comparison and recalibration.

This is a search heuristic only. A high score is not a rank lower bound, and a
singular reduction modulo one table prime does not make the rational
specialization singular over `Q`.

## Reproduce

Run the unit tests and actual-family smoke benchmark:

```bash
.venv/bin/python elkies-k3/scripts/test_h92_q12o5867_rootless_nagao.py
```

Run a small complete sieve:

```bash
.venv/bin/python elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py \
  --numerator-bound 1000 \
  --denominator-bound 1000 \
  --height-bucket-width 100 \
  --keep-per-bucket 32,16,8 \
  --output artifacts/local/elkies-k3/q12o5867-rootless-nagao-h1000.json
```

For search-scale bounds, export the same complete projective tables once and
compile the C++ lookup loop:

```bash
.venv/bin/python elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py \
  --export-cpp-tables artifacts/local/elkies-k3/q12o5867-rootless-nagao-tables.txt \
  --tables-only

g++ -O3 -std=c++17 \
  elkies-k3/scripts/scan_h92_q12o5867_rootless_nagao.cpp \
  -o artifacts/local/elkies-k3/scan-q12o5867-rootless-nagao

artifacts/local/elkies-k3/scan-q12o5867-rootless-nagao \
  artifacts/local/elkies-k3/q12o5867-rootless-nagao-tables.txt \
  10000 10000 100 32,16,8 1000 \
  artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h10000.json
```

The C++ arguments are `TABLE NUMERATOR_BOUND DENOMINATOR_BOUND
HEIGHT_BUCKET_WIDTH KEEP_PER_BLOCK FINALISTS OUTPUT [PARAMETER_SCALE]`. The
optional scale defaults to one. It uses the same signed canonical enumeration
and cumulative-score Pareto comparison as Python. The test suite compiles it
and checks exact integer block scores and complete survivor order against the
Python implementation.

## Initial bounded scans

On 2026-08-26, the C++ `H=1000` replay exactly matched all 24 Python survivors,
their complete ordering, and every integer block score. The Python scan took
4.02 seconds and the C++ scan 0.078 seconds on the same host, a roughly 52-fold
hot-loop speedup.

The first `H=10000` C++ scan visited 121,589,944 primitive projective
parameters in 8.90 seconds. Its stage sizes were `121589944 -> 478 -> 244 ->
191`. The first five heuristic finalists were:

| `u=a/b` | height | total score | good/bad table primes |
|---:|---:|---:|---:|
| `677/3402` | 3402 | 17.750962864310 | 25/0 |
| `-267/847` | 847 | 17.734848257487 | 25/0 |
| `-5954/7203` | 7203 | 17.654935023846 | 25/0 |
| `1002/437` | 1002 | 17.654717179183 | 25/0 |
| `-2631/9481` | 9481 | 17.645021113227 | 25/0 |

The ignored raw outputs are
`../artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h1000.json` and
`../artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h10000.json`.
These survivors have not been specialized, minimized, or tested for quotient
rank; they do not pass the fifteen-direction promotion gate merely by appearing
in this table.

## Disjoint-prime rerank

Rerank the fixed 191-candidate `H=10000` population on primes from 199 through
499, disjoint from all discovery blocks:

```bash
.venv/bin/python \
  elkies-k3/scripts/rerank_h92_q12o5867_rootless_nagao_holdout.py
```

There are 50 requested primes. The exact rational short model is unusable at
241 and 293 because their coefficient denominators are not invertible, leaving
48 complete projective tables. The output stores every `P^1(F_p)` row and the
same deterministic `10^-12` score units as the discovery sieve.

The rerank changes the extreme tail substantially. The holdout leader is
`89/1320` (discovery rank 56); the discovery leader `677/3402` falls to
holdout rank 121. There is no overlap between the discovery and holdout top
10, and only `-7801/1463` occurs in both top 25 lists. Ranking candidates by
the worse of their discovery and holdout ranks, then by rank sum, gives:

| `u=a/b` | discovery rank | holdout rank | discovery score | holdout score |
|---:|---:|---:|---:|---:|
| `-7801/1463` | 11 | 7 | 17.416321216952 | 17.278560415646 |
| `-5954/7203` | 3 | 30 | 17.654935023846 | 15.798511026100 |
| `278/8301` | 17 | 32 | 17.189033813320 | 15.734474481043 |
| `-119/30` | 34 | 20 | 16.922544522781 | 16.039278390398 |
| `601/418` | 35 | 2 | 16.914670584874 | 18.152814221295 |

The ignored raw result is
`../artifacts/local/elkies-k3/q12o5867-rootless-nagao-h10000-holdout-p199-499.json`.
This held-out agreement is a robustness filter for a heuristic score. It is
not evidence for a new Mordell--Weil direction or for rank 32.

## Projective CRT/Gauss constructor

The bounded constructor escapes rectangular height boxes by beam-combining the
six strongest local symbols at each of the 25 discovery primes. Finite symbols
impose `a=r*b (mod p)`; an infinity symbol imposes `b=0 (mod p)`. Each mixed
projective congruence lattice is Gauss-reduced exactly before short basis
combinations are enumerated and reranked on the 48 heldout primes.

This frontier is intentionally heuristic, not a sieve in the mathematical
sense.  Partial rational-representative height can decrease after another CRT
constraint is added; the exact width-one counterexample in
[`../elliptic-curves/tests/test_crt_lattice.py`](../elliptic-curves/tests/test_crt_lattice.py)
therefore applies to this pruning pattern as well.  A parameter omitted from
the finite beam is untested, not excluded.

```bash
.venv/bin/python \
  elkies-k3/scripts/construct_h92_q12o5867_rootless_nagao_crt.py \
  --output artifacts/local/elkies-k3/q12o5867-rootless-nagao-crt-gauss.json
```

The default bounded run retains a 512-state beam and constructs 49,197 novel
parameters after exact deduplication against all 191 `H=10000` survivors. Their
heights range from `46608240423243634438573` to
`115546775872581831137164928315`. The population hash is
`25a39812e6ef4b657795a11cc4a1fbb019cd298ff7ed22ec29b7bd89beaf59b6`.

The strongest balanced discovery/holdout candidate is

```text
13620847610649119692320505 / 11499364581254575952531276
discovery score 20.041789412616; holdout score 18.823665186417
```

The largest heldout score is `19.691437765590`, at

```text
-440362423111512637494908834 / 247919936030657359551312069
discovery score 20.021425441509
```

Both improve materially on the bounded-box score tails, but this is selection
on local traces. Neither score is rank evidence or passes the required
fifteen-direction quotient test.

A second independent run adds the 18 usable primes in `[199,313]` to the 25
original construction primes, then validates on all 73 primes in `[503,997]`:

```bash
.venv/bin/python \
  elkies-k3/scripts/construct_h92_q12o5867_rootless_nagao_crt.py \
  --extra-construction-prime-min 199 \
  --extra-construction-prime-max 313 \
  --validation-prime-min 503 \
  --validation-prime-max 997 \
  --specialize-top 3 \
  --output artifacts/local/elkies-k3/q12o5867-rootless-nagao-crt-gauss-c199-313-v503-997.json
```

This constructs 45,753 novel parameters. Its population hash is
`32d115e4488ee86906939d075cd66735dc75188cb9edd595e9cdd517392e19e1`.
Equal construction score/good/bad signatures receive the same competition
rank. The leading parameter is rank 1 on both construction and fresh
validation:

```text
-641700991247763976719358613776538997953485850041
--------------------------------------------------
 115855879811331922221886372780501938285356394515

construction score 31.653720311480; validation score 18.303316935069
```

All 43 construction and 73 validation reductions used for this parameter are
good. The artifact exactly specializes the top three balanced parameters,
checks nonsingularity and the `c4^3-c6^2=1728*Delta` identity, and supplies
denominator-cleared integral short models. It makes no minimality claim and
launches no ratpoints or other point search.

The default three prime blocks contain 25 usable primes between 19 and 197.
Custom semicolon-separated blocks accept explicit primes and inclusive ranges,
for example:

```bash
--prime-blocks '19-83;89-151;157-197'
```

The JSON output pins the input hash, requested and usable primes, rejected
model-denominator primes, every local table, all stage sizes and timings, and
the final ordered candidates. Promotion must separately construct and
minimize a specialization, evaluate and saturate the generic rank-17 subgroup,
and establish at least fifteen independent quotient directions.

## Exact point-search backends

This section records historical bounded searches. Under the Elkies 2026
descent-first policy these backends are no longer pre-descent discovery tools.
Every live point-search command requires a completed unconditional gate for
the same parameter and global minimal model. The specialization adapter itself
remains ungated because it performs no point search.

The specialization adapter and gate-protected CPU entry points are:

```bash
.venv/bin/python elliptic-curves/scripts/specialize_q12o5867_candidate.py \
  --a NUMERATOR --b DENOMINATOR --output SPECIALIZATION.json

.venv/bin/python elliptic-curves/scripts/probe_q12o5867_ratpoints.py \
  --input SPECIALIZATION.json --residual-selmer-gate GATE.json \
  --height 1000000000 \
  --denominator-bound 100000 --timeout 120 --output PROBE.json

.venv/bin/python elliptic-curves/scripts/probe_q12o5867_section_charts.py \
  --input SPECIALIZATION.json --residual-selmer-gate GATE.json \
  --pair-mode all --include-multiplicative \
  --height 100000 --denominator-bound 100000 \
  --per-chart-timeout 20 --output CHARTS.json
```

The first search completes the square exactly. The second searches all 136
charts `x=x_i+(x_j-x_i)X` and all 17 charts `x=x_i*X`, maps every result back
to the global minimal model, verifies it exactly, and immediately measures its
finite-quotient escape from the ordered generic rank-17 subgroup. Eight
complete 153-chart CPU runs through `H=D=10^4` and `H=D=10^5` rediscovered all
17 baseline abscissae and no novel point.

A CUDA widening used Samuel Li's public `ratpoints_gpu` commit
`904072e71a4f3b2896c4d6e88d31943166d05e38`, built locally with CUDA 13.3.
The executable SHA-256 was
`39361d73c83dcd2ab7125a1f22836cf15514302b4e152138962aba010078bcfa`.
The external tool was prepared outside the repository with:

```bash
git clone https://github.com/wgxli/ratpoints-gpu.git /tmp/q12-ratpoints-gpu
git -C /tmp/q12-ratpoints-gpu checkout 904072e71a4f3b2896c4d6e88d31943166d05e38
make -C /tmp/q12-ratpoints-gpu -j8 NVCC=/usr/local/cuda-13.3/bin/nvcc
```

On the leading fresh-validation CRT parameter, all 153 normalized charts were
searched through `H=D=500000`:

```bash
.venv/bin/python elliptic-curves/scripts/probe_q12o5867_section_charts.py \
  --input artifacts/local/elliptic-curves/q12o5867-specializations/q12o5867-specialization-crt3-v1.json \
  --residual-selmer-gate PASSING_GATE_FOR_THE_SAME_FIBRE.json \
  --pair-mode all --include-multiplicative \
  --height 500000 --denominator-bound 500000 \
  --per-chart-timeout 30 \
  --ratpoints /tmp/q12-ratpoints-gpu/ratpoints_gpu \
  --ratpoints-library /usr/local/cuda/targets/x86_64-linux/lib \
  --output artifacts/local/elliptic-curves/q12o5867-section-chart-probes/q12o5867-section-charts-crt3-v1-gpu-all-h500k-d500k.json
```

That run covered 25 times as many projective sites per chart as the `H=D=10^5`
stage, completed all charts without timeout, again recovered exactly the 17
baseline abscissae, and found no novel point. The result is a bounded search,
not a rank upper bound. Its JSON SHA-256 is
`f59ce596bb706f79027e24fd92eda51df961c34fb80ac78c5e6315a23086fad6`.

## Skew projective chart and fresh-prime ensemble

The inherited pencil coordinate is strongly unbalanced: after taking decimal
absolute values, the coefficients of both binary forms decrease by about
`7.7` orders of magnitude per parameter degree. The integer scale `50874487`
is the nearest prime above the elementary minimax balance point `10^7.7065`.
This observation does not change the family; it supplies the projective chart

```text
u = 50874487*v.
```

The C++ scanner now accepts this optional final scale argument, requires it to
be invertible at every table prime, scores the exact old parameter, and keeps
height buckets in `v`. Its regression test independently rescores every
returned old projective pair in Python. The skew scan is reproduced by:

```bash
g++ -O3 -std=c++17 \
  elkies-k3/scripts/scan_h92_q12o5867_rootless_nagao.cpp \
  -o artifacts/local/elkies-k3/scan-q12o5867-rootless-nagao-skew

artifacts/local/elkies-k3/scan-q12o5867-rootless-nagao-skew \
  artifacts/local/elkies-k3/q12o5867-rootless-nagao-tables.txt \
  10000 10000 100 32,16,8 10000 \
  artifacts/local/elkies-k3/q12o5867-rootless-nagao-skew-s50874487-h10000.json \
  50874487
```

It scores another 121,589,944 primitive chart parameters and retains 196.
The local artifact SHA-256 is
`e9892257b5c065beaad9130ce54a288f5a2f9cddd8db7e5beea7a6f24dbb5537`.
The same 48-prime heldout rerank has SHA-256
`345aee4e798060c55017ee88398bb78469c64be73ec49c9b8223dd9a02fc7aa3`.

Discovery/CRT selection was then replaced by a genuinely fresh-prime
ensemble. The reranker reconstructs the complete CRT populations and checks
their pinned hashes, unites all bounded populations exactly, centers and
standardizes the local contribution separately on every projective prime
table, and uses six round-robin prime blocks. Selection uses the 71 usable
primes in `[1009,1499]`; a fixed shortlist is confirmed on 63 disjoint usable
primes in `[1511,1999]`. The enlarged run is:

```bash
.venv/bin/python \
  elkies-k3/scripts/rerank_h92_q12o5867_rootless_nagao_ensemble.py \
  --additional-box \
    artifacts/local/elkies-k3/q12o5867-rootless-nagao-skew-s50874487-h10000.json \
  --include-specialized \
  --output \
    artifacts/generated-results/elliptic-curves/q12o5867-fresh-prime-ensemble-with-skew-shortlist-v1.json
```

The exact union has 95,337 parameters and ordered-population SHA-256
`9808a6d972e523ea0e3dfa354bd13889808ba15ea5873e48eae1eb19dc910c33`.
Eleven of 68 fixed selections have both positive confirmation total `z` and a
positive confirmation one-standard-error block lower score. The generated
artifact SHA-256 is
`70a3c7b5aab15ef4beaba08984eec540efca437e203e8cd4c8678b67f7cda206`.
The strongest three cost-relevant survivors are:

| `u` | source | selection `z` / LCB | confirmation `z` / LCB |
|---:|:---|---:|---:|
| `-145513755362346429303616429/240417263878420561750868108` | first CRT | `2.948 / 1.012` | `2.860 / 0.803` |
| `-129373820441/7819` | skew chart | `3.033 / 0.904` | `1.934 / 0.468` |
| `-59/101` | original box | `2.403 / 0.736` | `1.370 / 0.203` |

A subsequent four-scale widening repeated the same complete box at prime
scales `30000001`, `50000017`, `70000027`, and `100000007`. The five skew
boxes retain `196+190+188+211+202` rows. Their union with the original box and
both CRT populations has 96,127 unique parameters, one cross-source duplicate,
and ordered-population SHA-256
`b0ffb3301851a23421b6b4e829ed71410a03996c451d3fe8a967838fea2dcd20`.
The same two independent prime stages leave 17 positive confirmation-gate
rows; the stored top 16 are in
`../artifacts/generated-results/elliptic-curves/q12o5867-fresh-prime-ensemble-multiskew-shortlist-v1.json`,
SHA-256
`43e7fcb83158c62311bad596bdd07bda818c0d304070fdd9640719569ab38cf0`.
The best newly added cheap row is
`u=-135270004509/3260`, with selection `z/LCB=1.653/0.480` and independent
confirmation `1.833/0.493`. It also specializes with exact generic lower bound
17; its direct `10^9/10^5` ratpoints probe returned no finite abscissa.

The normalized slope search now also accepts an explicit ratpoints-compatible
CPU/CUDA executable. It clears each exact quartic by a square, reconstructs
both ordinate signs from every returned abscissa, verifies the transport, and
classifies exact signed sums of two listed basis sections before the quotient
gate. On each of the two cheap skew leaders, all 17 distinct-base charts were
searched through `H=D=10^6`. They returned respectively 28 and 29 exact points
beyond the literal displayed basis. Every one is exactly a signed pair sum of
the 17 generic sections, so none is a new Mordell--Weil direction. The two
self-consistent local artifact SHA-256 values are
`54cf679f7207d519298cc8361e4bb87a17aedfbeaae31b7018de92cc0b095564`
for `u=-129373820441/7819` and
`966bf973bc5897101345b21cb7872b6aec749ebe945d0403d991d2e58ecf2285`
for `u=-135270004509/3260`. This exact relation classification is stronger than
the otherwise one-sided finite-quotient non-escape result.

The same complete 17-chart `H=D=10^6` slope pass on the strongest arithmetic
ensemble survivor
`u=-145513755362346429303616429/240417263878420561750868108`
returned 29 points beyond the displayed basis; again all 29 are exact signed
pair sums. Its local artifact SHA-256 is
`9320ebf9025539d8829b0b5819446e0b3919874a2fad82ef36482ef8c1482242`.
The next two confirmed CRT rows,
`-380745322650012819400595381/231428974777023315028225950` and
`2562422866289678369255167/39602073966699250017405709`, likewise return only
exact signed pair sums (`28/28` and `29/29`). Their artifact SHA-256 values are
`88fa1db25a6b95c51a04c2dbc580790bcbf41cafca593952148b29ede40daaac`
and
`0cb09102a09f9171956a551e4c1d0a36702d11f5079f0d73b0233e86ce438da1`.

All three specialize exactly and retain the certified ordered generic rank-17
subgroup. Direct ratpoints searches through numerator height `10^9` and
denominator `10^5` found no finite abscissa on any of them. Four normalized
slope-quartic charts on `-59/101`, through numerator `10^8` and denominator
`10^4`, likewise produced no new point. These are bounded negative searches,
not dependence or rank upper bounds.

The earlier discovery-score promotions do not survive this calibration well:
for example the first 43-prime CRT construction leader has fresh standardized
`z=-0.268` and block LCB `-0.376`. Further expensive promotion should use the
confirmed ensemble rather than the discovery extreme.

## Discriminant-root projective charts

Six exact `PGL_2(Q)` charts center bounded boxes on distinct real-root
intervals of the binary discriminant.  Each chart scans the same
121,589,944 primitive pairs as the original `H=10000` box.  The normalized
matrices `[alpha,beta,gamma,delta]`, acting by
`u=(alpha*v+beta)/(gamma*v+delta)`, are:

```text
[28760230122, 31143670635,  -500,  -500]
[54857896711, 57520460246, -1000, -1000]
[41518313879, 44164218975, -1000, -1000]
[17207900084, 19267286851,  -500,  -500]
[ 7637406659,-4721226972,   1000,  1000]
[57125258651, 30687267094,  1000,  1000]
```

Their complete survivor counts are `190,183,195,190,204,202`.  The 1,164
new survivors are pairwise disjoint; against the five scalar-chart survivor
sets, only the common projective point `u=infinity` overlaps.  The audit
exactly verifies 18 distinct discriminant real-root sign-change brackets and
pins every chart matrix, scan hash, population, and confirmed candidate in
`../artifacts/generated-results/elliptic-curves/q12o5867-binary-gl2-chart-audit-v1.json`,
SHA-256
`c956e28fe42e9964f5f40f190009ea802804e890e02f3545c97ec07e8b149696`.

Adding these charts produces a 97,290-parameter exact union.  The same
disjoint selection and confirmation prime ranges leave 31 positive
confirmation-gate rows, 14 supplied by the new charts.  The resulting
ensemble is
`../artifacts/generated-results/elliptic-curves/q12o5867-fresh-prime-ensemble-gl2-shortlist-v1.json`,
SHA-256
`a2367dee7be4629a3875943ef96f451e067718a1394eb8388e57bd1d829e4207`.
The strongest new row is

```text
u = -19629011813143/307000
selection z/LCB = 1.368/0.349
confirmation z/LCB = 3.461/1.125
```

The cheapest confirmed new row is `u=-7509891247/200`, with selection
`z/LCB=2.046/0.522` and confirmation `1.414/0.182`.  Both specialize
exactly and retain the ordered generic rank-17 subgroup.  Direct ratpoints
searches through numerator height `10^9` and denominator `10^5` returned no
finite abscissa on either curve.  On the strongest row, the first complete
17-chart normalized-slope search through `H=D=10^6` returned 28 exact points
beyond the displayed basis, all 28 exactly signed sums of two generic
sections.  Its local artifact SHA-256 is
`7a23706055250dbc79c4e1e99d096f05d1f29d4c8988b0a98864da3a84375c2b`.
Thus these completed searches certify no rank jump; they are bounded negative
screens, not rank upper bounds.

## Exact PARI two-cover backend

`../elliptic-curves/scripts/probe_q12o5867_pari_two_cover.py` is a bounded,
owned-process descent backend. It now refuses to start without a passing
same-parameter, same-minimal-model residual-Selmer gate. It validates the exact specialization, runs
PARI `ell2cover` under wall-time, RSS, and stack caps, searches any returned
quartics with `hyperellratpoints`, replays the covering maps exactly, and sends
every mapped point directly to the finite-quotient gain gate.

The live toy cover/map regression passes, but the real resource audit is
negative. The 835-digit `a6` specialization `u=-267/847` remained inside
`ell2cover` for 30 seconds without returning a cover. Even projective
infinity, whose minimal `a6` has only 207 digits, exhausted a 16 GB PARI stack
before returning a cover; peak observed RSS was `16013324288` bytes. The
corresponding local artifact SHA-256 values are
`eb1a01b187e2f3f6e93816ebe6d4b00d2a243277ef762b7cfd5f3ea80a4a7079`
and
`37878490222419581fa900c32f7bffdcc79d2559cbc86a5e243fdd0f18afeaa5`.
Thus generic PARI descent is operational but not viable for these fibres on
the present 32 GB host. A field-specific BNF-free descent or a higher-memory
Magma computation remains the genuine descent route; timeout and stack
exhaustion provide no rank information.

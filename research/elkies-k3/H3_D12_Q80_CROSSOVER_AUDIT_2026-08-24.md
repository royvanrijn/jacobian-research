# H3 D12 to Q80 crossover audit — 2026-08-24

## Result

There is no cheap direct crossover from the current exact H3 `D12/MW5`
equation frame into the retained Q80 corridor.

The hoped-for target

```text
D12 -> Q80-LQ7-A1/MW16 --final q6--> rootless/MW17
```

has exact intersection degree

```text
F_D12 . F_candidate1 = 370213639961146392704841338.
```

In the current D12 marking its cheapest shortest-section decomposition has

```text
P.O = 108545315264050223557720445635678277043049764366553849,
height = 217090630528100447115440891271356554086099528733107699,
local D12 correction = 3,
vertical support = all 12 simple D12 components.
```

Thus candidate1 is not merely non-cheap: it is unusable as an equation-level
crossover.  The small CM24 marking `-P3`, with `P.O=0`, is entirely a
specialization/re-chambering collapse and cannot be used as its generic H3
D12 section.

The best stage in the complete audit is the common-prefix `D7+D5/MW5` frame
`Q80-P2`, but it still has

```text
d = 16328023738263177,
P.O = 211141449796038353605555467406487.
```

Among the eight retained low-q suffix stages, the best is the first `escape`
child `Q80-LQ1-D7D4`, with

```text
d = 34707400601887865,
P.O = 954004515294656151571507048228790.
```

Consequently no Q80 stage offers a credible replacement for the existing

```text
D12 -> A11 -> 2A5 -> ... -> A1 -> MW17
```

suffix.  The already-solved Q80 terminal q6 compiler remains useful as a
compiler pattern and CM24 regression, not as a direct H3 D12 shortcut.

## Exact transport

Every comparison uses one explicit pinned Neron--Severi marking:

```text
retained Q80 stage
  -> initial Q80 frame
  -> pinned R17 frame
  -> current equation-side H3 D12 frame.
```

The Q80-to-pinned map is
`data/fibrations/kumar_q80_rootless_target_to_q80_ns_transport.txt`.  The
pinned-to-D12 map is obtained by composing the two exact matrices in the
current `q24-equation-d13-to-pinned-r17.json` certificate.  This is a literal
integral transport.  No `qfisom`, fresh neighbor search, ADE-label match, or
CM24 class is used to identify a generic fibre.

The audit includes the Q80 start, both common-prefix q4 children, and all
eight retained low-q stages through the new rootless endpoint: eleven fibre
classes in total.

## Complete D12 score table

Here `c` is the exact D12 local correction.  `v` is the number of nonzero
simple D12 components in the cheapest shortest-section decomposition; it is
12 in every row.  The pole triple is `(deg Z, max deg X, max deg Y)` on the
current minimal D12 equation.

The RR column is a planning estimate, not an executed resolved-RR dimension:

```text
2 + 2(P.O) + connected vertical layers.
```

This connected-layer rule reproduces the executed orbit42 ambient exactly:
`2+2*3+1=9`.  A genuine crossover would still require resolved local modules
and a two-dimensional kernel.

| Q80 stage | generic frame | `d` | `P.O` | MW height | `c` | `v` | expected RR ambient | pole degrees `(Z,X,Y)` |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `Q80-START` | `E6+D5+A3/MW3` | `492925390917124368` | `192428162218790780838702471346573659` | `384856324437581561677404942693147319` | 3 | 12 | `384856324437581563172687111654727084` | `(192428162218790780838702471346573659, 384856324437581561677404942693147322, 577284486656372342516107414039720983)` |
| `Q80-P1` | `D9+A4/MW4` | `122519354993495723` | `11888187788007171357696926270156735` | `23776375576014342715393852540313474` | 0 | 12 | `23776375576014343087054574933998526` | `(11888187788007171357696926270156735, 23776375576014342715393852540313474, 35664563364021514073090778810470211)` |
| `Q80-P2` | `D7+D5/MW5` | `16328023738263177` | `211141449796038353605555467406487` | `422282899592076707211110934812975` | 3 | 12 | `422282899592076756741938926409014` | `(211141449796038353605555467406487, 422282899592076707211110934812978, 633424349388115060816666402219467)` |
| `Q80-LQ1-D7D4` (`escape`) | `D7+D4/MW6` | `34707400601887865` | `954004515294656151571507048228790` | `1908009030589312303143014096457584` | 0 | 12 | `1908009030589312408427420636056672` | `(954004515294656151571507048228790, 1908009030589312303143014096457584, 2862013545883968454714521144686376)` |
| `Q80-LQ2-A6A4` (`orbit424`) | `A6+A4/MW7` | `189408214646705618` | `28412135291541929039551662316710043` | `56824270583083858079103324633420087` | 3 | 12 | `56824270583083858653670452305859370` | `(28412135291541929039551662316710043, 56824270583083858079103324633420090, 85236405874625787118654986950130135)` |
| `Q80-LQ3-A6A3` (`orbit1222`) | `A6+A3/MW8` | `2785831137539739596` | `6146329370162830518246866939191358275` | `12292658740325661036493733878382716551` | 3 | 12 | `12292658740325661044944513047925416349` | `(6146329370162830518246866939191358275, 12292658740325661036493733878382716554, 18438988110488491554740600817574074831)` |
| `Q80-LQ4-A4A2A1` (`q6_7774`) | `A4+A2+A1/MW10` | `1574227338760046874097` | `1962642294983273936940622262285196209473176` | `3925284589966547873881244524570392418946356` | 0 | 12 | `3925284589966547873886019920833019722917164` | `(1962642294983273936940622262285196209473176, 3925284589966547873881244524570392418946356, 5887926884949821810821866786855588628419534)` |
| `Q80-LQ5-A3A2` (`q4_1938`) | `A3+A2/MW12` | `49375786943066877833412` | `1930786768160600427726829376352463361158429919` | `3861573536321200855453658752704926722316859839` | 3 | 12 | `3861573536321200855453808533452997026259874072` | `(1930786768160600427726829376352463361158429919, 3861573536321200855453658752704926722316859842, 5792360304481801283180488129057390083475289763)` |
| `Q80-LQ6-4A1` (`q4_6855`) | `4A1/MW13` | `1201925838129816927955115` | `1144093704726714833878682602110381758656249141058` | `2288187409453429667757365204220763517312498282117` | 3 | 12 | `2288187409453429667757368850245714323476056749991` | `(1144093704726714833878682602110381758656249141058, 2288187409453429667757365204220763517312498282120, 3432281114180144501636047806331145275968747423180)` |
| `Q80-LQ7-A1` (`candidate1`) | `A1/MW16` | `370213639961146392704841338` | `108545315264050223557720445635678277043049764366553849` | `217090630528100447115440891271356554086099528733107699` | 3 | 12 | `217090630528100447115440892394394370235718074930638021` | `(108545315264050223557720445635678277043049764366553849, 217090630528100447115440891271356554086099528733107702, 325635945792150670673161336907034831129149293099661553)` |
| `Q80-LQ8-ROOTLESS` | `rootless/MW17` | `54769401772670506721527679106` | `2375648853817871228249584803714803380511769640842682701216` | `4751297707635742456499169607429606761023539281685365402436` | 0 | 12 | `4751297707635742456499169607595748963037101387332475074681` | `(2375648853817871228249584803714803380511769640842682701216, 4751297707635742456499169607429606761023539281685365402436, 7126946561453613684748754411144410141535308922528048103654)` |

The exact ranking by `(d, expected ambient, P.O, vertical layers)` is

```text
Q80-P2
< Q80-LQ1-D7D4
< Q80-P1
< Q80-LQ2-A6A4
< Q80-START
< Q80-LQ3-A6A3
< Q80-LQ4-A4A2A1
< Q80-LQ5-A3A2
< Q80-LQ6-4A1
< Q80-LQ7-A1
< Q80-LQ8-ROOTLESS.
```

## CM24 section complexity

The specialized scores are useful only as compiler regressions.  Their
smallness does not survive generic transport to H3 D12.

| retained stage | CM24 horizontal | CM24 `P.O` / height | CM24 module complexity | status |
|---|---|---|---|---|
| `escape` | old q12 section | `0 / 3` | E6 exact jet | exact over `QQ(sqrt(-6))` |
| `orbit424` | rational 2-torsion chord | `0 / 0` | torsion chord | exact over `QQ(sqrt(-6))` |
| `orbit1222` | saturated A7 chord | `1 / 25/8` | raw saturated A7 module | exact `QQ(sqrt(-3))` parent |
| `q6_7774` | `P3` (MW word L1 `1`) | `0 / 8/7` | A1 plus connected A6 quotient | exact over `QQ(sqrt(-3))` |
| `q4_1938` | `-P1+P2+2P3` (L1 `4`) | `1 / 12/5` | connected A4 quotient | exact over `QQ(sqrt(-3))` |
| `q4_6855` | `2P1` (L1 `2`) | `0 / 3/5` | A3 middle-double | exact over `QQ(sqrt(-3))` |
| `candidate1` | `-P3` (L1 `1`) | `0 / 3/4` | whole A3 plus D4 outer complement | exact over `QQ(sqrt(-3))` |
| final q6 | `P2-P3` (L1 `2`) | `0 / 1` | whole A4 plus A5 residue `+/-4` | exact over `QQ(sqrt(-3))` |

## Reproduction

On this machine, from the repository root:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_h3_d12_q80_crossovers.sage
```

Expected terminal line:

```text
H3D12Q80_RESULT|stages=11|best=Q80-P2|best_d=16328023738263177|candidate1_d=370213639961146392704841338|candidate1_P.O=108545315264050223557720445635678277043049764366553849|candidate1_ambient=217090630528100447115440892394394370235718074930638021|artifact=artifacts/generated-results/h3-d12-q80-crossover-scores.json|status=PASS_EXACT_H3_D12_Q80_CROSSOVER_AUDIT
```

The generated JSON records every transported class, both composed integral
matrices, input hashes, all exact scores, the CM24 comparison, and the RR
estimate boundary.  It is reproducible output rather than a pinned source
artifact.

For this run its SHA-256 is
`5a2c4cb63a54b2f826840be6f51ec7e43de2e7935a331f9e1cc23c37e17b16c1`.
The checker source SHA-256 recorded by `MATH_STATUS.json` is
`902061b8cebf8ee10ce2990a04c3213d722890fbd5cc899b6db2aff1077d4faf`.

## Claim boundary

This closes the requested cheap crossover test.  It proves exact lattice
transport, intersection degrees, D12 MW heights/corrections, shortest-section
pole bounds, and vertical profiles.  It does not execute any new resolved
Riemann--Roch pencil, construct a Q80-stage equation over the generic H3
family, or identify the CM24 small sections with the enormous generic D12
sections.  The negative conclusion is a route-cost conclusion, not a theorem
that no other unsearched fibration can give a cheap crossover.

<!-- status-consumer: EC-K3-H3-D12-Q80-CROSSOVER-AUDIT 34f8f8038e591f00 -->

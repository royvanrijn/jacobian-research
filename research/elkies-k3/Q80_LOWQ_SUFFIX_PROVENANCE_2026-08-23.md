# Q80 low-q suffix provenance and backtracking map — 2026-08-23

## Purpose

This note records how the retained Q80 low-q suffix was actually produced. It is intended to prevent a recurring reconstruction mistake: names such as `7774`, `1938`, `6855`, and `candidate1` are **search-result labels**, not intrinsic names of sections or elliptic fibrations.

The generic lattice search selected a sequence of child frames first. The CM24 horizontal labels (`P3`, `-P1+P2+2P3`, `2P1`, `-P3`, `P2-P3`) were recovered **later**, after transporting those already-selected generic divisors into the CM24 specialization.

This distinction matters for characteristic-zero reconstruction. We should backtrack the selected search objects and their transports rather than repeatedly rediscover the final horizontal from node incidence.

For the repository-wide human-readable route map, including the H3 source route, H2 comparison route, reverse MW17 backtracks, and the role of Q80, see [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md). For the completed final-q6 reconstruction, see [`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).

## Canonical names for the Q80 suffix

New notes should identify a stage by **route + sequence + child geometry**. Historical search labels remain aliases for locating artifacts.

| canonical id | historical alias | neighbor | child |
|---|---|---:|---|
| `Q80-LQ1-D7D4` | `escape` | `q6` | `D7+D4/MW6` |
| `Q80-LQ2-A6A4` | `orbit424` | `q4` | `A6+A4/MW7` |
| `Q80-LQ3-A6A3` | `orbit1222` | `q4` | `A6+A3/MW8` |
| `Q80-LQ4-A4A2A1` | `q6_7774` | `q6` | `A4+A2+A1/MW10` |
| `Q80-LQ5-A3A2` | `q4_1938` | `q4` | `A3+A2/MW12` |
| `Q80-LQ6-4A1` | `q4_6855` | `q4` | `4A1/MW13` |
| `Q80-LQ7-A1` | `candidate1` | `q4` | `A1/MW16` |
| `Q80-LQ8-ROOTLESS` | `final q6` | `q6` | rootless/MW17 |

The two common-prefix children may analogously be called `Q80-P1-D9A4` and `Q80-P2-D7D5`.

The old ids are still valuable as **provenance**: `7774` is a candidate id, `1938` and `6855` are shell ids, and `candidate1` distinguishes one of three pairwise non-isometric A1/MW16 children. They should not be used as if they were mathematical names.

## 1. Where the alternate search really starts

The common Q80 prefix is

```text
E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5.
```

The low-q alternate search starts at this `D7+D5/MW5` frame.

The first important observation was that the old productive q12 MW class had a much cheaper representative in the **same MW coset**. Searching root corrections rather than changing the MW direction produced the q6 escape

```text
D7+D5/MW5
 --q6 (2,3)--> D7+D4/MW6.
```

The retained continuation then used q4 Weyl-orbit searches

```text
D7+D4/MW6
 --q4 orbit 424-->  A6+A4/MW7
 --q4 orbit 1222--> A6+A3/MW8.
```

Thus `orbit424` and `orbit1222` are already search/orbit labels. The later numeric names continue the same convention.

## 2. Exact retained search lineage

### `Q80-LQ4-A4A2A1` — historical `q6_7774`

Source frame: `A6+A3/MW8`.

Scoring artifact:

```text
data/fibrations/kumar_q80_a6a3_q6_chamber_scores.tsv
```

The retained row has

```text
candidate_id = 7774
q = 6, (a,b)=(2,3)
child = A4+A2+A1/MW10
root data = (7,28,30)
D.F = 2
P.O = 3
MW height = 219/28
vertical support = one fibre / two simple components
```

So `7774` is literally the q6 search candidate ID.

Retained vector:

```text
85,2699,1257,7718,3756,-41,3077,-4614,-6615,6032,2584,-1678,121,-736,-913,1,1165
```

### `Q80-LQ5-A3A2` — historical `q4_1938`

Source frame: the selected `Q80-LQ4-A4A2A1` child.

Scoring artifact:

```text
data/fibrations/kumar_q80_7774_q4_rank5_scores.tsv
```

The retained row is

```text
rank = 5
shell_id = 1938
eligible_id = 1935
q = 4, (a,b)=(2,2)
child = A3+A2/MW12
root data = (5,18,12)
D.F = 2
P.O = 1
MW height = 47/10
vertical support = one fibre / two components
```

Important: `1938` is a shell ID, **not the row rank**. There are earlier rows in the score table with the same root rank. The retained `1938` row was equation-friendly; its small connected vertical support became useful later.

Retained vector:

```text
-94,-1003,5298,4977,-1431,-1440,100,1,-1632,1893,1634,-1264,-4175,2248,-3111,1561,2842
```

### `Q80-LQ6-4A1` — historical `q4_6855`

Source frame: the selected `Q80-LQ5-A3A2` child.

Scoring artifact:

```text
data/fibrations/kumar_q80_1938_q4_4a1_scores.tsv
```

The retained row is the first-ranked `4A1/MW13` child:

```text
rank = 1
shell_id = 6855
eligible_id = 6853
q = 4, (a,b)=(2,2)
child = 4A1/MW13
root data = (4,8,16)
D.F = 2
P.O = 2
MW height = 19/3
vertical correction = 0
```

This is the origin of the name `6855`.

Retained vector:

```text
30693,-339,-2534,45446,10413,16390,-11527,5970,-18424,4193,21146,11296,25035,17925,-6032,4304,7717
```

### `Q80-LQ7-A1` — historical `candidate1`

Source frame: the selected `Q80-LQ6-4A1` child.

The balanced q4 shell on this new `4A1/MW13` frame has exactly three `A1/MW16` children. Exact integral-form tests show that the three are pairwise non-isometric and none is the older canonical-route `A1/MW16` frame.

The first retained child was historically called `candidate1`:

```text
q = 4, (a,b)=(2,2)
child = A1/MW16
D.F = 2
P.O = 2
MW height = 6
vertical correction = 0
```

Retained vector:

```text
21,671,-20182,-10366,27727,30558,5582,20831,-10195,-19691,6086,10389,20928,18651,16123,15473,-11496
```

The important historical point is that this child was **not** selected because of its later CM24 equation. It was a generic lattice-search child first.

### `Q80-LQ8-ROOTLESS` — terminal q6

Source frame: `Q80-LQ7-A1`.

Search/replay script:

```text
scripts/search_q80_new_lowq_final_q6_rootless.py
```

The rank-16 MW quotient made a naive radius-12 enumeration too large. The search decomposed q6 candidates by the A1 root pairing. The older canonical successful terminal q6 had A1 pairing `p=1`, so the new search used that as a guide.

The cheap `p=4` shell was exhausted first. The exact rational LDL / Fincke-Pohst streamer then searched the `p=1` shell and found the retained rootless child:

```text
q = 6, (a,b)=(2,3)
A1 pairing = 1
child = rootless/MW17
det = 948
```

Retained vector:

```text
-44717,-282065,63356,564493,-98198,249323,239104,-1054,-22328,-389456,-231271,-641746,-570362,-123785,227276,-186445,89497
```

This rootless hit is what promoted the historical `candidate1` from one of three `A1/MW16` possibilities to the retained suffix.

## 3. The actual generic suffix

The retained search-object lineage is therefore

```text
D7+D5/MW5
  -> Q80-LQ1-D7D4       (q6; old alias: escape)
D7+D4/MW6
  -> Q80-LQ2-A6A4       (q4; old alias: orbit424)
A6+A4/MW7
  -> Q80-LQ3-A6A3       (q4; old alias: orbit1222)
A6+A3/MW8
  -> Q80-LQ4-A4A2A1     (q6; old candidate_id: 7774)
A4+A2+A1/MW10
  -> Q80-LQ5-A3A2       (q4; old shell_id: 1938)
A3+A2/MW12
  -> Q80-LQ6-4A1        (q4; old shell_id: 6855)
4A1/MW13
  -> Q80-LQ7-A1         (q4; old alias: candidate1)
A1/MW16
  -> Q80-LQ8-ROOTLESS   (q6; p=1 shell hit)
rootless/MW17.
```

Machine-readable replay vectors are collected in

```text
data/fibrations/kumar_q80_new_lowq_rootless_path.tsv
data/fibrations/kumar_q80_new_lowq_rootless_geometry.tsv
```

and replayed by

```text
scripts/verify_q80_new_lowq_rootless_geometry.py
scripts/search_q80_new_lowq_final_q6_rootless.py
```

## 4. CM24 equation labels are a later layer

Only after the generic suffix had been selected was each retained divisor transported to CM24 and re-chambered. The special horizontals became

```text
Q80-LQ4-A4A2A1: P3
Q80-LQ5-A3A2:   -P1+P2+2P3
Q80-LQ6-4A1:    2P1
Q80-LQ7-A1:     -P3
Q80-LQ8-ROOTLESS generic divisor: P2-P3 at CM24
```

with large changes in `P.O`, height, fibre twist, and vertical support under specialization.

Therefore `P2-P3` is a **CM24 marking of the already-selected final generic divisor**. It is a regression/identification certificate, not the generic characteristic-zero definition of that section.

This explains why trying to lift the final GF(73) polynomial section directly can become singular or non-transverse without contradicting the generic rootless lattice certificate.

## 5. What the characteristic-zero compiler had lost

The exact characteristic-zero lift through `Q80-LQ7-A1` compiled the q4 pencil, computed its binary-quartic Jacobian, classified the child, and persisted the selected minimal Weierstrass coefficients `A,B`.

What it did **not** persist was a marked MW basis or enough rational-point transport through the quartic-to-Jacobian conversion. That turned the terminal q6 into an apparent fresh section-reconstruction problem even though the generic lattice search had already selected the divisor.

The late repair did not require retraining/reselecting the lattice route. It required reconstructing the lost marking information cheaply enough to identify the already-selected section.

## 6. Exact transport audit

`scripts/trace_q80_candidate1_marked_transport.sage` replays the retained suffix while keeping every unimodular neighbour transport instead of discarding it.

On the selected generic candidate1 frame the terminal divisor chamber-reduces as

```text
D_reduced = O + P - F
```

with exact section class

```text
P = (5,1,-44718,-282065,63356,564493,-98198,249323,239104,-1054,
     -22328,-389456,-231271,-641746,-570362,-123785,227276,-186445,89497)
```

and checks

```text
P^2 = -2
P.F = 1
P.O = 4.
```

The composed candidate1-to-`A6+A3` and candidate1-to-source transforms have determinants `-1` and `+1`, respectively. The same curve becomes a very large multisection in the earlier fibrations. This is the decisive reason not to expect the final section to become a small earlier MW basis vector by raw backward transport.

The transport audit therefore confirms both the search provenance and the coordinate convention while also showing that specialization, not generic back-transport, is where the small `P1,P2,P3` labels arise.

## 7. Successful final-horizontal reconstruction

The direct final section has only three node hits (`I3,I4,I6`). In `P.O=0` degree bounds this leaves a two-parameter x-coordinate and, at the historical prime `73`, the pinned point is singular/non-transverse in those coordinates. The old two-parameter resultant and digit-by-digit `73`-adic routes were therefore computationally poor construction methods.

The successful route is `scripts/recover_q80_final_q6_via_basis_sections.sage`:

1. reconstruct exact `P.O=0` sections that meet four or five reducible nodes, where interpolation is zero- or one-parameter;
2. reduce those exact sections modulo `73` in the transported candidate1 gauge;
3. test their exact elliptic-curve differences against the historical `P2-P3` modular point;
4. compute the matching difference over `QQ(sqrt(-3))(W)` exactly.

This gives the terminal horizontal exactly as an MW-basis difference. The GF73 point is used only to identify the correct pair/sign.

The resulting exact horizontal satisfies

```text
H = P2-P3  (historical CM24 marking)
P.O = 0
height = 1
hits = I3,I4,I6.
```

## 8. Exact RR and child closure

`scripts/certify_q80_final_q6_char0_rr_from_basis.sage` loads that exact horizontal and certifies the resolved terminal pencil over `QQ(sqrt(-3))`:

```text
ambient = 4
whole A4 quotient rank = 1
connected A5 quotient rank = 1
condition rank = 2
kernel dimension = 2
h0(D) = 2.
```

The exact A4 row reduces to the transported historical whole-A4 condition. The A5 row is recovered by a leading-jet toric calculation at the exact `I6` fibre and reduces to the transported `+/-4` quotient line.

`scripts/compile_q80_final_q6_char0_child.sage` then compiles the exact degree-four binary quartic and classifies the characteristic-zero CM24 child as

```text
finite fibres = 4 I3 + I4 + I6 + 2 I1
infinity = smooth
root lattice = 4A2+A3+A5
root rank = 16
root determinant = 1944
Euler number = 24
MW rank = 2.
```

The pinned repository certificate/model are

```text
data/fibrations/q80-final-q6-char0/Q80_CHAR0_FINAL_Q6_CERTIFICATE.md
data/fibrations/q80-final-q6-char0/q80_char0_final_q6_child.sage
```

This closes the characteristic-zero specialization shadow of the generic terminal neighbour. The generic child remains the independently certified rootless `MW17` frame with determinant `948`.

## 9. Superseded diagnostics

The following are retained only as diagnostics/history and should not be described as the active terminal construction:

- the direct two-parameter final-horizontal resultant/Groebner elimination;
- digit-by-digit `73`-adic lifting in the singular residue disk;
- local-73 singularity/tangent probes;
- forcing `P2-P3` as a generic section definition rather than using it as specialization regression data.

## Working rule

**Search provenance first, easy markings second, specialization regression last.**

For this suffix, the generic selected divisor is defined by the exact neighbor-search lineage. The terminal characteristic-zero equation is now reconstructed exactly, and the final-q6 marking gate should be treated as closed.
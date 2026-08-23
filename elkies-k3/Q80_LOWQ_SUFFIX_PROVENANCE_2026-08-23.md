# Q80 low-q suffix provenance and backtracking map — 2026-08-23

## Purpose

This note records how the retained Q80 low-q suffix was actually produced. It is intended to prevent a recurring reconstruction mistake: names such as `7774`, `1938`, `6855`, and `candidate1` are **search-result labels**, not intrinsic names of sections or elliptic fibrations.

The generic lattice search selected a sequence of child frames first. The CM24 horizontal labels (`P3`, `-P1+P2+2P3`, `2P1`, `-P3`, `P2-P3`) were recovered **later**, after transporting those already-selected generic divisors into the CM24 specialization.

This distinction matters for characteristic-zero reconstruction. We should backtrack the selected search objects and their transports rather than repeatedly rediscover the final horizontal from node incidence.

For the repository-wide human-readable route map, including the H3 source route, H2 comparison route, reverse MW17 backtracks, and the role of Q80, see [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md).

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

Therefore `P2-P3` is a **CM24 marking of the already-selected final generic divisor**. It should be used as a regression certificate, not assumed to be a generic characteristic-zero definition of that section.

This explains why trying to lift the final GF(73) polynomial section directly can become singular or non-transverse without contradicting the generic rootless lattice certificate.

## 5. What was lost in the current characteristic-zero equation pipeline

The exact characteristic-zero lift through `Q80-LQ7-A1` currently compiles the q4 pencil, computes its binary-quartic Jacobian, classifies the child, and persists the selected minimal Weierstrass coefficients `A,B`.

However, the compiler does **not** persist a marked MW basis or transport a sufficient set of rational points through the quartic-to-Jacobian conversion. In particular, `compile_q80_q4_candidate1_char0_family.sage` writes the selected parent model but not the new child MW marking.

That is the likely reason the final step has turned into a fresh section-reconstruction problem. The search/lattice chain itself already knows which divisor is intended; the equation conversion discarded the marking needed to identify it cheaply.

## 6. Backtracking plan

The preferred next approach is to recover the marking from the search lineage rather than continue resolving the singular GF(73) section scheme.

1. Reconstruct every retained child frame together with the exact unimodular neighbor transport for

   ```text
   Q80-LQ4-A4A2A1
     -> Q80-LQ5-A3A2
     -> Q80-LQ6-4A1
     -> Q80-LQ7-A1
     -> Q80-LQ8-ROOTLESS.
   ```

2. Compose those transports backward until reaching the earliest equation frame for which explicit marked sections are still available.

3. Express the final q6 horizontal/section class in that earlier marked frame.

4. At equation level, transport enough rational points through each binary-quartic/Jacobian conversion to preserve the required MW direction. If direct point transport is awkward, reconstruct a low-height basis on the child while the parent marking is still available and record the basis change.

5. Form the final horizontal by exact group law on the marked `Q80-LQ7-A1` parent.

6. Use the historical GF(73) `P2-P3`, profile `(0,2,2,0,4)`, A4 row, and A5 residue `+/-4` only as regression checks.

## 7. Current diagnostic status

The direct final-section reconstruction over the selected exact `Q80-LQ7-A1` parent found:

```text
historical GF(73) seed: (c0,c1)=(67,8)
reduced two-parameter Jacobian rank at 73: 1
q6 leading square coefficient: 1 (nonzero)
tangent direction: (17,1)
local transverse elimination obstruction order: 2
```

So the historical point is genuinely non-transverse in the special fibre; this is not merely a denominator-clearing artifact. Digit-by-digit 73-adic lifting is correspondingly expensive and is not the preferred construction route.

## Working rule

**Search provenance first, equation reconstruction second.**

For this suffix, the generic selected divisor is defined by the exact neighbor-search lineage. CM24 supplies a powerful marking/regression scaffold, but it should not be forced to carry the entire characteristic-zero reconstruction when the selected lattice transports can recover the intended object directly.

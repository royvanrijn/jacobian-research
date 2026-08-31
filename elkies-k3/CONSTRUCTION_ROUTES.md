# Elkies rank-17 K3: named construction and comparison routes

This note is the human-readable map of the reconstruction programme. It explains which
paths were discovered by working **backwards** from the recovered rank-17
Mordell--Weil lattice, which polarization is the actual source, which forward corridor
is currently selected, and where the H2/Q80 work fits.

It is a navigation and naming document, not a replacement for the repository status
ledger. `MATH_STATUS.json` remains the status authority; linked notes and exact verifiers
remain the proof/certificate sources. The script-level history and failed-route ledger
are in [`scripts/README.md`](scripts/README.md) and
[`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md).

## The one-picture summary

The project began at the right-hand side of this diagram and worked backwards:

```text
published level-474 Shimura point
          |
          v
principally polarized QM abelian surface
          |
          v
Dolgachev--Kumar K3
          |
          v
H3 Kumar source family: E7+E8 / MW2
          |
          |  selected H3 degree-two corridor
          v
E8+E6/MW3 -> D13/MW4 -> ... -> rootless/MW17
                                      ^
                                      |
                    recovered 17x17 MW lattice
                                      |
                 ---------------------+---------------------
                 |                    |                    |
                 v                    v                    v
        low-q MW2 backtrack    E6 backtrack        other lattice searches
          (q25,4,4,4)          (q90,4,4)

H2 Kumar comparison: E7+E8 / MW2
          |
          | q=80
          v
Q80: E6+D5+A3 / MW3
          |
          | Q80 Low-q Compiler Route
          v
        ... -> rootless/MW17
```

The lines meeting `rootless/MW17` in this picture are not one concatenated
historical path.  The lower arrows out of R17 are reverse-discovery searches;
the upper H3 arrows are the forward source construction selected only after
the source audit.

The crucial distinctions are:

- the **17x17 lattice is the endpoint constraint** from which the reverse searches
  started; its identification with the rootless endpoint of the selected H3
  corridor is now an explicit determinant-one isometry;
- the **H3 Kumar fibration is the actual source polarization** recovered from the
  source geometry;
- the displayed **H3 degree-two corridor is one exact source-to-target path**, not a
  theorem that it is shortest, globally optimal, or easiest to compile;
- **Q80 is a secondary fibration on the same determinant-948 K3 family**, found through
  the H2 comparison polarization. It is extremely useful as an equation/compiler route,
  but it is not the historical/source entrance.

## Target vocabulary: do not use `target` by itself

Several different objects were historically called the target.  They answer
different questions and are not interchangeable.  New notes, scripts, and
artifact fields should use one of the qualified names below.

| qualified name | precise meaning | role in the present programme |
|---|---|---|
| **recovered endpoint** or **R17** | the pinned rootless rank-17 positive frame `data/lattice/rank17_gram.txt` | input constraint for reverse searches; final lattice endpoint of the H3 and generic Q80 routes |
| **H3 source** | the level-474 Kumar `E7+E8/MW2` polarization with height Gram `[[21/2,3],[3,46]]` | geometric starting point selected by the source audit |
| **selected next lattice child** | the already-certified child frame of the retained neighbor | q12/orbit5867 rootless/MW17, integrally marked to pinned R17 |
| **selected next equation child** | the characteristic-zero Weierstrass model constructed for the retained child frame | the exact q12/orbit5867 rootless model, source-identified with geometric Picard rank 19 and full saturated MW lattice R17 |
| **specialization endpoint** | the rational member eventually identified with curve 273 | sought only after the generic rootless H3 family and its sections are explicit |
| **CM24 regression endpoint** | the specialized Q80 terminal model `4A2+A3+A5/MW2` | a compiler regression on a Picard-rank-jumping specialization, not R17 |

For the active H3 equation step, the parent/operation/child triple is therefore

```text
parent equation         = exact P1229-pointed 4A1/MW13 model
selected operation      = q12 degree-two pencil, orbit 5867
selected equation child = rootless/MW17 model marked to pinned R17
route endpoint          = pinned R17
fallback endpoint edge  = q12/orbit4484 (also lattice-certified)
```

The q8 Riemann--Roch compiler has constructed and marked its `4A1/MW13`
child. A complete good-prime polynomial shell showed that the nominal
four-branch cost word is not equation-effective. The corrected exact word
`499+500+69+511-489+933-913` constructs the q12 horizontal, and the direct
`22 -> 2` smooth Riemann--Roch calculation gives the rootless model over QQ.
An exhaustive p=131 polynomial shell and regular Hensel lifts then give
seventeen exact QQ sections whose determinant-948 height Gram is integrally
pinned to R17. An exact point at the old `v=0` I2 support identifies the
quartic with its Jacobian and hence with the pointed H3 forward chain. Exact
counts at good primes 131 and 137 give incompatible rank-20 reduction
discriminant square classes, proving geometric Picard rank 19. The unique
possible index-two enlargement has odd norm 73, so the full geometric MW
lattice is saturated R17 of exact rank 17 and trivial torsion.

<!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->

The reverse endpoint identification and all inverse stage transports are
certified in
[`RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md`](RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md).
That audit also separates the pinned dominant D13 lattice representative from
the component-nef D13 equation representative. They have an exact NS bridge,
but the bridge changes the embedded `U`, so their marked fibers and zeros must
be transported rather than identified by the common `D13/MW4` label.

## Two directions, two questions

The reverse and forward calculations live in the same determinant-948
Neron--Severi isometry class, but they were chosen for different purposes.

| direction | question being answered | selection rule | what it does not identify |
|---|---|---|---|
| **reverse-discovery backtracks from R17** | Can the 17-section frame be replaced by a fibration with fewer sections and more reducible fibres? | prefer exact primitive neighbors with low q, lower MW rank, useful root systems, and transportable divisor classes | the level-474 source polarization |
| **forward H3 Selected Degree-Two Corridor** | Starting from the H3 source selected by the Shimura/Humbert geometry, which certified neighbor chain returns to R17? | keep the pinned H3 marking; among nef degree-two presentations, retain rank-growing children that remain explicitly transportable | a proof that a reverse low-MW endpoint was the source, or that this corridor is optimal or cheapest to compile |

Thus the Low-q MW2 and E6 paths are discovery backtracks.  They helped reveal
the geometry but are not prefixes of the H3 Selected Degree-Two Corridor.
Conversely, the H3 q6/q8/q24 chain is the retained forward corridor even
though R17 was known before H3 was recovered; alternative q24 children,
lateral moves, and higher-q exits may still be cheaper equation routes.

## Canonical names

Use these names in new notes and code comments. Older names remain useful as search
provenance but should not be promoted to object names.

| canonical name | role | starts at | ends at |
|---|---|---|---|
| **H3 Source Family** | recovered source polarization and equation family | level-474 `H21 cap H92` curve | H3 `E7+E8/MW2` |
| **H3 Selected Degree-Two Corridor** | one certified lattice/chamber source-to-target path | H3 `E7+E8/MW2` | rootless `MW17` |
| **H3 Equation Route** | characteristic-zero realization of the selected corridor and its physical suffix | H3 source equation | currently exact through q4/orbit164 `2A3+2A1/MW9` |
| **Q80 Low-q Compiler Route** | secondary independent lattice/compiler route | Q80 `E6+D5+A3/MW3` | rootless `MW17` at generic lattice level |
| **Q80 CM24 Regression Route** | specialization scaffold for Q80 | CM24 specialization of Q80 | CM24 terminal `4A2+A3+A5/MW2` |
| **Low-q MW2 Backtrack** | reverse lattice ancestry discovered from the 17x17 target | rootless `MW17` | `E6+D4+2A2+A1/MW2` |
| **E6 Backtrack** | earlier reverse lattice ancestry | rootless `MW17` | `E6+2A3+2A1/MW3` |
| **H2 Symmetry Comparison** | comparison polarization selected by extra Atkin--Lehner symmetry | H2 `E7+E8/MW2` | Q60/Q80 and other comparison fibrations |
| **H2 Minimal-MW Comparison** | rank-minimizing comparison, not a transported source route | H2 `E7+E8/MW2` | `E8+D7+A1/MW1` |

The source and comparison labels `H2` and `H3` refer to the compatible Kumar height
lattices

```text
H2 = [4       0]
     [0   237/2]

H3 = [21/2   3]
     [   3  46].
```

H2 is the symmetry-selected comparison frame carrying the extra `w2=w237` involution.
The exact `H21 intersect H92` reconstruction identifies H3, not H2, as the source
polarization for the level-474 curve. See
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).

## 1. What we actually did backwards from MW17

The recovered endpoint is the positive-definite rank-17 Mordell--Weil lattice in

```text
data/lattice/rank17_gram.txt
```

with determinant `948`. Because the final fibration is rootless, its Neron--Severi
lattice is represented as

```text
NS = U + (-MW17).
```

For a positive frame `M`, the reverse neighbor searches look for a primitive isotropic
class

```text
f=(a,b,v),
ab=q,
v*M*v^t=2q.
```

Changing the embedded copy of `U` changes the elliptic fibration while keeping the same
K3 surface. This was the key reduction: instead of solving directly for seventeen
sections, search for another fibration on the same surface with a much smaller
Mordell--Weil rank and a larger reducible root system.

### Low-q MW2 Backtrack

This is the cleanest fully transported reverse ancestry found directly from the 17x17
target:

```text
rootless/MW17
 --q25--> A3+7A1/MW7
 --q4 --> D4+A3+2A2+2A1/MW4
 --q4 --> A5+D4+2A2+A1/MW3
 --q4 --> E6+D4+2A2+A1/MW2.
```

The path, exact witnesses, section profiles, and inverse divisor classes are recorded in
[`MW2_FIBRATION_PATH_2026-08-21.md`](MW2_FIBRATION_PATH_2026-08-21.md) and
[`MW2_RANK17_TRANSPORT_2026-08-21.md`](MW2_RANK17_TRANSPORT_2026-08-21.md).

This route answered the original computational question -- *can the 17-section endpoint
be rewritten as a much smaller-section fibration?* -- but it did not identify the true
source construction.

### E6 Backtrack

An earlier reverse route was

```text
rootless/MW17
 --q90--> MW7
 --q4 --> MW4
 --q4 --> E6+2A3+2A1/MW3.
```

It remains an exact Neron--Severi transport and useful provenance certificate; see
[`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md). The old split Weierstrass chart attached
to its terminal Kodaira type was not obtained by geometrically executing this chain, so
that old equation search is not a source reconstruction.

These reverse paths are therefore best called **backtracks**, not constructions.

## 2. How the true beginning was recovered

The source audit changed the direction of the project. The primary-source construction
is

```text
principally polarized QM abelian surface
  -> Dolgachev--Kumar K3 surface
  -> canonical E7+E8 elliptic fibration, MW rank 2
  -> elliptic-neighbor transformations
  -> rootless elliptic fibration, MW rank 17.
```

Exact binary-height-form and discriminant-glue classification leaves three compatible
`E7+E8/MW2` Kumar frames. H2 is distinguished by an extra Atkin--Lehner symmetry, but
the exact Humbert calculation at CM24 shows that the level-474 source lies on the **H3**
polarization, the `H21 intersect H92` component. Its normalization is birational to the
published genus-two model

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

That is why the programme is now being rebuilt **forward from H3**.

## 3. H3 Selected Degree-Two Corridor

Stable stage identifiers use

```text
H3-<sequence>-<child ADE>
```

rather than a search candidate number. They identify the selected corridor; they do not
assert route optimality.

| id | neighbor | child |
|---|---:|---|
| `H3-01-E8E6` | `q6` | `E8+E6/MW3` |
| `H3-02-D13` | `q8` | `D13/MW4` |
| `H3-03-D12` | `q24` | `D12/MW5` |
| `H3-04-A11` | `q6` | `A11/MW6` |
| `H3-05-2A5` | `q8` | `2A5/MW7` |
| `H3-06-3A3` | `q4` | `3A3/MW8` |
| `H3-07-A3-2A2` | `q4` | `A3+2A2/MW10` |
| `H3-08-5A1` | `q4` | `5A1/MW12` |
| `H3-09-4A1` | `q4` | `4A1/MW13` |
| `H3-10-3A1` | `q4` | `3A1/MW14` |
| `H3-11-2A1` | `q4` | `2A1/MW15` |
| `H3-12-A1` | `q4` | `A1/MW16` |
| `H3-13-ROOTLESS` | `q6` | rootless `MW17` |

The selected H3 lattice chain is

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4
 --q24--> D12/MW5
 --q6 --> A11/MW6
 --q8 --> 2A5/MW7
 --q4 --> 3A3/MW8
 --q4 --> A3+2A2/MW10
 --q4 --> 5A1/MW12
 --q4 --> 4A1/MW13
 --q4 --> 3A1/MW14
 --q4 --> 2A1/MW15
 --q4 --> A1/MW16
 --q6 --> rootless/MW17.
```

At the lattice/chamber level this corridor is exact, its rootless endpoint is
explicitly isometric to pinned R17, and every selected arrow is a nef
degree-two pencil. At the characteristic-zero equation level the source family is
explicit and the first two neighbors are exact through `H3-02-D13`. The active selected
equation frontier is `D13/MW4 --q24--> D12/MW5`.

The corridor arose from small-q, old-fibre-degree-two, rank-growing searches. At D13,
the stated proper presentations through q23 were checked and q24 is the first rank-growing
shell, but three q24 orbits give `D12/MW5`; orbit 85 was selected. The later suffix uses
deterministic first hits from that one frame. Therefore alternative q24 children,
lateral moves and larger-q exits remain legitimate candidates for an easier equation
route.

### Why these H3 arrows were retained

The sequence was not selected merely because its ADE labels looked attractive.
The stored search and chamber certificates give the following decision
ledger:

1. The source audit selects H3 from the level-474 `H21 intersect H92`
   geometry.  It is not selected by reversing either low-MW backtrack.
2. From H3, the q6 class gives a nef old-fibre-degree-two
   `E8+E6/MW3` child.  Its q8 continuation is again degree two and gives the
   pinned `D13/MW4` frame.  These are the two equation-level neighbors already
   executed over `QQ`.
3. From D13, every proper presentation through q=20 has MW rank at most four;
   q=21 and q=22 also do not grow, and q=23 has no proper factor presentation.
   The first rank growth is the q24 factorization `(12,2)`.  Exactly three
   primitive orbits give `D12/MW5`; orbit 85 is the retained nef degree-two
   representative.
4. The next retained rank-growing degree-two choices are D12 orbit 42 at q6,
   A11 orbit 922 at q8, and `2A5` orbit 472 at q4.  The compact `3A3` orbit 323
   then jumps to `A3+2A2/MW10` at q4.
5. From MW10, the stored q4 first hits successively remove the remaining root
   rank: orbits 207, 52, 114, 498, and 981 reach `A1/MW16` while preserving
   nefness and old-fibre degree two.
6. At `A1/MW16`, a bounded q4 search found no rootless child.  The streamed q6
   search found exact orbit 2247, whose certified child is R17.  This proves
   the selected final arrow, not an exhaustive minimal-q theorem for the full
   q6 shell.

The detailed counts, witnesses, and chamber arguments remain canonical in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md), and the eleven-step
integral replay remains
[`scripts/verify_h3_d13_to_mw17_path.sage`](scripts/verify_h3_d13_to_mw17_path.sage).

The detailed source and lattice certificates are in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md); the corrected q8 equation
certificate and supersession boundary are summarized in
[`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md). Selection reasons and
failed alternatives are in
[`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md).

## 4. H2 Symmetry Comparison -- why Q80 exists

Before the H3 source marking was identified, H2 was the natural candidate because

```text
H2 = diag(4,237/2)
```

has the extra `w2=w237` symmetry. It is still a valid Kumar polarization on the same
determinant-948 genus and remains extremely useful for comparison and equation
engineering.

The first neighbor forced to use the level-79 direction occurs at `q=60` and gives an
`E8+E6/MW3` comparison frame. A bounded CM-stability search through `q=80` then found
the much better deformation chart

```text
H2 E7+E8/MW2
 --q80--> Q80 = E6+D5+A3/MW3.
```

At the CM24 anchor Q80 gains only one root and keeps all three generic MW directions,
which is why it became the preferred compiler laboratory.

The name **Q80** should mean this specific q80 H2 neighbor / ambient chart; it should not
be used as a name for the source family.

## 5. Q80 Low-q Compiler Route

The common Q80 prefix is

```text
Q80 E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5.
```

From there the retained low-q route is:

| canonical id | old/search alias | neighbor | child |
|---|---|---:|---|
| `Q80-LQ1-D7D4` | `escape` | `q6` | `D7+D4/MW6` |
| `Q80-LQ2-A6A4` | `orbit424` | `q4` | `A6+A4/MW7` |
| `Q80-LQ3-A6A3` | `orbit1222` | `q4` | `A6+A3/MW8` |
| `Q80-LQ4-A4A2A1` | `q6_7774` | `q6` | `A4+A2+A1/MW10` |
| `Q80-LQ5-A3A2` | `q4_1938` | `q4` | `A3+A2/MW12` |
| `Q80-LQ6-4A1` | `q4_6855` | `q4` | `4A1/MW13` |
| `Q80-LQ7-A1` | `candidate1` | `q4` | `A1/MW16` |
| `Q80-LQ8-ROOTLESS` | `final q6` | `q6` | rootless `MW17` |

The old labels are intentionally preserved as aliases because they identify search
artifacts. They are **not** intrinsic names:

- `7774` is a candidate id;
- `1938` and `6855` are shell ids;
- `candidate1` was one of three non-isometric A1 children;
- its later CM24 section label was recovered only after the generic lattice child had
  already been selected.

The complete generic route is therefore

```text
Q80 E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5
 --q6--> D7+D4/MW6
 --q4--> A6+A4/MW7
 --q4--> A6+A3/MW8
 --q6--> A4+A2+A1/MW10
 --q4--> A3+A2/MW12
 --q4--> 4A1/MW13
 --q4--> A1/MW16
 --q6--> rootless/MW17.
```

At the generic lattice level the full route is certified. The terminal generic divisor
is also pinned as an exact section class on the selected candidate1 frame by replaying
the retained neighbour transports.

The characteristic-zero CM24 equation shadow is exact through the terminal q6. The
former marking-loss obstruction was resolved by reconstructing easier high-incidence
`P.O=0` sections on the exact candidate1 parent and identifying the required exact
group-law difference by reduction to the historical `P2-P3` modular point. This avoids
the singular direct three-node Hensel/resultant problem.

The exact terminal resolved RR pencil has ambient dimension 4, condition rank 2, kernel
dimension 2, and `h0(D)=2`. Its exact degree-four binary quartic compiles to the
specialized child `4A2+A3+A5/MW2` with smooth infinity, root rank 16, root determinant
1944, and Euler number 24.

The search provenance and successful marking reconstruction are in
[`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md).
The final proof summary and reproduction commands are in
[`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).

## 6. Q80 CM24 Regression Route

CM24 is a specialization of the Q80 family with larger Picard rank. It is a
regression/compiler scaffold, not the generic construction.

After the `Q80-LQ3-A6A3` stage the specialized route is

```text
2A6+3A1/MW3
 --q6--> A5+2A4+2A1/MW3
 --q4--> 2A4+2A3+A1/MW3
 --q4--> A1+2A3+2D4/MW3
 --q4--> A1+A2+A3+A4+A5/MW3
 --q6--> 4A2+A3+A5/MW2.
```

The historical specialized horizontal labels are

```text
Q80-LQ4-A4A2A1 : P3
Q80-LQ5-A3A2   : -P1+P2+2P3
Q80-LQ6-4A1    : 2P1
Q80-LQ7-A1     : -P3
Q80-LQ8-ROOTLESS generic divisor : P2-P3 at CM24.
```

Those identities are specialization data. In particular `P2-P3` must not be mistaken
for a generic characteristic-zero definition of the final q6 horizontal.

The CM24 equation corridor was first closed over `GF(73)` and is now exact in
characteristic zero over `QQ(sqrt(-3))` through the terminal q6. Its terminal frame is
not rootless because the specialization has extra algebraic classes. See
[`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) and
[`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).

## 7. Other useful comparison paths

### H2 Minimal-MW Comparison

From H2 an exact q4 neighbor gives

```text
H2 E7+E8/MW2 --q4--> E8+D7+A1/MW1.
```

MW rank zero is impossible on this determinant-948 NS, so MW1 is rank-minimal. This is
useful structurally, but it is not connected to the rootless 17x17 target by a pinned
explicit neighbor transport and its remaining section has pole order 58. It is therefore
a comparison endpoint, not the active construction route. See
[`MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md`](MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md).

### Q60 Compact Comparison

The H2 q60 `E8+E6/MW3` frame is the first compact neighbor forced to use the level-79
direction. It was important in understanding pole reduction and CM boundaries, but Q80
is the more CM-stable deformation chart and H3 is the actual source polarization. Call
this the **Q60 Compact Comparison**, not the source route.

## 8. What is complete, and what is not

| route/family | lattice / chamber | characteristic-zero equations | intended use |
|---|---|---|---|
| H3 Source Family | source polarization and marked family exact | exact `E7+E8/MW2` source | canonical beginning |
| H3 Selected Degree-Two Corridor | complete to pinned R17, with lossless inverse NS transport back to H3 | exact through q8/orbit376 and q12/orbit5867, including source identity and full saturated R17 | completed construction proof; arithmetic base-change search is now primary |
| Q80 Low-q Compiler Route | complete to rootless/MW17 | terminal CM24 shadow exact; generic rootless endpoint remains a lattice certificate | secondary construction/compiler |
| Q80 CM24 Regression Route | specialization lattice/equation corridor complete | complete through `4A2+A3+A5/MW2` | regression and module discovery |
| Low-q MW2 Backtrack | exact transport from MW17 to MW2 | endpoint rational model exists, but this is not the source route | reverse provenance / alternate reconstruction |
| E6 Backtrack | exact transport from MW17 to MW3 | old guessed split chart not a transported construction | historical reverse provenance |
| H2 Minimal-MW Comparison | exact H2->MW1 neighbor | not a target construction | structural comparison |
<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 2c6a2a36699933ab -->

## 9. Recommended language in future notes

Use:

- **rootless MW17 target** or `R17` for the recovered 17x17 endpoint;
- **H3 Source Family** for the true Kumar-source polarization/equation family;
- **H3 Selected Degree-Two Corridor** for the currently certified lattice path;
- **H3 Equation Route** for the characteristic-zero realization, currently exact
  through the C8-pointed q4/orbit164 child;
- **H2 Symmetry Comparison** for the `diag(4,237/2)` polarization;
- **Q80 Low-q Compiler Route** for the secondary rootless lattice path;
- **Q80 CM24 Regression Route** for its specialized equation scaffold;
- **Low-q MW2 Backtrack** and **E6 Backtrack** for reverse-discovered lattice ancestries.

For individual Q80 suffix stages, prefer `Q80-LQ*-<child>` identifiers and put `7774`,
`1938`, `6855`, `candidate1`, etc. in parentheses as historical aliases only.

The conceptual rule is:

> **Reverse search found useful alternative fibrations on the R17 surface;
> the independent source audit selected the H3 Source Family; the H3 Selected
> Degree-Two Corridor proves one forward route, not the shortest or cheapest
> one. Q80 is a separate compiler route, and its CM24 endpoint is a
> specialization regression rather than R17.**

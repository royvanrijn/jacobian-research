# Elkies rank-17 K3: canonical construction routes

This note is the human-readable map of the reconstruction programme.  It
explains which paths were discovered by working **backwards** from the recovered
rank-17 Mordell--Weil lattice, which path is the actual source construction,
and where the H2/Q80 work fits.

It is a navigation and naming document, not a replacement for the repository
status ledger.  `MATH_STATUS.json` remains the status authority; the linked
notes and exact verifiers remain the proof/certificate sources.

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
H3 Kumar source: E7+E8 / MW2
          |
          |  H3 Source Route  (PRIMARY)
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
          | Q80 Low-q Compiler Route  (SECONDARY)
          v
        ... -> rootless/MW17
```

The crucial distinction is:

- the **17x17 lattice is the endpoint constraint** from which the reverse
  searches started;
- the **H3 Kumar fibration is the actual source polarization** recovered from
  the source geometry;
- **Q80 is a secondary fibration on the same determinant-948 K3 family**, found
  through the H2 comparison polarization.  It is extremely useful as an
  equation/compiler route, but it is not the historical/source entrance.

## Canonical route names

Use these names in new notes and code comments.  Older names remain useful as
search provenance but should not be promoted to object names.

| canonical name | role | starts at | ends at |
|---|---|---|---|
| **H3 Source Route** | primary source-to-target construction | H3 `E7+E8/MW2` | rootless `MW17` |
| **Q80 Low-q Compiler Route** | secondary independent equation/compiler route | Q80 `E6+D5+A3/MW3` | rootless `MW17` |
| **Q80 CM24 Regression Route** | specialization scaffold for the Q80 route | CM24 specialization of Q80 | CM24 terminal `4A2+A3+A5/MW2` |
| **Low-q MW2 Backtrack** | reverse lattice ancestry discovered from the 17x17 target | rootless `MW17` | `E6+D4+2A2+A1/MW2` |
| **E6 Backtrack** | earlier reverse lattice ancestry | rootless `MW17` | `E6+2A3+2A1/MW3` |
| **H2 Symmetry Comparison Route** | comparison polarization selected by the extra Atkin--Lehner symmetry | H2 `E7+E8/MW2` | Q60/Q80 and other comparison fibrations |
| **H2 Minimal-MW Comparison** | rank-minimizing comparison, not a transported source route | H2 `E7+E8/MW2` | `E8+D7+A1/MW1` |

The source and comparison labels `H2` and `H3` refer to the compatible Kumar
height lattices

```text
H2 = [4       0]
     [0   237/2]

H3 = [21/2   3]
     [   3  46].
```

H2 is the symmetry-selected comparison frame carrying the extra `w2=w237`
involution.  The exact `H21 intersect H92` reconstruction identifies H3, not
H2, as the source polarization for the level-474 curve.  See
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).

## 1. What we actually did backwards from MW17

The recovered endpoint is the positive-definite rank-17 Mordell--Weil lattice
in

```text
data/lattice/rank17_gram.txt
```

with determinant `948`.  Because the final fibration is rootless, its
Neron--Severi lattice is represented as

```text
NS = U + (-MW17).
```

For a positive frame `M`, the reverse neighbor searches look for a primitive
isotropic class

```text
f=(a,b,v),
ab=q,
v*M*v^t=2q.
```

Changing the embedded copy of `U` changes the elliptic fibration while keeping
the same K3 surface.  This was the key reduction: instead of solving directly
for seventeen sections, search for another fibration on the same surface with
a much smaller Mordell--Weil rank and a larger reducible root system.

### Low-q MW2 Backtrack

This is the cleanest fully transported reverse ancestry found directly from
the 17x17 target:

```text
rootless/MW17
 --q25--> A3+7A1/MW7
 --q4 --> D4+A3+2A2+2A1/MW4
 --q4 --> A5+D4+2A2+A1/MW3
 --q4 --> E6+D4+2A2+A1/MW2.
```

The path, exact witnesses, section profiles, and inverse divisor classes are
recorded in
[`MW2_FIBRATION_PATH_2026-08-21.md`](MW2_FIBRATION_PATH_2026-08-21.md) and
[`MW2_RANK17_TRANSPORT_2026-08-21.md`](MW2_RANK17_TRANSPORT_2026-08-21.md).

This route answered the original computational question -- *can the 17-section
endpoint be rewritten as a much smaller-section fibration?* -- but it did not
identify the true source construction.

### E6 Backtrack

An earlier reverse route was

```text
rootless/MW17
 --q90--> MW7
 --q4 --> MW4
 --q4 --> E6+2A3+2A1/MW3.
```

It remains an exact Neron--Severi transport and useful provenance certificate;
see [`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).  The old split Weierstrass
chart attached to its terminal Kodaira type was not obtained by geometrically
executing this chain, so that old equation search is not a source
reconstruction.

These reverse paths are therefore best called **backtracks**, not
constructions.

## 2. How the true beginning was recovered

The source audit changed the direction of the project.  The primary-source
construction is

```text
principally polarized QM abelian surface
  -> Dolgachev--Kumar K3 surface
  -> canonical E7+E8 elliptic fibration, MW rank 2
  -> elliptic-neighbor transformations
  -> rootless elliptic fibration, MW rank 17.
```

Exact binary-height-form and discriminant-glue classification leaves three
compatible `E7+E8/MW2` Kumar frames.  H2 is distinguished by an extra
Atkin--Lehner symmetry, but the exact Humbert calculation at CM24 shows that
the level-474 source lies on the **H3** polarization, the `H21 intersect H92`
component.  Its normalization is birational to the published genus-two model

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

That is why the programme is now being rebuilt **forward from H3**.

## 3. H3 Source Route -- the canonical construction

Canonical stage identifiers use

```text
H3-<sequence>-<child ADE>
```

rather than a search candidate number.

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

So the canonical H3 lattice chain is

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

At the lattice/chamber level this route is exact and every selected arrow is a
nef degree-two pencil.  At the characteristic-zero equation level the source
family is explicit and the first two neighbors are exact through
`H3-02-D13`.  The active equation frontier is therefore **D13/MW4 forward**,
not the recovery of the source or of the first q8 child.

The detailed source and lattice certificates are in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md); the corrected q8 equation
certificate and supersession boundary are summarized in
[`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

## 4. H2 Symmetry Comparison Route -- why Q80 exists

Before the H3 source marking was identified, H2 was the natural candidate
because

```text
H2 = diag(4,237/2)
```

has the extra `w2=w237` symmetry.  It is still a valid Kumar polarization on
the same determinant-948 genus and remains extremely useful for comparison
and equation engineering.

The first neighbor forced to use the level-79 direction occurs at `q=60` and
gives an `E8+E6/MW3` comparison frame.  A bounded CM-stability search through
`q=80` then found the much better deformation chart

```text
H2 E7+E8/MW2
 --q80--> Q80 = E6+D5+A3/MW3.
```

At the CM24 anchor Q80 gains only one root and keeps all three generic MW
directions, which is why it became the preferred compiler laboratory.

The name **Q80** should mean this specific `q=80` H2 neighbor / ambient chart;
it should not be used as a name for the source family.

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

The old labels are intentionally preserved as aliases because they identify
search artifacts.  They are **not** intrinsic names:

- `7774` is a candidate id;
- `1938` and `6855` are shell ids;
- `candidate1` was one of three non-isometric A1 children;
- its later CM24 section label was recovered only after the generic lattice
  child had already been selected.

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

At the generic lattice level the full route is certified.  The
characteristic-zero equation pipeline has progressed through the
`Q80-LQ7-A1` child: it can compile the penultimate q4 pencil, form and classify
its binary-quartic Jacobian, and persist the A1/MW16 parent for the final
neighbor.  The current final-q6 obstruction is **marking loss**, not discovery
of the lattice path: the equation conversion did not persist enough MW point
transport to identify the required final horizontal cheaply.

The exact search provenance and current backtracking repair plan are in
[`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md).

## 6. Q80 CM24 Regression Route

CM24 is a specialization of the Q80 family with larger Picard rank.  It is a
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

Those identities are specialization data.  In particular `P2-P3` must not be
mistaken for a generic characteristic-zero definition of the final q6
horizontal.

The CM24 equation corridor is exact through the final q6 over `GF(73)`.  Its
terminal frame is not rootless because the specialization has extra algebraic
classes.  See [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md).

## 7. Other useful comparison paths

### H2 Minimal-MW Comparison

From H2 an exact `q=4` neighbor gives

```text
H2 E7+E8/MW2 --q4--> E8+D7+A1/MW1.
```

MW rank zero is impossible on this determinant-948 NS, so MW1 is
rank-minimal.  This is useful structurally, but it is not connected to the
rootless 17x17 target by a pinned explicit neighbor transport and its remaining
section has pole order 58.  It is therefore a comparison endpoint, not the
active construction route.  See
[`MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md`](MW1_OPTIMAL_FIBRATION_PATH_2026-08-21.md).

### Q60 comparison

The H2 `q=60` `E8+E6/MW3` frame is the first compact neighbor forced to use
the level-79 direction.  It was important in understanding pole reduction and
CM boundaries, but Q80 is the more CM-stable deformation chart and H3 is the
actual source polarization.  Call this the **Q60 Compact Comparison**, not the
source route.

## 8. What is complete, and what is not

| route | lattice / chamber | characteristic-zero equations | intended use |
|---|---|---|---|
| H3 Source Route | complete to rootless/MW17 | exact through D13/MW4 | primary construction |
| Q80 Low-q Compiler Route | complete to rootless/MW17 | exact pipeline through A1/MW16 parent; final q6 marking still to preserve/recover | secondary construction/compiler |
| Q80 CM24 Regression Route | specialization lattice/equation corridor complete | complete over `GF(73)` | regression and module discovery |
| Low-q MW2 Backtrack | exact transport from MW17 to MW2 | endpoint rational model exists, but this is not the source route | reverse provenance / alternate reconstruction |
| E6 Backtrack | exact transport from MW17 to MW3 | old guessed split chart not a transported construction | historical reverse provenance |
| H2 Minimal-MW Comparison | exact H2->MW1 neighbor | not a target construction | structural comparison |

## 9. Recommended language in future notes

Use:

- **rootless MW17 target** or `R17` for the recovered 17x17 endpoint;
- **H3 Source Route** for the true Kumar-source construction;
- **H2 Symmetry Comparison** for the `diag(4,237/2)` polarization;
- **Q80 Low-q Compiler Route** for the secondary rootless path;
- **Q80 CM24 Regression Route** for its specialized equation scaffold;
- **Low-q MW2 Backtrack** and **E6 Backtrack** for reverse-discovered lattice
  ancestries.

For individual Q80 suffix stages, prefer `Q80-LQ*-<child>` identifiers and put
`7774`, `1938`, `6855`, `candidate1`, etc. in parentheses as historical
aliases only.

The conceptual rule is simple:

> **Backtracking discovered the geometry; source recovery selected H3; forward
> equation reconstruction is now replaying the H3 route, with Q80 as an
> independent compiler path to the same rootless target.**

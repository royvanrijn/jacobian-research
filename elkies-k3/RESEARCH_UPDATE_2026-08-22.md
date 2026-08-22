# Research update — 2026-08-22

## H3 q8 re-audit changes the active equation priority

The first H3 neighbour remains exact:

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3.
```

Its complete resolved RR cover, two-dimensional q6 pencil, minimal child equation, fibre classification, MW height lattice, and Neron--Severi discriminant check remain re-authorized.

The second q8 equation-level hop is **paused**. A new regression audit proves that the exact rational child section used by `derive_h92_q6_child_q8_marking.sage` cannot have the MW coordinate assigned to it.

The claimed coordinate is

```text
(-2,-2,0),
```

which has height `24` in the certified child height Gram

```text
[[8/3,1/3,-1],
 [1/3,8/3,1],
 [-1,1,46]].
```

But the exact rational functions for that section meet the standard zero section transversely over a squarefree degree-46 divisor of smooth fibres. Hence `S.O=46`. Shioda's self-height formula, even using an intentionally loose maximal `E8+E6` correction of `36`, forces height at least `60`.

Reproduce the contradiction with

```bash
sage -python elkies-k3/scripts/audit_h92_q6_child_q8_marking_height.sage
```

The authoritative H3 q8 checkpoint is now [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

### What remains exact on H3 q8

The q8 **lattice** classification remains valid. A second audit also resolves the degree-18/degree-16 representative confusion:

```text
dominant D13 representative
  --122 finite old-fibre root reflections-->
degree-18 historical ambient class,

classifier-nef representative
  --120 finite old-fibre root reflections-->
degree-16 chamber class.
```

The degree-16 terminal vector is

```text
(22,16,-14,-20,-27,-40,-33,-26,-18,-4,-5,-7,-10,-8,-6,-4,-2,8,0),
```

exactly the class independently found by the later horizontal/fibre chamber reduction. This re-authorizes the degree-16 **lattice representative**, but not the hand-translated `q6^8` local compiler built afterward.

Reproduce with

```bash
sage -python elkies-k3/scripts/audit_h92_q8_representative_selection.sage
```

### Retracted/conditional q8 equation work

Until the point-to-MW bridge is repaired, treat as conditional diagnostics rather than q8 certificates:

- the q8 child marking built from the rational section `S`;
- child-side component-nef/nef/dominant chord, finite, q-frame, infinity, global-intersection and branch constructions that depend on that marking;
- the historical degree-18 `true1600` pipeline as the RR system of the final q8 moving divisor;
- the experimental degree-16 `corrected1278` q6^8 local compiler and its modular survivor sequence.

The next H3 task is not another local rank probe. It is to compute the actual canonical heights and pairwise height pairings of the exact child points `old_zero`, `affine_E7`, `E7_7`, their two differences, and the selected group combination, then match those rational points to the pinned MW lattice independently of the NS labels.

## Q80 is now the live equation-construction route

The Q80 route is independent of the H3 marking failure. Its generic lattice corridor remains certified all the way to a new rootless `MW17` frame, and its equation-level CM24 corridor has advanced through q6_7774 and q4_1938, with later low-q stages continuing toward the rootless endpoint.

The retained generic continuation is

```text
D7+D5/MW5
 --q6 (2,3)--> D7+D4/MW6
 --q4 (2,2)--> A6+A4/MW7
 --q4 (2,2)--> A6+A3/MW8
 --q6 (2,3)--> A4+A2+A1/MW10
 --q4 (2,2)--> A3+A2/MW12
 --q4 (2,2)--> 4A1/MW13
 --q4 (2,2)--> A1/MW16
 --q6 (2,3)--> rootless/MW17.
```

All retained new divisors chamber-reduce to old-fibre degree `2`.

Detailed Q80 status remains in:

- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md)
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv)
- [`Q80_LOW_Q_ALTERNATE_2026-08-22.md`](Q80_LOW_Q_ALTERNATE_2026-08-22.md)
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md)

## Reusable compiler lesson

The Q80 equation work reinforces a rule that also matters when H3 q8 resumes: a connected vertical ADE correction must be compiled as a **single resolved quotient-line/module condition**. Treating listed exceptional components as independent evaluation rows can create spurious rank and false RR obstructions.

## Current execution order

1. Continue the Q80 equation corridor toward the rootless endpoint.
2. In parallel, repair the H3 q8 rational-point/MW bridge by independent height and component calculations.
3. Only after that bridge is certified, rebuild the H3 q8 local modules from the corrected rational section or derive the degree-16 source-chamber modules directly.
4. Do not spend more time on `true1600`, `corrected1278`, or additional H3 q8 local jet screens in their current form.

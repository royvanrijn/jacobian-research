# H3 q=8 current frontier

Status: 2026-08-22, after the point-to-MW marking re-audit.

The active checkpoint is [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md). It supersedes the active q8 conclusions in the earlier module-intersection notes and in the q8 sections of `BISECTION_COLLISION_SEARCH.md` wherever they conflict.

## Objective

The lattice target is unchanged:

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

The first q6 hop is still exact. The second q8 equation-level hop is **not currently certified**.

## Re-authorized exact facts

The following remain valid:

- the exact first q6 neighbour and its complete resolved RR cover;
- the globally minimal q6 child with `E8+E6/MW3`;
- the child MW height Gram

  ```text
  [[8/3,1/3,-1],
   [1/3,8/3,1],
   [-1,1,46]];
  ```

- the q8 lattice orbit classification and abstract D13/MW4 target;
- the independently chamber-certified q8 classifier-nef representative;
- the source-side representative reconciliation:

  ```text
  dominant D13 representative --122 old-fibre root reflections--> degree-18 class
  classifier-nef representative --120 old-fibre root reflections--> degree-16 class
  ```

- the degree-16 terminal vector

  ```text
  (22,16,-14,-20,-27,-40,-33,-26,-18,-4,-5,-7,-10,-8,-6,-4,-2,8,0),
  ```

  which is exactly the later experimental chamber-reduced class;
- the exact rational child points and section formulas as rational points on the q6 child, without the withdrawn MW/NS identification;
- the independent Q80 programme.

Reproduce the representative reconciliation with

```bash
sage -python elkies-k3/scripts/audit_h92_q8_representative_selection.sage
```

## Confirmed invalid bridge

`derive_h92_q6_child_q8_marking.sage` assigns the constructed rational section `S` the MW coordinate

```text
(-2,-2,0),
```

which has height `24` in the certified child MW lattice.

But the exact rational functions for `S` have a square/cube denominator root `h` of degree `46`, coprime to the discriminant, and the zero-section parameter `z=-x/y` has a simple zero at every root of `h`. Hence

```text
S.O = 46
```

on smooth fibres.

Shioda's height formula on the elliptic K3 gives

```text
<S,S> = 4 + 2(S.O) - sum_v contr_v(S).
```

Even using the deliberately loose bound

```text
sum_v contr_v(S) <= maxdiag(E8^-1)+maxdiag(E6^-1) = 30+6 = 36,
```

the rational section must have height at least

```text
4 + 2*46 - 36 = 60,
```

contradicting `24`.

Reproduce this with

```bash
sage -python elkies-k3/scripts/audit_h92_q6_child_q8_marking_height.sage
```

Expected status:

```text
PASS_CONFIRMED_MARKING_BRIDGE_CONTRADICTION
```

## Retracted / conditional q8 equation work

Until the rational point-to-MW bridge is repaired, do **not** treat the following as exact q8-neighbour results:

- the child q8 marking `relative_child_section_MW_coordinates=(-2,-2,0)`;
- child-side component-nef/nef/dominant q8 chord, collision, finite, q-frame, fractional, infinity, global-intersection, or branch conclusions that depend on that marking;
- the degree-18 `true1600` source pipeline as the RR system of the final q8 moving divisor;
- the hand-converted degree-16 `corrected1278` q6^8 local compiler and its `1278 -> 14 -> 7` survivor sequence.

The modular ranks from those experiments remain useful diagnostics for the declared inputs, but they are no longer geometry certificates for the target q8 pencil.

## What the recent degree-16 work did establish

The degree-18 class historically called `source-nef` was not the classifier-nef representative. A bounded effective-root audit found a degree-two `(-2)` wall with negative pairing, and complete alternating reduction reached the exact degree-16 vector above. The representative audit independently confirms that this vector is the finite-root reduction of the classifier's nef representative.

So retain the **degree-16 lattice chamber result**, but discard the inference

```text
degree 16 => replace every q6^9 local module by a q6^8 module.
```

That local-module step was never proved and produced contradictory downstream behaviour.

## Next exact gate

Do not run more q8 local rank probes yet.

First repair the child point-to-MW bridge. On the final minimal q6 child, independently compute the heights, pairings, zero intersections, and reducible-fibre component corrections of

```text
old_zero,
affine_E7 point,
E7_7 point,
E7_7-old_zero,
E7_7-affine_E7,
2*(E7_7-old_zero)+2*(E7_7-affine_E7).
```

Then match their actual height Gram to the pinned rank-three MW lattice. Only after that match should the q8 lattice coordinate `(0,-2,0)` be converted into an actual rational child section.

If that bridge is repaired, rebuild the q8 local modules from the corrected section. If it is not, derive the degree-16 source-chamber local modules directly from the resolved divisor, with no q6-power shortcut.

## Execution priority

The H3 lattice route remains valid and attractive, but its q8 equation hop is paused at this marking repair. The Q80 route is therefore the live equation-construction path for now.

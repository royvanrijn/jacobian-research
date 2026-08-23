# H3 q24 D12 suffix ledger

Date: 2026-08-23

Status: **q24 equation-side marking is pinned to the R17-directed lattice corridor; the
fresh first-hit D12 -> A1 replay is retained as an alternate branch.**

## 1. What was originally found

Commit `bfb96bf526b0b289e12d531e7aff697cae5f39c6` added a fresh replay starting from a
D12 frame exported from the equation-side q24 divisor `D24eq`:

```text
scripts/export_h92_q24_native_d12_frame.sage
scripts/run_h92_q24_native_suffix_to_a1.py
```

That replay found the valid degree-two sequence

```text
D12
 --q6 orbit 55--> A11
 --q8 orbit 1960--> 2A5
 --q4 orbit 261--> 3A3
 --q4 orbit 567--> A3+2A2
 --q4 orbit 3822--> 5A1
 --q4 orbit 8598--> 4A1
 --q4 orbit 13487--> 3A1
 --q4 orbit 18368--> 2A1
 --q4 orbit 17593--> A1.
```

Its terminal status remains a correct statement:

```text
PASS_Q24_NATIVE_D12_TO_A1
```

It proves that this D12 coordinate frame has a native degree-two path to an A1/MW16
fibration.  It does **not** prove that this is the A1 lying on the selected route to the
recovered R17 endpoint.

## 2. Why that branch is not the preferred R17 suffix

Two exact diagnostics exposed the distinction.

First, the fresh A1 admits a q6 rootless child (orbit 504), but that rootless
determinant-948 frame is not integrally isometric to pinned R17.  Exhausting the initial
10,000-vector capped q6 sample found no pinned-R17 child.

Second, pulling the pinned R17 fibration itself into that fresh A1 coordinates gives an
enormous old-fibre degree rather than degree two.  Therefore that A1 is a genuinely
different elliptic fibration on the same K3 surface.

This is a useful warning: equal `ADE/MW` labels do not determine the embedded elliptic
`U`, and repeated deterministic first-hit searches can leave a previously selected
marked corridor while preserving the same sequence of root ranks.

## 3. The marking issue at D13

A temporary diagnostic compared the raw H3-coordinate q24 fibre rays and appeared to
show divergence already at D12.  That interpretation was too strong.

The selected historical/pinned D13 search frame and the current equation D13 frame are
related by a nontrivial exact NS marking change.  That change does not preserve the
literal embedded standard `U`.

The current equation-frame transport stack is authoritative:

```text
pinned/dominant D13
  -> nef/component marking
  -> source-H3 ambient
  -> physical component reflections
  -> q6 Weyl transport
  -> q6 Eichler translation
  -> current equation D13.
```

`compare_h92_three_q24_d12_current_equation_profiles.sage` independently transports all
three q24 D12 candidates through this stack.  Orbit 85 reproduces the pinned equation
profile

```text
MW=(-2,1,-1,1)
height=52
correction=0
P.O=24
vertical_F=-7
root_L1=69
support=13
```

so orbit 85 remains the selected q24 child in the current equation marking.

## 4. Exact R17-directed closeout

The authoritative closeout is now

```text
scripts/certify_h92_q24_equation_d13_to_pinned_r17.sage
```

with terminal status

```text
PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH
```

It reconstructs the exact 19x19 coordinate map from the current equation D13 marking to
the historical selected D13 marking, transports q24 orbit 85 into current equation
coordinates, and deliberately chooses the resulting D12 basis as the already-certified
R17-directed D12 frame.

From there the exact suffix is:

| step | q | orbit | child | MW |
|---:|---:|---:|---|---:|
| 1 | 6 | 42 | `A11` | 6 |
| 2 | 8 | 922 | `2A5` | 7 |
| 3 | 4 | 472 | `3A3` | 8 |
| 4 | 4 | 323 | `A3+2A2` | 10 |
| 5 | 4 | 207 | `5A1` | 12 |
| 6 | 4 | 52 | `4A1` | 13 |
| 7 | 4 | 114 | `3A1` | 14 |
| 8 | 4 | 498 | `2A1` | 15 |
| 9 | 4 | 981 | `A1` | 16 |
| 10 | 6 | 2247 | rootless | 17 |

The final q6 is therefore explicitly

```text
A1/MW16 --q6 orbit 2247, old-fibre degree 2--> rootless/MW17.
```

The resulting rootless frame is then identified with
`data/lattice/rank17_gram.txt` by the stored determinant-one endpoint isometry.

Successful terminal lines include:

```text
Q24R17MAP_FINAL_Q6|parent=A1/MW16|q=6|orbit=2247|degree=2|child=rootless/MW17|status=PASS
Q24R17MAP_PINNED|...|pinned_R17=1|status=PASS
Q24R17MAP_RESULT|...|status=PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH
```

## 5. Proof boundary

This closes the route at the **integral lattice / Neron--Severi / marking** level from
the current equation D13 frame all the way to pinned R17.

It does not execute the downstream characteristic-zero genus-one pencils and
Weierstrass/Jacobian transformations.

The current equation programme is therefore:

```text
exact D13 equation
  -- q24 orbit85 -->
D12 equation
  -- q6 -->
A11
  -- q8 -->
2A5
  -- q4 ... -->
A1
  -- q6 orbit2247 -->
rootless R17 equation.
```

The pointed q24 D12 work already recovers useful modular data for the next q6 step:
`recover_h92_q24_pointed_zero_pole_sections.sage` reconstructs an explicit A11 target
section modulo a good prime as a group-law combination of easy `P.O=0` sections.
Its own proof boundary states that the resolved q6 RR pencil and A11 child equation
remain to be compiled.

## 6. Preferred-route consequence

The preferred route is the transported R17-directed suffix described in
`H3_PREFERRED_PATH.md`.

The fresh first-hit branch
`55,1960,261,567,3822,8598,13487,18368,17593` remains valuable as an alternate
fibration/compiler experiment, but it must not be called the canonical R17 suffix.

q32 is likewise retained as alternate/regression work and as a source of reusable
pointed-quartic/spinor-marking techniques.

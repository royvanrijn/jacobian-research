# H3 preferred construction path

Status: **preferred route and proof-boundary note**.

This file records the route that should be treated as the default H3 construction path
after the exact equation-D13-to-pinned-R17 lattice/marking closeout on 2026-08-23.

The crucial distinction is now explicit:

- the **current equation D13 marking** is not the same literal embedded `U` as the older
  pinned/dominant D13 search frame;
- q24 orbit 85 nevertheless transports exactly into the current equation D13 frame;
- choosing the q24 child basis by that exact marking transport gives the
  **R17-directed D12 frame**;
- the certified suffix then replays literally through the final q6 orbit 2247 to the
  original pinned 17x17 R17 lattice.

The fresh first-hit suffix produced directly from a separately minimized `D24eq` D12
frame (orbits `55,1960,261,567,3822,8598,13487,18368,17593`) is a valid alternate
lattice branch, but it is **not** the preferred R17-directed suffix.

## Preferred route

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4                         [current exact equation D13 marking]
 --q24 orbit 85 --> D12/MW5              [transported R17-directed D12 basis]
 --q6 orbit 42 --> A11/MW6
 --q8 orbit 922 --> 2A5/MW7
 --q4 orbit 472 --> 3A3/MW8
 --q4 orbit 323 --> A3+2A2/MW10
 --q4 orbit 207 --> 5A1/MW12
 --q4 orbit 52 --> 4A1/MW13
 --q4 orbit 114 --> 3A1/MW14
 --q4 orbit 498 --> 2A1/MW15
 --q4 orbit 981 --> A1/MW16
 --q6 orbit 2247 --> rootless/MW17 = pinned R17
```

All selected lattice arrows have old-fibre degree two.  The exact closeout script

```text
scripts/certify_h92_q24_equation_d13_to_pinned_r17.sage
```

reconstructs the full 19-dimensional marking change from the current equation D13 frame
to the historical selected D13 frame, transports q24 orbit 85 into equation coordinates,
replays the complete suffix, and verifies a determinant-one NS transport to
`data/lattice/rank17_gram.txt`.

Expected terminal status:

```text
PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH
```

## What is now closed

At the **lattice / Neron--Severi / marking** level, the route from the current exact
equation D13 frame to the original recovered R17 endpoint is now closed.

In particular the previously uncertain final step is explicit:

```text
A1/MW16 --q6 orbit 2247, old-fibre degree 2--> rootless/MW17 = pinned R17.
```

The successful closeout also resolves the apparent contradiction from the fresh
`D24eq -> ... -> A1` first-hit replay.  That replay changed coordinate/minimization
choices and then selected new first hits at every stage.  Matching ADE/MW labels was
not enough to keep it on the R17-directed marked corridor.

## Proof boundary

Do not promote this lattice/marking closeout to an equation-level claim.

Current equation status is:

1. **H3 source -> E8+E6 -> D13:** exact characteristic-zero equations and markings.
2. **D13 -> D12 q24:** the selected orbit-85 divisor and its equation-frame marking are
   pinned exactly, and modular/pointed D12 compiler work is available; the complete
   characteristic-zero D12 Weierstrass/Jacobian + marked-section certificate still
   needs to be finalized unless a later dedicated QQ certificate supersedes this note.
3. **D12 -> R17 suffix:** exact integral neighbor/marking transport is fully pinned, but
   the downstream characteristic-zero Weierstrass pencils have not all been executed.
4. The pointed D12 work has already recovered an explicit modular A11 target section
   from easy zero-pole sections; the next concrete compiler target is the q6
   `D12 -> A11` resolved RR pencil and child equation.

## Route priority

Use this order unless new evidence changes the equation cost materially:

1. finish the characteristic-zero q24 `D13 -> D12` child model in the pinned
   R17-directed marking;
2. compile the q6 `D12 -> A11` equation using the pointed/zero-pole section machinery;
3. continue the exact transported suffix
   `q8, q4, q4, q4, q4, q4, q4, q4`;
4. execute the final q6 orbit 2247 at equation level and verify the resulting generic
   rootless family carries the recovered rank-17 MW lattice;
5. only then perform the specialization/identification with curve 273.

q32 remains useful as an alternate compiler and marking-recovery laboratory, but it is
not the preferred H3 route.

## Alternate first-hit branch

The script

```text
scripts/run_h92_q24_native_suffix_to_a1.py
```

found the valid first-hit sequence

```text
55, 1960, 261, 567, 3822, 8598, 13487, 18368, 17593
```

from its minimized q24 D12 frame.  Its endpoint A1 is not the R17-directed A1; a direct
pullback of pinned R17 into that A1 has enormous degree.  Keep this branch as a useful
alternate lattice/compiler experiment, not as the canonical R17 suffix.

For the closeout details see `H3_Q24_R17_CLOSEOUT_2026-08-23.md`.  For broader route
provenance see `CONSTRUCTION_ROUTES.md` and `SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`.

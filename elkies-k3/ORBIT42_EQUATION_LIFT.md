# H3 orbit42 equation lift

This note records the exact equation-level boundary for the selected

```text
D12/MW5 --q6 orbit42--> A11/MW6
```

continuation.  It is the canonical note for the orbit42 milestones below.
The complete lattice/chamber corridor to the pinned rootless/MW17 frame is
already exact; the equation pencils after D12 are not.

## Exact parent and orbit42 profile

The parent is the characteristic-zero q24/orbit85 child certified by
`EC-K3-H3-Q24-QQ-D12`.  In its current equation frame the selected orbit42
class has

```text
mw = (-1,0,-1,-1,0)
height = 7
local correction = 3
P.O = 3
fibre twist = 0
D42 = O + P + V.
```

The correction-one, `P.O=2`, and zero-pole/dx4 target-section assumptions are
obsolete for this class.

## Completed exact milestones

`preflight_h92_q24_orbit42_component_valuation_qq.sage` binds the corrected
profile to the exact D12 equation and certifies the local I8* Weierstrass
orders `(2,3,14)`.  It does not construct the Riemann--Roch pencil.

`map_h92_q24_orbit42_i8star_physical_components_qq.sage` resolves the I8*
singularity over QQ, records eleven blow-up centres and twelve physical
components, and matches their graph and fibre multiplicities to the abstract
D12 root frame.  Exactly two spinor-arm orientations survive.  In the first,
the marked section meets `C11`; in the second it meets `C10`.  This is an exact
physical marking, not an A11 child equation.

The separate zero-pole audit reconstructs nine signed pairs, hence eighteen
explicit rational zero-pole sections on the exact D12 model.  Each section
satisfies the characteristic-zero Weierstrass identity and its pinned mod-53
regression.  These are the identity-class sections.  The selected orbit42
class is the missing nontrivial discriminant class, so this exact computation
closes the easy zero-pole shortcut without proving that no other construction
of the target section exists.

The proposed fast-q6 transport has also been closed exactly.  The named
unimodular equation-D13 to ambient H3-NS matrix passes its Gram and q24-fibre
round trips, but gives q6 degrees `435` for `O12` and `703` for `P42` (their
q8 degrees remain `30` and `48`).  They are high-degree q6 multisections, not
q6 rational points.  `run_h92_q24_orbit42_fast_parallel.py` is now a clean
rejection gate and never launches its obsolete transport/orientation jobs.

A second shortcut audit observes that twice the R3-zero orbit42 MW target is
an L1-six combination of the eighteen abstract identity-shell vectors.  At
the pinned prime `100003`, the exact equation points match that shell in eight
pointed-anchor-compatible ways.  Four rational degree-three relative halves
survive, but every chord branch polynomial is squarefree of degree `18`; none
has an A11 fibre.  This is a modular rejection of the rational-halving
shortcut, not a characteristic-zero non-existence theorem.

## Active gate

The next gate is the resolved-surface Riemann--Roch trivialization for
`D42=O+P+V`, carried out for both surviving physical orientations.  Success
requires a two-dimensional exact kernel, a compiled degree-two pencil, and an
exact child classification `A11/MW6`.  Only after that certificate exists does
the equation route advance to

```text
A11/MW6
 --q8 orbit922--> 2A5/MW7
 --q4 orbit472--> 3A3/MW8
 --q4 orbit323--> A3+2A2/MW10
 --q4 orbit207--> 5A1/MW12
 --q4 orbit52 --> 4A1/MW13
 --q4 orbit114--> 3A1/MW14
 --q4 orbit498--> 2A1/MW15
 --q4 orbit981--> A1/MW16
 --q6 orbit2247--> rootless/MW17
 -> pinned R17.
```

The operational stage ledger and version-locked launch commands live in
[`scripts/success-path/`](scripts/success-path/).  That ledger is a workflow
index; `MATH_STATUS.json` remains the sole authority for mathematical status.

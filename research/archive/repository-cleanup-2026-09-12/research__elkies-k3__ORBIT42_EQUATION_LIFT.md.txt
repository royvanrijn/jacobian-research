# H3 orbit42 equation lift

This note records the exact equation-level boundary for the selected

```text
D12/MW5 --q6 orbit42--> A11/MW6
```

continuation.  It is the canonical note for the orbit42 milestones below.
The complete lattice/chamber corridor to the pinned rootless/MW17 frame is
already exact. The orbit42 equation edge is exact. The later physical equation
route is now also complete through q8/orbit376 and q12/orbit5867 to the
source-identified, saturated R17 endpoint; see
[`PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md`](PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md).
The detailed milestones below retain the proof boundary as it stood while the
orbit42 edge itself was being closed.

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
class is in the nontrivial discriminant class, so these eighteen sections alone
do not generate it.

The two missing spinor-class zero-pole sections are now also exact.  On the
same R3-zero D12 model they have `deg(x)=3`, constant `y`, and are negatives
of one another.  The characteristic-zero solve has a unique cubic `x`, its
Weierstrass right-hand side is a rational constant square, and the two signs
reduce to the pinned residues

```text
x = 49445 u^3 + 5652 u^2 + 88005 u + 74701,
y = 45430 or 54573                         (mod 100003).
```

Thus the full twenty-point abstract zero-pole shell is represented exactly.
This supplies a construction aid and repairs the earlier incomplete spinor
enumeration; it does not replace the resolved line-bundle trivialization.
The reproducing command is

```bash
sage -python elkies-k3/scripts/lift_h92_q24_orbit42_spinor_zero_pole_sections_qq.sage
```

and its exact output is
`artifacts/local/elkies-k3/q24-orbit42-spinor-zero-pole-sections-qq.json`.

Combining this pair with the eighteen exact identity sections gives four
exact marking-compatible orbit42 candidates on the R3-zero model.  After the
degree-two zero change each has `P.O=3`, projective degrees `(9,9,3)`, and an
exact characteristic-zero Weierstrass identity.  Their reductions equal the
complete four-point mod-100003 census.  The coefficient height is substantial
(about 78,700 bits), so the equations live only in the generated artifact.
The command is

```bash
sage -python elkies-k3/scripts/construct_h92_q24_orbit42_exact_section_candidates_qq.sage
```

and the output is
`artifacts/local/elkies-k3/q24-orbit42-exact-section-candidates-qq.json`.
The points are exact; their four-way orbit42 marking is still pinned through
the modular shell isometry.  The resolved RR calculation must select the
physical orientation and provide the exact line-bundle kernel.

## Exact resolved-RR lift and A11 child

`lift_h92_q24_orbit42_resolved_rr_qq.sage` transports the selected exact
section to the minimal D12 equation and compiles the C10 physical orientation.
The smooth collision conditions give `9 -> 3`, and the weighted C01 valuation
with weights `(2,2,3)` gives `3 -> 2`.  Thus the exact resolved line-bundle
calculation has

```text
ambient = 9
smooth collision rank = 6
post-collision dimension = 3
resolved rank = 1
kernel = h0(D42) = 2.
```

The resulting binary quartic has degree four.  Its exact minimized Jacobian
has `deg(A,B,Delta)=(8,12,24)`, one `I12` fibre and twelve geometrically nodal
`I1` fibres.  Hence its root lattice is A11, its Euler number is 24, and its
MW rank is 6 under the repository's rank-19 marking.  The opposite C11
orientation has one additional independent C10-arm condition modulo the
pinned good prime and is rejected.

Reproduce the exact equation edge with

```bash
sage -python elkies-k3/scripts/lift_h92_q24_orbit42_resolved_rr_qq.sage
```

The output is
`artifacts/local/elkies-k3/q24-d12-to-a11-orbit42-resolved-rr-qq.json`, with
status `PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR`.

The equation-side child marking is certified separately by restricting the
exact pencil to all eighteen exact identity-shell curves at the pinned good
prime `100003`.  Their ordered new-fibre degrees are

```text
6,4,2,8,7,3,9,1,2,8,6,4,7,3,3,7,9,1.
```

Comparison with every A11 neighbor and all eight pointed shell isometries has
exactly two matches: orbit64/mapping7 in the selected C10 orientation and
orbit65/mapping6 in the spinor-conjugate C11 orientation.  Thus the concrete
child frame for the next edge is
`artifacts/local/elkies-k3/q24-downstream-lift/d12-c10a-zero-q6-frames/q6-o0064-r11-n132-d12-ad4a027cb197.txt`.
Reproduce this marking with

```bash
sage -python elkies-k3/scripts/certify_h92_q24_orbit42_a11_equation_marking.sage
```

Its output is
`artifacts/local/elkies-k3/q24-a11-equation-marking-orbit64-mod100003.json`.
The equation, quartic, Jacobian, and fibre classification are exact over QQ;
the shell-to-lattice orientation identification retains the stated pinned
good-reduction boundary.

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

The orbit42 gate is closed.  The next gate is the exact equation lift from the
selected orbit64 A11 frame through the q8 neighbor whose child is `2A5/MW7`.
The historical orbit number must not be copied into this frame.  The exact
construction transport

```bash
sage -python elkies-k3/scripts/certify_h92_q24_a11_q8_construction_fingerprint.sage
```

uses the stronger historical fingerprint

```text
new fibre = O + P - 2F,   vertical-root correction = 0
```

and splits the frame isometry into the `A11` root chain, the rank-six MW
quotient, and the integral glue condition.  Of the 32 root/MW combinations,
eight respect the glue and four land in the nef q8 search; the latter give
exactly equation-side orbits `12` and `2162`, both with `2A5/MW7` children.
The working target is equation orbit `12`, selected by the explicit minimum
MW-coordinate L1 convention (`2` rather than `3`).  This fixes a construction-
compatible divisor but does not yet construct its section function or pencil.

For the orbit-12 common local profile, the component-six section builder now
supports recursive infinity elimination.  At the preferred discovery prime
`p=43`, eliminating eight of the 18 `R` coefficients gives a 26-variable,
28-equation, 834 KB sparse system with 16 infinity branches.  A bounded msolve
benchmark completed its degree-10 matrix and reached degree 11 before being
stopped; this is a search-performance observation, not a section or
non-existence result.

<!-- status-consumer: EC-K3-H3-A11-Q8-CONSTRUCTION-TARGET -->

The construction audit one level deeper changes the preferred section target.
Under the pinned equation marking, the eighteen exact identity-shell points
span rank five but index five in the saturated hyperplane with sixth MW
coordinate zero.  The earlier pole-order-four vector is in the wrong coset.
An exhaustive exact height/correction enumeration proves that the selected
orbit-12 coset has no section with `P.O <= 4` and exactly two with `P.O=5`.
The preferred one is

```text
M = (1,0,0,0,0,1),  height = 47/4,  correction = 9/4,
P12 = M + S6 - 2*S2 - 2*S8.
```

Here `S_i` denotes the exact identity-shell point in zero-based equation
order.  The equality is literal in the pinned A11 MW marking.  A new
free-infinity chart covers every leading smooth-fibre value at once; over
`GF(100003)` it has 36 variables, 37 equations and maximum degree six.  A
bounded four-thread msolve mode-42 run completed two degree-eight reductions,
entered a third `257692 x 2857438` matrix and stopped at 600 seconds (about
13 GB resident).  It found no section, and this bounded run is not a
non-existence result.  The exact next gate is characteristic-zero coordinates
for `M`, followed by the displayed group-law reconstruction of `P12`.

The degree-one alternative has now been exhausted exactly, without a section
ansatz.  Möbius inversion of the resolved pencil and binary-quartic
covariants transport identity-shell indices 7 and 17 and spinor index 0 to
three exact characteristic-zero A11 points.  The pointed-quartic replay gives
the opposite point as the negative covariant point in all three cases.  The
selected unimodular D12-to-A11 transition then proves that the sixth child MW
coordinate is exactly the fifth parent MW coordinate.  All eighteen identity
vectors and both spinor vectors have that parent coordinate zero, so none of
these exact degree-one transports supplies the missing direction.  Switching
to the other construction-compatible q8 orbit 2162 only changes the required
sixth coordinate from `+1` to `-1`.

The smallest primitive parent carrier is the D12 vector
`(0,0,0,0,1)`, with height 12, correction 0 and `P.O=4`; it is not in the
exact zero-pole shell.  A proposed quintic shortcut was rejected by exact
marking replay: it applied the orbit64 transition for the selected `R3`-zero
D12 frame to coordinates taken in the distinct `A0`-zero frame.  In the
compatible selected marking, `close_P24` has A11 degree 46 and MW vector
`(33,-77,31,-38,7,1)`, while `oldI9_A0` has degree 4 and vector
`(2,-6,3,-3,1,0)`; the claimed group word does not equal `M` and the
index-five shell is not saturated.

An exact audit of all forty stored explicit (-2)-classes finds no positive
single-carrier replacement.  The first positive subset is
`AJ(H3_simple_2)-2*AJ(H3_simple_8)` plus shell points, requiring degree-40
and degree-44 traces, so it is a correctness fallback rather than the active
construction.  An exhaustive equation-cost scan also tested lateral A11
neighbours: the apparent best candidates 849 and 591 fail the full nefness
gate, while the first passing candidates lead to different ADE types or lack
a certified continuation to pinned R17.

The active no-large-elimination direction therefore keeps equation orbit12
and exploits its split `I12` fibre directly.  With the formal nodal centre
`c(s)`, the identity

```text
y^2 = (x-c)^2*(x+2*c) + g(s),   ord_s(g)=12
```

and the component-3 substitution `x-c=s^3*Q/Z^2`, `y=s^3*R/Z^3`
turn the local section condition into a norm/Pell-style coefficient
recurrence.  Fixed-infinity generation reduces the benchmark to 16 variables
and 18 equations, but naive full substitution raises total degree to 36;
that is evidence to use bidirectional coefficient recurrences or resolved
linear Riemann--Roch, not a large Gröbner calculation.  Modular systems remain
bounded discovery aids only.

<!-- status-consumer: EC-K3-H3-A11-Q8-TARGET-COSET-BRIDGE 8d17ab150a7e3567 -->

That resolved route has now closed the q8 edge exactly.  The component-3
chart first reconstructs the small residual `R=P12-M` over `QQ`; exact
fraction-free group law then gives `P12` and `D=P12-O_pinned` without lifting
the former 36-variable difference chart.  The final horizontal has projective
degrees `(16,24,6)`, `P.O=6`, and central `I12` depth 6.

For `O+D-2F`, use

```text
a=AA/Z^2, deg(AA)<=10;  b=BB/Z, deg(BB)<=2.
```

This is a complete 14-dimensional ambient.  The congruence
`AA*X=BB*Y mod Z^2` has rank 12.  Inverting `X mod Z^2` reduces its exact
solution to `AA=BB*Y/X mod Z^2` and the single condition that the degree-11
coefficient vanish, leaving `h0=2`.  No Groebner basis is used.
Fraction-free chord elimination gives `Z^6` times a quartic.  Its globally
minimal Jacobian has degrees `(8,12,24)`, fibres `2I6+12I1`, root lattice
`2A5`, Euler number 24, and MW rank 7 in the fixed rank-19 marking.

The equation marking is exact as well.  At the old `I12` value the horizontal
specializes to the node.  If `N=AA1-T*AA0` and `Db=T*BB0-BB1`, the affine
component has normalized quartic ordinate

```text
w_affine=(N^2-3*X*Db^2)/Z^3.
```

The opposite sign is the only other degree-one old-fibre curve,
`old_A11_component_9`, and the pointed-quartic map sends it to infinity as
the selected child zero.  The ten degree-zero curves form the two physical
chains `(0,3,4,5,10)` and `(1,2,6,7,8)`.  The full equation-A11 to
component-9-zero child transport and its inverse both have determinant `-1`
and preserve the NS Gram form exactly.

The primary outputs are
`artifacts/local/elkies-k3/q24-a11-to-2a5-q8-resolved-rr-qq.json` and
`artifacts/local/elkies-k3/q24-a11-to-2a5-q8-equation-marking-qq.json`, with
statuses `PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR` and
`PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING`.

<!-- status-consumer: EC-K3-H3-A11-Q8-QQ-2A5 b7aaf4bf483eac68 -->

The equation route was subsequently re-audited in the physical component
chamber. The q6/orbit1307 and changed-zero alternatives below retain exact
lattice information but are no longer the active equation continuation. The
current physical route is

```text
2A5/MW7
 --q4 orbit208--> 3A3/MW8
 --q4 orbit1584--> D4+A3+3A1/MW7
 --q4 orbit164--> 2A3+2A1/MW9
 --q8 orbit376--> 4A1/MW13
 --q12 orbit4484--> rootless/MW17
 -> pinned R17.
```

The first three arrows have exact characteristic-zero equations and effective
zeros. The q8/orbit376 and q12/orbit4484 arrows are exact marked-lattice
edges but remain equation-open. See
[`PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md`](PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md).

<!-- status-consumer: EC-K3-H3-Q24-O42-QQ-A11 ffa4308117c55056 -->
<!-- status-consumer: EC-K3-H3-A11-Q8-CONSTRUCTION-TARGET c892eec88af45f08 -->
<!-- status-consumer: EC-K3-H3-A11-R17-Q6O1307-PROMOTED-LATTICE-ROUTE 9ee5630063324558 -->
<!-- status-consumer: EC-K3-H3-A11-R17-Q4O230-Q6O1315-PROMOTED-LATTICE-ROUTE a8889fef54ee3b47 -->
<!-- status-consumer: EC-K3-H3-FIRST-Q8-Q4O11-PROMOTED-LATTICE-ROUTE ba6ee6488fce2411 -->

The operational stage ledger and version-locked launch commands live in
[`scripts/success-path/`](scripts/success-path/).  That ledger is a workflow
index; `MATH_STATUS.json` remains the sole authority for mathematical status.

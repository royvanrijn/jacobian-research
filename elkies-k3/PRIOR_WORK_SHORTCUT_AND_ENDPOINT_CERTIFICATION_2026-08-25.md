# Prior-work shortcut and endpoint-certification directive

Date: 2026-08-25.

## Outcome

The historical construction already establishes that the determinant-948
rank-19 K3 admits a rootless elliptic fibration of Mordell--Weil rank 17 over
`QQ`.  Public sources also describe the construction method: start from an
easy elliptic model, pass through a chain of 2-neighbours (and occasionally a
3-neighbour), discover coefficients modulo a small prime, lift p-adically,
recognize rational coefficients by lattice reduction, and verify the resulting
identities exactly.

This is substantially the method being reconstructed here.  The missing
historical deliverables are concrete rather than conceptual: no publicly
readable primary source inspected by the repository supplies the final
Weierstrass coefficients, a rank-17 section basis, the actual neighbour chain,
or its coordinate transformations.  See
[`../elliptic-curves/ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md`](../elliptic-curves/ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md).

The consequence for current work is a change of emphasis.  The remaining
neighbour chain is a **construction scaffold**, not a requirement to build a
maximal theorem package at every intermediate fibration.  Exact work that is
needed to reach and certify the final equation remains mandatory.  Exhaustive
intermediate orbit classification, alternate-route optimization, and complete
MW/physical-component recovery should be deferred unless they directly unlock
the next equation or are needed by the final proof.

This directive does not weaken `MATH_STATUS.json`: no endpoint claim is
promoted until the final certificate and its replay exist.

## What prior work contributes

Elkies's 2007 lecture report states all of the following:

- the non-CM rational point on the level `6*79` Shimura curve produces a
  rank-17 elliptic K3 over `QQ(t)`;
- the desired equation is reached from an easier elliptic model through
  2-neighbours and occasional 3-neighbours;
- modular search, p-adic Newton lifting, lattice reduction, and exact
  substitution are the intended computational pipeline;
- the rootless fibration has 1,311 pairs of polynomial sections of height four.

Primary source:
[N. Elkies, *Three lectures on elliptic surfaces and curves of high rank*](https://arxiv.org/html/0709.2908).
The Shimura/Kumar source construction is described further in
[N. Elkies, *Shimura curve computations via K3 surfaces of Neron--Severi rank at least 19*](https://arxiv.org/html/0802.1301)
and the explicit genus-two-to-`E8+E7` model comes from
[A. Kumar, *K3 surfaces associated with curves of genus two*](https://arxiv.org/abs/math/0701669).

These sources validate the overall route architecture but do not identify the
historical intermediate markings.  They therefore cannot replace the two
remaining equation lifts unless the missing original data are recovered.

## Minimum construction record for an intermediate edge

For an intermediate edge used only to continue the construction, retain the
smallest exact package that makes the next edge reproducible:

1. the exact parent equation and selected primitive pencil;
2. the exact two-dimensional Riemann--Roch space, or an equivalent exact
   function-field construction of the pencil;
3. the binary quartic/genus-one model, its Jacobian, and exact identities;
4. a rational origin and enough component/curve pointing to define the next
   selected divisor unambiguously;
5. the exact parent-to-child rational maps or the minimum lossless class
   transport required for the continuation;
6. the reproducing command, software assumptions, and output hash.

Do **not** block construction merely to complete any of the following when it
is not used by the next edge or final endpoint proof:

- a complete MW basis at the intermediate fibration;
- exhaustive automorphism or neighbour-orbit classification;
- proof that the route or compiler score is globally optimal;
- alternate suffixes and changed-zero loops;
- full pointing of every old curve and every physical component;
- a general theorem extracted from a single successful edge.

If an omitted marking later becomes necessary, recover that specific marking
then.  Never replace an exact requirement by an ADE/MW label alone.

## Direct R17 endpoint gate

Once a rootless equation is obtained, certify the endpoint directly rather
than replaying every optional intermediate invariant.  The final package must
contain:

1. an exact elliptic K3 Weierstrass equation over `QQ` and an exact minimal
   fibre/Euler audit;
2. proof that every singular fibre is irreducible, so the fibre-root lattice is
   zero;
3. seventeen rational sections, each checked by literal substitution;
4. their exact canonical-height Gram matrix and an integral isometry to
   [`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt), including
   determinant `948`, torsion, and saturation checks;
5. proof that the geometric Picard number is `19`, either by exact
   identification with the already-certified non-CM H3 surface or by an
   independent good-reduction upper bound;
6. an exact source-identity certificate: composed birational maps, a lossless
   NS isometry, or another exact invariant package that proves this is the
   intended level-474 surface rather than merely a K3 with similar numerical
   data;
7. pinned replay commands, software versions, and hashes, followed only then
   by the `MATH_STATUS.json` update and regenerated `STATUS.md`.

Items 2--5 and Shioda--Tate then give geometric and arithmetic Mordell--Weil
rank exactly `17`.  The general justification is recorded in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md),
especially Theorems A and F.

## Section recovery at the endpoint

Do not transport a full MW basis through every remaining quartic/Jacobian
conversion unless that is cheaper than direct recovery.  On a rootless
elliptic K3, height-four sections are polynomial sections.  The known count of
1,311 pairs gives a strong regression fingerprint and makes the historical
finite-field/Hensel strategy natural:

1. solve the bounded polynomial-section system over one or more suitable good
   primes;
2. organize solutions using the pinned R17 short-vector graph and select a
   rank-17 Gram basis;
3. Hensel/Newton lift only that basis or a small generating shell;
4. recognize rational coefficients and certify over `QQ` by substitution and
   exact height pairing;
5. verify saturation and, as a regression rather than a prerequisite, recover
   the full 1,311-pair shell if affordable.

This is preferable to a fresh unconstrained 17-section nonlinear solve and may
be preferable to transporting seventeen sections through all remaining edges.

## Current routing instruction

The q4/orbit1584 and q4/orbit164 equation edges are now certified exactly as
described below.  Preserve those results narrowly.  The intended construction
continuation is:

```text
... --q4/orbit1584--> D4+A3+3A1/MW7
    --q4/orbit164--> 2A3+2A1/MW9
    --q8/orbit376--> 4A1/MW13
    --q12/orbit5867--> rootless R17/MW17.
```

Until this route fails an exact construction gate:

- prioritize q8/orbit376, then the preferred optional q12/orbit5867 equation
  lift; retain q12/orbit4484 as the certified fallback;
- carry only the markings needed for the next divisor and final source
  identity;
- park q323 work, alternate suffix searches, changed-zero optimization, and
  global score improvement unless they directly unblock one of those lifts;
- after the rootless equation appears, switch immediately to the direct R17
  endpoint gate above.

The separate fixed-ADE audit corrected a false positive in the parked
q4/orbit323 work. Branch 33 inverse matched the target profile and three
inherited intersections only in isolation; it occurs in neither complete
signed eleven-section graph solution. A resolved summand replay reconstructs
the same excluded branch exactly. Exhaustive q4 enumeration then supplies an
exact QQ replacement: lifted branch 16 with the q4/orbit1599 trace has a
`5 -> 3 -> 2` RR plane and `I4+2I3+14I1`, hence `A3+2A2/MW10`. This restores
the fixed-ADE equation step without a Groebner basis. The two complete signed
graph solutions agree on branch 16 itself: sign `+1` and the same raw NS
class. In physical coordinates the certified lattice fibre satisfies
`D=O+P+V`, with `V` supported in the old root lattice and zero MW tail.
Therefore this edge has no pinned-suffix marking ambiguity; the remaining
twofold choice affects only an unused lifted branch.  This statement concerns
the q4/orbit1599 edge itself: its resulting marked `A3+2A2` frame is not the
stored canonical frame that feeds q4/orbit207, so it remains an optional
fixed-ADE route rather than a replacement for either the fixed corridor or the
lower-cost q4/orbit1584 route.

## Corrected fixed-corridor q4/orbit323 equation on 2026-08-26

The complete eleven-section graph excluded the target from that particular
rank-seven lifted subgroup; it did not prove that the marked target section was
absent.  The missing eighth Mordell--Weil direction is supplied by the already
exact inherited `C7` section.  In the correctly oriented marked quotient,

```text
2*T = P8 + 2*P18 + P33 - 2*C7  modulo the trivial lattice.
```

The opposite global orientation has `P.O=3`; the displayed orientation has the
certified target invariant `P.O=1`.  Its duplication quartic factors as one
linear and one cubic factor.  The unique marked rational half has compact
degrees `x=(6,2)`, `y=(8,3)` and maximum coefficient sizes 94 and 135 bits.
The doubled NS relation is exact, and the trivial lattice has Smith invariants
all equal to one, ruling out a torsion ambiguity.  This construction takes a
fraction of a second and uses neither a new modular shell nor a Groebner basis.

For the corrected divisor, `D=O+T+V` with

```text
V = second_old_I6_I4_component_2 + second_old_I6_I4_component_3.
```

Thus the relevant connected resolved trace is at compact `t=0`.  The exact
Riemann--Roch calculation is `5 -> 3 -> 2`; the chord radicand is two squared
linear factors times a quartic.  Binary-quartic invariants give a minimal
Jacobian with `I4+2I3+14I1`, Euler number 24, root data `(7,24,36)`, and
`A3+2A2/MW10`.  This closes the equation/RR/Jacobian portion of the prescribed
q4/orbit323 edge.

The child-zero gate is also now exact.  A split-`I4` toric arc selects
`old_A11_component_2` as `W=+L0(u)` in the stored positive-square-root
normalization; direct pointing gives the stored child invariants with exact
`81/729` scaling and constructs the opposite branch over `QQ(u)`.  Reflecting
the entire remaining fixed suffix in the old-zero wall, then replaying nine
vertical component reflections, preserves every adjacent degree-two pairing
and converts the old q4/orbit207 class into a physical q12 edge to `5A1/MW12`.
Its marked horizontal has `P.O=10`, height `65/3`, and
`D=O+P-4F`, with no vertical residual.  The edge and its finite horizontal-wall
gate are lattice-certified; its horizontal equation is the remaining gate.

For that equation gate, the q12/orbit5867 polynomial-shell method has been
reused instead of starting a degree-16 Abel inversion.  The fibre-safe prime
61 gives a complete signed polynomial `P.O=0` shell of size 602, with 120
ordinary rank-12 Hensel candidates.  Prime 31 is explicitly rejected because
the two marked `I3` valuations become `(4,4)`.  Resolved component naming,
target-coset selection, and characteristic-zero lifting remain open; this
bounded modular shell is not yet the q12 equation.

## Exact physical-suffix closeout on 2026-08-25

The first two edges of the route above have passed their characteristic-zero
construction gates.

For q4/orbit1584, the exact physical marked-edge certificate is combined with
an exact `QQ` Riemann--Roch calculation.  The unique resolved double-branch
jet gives a two-dimensional section space and a degree-four binary quartic.
Its minimal Jacobian has finite fibres `I4+3I2+8I1`, an `I0*` fibre at
infinity, Euler number 24, and hence root lattice `D4+A3+3A1` and MW rank 7
under the already-certified Picard-rank-19 source identity.  The
second-I6-affine component is then attached as an exact rational zero.  Exact
specialization at the three `I2` fibres distinguishes the old `C0` and `C4`
branches and supplies the `C0` marking required by q4/orbit164.

For q4/orbit164, the selected horizontal meets two identity components of
finite `I2` fibres.  Their two exact rational double-branch values determine
the linear interpolation used in

```text
u = (m-l(T))/((T-r0)(T-r1)).
```

Substitution removes the denominator square and realizes the complete
`4 -> 2 -> 2` Riemann--Roch calculation without a Groebner basis.  The
resulting degree-four quartic has a minimal Jacobian with finite fibres
`I4+2I2+12I1`, an `I4` fibre at infinity, Euler number 24, root lattice
`2A3+2A1`, and MW rank 9.  The unique finite-`I4` identity branch is exactly
old `C8`; the opposite quartic sign reconstructs the contracted node.
Pointing at `C8` reproduces the stored short Jacobian and supplies the source
zero for q8/orbit376.

The inherited-`P1` no-Groebner construction has also passed its first complete
modular gate.  Certified pointed maps through q4/orbits 208, 1584, and 164
give base degrees `3 -> 6 -> 7`; fibrewise Abel reduction uses a `7 x 8`
kernel in `L(8O)`.  Over `GF(131)` the resulting 122 good fibres interpolate
an exact q4/orbit164 section with coordinate degrees `x=(32,28)` and
`y=(47,42)` (inside the predicted numerator bound 48), verified on the
Weierstrass equation with independent holdouts.
The marked-tail calculation gives

```text
H = T - 8*C1 + 5*C2 - 2*C3 - 7*C7 + R,
R = (0,0,0,0,-1,-1,2,0,0),
```

where `T` is the Abel trace and `H` is the q8/orbit376 horizontal.  The named
degree-one sections span only the first four MW coordinates.  An exact
fourfold pole-growth audit now gives the equation-side basis `B0,...,B7` the
correct height Gram, of determinant `459/8`.  Of its 16 integral embeddings
in the C8-pointed marked MW9 lattice, eight survive the valid first-seven
component profiles.  All eight contain `R`, with different basis words, but
the complete modular Abel trace selects exactly one of them by the certified
q8 pole fingerprint.  The resulting identity over `GF(131)(t)` is

```text
H = T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7,
deg(x_H)=(12,8), deg(y_H)=(18,12), P.O(H)=4.
```

This also corrects an unresolved-node error: `B7` has an odd infinity-`I4`
label, not the former coarse label 2, and the exact section
`N=2*B0+B5+B7` has height `13/4`, not `3`.  A single local resolved `I4`
chart is not needed for the modular selection: pole growth already chooses
the marked embedding.

The horizontal itself has now been reconstructed over `QQ(t)` without first
lifting the larger Abel trace.  Twenty-two exact good-prime horizontals give
a 566-bit CRT modulus.  Independent coefficientwise rational reconstruction
still leaves 12 unresolved coefficients because monic normalization hides a
large common scale.  A simultaneous projective LLL reconstruction of the
22- and 32-entry coordinate vectors instead selects its first candidate pair,
with primitive-vector maxima 363 bits for `x` and 526 bits for `y`.  Literal
substitution proves

```text
y_H^2 = x_H^3 + A*x_H + B,
deg(x_H)=(12,8), deg(y_H)=(18,12), P.O(H)=4
```

over `QQ(t)`, and coefficient reduction reproduces every one of the 22 input
prime artifacts.  Doubling twice gives fourfold `x` pole degree 172 and hence
canonical height `(4+172)/16=11`, exactly matching the marked q8 class.  The
calculation uses only fibrewise `7 x 8` kernels,
univariate interpolation, affine elliptic group law, CRT, and small LLL; no
Groebner basis is used.

The q8 Riemann--Roch compilation is also exact. Writing
the recovered horizontal as `x_H=X/Z^2`, `y_H=Y/Z^3` with
`deg(X,Y,Z)=(12,18,4)`, the saturated chord ambient is

```text
a=AA/Z^2, deg(AA)<=8;  b=BB/Z, deg(BB)<=2.
```

It has dimension 12.  The collision congruence
`AA*X=BB*Y mod Z^2` has rank 8 and leaves an exact four-plane.  Inverting `X`
modulo `Z^2` constructs the plane directly; an independent full `8 x 12`
coefficient matrix and all 22 CRT primes verify it. The vertical divisor is
the union of the two complete nonidentity `A3` chains in the old `I4` fibres.
Split-toric resolution produces one quotient row from each chain; the rows
are independent, so the resolved plane has `h0=2`. Exact division by `Z^4`
leaves a quartic whose minimal Jacobian has `4I2+16I1`, Euler number 24, and
root data `(4,8,16)`, hence `4A1/MW13`. The selected marked embedding fixes
the finite `T=0` I4 chain as old components `6`, missing, `5`; the exact `B0`
tangent selects old component 6 as the initial pointed origin. The
lattice-preferred `P1229` is instead the nonidentity component of the finite
I2 at `T=25281/168246841`. A rational arc on its exceptional conic selects
quartic sign `-1`; direct pointing verifies `81*A_pointed=A_child` and
`729*B_pointed=B_child`. This closes the q8 equation-to-NS zero marking
without a group-law translation and without a Groebner basis.

The exact marked lattice continuation is

```text
3A3/MW8
  --q4/orbit1584--> D4+A3+3A1/MW7
  --q4/orbit164--> 2A3+2A1/MW9
  --q8/orbit376--> 4A1/MW13
  --q12/orbit5867--> rootless/MW17.
```

The q12/orbit5867 edge is the preferred optional final compiler target after
q8/orbit376. It is fully marked and pinned to the same R17 lattice as the
q12/orbit4484 fallback. Its nominal optimized four-`P.O=0` lattice word has
q4/orbit164 parent degrees `(3,2,1,2)` and parent `a-b` values `(2,2,1,1)`,
lowering the planning totals from `(10,8)` for orbit4484 to `(8,6)`. The
canonical comparison is
[`../artifacts/generated-results/elkies-k3-h3-q4o1584-route-optimization-handoff.json`](../artifacts/generated-results/elkies-k3-h3-q4o1584-route-optimization-handoff.json).
Its durable identity gate is the certificate label `q12o5867_after_q8o376`
and exact fibre fingerprint
`d676cab5918a08add2f743081cf76932d02adc9e18ab9c808dc05760fe0157dd`.
The index 30357 in the expanded bounded frontier is sample-local and does not
replace the promoted orbit5867 label.

The nominal word is not equation-effective: complete polynomial `P.O=0`
shells at primes 83, 89, 137, and 151 contain no Q1 candidate. This is a
four-prime negative experiment, not a characteristic-zero nonexistence proof.
The actual exact compiler word is
`499+500+69+511-489+933-913`, corresponding to new classes
`499+500+69+511` and equation-coordinate correction cosets
`-489+933-913`. Its parent degrees are `(4,2,1,5)` and parent `a-b` values are
`(4,2,1,4)`. Ordinary regular Hensel lifts and exact group law construct the
four required sections over QQ without elimination.

The resulting horizontal has x numerator/denominator degrees `(24,20)`, y
degrees `(35,30)`, `P.O=10`, q8-I2 profile `(1,1,1,1)`, height 22, and the
exact marked q12 tail. For `D=O+P-4F`, the smooth chord ambient has dimension
22 and the exact congruence has rank 20, hence `h0=2`, with no vertical rows.
After removing the squared degree-10 factor, the binary quartic gives a
minimal Jacobian of degrees `(8,12,24)`. Its discriminant is squarefree of
degree 24 and infinity is smooth, so the geometric fibres are `24I1`, Euler
number 24, and root rank zero. No Groebner basis or surface elimination is
used.

The direct endpoint basis is also exact. At the smallest good prime 131, an
exhaustive generated-C polynomial shell has 2,622 signed sections and every
13-by-12 coefficient Jacobian has rank 12. Exact modular intersections select
seventeen sections with height determinant 948 and determinant-minus-one
transport to pinned R17. All seventeen lift independently over QQ at 4,096
p-adic digits; their maximum coefficient sizes lie between 7,895 and 7,933
bits. Literal QQ substitution and exact finite/infinity intersections recover
the same determinant-948 Gram. Torsion is trivial because rootlessness gives
height `4+2(P.O)>0` for every nonzero section. Hence MW rank at least 17 is
unconditional.

The remaining endpoint gates admit a short exact closeout. Write the stored
q12 binary quartic as `w^2=q_u(v)`, where `v` is the q8 base and `u` is the
new base. The old finite `I2` support `v=0` gives
`q_u(0)=W_0(u)^2` with `W_0` an exact quadratic in `QQ[u]`. Thus
`(v,w)=(0,W_0(u))` is a `QQ(u)`-point on the genus-one pencil. Applying the
standard nonbranch pointed-quartic construction gives generalized
Weierstrass coefficients whose short invariants satisfy

```text
81*A_pointed = A_endpoint,
729*B_pointed = B_endpoint.
```

Equivalently, `x_endpoint=9*x_pointed` and
`y_endpoint=27*y_pointed`. The terminal Jacobian is therefore the same
elliptic K3 as the q12 pencil, not a nontrivial torsor. Composing this exact
pointing with the already-pointed forward chain proves identity with the
certified H3 source.

For the Picard upper bound, exact point counts at two good primes are:

| `p` | `#X(F_p)` | `#X(F_{p^2})` | residual trace `s1` | pair trace `r` | reduction NS discriminant square class |
|---:|---:|---:|---:|---:|---:|
| 131 | 19,308 | 294,853,764 | -343 | -212 | `-23,700 ~ -948` |
| 137 | 21,036 | 352,653,204 | -337 | -200 | `-35,076 = -948*37` |

The nineteen explicit characteristic-zero divisor classes contribute `19p`
and `19p^2` to the first two `H^2` traces. The residual cubic has one
eigenvalue `-p` at both primes; the other two have trace `r`. Since `r/p` is
not in `{-2,-1,0,1,2}`, that pair is not `p` times roots of unity. The Tate
theorem for K3 surfaces in odd characteristic therefore gives geometric
Picard rank 20 at both reductions. Over `F_{p^2}`, Artin--Tate gives the
displayed discriminant square class `r^2-4p^2`. Their ratio is `37/25`, not a
rational square. If the characteristic-zero Picard rank were 20, its
Neron--Severi discriminant would specialize to the same square class at both
rank-20 reductions. Hence `rho<=19`; the nineteen explicit classes give
`rho=19`.

Shioda--Tate now gives exact MW rank 17. To prove saturation, let `L` be the
displayed determinant-948 section lattice. If `L` had proper index `n` in the
full integral MW lattice, then `n^2` would divide 948, so `n=2`. The mod-2
radical of the pinned Gram is one-dimensional, represented by

```text
(1,0,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0)/2.
```

Its norm is 73. But on a rootless elliptic K3 the full MW height lattice is
even integral: there are no fibre corrections and every self-height is
`4+2(P.O)`. The unique possible index-two enlargement is odd, so it cannot
occur. The seventeen sections therefore form the full saturated geometric MW
lattice, integrally isometric to pinned R17, with exact rank 17, trivial
torsion, and determinant 948.

Every edge in this displayed route is primitive, nef, of old-fibre degree
two, and has mutually inverse unimodular Neron--Severi transports.  The final
rootless positive frame is integrally isometric to
[`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt).  This is an
exact marked-lattice and endpoint-isometry statement. Together with the
reconstruction above it supplies the exact P1229-pointed q8 horizontal,
resolved pencil, quartic, `4A1/MW13` Jacobian, q12 horizontal, and rootless
q12/orbit5867 equation and seventeen exact endpoint sections with a pinned
determinant-948 height Gram. The exact quartic pointing, two-prime
van-Luijk/Artin--Tate calculation, and even-overlattice audit above close the
source-identity, Picard-rank, and saturation gates.

The exact equation replays are:

```bash
sage -python elkies-k3/scripts/certify_h92_q4o208_physical_q4o1584_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_equation_marking_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_physical_q4o164_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o164_c8_equation_marking_qq.sage
sage -python elkies-k3/scripts/reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage
sage -python elkies-k3/scripts/compile_h92_q4o164_q8o376_smooth_rr_qq.sage
sage -python elkies-k3/scripts/construct_h92_q12o5867_target_horizontal_qq.sage
sage -python elkies-k3/scripts/compile_h92_q12o5867_smooth_rr_qq.sage
sage -python elkies-k3/scripts/construct_h92_q12o5867_rootless_p0_shell_mod131.sage
python3 elkies-k3/scripts/select_h92_q12o5867_rootless_mod131_basis.py
sage -python elkies-k3/scripts/lift_h92_q12o5867_rootless_selected_basis_qq.sage
sage -python elkies-k3/scripts/certify_h92_q12o5867_rootless_height_basis_qq.sage
sage -python elkies-k3/scripts/certify_h92_q12o5867_endpoint_qq.sage
```

Their terminal statuses are respectively
`PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN`,
`PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING`,
`PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN`, and
`PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING`.

This closeout completes the q12/orbit5867 equation, proves its exact identity
with the H3 source, establishes geometric Picard rank 19, and identifies the
full saturated geometric Mordell--Weil lattice with pinned R17. The
q12/orbit4484 lattice edge remains a certified but unnecessary fallback.

<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-R17-BASIS a2097150acf00645 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->

## Highest-value external shortcut

The single highest-value external action is recovery of Elkies's original
data.  A focused request should ask for any surviving copy of:

- the rank-17 Weierstrass coefficients;
- the 17 section coordinates or an MW basis;
- the 2/3-neighbour sequence and coordinate transformations;
- the quadratic base-change map and eighteenth section;
- the specialization parameters used for the rank-28 or rank-29 searches;
- Magma, PARI/GP, Sage, or handwritten computation files from the 2006--2007
  construction.

No external message should be sent without user authorization.  If any such
data are recovered, import them with provenance and hashes and first test them
against determinant `948`, the 1,311 height-four pairs, and the level-474
source marking before changing the construction route.

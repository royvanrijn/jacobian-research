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
Groebner basis is used.  The remaining q8 gate is the resolved
Riemann--Roch compilation from this exact horizontal to the `4A1/MW13`
child equation.

The exact marked lattice continuation is

```text
3A3/MW8
  --q4/orbit1584--> D4+A3+3A1/MW7
  --q4/orbit164--> 2A3+2A1/MW9
  --q8/orbit376--> 4A1/MW13
  --q12/orbit5867--> rootless/MW17.
```

The q12/orbit5867 edge is the preferred optional final compiler target after
q8/orbit376.  It is fully marked and pinned to the same R17 lattice as the
q12/orbit4484 fallback.  Its optimized four-`P.O=0` word has q4/orbit164
parent degrees `(3,2,1,2)` and parent `a-b` values `(2,2,1,1)`, lowering the
corresponding totals from `(10,8)` for orbit4484 to `(8,6)`.  This is an exact
lattice/Mordell--Weil word and a compiler-cost improvement, not yet a
characteristic-zero equation certificate.  The canonical comparison is
[`../artifacts/generated-results/elkies-k3-h3-q4o1584-route-optimization-handoff.json`](../artifacts/generated-results/elkies-k3-h3-q4o1584-route-optimization-handoff.json).

Every edge in this displayed route is primitive, nef, of old-fibre degree
two, and has mutually inverse unimodular Neron--Severi transports.  The final
rootless positive frame is integrally isometric to
[`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt).  This is an
exact marked-lattice and endpoint-isometry statement.  Together with the
reconstruction above it now supplies the exact q8 horizontal, but it does
**not** yet supply the q8/orbit376 child equation, a q12/orbit5867 (or
fallback q12/orbit4484) equation, seventeen endpoint sections, or the direct
endpoint certificate.

The exact equation replays are:

```bash
sage -python elkies-k3/scripts/certify_h92_q4o208_physical_q4o1584_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_equation_marking_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o1584_physical_q4o164_rr_qq.sage
sage -python elkies-k3/scripts/certify_h92_q4o164_c8_equation_marking_qq.sage
sage -python elkies-k3/scripts/reconstruct_h92_q4o164_q8_horizontal_crt_qq.sage
```

Their terminal statuses are respectively
`PASS_EXACT_QQ_PHYSICAL_Q4O1584_D4_A3_3A1_RR_AND_JACOBIAN`,
`PASS_EXACT_QQ_Q4O1584_SECOND_AFFINE_POINTING_AND_C0_MARKING`,
`PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN`, and
`PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING`.

This closeout strengthens the directive rather than changing it: the active
equation gate is now the q8/orbit376 resolved RR compiler, followed preferably
by q12/orbit5867 and the direct R17 endpoint package; q12/orbit4484 remains the
certified fallback.

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

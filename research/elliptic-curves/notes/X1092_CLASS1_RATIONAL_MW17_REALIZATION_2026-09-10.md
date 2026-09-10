# X1092 class 1: an exact rational MW17 realization

**Current operation:** [ordinary prospective search is running](CLASS1_PROSPECTIVE_ORDINARY_SEARCH_2026-09-10.md)
on the byte-frozen parent. The strict-class prerequisite in the original
realization handoff below has been removed for ordinary search. The separate
arithmetic construction intake remains UNKNOWN; other parent classes are paused.

**Completed:** a primitive rational marked nef `U`, an exact rational elliptic
equation and birational maps to the recovered source K3, and seventeen rational
generic sections with a saturated height Gram of determinant **1092**. The
generic rank is **exactly 17**. The new frame is census class **1**, integrally
nonisometric to curve302's class **6**.

**Still UNKNOWN:** construction of a new arithmetic strict class modulo the
full inherited image on this marked parent. Four bounded generic carrier
constructions passed; none is asserted to be an arithmetic strict class.
The strict preflight fails closed. **No parameter panel was commissioned or
run.** There is no new specialized rank or conductor claim.

## Exact transport and the prescribed witness

The starting certificate is the
[complete census](DET1092_ROOTLESS_J2_CENSUS_2026-09-10.md), with the retained
shared-core witness: known embedding **25**, new embedding **33**, sixth index
**155**, common primitive rank-16 core of determinant **4100**. Its exact
common-core bases, integral inclusions and Gram are retained and checked in the
[marking packet](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_marking_v1.json).

A low-degree transport of the requested *frame type* was found among the
already retained generic degree-two quotient representatives: mask **109158**.
In the rational source marking `NS=U + (-G302)`, put

```
w = (-1,1,0,-1,-1,-1,0,1,1,0,0,0,0,0,-1,1,1)
w G302 w^t = 12
D = (3,2,w) = O_old + P_w - F_old
O = (-1,1,0,...,0) = O_old.
```

For the row matrix `T=[D; D+O; W]`, where `W` is the saturated orthogonal
kernel, the packet contains the full exact identity

```
det(T) = +/-1
T NS T^t = U + (-H).
```

It also retains integral unimodular isometries from both the class-1 census
representative and **embedding 33** into `H`, and the resulting 17 by 19
embedding-33 frame rows in the rational source NS. Stacking those rows beneath
`D,D+O` is again unimodular. The requested complement is therefore transported
exactly; this is **not a claim that the particular common-core identification
with embedding 25 extends to this chosen transport**. No such extra condition
is needed to realize the requested J2 class.

Class 1 has minimum 4, 2,436 minimum vectors and automorphism order 2 in the
census. The isometry checker also rejects an isometry to the known class-6
frame. This realizes one new J2 type on the same rational surface. It does
not classify J1/surface-automorphism orbits or construct a different K3 surface.

## Effectivity, nefness and the rational zero

`D^2=0`, `D.O=1`, and `D.F_old=2`. Primitivity follows already from `D.O=1`.
K3 Riemann--Roch makes one of `D,-D` effective; the latter has negative pairing
with the known nef old fibre, so **D is effective**.

If an irreducible effective curve has negative intersection with effective D,
it is a negative-square component of D, hence a `(-2)` curve and has old fibre
degree at most 2. The old rootless frame excludes vertical `(-2)` components.
The checker exhausts the negative-intersection ellipsoids for old degrees
1 and 2, including both signs and the integrality congruence for the first
U coordinate. The homogenized ellipsoids return respectively 2 and 4 short
vectors; **none yields a negative wall**. Thus D is nef.

The zero is the *existing effective rational section* O, not an abstract root
whose effectivity or rationality remains to be established. A primitive nef
isotropic divisor on this K3 gives the elliptic pencil, and `D.O=1` gives its
rational zero. The equation below makes the pencil explicit over Q.

## Equation, maps and sections

The generic trace is built only from the seventeen source sections. It has
`x(P_w)=Nx/h^2`, `y(P_w)=Ny/h^3`, with `deg h=4`. The existing exact
degree-two-neighbour construction is reused: solve

```
a Nx - b Ny = 0 mod h^2,  deg a <= 7, deg b <= 1.
```

The constraint matrix has rank 8 and kernel dimension 2. Its basis `(a0,b0),
(a1,b1)` gives the rational pencil coordinate

```
u = [a1 (x h^2-Nx) + b1 (y h^3+Ny)]
  / [a0 (x h^2-Nx) + b0 (y h^3+Ny)].
```

The residual chord gives a smooth quartic with a rational point coming from
O. A **pointed birational quartic conversion** supplies the elliptic equation;
the degree-four covariant map to a Jacobian is not used as a birational map.
The [equation packet](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_equation_v1.json)
contains every rational coefficient and both map directions. Their equations
and inverse compositions are checked exactly in the quartic function field.

For practical loading use the
[polynomial parent](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_compact_parent_v1.json):

```
y^2 = x^3 + A(s) x + B(s),  deg A <= 8, deg B <= 12.
```

Its JSON coefficient arrays are exact rationals in ascending degree, with
all seventeen points and their height Gram. The base change and Weierstrass
scaling back to the pointed equation are explicit and replayed. Coefficient
normalization reduces the maximum numerator-plus-denominator bit measure from
11,921 to 1,347. This is a bounded chart improvement, not a minimality claim.

The complete old degree-one section window through height 12 has 213
candidates but spans only rank 16 in the child frame. One generically selected
rational bisection closes the gap. Its source NS class is `(2,2,v)`, where

```
v = (0,0,-1,-1,-1,-2,1,1,0,0,1,-1,-1,2,0,0,1),  v^2=10.
```

The Euclidean residual line through `P_{-v}` constructs this curve exactly.
Its intersection with D is one. Solving its incidence with the new pencil
gives a rational section without using a specialized point. Sixteen old
sections and this bisection have **index one** in H. Every transported point
is checked on the equation and against its source curve. The
[section packet](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_sections_v1.json)
contains the coordinates, source words, integral frame coordinates and full
positive-definite height Gram, of determinant 1092.

Since the new complement is rootless, there are no reducible-fibre height
corrections. Exact divisor projection therefore computes the Shioda Gram.
The seventeen independent rational sections give the lower bound. The
already certified geometric Picard rank 19 and Shioda--Tate give the upper
bound 17. Thus the generic rank is exactly 17, and the displayed group is
saturated in the geometric MW lattice.

## Carrier and strict preflight: keep the gates distinct

A first bounded window of norm-ten sums/differences of two basis sections had
68 candidates and no coset-minimum-ten survivors. This says nothing about the
complete carrier atlas. The separately retained triple window has 646
norm-ten candidates; its first four hash-ordered representatives passing
exact minimum-ten coset checks were frozen before construction.

All four cold trace constructions and Sage-free Euclidean carrier checks
passed. A separate Sage replay checks the trace and both polynomial map
identities. The four recorded in-process cold totals sum to approximately
0.313 seconds; this excludes Sage startup and all parent derivation work and
**is not an end-to-end efficiency claim**. No rational parameter was tested.
These are four exact genus-zero carrier equations, not evidence of a heavy
specialization tail or new strict classes.

The
[strict preflight](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_strict_preflight_v1.json)
binds the new cubic and inherited seventeen-section basis, and runs the
existing fail-closed dispatcher. All four arithmetic gates remain UNKNOWN:

1. nonzero arithmetic strict class constructed modulo the full inherited image;
2. its cover soluble;
3. a new rational point recovered;
4. an independent quotient direction certified.

The current theory task confirmed that the successful strict solvers take
pinned arithmetic-field inputs and do not provide a certified general
parameter-to-class constructor. The complete inherited everywhere-even
half-ideal image, strict character space and ordinary-unramified character
space for this new parent have **not** been computed. In particular,
curve302's zero generic-strict dimension must not be copied here. A genus-zero
NS carrier does not certify an arithmetic strict class.

**Next:** construct a principal-square identity in the correctly marked cubic
algebra, together with a character proving novelty modulo the full inherited
image, then attempt its cover under a declared bound. The matched parameter
panel remains blocked at this arithmetic gate. There is no reason to switch
to another tied frame because class 1 has successfully passed all realization
gates; the remaining gap is the prospective arithmetic constructor.

## Replay and bounds

Producer:
`research/elkies-k3/scripts/realize_x1092_class1.sage`.
Run with `sage -python`, one stage at a time: `marking`, `trace`, `equation`,
`section_plan`, `sections`, `normalize`, `compact`, `carrier_plan`,
`carrier --index 0` through `3`, then `strict_preflight`. These were run under
120-second caps for the lattice/chart stages, 600 seconds for trace/equation/
section compilation and 180 seconds per carrier. No stage needed its cap.
The producer checkpoints locally between stages; the final proof packets
are portable JSON. The two-section empty-window packet is historical input
selection evidence; the active producer generates the triple-window protocol.

Read-only replay, requiring no local producer checkpoints:

```sh
timeout 180 sage -python research/elkies-k3/scripts/verify_x1092_class1_realization.sage
```

The separately written checker reconstructs the source points, lattice
transport, nef test, rational equation, inverse maps, all section transports,
height Gram, compact chart and four carrier identities from JSON. It reuses
Sage/PARI and is not a different-engine or formal proof. The
[manifest](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_manifest_v1.json)
binds the producer, checker and result packets. A deliberately altered frame
Gram is rejected by the checker. Local failed development attempts
and logs remain under `research/artifacts/local/elkies-k3/x1092-class1-realization-v1/`;
only passing final packets support the claims above.

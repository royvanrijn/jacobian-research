# Determinant 1092: a finite-field certificate for bisection non-incidence

The subsequent [uniform reduction theorem](DET1092_UNIFORM_SMOOTH_BISECTION_REDUCTION_2026-09-09.md)
proves all40,917 minimum conics stay smooth at149/151 and supplies a four-chart
remedy for pole failures. The original45-trial protocol and outcomes below
are unchanged; no complete atlas incidence test has been run.
The [Euclidean closure](DET1092_EUCLIDEAN_BISECTION_FORMULA_2026-09-09.md)
subsequently resolves all three rational survivors as nonsplit and closes
the fifteen-orbit panel, retaining these original modular UNKNOWN outcomes.

## Result and limits

**New deduction; independently verified application.** An exact modular gate
can exclude an entire smooth rational bisection translation orbit over a
specified rational parent parameter, without constructing its rational trace
point, rational Riemann--Roch kernel, conic equation or parametrization. A
maximal pole-degree test is essential: it justifies lifting the calculation
from the finite function field. Failed gates remain **UNKNOWN**.

On the frozen twelve-new-orbit panel at302, this excludes nine orbits. Three
remain UNKNOWN. Two of three old conic regressions are also excluded by these
primes; the third,47755, is already excluded by its old5-adic certificate but
survives this particular panel. All45 trial outcomes are retained.

The new exclusions also hold on explicit rational residue neighbourhoods of
302. They are arithmetic non-incidence statements, **not search-height
failures**. They do not locate302 relative to the full40,917-orbit atlas,
construct a seed, or upper-bound the rank of any elliptic fibre. No full-atlas
run, extra prime, point search, parameter sweep, later exceptional point,
production mutation or detached job was used.

## 1. Characteristic-zero geometric input

**Verified application of completed inputs.** Retain the
[completed degree-two census](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md)
and its [rationality theorem](DET1092_RATIONAL_BISECTION_INDEX_AND_ODD_DIVISOR_CONSTRUCTION_2026-09-09.md).
For a saved minimum word of norm10, the parent contains the unique smooth
rational curve

\[
 C_w=2O+4F+\phi(w),\qquad C_w^2=-2,\quad C_w.O=0,\quad C_w.F=2.
\]

The parent has24 irreducible singular fibres and holomorphic Euler
characteristic2. The trace section and residual class satisfy

\[
 P=P_{-w}=O+5F-\phi(w),\qquad P.O=3,\qquad C_w+P=3O+9F. \tag{1}
\]

**Established literature.** The section-class and height identities are
standard elliptic-surface intersection theory; see
[Schuett--Shioda, sections6 and11](https://arxiv.org/pdf/0907.0298).
Here all reducible-fibre correction terms vanish. The census and rationality
proofs are dependencies, not re-proved by the modular checker.

Use the original parameter `t`, with302 at `t=0`, and the short equation

\[
 X=x+b_2/12,\quad Y=y+(a_1x+a_3)/2,\qquad
 Y^2=X^3+A(t)X+B(t),\quad A=-c_4/48,\ B=-c_6/864.
\]

Here `deg A<=8`, `deg B<=12`. With monic finite pole polynomial `h`, (1)
gives

\[
 X(P)=N_x/h^2,\quad Y(P)=N_y/h^3,\qquad
 \deg h\le3,\quad\deg N_x\le10,\quad\deg N_y\le15.       \tag{2}
\]

To include the infinity case correctly, let `d=deg h` and `e=3-d` be the
intersection multiplicity with `O` at infinity. The minimal coordinates
there are `t^-4 X,t^-6 Y`. Their pole orders are at most `2e,3e`.
Thus the bounds on the numerator degrees are `4+2d+2e=10` and
`6+3d+3e=15`. Equality of these bounds is not assumed.

A section of `H0(3O+9F)` is

\[
 f_0(t)+f_1(t)X+f_2(t)Y,
 \qquad (\deg f_0,\deg f_1,\deg f_2)\le(9,5,3).           \tag{3}
\]

Vanishing on `P`, after multiplication by `h^3`, gives a19-by20 matrix:
its columns are the coefficients of

\[
 t^j h^3\ (0\le j\le9),\quad
 t^j N_xh\ (0\le j\le5),\quad
 t^j N_y\ (0\le j\le3),                                 \tag{4}
\]

in degrees0 through18. The divisor (3) contains `P`; its residual divisor
has class `C_w`. Uniqueness of the effective irreducible `(-2)` curve makes
that residual divisor exactly `C_w`, not an arbitrary quadratic cover.

## 2. Maximal-pole lifting lemma

**New deduction from elementary valuation and linear algebra.** Let `p>3`
be a prime with `A,B` coefficientwise `p`-integral and parent discriminant
nonzero modulo `p` as a polynomial. Assume the input generic sections have
the indicated reductions. Compute `-sum w_i S_i` entirely over `Fp(t)`.
Require that it is nonzero and its reduced `X` denominator has degree6.

These tests imply that the monic characteristic-zero `h` in (2) has degree3
and is coefficientwise `p`-integral, and that its `N_x,N_y` are integral too,
with reductions exactly those computed in the finite function field.

### Proof of the pole gate

Use the Gauss valuation on `Q_p(t)`, whose valuation ring has uniformizer
`p` and residue field `Fp(t)`. The integral short equation has unit
discriminant in this valuation ring, so its projective elliptic curve is
smooth and proper. Reduction is a group homomorphism. Therefore the modular
word is the reduction of the rational section `P`, whether or not the
intermediate affine addition formulas have inconvenient denominators.

The modular trace is affine, hence `X(P),Y(P)` belong to the Gauss valuation
ring. Write `X(P)=n/d` in reduced form, with `deg d<=6`, and scale the pair
to have integral coefficients and at least one unit coefficient. Affine
reduction forces the denominator to have a unit coefficient. Its reduction
may initially have factors in common with the reduced numerator.

But the reduced denominator of the quotient has degree6. Since `deg d<=6`,
the denominator's leading coefficient is a unit and there can have been no
positive-degree cancellation. Normalizing it to be monic is therefore
integral and commutes with reduction. In characteristic zero this monic
denominator is `h^2`. Gauss multiplicativity implies that monic `h` is
integral; its reduction is the monic square root of the saved denominator.
Finally `N_x=X(P)h^2` and `N_y=Y(P)h^3` are polynomials in the Gauss valuation
ring, hence coefficientwise integral. This proves the claim.

A lower modular denominator degree does **not** prove that this lifting
fails. It simply loses this certificate and is rejected. This is why the
degree gate must not be silently omitted in a faster implementation.

### Proof of the kernel gate

Require a nonzero19-by19 minor of (4) modulo `p`, and save both its columns
and determinant. The corresponding characteristic-zero matrix is integral
by the pole lemma. Its minor is a unit. Solving against that invertible
block proves that its integral kernel is free of rank one and reduces onto
the full modular kernel. In particular, a saved nonzero modular null vector
lifts, after a unit scaling, to the rational line (3).

This argument needs no rational nullspace solve. A rank drop gives UNKNOWN.
The elementary proof, not the CAS finite-field calculation alone, is what
permits a characteristic-zero conclusion. It is written mathematics, not
formal verification or external review.

## 3. Target obstruction and neighbourhood corollary

**New deduction.** At a `p`-integral target `t0`, also require

\[
 h(t_0)\ne0,\quad\Delta_E(t_0)\ne0,\quad f_2(t_0)\ne0
 \pmod p.                                                   \tag{5}
\]

Eliminate `Y` from (3) and divide out the known trace root:

\[
 (f_0+f_1X)^2-f_2^2(X^3+A X+B)=(X-X(P))\,q_{w,t_0}(X).     \tag{6}
\]

All coefficients here are evaluated at `t0`. The quadratic has integral
coefficients and leading coefficient `-f_2(t0)^2`, a unit. If its reduced
discriminant is a nonzero nonsquare, it has no `Q_p` root: any such root
would be integral by the unit-leading-coefficient valuation argument and
would reduce to an `Fp` root. Consequently `C_w` has no rational point over
`t0`. No point at `O` is missed: `C_w.O=0`; also the specialized line has
nonzero `Y` coefficient. Rational generic translations preserve `t0` and
act as automorphisms of its smooth elliptic fibre. The **entire translation
orbit** is excluded.

All conditions and the quadratic reduction remain unchanged on

\[
 t\in t_0+p\mathbf Z_p.                                    \tag{7}
\]

Thus a single certificate excludes every rational parameter in that
residue ball, not only the tested parameter. No continuity estimate or
additional parameter sampling is needed. This is a cover-specific local
obstruction, not a criterion for dependence of a split point or a claim
about Selmer/Sha classes.

## 4. Frozen panel and independent replay

**Verified application.** The preconstruction rule takes the first twelve
`category=rational` rows by numeric saved mask, plus exactly three old map
regressions8044,47755,103186. The primes are149,151,157, with no replacements.
The word/prime choices use generic inputs only;302 is an explicit
retrospective target for incidence evaluation. Numeric-mask selection is
coordinate dependent, not a basis-invariant scheduler. The lemma itself
concerns the intrinsic curve, section and pole divisor.

One certificate for each new excluded orbit is:

| Orbit | Prime | Reduced discriminant |
| ---: | ---: | ---: |
|87|151|135|
|103|149|137|
|109|149|41|
|110|149|56|
|115|149|43|
|117|151|63|
|121|149|12|
|122|151|71|
|123|149|93|

For example, orbit87 at151 has residual polynomial

\[
 63X^2+112X+81,\qquad\operatorname{disc}=135\pmod{151},
 \qquad135^{75}\equiv-1\pmod{151}.
\]

The [orbit87 record](../../artifacts/generated-results/elliptic-curves/det1092_modular_bisection_direct_v1/orbit-87.json)
also supplies the generic word, degree-three pole polynomial, two numerators,
20 kernel coefficients and a19-column minor of determinant64 at151. These
make the small discriminant assertion a verifiable certificate for the
actual characteristic-zero curve rather than an unattached square test.

The new UNKNOWN masks are61,107,111. The three old regression orbits give
two exclusions at the frozen primes:8044 and103186. Orbit47755 is UNKNOWN
here despite its old exact obstruction at5. No retuning was done to hide
that conservative failure.

Across all45 trials the statuses are:

- 20 nonsquare exclusions;
- 20 split or repeated reductions, all UNKNOWN over `Q`;
- 4 maximal-pole-gate failures, UNKNOWN;
- 1 bad target/trace gate, UNKNOWN.

The table and the two successful regression exclusions use only149 and151.
Therefore all these eleven orbits miss every rational parameter
`t=22499*a/b` with integers `a,b`, `b!=0`, `gcd(b,22499)=1`. This explicit
neighbourhood statement is a corollary of (7), not a new parameter campaign.

The first implementation constructed rational traces but not rational RR
kernels; its internal runtime was10.833 seconds. The direct finite-field
implementation used the same panel and took0.291 seconds internally. These
are small-sample timings, not a forecast for a complete atlas calculation.
Each process was capped at25 seconds.

The independent checker uses reversed-order, manually implemented elliptic
addition instead of `EllipticCurve` group arithmetic. It verifies the saved
nonzero minors and kernel equations without solving a nullspace, divides
the cubic by the trace factor, and tests every residue for quadratic roots.
All45 exposures replay. The40 admissible modular RR outcomes also match the
rational-trace prototype; eight nonzero old conic squareclasses agree.
Two successive runs give the same report. An initial development run failed
on coercing rational coefficient strings directly into `Fp`; explicit `QQ`
parsing fixed that checker-only input conversion before its report was sealed.

## 5. Consequence for theorem hunting

**Established application / explanatory boundary.** All40,917 conics being
rational says there are many rational source curves on the common parent.
It does not say that any particular fibre, including302, meets them
rationally. The new certificates establish arithmetic nonsplitting for
nine further complete orbits. Searching taller translates within those
orbits cannot change this. This statement does not exclude general pointed
quartics or the known positive-genus seed carrier.

The [complete-atlas nonuniversality theorem](DET1092_NORM8_DEGENERATION_AND_SMOOTH_ATLAS_OBSTRUCTION_2026-09-09.md)
already shows that smooth rational bisections cannot account for all rank
jumps. The more specific question whether302 lies outside *all* of them is
still UNKNOWN. This prototype supplies a fail-closed way to approach that
question, not its answer. A complete run would be a separately scoped,
checkpointed proof campaign; none has been launched. Survivors would still
require exact incidence decisions, and a split point would require an
independence certificate.

## 6. Evidence and reproduction

Immutable packets:

- [Direct finite-field protocol](../../artifacts/generated-results/elliptic-curves/det1092_modular_bisection_direct_v1/protocol.json),
  [complete outcomes](../../artifacts/generated-results/elliptic-curves/det1092_modular_bisection_direct_v1/summary.json),
  [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_modular_bisection_direct_v1/independent-replay.json).
- [Rational-trace prototype](../../artifacts/generated-results/elliptic-curves/det1092_modular_bisection_exclusion_v1/summary.json),
  retained as a distinct implementation and regression, not overwritten.

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_modular_bisection_direct.sage
```

Checker SHA256:
`c13500b14be30eb60804e3ca49db49ba44d7a9d456854855e8fbd9c66f5289e7`.
Independent report SHA256:
`1c08e334d1766c49ab5456d0f307397efb7acdecdec325fbf53f4eff18c36946`.
No conjectural seed-incidence discriminator or amplification theorem is
claimed by this packet.

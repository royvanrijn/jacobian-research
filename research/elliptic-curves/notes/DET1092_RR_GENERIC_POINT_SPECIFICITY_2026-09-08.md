# Generic elliptic points can give non-generic RR Jacobian classes

## Exact result

**New deduction, independently verified application.** On all eight frozen
V4 fibres, and on a separate generic-point control on302, a point known to
belong to the original elliptic MW17 gives a rational Jacobian class
independent of the full rank17 generic Picard image on its RR member.

Thus the distinction certified for the first302 unlock is real, but is
**not sufficient to certify that the marked elliptic point is an extra
Mordell--Weil direction**. The control points are exactly generic basis
section0; their classes modulo elliptic MW17 are zero, not merely undetected.

The [completed panel certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_generic_point_controls_v2/panel-replay.json)
contains the symbolic construction, all nine independent replays, the finite
exposure totals and a post-arithmetic snapshot of the V4 terminal certificates.
The eight V4 outcomes are independently verified bounded nulls. They are
**not proofs of rank exactly17** or absence of other exceptional points.

## Frozen source-only assignment

**Verified application.** Keep the same generic RR net `A,B` and its
historically calibrated generic centre. Choose `S0`, the first section in
the pinned generic17 basis, before any control arithmetic. Define

\[
r_0(t)=-\frac{B_0(t)+B_1(t)x_{S_0}(t)+B_2(t)y_{S_0}(t)}
                  {A_0(t)+A_1(t)x_{S_0}(t)+A_2(t)y_{S_0}(t)}.
\]

The explicitly recorded coprime numerator and denominator define a map
`P1_t -> P1_u` of **degree7**. The universal control family is

\[
\boxed{\mathcal C_t:s^2=cq(T;r_0(t),0).}
\]

Here `t` indexes the original elliptic fibre, while `T` is the coordinate
along the RR curve. This is a symbolic rational-function construction, not
an interpolation of the nine cases. The certificate also expands the
rational function `s0(t)` for its marked point

\[
P(t)=(T=t,s=s_0(t)).
\]

With `f_i(T;t)=B_i(T)+r0(t)A_i(T)` and
`m_num=-f1+a1*f2/2`, its formula is

\[
s_0(t)=
\left.
\frac{(2x_{S_0}^{\rm short}+c_x)f_2^2-m_{\rm num}^2}{h^3}
\right|_{T=t}.
\]

The replay verifies the identity `s0(t)^2=c*q(t;r0(t),0)` over `Q(t)` and
the complete rational map back to **exactly `S0(t)`**, not to an unknown
elliptic point. The inherited Jacobian basepoint is still
`P0(t)=C_t intersect O`, provided by `D.O=1`.

No exceptional point, winning chart coordinate, rank score or search outcome
enters this assignment. The coefficient data of the generic net are retained
from the calibrated centre; this does not claim that the original centre
was historically selected blindly. The chosen generic basis marking is
fixed, and no basis-invariance claim is needed for this specificity control.

**Excluded parameters.** The family and maps are used off the poles of the
displayed rational functions and off the sextic discriminant locus. The
finite proof additionally requires nonbranch, finite marked/base points and
the displayed section-divisor denominators to be invertible. All these
conditions are checked exactly at every panel parameter; none was replaced.

## Matched invariant, different elliptic interpretation

**Verified application.** For each member let `H` be the full normalized
restriction image of the original K3 Picard group, and let
`xi=[P(t)-P0(t)]`. The proof uses all17 generic section restrictions and the
fibre class, not only the44 inherited rational intersections. The same
finite-character and theta-characteristic proof as in the
[first-unlock certificate](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md)
gives

\[
\operatorname{rank}H=17,\qquad
\operatorname{rank}\langle H,\xi\rangle=18
\]

on all nine members. The whole Jacobian rank is only bounded below by18.

| Marked point | Class in elliptic fibre modulo MW17 | Jacobian class modulo rational span of `H` | Sha image of its true2-Kummer class |
|---|---|---|---|
| Known first302 unlock, on its retrospective member | Independent extra direction | Independent | Zero |
| Generic section0 on each of the eight V4 controls | Exactly zero | Independent | Zero |
| Generic section0 on302, on its source-only control member | Exactly zero | Independent | Zero |

The two302 members are intentionally different: one is conditioned on the
known unlock; the other follows the same generic-input rule as the controls.
This is a **specificity control**, not a claim to have prospectively selected
the successful302 RR member.

## Complete finite exposure and Selmer boundary

**Verified application.** Every row uses exactly the same64 frozen primes,
without refill after bad reduction or a nonunit divisor.

| Case | Valid blocks | Retained skips | Fake-image ranks before/after | Certified subgroup ranks before/after |
|---|---:|---:|---|---|
| scale-0131232 |38|26|16 /17|17 /18|
| scale-0257585 |39|25|16 /17|17 /18|
| scale-0487239 |42|22|16 /17|17 /18|
| scale-0177036 |39|25|16 /17|17 /18|
| scale-0043332 |42|22|16 /17|17 /18|
| scale-0590501 |41|23|16 /17|17 /18|
| scale-0290097 |41|23|16 /17|17 /18|
| scale-0748009 |41|23|16 /17|17 /18|
| 302 generic-point control |38|26|16 /17|17 /18|

Across576 prime attempts there are361 valid blocks,94 denominator skips,
57 bad-sextic skips and64 nonunit-divisor skips. All are retained.
Independent replay uses norms in actual finite extension fields rather
than the constructor's polynomial modular exponentiation.

**Established interpretation and verified application.** Each Jacobian has
an explicitly certified18-dimensional **rational Kummer subspace** inside
its true2-Selmer group. The marked class is outside the generic-image Kummer
subgroup, but is represented by a rational point. Its Sha image is zero;
the Cassels--Tate pairing restricted to this known rational subspace is zero.
These follow from the rational-point/Selmer/Sha sequence, not from treating
finite local survival as a global point. For the relation between even-degree
fake descent and the true groups, see
[Poonen--Schaefer, sections11--13](https://math.mit.edu/~poonen/papers/descent.pdf).

The full2-Selmer dimensions, additional residual Selmer classes, nonzero Sha
classes and full Cassels--Tate matrices are **not computed**. None is asserted
absent. The controls' elliptic relative Kummer class is exactly zero because
the marked point is `S0`. The earlier first302 elliptic class remains
non-strict modulo MW17 for its frozen bad-place definition. No unidentified
map to an MW16 cubic ideal class, or uncomputed Jacobian strict kernel, is
being substituted for these distinct statements.

## Why this happens, already over the function field

**New deduction.** The marked section and all generic Picard-image classes
exist over `Q(t)`. A relation over that field would specialize to a relation
at the verified smooth first control. Therefore

\[
\operatorname{rank}\operatorname{Jac}(\mathcal C_t)(\mathbf Q(t))\ge18,
\qquad \operatorname{rank}E_t(\mathbf Q(t))=17.
\]

The second equality is the already proved full geometric/arithmetic rank of
the original parent; there is no new elliptic direction in this construction.

The degree7 cover explains the discrepancy. A generic RR member intersects
`S0` in a degree7 divisor. After the base change `u=r0(t)`, one of its points
becomes individually rational. Restricting the whole divisor and marking
one of its points are different operations. Precisely,

\[
\operatorname{Tr}_{\mathbf Q(t)/\mathbf Q(u)}\xi(t)=\Phi(S_0),
\qquad u=r_0(t).
\]

Indeed the seven conjugates are exactly the zeros of
`num(r0)(T)-u*den(r0)(T)`, the previously verified intersection divisor.
Thus `7*xi-Phi(S0)` has trace zero. A new class relative to the constant
Picard image can come from resolving this **generic degree7 intersection**,
without resolving any exceptional elliptic class.

This is the exact obstruction to promoting the Jacobian distinction alone
to an elliptic-extra-point certificate. It does not rule out a more refined
arithmetic discriminator, nor establish exact ranks of the V4 fibres.

## Reproduction and immutable failures

The v1 protocol was frozen before any local test, then stopped on the first
case because Sage represented the exact quotient `A2/h` in a fraction field.
The retained v2 repair only coerces this known polynomial before indexing
its constant coefficient. Independent replay checks that the selector,
cases, primes and limits agree exactly between protocols.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_generic_control_panel.sage
```

Individual constructors and replays took about one second each. Each had
a25-second cap and per-prime checkpoints; the complete panel replay also
has a25-second cap. No point search, class-group calculation, global Selmer
campaign, or V3/V4 modification was performed. V4 outcome certificates were
read only after all control arithmetic had passed.

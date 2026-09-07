# 302 constructor transfer: two bounded attempts

No additional302 strict class was constructed in this checkpoint. The
[positive reference result](CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md)
is unchanged: two point-free classes, two new ideal-class2-torsion factors,
and two rational lifts verified retrospectively. The302 work below uses
only the equation and marked generic sections and supplies bounded
constructor evidence, not another general capacity exclusion.

## Arithmetic units

The unit attempt starts from the unit ideal and the six generic ordinary
unramified half ideals. It applies81 fixed nonnegative archimedean reduction
directions to each. Repeated reduced ideals would provide units through
the ratio of their exact principal multipliers; a nontrivial unit would
then be tested for a strict generic correction.

The [independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_arithmetic_unit_class_verification_v1.json)
checks all567 ideal identities. The81 repeated-ideal pairs give only
units `+/-1`. Thus this finite run supplied no new norm-one unit, and no
additional class. It does not compute the unit group. The fixed attempt
is closed; no larger automatic restart follows from this protocol.

## Explicit rational bisections

The new geometric constructor uses the quartic and its fifteen generic
lines. An intersecting pair of those lines gives a residual plane conic,
which is a generic section. Choose three further generic lines, disjoint
from the first pair, with incidence graph a three-vertex path. There are
266 such configurations, frozen before testing their zero fibres.

A unique smooth quadric containing the conic and those three lines has a
residual curve of bidegree `(2,1)`. When irreducible, this is a rational
twisted cubic. Its two-to-one map to the actual parent parameter
`t=W/L` is computed exactly:

\[
 t=\frac{a s_0^2+b s_0s_1+c s_1^2}
          {d s_0^2+e s_0s_1+f s_1^2}.
\]

Thus this particular bisection fibre splits rationally exactly when

\[
 \Delta_D(t)=(b-et)^2-4(a-dt)(c-ft)\in\mathbf Q^2,
\]

with repeated fibres handled separately. The numerator and denominator
are coprime binary quadratics, so this includes the projective root at
infinity. These are explicit **solubility conditions for specified
bisections**, not rank predictors.

The [construction artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_conic_triple_bisection_v1.json)
records the quadric, twisted-cubic parametrization and quadratic polynomial
Delta for every retained bisection. The
[separate verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_conic_triple_bisection_verification_v1.json)
replays the curve equations, degree-two maps and discriminants.

| Fixed catalogue outcome | Count |
|---|---:|
| Configurations |266|
| Duplicate quadrics |74|
| Distinct quadrics |192|
| Irreducible rational twisted-cubic bisections |178|
| Reducible residuals |14|
| Rationally split bisections at zero |0|

Every Delta_D(0) is a nonzero nonsquare. The178 zero fibres represent114
distinct quadratic fields, none equal to a nonsplit field from the previous
four-line catalogue. Thus this is a different construction, but it still
does not lift rationally at302.

Each of the14 reducible residuals is verified to be a rational line plus
a rational conic, both mapping with degree one to t. They are generic
sections and hence belong to the certified full saturated generic basis.
They do not supply an additional specialization direction.

## Consequence for the next step

The successful reference constructor required new global principal-ideal
dependencies. Neither the fixed unit route nor this rational-bisection
catalogue supplied one on302. Increasing their budgets is not justified by
the present output. The additional302 class and its class-creation event
remain the primary missing objects.

These computations do not explain the large jump by visibility, prove that
other covers fail, or exclude fibre-specific class-group mechanisms. The
new quadratic conditions are retained as exact retrospective controls. No
exceptional points, active search charts or new parameters were inputs.

The [checkpoint archive](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_transfer_attempts_evidence_v1.zip)
and [manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_transfer_attempts_evidence_v1.json)
preserve the terminal calculations. Restore into an empty replay checkout.
Cheap checks are:

```sh
sage -python elliptic-curves/rank-jump/verify_curve302_unit_attempt.py check
sage -python elliptic-curves/rank-jump/verify_curve302_conic_triple_bisection.py check
```

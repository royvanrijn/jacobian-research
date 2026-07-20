# Automatically clustered obstruction report

## One-divisor pole matrices

- degree `0`: rank `56`, earliest pole `(1, 0)`, representatives `['0']`
- degree `1`: rank `56`, earliest pole `(1, 0)`, representatives `['v']`
- degree `2`: rank `60`, earliest pole `(1, 0)`, representatives `['v**2', 'v**2 + v']`
- degree `3`: rank `70`, earliest pole `(1, 0)`, representatives `['v**3', 'v**3 + v**2 + v']`
- degree `4`: rank `80`, earliest pole `(1, 0)`, representatives `['v**4', 'v**4 + v**3 + v**2 + v']`

No tested prime changed these rational ranks.

## Collision-normalized F4 families

- `triangle_d3` (10 terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.
- `triangle_d4` (15 terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.
- `x_heavy_15` (15 terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.
- `y_heavy_15` (15 terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.
- `balanced_15` (15 terms/coordinate): exact `Q` basis `[1]`; all saturated modular runs empty.

## Two-divisor toric families

- chart `['y/x', 'x**2/y']`, 10 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.
- chart `['y/x', 'x**2/y']`, 15 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.
- chart `['y/x', 'x/y**2']`, 10 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.
- chart `['y/x', 'x/y**2']`, 15 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.
- chart `['y/x**2', 'x/y']`, 10 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.
- chart `['y/x**2', 'x/y']`, 15 terms/coordinate: exact `Q` basis `[1]`; three modular runs empty.

## Translated two-divisor families

- `diamond_r2`: support 13, cancellation kernel 3, pole rank 10; exact `Q` basis `[1]` and three modular runs empty.
- `box_r2`: support 25, cancellation kernel 6, pole rank 19; exact `Q` basis `[1]` and three modular runs empty.
- `asymmetric_30`: support 30, cancellation kernel 8, pole rank 22; exact `Q` basis `[1]` and three modular runs empty.
- `shifted_42`: support 42, cancellation kernel 9, pole rank 33; exact `Q` basis `[1]` and three modular runs empty.
- `box_r3`: support 49, cancellation kernel 10, pole rank 39; exact `Q` basis `[1]` and three modular runs empty.

## Newton-directed unequal-coordinate families

- `total_ratio_2_3`: kernel dimensions 18/34; recorded modular F4 timeout; no exact emptiness claim.
- `weighted_3_2_ratio_2_3`: kernel dimensions 13/21; exact collision-normalized `Q` basis `[1]`; three modular runs empty.
- `weighted_2_3_ratio_2_3`: kernel dimensions 13/21; exact collision-normalized `Q` basis `[1]`; three modular runs empty.
- `corner_2_7_square_cube`: kernel dimensions 11/13; exact collision-normalized `Q` basis `[1]`; three modular runs empty.
- `corner_7_2_square_cube`: kernel dimensions 11/13; exact collision-normalized `Q` basis `[1]`; three modular runs empty.

## Recurring reason

Pole cancellation alone has a large kernel and is not the obstruction. In every exactly completed sparse family, collision normalization turns the ideal into the unit ideal; all Stage D kernels still contain the inverse-chart Keller control. The dense total-degree Stage D control is explicitly unclassified after a scaling timeout.

## Reduced degree-three contradiction core

The original 23 equations reduce greedily to 19. The constant Jacobian coefficient and `Py(0)=0` are redundant because the retained derivative normalization fixes the constant term. Retained labels:

`jac_x4_y0`, `jac_x3_y0`, `jac_x2_y1`, `jac_x2_y0`, `jac_x1_y3`, `jac_x1_y2`, `jac_x1_y1`, `jac_x1_y0`, `jac_x0_y4`, `jac_x0_y3`, `jac_x0_y2`, `jac_x0_y1`, `P_at_0`, `Q_at_0`, `P_at_1`, `Q_at_1`, `Px_at_0_minus_1`, `Qx_at_0`, `Qy_at_0_minus_1`.

The retained system is exactly empty over `Q`. It is greedily reduced, not claimed inclusion-minimal, because single-deletion ideals become costly positive-dimensional systems.

## Published `(9,27)` regression

The nine-equation Section 5 elimination reproduces equation (5.9) over `Q`: terminal relation found = `True`. The subsequent valuation split gives the published order contradiction `70/66/72` and degree contradiction `76/78/75`.
This assumes the paper's prior reductions; it is not an independent re-proof of every cited proposition.

## Reduced `(72,108)` target

- broad `case_with_vertical_vertices`: 186 variables, 302 equations; not launched.
- broad `case_without_vertical_vertices`: 72 variables, 92 equations; not launched.
- experimental edge-power `case_with_vertical_vertices`: 163 variables, 278 equations; not launched.
- experimental edge-power `case_without_vertical_vertices`: 49 variables, 88 equations; timeout `30s` (complexity only).

No `(72,108)` survivor or exclusion has been obtained. Original affine collision constraints are not valid in these reduced coordinates unless explicitly transported.

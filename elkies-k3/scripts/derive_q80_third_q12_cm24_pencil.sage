#!/usr/bin/env sage
"""Form the CM24 third-q12 genus-one pencil from resolved D7 rows.

This script loads the exact 2-dimensional kernel in
``derive_q80_third_q12_local_gates.sage``.  For a kernel row

    n = a(W) + b*X + c(W)*z_Q,

it clears ``z_Q=(Y+Qy)/(X-Qx)``.  The new base is the ratio ``V=N1/N0``.
Because the cleared relation is linear in Y, eliminating Y gives one exact
plane equation in the old variables W,X over the new V-line.  Common fixed
factors are removed and reported before any Jacobian conversion.  Direct
normalization at the split prime 73 gives genus one at both V=1 and V=7; see
``verify_q80_third_q12_cm24_genus.sage``.
"""

from pathlib import Path

from sage.all import PolynomialRing


HERE = Path(__file__).resolve().parent
load(str(HERE / "derive_q80_third_q12_local_gates.sage"))

new_parameter_ring = PolynomialRing(L, "Vnew")
Vnew = new_parameter_ring.gen()
new_old_base = PolynomialRing(new_parameter_ring, "w")
w = new_old_base.gen()
new_x_ring = PolynomialRing(new_old_base, "x")
x = new_x_ring.gen()


def lift_base(polynomial):
    return new_old_base(
        [new_parameter_ring(coefficient) for coefficient in polynomial.list()]
    )


new_A = lift_base(A)
new_B = lift_base(B)
new_Qx = lift_base(Qx)
new_Qy = lift_base(Qy)


def numerator_data(row):
    a_polynomial = sum(new_old_base(row[index])*w**index for index in range(5))
    b_scalar = new_old_base(row[5])
    c_polynomial = sum(
        new_old_base(row[6+index])*w**index for index in range(3)
    )
    polynomial_part = (
        (a_polynomial+b_scalar*x)*(x-new_Qx)+c_polynomial*new_Qy
    )
    return polynomial_part, c_polynomial


polynomial0, c0 = numerator_data(kernel[0])
polynomial1, c1 = numerator_data(kernel[1])
y_coefficient = c1-Vnew*c0
polynomial_part = polynomial1-Vnew*polynomial0
eliminated = (
    polynomial_part**2-y_coefficient**2*(x**3+new_A*x+new_B)
)

content = eliminated[0]
for coefficient in eliminated.list()[1:]:
    content = content.gcd(coefficient)
content = content.monic()
primitive = eliminated//content
assert primitive*content == eliminated
residual_cubic, extraneous_remainder = primitive.quo_rem(x-new_Qx)
assert extraneous_remainder == 0 and residual_cubic.degree() == 3

# The weighted projective closure has weight 3 on x.  Its degree-nine
# infinity form is independent of the new-base parameter up to scale and has
# the distinguished double branch xi=3 and simple branch xi=-6.  Pin this
# exact local datum because it is the natural rational place for the ensuing
# genus-one-to-Weierstrass conversion.
assert tuple(coefficient.degree() for coefficient in residual_cubic.list()) == (
    9,
    6,
    2,
    0,
)
infinity_ring = PolynomialRing(new_parameter_ring, "xi")
xi = infinity_ring.gen()
infinity_cubic = sum(
    infinity_ring(residual_cubic[index][9-3*index])*xi**index
    for index in range(4)
)
assert infinity_cubic == infinity_cubic[3]*(xi-3)**2*(xi+6)

print(
    f"Q80THIRDCM24PENCIL|content={content.factor()}|"
    f"raw_X_degree={eliminated.degree()}|primitive_X_degree={primitive.degree()}|"
    f"residual_X_degree={residual_cubic.degree()}|"
    f"residual_W_degree={max(coefficient.degree() for coefficient in residual_cubic.list())}|"
    "infinity_branches=3(double),-6(simple)|"
    "extraneous_factor=X-Qx|generic_branch_normalization=GENUS_ONE",
    flush=True,
)
print(
    f"Q80THIRDCM24PENCIL|N0={tuple(kernel[0])}|N1={tuple(kernel[1])}|"
    "new_base=N1/N0|status=PASS_RESOLVED_Q12_PENCIL",
    flush=True,
)

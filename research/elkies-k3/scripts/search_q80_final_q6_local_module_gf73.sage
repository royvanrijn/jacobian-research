#!/usr/bin/env sage
"""Solve the corrected final-q6 5D-to-2D local module over GF(73).

For the true specialized horizontal section S, put L=O+S+2F.  The ambient
space H0(L) has the Smith-saturated basis ``1,s,s^2,qsat,s*qsat``, where
q0 is the residue-corrected chord and

  qsat=(q0+63)/(s-27).

The Smith invariants of the numerator module are ``(1,s-27)``; the older
``q0/(s-2),q0`` trial basis was not saturated.  The desired divisor is

  D=L-R_A1-(R5+R8)_A3-(R11+R13+R16)_D4.

Each disjoint fiber therefore contributes one evaluation row.  We recover
the possible row values directly from the finite I2/I4 and infinite I0*
limits of the chord branch, take every 3x5 kernel, and test the induced
fractional-linear q0 gauge exactly.  The lattice target is
``A1+2A3+2A4``, root data ``(15,66,800)``.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, FunctionField, PolynomialRing, gcd, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--section-index", type=int, default=0)
parser.add_argument("--arbitrary-a3-row", action="store_true")
parser.add_argument("--arbitrary-infinity-row", action="store_true")
parser.add_argument("--i2-index", type=int, choices=(0, 1, 2))
parser.add_argument("--i4-index", type=int, choices=(0, 1, 2))
parser.add_argument("--i4-value-index", type=int)
parser.add_argument("--write-artifact", action="store_true")
parser.add_argument(
    "--i4-gate",
    choices=("value", "derivative", "bzero", "azero"),
    default="value",
)
arguments = parser.parse_args()
FIFTH = ROOT / "artifacts/generated-results/q80-deforming-fifth-pair23-gf73.json"
Q6 = ROOT / "artifacts/generated-results/q80-deforming-fifth-q6-horizontal-candidates-gf73.json"
KNOWN_HASHES = {
    FIFTH: "23fc49bce2618a6d3c5f5e18ded34b4ffbee220be83523ae250bf7774a91db14",
    Q6: "fcd61f89daab0a68785a006e6b10dc3829b1c30c24243b67a4e1b80c7d6e6e09",
}
payloads = {}
for path, expected in KNOWN_HASHES.items():
    content = path.read_bytes()
    assert hashlib.sha256(content).hexdigest() == expected
    payloads[path] = json.loads(content)
fifth = payloads[FIFTH]
section = payloads[Q6]["q6_candidates"][arguments.section_index]
assert section["I2_labels"] == [0, 0, 0]
assert section["I0star_correction"] == "1"

finite = GF(73, impl="modn")
parameter = FunctionField(finite, "R")
R = parameter.gen()
s_ring = PolynomialRing(parameter, "s")
s = s_ring.gen()
s_field = s_ring.fraction_field()
A = s_ring(fifth["A_coefficients_low_to_high"])
B = s_ring(fifth["B_coefficients_low_to_high"])


def polynomial(coefficients):
    return s_ring([parameter(value) for value in coefficients])


X = s_field(polynomial(section["X_numerator_coefficients_low_to_high"]))/s_field(
    polynomial(section["X_denominator_coefficients_low_to_high"])
)
Y = s_field(polynomial(section["Y_numerator_coefficients_low_to_high"]))/s_field(
    polynomial(section["Y_denominator_coefficients_low_to_high"])
)
assert Y**2 == X**3+s_field(A)*X+s_field(B)


def branch_for_raw_chord(raw_chord):
    return raw_chord**4-6*X*raw_chord**2-8*Y*raw_chord-3*X**2-4*s_field(A)


def radical_degree(value):
    numerator = value.numerator()
    return int((numerator//gcd(numerator, numerator.derivative())).degree())


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return n-1, n*(n-1), n, n, f"I{n}"
    if ord_delta == 2:
        return 0, 0, 1, 2, "II"
    if ord_delta == 3:
        return 1, 2, 2, 3, "III"
    if ord_delta == 4:
        return 2, 6, 3, 4, "IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, 24, 4, 6, "I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta-6)
        rank = n+4
        return rank, 2*rank*(rank-1), 4, n+6, f"I{n}*"
    if ord_delta == 8:
        return 6, 72, 3, 8, "IV*"
    if ord_delta == 9:
        return 7, 126, 2, 9, "III*"
    if ord_delta == 10:
        return 8, 240, 1, 10, "II*"
    raise ArithmeticError((ord_a, ord_b, ord_delta))


def classify_quartic(quartic, twist):
    """Return the unit-corrected Jacobian and placewise minimal root data."""
    coefficients = list(quartic.list())+[parameter(0)]*5
    e, d_coefficient, c, b, a = coefficients[:5]
    invariant_I = 12*a*e-3*b*d_coefficient+c**2
    invariant_J = (
        72*a*c*e+9*b*c*d_coefficient-27*a*d_coefficient**2
        -27*b**2*e-2*c**3
    )
    jacobian_A = twist**2*(-27*invariant_I)
    jacobian_B = twist**3*(-27*invariant_J)
    delta_core = twist**6*(4*invariant_I**3-invariant_J**2)
    root_rank = root_count = 0
    root_determinant = 1
    euler = 0
    finite_signature = []
    finite_scalings = []
    minimalizing_unit = parameter(1)
    finite_factors = set()
    for value in (jacobian_A, jacobian_B, delta_core):
        for polynomial_value in (value.numerator(), value.denominator()):
            finite_factors.update(factor for factor, _ in polynomial_value.factor())
    for factor in sorted(finite_factors, key=str):
        raw_orders = (
            int(jacobian_A.valuation(factor)),
            int(jacobian_B.valuation(factor)),
            int(delta_core.valuation(factor)),
        )
        scaling = min(
            raw_orders[0]//4, raw_orders[1]//6, raw_orders[2]//12
        )
        orders = (
            raw_orders[0]-4*scaling,
            raw_orders[1]-6*scaling,
            raw_orders[2]-12*scaling,
        )
        if scaling:
            finite_scalings.append((str(factor), scaling))
            minimalizing_unit *= parameter(factor)**(-scaling)
        if orders[2] == 0:
            continue
        data = kodaira_data(*orders)
        rank, count, determinant, local_euler, kind = data
        degree = int(factor.degree())
        root_rank += degree*rank
        root_count += degree*count
        root_determinant *= determinant**degree
        euler += degree*local_euler
        finite_signature.append(
            (str(factor), degree, list(raw_orders), scaling, list(orders), kind)
        )
    raw_infinity_orders = tuple(
        int(value.denominator().degree()-value.numerator().degree())
        for value in (jacobian_A, jacobian_B, delta_core)
    )
    infinity_scaling = min(
        raw_infinity_orders[0]//4,
        raw_infinity_orders[1]//6,
        raw_infinity_orders[2]//12,
    )
    infinity_orders = (
        raw_infinity_orders[0]-4*infinity_scaling,
        raw_infinity_orders[1]-6*infinity_scaling,
        raw_infinity_orders[2]-12*infinity_scaling,
    )
    infinity_kind = "smooth"
    if infinity_orders[2] > 0:
        rank, count, determinant, local_euler, infinity_kind = kodaira_data(
            *infinity_orders
        )
        root_rank += rank
        root_count += count
        root_determinant *= determinant
        euler += local_euler
    assert euler == 24
    minimal_A = jacobian_A*minimalizing_unit**4
    minimal_B = jacobian_B*minimalizing_unit**6
    minimal_delta = delta_core*minimalizing_unit**12
    assert minimal_A.denominator() == 1
    assert minimal_B.denominator() == 1
    assert minimal_delta.denominator() == 1
    minimal_A = minimal_A.numerator()
    minimal_B = minimal_B.numerator()
    minimal_delta = minimal_delta.numerator()
    assert minimal_A.degree() <= 8
    assert minimal_B.degree() <= 12
    assert minimal_delta.degree() <= 24
    return {
        "status": "classified",
        "finite_scalings": [list(row) for row in finite_scalings],
        "finite_signature": [list(row) for row in finite_signature],
        "raw_infinity_orders": list(raw_infinity_orders),
        "infinity_scaling": infinity_scaling,
        "infinity_orders": list(infinity_orders),
        "infinity_fiber": infinity_kind,
        "root_data": [root_rank, root_count, root_determinant],
        "CM24_MW_rank": 18-root_rank,
        "A_rational": str(jacobian_A),
        "B_rational": str(jacobian_B),
        "Delta_rational": str(delta_core),
        "A_coefficients_low_to_high": list(map(int, minimal_A.list())),
        "B_coefficients_low_to_high": list(map(int, minimal_B.list())),
        "Delta_coefficients_low_to_high": list(map(int, minimal_delta.list())),
    }


# Cancel the unique finite pole and then the s^2,s principal part at infinity.
x_denominator_factors = tuple(X.denominator().factor())
assert len(x_denominator_factors) == 1
pole_factor, pole_exponent = x_denominator_factors[0]
assert pole_factor.degree() == 1 and pole_exponent == 2
pole_center = finite(-pole_factor[0]/pole_factor[1])
linear = s_field(s-pole_center)
best = None
for residue in finite:
    raw = s_field(R)+s_field(residue)/linear
    degree = radical_degree(branch_for_raw_chord(raw))
    row = (degree, int(residue), (), raw)
    if best is None or degree < best[0]:
        best = row
for power in (2, 1):
    refined = None
    for coefficient in finite:
        raw = best[3]+s_field(coefficient)*s**power
        degree = radical_degree(branch_for_raw_chord(raw))
        row = (degree, best[1], best[2]+(int(coefficient),), raw)
        if refined is None or degree < refined[0]:
            refined = row
    best = refined
ambient_branch_degree, pole_residue, polynomial_part, raw_chord = best
h = raw_chord-s_field(R)
assert len(polynomial_part) == 2

# Exact Smith saturation certificate for the numerator module generated by
# 1 and q0.  The first invariant is one and the gcd of the 2x2 minors is H;
# adjoining (q0+63)/H saturates the missing factor.
H = s_ring(s-pole_center)
assert X.denominator() == H**2 and Y.denominator() == H**3
Nx = X.numerator()
Ny = Y.numerator()
K, K_remainder = (Ny+finite(pole_residue)*Nx).quo_rem(H)
assert K_remainder == 0
smith_numerator_module = matrix(s_ring, [
    [-Nx, K],
    [H**2, -finite(pole_residue)*H],
    [0, H**2],
])
assert gcd(smith_numerator_module.list()).monic() == 1
assert gcd(smith_numerator_module.minors(2)).monic() == H.monic()
saturated_numerator_column = []
for row_index in range(3):
    quotient, remainder = (
        smith_numerator_module[row_index, 1]
        + finite(63)*smith_numerator_module[row_index, 0]
    ).quo_rem(H)
    assert remainder == 0
    saturated_numerator_column.append(quotient)

# Rebuild the translated quartic over GF(73)(s), with q the ambient chord q0.
old_base = FunctionField(finite, "t")
t = old_base.gen()
q_ring = PolynomialRing(old_base, "q")
q = q_ring.gen()


def to_old_base(value):
    numerator = value.numerator()
    denominator = value.denominator()
    numerator_t = sum(finite(c)*t**i for i, c in enumerate(numerator.list()))
    denominator_t = sum(finite(c)*t**i for i, c in enumerate(denominator.list()))
    return old_base(numerator_t/denominator_t)


X_t = to_old_base(X)
Y_t = to_old_base(Y)
A_t = to_old_base(s_field(A))
h_t = to_old_base(h)
translated = (
    (q+h_t)**4-6*X_t*(q+h_t)**2-8*Y_t*(q+h_t)
    -3*X_t**2-4*A_t
)


def evaluate_at(value, point):
    return value.numerator()(point)/value.denominator()(point)


delta = -finite(16)*(4*A**3+27*B**2)
i2_roots = tuple(
    finite(-factor[0]/factor[1])
    for factor, exponent in delta.factor() if exponent == 2
)
i4_roots = tuple(
    finite(-factor[0]/factor[1])
    for factor, exponent in delta.factor() if exponent == 4
)
assert len(i2_roots) == len(i4_roots) == 3


def finite_branch_values(point):
    specialized = PolynomialRing(finite, "z")([
        evaluate_at(coefficient, point) for coefficient in translated.list()
    ])
    values = []
    for factor, _ in specialized.factor():
        if factor.degree() == 1:
            values.append(finite(-factor[0]/factor[1]))
    return tuple(sorted(set(values), key=int))


def finite_branch_factorization(point):
    specialized = PolynomialRing(finite, "z")([
        evaluate_at(coefficient, point) for coefficient in translated.list()
    ])
    return tuple(
        (int(-factor[0]/factor[1]), int(exponent))
        for factor, exponent in specialized.factor()
        if factor.degree() == 1
    )


i2_values = {int(root): finite_branch_values(root) for root in i2_roots}
i4_values = {int(root): finite_branch_values(root) for root in i4_roots}
i4_factorizations = {
    int(root): finite_branch_factorization(root) for root in i4_roots
}

# First infinitely-near q0 jets at each I4 node.  The endpoint subtraction
# can impose the derivative of a+b*q0 rather than its value on the contracted
# nodal model, so retain both exact possibilities.
common_denominator = lcm([
    coefficient.denominator() for coefficient in translated.list()
])
local_ring = PolynomialRing(finite, names=("eps", "jet"))
eps, jet = local_ring.gens()


def local_jets(point, value):
    expression = local_ring.zero()
    for q_degree, coefficient in enumerate(translated.list()):
        scaled = coefficient*common_denominator
        assert scaled.denominator() == 1
        polynomial_scaled = scaled.numerator()
        local_coefficient = sum(
            finite(c)*(finite(point)+eps)**degree
            for degree, c in enumerate(polynomial_scaled.list())
        )
        expression += local_coefficient*(finite(value)+jet*eps)**q_degree
    minimum_eps = min(exponents[0] for exponents in expression.dict())
    leading_jet = PolynomialRing(finite, "j")([
        sum(
            coefficient
            for (eps_degree, jet_degree), coefficient in expression.dict().items()
            if eps_degree == minimum_eps and jet_degree == degree
        )
        for degree in range(expression.degree(jet)+1)
    ])
    return tuple(sorted({root for root, _ in leading_jet.roots()}, key=int))


i4_jets = {
    int(root): {
        int(value): local_jets(root, value)
        for value in i4_values[int(root)]
    }
    for root in i4_roots
}

# On an I2 nonidentity component the old coordinates contract to the nodal
# point.  Since S is on the identity component, q0 has the constant value
# (Y_S)/(x_node-X_S)-h there; this is the actual A1 restriction row.
x_ring = PolynomialRing(finite, "x")
x_node_variable = x_ring.gen()


def evaluate_surface(value, point):
    return finite(value.numerator()(point)/value.denominator()(point))


def node_x(point):
    cubic = (
        x_node_variable**3+finite(A(point))*x_node_variable+finite(B(point))
    )
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return finite(-common[0]/common[1])


i2_component_values = {}
for root in i2_roots:
    node = node_x(root)
    value = (
        evaluate_surface(Y, root)/(node-evaluate_surface(X, root))
        - evaluate_surface(h, root)
    )
    i2_component_values[int(root)] = (finite(value),)

# At I0*, q0 has half-integral leading order.  In the ramified chart
# z=1/s=tau^2, q0=tau^-3*w and the leading quartic is
# w^4+60*w^2+65.  Smith saturation compensates this at s=27; in the basis
# (1,s,s^2,qsat,s*qsat) the transported outer-component gate is a2=0.
coefficient_poles = []
for degree, coefficient in enumerate(translated.list()):
    pole = coefficient.numerator().degree()-coefficient.denominator().degree()
    coefficient_poles.append((degree, pole))
assert tuple(coefficient_poles) == ((0, 6), (1, 4), (2, 3), (3, -1), (4, 0))
w_ring = PolynomialRing(finite, "w")
w = w_ring.gen()
infinity_polynomial = w**4+60*w**2+65
assert infinity_polynomial == (w+5)*(w+34)*(w+39)*(w+68)


def kernel_function(row):
    a0, a1, a2, b0, b1 = map(finite, row)
    qsat_coefficient = s_field(b0+b1*s)/s_field(s-27)
    return (
        s_field(a0+a1*s+a2*s**2)+63*qsat_coefficient,
        qsat_coefficient,
    )


# A generic FunctionField(GF(73))(s) substitution is much more expensive than
# the same calculation after specializing the prospective new base R.  Any
# square factor over GF(73)(R) remains square at a good specialization, so a
# specialized radical degree above eight safely excludes an odd quartic.  Do
# this before constructing the symbolic cover; it makes the exhaustive P^2
# quotient search practical.
fast_ring = PolynomialRing(finite, "u")
u = fast_ring.gen()
fast_field = fast_ring.fraction_field()


def to_fast(value):
    numerator = value.numerator()
    denominator = value.denominator()
    return fast_field(
        fast_ring([finite(coefficient) for coefficient in numerator.list()])
        / fast_ring([finite(coefficient) for coefficient in denominator.list()])
    )


X_fast = to_fast(X)
Y_fast = to_fast(Y)
A_fast = fast_field(fast_ring([finite(c) for c in A.list()]))
h_fast = to_fast(h)


def fast_kernel_function(row):
    a0, a1, a2, b0, b1 = map(finite, row)
    qsat_coefficient = fast_field(b0+b1*u)/fast_field(u-27)
    return (
        fast_field(a0+a1*u+a2*u**2)+63*qsat_coefficient,
        qsat_coefficient,
    )


def fast_kernel_possible(kernel):
    (a0, b0), (a1, b1) = map(fast_kernel_function, kernel.rows())
    for new_base_value in (2, 3, 5, 7, 11):
        denominator = finite(new_base_value)*b0-b1
        if denominator == 0:
            continue
        q_value = (a1-finite(new_base_value)*a0)/denominator+h_fast
        cover = (
            q_value**4-6*X_fast*q_value**2-8*Y_fast*q_value
            -3*X_fast**2-4*A_fast
        )
        square_class = cover.numerator()*cover.denominator()
        radical = square_class//gcd(square_class, square_class.derivative())
        return int(radical.degree()) <= 8
    return False


def cheap_radical_degree(polynomial_value):
    specialized_ring = PolynomialRing(finite, "z")
    for value in (2, 3, 5, 7, 11):
        coefficients = []
        for coefficient in polynomial_value.list():
            if coefficient.denominator()(finite(value)) == 0:
                break
            coefficients.append(
                coefficient.numerator()(finite(value))
                / coefficient.denominator()(finite(value))
            )
        else:
            specialized = specialized_ring(coefficients)
            radical = specialized//gcd(specialized, specialized.derivative())
            return int(radical.degree())
    return None


def cover_factorization_for_kernel(kernel):
    (a0, b0), (a1, b1) = map(kernel_function, kernel.rows())
    # If this determinant vanishes then the quotient is a rational function
    # of the old base s alone.  Cross multiplication can otherwise leave a
    # spurious genus-one fixed factor (the a2=0 candidate does exactly this).
    if a0*b1-a1*b0 == 0:
        return None
    if not fast_kernel_possible(kernel):
        return None
    denominator = s_field(R)*b0-b1
    if denominator == 0:
        return None
    q0 = (a1-s_field(R)*a0)/denominator
    cover = branch_for_raw_chord(q0+h)
    square_class = cover.numerator()*cover.denominator()
    cheap_degree = cheap_radical_degree(square_class)
    if cheap_degree is not None and cheap_degree > 8:
        return None
    factorization = square_class.factor()
    odd_degree = sum(
        int(factor.degree())
        for factor, exponent in factorization
        if int(exponent) % 2
    )
    if odd_degree != 4:
        return None
    return factorization


if arguments.arbitrary_infinity_row:
    chosen_i2_indices = (
        range(3) if arguments.i2_index is None else (arguments.i2_index,)
    )
    chosen_i4_indices = (
        range(3) if arguments.i4_index is None else (arguments.i4_index,)
    )
    infinity_hits = []
    infinity_tests = 0
    old_base_collapses = 0
    for i2_index in chosen_i2_indices:
        i2_root = i2_roots[i2_index]
        i2_value = i2_component_values[int(i2_root)][0]
        row_i2 = [
            1, i2_root, i2_root**2,
            (i2_value+63)/(i2_root-27),
            i2_root*(i2_value+63)/(i2_root-27),
        ]
        for i4_index in chosen_i4_indices:
            i4_root = i4_roots[i4_index]
            chosen_values = i4_values[int(i4_root)]
            if arguments.i4_gate in ("bzero", "azero"):
                chosen_values = (None,)
            if arguments.i4_value_index is not None:
                if arguments.i4_value_index >= len(chosen_values):
                    continue
                chosen_values = (chosen_values[arguments.i4_value_index],)
            for i4_value in chosen_values:
                if arguments.i4_gate == "derivative":
                    derivatives = i4_jets[int(i4_root)][int(i4_value)]
                else:
                    derivatives = (None,)
                for i4_derivative in derivatives:
                    if arguments.i4_gate == "bzero":
                        row_i4 = [0, 0, 0, 1, i4_root]
                    elif arguments.i4_gate == "azero":
                        row_i4 = [1, i4_root, i4_root**2, 0, 0]
                    elif arguments.i4_gate == "derivative":
                        row_i4 = [
                            0, 1, 2*i4_root,
                            (
                                i4_derivative/(i4_root-27)
                                - (i4_value+63)/(i4_root-27)**2
                            ),
                            (
                                (i4_value+63)/(i4_root-27)
                                + i4_root*(
                                    i4_derivative/(i4_root-27)
                                    - (i4_value+63)/(i4_root-27)**2
                                )
                            ),
                        ]
                    else:
                        row_i4 = [
                            1, i4_root, i4_root**2,
                            (i4_value+63)/(i4_root-27),
                            i4_root*(i4_value+63)/(i4_root-27),
                        ]
                    ambient_kernel = matrix(
                        finite, [row_i2, row_i4]
                    ).right_kernel_matrix()
                    assert ambient_kernel.nrows() == 3
                    projective_rows = (
                        [(finite(1), left, right) for left in finite for right in finite]
                        + [(finite(0), finite(1), right) for right in finite]
                        + [(finite(0), finite(0), finite(1))]
                    )
                    assert len(projective_rows) == 73**2+73+1
                    for functional in projective_rows:
                        subkernel = matrix(
                            finite, [functional]
                        ).right_kernel_matrix()
                        kernel = subkernel*ambient_kernel
                        assert kernel.nrows() == 2
                        (aa, ba), (ab, bb) = map(kernel_function, kernel.rows())
                        if aa*bb-ab*ba == 0:
                            old_base_collapses += 1
                            continue
                        infinity_tests += 1
                        factorization = cover_factorization_for_kernel(kernel)
                        if factorization is None:
                            continue
                        original_row = matrix(finite, [functional])*ambient_kernel
                        record = {
                        "I2_index": i2_index,
                        "I2_root": int(i2_root),
                        "I2_value": int(i2_value),
                        "I4_index": i4_index,
                        "I4_root": int(i4_root),
                        "I4_value": (
                            None if i4_value is None else int(i4_value)
                        ),
                        "I4_gate": arguments.i4_gate,
                        "I4_derivative": (
                            None if i4_derivative is None else int(i4_derivative)
                        ),
                        "quotient_functional": list(map(int, functional)),
                        "infinity_row": list(map(int, original_row.row(0))),
                        "kernel_rows": [
                            list(map(int, row)) for row in kernel.rows()
                        ],
                        "factor_degrees_exponents": [
                            [int(factor.degree()), int(exponent)]
                            for factor, exponent in factorization
                        ],
                        }
                        infinity_hits.append(record)
                        print(
                        "Q80FINALQ6MODULEINF|"
                        f"I2={(i2_index, int(i2_root), int(i2_value))}|"
                        f"I4={(i4_index, int(i4_root), i4_value, arguments.i4_gate, i4_derivative)}|"
                        f"functional={tuple(map(int,functional))}|"
                        f"row={tuple(record['infinity_row'])}|"
                        f"kernel={tuple(tuple(map(int,row)) for row in kernel.rows())}|"
                        f"factors={tuple(tuple(item) for item in record['factor_degrees_exponents'])}|"
                        "status=GENUS_ONE_HIT",
                            flush=True,
                        )
    print(
        "Q80FINALQ6MODULEINF|"
        f"i2_index={arguments.i2_index}|i4_index={arguments.i4_index}|"
        f"i4_value_index={arguments.i4_value_index}|i4_gate={arguments.i4_gate}|"
        f"tests={infinity_tests}|"
        f"old_base_collapses={old_base_collapses}|"
        f"genus_one_hits={len(infinity_hits)}|"
        "status=PASS_EXACT_PROJECTIVE_INFINITY_SEARCH",
        flush=True,
    )
    raise SystemExit(0)


if arguments.arbitrary_a3_row:
    chosen_i2_indices = (
        range(3) if arguments.i2_index is None else (arguments.i2_index,)
    )
    arbitrary_hits = []
    arbitrary_tests = 0
    for i2_index in chosen_i2_indices:
        i2_root = i2_roots[i2_index]
        i2_value = i2_component_values[int(i2_root)][0]
        row_i2 = [
            1, i2_root, i2_root**2,
            (i2_value+63)/(i2_root-27),
            i2_root*(i2_value+63)/(i2_root-27),
        ]
        row_infinity = [0, 0, 1, 0, 0]
        ambient_kernel = matrix(finite, [row_i2, row_infinity]).right_kernel_matrix()
        assert ambient_kernel.nrows() == 3
        projective_rows = (
            [(finite(1), left, right) for left in finite for right in finite]
            + [(finite(0), finite(1), right) for right in finite]
            + [(finite(0), finite(0), finite(1))]
        )
        assert len(projective_rows) == 73**2+73+1
        for functional in projective_rows:
            subkernel = matrix(finite, [functional]).right_kernel_matrix()
            kernel = subkernel*ambient_kernel
            assert kernel.nrows() == 2
            arbitrary_tests += 1
            factorization = cover_factorization_for_kernel(kernel)
            if factorization is None:
                continue
            row = {
                "I2_index": i2_index,
                "I2_root": int(i2_root),
                "I2_value": int(i2_value),
                "quotient_functional": list(map(int, functional)),
                "kernel_rows": [list(map(int, row)) for row in kernel.rows()],
                "factor_degrees_exponents": [
                    [int(factor.degree()), int(exponent)]
                    for factor, exponent in factorization
                ],
            }
            arbitrary_hits.append(row)
            print(
                "Q80FINALQ6MODULEARB|"
                f"I2={(i2_index, int(i2_root), int(i2_value))}|"
                f"functional={tuple(map(int,functional))}|"
                f"kernel={tuple(tuple(map(int,row)) for row in kernel.rows())}|"
                f"factors={tuple(tuple(item) for item in row['factor_degrees_exponents'])}|"
                "status=GENUS_ONE_HIT",
                flush=True,
            )
    print(
        "Q80FINALQ6MODULEARB|"
        f"i2_index={arguments.i2_index}|tests={arbitrary_tests}|"
        f"genus_one_hits={len(arbitrary_hits)}|"
        "status=PASS_EXACT_PROJECTIVE_ROW_SEARCH",
        flush=True,
    )
    raise SystemExit(0)


hits = []
tests = 0
generic_factorizations = 0
for i2_root in i2_roots:
    for i2_value in i2_component_values[int(i2_root)]:
        row_i2 = [
            1, i2_root, i2_root**2,
            (i2_value+63)/(i2_root-27),
            i2_root*(i2_value+63)/(i2_root-27),
        ]
        for i4_root in i4_roots:
            i4_gates = [
                ("value", value, None)
                for value in i4_values[int(i4_root)]
            ]
            i4_gates += [
                ("derivative", value, derivative)
                for value in i4_values[int(i4_root)]
                for derivative in i4_jets[int(i4_root)][int(value)]
            ]
            i4_gates += [
                ("mixed", value, (derivative, int(scalar)))
                for value in i4_values[int(i4_root)]
                for derivative in i4_jets[int(i4_root)][int(value)]
                for scalar in finite
            ]
            i4_gates += [("b(c)=0", None, None), ("a(c)=0", None, None)]
            for i4_gate, i4_value, i4_derivative in i4_gates:
                if i4_gate == "b(c)=0":
                    row_i4 = [0, 0, 0, 1, i4_root]
                elif i4_gate == "a(c)=0":
                    row_i4 = [1, i4_root, i4_root**2, 0, 0]
                elif i4_gate == "derivative":
                    row_i4 = [
                        0,
                        1,
                        2*i4_root,
                        (
                            i4_derivative/(i4_root-27)
                            - (i4_value+63)/(i4_root-27)**2
                        ),
                        (
                            (i4_value+63)/(i4_root-27)
                            + i4_root*(
                                i4_derivative/(i4_root-27)
                                - (i4_value+63)/(i4_root-27)**2
                            )
                        ),
                    ]
                elif i4_gate == "mixed":
                    derivative, scalar = i4_derivative
                    value_row = [
                        1, i4_root, i4_root**2,
                        (i4_value+63)/(i4_root-27),
                        i4_root*(i4_value+63)/(i4_root-27),
                    ]
                    derivative_row = [
                        0, 1, 2*i4_root,
                        derivative/(i4_root-27)-(i4_value+63)/(i4_root-27)**2,
                        (
                            (i4_value+63)/(i4_root-27)
                            + i4_root*(
                                derivative/(i4_root-27)
                                - (i4_value+63)/(i4_root-27)**2
                            )
                        ),
                    ]
                    row_i4 = [
                        left+finite(scalar)*right
                        for left, right in zip(derivative_row, value_row)
                    ]
                else:
                    row_i4 = [
                        1, i4_root, i4_root**2,
                        (i4_value+63)/(i4_root-27),
                        i4_root*(i4_value+63)/(i4_root-27),
                    ]
                for infinity_value in ("a2=0",):
                    row_infinity = [0, 0, 1, 0, 0]
                    constraints = matrix(finite, [row_i2, row_i4, row_infinity])
                    if constraints.rank() != 3:
                        continue
                    kernel = constraints.right_kernel_matrix()
                    assert kernel.nrows() == 2
                    (a0, b0), (a1, b1) = map(kernel_function, kernel.rows())
                    denominator = s_field(R)*b0-b1
                    if denominator == 0:
                        continue
                    q0 = (a1-s_field(R)*a0)/denominator
                    old_chord = q0+h
                    cover = branch_for_raw_chord(old_chord)
                    square_class = cover.numerator()*cover.denominator()
                    tests += 1
                    cheap_degree = cheap_radical_degree(square_class)
                    if cheap_degree is not None and cheap_degree > 8:
                        continue
                    factorization = square_class.factor()
                    generic_factorizations += 1
                    odd_part = s_ring.one()
                    for factor, exponent in factorization:
                        if int(exponent) % 2:
                            odd_part *= factor
                    odd_degree = int(odd_part.degree())
                    if odd_degree != 4:
                        continue
                    classification = classify_quartic(
                        odd_part.monic(), parameter(factorization.unit())
                    )
                    record = {
                        "I2_root": int(i2_root),
                        "I2_value": int(i2_value),
                        "I4_root": int(i4_root),
                        "I4_gate": i4_gate,
                        "I4_value": None if i4_value is None else int(i4_value),
                        "I4_derivative": (
                            None if i4_derivative is None
                            else list(i4_derivative) if isinstance(i4_derivative, tuple)
                            else int(i4_derivative)
                        ),
                        "I0star_gate": infinity_value,
                        "constraint_rows": [list(map(int, row)) for row in constraints.rows()],
                        "kernel_rows": [list(map(int, row)) for row in kernel.rows()],
                        "q0_in_s_R": str(q0),
                        "factor_degrees_exponents": [
                            [int(factor.degree()), int(exponent)]
                            for factor, exponent in factorization
                        ],
                        "classification": classification,
                    }
                    hits.append(record)
                    print(
                        "Q80FINALQ6MODULE|"
                        f"I2={(int(i2_root), int(i2_value))}|"
                        f"I4={(int(i4_root), i4_gate, i4_value, i4_derivative)}|"
                        f"I0star={infinity_value}|kernel={tuple(tuple(map(int,row)) for row in kernel.rows())}|"
                        f"factors={tuple(tuple(row) for row in record['factor_degrees_exponents'])}|"
                        f"root_data={tuple(classification['root_data'])}|"
                        "status=GENUS_ONE_HIT",
                        flush=True,
                    )

print(
    "Q80FINALQ6MODULE|"
        f"section={arguments.section_index}|pole_center={int(pole_center)}|pole_residue={pole_residue}|"
    f"polynomial_part_s2_s={polynomial_part}|ambient_branch_degree={ambient_branch_degree}|"
    f"I2_branch_values={i2_values}|I2_component_values={i2_component_values}|"
    f"I4_values={i4_values}|I4_factorizations={i4_factorizations}|I4_jets={i4_jets}|"
    f"infinity_polynomial={infinity_polynomial}|infinity_gate=a2=0|"
    f"tests={tests}|generic_factorizations={generic_factorizations}|"
    f"genus_one_hits={len(hits)}|"
    "status=PASS_EXACT_LOCAL_MODULE_SEARCH",
    flush=True,
)

pinned_hits = [
    row for row in hits
    if (
        row["I2_root"], row["I2_value"], row["I4_root"],
        row["I4_gate"], row["I4_value"], row["I0star_gate"]
    ) == (72, 6, 64, "value", 3, "a2=0")
]
assert len(pinned_hits) == 1
assert pinned_hits[0]["kernel_rows"] == [[1, 0, 0, 41, 48], [0, 1, 0, 6, 72]]
assert pinned_hits[0]["classification"]["root_data"] == [15, 66, 800]

if arguments.write_artifact:
    artifact = {
        "schema": "q80-final-q6-saturated-module-gf73-v1",
        "status": "exact_finite_field_certificate",
        "prime": 73,
        "section_index": arguments.section_index,
        "source_artifacts": [
            {"path": str(path.relative_to(ROOT)), "sha256": digest}
            for path, digest in KNOWN_HASHES.items()
        ],
        "smith_certificate": {
            "H": str(H),
            "invariant_factors": ["1", str(H.monic())],
            "numerator_module_rows": [
                [str(entry) for entry in row]
                for row in smith_numerator_module.rows()
            ],
            "saturated_generator": "qsat=(q0+63)/(s-27)",
            "saturated_numerator_column": [
                str(entry) for entry in saturated_numerator_column
            ],
        },
        "ambient_basis": ["1", "s", "s^2", "qsat", "s*qsat"],
        "pinned_constraints": {
            "I2": {"root": 72, "q0_value": 6},
            "I4": {"root": 64, "q0_value": 3},
            "I0star_row": [0, 0, 1, 0, 0],
        },
        "selected_hit": pinned_hits[0],
        "target_root_data": [15, 66, 800],
        "rank_claim": None,
        "reproduce": (
            "sage elkies-k3/scripts/search_q80_final_q6_local_module_gf73.sage "
            "--section-index=0 --write-artifact"
        ),
    }
    output_path = ROOT / (
        "artifacts/generated-results/q80-final-q6-saturated-module-gf73.json"
    )
    encoded = json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n"
    output_path.write_text(encoded)
    print(
        "Q80FINALQ6MODULE|"
        f"artifact={output_path}|"
        f"sha256={hashlib.sha256(encoded.encode()).hexdigest()}|"
        "target=(15,66,800)|status=PASS_ARTIFACT_WRITE",
        flush=True,
    )

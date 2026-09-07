#!/usr/bin/env sage
"""Certify the two q=80 rank-19 tangent branches over GF(7).

The seed is the unique normalized q=80 one-A1 surface over GF(7) carrying
the two polynomial section directions that continue generically.  The third
direction is the one-pole class selected by exact function-field pole gates in
``verify_q80_cm24_seed_gf7.sage``.

There are 46 variables: four surface parameters, two polynomial sections,
and one section with an unknown simple pole.  Their Weierstrass identities
give 13+13+19=45 equations.  After resolving the forced component-zero
coordinates, the exact Jacobian has nullity two rather than one.  The exact
quadratic obstruction splits into two lines, and both normalized directions
lift uniquely through the displayed finite formal order.  This is a bounded
formal calculation, not a proof that either branch algebraizes over QQ.
"""

from sage.all import GF, Matrix, PolynomialRing, PowerSeriesRing, binomial, vector


PROTOCOL = "Q80RANK19GF7"
field = GF(7)

names = (
    ("d", "p", "q", "e")
    + tuple(f"x1{index}" for index in range(5))
    + tuple(f"y1{index}" for index in range(7))
    + tuple(f"x2{index}" for index in range(5))
    + tuple(f"y2{index}" for index in range(7))
    + ("lam",)
    + tuple(f"n{index}" for index in range(7))
    + tuple(f"m{index}" for index in range(10))
)
parameters = PolynomialRing(field, names=names)
variables = parameters.gens()
parameter_field = parameters.fraction_field()
variables_k = tuple(map(parameter_field, variables))

(
    d,
    p,
    q,
    e,
    *section_variables,
) = variables_k
x1 = section_variables[0:5]
y1 = section_variables[5:12]
x2 = section_variables[12:17]
y2 = section_variables[17:24]
lam = section_variables[24]
n = section_variables[25:32]
m = section_variables[32:42]
assert len(variables) == 46 and len(m) == 10

polynomials = PolynomialRing(parameter_field, "T")
T = polynomials.gen()
series = PolynomialRing(parameter_field, "s")
s = series.gen()

r = -3 * d**2 + 3 - p - q
A = T**2 * (-3 + p*T + q*T**2 + r*T**3)

A_at_one = series(A(T=1+s))
u = (A_at_one + 3*d**2) / (-3*d**2)
branch = 2*d**3 * (
    1 + parameter_field(3)/2*u + parameter_field(3)/8*u**2
    - parameter_field(1)/16*u**3
)
branch_jets = vector(parameter_field, [branch[index] for index in range(4)])
jet_matrix = Matrix(
    parameter_field,
    4,
    4,
    lambda row, column: series((1+s)**(4+column))[row],
)
fixed_jets = vector(
    parameter_field,
    [series(2*(1+s)**3 + e*(1+s)**8)[index] for index in range(4)],
)
b1, b2, b3, b4 = jet_matrix.solve_right(branch_jets-fixed_jets)
B = T**3 * (2 + b1*T + b2*T**2 + b3*T**3 + b4*T**4 + e*T**5)


def coordinate_polynomial(coefficients):
    return sum(coefficients[index]*T**index for index in range(len(coefficients)))


X1 = coordinate_polynomial(x1)
Y1 = coordinate_polynomial(y1)
X2 = coordinate_polynomial(x2)
Y2 = coordinate_polynomial(y2)
N = coordinate_polynomial(n)
M = coordinate_polynomial(m)
h = T-lam

residuals = (
    Y1**2-X1**3-A*X1-B,
    Y2**2-X2**3-A*X2-B,
    M**2-N**3-A*N*h**4-B*h**6,
)
equations = tuple(
    residual[index]
    for residual, expected_degree in zip(residuals, (12, 12, 18))
    for index in range(expected_degree+1)
)
assert len(equations) == 45

# Both continuing nonidentity directions pass through the I4 node at T=1.
# These equations select their strict-transform branches; on the singular
# Weierstrass model their first-order content is not recovered from the raw
# section identity alone.
component_equations = (
    X1(T=1)-d,
    Y1(T=1),
    N(T=1)-d*(1-lam)**2,
    M(T=1),
    # At T=0 the D5 first exceptional cubic is (X-1)^2*(X+2).
    # P1 follows its double-root branch.
    x1[1]-1,
    # First exceptional conics at the I4 node.  Squaring retains both split
    # orientations while removing the singular tangent of the raw node chart.
    (6*d*Y1.derivative()(T=1))**2
    - 3*d*(6*d*X1.derivative()(T=1)+A.derivative()(T=1))**2,
    (
        6*d*(M/h**3).derivative()(T=1)
    )**2
    - 3*d*(
        6*d*(N/h**2).derivative()(T=1)+A.derivative()(T=1)
    )**2,
)
resolved_equations = equations + component_equations

seed_values = (
    (3, 4, 3, 2)
    + (0, 1, 2, 0, 0)
    + (0, 0, 2, 1, 4, 0, 0)
    + (1, 2, 2, 6, 4)
    + (1, 3, 3, 2, 3, 3, 6)
    + (6,)
    + (1, 4, 3, 0, 6, 1, 4)
    # The selected third direction is the negative of scanner candidate 3.
    + (6, 1, 5, 2, 4, 3, 6, 5, 4, 6)
)
assert len(seed_values) == len(variables)
seed = {variable: field(value) for variable, value in zip(variables, seed_values)}

evaluations = tuple(field(equation.subs(seed)) for equation in equations)
assert not any(evaluations)
assert not any(field(equation.subs(seed)) for equation in component_equations)

jacobian = Matrix(
    field,
    [
        [field(equation.derivative(variable).subs(seed)) for variable in variables]
        for equation in resolved_equations
    ],
)
raw_rank = jacobian.rank()
print(
    f"{PROTOCOL}|stage=raw_jacobian|rank={raw_rank}"
    f"|nullity={jacobian.right_kernel().dimension()}",
    flush=True,
)

# The nonidentity D5 and E6 component chart for P1 is obtained by the local
# divisibilities ord_0(X1,Y1)>=(1,2) and ord_infinity(X1bar,Y1bar)>=(2,2).
# These are precisely the seven zero coordinates below.  Leaving them in a
# dense Weierstrass section ansatz produces the familiar singular section
# scheme at the rational double points.
component_zero_names = {
    "x10", "y10", "y11", "x13", "x14", "y15", "y16",
}
active_columns = [
    index for index, name in enumerate(names) if name not in component_zero_names
]
active_names = tuple(names[index] for index in active_columns)
resolved_jacobian = jacobian.matrix_from_columns(active_columns)
rank = resolved_jacobian.rank()
tangent = resolved_jacobian.right_kernel_matrix()
print(
    f"{PROTOCOL}|stage=resolved_jacobian|rank={rank}|nullity={tangent.nrows()}",
    flush=True,
)
assert len(active_names) == 39
assert rank == 37 and tangent.nrows() == 2
for diagnostic_index, diagnostic_row in enumerate(tangent.rows(), 1):
    support = tuple(
        active_names[index]
        for index, value in enumerate(diagnostic_row)
        if value
    )
    print(
        f"{PROTOCOL}|stage=tangent|row={diagnostic_index}"
        f"|support={','.join(support)}|values="
        + ",".join(map(str, map(int, diagnostic_row))),
        flush=True,
    )

# Compute the quadratic Kuranishi obstruction in the two tangent variables.
# Work in the exact truncated algebra F_7[tau0,tau1]/(tau0,tau1)^3 so the
# rational coefficient functions are expanded without choosing integer lifts.
tangent_cone_ring = PolynomialRing(field, names=("tau0", "tau1"))
tau0, tau1 = tangent_cone_ring.gens()
tangent_cone_quotient = tangent_cone_ring.quotient(
    (tau0**3, tau0**2*tau1, tau0*tau1**2, tau1**3),
    names=("T0", "T1"),
)
T0, T1 = tangent_cone_quotient.gens()
active_lookup = {column: index for index, column in enumerate(active_columns)}
cone_images = []
for column, variable in enumerate(variables):
    value = tangent_cone_quotient(field(seed[variable]))
    if column in active_lookup:
        active_index = active_lookup[column]
        value += tangent_cone_quotient(tangent[0, active_index])*T0
        value += tangent_cone_quotient(tangent[1, active_index])*T1
    cone_images.append(value)
cone_map = parameters.hom(cone_images, tangent_cone_quotient)
cone_values = vector(
    tangent_cone_quotient,
    [
        cone_map(parameters(equation.numerator()))
        * cone_map(parameters(equation.denominator())).inverse_of_unit()
        for equation in resolved_equations
    ],
)
left_kernel = resolved_jacobian.left_kernel_matrix()
cone_obstructions = left_kernel.change_ring(tangent_cone_quotient)*cone_values
cone_polynomials = []
for obstruction in cone_obstructions:
    polynomial = tangent_cone_ring(obstruction.lift())
    if polynomial and polynomial not in cone_polynomials:
        cone_polynomials.append(polynomial)
cone_ideal = tangent_cone_ring.ideal(cone_polynomials)
cone_basis = tuple(cone_ideal.groebner_basis())
expected_cone = tau0**2-tau0*tau1+tau1**2
assert cone_basis == (expected_cone,)
assert expected_cone == (tau0+2*tau1)*(tau0-3*tau1)
print(
    f"{PROTOCOL}|stage=quadratic_cone|equation={expected_cone}"
    f"|factorization={expected_cone.factor()}",
    flush=True,
)

# Lift both tangent lines in F_7[[h]], fixing p=4+h to remove formal
# reparametrization.  At each order the next-order compatibility condition
# selects exactly one of the seven remaining kernel corrections.
equation_numerators = tuple(
    parameters(equation.numerator()) for equation in resolved_equations
)
equation_denominators = tuple(
    parameters(equation.denominator()) for equation in resolved_equations
)
seed_active = vector(
    field, [field(seed[variables[column]]) for column in active_columns]
)
p_column = active_names.index("p")


formal_series_ring = PowerSeriesRing(field, "h", default_prec=512)
finite_jet_matrix_inverse = Matrix(
    field, 4, 4, lambda row, column: binomial(4+column, row)
).inverse()


def _poly_add(left, right):
    length = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else formal_series_ring.zero())
        +(right[index] if index < len(right) else formal_series_ring.zero())
        for index in range(length)
    ]


def _poly_mul(left, right, truncate=None):
    length = len(left)+len(right)-1
    if truncate is not None:
        length = min(length, truncate)
    answer = [formal_series_ring.zero() for _ in range(length)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i+j < length:
                answer[i+j] += left_value*right_value
    return answer


def _poly_pow(value, exponent, truncate=None):
    answer = [formal_series_ring.one()]
    for _ in range(exponent):
        answer = _poly_mul(answer, value, truncate=truncate)
    return answer


def _poly_subtract(first, *others):
    answer = first
    for other in others:
        answer = _poly_add(answer, [-value for value in other])
    return answer


def _evaluate_formal_series(values):
    D, P, Q, E = values[:4]
    remaining = values[4:]
    X1s = list(remaining[0:5])
    Y1s = list(remaining[5:12])
    X2s = list(remaining[12:17])
    Y2s = list(remaining[17:24])
    LAMs = remaining[24]
    Ns = list(remaining[25:32])
    Ms = list(remaining[32:42])
    R = -3*D**2+3-P-Q
    As = [formal_series_ring.zero(), formal_series_ring.zero(), formal_series_ring(-3), P, Q, R]
    A_at_one = [
        sum(
            (formal_series_ring(binomial(index, derivative))*As[index]
             for index in range(derivative, 6)),
            formal_series_ring.zero(),
        )
        for derivative in range(4)
    ]
    denominator = -3*D**2
    u_values = list(A_at_one)
    u_values[0] += 3*D**2
    u_values = [value/denominator for value in u_values]
    branch_values = [formal_series_ring.one()]+[formal_series_ring.zero()]*3
    for scalar, exponent in ((field(3)/2, 1), (field(3)/8, 2), (-field(1)/16, 3)):
        branch_values = _poly_add(
            branch_values,
            [scalar*value for value in _poly_pow(u_values, exponent, truncate=4)],
        )
    branch_values = [2*D**3*value for value in branch_values]
    fixed_values = [
        formal_series_ring(2*binomial(3, derivative))+E*binomial(8, derivative)
        for derivative in range(4)
    ]
    difference = [branch_values[index]-fixed_values[index] for index in range(4)]
    middle_B = [
        sum(
            (formal_series_ring(finite_jet_matrix_inverse[row, column])*difference[column]
             for column in range(4)),
            formal_series_ring.zero(),
        )
        for row in range(4)
    ]
    Bs = [formal_series_ring.zero()]*3+[formal_series_ring(2)]+middle_B+[E]
    Hs = [-LAMs, formal_series_ring.one()]
    residual_polynomials = (
        _poly_subtract(_poly_pow(Y1s, 2), _poly_pow(X1s, 3), _poly_mul(As, X1s), Bs),
        _poly_subtract(_poly_pow(Y2s, 2), _poly_pow(X2s, 3), _poly_mul(As, X2s), Bs),
        _poly_subtract(
            _poly_pow(Ms, 2), _poly_pow(Ns, 3),
            _poly_mul(_poly_mul(As, Ns), _poly_pow(Hs, 4)),
            _poly_mul(Bs, _poly_pow(Hs, 6)),
        ),
    )
    result = []
    for residual, expected_degree in zip(residual_polynomials, (12, 12, 18)):
        residual += [formal_series_ring.zero()]*(expected_degree+1-len(residual))
        result.extend(residual[:expected_degree+1])
    value_at_one = lambda values: sum(values, formal_series_ring.zero())
    derivative_at_one = lambda values: sum(
        (formal_series_ring(index)*values[index] for index in range(1, len(values))),
        formal_series_ring.zero(),
    )
    X1_one = value_at_one(X1s)
    Y1_one = value_at_one(Y1s)
    N_one = value_at_one(Ns)
    M_one = value_at_one(Ms)
    A_prime = derivative_at_one(As)
    X1_prime = derivative_at_one(X1s)
    Y1_prime = derivative_at_one(Y1s)
    H_one = 1-LAMs
    M_affine_prime = derivative_at_one(Ms)/H_one**3-3*M_one/H_one**4
    N_affine_prime = derivative_at_one(Ns)/H_one**2-2*N_one/H_one**3
    result.extend((
        X1_one-D, Y1_one, N_one-D*H_one**2, M_one, X1s[1]-1,
        (6*D*Y1_prime)**2-3*D*(6*D*X1_prime+A_prime)**2,
        (6*D*M_affine_prime)**2-3*D*(6*D*N_affine_prime+A_prime)**2,
    ))
    assert len(result) == len(resolved_equations) == 52
    return result


def residual_coefficient(coefficients, order):
    images = []
    for column, variable in enumerate(variables):
        if column in active_lookup:
            active_index = active_lookup[column]
            value = formal_series_ring(
                [coefficient[active_index] for coefficient in coefficients]
            ).add_bigoh(order+1)
        else:
            value = formal_series_ring(field(seed[variable])).add_bigoh(order+1)
        images.append(value)
    return vector(field, [value[order] for value in _evaluate_formal_series(images)])


formal_order = 7
surface_names = ("d", "p", "q", "e")
surface_columns = tuple(active_names.index(name) for name in surface_names)
branch_records = []
for slope in (5, 3):
    first_order = slope*tangent[0]+tangent[1]
    assert first_order[p_column] == 1
    coefficients = [seed_active, vector(field, first_order)]
    for order in range(2, formal_order):
        right_hand_side = -residual_coefficient(coefficients, order)
        assert not any(left_kernel*right_hand_side)
        particular = resolved_jacobian.solve_right(right_hand_side)
        survivors = []
        for coefficient0 in field:
            for coefficient1 in field:
                correction = (
                    particular
                    + coefficient0*tangent[0]
                    + coefficient1*tangent[1]
                )
                if correction[p_column]:
                    continue
                trial = coefficients+[vector(field, correction)]
                next_right_hand_side = -residual_coefficient(trial, order+1)
                if not any(left_kernel*next_right_hand_side):
                    survivors.append(trial)
        assert len(survivors) == 1
        coefficients = survivors[0]
    surface_series = tuple(
        tuple(int(coefficient[column]) for coefficient in coefficients)
        for column in surface_columns
    )
    branch_records.append((slope, surface_series))
    print(
        f"{PROTOCOL}|stage=formal_branch|line=tau0:{slope},tau1:1"
        f"|surface_tangent="
        + ",".join(str(series[1]) for series in surface_series)
        + "|series="
        + ";".join(
            f"{name}:{','.join(map(str, series))}"
            for name, series in zip(surface_names, surface_series)
        ),
        flush=True,
    )

assert tuple(record[1][0][1] for record in branch_records) == (5, 3)
print(
    f"{PROTOCOL}|raw_variables={len(variables)}|resolved_variables={len(active_names)}"
    f"|raw_equations={len(equations)}|resolved_equations={len(resolved_equations)}"
    f"|jacobian_rank={rank}|tangent_dimension=2|formal_order={formal_order}"
    f"|status=PASS",
    flush=True,
)

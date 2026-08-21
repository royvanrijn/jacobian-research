#!/usr/bin/env sage
"""Verify the historical standard-P1 q=80 seed and its obstructions mod 11.

The P2 component gate was later found to be wrong at the residual I2, so this
script is a reproducible negative-branch diagnostic.  The variables encode the I1*+I4+IV* ambient, the extra I2 node, and the
transported generic MW basis with pole profile (0,0,1).  All equations are
polynomial.  Full column rank of the overdetermined Jacobian certifies that a
square subsystem can be selected for multivariate Hensel lifting.
"""

from itertools import product

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, rational_reconstruction, vector


field = GF(11)
names = (
    "d", "p", "q", "r", "b1", "b2", "b3", "b4", "e", "rho", "w",
    "a", "b", "u0", "u1", "u2",
    "k2", "l2", "x21", "x22", "x23",
    "y21", "y22", "y23", "y24", "y25",
    "z30", "z31", "z32", "h3", "l3",
    "x31", "x32", "x33", "x34", "x35", "x36", "x37", "x38", "x39",
    "y31", "y32", "y33", "y34", "y35", "y36", "y37", "y38",
    "y39", "y310", "y311", "y312", "y313", "y314",
)
parameters = PolynomialRing(ZZ, names=names)
(
    d, p, q, r, b1, b2, b3, b4, e, rho, w,
    a, b, u0, u1, u2,
    k2, l2, x21, x22, x23,
    y21, y22, y23, y24, y25,
    z30, z31, z32, h3, l3,
    x31, x32, x33, x34, x35, x36, x37, x38, x39,
    y31, y32, y33, y34, y35, y36, y37, y38,
    y39, y310, y311, y312, y313, y314,
) = parameters.gens()
polynomials = PolynomialRing(parameters, "T")
T = polynomials.gen()

A = T**2 * (-3 + p*T + q*T**2 + r*T**3)
B = T**3 * (2 + b1*T + b2*T**2 + b3*T**3 + b4*T**4 + e*T**5)
discriminant = 4*A**3 + 27*B**2

equations = []
equation_tags = []


def add(tag, values):
    for index, value in enumerate(values):
        value = parameters(value)
        if value:
            equations.append(value)
            equation_tags.append(f"{tag}:{index}")


# The normalized I4 branch at T=1.  The first two equations choose its sign;
# the first three nonconstant discriminant jets then impose exact order >=4.
add("ambient", (A(1) + 3*d**2, B(1) - 2*d**3))
add("I4", (discriminant.derivative(j)(1) for j in range(1, 4)))

# The new I2 at rho: w is the double root of the singular cubic and the
# discriminant derivative supplies the second vanishing order.
add(
    "I2",
    (
        w**3 + A(rho)*w + B(rho),
        3*w**2 + A(rho),
        discriminant.derivative()(rho),
    ),
)

# P1: nonidentity at D5, A3, and E6.  The visible incidences force x(0)=0,
# x(1)=d, deg(x)<=3 and y divisible by T^2(T-1), deg(y)<=5.
x1 = T * (a + b*T + (d-a-b)*T**2)
y1 = T**2 * (T-1) * (u0 + u1*T + u2*T**2)
add("P1", (y1**2 - x1**3 - A*x1 - B).list())
add("P1_I2_incidence", (x1(rho) - w, y1(rho)))

# P2: identity components at zero and infinity are parametrized by k2,l2.
x2 = k2**2 + x21*T + x22*T**2 + x23*T**3 + l2**2*T**4
y2 = k2**3 + y21*T + y22*T**2 + y23*T**3 + y24*T**4 + y25*T**5 + l2**3*T**6
add("P2", (y2**2 - x2**3 - A*x2 - B).list())

# P3 has P.O=3 generically.  At CM24 a quadratic factor cancels from Z,
# X, and Y with multiplicities 1,2,3, leaving the visible P.O=1 section.
Z3 = T**3 + z32*T**2 + z31*T + z30
X3 = (
    h3**2*z30**2
    + x31*T + x32*T**2 + x33*T**3 + x34*T**4 + x35*T**5
    + x36*T**6 + x37*T**7 + x38*T**8 + x39*T**9
    + l3**2*T**10
)
Y3 = (
    h3**3*z30**3
    + y31*T + y32*T**2 + y33*T**3 + y34*T**4 + y35*T**5
    + y36*T**6 + y37*T**7 + y38*T**8 + y39*T**9
    + y310*T**10 + y311*T**11 + y312*T**12 + y313*T**13
    + y314*T**14 + l3**3*T**15
)
add("P3", (Y3**2 - X3**3 - A*X3*Z3**4 - B*Z3**6).list())
add(
    "P3_incidence",
    (
        X3(1) - d*Z3(1)**2,
        Y3(1),
        X3(rho) - w*Z3(rho)**2,
        Y3(rho),
    ),
)

# Component-correct survivor found by the compiled full mod-11 scan.
seed_values = (
    3, 2, 3, 4, 1, 0, 10, 10, 9, 9, 6,
    1, 2, 6, 3, 0,
    8, 1, 6, 4, 7,
    6, 0, 3, 1, 5,
    6, 1, 6, 7, 5,
    7, 2, 2, 0, 5, 2, 10, 3, 0,
    1, 6, 2, 5, 3, 3, 1, 5, 5, 9, 8, 5, 6, 0,
)
assert len(seed_values) == len(names)
seed = {generator: ZZ(value) for generator, value in zip(parameters.gens(), seed_values)}
residuals = vector(ZZ, [equation.subs(seed) for equation in equations])
failures = [
    (tag, value % 11)
    for tag, value in zip(equation_tags, residuals)
    if value % 11
]
assert not failures, failures

jacobian = matrix(
    field,
    [
        [field(equation.derivative(generator).subs(seed)) for generator in parameters.gens()]
        for equation in equations
    ],
)
rank = jacobian.rank()
pivot_rows = tuple(jacobian.transpose().pivots())
assert len(pivot_rows) == rank

# Remove the CM-only I2 equations/incidences and rho,w to inspect the actual
# one-dimensional generic q=80 rank-19 deformation system at the same seed.
generic_row_indices = [
    index
    for index, tag in enumerate(equation_tags)
    if not tag.startswith("I2:")
    and not tag.startswith("P1_I2_incidence:")
    and tag not in ("P3_incidence:2", "P3_incidence:3")
]
generic_column_indices = [index for index in range(len(names)) if index not in (9, 10)]
generic_jacobian = jacobian.matrix_from_rows_and_columns(
    generic_row_indices, generic_column_indices
)
generic_kernel = generic_jacobian.right_kernel_matrix()
generic_surface_projection = generic_kernel.matrix_from_columns(range(9))
print(
    f"Q80GENERICSEED|variables={len(generic_column_indices)}|"
    f"equations={len(generic_row_indices)}|rank={generic_jacobian.rank()}|"
    f"kernel={generic_kernel.nrows()}|"
    f"surface_tangent_dimension={generic_surface_projection.rank()}",
    flush=True,
)

generic_left_kernel = generic_jacobian.left_kernel_matrix()
generic_residuals = vector(ZZ, [residuals[index] for index in generic_row_indices])
generic_first_rhs = vector(
    field, [field(-(value // 11)) for value in generic_residuals]
)
generic_first_liftable = not any(generic_left_kernel * generic_first_rhs)
print(
    f"Q80GENERICLIFT|modulus=121|liftable={ZZ(generic_first_liftable)}",
    flush=True,
)
if generic_first_liftable:
    generic_first_delta = generic_jacobian.solve_right(generic_first_rhs)
    independent_surface_rows = generic_surface_projection.transpose().pivots()
    assert len(independent_surface_rows) == 3
    pure_coefficients = generic_surface_projection.left_kernel_matrix()
    assert pure_coefficients.nrows() == 8
    adapted_generic_kernel = matrix(
        field,
        [generic_kernel[index] for index in independent_surface_rows]
    ).stack(pure_coefficients * generic_kernel)
    assert adapted_generic_kernel.rank() == 11

    def expand_generic_row(row):
        expanded = vector(field, [0] * len(names))
        for index, column in enumerate(generic_column_indices):
            expanded[column] = row[index]
        return expanded

    generic_kernel_full = matrix(
        field, [expand_generic_row(row) for row in adapted_generic_kernel.rows()]
    )
    generic_delta_full = expand_generic_row(generic_first_delta)
    generic_tangent_names = tuple(f"g{index + 1}" for index in range(11))
    generic_tangent_ZZ = PolynomialRing(ZZ, names=generic_tangent_names)
    generic_tangent_F = PolynomialRing(field, names=generic_tangent_names)
    generic_tangent_ZZ_variables = generic_tangent_ZZ.gens()
    generic_first_lift = vector(ZZ, seed_values) + 11 * vector(
        ZZ, [ZZ(value) for value in generic_delta_full]
    )
    generic_tangent_substitution = {}
    for column, generator in enumerate(parameters.gens()):
        value = generic_tangent_ZZ(generic_first_lift[column])
        for row, tangent in enumerate(generic_kernel_full.rows()):
            value += 11 * ZZ(tangent[column]) * generic_tangent_ZZ_variables[row]
        generic_tangent_substitution[generator] = value
    generic_second_values = []
    for equation_index in generic_row_indices:
        expanded = generic_tangent_ZZ(
            equations[equation_index].subs(generic_tangent_substitution)
        )
        assert all(
            coefficient % 121 == 0 for coefficient in expanded.coefficients()
        )
        divided = generic_tangent_ZZ(
            {
                monomial: coefficient // 121
                for monomial, coefficient in expanded.dict().items()
            }
        )
        generic_second_values.append(generic_tangent_F(divided))
    generic_obstructions = []
    for left_row in generic_left_kernel.rows():
        obstruction = -sum(
            generic_tangent_F(left_row[index]) * generic_second_values[index]
            for index in range(len(generic_row_indices))
        )
        if obstruction:
            generic_obstructions.append(obstruction)
    generic_obstruction_ideal = generic_tangent_F.ideal(generic_obstructions)
    generic_obstruction_groebner = generic_obstruction_ideal.groebner_basis()
    print(
        f"Q80GENERICOBSTRUCTION|variables=11|equations={len(generic_obstructions)}|"
        f"dimension={generic_obstruction_ideal.dimension()}|"
        f"surface_coordinates={generic_tangent_names[:3]}|"
        f"pure_coordinates={generic_tangent_names[3:]}",
        flush=True,
    )
    print(
        "Q80GENERICOBSTRUCTION|leading_monomials="
        + ";".join(str(polynomial.lm()) for polynomial in generic_obstruction_groebner),
        flush=True,
    )
    generic_tangent_variables = generic_tangent_F.gens()
    generic_surface_elimination = generic_obstruction_ideal.elimination_ideal(
        generic_tangent_variables[3:]
    )
    generic_surface_ring = PolynomialRing(
        field, names=generic_tangent_names[:3]
    )
    generic_surface_equations = [
        generic_surface_ring(generator)
        for generator in generic_surface_elimination.gens()
        if generator
    ]
    generic_surface_ideal = generic_surface_ring.ideal(generic_surface_equations)
    print(
        f"Q80GENERICSURFACE|equations={len(generic_surface_equations)}|"
        f"dimension={generic_surface_ideal.dimension()}|"
        "ideal=" + ";".join(map(str, generic_surface_equations)),
        flush=True,
    )
    generic_slice_points = []
    generic_slice_point = None
    generic_slice_trials = 0
    slice_ring = PolynomialRing(field, names=generic_tangent_names[:5])
    slice_variables = slice_ring.gens()
    for trial in range(128):
        generic_slice_trials += 1
        fixed_values = []
        state = ZZ(trial + 1)
        for _ in range(6):
            state = (1103515245 * state + 12345) % (2**31)
            fixed_values.append(field(state % 11))
        fixed_substitution = dict(
            zip(generic_tangent_variables[5:], fixed_values)
        )
        slice_polynomials = [
            slice_ring(polynomial.subs(fixed_substitution))
            for polynomial in generic_obstructions
        ]
        slice_ideal = slice_ring.ideal(slice_polynomials)
        if slice_ideal.dimension() != 0:
            continue
        for point in slice_ideal.variety(ring=field):
            coordinates = [
                point.get(variable, field.zero()) for variable in slice_variables
            ] + fixed_values
            if any(coordinates[:3]):
                coordinate_tuple = tuple(coordinates)
                if coordinate_tuple not in generic_slice_points:
                    generic_slice_points.append(coordinate_tuple)
                if len(generic_slice_points) >= 64:
                    break
        if len(generic_slice_points) >= 64:
            break
    if generic_slice_points:
        generic_slice_point = generic_slice_points[0]
    print(
        f"Q80GENERICTANGENT|slice_trials={generic_slice_trials}|"
        f"points={len(generic_slice_points)}|first="
        + (
            ",".join(map(str, generic_slice_point))
            if generic_slice_point is not None else ""
        ),
        flush=True,
    )
    if generic_slice_point is not None:
        selected_generic_delta = vector(field, generic_first_delta)
        for coefficient, tangent in zip(
            generic_slice_point, adapted_generic_kernel.rows()
        ):
            selected_generic_delta += coefficient * tangent
        reduced_seed = vector(
            ZZ, [seed_values[index] for index in generic_column_indices]
        )
        generic_current = reduced_seed + 11 * vector(
            ZZ, [ZZ(value) for value in selected_generic_delta]
        )

        def generic_unfixed_values(point):
            full_point = vector(ZZ, seed_values)
            for index, column in enumerate(generic_column_indices):
                full_point[column] = point[index]
            substitution = dict(zip(parameters.gens(), full_point))
            return vector(ZZ, [
                equations[index].subs(substitution)
                for index in generic_row_indices
            ])

        def advance_unfixed_generic(current, exponent):
            modulus = 11**exponent
            values = generic_unfixed_values(current)
            rhs = vector(field, [field(-(value // modulus)) for value in values])
            if any(generic_left_kernel * rhs):
                return None
            particular = generic_jacobian.solve_right(rhs)
            next_modulus = modulus * 11

            def obstruction(coefficients):
                delta = vector(field, particular)
                for coefficient, tangent in zip(
                    coefficients, adapted_generic_kernel.rows()
                ):
                    delta += coefficient * tangent
                candidate = current + modulus * vector(
                    ZZ, [ZZ(value) for value in delta]
                )
                candidate_values = generic_unfixed_values(candidate)
                assert all(value % next_modulus == 0 for value in candidate_values)
                next_rhs = vector(field, [
                    field(-(value // next_modulus)) for value in candidate_values
                ])
                return generic_left_kernel * next_rhs, candidate

            zero = vector(field, [0] * adapted_generic_kernel.nrows())
            constant, _ = obstruction(zero)
            columns = []
            for index in range(adapted_generic_kernel.nrows()):
                unit = vector(field, [0] * adapted_generic_kernel.nrows())
                unit[index] = 1
                value, _ = obstruction(unit)
                columns.append(value - constant)
            compatibility = matrix(field, columns).transpose()
            try:
                coefficients = compatibility.solve_right(-constant)
            except ValueError:
                return None
            final, candidate = obstruction(coefficients)
            assert not any(final)
            return vector(ZZ, [ZZ(value % next_modulus) for value in candidate])

        selected_surviving_point = None
        best_unfixed_exponent = 2
        for candidate_point in generic_slice_points:
            candidate_delta = vector(field, generic_first_delta)
            for coefficient, tangent in zip(
                candidate_point, adapted_generic_kernel.rows()
            ):
                candidate_delta += coefficient * tangent
            candidate_current = reduced_seed + 11 * vector(
                ZZ, [ZZ(value) for value in candidate_delta]
            )
            unfixed_trial = candidate_current
            unfixed_exponent = 2
            while unfixed_trial is not None and unfixed_exponent < 7:
                unfixed_trial = advance_unfixed_generic(
                    unfixed_trial, unfixed_exponent
                )
                unfixed_exponent += 1
            best_unfixed_exponent = max(
                best_unfixed_exponent,
                unfixed_exponent if unfixed_trial is not None else unfixed_exponent - 1,
            )
            if unfixed_trial is not None:
                generic_slice_point = candidate_point
                selected_generic_delta = candidate_delta
                generic_current = candidate_current
                selected_surviving_point = candidate_point
                break
        print(
            f"Q80GENERICUNFIXED|survives_to_exponent="
            f"{best_unfixed_exponent}|selected="
            + (
                ",".join(map(str, selected_surviving_point))
                if selected_surviving_point is not None else ""
            ),
            flush=True,
        )
        if selected_surviving_point is None:
            changing_surface_columns = [
                column
                for column in range(9)
                if selected_generic_delta[column] != 0
            ]
        else:
            changing_surface_columns = [
                column
                for column in range(9)
                if selected_generic_delta[column] != 0
            ]
        assert changing_surface_columns
        print(
            "Q80GENERICTANGENT|surface_delta="
            + ",".join(map(str, selected_generic_delta[:9])),
            flush=True,
        )
        parameter_column = changing_surface_columns[2]
        parameter_target = ZZ(generic_current[parameter_column])
        selector = vector(field, [0] * len(generic_column_indices))
        selector[parameter_column] = 1
        fixed_jacobian = generic_jacobian.stack(matrix(field, [selector]))
        fixed_left_kernel = fixed_jacobian.left_kernel_matrix()
        fixed_kernel = fixed_jacobian.right_kernel_matrix()

        def generic_exact_values(point):
            full_point = vector(ZZ, seed_values)
            for index, column in enumerate(generic_column_indices):
                full_point[column] = point[index]
            substitution = dict(zip(parameters.gens(), full_point))
            values = [
                equations[index].subs(substitution)
                for index in generic_row_indices
            ]
            values.append(point[parameter_column] - parameter_target)
            return vector(ZZ, values)

        def advance_fixed_generic(current, exponent):
            modulus = 11**exponent
            values = generic_exact_values(current)
            assert all(value % modulus == 0 for value in values)
            rhs = vector(field, [field(-(value // modulus)) for value in values])
            if any(fixed_left_kernel * rhs):
                return None
            particular = fixed_jacobian.solve_right(rhs)
            next_modulus = modulus * 11

            def obstruction(coefficients):
                delta = vector(field, particular)
                for coefficient, tangent in zip(coefficients, fixed_kernel.rows()):
                    delta += coefficient * tangent
                candidate = current + modulus * vector(
                    ZZ, [ZZ(value) for value in delta]
                )
                candidate_values = generic_exact_values(candidate)
                assert all(value % next_modulus == 0 for value in candidate_values)
                next_rhs = vector(
                    field,
                    [field(-(value // next_modulus)) for value in candidate_values],
                )
                return fixed_left_kernel * next_rhs, candidate

            zero = vector(field, [0] * fixed_kernel.nrows())
            constant, _ = obstruction(zero)
            columns = []
            for index in range(fixed_kernel.nrows()):
                unit = vector(field, [0] * fixed_kernel.nrows())
                unit[index] = 1
                value, _ = obstruction(unit)
                columns.append(value - constant)
            compatibility = matrix(field, columns).transpose()
            try:
                coefficients = compatibility.solve_right(-constant)
            except ValueError:
                return None
            final, candidate = obstruction(coefficients)
            assert not any(final)
            return vector(
                ZZ, [ZZ(value % next_modulus) for value in candidate]
            )

        generic_exponent = 2
        while generic_exponent < 20:
            advanced = advance_fixed_generic(generic_current, generic_exponent)
            if advanced is None:
                break
            generic_current = advanced
            generic_exponent += 1
        generic_modulus = 11**generic_exponent
        print(
            f"Q80GENERICPADIC|parameter={names[generic_column_indices[parameter_column]]}|"
            f"target={parameter_target}|exponent={generic_exponent}|"
            f"surface_residues={','.join(map(str, generic_current[:9]))}",
            flush=True,
        )
        generic_reconstructed = []
        generic_reconstruction_failed = False
        for residue in generic_current:
            try:
                generic_reconstructed.append(
                    QQ(rational_reconstruction(residue, generic_modulus))
                )
            except Exception:
                generic_reconstruction_failed = True
                break
        generic_exact_pass = False
        if not generic_reconstruction_failed:
            generic_exact_pass = all(
                value == 0 for value in generic_exact_values(generic_reconstructed)
            )
        print(
            f"Q80GENERICRECOGNIZE|all_coordinates="
            f"{ZZ(not generic_reconstruction_failed)}|"
            f"exact_pass={ZZ(generic_exact_pass)}|"
            "surface="
            + (
                ",".join(map(str, generic_reconstructed[:9]))
                if not generic_reconstruction_failed else ""
            ),
            flush=True,
        )

print(
    f"Q80CM24SEED|prime=11|variables={len(names)}|equations={len(equations)}|"
    f"jacobian_rank={rank}",
    flush=True,
)
print(
    "Q80CM24SEED|pivot_equations="
    + ",".join(equation_tags[index] for index in pivot_rows),
    flush=True,
)
kernel = jacobian.right_kernel_matrix()
for index, tangent in enumerate(kernel.rows(), 1):
    support = [names[column] for column, value in enumerate(tangent) if value]
    surface_support = [name for name in support if name in names[:11]]
    print(
        f"Q80CM24SEED|tangent={index}|surface_support={','.join(surface_support)}|"
        f"full_support={','.join(support)}",
        flush=True,
    )
if rank == len(names):
    print("Q80CM24SEED|hensel_regular=1|status=PASS", flush=True)
else:
    print("Q80CM24SEED|hensel_regular=0|status=SINGULAR_SEED", flush=True)

# First singular-Hensel step.  The seed lifts to 11^2 exactly when
# J*delta=-F(seed)/11 is consistent over GF(11).
first_rhs = vector(field, [field(-(value // 11)) for value in residuals])
left_kernel = jacobian.left_kernel_matrix()
first_obstruction = left_kernel * first_rhs
first_liftable = not any(first_obstruction)
print(
    f"Q80CM24LIFT|modulus=121|liftable={ZZ(first_liftable)}|"
    f"left_obstruction={tuple(first_obstruction)}",
    flush=True,
)

if first_liftable:
    first_delta = jacobian.solve_right(first_rhs)
    assert jacobian * first_delta == first_rhs
    surface_kernel = [
        tangent
        for tangent in kernel.rows()
        if any(tangent[column] for column in range(11))
    ]
    assert len(surface_kernel) == 3
    seed_vector = vector(ZZ, seed_values)
    second_order_survivors = []
    for coefficients in product(range(11), repeat=3):
        delta = vector(field, first_delta)
        for coefficient, tangent in zip(coefficients, surface_kernel):
            delta += field(coefficient) * tangent
        lifted = seed_vector + 11 * vector(ZZ, [ZZ(value) for value in delta])
        lifted_substitution = dict(zip(parameters.gens(), lifted))
        values = vector(ZZ, [
            equation.subs(lifted_substitution) for equation in equations
        ])
        assert all(value % 121 == 0 for value in values)
        second_rhs = vector(field, [field(-(value // 121)) for value in values])
        if not any(left_kernel * second_rhs):
            second_order_survivors.append((coefficients, lifted, second_rhs))
    print(
        f"Q80CM24LIFT|modulus=1331|surface_tuples_tested={11**3}|"
        f"survivors={len(second_order_survivors)}|"
        "coefficients="
        + ";".join(
            ",".join(map(str, coefficients))
            for coefficients, _, _ in second_order_survivors[:32]
        ),
        flush=True,
    )

    # Compute the full second-order obstruction ideal in all nine tangent
    # coordinates.  Substitution by x_1+11*K*c is coefficientwise divisible
    # by 11^2; after division and reduction, the Jacobian left-kernel gives
    # exact quadratic compatibility equations over GF(11).
    tangent_names = tuple(f"c{index + 1}" for index in range(kernel.nrows()))
    tangent_integer_ring = PolynomialRing(ZZ, names=tangent_names)
    tangent_integer_variables = tangent_integer_ring.gens()
    tangent_field_ring = PolynomialRing(field, names=tangent_names)
    tangent_field_variables = tangent_field_ring.gens()
    first_lift = seed_vector + 11 * vector(
        ZZ, [ZZ(value) for value in first_delta]
    )
    tangent_substitution = {}
    for column, generator in enumerate(parameters.gens()):
        value = tangent_integer_ring(first_lift[column])
        for row, tangent in enumerate(kernel.rows()):
            value += 11 * ZZ(tangent[column]) * tangent_integer_variables[row]
        tangent_substitution[generator] = value

    second_order_values = []
    for equation in equations:
        expanded = tangent_integer_ring(equation.subs(tangent_substitution))
        assert all(coefficient % 121 == 0 for coefficient in expanded.coefficients())
        divided = tangent_integer_ring(
            {
                monomial: coefficient // 121
                for monomial, coefficient in expanded.dict().items()
            }
        )
        second_order_values.append(tangent_field_ring(divided))
    obstruction_polynomials = []
    for left_row in left_kernel.rows():
        obstruction = -sum(
            tangent_field_ring(left_row[index]) * second_order_values[index]
            for index in range(len(equations))
        )
        if obstruction:
            obstruction_polynomials.append(obstruction)
    obstruction_ideal = tangent_field_ring.ideal(obstruction_polynomials)
    print(
        f"Q80CM24OBSTRUCTION|variables={len(tangent_names)}|"
        f"equations={len(obstruction_polynomials)}|"
        f"dimension={obstruction_ideal.dimension()}",
        flush=True,
    )
    obstruction_groebner = obstruction_ideal.groebner_basis()
    print(
        "Q80CM24OBSTRUCTION|basis_size="
        + str(len(obstruction_groebner))
        + "|leading_monomials="
        + ";".join(str(polynomial.lm()) for polynomial in obstruction_groebner),
        flush=True,
    )
    def exact_values(point):
        substitution = dict(zip(parameters.gens(), point))
        return vector(ZZ, [equation.subs(substitution) for equation in equations])

    def advance_with_lookahead(current, exponent):
        """Lift one digit while choosing kernel coordinates for the next."""
        modulus = 11**exponent
        values = exact_values(current)
        assert all(value % modulus == 0 for value in values)
        rhs = vector(field, [field(-(value // modulus)) for value in values])
        if any(left_kernel * rhs):
            return None
        particular = jacobian.solve_right(rhs)
        next_modulus = modulus * 11

        def lookahead_obstruction(kernel_coefficients):
            delta = vector(field, particular)
            for coefficient, tangent in zip(kernel_coefficients, kernel.rows()):
                delta += coefficient * tangent
            candidate = current + modulus * vector(
                ZZ, [ZZ(value) for value in delta]
            )
            candidate_values = exact_values(candidate)
            assert all(value % next_modulus == 0 for value in candidate_values)
            next_rhs = vector(
                field,
                [field(-(value // next_modulus)) for value in candidate_values],
            )
            return left_kernel * next_rhs, candidate

        zero_coefficients = vector(field, [0] * kernel.nrows())
        constant_obstruction, _ = lookahead_obstruction(zero_coefficients)
        columns = []
        for index in range(kernel.nrows()):
            unit = vector(field, [0] * kernel.nrows())
            unit[index] = 1
            obstruction, _ = lookahead_obstruction(unit)
            columns.append(obstruction - constant_obstruction)
        lookahead_matrix = matrix(field, columns).transpose()
        try:
            kernel_coefficients = lookahead_matrix.solve_right(-constant_obstruction)
        except ValueError:
            return None
        final_obstruction, candidate = lookahead_obstruction(kernel_coefficients)
        assert not any(final_obstruction)
        return vector(
            ZZ, [ZZ(value % next_modulus) for value in candidate]
        )

    sparse_obstruction_points = []
    selected_branch = None
    for pure_index in range(3, 9):
        for values in product(range(11), repeat=4):
            point = [field.zero()] * 9
            point[:3] = map(field, values[:3])
            point[pure_index] = field(values[3])
            substitution = dict(zip(tangent_field_variables, point))
            if all(
                polynomial.subs(substitution) == 0
                for polynomial in obstruction_polynomials
            ):
                sparse_obstruction_points.append(tuple(point))
                selected_first_delta = vector(field, first_delta)
                for coefficient, tangent in zip(point, kernel.rows()):
                    selected_first_delta += coefficient * tangent
                candidate_mod_121 = seed_vector + 11 * vector(
                    ZZ, [ZZ(value) for value in selected_first_delta]
                )
                trial = candidate_mod_121
                trial_exponent = 2
                while trial is not None and trial_exponent < 7:
                    trial = advance_with_lookahead(trial, trial_exponent)
                    trial_exponent += 1
                if trial is not None:
                    selected_branch = (tuple(point), trial, trial_exponent)
                    break
        if selected_branch is not None:
            break
    print(
        f"Q80CM24OBSTRUCTION|one_pure_coordinate_search={6 * 11**4}|"
        f"points_tested={len(sparse_obstruction_points)}|"
        "selected="
        + (
            ",".join(map(str, selected_branch[0]))
            if selected_branch is not None else ""
        ),
        flush=True,
    )

    # Once the quadratic gate has been crossed, one-digit look-ahead is
    # linear in the nine kernel parameters: Hessian terms acquire an extra
    # factor of 11.  Continue the selected branch to high p-adic precision.
    if selected_branch is not None:
        _, current, exponent = selected_branch
        while exponent < 20:
            advanced = advance_with_lookahead(current, exponent)
            assert advanced is not None
            current = advanced
            exponent += 1
        final_modulus = 11**exponent
        print(
            f"Q80CM24PADIC|prime=11|exponent={exponent}|modulus={final_modulus}|"
            "surface_residues="
            + ",".join(map(str, current[:11])),
            flush=True,
        )

        reconstructed = []
        reconstruction_failed = False
        for residue in current:
            try:
                reconstructed.append(QQ(rational_reconstruction(residue, final_modulus)))
            except Exception:
                reconstruction_failed = True
                break
        exact_pass = False
        if not reconstruction_failed:
            exact_substitution = dict(zip(parameters.gens(), reconstructed))
            exact_pass = all(
                equation.subs(exact_substitution) == 0 for equation in equations
            )
        print(
            f"Q80CM24RECOGNIZE|all_coordinates={ZZ(not reconstruction_failed)}|"
            f"exact_pass={ZZ(exact_pass)}|"
            "surface="
            + (
                ",".join(map(str, reconstructed[:11]))
                if not reconstruction_failed else ""
            ),
            flush=True,
        )

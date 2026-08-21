#!/usr/bin/env sage
"""Analyze a bounded split-prime implicit ideal for a q=80 CM24 branch.

The input is produced by ``extend_q80_rank19_branches_gf7.sage``.  Its
generators are relations on a finite formal jet, so every conclusion here is
about the resulting candidate ideal; it is not by itself a characteristic-zero
algebraization certificate.
"""

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, PowerSeriesRing, ZZ, vector


parser = argparse.ArgumentParser()
parser.add_argument(
    "--input",
    default="artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json",
)
parser.add_argument(
    "--check-prime",
    action="store_true",
    help="run the expensive full projective primality check",
)
parser.add_argument(
    "--check-genus",
    action="store_true",
    help="run the expensive black-box projective genus check",
)
parser.add_argument(
    "--plane-pair",
    nargs=2,
    metavar=("X", "Y"),
    help="eliminate to the indicated affine coordinate pair and test its plane model",
)
parser.add_argument("--print-plane", action="store_true")
parser.add_argument(
    "--plane-delta",
    action="store_true",
    help="compute local delta invariants on the selected plane projection",
)
parser.add_argument(
    "--plane-normalize",
    action="store_true",
    help="run positive-characteristic normalization on the selected affine plane model",
)
parser.add_argument(
    "--plane-projective-normalize",
    action="store_true",
    help="normalize the homogeneous coordinate ring of the selected plane model",
)
parser.add_argument(
    "--print-normalization-generators",
    action="store_true",
    help="include large affine/projective normalization module generators in output",
)
parser.add_argument(
    "--projective-parameter",
    action="store_true",
    help="recover a degree-one parameter from the normalized rational normal curve",
)
parser.add_argument(
    "--parameter-output",
    help="write the recovered modular surface parameterization to JSON",
)
parser.add_argument(
    "--normalization-ratios",
    action="store_true",
    help="search ratios of linear normalization generators for a degree-one parameter",
)
parser.add_argument(
    "--ratio-degree",
    type=int,
    default=1,
    help="maximum monomial degree used by --normalization-ratios",
)
parser.add_argument(
    "--search-normalization-units",
    action="store_true",
    help="exhaustively search normalized linear forms for a nonconstant unit",
)
parser.add_argument(
    "--split-normalization-unit",
    action="store_true",
    help="split the quadratic map supplied by a discovered degree-two unit",
)
parser.add_argument(
    "--coordinate-fibers",
    action="store_true",
    help="report affine fiber lengths over every GF(7) value of each coordinate",
)
parser.add_argument(
    "--exclude-parameter",
    help="verify that an exact GF(7) surface parameter does not lie in this ideal",
)
arguments = parser.parse_args()

input_path = Path(arguments.input)
raw_bytes = input_path.read_bytes()
artifact = json.loads(raw_bytes)
schema = artifact.get("schema")
if schema == "q80-cm24-formal-branch-ideal-v1":
    field_text = artifact.get("field", "")
    if not (field_text.startswith("GF(") and field_text.endswith(")")):
        raise ValueError("the branch ideal does not declare a prime field")
    prime = ZZ(field_text[3:-1])
    variable_names = tuple(artifact["affine_variables"])
    relation_rows = None
elif (
    schema == "q80-cm24-split-prime-formal-branch-v1"
    and artifact.get("kind") == "canonical_centered_relation_space"
):
    prime = ZZ(artifact["prime"])
    variable_names = tuple(artifact["centered_variables"])
    relation_rows = artifact["rref_basis"]
else:
    raise ValueError("unexpected q80 branch ideal schema")
if not prime.is_prime():
    raise ValueError("the branch ideal characteristic is not prime")
finite = GF(prime)


def stage(name, started, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(
        f"Q80RANK19IDEAL|stage={name}|seconds={time.monotonic()-started:.3f}"
        + (f"|{payload}" if payload else ""),
        flush=True,
    )


started = time.monotonic()
plane_geometric_genus = None
affine_ring = PolynomialRing(
    finite, names=variable_names, order="degrevlex"
)
if relation_rows is None:
    affine_generators = tuple(
        affine_ring(polynomial) for polynomial in artifact["affine_generators"]
    )
else:
    monomials = tuple(
        affine_ring.monomial(*map(ZZ, exponents))
        for exponents in artifact["monomial_exponents"]
    )
    affine_generators = tuple(
        sum(
            (finite(coefficient)*monomial
             for coefficient, monomial in zip(row, monomials)),
            affine_ring.zero(),
        )
        for row in relation_rows
    )
affine_ideal = affine_ring.ideal(affine_generators)
stage(
    "affine",
    started,
    generators=len(affine_generators),
    dimension=affine_ideal.dimension(),
)

if arguments.exclude_parameter:
    started = time.monotonic()
    excluded = json.loads(Path(arguments.exclude_parameter).read_text())
    parameter_ring = PolynomialRing(finite, "t")
    parameter_field = parameter_ring.fraction_field()
    images = tuple(
        parameter_field(excluded["functions"][name]["value"])
        for name in artifact["affine_variables"]
    )
    residuals = tuple(
        generator(*images) for generator in affine_generators
    )
    stage(
        "exclude_parameter",
        started,
        input=arguments.exclude_parameter,
        generators=len(residuals),
        vanishing=sum(not residual for residual in residuals),
        excluded=int(all(residual for residual in residuals)),
    )

if arguments.coordinate_fibers:
    started = time.monotonic()
    for variable in affine_ring.gens():
        lengths = []
        for value in finite:
            fiber_ideal = affine_ideal+affine_ring.ideal(variable-value)
            lengths.append(
                -1 if fiber_ideal.dimension() != 0 else fiber_ideal.vector_space_dimension()
            )
        stage(
            "coordinate_fibers",
            started,
            coordinate=variable,
            lengths=tuple(lengths),
            maximum=max(lengths),
        )

if arguments.plane_pair:
    keep_names = tuple(arguments.plane_pair)
    if len(set(keep_names)) != 2 or not set(keep_names) <= set(affine_ring.variable_names()):
        raise ValueError("--plane-pair must name two distinct affine variables")
    eliminated_variables = tuple(
        variable
        for variable in affine_ring.gens()
        if str(variable) not in keep_names
    )
    started = time.monotonic()
    elimination_ideal = affine_ideal.elimination_ideal(eliminated_variables)
    plane_ring = PolynomialRing(finite, names=keep_names, order="degrevlex")
    plane_generators = tuple(
        plane_ring(str(polynomial)) for polynomial in elimination_ideal.gens()
    )
    plane_ideal = plane_ring.ideal(plane_generators)
    factor_degrees = ()
    irreducible = None
    if len(plane_generators) == 1:
        factorization = plane_generators[0].factor()
        factor_degrees = tuple(
            (factor.total_degree(), multiplicity)
            for factor, multiplicity in factorization
        )
        irreducible = len(factorization) == 1 and factorization[0][1] == 1
    stage(
        "plane_projection",
        started,
        coordinates=",".join(keep_names),
        generators=len(plane_generators),
        dimension=plane_ideal.dimension(),
        degrees=",".join(str(polynomial.total_degree()) for polynomial in plane_generators),
        bidegrees=",".join(
            f"{polynomial.degree(plane_ring.gen(0))}:{polynomial.degree(plane_ring.gen(1))}"
            for polynomial in plane_generators
        ),
        irreducible="NA" if irreducible is None else int(irreducible),
        factor_degrees=factor_degrees,
    )
    if arguments.print_plane:
        print(
            "Q80RANK19IDEAL|stage=plane_equations|"
            f"coordinates={','.join(keep_names)}|equations={plane_generators}",
            flush=True,
        )
    if arguments.plane_delta:
        if len(plane_generators) != 1 or not irreducible:
            raise ValueError("local delta analysis requires one irreducible plane equation")
        from sage.libs.singular.function_factory import ff

        delta_loc = ff.normal__lib.deltaLoc

        def chart_delta(polynomial, extra_equations, label):
            chart_ring = polynomial.parent()
            singular_ideal = chart_ring.ideal(
                [polynomial]
                + [polynomial.derivative(variable) for variable in chart_ring.gens()]
                + list(extra_equations)
            )
            if singular_ideal == chart_ring.ideal(1):
                stage(label, time.monotonic(), components=0, delta=0)
                return 0
            component_started = time.monotonic()
            components = singular_ideal.minimal_associated_primes()
            values = tuple(tuple(map(ZZ, delta_loc(polynomial, component))) for component in components)
            stage(
                label,
                component_started,
                components=len(components),
                singular_scheme_degree=singular_ideal.vector_space_dimension(),
                invariants=values,
                delta=sum(value[0] for value in values),
            )
            return sum(value[0] for value in values)

        delta_started = time.monotonic()
        affine_delta = chart_delta(plane_generators[0], (), "plane_delta_affine")

        projective_ring = PolynomialRing(
            finite, names=("X", "Y", "Z"), order="degrevlex"
        )
        X, Y, Z = projective_ring.gens()
        plane_degree = plane_generators[0].total_degree()
        projective_polynomial = projective_ring.zero()
        for exponents, coefficient in plane_generators[0].dict().items():
            projective_polynomial += (
                coefficient
                * X**exponents[0]
                * Y**exponents[1]
                * Z**(plane_degree-sum(exponents))
            )

        infinity_x_ring = PolynomialRing(
            finite, names=("y", "z"), order="degrevlex"
        )
        y, z = infinity_x_ring.gens()
        infinity_x_polynomial = infinity_x_ring(
            projective_polynomial.subs({X: 1, Y: y, Z: z})
        )
        infinity_delta = chart_delta(
            infinity_x_polynomial, (z,), "plane_delta_infinity_X"
        )

        infinity_y_ring = PolynomialRing(
            finite, names=("x", "z"), order="degrevlex"
        )
        x, z = infinity_y_ring.gens()
        infinity_y_polynomial = infinity_y_ring(
            projective_polynomial.subs({X: x, Y: 1, Z: z})
        )
        origin = infinity_y_ring.ideal(x, z)
        if all(polynomial in origin for polynomial in (
            infinity_y_polynomial,
            infinity_y_polynomial.derivative(x),
            infinity_y_polynomial.derivative(z),
        )):
            infinity_delta += chart_delta(
                infinity_y_polynomial, (x, z), "plane_delta_infinity_Y_extra"
            )
        total_delta = affine_delta+infinity_delta
        plane_arithmetic_genus = (plane_degree-1)*(plane_degree-2)//2
        plane_geometric_genus = plane_arithmetic_genus-total_delta
        stage(
            "plane_delta_total",
            delta_started,
            arithmetic_genus=plane_arithmetic_genus,
            delta=total_delta,
            geometric_genus=plane_geometric_genus,
        )
    if arguments.plane_normalize:
        if len(plane_generators) != 1 or not irreducible:
            raise ValueError("plane normalization requires one irreducible equation")
        from sage.libs.singular.function_factory import ff

        started = time.monotonic()
        normalization = ff.normal__lib.normalP(plane_ideal)
        module_generators = tuple(normalization[0][0])
        affine_normalization_data = dict(
            components=len(normalization[0]),
            delta_data=normalization[1],
        )
        if arguments.print_normalization_generators:
            affine_normalization_data["module_generators"] = normalization[0]
        stage(
            "plane_affine_normalization",
            started,
            **affine_normalization_data,
        )
        denominator = module_generators[-1]
        numerators = module_generators[:-1]
        normalization_ring = PolynomialRing(
            finite,
            names=keep_names+tuple(f"U{index}" for index in range(len(numerators))),
            order="degrevlex",
        )
        normalization_variables = normalization_ring.gens()

        def embed_plane(polynomial):
            result = normalization_ring.zero()
            for exponents, coefficient in polynomial.dict().items():
                monomial = normalization_ring.one()
                for variable, exponent in zip(normalization_variables[:2], exponents):
                    monomial *= variable**exponent
                result += coefficient*monomial
            return result

        embedded_denominator = embed_plane(denominator)
        presentation = normalization_ring.ideal(
            [embed_plane(plane_generators[0])]
            + [
                embedded_denominator*normalization_variables[2+index]
                - embed_plane(numerator)
                for index, numerator in enumerate(numerators)
            ]
        )
        presentation, presentation_saturation = presentation.saturation(
            normalization_ring.ideal(embedded_denominator)
        )
        stage(
            "plane_normalization_presentation",
            time.monotonic(),
            variables=",".join(normalization_ring.variable_names()),
            dimension=presentation.dimension(),
            saturation_exponent=presentation_saturation,
            groebner_generators=len(presentation.groebner_basis()),
        )
        for variable in normalization_variables:
            lengths = []
            for value in finite:
                fiber = presentation+normalization_ring.ideal(variable-value)
                lengths.append(
                    -1 if fiber.dimension() != 0 else fiber.vector_space_dimension()
                )
            stage(
                "plane_normalization_fibers",
                time.monotonic(),
                coordinate=variable,
                lengths=tuple(lengths),
                maximum=max(lengths),
            )
        if arguments.normalization_ratios:
            ratio_basis = tuple(
                monomial
                for degree in range(arguments.ratio_degree+1)
                for monomial in normalization_ring.monomials_of_degree(degree)
            )
            ratio_names = tuple(map(str, ratio_basis))
            for numerator_index in range(len(ratio_basis)):
                for denominator_index in range(len(ratio_basis)):
                    if numerator_index == denominator_index:
                        continue
                    numerator = ratio_basis[numerator_index]
                    ratio_denominator = ratio_basis[denominator_index]
                    lengths = []
                    for value in finite:
                        fiber = presentation+normalization_ring.ideal(
                            numerator-value*ratio_denominator
                        )
                        fiber, _ = fiber.saturation(
                            normalization_ring.ideal(ratio_denominator)
                        )
                        lengths.append(
                            -1 if fiber.dimension() != 0 else fiber.vector_space_dimension()
                        )
                    if max(lengths) <= 2:
                        stage(
                            "plane_normalization_ratio",
                            time.monotonic(),
                            ratio=f"{ratio_names[numerator_index]}/{ratio_names[denominator_index]}",
                            lengths=tuple(lengths),
                            maximum=max(lengths),
                        )
        if arguments.search_normalization_units:
            unit_basis = (normalization_ring.one(),)+normalization_variables
            unit_names = ("1",)+tuple(map(str, normalization_variables))
            searched = 0
            unit_hits = 0
            parameter_hits = 0
            degree_two_units = []
            search_started = time.monotonic()
            for pivot in range(len(unit_basis)):
                for tail in itertools.product(
                    range(finite.order()), repeat=len(unit_basis)-pivot-1
                ):
                    coefficients = (0,)*pivot+(1,)+tail
                    candidate = sum(
                        finite(coefficient)*basis_element
                        for coefficient, basis_element in zip(coefficients, unit_basis)
                    )
                    searched += 1
                    if any(
                        candidate-finite(value) in presentation
                        for value in range(finite.order())
                    ):
                        continue
                    if not (presentation+normalization_ring.ideal(candidate)).is_one():
                        continue
                    unit_hits += 1
                    lengths = []
                    for value in finite:
                        fiber = presentation+normalization_ring.ideal(candidate-value)
                        lengths.append(
                            0 if fiber.is_one() else fiber.vector_space_dimension()
                        )
                    degree = max(lengths)
                    if degree == 1:
                        parameter_hits += 1
                    if degree == 2:
                        degree_two_units.append((candidate, tuple(lengths)))
                    if degree <= 2:
                        stage(
                            "plane_normalization_unit",
                            time.monotonic(),
                            coefficients=coefficients,
                            basis=unit_names,
                            function=candidate,
                            lengths=tuple(lengths),
                            degree=degree,
                        )
            stage(
                "plane_normalization_unit_search",
                search_started,
                searched=searched,
                units=unit_hits,
                degree_one=parameter_hits,
            )
            if arguments.split_normalization_unit:
                if not degree_two_units:
                    raise ValueError("no degree-two normalization unit was found")
                for unit_index, (unit, lengths) in enumerate(degree_two_units):
                    branch_values = tuple(
                        finite(index) for index, length in enumerate(lengths) if length == 1
                    )
                    if len(branch_values) != 2:
                        raise ValueError("degree-two unit does not have two rational branch values")
                    split_ring = PolynomialRing(
                        finite,
                        names=normalization_ring.variable_names()+("T",),
                        order="degrevlex",
                    )
                    split_variables = split_ring.gens()

                    def embed_normalization(polynomial):
                        result = split_ring.zero()
                        for exponents, coefficient in polynomial.dict().items():
                            monomial = split_ring.one()
                            for variable, exponent in zip(split_variables, exponents):
                                monomial *= variable**exponent
                            result += coefficient*monomial
                        return result

                    split_base = split_ring.ideal(
                        tuple(embed_normalization(polynomial) for polynomial in presentation.gens())
                    )
                    embedded_unit = embed_normalization(unit)
                    T = split_variables[-1]
                    for scalar in (finite(1), finite(3)):
                        quadratic = (
                            (embedded_unit-branch_values[1])*T**2
                            - scalar*(embedded_unit-branch_values[0])
                        )
                        split_ideal = split_base+split_ring.ideal(quadratic)
                        split_started = time.monotonic()
                        components = split_ideal.minimal_associated_primes()
                        stage(
                            "plane_normalization_unit_split",
                            split_started,
                            unit=unit_index,
                            branch_values=branch_values,
                            scalar=scalar,
                            components=len(components),
                            component_generators=tuple(tuple(component.gens()) for component in components),
                        )
    if arguments.plane_projective_normalize:
        if len(plane_generators) != 1 or not irreducible:
            raise ValueError("projective normalization requires one irreducible equation")
        from sage.libs.singular.function_factory import ff

        started = time.monotonic()
        projective_plane_ideal = plane_ideal.homogenize()
        projective_normalization = ff.normal__lib.normalP(
            projective_plane_ideal, "isPrim"
        )
        projective_module = tuple(projective_normalization[0][0])
        projective_normalization_data = dict(
            variables=projective_plane_ideal.ring().variable_names(),
            module_generator_count=len(projective_module),
            module_generator_degrees=tuple(
                polynomial.total_degree() for polynomial in projective_module
            ),
            delta_data=projective_normalization[1],
        )
        if arguments.print_normalization_generators:
            projective_normalization_data["module_generators"] = projective_module
        stage(
            "plane_projective_normalization",
            started,
            **projective_normalization_data,
        )
        projective_plane_ring = projective_plane_ideal.ring()
        projective_denominator = projective_module[-1]
        projective_numerators = projective_module[:-1]
        normalized_ring = PolynomialRing(
            finite,
            names=projective_plane_ring.variable_names()
            + tuple(f"V{index}" for index in range(len(projective_numerators))),
            order="degrevlex",
        )
        normalized_variables = normalized_ring.gens()

        def embed_projective(polynomial):
            result = normalized_ring.zero()
            for exponents, coefficient in polynomial.dict().items():
                monomial = normalized_ring.one()
                for variable, exponent in zip(normalized_variables[:3], exponents):
                    monomial *= variable**exponent
                result += coefficient*monomial
            return result

        embedded_projective_denominator = embed_projective(projective_denominator)
        normalized_presentation = normalized_ring.ideal(
            [embed_projective(projective_plane_ideal.gens()[0])]
            + [
                embedded_projective_denominator*normalized_variables[3+index]
                - embed_projective(numerator)
                for index, numerator in enumerate(projective_numerators)
            ]
        )
        normalized_presentation, normalized_saturation = normalized_presentation.saturation(
            normalized_ring.ideal(embedded_projective_denominator)
        )
        normalized_hilbert_polynomial = normalized_presentation.hilbert_polynomial(
            algorithm="singular"
        )
        stage(
            "plane_projective_normalization_presentation",
            time.monotonic(),
            variables=normalized_ring.variable_names(),
            dimension=normalized_presentation.dimension()-1,
            saturation_exponent=normalized_saturation,
            hilbert_polynomial=normalized_hilbert_polynomial,
            groebner_generators=len(normalized_presentation.groebner_basis()),
        )
        parameter_hits = []
        for numerator_index, numerator in enumerate(normalized_variables):
            for denominator_index, denominator in enumerate(normalized_variables):
                if numerator_index == denominator_index:
                    continue
                lengths = []
                for value in finite:
                    fiber = normalized_presentation+normalized_ring.ideal(
                        denominator-1, numerator-value
                    )
                    lengths.append(
                        0 if fiber.is_one() else fiber.vector_space_dimension()
                    )
                degree = max(lengths)
                if degree <= 2:
                    parameter_hits.append((numerator_index, denominator_index, degree))
                    stage(
                        "plane_projective_normalization_ratio",
                        time.monotonic(),
                        ratio=f"{numerator}/{denominator}",
                        lengths=tuple(lengths),
                        degree=degree,
                    )
        stage(
            "plane_projective_normalization_ratio_search",
            time.monotonic(),
            tested=len(normalized_variables)*(len(normalized_variables)-1),
            degree_at_most_two=len(parameter_hits),
            degree_one=sum(degree == 1 for _, _, degree in parameter_hits),
        )
        if arguments.projective_parameter:
            precision = 40
            series_ring = PowerSeriesRing(finite, "s", default_prec=precision)
            s = series_ring.gen()

            def evaluate_plane_series(polynomial, d_series, q_series):
                return sum(
                    coefficient*d_series**exponents[0]*q_series**exponents[1]
                    for exponents, coefficient in polynomial.dict().items()
                )

            plane_polynomial = plane_generators[0]
            q_series = series_ring.zero()
            plane_q_derivative = plane_polynomial.derivative(plane_ring.gen(1))
            for _ in range(8):
                residual = evaluate_plane_series(plane_polynomial, s, q_series)
                derivative = evaluate_plane_series(plane_q_derivative, s, q_series)
                q_series = (q_series-residual/derivative).add_bigoh(precision)
            assert not evaluate_plane_series(plane_polynomial, s, q_series)

            def evaluate_projective_series(polynomial):
                return sum(
                    coefficient*s**exponents[0]*q_series**exponents[1]
                    for exponents, coefficient in polynomial.dict().items()
                )

            denominator_series = evaluate_projective_series(projective_denominator)
            normalized_coordinate_series = [s, q_series, series_ring.one()]
            for numerator in projective_numerators:
                numerator_series = evaluate_projective_series(numerator)
                denominator_valuation = denominator_series.valuation()
                numerator_valuation = numerator_series.valuation()
                if numerator_valuation < denominator_valuation:
                    raise AssertionError("normalization coordinate has a pole at CM24")
                normalized_coordinate_series.append(
                    (
                        numerator_series.shift(-denominator_valuation)
                        / denominator_series.shift(-denominator_valuation)
                    ).add_bigoh(precision-denominator_valuation)
                )
            normalized_coordinate_series = tuple(normalized_coordinate_series)
            jet_matrix = Matrix(
                finite,
                [
                    [coordinate[order] for coordinate in normalized_coordinate_series]
                    for order in range(9)
                ],
            )
            if jet_matrix.rank() != 9:
                raise AssertionError("normalized degree-eight coordinates lack a full CM jet")
            order_eight_kernel = jet_matrix.matrix_from_rows(range(8)).right_kernel_matrix()
            order_seven_kernel = jet_matrix.matrix_from_rows(range(7)).right_kernel_matrix()
            if order_eight_kernel.nrows() != 1 or order_seven_kernel.nrows() != 2:
                raise AssertionError("unexpected osculating flag dimensions")
            high_vector = vector(finite, order_eight_kernel.row(0))
            low_vector = next(
                vector(finite, row)
                for row in order_seven_kernel.rows()
                if Matrix(finite, [high_vector, row]).rank() == 2
            )

            def normalized_vector(row):
                pivot = next(value for value in row if value)
                return row/pivot

            high_vector = normalized_vector(high_vector)
            low_vector = normalized_vector(low_vector)
            high_form = sum(
                coefficient*variable
                for coefficient, variable in zip(high_vector, normalized_variables)
            )
            low_form = sum(
                coefficient*variable
                for coefficient, variable in zip(low_vector, normalized_variables)
            )
            high_series = sum(
                coefficient*coordinate
                for coefficient, coordinate in zip(
                    high_vector, normalized_coordinate_series
                )
            )
            low_series = sum(
                coefficient*coordinate
                for coefficient, coordinate in zip(
                    low_vector, normalized_coordinate_series
                )
            )
            if (high_series.valuation(), low_series.valuation()) != (8, 7):
                raise AssertionError("osculating forms do not have orders eight and seven")
            parameter_series = (high_series/low_series).add_bigoh(precision-7)
            inverse_series = parameter_series.reverse().add_bigoh(precision-8)
            d_in_parameter = inverse_series
            q_in_parameter = q_series(inverse_series).add_bigoh(precision-8)

            def pade_candidate(series, numerator_degree, denominator_degree):
                coefficient_count = series.prec()
                coefficients = tuple(series[index] for index in range(coefficient_count))
                if denominator_degree == 0:
                    if any(coefficients[numerator_degree+1:]):
                        return None
                    return coefficients[:numerator_degree+1], (finite.one(),)
                rows = tuple(range(numerator_degree+1, coefficient_count))
                coefficient_matrix = Matrix(
                    finite,
                    [
                        [
                            coefficients[index-offset] if index >= offset else finite.zero()
                            for offset in range(1, denominator_degree+1)
                        ]
                        for index in rows
                    ],
                )
                target = vector(
                    finite, [-coefficients[index] for index in rows]
                )
                if coefficient_matrix.rank() != denominator_degree:
                    return None
                if target not in coefficient_matrix.column_space():
                    return None
                denominator = (finite.one(),)+tuple(
                    coefficient_matrix.solve_right(target)
                )
                numerator = tuple(
                    sum(
                        denominator[offset]*coefficients[index-offset]
                        for offset in range(min(denominator_degree, index)+1)
                    )
                    for index in range(numerator_degree+1)
                )
                return numerator, denominator

            parameter_polynomial_ring = PolynomialRing(finite, "t")
            t = parameter_polynomial_ring.gen()
            parameter_field = parameter_polynomial_ring.fraction_field()

            def minimal_pade(series, maximum_degree):
                for total_degree in range(2*maximum_degree+1):
                    candidates = []
                    for denominator_degree in range(maximum_degree+1):
                        numerator_degree = total_degree-denominator_degree
                        if 0 <= numerator_degree <= maximum_degree:
                            candidate = pade_candidate(
                                series, numerator_degree, denominator_degree
                            )
                            if candidate is not None:
                                candidates.append((numerator_degree, denominator_degree, candidate))
                    if len(candidates) == 1:
                        numerator_degree, denominator_degree, (numerator, denominator) = candidates[0]
                        numerator_polynomial = sum(
                            coefficient*t**index for index, coefficient in enumerate(numerator)
                        )
                        denominator_polynomial = sum(
                            coefficient*t**index for index, coefficient in enumerate(denominator)
                        )
                        return (
                            numerator_degree,
                            denominator_degree,
                            parameter_field(numerator_polynomial/denominator_polynomial),
                        )
                return None

            d_candidate = minimal_pade(d_in_parameter, 8)
            q_candidate = minimal_pade(q_in_parameter, 8)
            if d_candidate is None or q_candidate is None:
                raise AssertionError("failed to reconstruct normalized plane coordinates")
            d_function = d_candidate[2]
            q_function = q_candidate[2]
            plane_residual = sum(
                parameter_field(coefficient)
                * d_function**exponents[0]
                * q_function**exponents[1]
                for exponents, coefficient in plane_polynomial.dict().items()
            )
            if plane_residual:
                raise AssertionError("reconstructed parameterization misses plane equation")
            remaining_ring = PolynomialRing(
                parameter_field, names=("P", "E"), order="lex"
            )
            remaining_P, remaining_E = remaining_ring.gens()
            substituted_generators = []
            for polynomial in affine_generators:
                substituted = remaining_ring.zero()
                for exponents, coefficient in polynomial.dict().items():
                    substituted += (
                        parameter_field(coefficient)
                        * d_function**exponents[0]
                        * remaining_P**exponents[1]
                        * q_function**exponents[2]
                        * remaining_E**exponents[3]
                    )
                substituted_generators.append(substituted)
            remaining_ideal = remaining_ring.ideal(substituted_generators)
            remaining_groebner = remaining_ideal.groebner_basis()

            def solved_coordinate(variable, other_variable):
                for relation in remaining_groebner:
                    if relation.degree(variable) == 1 and relation.degree(other_variable) == 0:
                        coefficient = relation.monomial_coefficient(variable)
                        constant = relation-variable*coefficient
                        if constant.degree() <= 0:
                            return parameter_field(-constant/coefficient)
                raise AssertionError(f"no linear solution for {variable}")

            p_function = solved_coordinate(remaining_P, remaining_E)
            e_function = solved_coordinate(remaining_E, remaining_P)
            if any(
                relation.subs({remaining_P: p_function, remaining_E: e_function})
                for relation in substituted_generators
            ):
                raise AssertionError("recovered P,E do not satisfy the modular ideal")

            def rational_degrees(function):
                return (
                    function.numerator().degree(),
                    function.denominator().degree(),
                )

            stage(
                "plane_projective_parameter",
                time.monotonic(),
                parameter=f"({high_form})/({low_form})",
                D_degrees=d_candidate[:2],
                D=d_function,
                Q_degrees=q_candidate[:2],
                Q=q_function,
                exact_plane_substitution=1,
            )
            stage(
                "surface_projective_parameter",
                time.monotonic(),
                D_degrees=rational_degrees(d_function),
                D=d_function,
                P_degrees=rational_degrees(p_function),
                P=p_function,
                Q_degrees=rational_degrees(q_function),
                Q=q_function,
                E_degrees=rational_degrees(e_function),
                E=e_function,
                exact_ideal_substitution=1,
            )
            if arguments.parameter_output:
                parameter_artifact = {
                    "schema": "q80-cm24-formal-branch-parameter-v1",
                    "status": "bounded_mod7_formal_evidence",
                    "field": f"GF({prime})",
                    "source_ideal": str(input_path),
                    "source_ideal_sha256": hashlib.sha256(raw_bytes).hexdigest(),
                    "characteristic_zero_slope": "8/87",
                    "slope_mod_7": "5:1",
                    "parameter_at_cm24": "t=0",
                    "centered_coordinates": {
                        "D": "d-3", "P": "p-4", "Q": "q-3", "E": "e-2"
                    },
                    "normalization_hilbert_polynomial": str(normalized_hilbert_polynomial),
                    "normalization_parameter": f"({high_form})/({low_form})",
                    "functions": {
                        "D": {"degrees": list(rational_degrees(d_function)), "value": str(d_function)},
                        "P": {"degrees": list(rational_degrees(p_function)), "value": str(p_function)},
                        "Q": {"degrees": list(rational_degrees(q_function)), "value": str(q_function)},
                        "E": {"degrees": list(rational_degrees(e_function)), "value": str(e_function)},
                    },
                    "verification": {
                        "plane_substitution": True,
                        "all_14_ideal_generators": True,
                    },
                    "caveat": (
                        "This parameterizes the bounded GF(7) formal-branch candidate ideal. "
                        "It is not yet a characteristic-zero family certificate."
                    ),
                }
                parameter_output_path = Path(arguments.parameter_output)
                parameter_output_path.parent.mkdir(parents=True, exist_ok=True)
                parameter_output_path.write_text(
                    json.dumps(
                        parameter_artifact, indent=2, sort_keys=True, default=int
                    )+"\n"
                )
                print(
                    "Q80RANK19IDEAL|stage=parameter_artifact|"
                    f"output={parameter_output_path}|status=WRITTEN",
                    flush=True,
                )
    if plane_ideal.dimension() == 1 and arguments.check_genus:
        started = time.monotonic()
        plane_genus = plane_ideal.genus()
        stage(
            "plane_normalization",
            started,
            coordinates=",".join(keep_names),
            geometric_genus=plane_genus,
        )

started = time.monotonic()
if "homogeneous_generators" in artifact:
    homogeneous_ring = PolynomialRing(
        finite, names=tuple(artifact["homogeneous_variables"]), order="degrevlex"
    )
    homogeneous_generators = tuple(
        homogeneous_ring(polynomial)
        for polynomial in artifact["homogeneous_generators"]
    )
    raw_projective_ideal = homogeneous_ring.ideal(homogeneous_generators)
    z = homogeneous_ring(artifact["homogeneous_variables"][0])
    projective_ideal, saturation_exponent = raw_projective_ideal.saturation(
        homogeneous_ring.ideal(z)
    )
else:
    # A generic relation-basis artifact records affine generators only.
    # Sage homogenizes a Groebner basis here, giving the projective closure
    # rather than merely the ideal generated by separately homogenized input.
    projective_ideal = affine_ideal.homogenize()
    homogeneous_generators = tuple(projective_ideal.gens())
    saturation_exponent = "implicit"
hilbert_polynomial = projective_ideal.hilbert_polynomial(algorithm="singular")
projective_dimension = projective_ideal.dimension()-1
if projective_dimension != 1:
    raise AssertionError(f"candidate projective ideal has dimension {projective_dimension}")
degree = ZZ(hilbert_polynomial[1])
arithmetic_genus = ZZ(1-hilbert_polynomial(0))
stage(
    "projective_closure",
    started,
    raw_generators=len(homogeneous_generators),
    saturation_exponent=saturation_exponent,
    dimension=projective_dimension,
    degree=degree,
    hilbert_polynomial=hilbert_polynomial,
    arithmetic_genus=arithmetic_genus,
)

is_prime = None
if arguments.check_prime:
    started = time.monotonic()
    is_prime = projective_ideal.is_prime()
    stage("prime", started, is_prime=int(is_prime))

geometric_genus = None
if arguments.check_genus:
    if is_prime is False:
        raise ValueError("refusing genus computation on a reducible candidate ideal")
    started = time.monotonic()
    geometric_genus = projective_ideal.genus()
    stage("normalization", started, geometric_genus=geometric_genus)

structural_signature = degree == 15 and arithmetic_genus == 13
print(
    "Q80RANK19IDEAL|"
    f"input={input_path}|sha256={hashlib.sha256(raw_bytes).hexdigest()}|"
    f"degree15_pa13={int(structural_signature)}|"
    f"projective_prime={'UNCHECKED' if is_prime is None else int(is_prime)}|"
    f"projective_genus={'UNCHECKED' if geometric_genus is None else geometric_genus}|"
    f"plane_genus={'UNCHECKED' if plane_geometric_genus is None else plane_geometric_genus}|"
    f"scope=bounded_mod{prime}_candidate_ideal|status=PASS",
    flush=True,
)

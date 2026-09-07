"""Exact certificates for reciprocal-boundary suspension links.

The routines in this module are deliberately algebraic.  They do not infer a
classification from numerical sampling: valuations use exact divisibility,
boundary residue polynomials use Groebner elimination, and spectral
compatibility is a gcd in a rational-function field.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from itertools import product
from typing import Mapping, Optional

import sympy as sp


@dataclass(frozen=True)
class BoundaryReconstructionCandidate:
    """A boundary projection and one proposed primitive reconstruction residue.

    ``base_coordinates`` are functions on the source.  ``base_symbols`` are
    independent symbols used for elimination.  The parametrization describes
    the normalization chart of the base boundary over ``QQ(boundary_parameter)``.
    """

    name: str
    source_variables: tuple[sp.Symbol, ...]
    boundary: sp.Expr
    base_coordinates: tuple[sp.Expr, sp.Expr]
    reconstruction: sp.Expr
    base_symbols: tuple[sp.Symbol, sp.Symbol]
    reconstruction_symbol: sp.Symbol
    boundary_parameter: sp.Symbol
    boundary_parametrization: Mapping[sp.Symbol, sp.Expr]
    reconstruction_is_primitive: bool = True


@dataclass(frozen=True)
class FullSteinCertificate:
    """The relative constant field of the generic boundary fiber.

    The calculation is independent of the proposed reconstruction residue.
    The primitive Jacobian LND is reduced to the generic boundary fiber, a
    local slice is constructed there, and invariantization produces a finite
    presentation of the full relative constant field.
    """

    generic_fiber_equations: tuple[sp.Expr, ...]
    primitive_lnd_coefficients: tuple[sp.Expr, ...]
    induced_lnd_coefficients: tuple[sp.Expr, ...]
    local_slice_numerator: sp.Expr
    local_slice_denominator: sp.Expr
    invariant_generators: tuple[sp.Expr, ...]
    invariant_relations: tuple[sp.Expr, ...]
    primitive_element: sp.Expr
    primitive_polynomial: sp.Expr
    degree: int
    geometric_component_count: int


@dataclass(frozen=True)
class BoundaryReconstructionCertificate:
    name: str
    boundary_parametrization_valid: bool
    elimination_relations: tuple[sp.Expr, ...]
    residue_polynomial: sp.Expr
    residue_degree: int
    full_stein: FullSteinCertificate
    stein_degree: int
    reconstruction_generates_full_stein: bool
    hidden_cover_detected: bool
    primitive_marking_agrees: bool


@dataclass(frozen=True)
class ReciprocalLinkCandidate:
    """Polynomial/rational data of an orientation-reversing reciprocal link."""

    name: str
    source_variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol]
    boundary: sp.Expr
    s: sp.Expr
    P: sp.Expr
    Q: sp.Expr
    m: int
    r: int
    coordinate_symbols: tuple[sp.Symbol, sp.Symbol, sp.Symbol]
    inverse_source: Optional[tuple[sp.Expr, sp.Expr, sp.Expr]]
    boundary_parameter: sp.Symbol
    boundary_parametrization: Mapping[sp.Symbol, sp.Expr]
    reconstruction_is_primitive: bool = True


@dataclass(frozen=True)
class SpectralObstructionCertificate:
    spectral_polynomial: sp.Expr
    transformed_constraint: sp.Expr
    common_factor: sp.Expr
    common_degree: int
    excludes_candidate: bool


@dataclass(frozen=True)
class ReciprocalLinkCertificate:
    name: str
    D: sp.Expr
    Y: sp.Expr
    x: sp.Expr
    y: sp.Expr
    B: sp.Expr
    valuations: tuple[tuple[str, int], ...]
    marked_valuation_pattern: bool
    reciprocal_identity: bool
    straightening_identity: bool
    polynomial_straightening: bool
    chart_jacobian: sp.Expr
    straightening_jacobian: sp.Expr
    chart_ratio: sp.Expr
    lnd_coefficients: tuple[sp.Expr, sp.Expr, sp.Expr]
    lnd_content: sp.Expr
    lnd_on_B: sp.Expr
    localized_coordinate_certificate: bool
    boundary: BoundaryReconstructionCertificate
    spectral: SpectralObstructionCertificate
    verdict: str


def _poly_order(expression: sp.Expr, prime: sp.Expr, variables) -> int:
    expression = sp.expand(expression)
    if expression == 0:
        raise ValueError("valuation of zero is not finite")
    current = sp.Poly(expression, *variables, domain="EX")
    divisor = sp.Poly(sp.expand(prime), *variables, domain="EX")
    order = 0
    while True:
        quotient, remainder = current.div(divisor)
        if not remainder.is_zero:
            return order
        order += 1
        current = quotient


def prime_valuation(expression: sp.Expr, prime: sp.Expr, variables) -> int:
    """Return the exact prime order of a rational expression."""
    numerator, denominator = sp.fraction(sp.cancel(expression))
    return _poly_order(numerator, prime, variables) - _poly_order(
        denominator, prime, variables
    )


def _is_polynomial(expression: sp.Expr, variables) -> bool:
    _, denominator = sp.fraction(sp.cancel(expression))
    return sp.Poly(denominator, *variables, domain="EX").total_degree() == 0


def jacobian_derivation(
    first: sp.Expr,
    second: sp.Expr,
    value: sp.Expr,
    variables,
) -> sp.Expr:
    return sp.factor(
        sp.det(sp.Matrix([first, second, value]).jacobian(variables))
    )


def _monic_polynomial_over_parameter(
    expressions: list[sp.Expr],
    variable: sp.Symbol,
    parameter: sp.Symbol,
) -> sp.Poly:
    field = sp.QQ.frac_field(parameter)
    polynomials: list[sp.Poly] = []
    for expression in expressions:
        numerator, _ = sp.fraction(sp.cancel(expression))
        if numerator == 0 or variable not in numerator.free_symbols:
            continue
        polynomials.append(sp.Poly(numerator, variable, domain=field))
    if not polynomials:
        raise ValueError("elimination produced no algebraic reconstruction relation")
    generator = reduce(sp.gcd, polynomials)
    if generator.degree() < 1:
        raise ValueError("reconstruction residue is transcendental in this chart")
    return generator.monic()


def _derivation(
    expression: sp.Expr,
    coefficients: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * sp.diff(expression, variable)
            for coefficient, variable in zip(coefficients, variables)
        )
    )


def _primitive_jacobian_lnd(
    candidate: BoundaryReconstructionCandidate,
) -> tuple[sp.Expr, ...]:
    coefficients = tuple(
        jacobian_derivation(
            candidate.base_coordinates[0],
            candidate.base_coordinates[1],
            variable,
            candidate.source_variables,
        )
        for variable in candidate.source_variables
    )
    nonzero = [coefficient for coefficient in coefficients if coefficient != 0]
    if not nonzero:
        raise ValueError("the boundary projection has zero Jacobian derivation")
    content = reduce(sp.gcd, nonzero)
    primitive = tuple(sp.cancel(coefficient / content) for coefficient in coefficients)
    if not all(_is_polynomial(value, candidate.source_variables) for value in primitive):
        raise ValueError("could not divide the Jacobian derivation by polynomial content")
    return primitive


def _quotient_degree(groebner_basis: sp.GroebnerBasis) -> int:
    """Dimension of a zero-dimensional quotient from its leading monomials."""
    if not groebner_basis.is_zero_dimensional:
        raise ValueError("the invariant algebra is not finite over the boundary field")
    variable_count = len(groebner_basis.gens)
    leading = [
        polynomial.LM(order=groebner_basis.order).exponents
        for polynomial in groebner_basis.polys
    ]
    bounds: list[Optional[int]] = [None] * variable_count
    for monomial in leading:
        support = [index for index, exponent in enumerate(monomial) if exponent]
        if len(support) == 1:
            index = support[0]
            exponent = monomial[index]
            if bounds[index] is None or exponent < bounds[index]:
                bounds[index] = exponent
    if any(bound is None for bound in bounds):
        raise ValueError("zero-dimensional leading ideal has no pure-power bound")
    standard = 0
    for monomial in product(*(range(int(bound)) for bound in bounds)):
        if not any(
            all(exponent >= divisor for exponent, divisor in zip(monomial, lead))
            for lead in leading
        ):
            standard += 1
    return standard


def _primitive_element_certificate(
    invariant_relations: tuple[sp.Expr, ...],
    invariant_symbols: tuple[sp.Symbol, ...],
    invariant_expressions: tuple[sp.Expr, ...],
    parameter: sp.Symbol,
    degree: int,
) -> tuple[sp.Expr, sp.Expr]:
    primitive_symbol = sp.Symbol("Z_Stein")
    coefficient_vectors = []
    for index in range(len(invariant_symbols)):
        coefficient_vectors.append(
            tuple(1 if index == position else 0 for position in range(len(invariant_symbols)))
        )
    # For two distinct embeddings, equality on
    # Z_0+s*Z_1+...+s^(n-1)*Z_(n-1) excludes at most n-1 values of s.
    # There are at most degree*(degree-1)/2 pairs, so testing one more integer
    # than this bound is a deterministic primitive-element search in
    # characteristic zero, not a heuristic cutoff.
    bad_seed_bound = (
        degree * (degree - 1) // 2 * max(1, len(invariant_symbols) - 1)
    )
    for seed in range(1, bad_seed_bound + 2):
        coefficient_vectors.append(
            tuple(seed**position for position in range(len(invariant_symbols)))
        )

    field = sp.QQ.frac_field(parameter)
    for coefficients in coefficient_vectors:
        combination = sum(
            coefficient * symbol
            for coefficient, symbol in zip(coefficients, invariant_symbols)
        )
        basis = sp.groebner(
            [*invariant_relations, primitive_symbol - combination],
            *invariant_symbols,
            primitive_symbol,
            order="lex",
            domain=field,
        )
        eliminated = [
            polynomial.as_expr()
            for polynomial in basis.polys
            if not (polynomial.as_expr().free_symbols & set(invariant_symbols))
        ]
        try:
            minimal = _monic_polynomial_over_parameter(
                eliminated, primitive_symbol, parameter
            )
        except ValueError:
            continue
        if minimal.degree() == degree:
            expression = sp.factor(
                sum(
                    coefficient * invariant
                    for coefficient, invariant in zip(coefficients, invariant_expressions)
                )
            )
            return expression, sp.factor(minimal.as_expr())
    raise ValueError("failed to find a primitive element for the full Stein field")


def classify_full_stein_degree(
    candidate: BoundaryReconstructionCandidate,
    max_lnd_order: int = 32,
) -> FullSteinCertificate:
    """Compute the full generic Stein field, independently of ``reconstruction``.

    On the generic boundary fiber the primitive Jacobian LND has a rational
    local slice.  The associated finite exponential projection sends every
    source coordinate to an invariant.  These images generate the relative
    constant field, so the dimension of their exact relation algebra is the
    full Stein degree.
    """
    source = candidate.source_variables
    parameter = candidate.boundary_parameter
    field = sp.QQ.frac_field(parameter)
    base_x, base_y = candidate.base_symbols
    try:
        parameter_x = candidate.boundary_parametrization[base_x]
        parameter_y = candidate.boundary_parametrization[base_y]
    except KeyError as error:
        raise ValueError("boundary parametrization must specify both base symbols") from error

    generic_equations = tuple(
        sp.factor(sp.fraction(sp.cancel(expression))[0])
        for expression in (
            candidate.boundary,
            candidate.base_coordinates[0] - parameter_x,
            candidate.base_coordinates[1] - parameter_y,
        )
    )
    fiber_basis = sp.groebner(
        generic_equations,
        *source,
        order="lex",
        domain=field,
    )
    if any(polynomial.as_expr() == 1 for polynomial in fiber_basis.polys):
        raise ValueError("the proposed boundary chart has empty generic fiber")

    def reduce_on_fiber(expression: sp.Expr) -> sp.Expr:
        _, remainder = fiber_basis.reduce(sp.cancel(expression))
        return sp.factor(remainder)

    primitive = _primitive_jacobian_lnd(candidate)
    induced = tuple(reduce_on_fiber(coefficient) for coefficient in primitive)
    if all(coefficient == 0 for coefficient in induced):
        raise ValueError("the primitive LND vanishes on the generic boundary fiber")

    best_derivatives: Optional[list[sp.Expr]] = None
    for variable in source:
        derivatives = [reduce_on_fiber(variable)]
        for _ in range(max_lnd_order):
            next_derivative = reduce_on_fiber(
                _derivation(derivatives[-1], primitive, source)
            )
            if next_derivative == 0:
                break
            derivatives.append(next_derivative)
        else:
            raise ValueError(
                f"LND iteration did not terminate within order {max_lnd_order}"
            )
        if len(derivatives) > 1 and (
            best_derivatives is None or len(derivatives) > len(best_derivatives)
        ):
            best_derivatives = derivatives
    if best_derivatives is None:
        raise ValueError("could not find a local slice on the generic boundary fiber")

    slice_numerator = best_derivatives[-2]
    slice_denominator = best_derivatives[-1]
    local_slice = sp.cancel(slice_numerator / slice_denominator)

    invariant_expressions: list[sp.Expr] = []
    for variable in source:
        terms = [reduce_on_fiber(variable)]
        for _ in range(max_lnd_order):
            derivative = reduce_on_fiber(_derivation(terms[-1], primitive, source))
            if derivative == 0:
                break
            terms.append(derivative)
        invariant = sum(
            (-local_slice) ** order * derivative / sp.factorial(order)
            for order, derivative in enumerate(terms)
        )
        invariant = sp.factor(sp.cancel(invariant))
        if reduce_on_fiber(_derivation(invariant, primitive, source)) != 0:
            raise ValueError("local-slice invariantization failed its exact check")
        invariant_expressions.append(invariant)

    invariant_symbols = tuple(
        sp.Symbol(f"Z_Stein_{index}") for index in range(len(source))
    )
    inverse_symbol = sp.Symbol("H_Stein_inverse")
    graph_equations = list(generic_equations)
    for symbol, invariant in zip(invariant_symbols, invariant_expressions):
        numerator, denominator = sp.fraction(sp.cancel(invariant))
        graph_equations.append(sp.expand(symbol * denominator - numerator))
    graph_equations.append(inverse_symbol * slice_denominator - 1)
    graph_basis = sp.groebner(
        graph_equations,
        *source,
        inverse_symbol,
        *invariant_symbols,
        order="lex",
        domain=field,
    )
    eliminated_variables = set(source) | {inverse_symbol}
    invariant_relations = tuple(
        sp.factor(polynomial.as_expr())
        for polynomial in graph_basis.polys
        if not (polynomial.as_expr().free_symbols & eliminated_variables)
    )
    if not invariant_relations:
        raise ValueError("elimination produced no full-Stein invariant relations")
    invariant_basis = sp.groebner(
        invariant_relations,
        *invariant_symbols,
        order="lex",
        domain=field,
    )
    degree = _quotient_degree(invariant_basis)
    primitive_element, primitive_polynomial = _primitive_element_certificate(
        invariant_relations,
        invariant_symbols,
        tuple(invariant_expressions),
        parameter,
        degree,
    )
    return FullSteinCertificate(
        generic_fiber_equations=generic_equations,
        primitive_lnd_coefficients=tuple(sp.factor(value) for value in primitive),
        induced_lnd_coefficients=induced,
        local_slice_numerator=sp.factor(slice_numerator),
        local_slice_denominator=sp.factor(slice_denominator),
        invariant_generators=tuple(invariant_expressions),
        invariant_relations=invariant_relations,
        primitive_element=primitive_element,
        primitive_polynomial=primitive_polynomial,
        degree=degree,
        geometric_component_count=degree,
    )


def classify_boundary_reconstruction(
    candidate: BoundaryReconstructionCandidate,
) -> BoundaryReconstructionCertificate:
    """Eliminate source coordinates and compute the residue/Stein degree."""
    source = candidate.source_variables
    base_x, base_y = candidate.base_symbols
    residue = candidate.reconstruction_symbol
    equations = [
        candidate.boundary,
        base_x - candidate.base_coordinates[0],
        base_y - candidate.base_coordinates[1],
        residue - candidate.reconstruction,
    ]
    basis = sp.groebner(
        equations,
        *source,
        base_x,
        base_y,
        residue,
        order="lex",
    )
    source_set = set(source)
    elimination_relations = tuple(
        sp.factor(polynomial.as_expr())
        for polynomial in basis.polys
        if not (polynomial.as_expr().free_symbols & source_set)
    )
    parametrized = [
        sp.cancel(relation.subs(candidate.boundary_parametrization))
        for relation in elimination_relations
    ]
    residue_polynomial = _monic_polynomial_over_parameter(
        parametrized,
        residue,
        candidate.boundary_parameter,
    )

    # Every eliminated base relation must vanish on the proposed boundary
    # normalization chart.
    base_only = [
        relation
        for relation in elimination_relations
        if residue not in relation.free_symbols
    ]
    parameter_valid = all(
        sp.cancel(relation.subs(candidate.boundary_parametrization)) == 0
        for relation in base_only
    )
    degree = residue_polynomial.degree()
    full_stein = classify_full_stein_degree(candidate)
    generates_full_stein = degree == full_stein.degree
    return BoundaryReconstructionCertificate(
        name=candidate.name,
        boundary_parametrization_valid=parameter_valid,
        elimination_relations=elimination_relations,
        residue_polynomial=sp.factor(residue_polynomial.as_expr()),
        residue_degree=degree,
        full_stein=full_stein,
        stein_degree=full_stein.degree,
        reconstruction_generates_full_stein=generates_full_stein,
        hidden_cover_detected=not generates_full_stein,
        primitive_marking_agrees=(
            candidate.reconstruction_is_primitive == generates_full_stein
        ),
    )


def spectral_polynomial(m: int, r: int, variable: sp.Symbol) -> sp.Poly:
    """Return J_(mr,r) using exact beta-integral coefficients."""
    if m < 1 or r < 1:
        raise ValueError("m and r must be positive")
    degree = m * r
    expression = sum(
        (-variable) ** j
        * sp.binomial(degree, j)
        * sp.factorial(r)
        * sp.factorial(j)
        / sp.factorial(r + j + 1)
        for j in range(degree + 1)
    )
    return sp.Poly(expression, variable, domain=sp.QQ)


def spectral_obstruction(
    residue_polynomial: sp.Expr,
    residue_symbol: sp.Symbol,
    boundary_parameter: sp.Symbol,
    boundary_y: sp.Expr,
    m: int,
    r: int,
) -> SpectralObstructionCertificate:
    """Compare a boundary minimal polynomial with the Keller leading equation."""
    q = sp.Symbol("q_spectral")
    spectral = spectral_polynomial(m, r, q)
    substituted = sp.cancel(
        spectral.as_expr().subs(q, residue_symbol / boundary_y ** (m + 1))
    )
    numerator, _ = sp.fraction(substituted)
    field = sp.QQ.frac_field(boundary_parameter)
    residue_poly = sp.Poly(residue_polynomial, residue_symbol, domain=field).monic()
    constraint_poly = sp.Poly(numerator, residue_symbol, domain=field).monic()
    common = sp.gcd(residue_poly, constraint_poly).monic()
    common_degree = common.degree()
    # Over the algebraically closed ground field the transformed spectral
    # polynomial is a product of linear factors
    # Z-q_i*boundary_y^(m+1).  Thus a primitive geometric residue extension
    # of degree >1 is excluded even if several conjugate factors remain
    # bundled over QQ(parameter).
    geometrically_excluded = residue_poly.degree() > 1 or common_degree == 0
    return SpectralObstructionCertificate(
        spectral_polynomial=sp.factor(spectral.as_expr()),
        transformed_constraint=sp.factor(constraint_poly.as_expr()),
        common_factor=sp.factor(common.as_expr()),
        common_degree=common_degree,
        excludes_candidate=geometrically_excluded,
    )


def _localized_inverse_certificate(candidate, x, y, B) -> bool:
    if candidate.inverse_source is None:
        return False
    coordinate_symbols = candidate.coordinate_symbols
    forward_substitution = dict(zip(coordinate_symbols, (x, y, B)))
    first_composition = all(
        sp.cancel(inverse.subs(forward_substitution) - source) == 0
        for inverse, source in zip(candidate.inverse_source, candidate.source_variables)
    )
    inverse_substitution = dict(zip(candidate.source_variables, candidate.inverse_source))
    second_composition = all(
        sp.cancel(forward.subs(inverse_substitution) - coordinate) == 0
        for forward, coordinate in zip((x, y, B), coordinate_symbols)
    )
    return first_composition and second_composition


def classify_reciprocal_link(
    candidate: ReciprocalLinkCandidate,
) -> ReciprocalLinkCertificate:
    """Return an exact certificate for one reciprocal-link candidate."""
    variables = candidate.source_variables
    A = candidate.boundary
    Y = sp.cancel(candidate.Q - candidate.s * candidate.P)
    D = sp.cancel(1 - candidate.s * Y**candidate.m)
    x = sp.cancel(A * candidate.s)
    y = Y
    B = sp.cancel(candidate.P / A)

    valuations = (
        ("s", prime_valuation(candidate.s, A, variables)),
        ("Y", prime_valuation(Y, A, variables)),
        ("P", prime_valuation(candidate.P, A, variables)),
        ("B", prime_valuation(B, A, variables)),
        ("D", prime_valuation(D, A, variables)),
    )
    polynomial_straightening = all(
        _is_polynomial(value, variables) for value in (x, y, B)
    )
    reciprocal_identity = sp.cancel(D - 1 / A) == 0
    straightening_identity = all(
        sp.cancel(value) == 0
        for value in (
            A - 1 - x * y**candidate.m,
            candidate.P - A * B,
            candidate.Q - y - x * B,
        )
    )

    chart_jacobian = sp.factor(
        sp.det(sp.Matrix([candidate.s, candidate.P, candidate.Q]).jacobian(variables))
    )
    straightening_jacobian = sp.factor(
        sp.det(sp.Matrix([x, y, B]).jacobian(variables))
    )
    chart_ratio = sp.factor(sp.cancel(chart_jacobian / straightening_jacobian))
    coefficients = tuple(
        jacobian_derivation(x, y, variable, variables) for variable in variables
    )
    nonzero_coefficients = [coefficient for coefficient in coefficients if coefficient != 0]
    content = reduce(sp.gcd, nonzero_coefficients) if nonzero_coefficients else sp.Integer(0)
    lnd_on_B = jacobian_derivation(x, y, B, variables)
    localized_certificate = _localized_inverse_certificate(candidate, x, y, B)

    base_x, base_y, residue = candidate.coordinate_symbols
    boundary_candidate = BoundaryReconstructionCandidate(
        name=candidate.name,
        source_variables=variables,
        boundary=A,
        base_coordinates=(x, y),
        reconstruction=B,
        base_symbols=(base_x, base_y),
        reconstruction_symbol=residue,
        boundary_parameter=candidate.boundary_parameter,
        boundary_parametrization=candidate.boundary_parametrization,
        reconstruction_is_primitive=candidate.reconstruction_is_primitive,
    )
    boundary_certificate = classify_boundary_reconstruction(boundary_candidate)
    boundary_y = sp.cancel(
        candidate.boundary_parametrization[base_y]
    )
    spectral_certificate = spectral_obstruction(
        boundary_certificate.residue_polynomial,
        residue,
        candidate.boundary_parameter,
        boundary_y,
        candidate.m,
        candidate.r,
    )

    marked_pattern = dict(valuations) == {
        "s": -1,
        "Y": 0,
        "P": 1,
        "B": 0,
        "D": -1,
    }
    reciprocal_valid = (
        marked_pattern
        and reciprocal_identity
        and straightening_identity
        and polynomial_straightening
        and localized_certificate
        and boundary_certificate.boundary_parametrization_valid
    )
    if not reciprocal_valid:
        verdict = "invalid_or_incomplete_reciprocal_certificate"
    elif spectral_certificate.excludes_candidate:
        verdict = "excluded_by_keller_boundary_spectrum"
    elif boundary_certificate.hidden_cover_detected:
        verdict = "hidden_stein_cover_detected"
    elif boundary_certificate.stein_degree == 1:
        verdict = "passes_marked_cancellation_boundary_prefilter"
    else:
        verdict = "requires_primitive_stein_analysis"

    return ReciprocalLinkCertificate(
        name=candidate.name,
        D=D,
        Y=Y,
        x=x,
        y=y,
        B=B,
        valuations=valuations,
        marked_valuation_pattern=marked_pattern,
        reciprocal_identity=reciprocal_identity,
        straightening_identity=straightening_identity,
        polynomial_straightening=polynomial_straightening,
        chart_jacobian=chart_jacobian,
        straightening_jacobian=straightening_jacobian,
        chart_ratio=chart_ratio,
        lnd_coefficients=coefficients,
        lnd_content=sp.factor(content),
        lnd_on_B=lnd_on_B,
        localized_coordinate_certificate=localized_certificate,
        boundary=boundary_certificate,
        spectral=spectral_certificate,
        verdict=verdict,
    )


def standard_cancellation_example() -> ReciprocalLinkCandidate:
    """Return the smallest complete reciprocal candidate, with q=3."""
    u, v, w = sp.symbols("u v w")
    X, Y, Z = sp.symbols("X Y Z")
    tau = sp.symbols("tau", nonzero=True)
    A = 1 + u * v
    B = A**2 * w + 3 * v**2
    return ReciprocalLinkCandidate(
        name="cancellation_(1,1,q=3)",
        source_variables=(u, v, w),
        boundary=A,
        s=u / A,
        P=A * B,
        Q=v + u * B,
        m=1,
        r=1,
        coordinate_symbols=(X, Y, Z),
        inverse_source=(X, Y, (Z - 3 * Y**2) / (1 + X * Y) ** 2),
        boundary_parameter=tau,
        boundary_parametrization={X: -1 / tau, Y: tau},
    )


def split_plinth_example(
    degree: int,
    *,
    hidden_reconstruction: bool = False,
    plinth_exponent: int = 2,
) -> BoundaryReconstructionCandidate:
    """Return ``u^n v+w^f`` with full boundary Stein degree ``f``."""
    if degree < 2:
        raise ValueError("the split-plinth degree must be at least two")
    if plinth_exponent < 1:
        raise ValueError("the plinth exponent must be positive")
    u, v, w = sp.symbols("u v w")
    X, Y, Z = sp.symbols("X Y Z")
    tau = sp.symbols("tau")
    invariant = u**plinth_exponent * v + w**degree
    reconstruction = invariant if hidden_reconstruction else w
    return BoundaryReconstructionCandidate(
        name=(
            f"split_plinth_hidden_degree_{degree}"
            if hidden_reconstruction
            else f"split_plinth_degree_{degree}"
        ),
        source_variables=(u, v, w),
        boundary=u,
        base_coordinates=(u, invariant),
        reconstruction=reconstruction,
        base_symbols=(X, Y),
        reconstruction_symbol=Z,
        boundary_parameter=tau,
        boundary_parametrization={X: 0, Y: tau},
        reconstruction_is_primitive=not hidden_reconstruction,
    )


def masuda_plinth_example() -> BoundaryReconstructionCandidate:
    """Return Masuda's standard degree-two plinth/Stein countermodel."""
    return split_plinth_example(2)


def masuda_hidden_cover_example() -> BoundaryReconstructionCandidate:
    """Return the Masuda fiber while deliberately displaying only its base residue.

    The proposed reconstruction has degree one, but the induced boundary LND
    has the hidden constant ``w`` with ``w**2=tau``.  This is the regression
    model showing why full Stein degree cannot be inferred from one residue.
    """
    return split_plinth_example(2, hidden_reconstruction=True)

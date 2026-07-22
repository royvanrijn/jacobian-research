"""Exact certificates for reciprocal-boundary suspension links.

The routines in this module are deliberately algebraic.  They do not infer a
classification from numerical sampling: valuations use exact divisibility,
boundary residue polynomials use Groebner elimination, and spectral
compatibility is a gcd in a rational-function field.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
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
class BoundaryReconstructionCertificate:
    name: str
    boundary_parametrization_valid: bool
    elimination_relations: tuple[sp.Expr, ...]
    residue_polynomial: sp.Expr
    residue_degree: int
    stein_degree: Optional[int]


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
    return BoundaryReconstructionCertificate(
        name=candidate.name,
        boundary_parametrization_valid=parameter_valid,
        elimination_relations=elimination_relations,
        residue_polynomial=sp.factor(residue_polynomial.as_expr()),
        residue_degree=degree,
        stein_degree=degree if candidate.reconstruction_is_primitive else None,
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


def masuda_plinth_example() -> BoundaryReconstructionCandidate:
    """Return the standard degree-two plinth/Stein countermodel."""
    u, v, w = sp.symbols("u v w")
    X, Y, Z = sp.symbols("X Y Z")
    tau = sp.symbols("tau")
    invariant = u**2 * v + w**2
    return BoundaryReconstructionCandidate(
        name="Masuda_degree_two_plinth",
        source_variables=(u, v, w),
        boundary=u,
        base_coordinates=(u, invariant),
        reconstruction=w,
        base_symbols=(X, Y),
        reconstruction_symbol=Z,
        boundary_parameter=tau,
        boundary_parametrization={X: 0, Y: tau},
    )

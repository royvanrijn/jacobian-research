"""Finite invariant pipeline for minimal-boundary package predicates.

The module deliberately separates two layers:

* ``FiniteNormalizationRecord`` is an exact, finite export of canonical
  normalization data (prime signatures, valuation rows, collision modules,
  and a marked quotient chart).
* ``extract_minimal_boundary`` uses only that export.  It never dispatches
  on a construction or fixture name.

The export is not computed from arbitrary polynomial maps here.  Constructing
it from a bare Keller map is the remaining normalization/frontend problem.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from functools import reduce
from math import gcd
from typing import Iterable

import sympy as sp


class Outcome(str, Enum):
    """Three-valued result of a finite predicate test."""

    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class PredicateResult:
    """One predicate result with a machine-readable certificate."""

    predicate: str
    outcome: Outcome
    reason: str
    certificate: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class PrimeRecord:
    """Intrinsic height-one signature from the canonical normalization."""

    label: str
    boundary: bool
    critical_fitting: bool
    quotient_image: str
    target_image: str
    ramification_index: int
    residue_degree: int
    different_exponent: int
    color: str
    completed_incidence: tuple[int, ...]
    automorphism_orbit_size: int = 1

    @property
    def signature(self) -> tuple[object, ...]:
        return (
            self.color,
            self.ramification_index,
            self.residue_degree,
            self.different_exponent,
            self.completed_incidence,
        )


@dataclass(frozen=True)
class LinkRecord:
    """Finite rank-one localization-link audit."""

    label: str
    normal_source: bool
    normal_target: bool
    source_ufd: bool
    target_ufd: bool
    scalar_source_units: bool
    scalar_target_units: bool
    source_parameter_prime: bool
    target_parameter_prime: bool
    unit_lattice_rank: int
    target_parameter_support_in_source: tuple[str, ...]
    source_parameter_support_in_target: tuple[str, ...]
    orientation_exponent: int
    forward_generator_valuations: tuple[int, ...] = ()
    reverse_generator_valuations: tuple[int, ...] = ()


@dataclass(frozen=True)
class LedgerPrime:
    """One graph-prime row in the signed relative-canonical ledger."""

    label: str
    exceptional: bool
    localization_boundary: bool
    ramification_qx: int
    pulled_back_core_ramification: int
    pulled_back_target_ramification: int

    @property
    def signed_coefficient(self) -> int:
        return (
            self.ramification_qx
            + self.pulled_back_core_ramification
            - self.pulled_back_target_ramification
        )

    @property
    def relevant(self) -> bool:
        return self.exceptional or self.localization_boundary or any(
            (
                self.ramification_qx,
                self.pulled_back_core_ramification,
                self.pulled_back_target_ramification,
            )
        )


@dataclass(frozen=True)
class CriticalCurveRecord:
    """Normalization and puncture-unit data of the selected critical curve."""

    geometrically_integral: bool
    smooth: bool
    genus: int
    parameter: sp.Symbol
    punctures: tuple[str, ...]
    unit_generators: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class CollisionRecord:
    """Closed-collision nilradical generation audit."""

    label: str
    nilradical_dimension: int
    conormal_generators: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class ConormalRecord:
    """Primitive conormal content and residue marking."""

    height_one_content_vectors: tuple[tuple[int, ...], ...]
    collisions: tuple[CollisionRecord, ...]
    residue_coefficient: sp.Expr


@dataclass(frozen=True)
class ChartRecord:
    """Marked quotient chart used only by the coefficient straightener."""

    source_variables: tuple[sp.Symbol, ...]
    target_expressions: tuple[sp.Expr, ...]
    controlled_divisor: sp.Expr
    controlled_exponent: int
    reciprocal_chart: bool
    fitting_support_rows: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class FiniteNormalizationRecord:
    """Canonical finite normalization exported to the invariant pipeline."""

    name: str
    primes: tuple[PrimeRecord, ...]
    links: tuple[LinkRecord, ...]
    ledger: tuple[LedgerPrime, ...]
    unrecorded_graph_primes: tuple[str, ...]
    critical_curve: CriticalCurveRecord
    conormal: ConormalRecord
    residue_parameter: sp.Symbol
    residue_generators: tuple[sp.Expr, ...]
    chart: ChartRecord | None


@dataclass(frozen=True)
class ChartCertificate:
    """Formula-blind output of the marked coefficient straightener."""

    mechanism: str
    core_jacobian_unit: sp.Expr
    extracted_mark: sp.Expr
    target_shear: sp.Expr
    mbp_chart: bool


@dataclass(frozen=True)
class ExtractionReport:
    """Complete output of one invariant extraction run."""

    record_name: str
    selected_prime: str | None
    puncture_rank: int | None
    valuation_rows: tuple[tuple[int, ...], ...]
    chart_certificate: ChartCertificate | None
    predicates: tuple[PredicateResult, ...]

    @property
    def mbpkg(self) -> Outcome:
        outcomes = {result.outcome for result in self.predicates}
        if Outcome.FAIL in outcomes:
            return Outcome.FAIL
        if Outcome.UNKNOWN in outcomes:
            return Outcome.UNKNOWN
        return Outcome.PASS


def _certificate(**entries: object) -> tuple[tuple[str, str], ...]:
    return tuple((key, str(value)) for key, value in sorted(entries.items()))


def _integer_rank(rows: Iterable[Iterable[int]]) -> int:
    rows = tuple(tuple(int(value) for value in row) for row in rows)
    return 0 if not rows else int(sp.Matrix(rows).rank())


def _primitive_integer_row(row: Iterable[int]) -> bool:
    values = tuple(abs(int(value)) for value in row if int(value) != 0)
    return bool(values) and reduce(gcd, values) == 1


def _polynomial_order(
    polynomial: sp.Poly, factor: sp.Poly
) -> tuple[int, sp.Poly]:
    if polynomial.is_zero:
        raise ValueError("valuation of the zero function is not finite")
    order = 0
    quotient = polynomial
    while quotient.rem(factor).is_zero:
        quotient = quotient.exquo(factor)
        order += 1
    return order, quotient


def _valuation_at_puncture(
    expression: sp.Expr, parameter: sp.Symbol, puncture: str
) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    numerator_poly = sp.Poly(numerator, parameter)
    denominator_poly = sp.Poly(denominator, parameter)
    if puncture == "infinity":
        return denominator_poly.degree() - numerator_poly.degree()
    if puncture == "zero":
        point = sp.Integer(0)
    elif puncture.startswith("finite:"):
        point = sp.sympify(puncture.removeprefix("finite:"))
    else:
        raise ValueError(f"unsupported puncture label {puncture!r}")
    factor = sp.Poly(parameter - point, parameter)
    numerator_order, _ = _polynomial_order(numerator_poly, factor)
    denominator_order, _ = _polynomial_order(denominator_poly, factor)
    return numerator_order - denominator_order


def _valuation_row(
    expression: sp.Expr,
    parameter: sp.Symbol,
    punctures: tuple[str, ...],
) -> tuple[int, ...]:
    return tuple(
        _valuation_at_puncture(expression, parameter, puncture)
        for puncture in punctures
    )


def select_critical_prime(
    record: FiniteNormalizationRecord,
) -> tuple[PrimeRecord | None, PredicateResult]:
    """Select the SCB prime from intrinsic signatures and automorphism orbits."""

    candidates = tuple(
        prime
        for prime in record.primes
        if prime.boundary and prime.critical_fitting and prime.quotient_image
    )
    fixed = tuple(
        prime for prime in candidates if prime.automorphism_orbit_size == 1
    )
    if len(fixed) != 1:
        return None, PredicateResult(
            "SCB",
            Outcome.FAIL,
            "the critical boundary signature is not uniquely fixed",
            _certificate(
                candidates=",".join(prime.label for prime in candidates),
                fixed=",".join(prime.label for prime in fixed),
            ),
        )
    selected = fixed[0]
    same_image = tuple(
        prime
        for prime in candidates
        if prime.quotient_image == selected.quotient_image
    )
    if len(same_image) != 1:
        return None, PredicateResult(
            "SCB",
            Outcome.FAIL,
            "more than one critical boundary lies over the selected core image",
            _certificate(
                image=selected.quotient_image,
                primes=",".join(prime.label for prime in same_image),
            ),
        )
    return selected, PredicateResult(
        "SCB",
        Outcome.PASS,
        "one fixed critical boundary signature is isolated",
        _certificate(
            prime=selected.label,
            signature=selected.signature,
            quotient_image=selected.quotient_image,
        ),
    )


def compute_punctures(
    curve: CriticalCurveRecord,
) -> tuple[int | None, tuple[tuple[int, ...], ...], PredicateResult]:
    """Compute puncture rank and validate the exported unit valuation lattice."""

    puncture_count = len(curve.punctures)
    try:
        rows = tuple(
            _valuation_row(unit, curve.parameter, curve.punctures)
            for unit in curve.unit_generators
        )
    except (ValueError, sp.PolynomialError) as error:
        return None, (), PredicateResult(
            "PR<=1", Outcome.FAIL, f"puncture valuation failed: {error}"
        )
    if any(sum(row) != 0 for row in rows):
        return None, rows, PredicateResult(
            "PR<=1",
            Outcome.FAIL,
            "a principal unit row violates the degree-zero relation",
        )
    rank = _integer_rank(rows)
    expected_rank = max(0, puncture_count - 1)
    complete = rank == expected_rank
    geometric = (
        curve.geometrically_integral and curve.smooth and curve.genus == 0
    )
    outcome = (
        Outcome.PASS
        if geometric and complete and rank <= 1
        else Outcome.FAIL
    )
    reason = (
        "smooth rational normalization has a complete rank-at-most-one "
        "puncture lattice"
        if outcome == Outcome.PASS
        else "critical normalization or puncture-unit lattice fails PR<=1"
    )
    return rank, rows, PredicateResult(
        "PR<=1",
        outcome,
        reason,
        _certificate(
            punctures=curve.punctures,
            computed_rank=rank,
            expected_rank=expected_rank,
            row_rank=rank,
        ),
    )


def test_saturation(links: tuple[LinkRecord, ...]) -> PredicateResult:
    """Test SAT, including rank-one orientation extracted from unit rows."""

    failures: list[str] = []
    orientations: list[int] = []
    for link in links:
        structural = all(
            (
                link.normal_source,
                link.normal_target,
                link.source_ufd,
                link.target_ufd,
                link.scalar_source_units,
                link.scalar_target_units,
                link.source_parameter_prime,
                link.target_parameter_prime,
            )
        )
        supports = (
            link.target_parameter_support_in_source == ("a",)
            and link.source_parameter_support_in_target == ("d",)
        )
        oriented = (
            link.unit_lattice_rank == 1
            and link.orientation_exponent in (-1, 1)
        )
        if not (structural and supports and oriented):
            failures.append(link.label)
        orientations.append(link.orientation_exponent)
    outcome = Outcome.FAIL if failures or not links else Outcome.PASS
    return PredicateResult(
        "SAT",
        outcome,
        (
            "every localization link is saturated with intrinsic orientation"
            if outcome == Outcome.PASS
            else "one or more localization links fail saturation"
        ),
        _certificate(
            failed=",".join(failures),
            orientations=tuple(orientations),
        ),
    )


def test_boundary_monotonicity(
    links: tuple[LinkRecord, ...],
) -> PredicateResult:
    """Test BM on the finite algebra-generator valuation lists."""

    failures: list[str] = []
    rows: list[tuple[int, ...]] = []
    for link in links:
        values = (
            link.forward_generator_valuations
            if link.orientation_exponent == 1
            else link.reverse_generator_valuations
        )
        rows.append(values)
        if any(value < 0 for value in values):
            failures.append(link.label)
    outcome = Outcome.FAIL if failures or not links else Outcome.PASS
    return PredicateResult(
        "BM",
        outcome,
        (
            "all positive-oriented chart generators have nonnegative boundary value"
            if outcome == Outcome.PASS
            else "a positive-oriented generator has a boundary pole"
        ),
        _certificate(failed=",".join(failures), generator_rows=tuple(rows)),
    )


def test_ledger(
    ledger: tuple[LedgerPrime, ...],
    unrecorded: tuple[str, ...],
) -> PredicateResult:
    """Test complete signed cancellation and irredundancy of graph primes."""

    nonzero = tuple(
        row.label for row in ledger if row.signed_coefficient != 0
    )
    irrelevant = tuple(row.label for row in ledger if not row.relevant)
    outcome = (
        Outcome.PASS
        if ledger and not nonzero and not irrelevant and not unrecorded
        else Outcome.FAIL
    )
    return PredicateResult(
        "LC",
        outcome,
        (
            "the recorded graph-prime ledger is balanced and irredundant"
            if outcome == Outcome.PASS
            else "the graph-prime ledger is unbalanced, incomplete, or redundant"
        ),
        _certificate(
            signed_rows=tuple(
                (row.label, row.signed_coefficient) for row in ledger
            ),
            nonzero=nonzero,
            irrelevant=irrelevant,
            unrecorded=unrecorded,
        ),
    )


def test_primitive_conormal(
    conormal: ConormalRecord,
    puncture_rank: int | None,
    curve: CriticalCurveRecord,
) -> PredicateResult:
    """Test height-one, closed-collision, and residue primitivity."""

    height_one = bool(conormal.height_one_content_vectors) and all(
        _primitive_integer_row(vector)
        for vector in conormal.height_one_content_vectors
    )
    collision_ranks = tuple(
        (
            collision.label,
            (
                0
                if not collision.conormal_generators
                else int(sp.Matrix(collision.conormal_generators).rank())
            ),
            collision.nilradical_dimension,
        )
        for collision in conormal.collisions
    )
    collision_failures = tuple(
        label
        for label, generated, nilradical in collision_ranks
        if generated != nilradical
    )
    try:
        residue_valuations = _valuation_row(
            conormal.residue_coefficient,
            curve.parameter,
            curve.punctures,
        )
    except (ValueError, sp.PolynomialError):
        residue_valuations = ()
    if puncture_rank == 0:
        numerator, denominator = sp.fraction(
            sp.cancel(conormal.residue_coefficient)
        )
        residue = (
            not sp.sympify(denominator).has(curve.parameter)
            and sp.Poly(numerator, curve.parameter).degree() == 1
        )
    elif puncture_rank == 1:
        residue = (
            len(residue_valuations) == len(curve.punctures)
            and sum(residue_valuations) == 0
            and _primitive_integer_row(residue_valuations)
        )
    else:
        residue = False
    outcome = (
        Outcome.PASS
        if height_one and not collision_failures and residue
        else Outcome.FAIL
    )
    return PredicateResult(
        "PC",
        outcome,
        (
            "the conormal class is primitive in height one and at collisions"
            if outcome == Outcome.PASS
            else "the conormal content, a collision, or its residue mark is imprimitive"
        ),
        _certificate(
            height_one_content=conormal.height_one_content_vectors,
            collision_ranks=collision_ranks,
            failed_collisions=collision_failures,
            residue_coefficient=conormal.residue_coefficient,
            residue_valuations=residue_valuations,
        ),
    )


def test_noncontraction(
    parameter: sp.Symbol,
    generators: tuple[sp.Expr, ...],
) -> PredicateResult:
    """Test NC by exact differentiation in the rational residue field k(t)."""

    varying = tuple(
        sp.factor(generator)
        for generator in generators
        if sp.factor(sp.diff(generator, parameter)) != 0
    )
    outcome = Outcome.PASS if varying else Outcome.FAIL
    return PredicateResult(
        "NC",
        outcome,
        (
            "the residue image has transcendence degree one"
            if outcome == Outcome.PASS
            else "every exported residue generator is constant"
        ),
        _certificate(
            parameter=parameter,
            varying_generators=varying,
        ),
    )


def _constant_quotient(numerator: sp.Expr, denominator: sp.Expr) -> sp.Expr | None:
    quotient = sp.cancel(numerator / denominator)
    if sp.denom(quotient) != 1 or quotient == 0:
        return None
    return quotient if not quotient.free_symbols else None


def _straighten_positive(chart: ChartRecord) -> ChartCertificate | None:
    if len(chart.source_variables) != 2 or len(chart.target_expressions) != 2:
        return None
    w, q = chart.source_variables
    first, second = chart.target_expressions
    if sp.expand(first - q) != 0:
        return None
    jacobian = sp.factor(
        sp.Matrix(chart.target_expressions)
        .jacobian(chart.source_variables)
        .det()
    )
    unit = _constant_quotient(
        jacobian, chart.controlled_divisor**chart.controlled_exponent
    )
    if unit is None or chart.controlled_exponent != 1:
        return None
    divisor_poly = sp.Poly(chart.controlled_divisor, q)
    if divisor_poly.degree() != 1:
        return None
    leading = divisor_poly.coeff_monomial(q)
    h = sp.factor(
        -divisor_poly.coeff_monomial(1) / leading
    )
    normalized = sp.expand(second / (-unit * leading))
    primitive = sp.integrate(h, w)
    shear = sp.factor(normalized - (w * q - primitive))
    if sp.diff(shear, w) != 0:
        return None
    return ChartCertificate(
        mechanism="positive-section",
        core_jacobian_unit=unit,
        extracted_mark=w,
        target_shear=shear,
        mbp_chart=True,
    )


def _fitting_support_rank(chart: ChartRecord) -> int:
    return _integer_rank(chart.fitting_support_rows)


def _straighten_reciprocal(chart: ChartRecord) -> ChartCertificate | None:
    if len(chart.source_variables) != 3 or len(chart.target_expressions) != 3:
        return None
    P, S, Q = chart.source_variables
    first, second, third = chart.target_expressions
    if sp.expand(first - P) != 0:
        return None
    jacobian = sp.factor(
        sp.Matrix(chart.target_expressions)
        .jacobian(chart.source_variables)
        .det()
    )
    unit = _constant_quotient(
        jacobian, chart.controlled_divisor**chart.controlled_exponent
    )
    if unit is None:
        return None

    # Cancellation-type integral core: the second target is the preserved
    # Q-coordinate and d(third)/dS is a unit times D^e.
    if sp.expand(second - Q) == 0:
        derivative_unit = _constant_quotient(
            sp.diff(third, S),
            chart.controlled_divisor**chart.controlled_exponent,
        )
        if derivative_unit is None:
            return None
        shear = sp.factor(third.subs(S, 0))
        return ChartCertificate(
            mechanism="reciprocal-integral",
            core_jacobian_unit=unit,
            extracted_mark=S,
            target_shear=shear,
            mbp_chart=_fitting_support_rank(chart) == 1,
        )

    # Marked-line incidence: change from Q to the second target B and read
    # C=Y-BX coefficientwise.  No family formula is supplied.
    beta = sp.factor(second - Q)
    if sp.diff(beta, Q) != 0:
        return None
    B = sp.Dummy("B")
    third_in_B = sp.expand(third.subs(Q, B - beta))
    polynomial = sp.Poly(third_in_B, B)
    if polynomial.degree() > 1:
        return None
    mark = sp.factor(-polynomial.coeff_monomial(B))
    base = sp.factor(polynomial.coeff_monomial(1))
    if sp.diff(mark, P) != 0 or sp.diff(mark, Q) != 0:
        return None
    if sp.expand(third_in_B - (base - B * mark)) != 0:
        return None
    return ChartCertificate(
        mechanism="quadratic-incidence",
        core_jacobian_unit=unit,
        extracted_mark=mark,
        target_shear=base,
        # Definition 1.9 currently names only weighted and cancellation
        # charts.  A quadratic-incidence extraction is recorded but is not
        # silently promoted to MBP chart-straightenability.
        mbp_chart=False,
    )


def attempt_chart_straightening(
    chart: ChartRecord | None,
    links: tuple[LinkRecord, ...],
) -> tuple[ChartCertificate | None, PredicateResult]:
    """Attempt CS by exact coefficient extraction, never by fixture name."""

    if chart is None:
        return None, PredicateResult(
            "CS",
            Outcome.UNKNOWN,
            "no marked quotient chart was exported",
        )
    orientations = tuple(link.orientation_exponent for link in links)
    positive = bool(orientations) and all(
        value == 1 for value in orientations
    )
    if positive == chart.reciprocal_chart:
        return None, PredicateResult(
            "CS",
            Outcome.FAIL,
            "the chart type disagrees with the intrinsic link orientation",
            _certificate(
                orientations=orientations,
                reciprocal_chart=chart.reciprocal_chart,
            ),
        )
    certificate = (
        _straighten_positive(chart)
        if positive
        else _straighten_reciprocal(chart)
    )
    if certificate is None:
        return None, PredicateResult(
            "CS",
            Outcome.FAIL,
            "the coefficient straightener found no admissible marked chart",
            _certificate(orientations=orientations),
        )
    outcome = Outcome.PASS if certificate.mbp_chart else Outcome.FAIL
    reason = (
        "coefficient extraction certifies an MBP chart"
        if outcome == Outcome.PASS
        else "a chart was extracted, but it lies outside the two MBP1 chart classes"
    )
    return certificate, PredicateResult(
        "CS",
        outcome,
        reason,
        _certificate(
            orientations=orientations,
            mechanism=certificate.mechanism,
            mark=certificate.extracted_mark,
            jacobian_unit=certificate.core_jacobian_unit,
        ),
    )


def extract_minimal_boundary(
    record: FiniteNormalizationRecord,
) -> ExtractionReport:
    """Run all eight MBP predicates on one finite normalization export."""

    selected, scb = select_critical_prime(record)
    puncture_rank, rows, pr = compute_punctures(record.critical_curve)
    sat = test_saturation(record.links)
    bm = test_boundary_monotonicity(record.links)
    lc = test_ledger(record.ledger, record.unrecorded_graph_primes)
    pc = test_primitive_conormal(
        record.conormal, puncture_rank, record.critical_curve
    )
    nc = test_noncontraction(
        record.residue_parameter, record.residue_generators
    )
    chart_certificate, cs = attempt_chart_straightening(
        record.chart, record.links
    )
    predicates = (scb, sat, bm, lc, pr, pc, nc, cs)
    return ExtractionReport(
        record_name=record.name,
        selected_prime=selected.label if selected is not None else None,
        puncture_rank=puncture_rank,
        valuation_rows=rows,
        chart_certificate=chart_certificate,
        predicates=predicates,
    )


def mutate_record(
    record: FiniteNormalizationRecord, **changes: object
) -> FiniteNormalizationRecord:
    """Typed convenience wrapper used by countermodel fixtures."""

    return replace(record, **changes)

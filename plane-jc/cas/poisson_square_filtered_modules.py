#!/usr/bin/env python3
"""Search-facing associated-stratum filter for the three-band Poisson box.

The exact primary data are certified by
``poisson_square_normalized_defect.sing``.  This module packages the eight
associated strata as normalized parametrizations and reduces proposed lower
coefficient equations on their dense charts before a larger ideal is built.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Iterable

import sympy as sp

from poisson_square_rigidity import standard_three_band_problem


@dataclass(frozen=True)
class FilteredPoissonStratum:
    """One associated stratum and its certified cyclic normal fibers."""

    name: str
    role: str
    normalized_dimension: int
    parameters: tuple[sp.Symbol, ...]
    nonzero_factors: tuple[sp.Expr, ...]
    substitution_items: tuple[tuple[sp.Symbol, sp.Expr], ...]
    d3_hilbert_vector: tuple[int, ...] | None = None
    d2_hilbert_vector: tuple[int, ...] | None = None
    d3_socle_dimension: int | None = None
    d2_socle_dimension: int | None = None
    d3_presentation: str | None = None
    d2_presentation: str | None = None

    @cached_property
    def substitution(self) -> dict[sp.Symbol, sp.Expr]:
        return dict(self.substitution_items)

    @property
    def d3_length(self) -> int | None:
        if self.d3_hilbert_vector is None:
            return None
        return sum(self.d3_hilbert_vector)

    @property
    def d2_length(self) -> int | None:
        if self.d2_hilbert_vector is None:
            return None
        return sum(self.d2_hilbert_vector)


@dataclass(frozen=True)
class LowerBandStratumAudit:
    """Effect of proposed lower equations on one dense associated chart."""

    stratum: str
    status: str
    restrictions: tuple[sp.Expr, ...]
    d3_hilbert_vector: tuple[int, ...] | None
    d2_hilbert_vector: tuple[int, ...] | None


def _coefficient_substitution(
    A: sp.Expr,
    B: sp.Expr,
    C: sp.Expr,
    D: sp.Expr,
    t: sp.Symbol,
) -> tuple[tuple[sp.Symbol, sp.Expr], ...]:
    items: list[tuple[sp.Symbol, sp.Expr]] = []
    for prefix, polynomial, degree in (
        ("p_3", A, 3),
        ("p_2", B, 4),
        ("q_2", C, 2),
        ("q_1", D, 3),
    ):
        expanded = sp.Poly(sp.expand(polynomial), t)
        items.extend(
            (
                sp.Symbol(f"{prefix}_{index}"),
                sp.factor(expanded.coeff_monomial(t**index)),
            )
            for index in range(degree + 1)
        )
    return tuple(items)


def associated_stratum_ledger() -> tuple[FilteredPoissonStratum, ...]:
    """Return normalized dense charts for all eight associated primes."""

    t = sp.Symbol("t")
    a, e, r = sp.symbols("a e r")
    alpha, delta, k, lam = sp.symbols("alpha delta k lambda")
    b, c = sp.symbols("b c")

    tangent_d = 1 + e * r
    h = t - r
    tangent = (
        a * h**3,
        -h / tangent_d - e * h**2 / (2 * tangent_d**2),
        -sp.Rational(3, 2) * a * tangent_d**2 * h**2,
        1 + e * t,
    )

    D_linear = 1 + delta * t
    B_linear = k * D_linear**2 + 1 / (2 * delta)
    component_c_zero = (alpha * D_linear**3, B_linear, 0, D_linear)
    component_a_zero = (0, B_linear, lam * B_linear, D_linear)
    core = (0, B_linear, 0, D_linear)
    constant_d = (alpha, b - t, 0, 1)
    constant_bc = (0, 1 / (2 * delta), c, D_linear)
    core_constant_d = (0, b - t, 0, 1)
    core_constant_bc = (0, 1 / (2 * delta), 0, D_linear)

    def record(
        name: str,
        role: str,
        dimension: int,
        parameters: tuple[sp.Symbol, ...],
        nonzero: tuple[sp.Expr, ...],
        bands: tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr],
        d3: tuple[int, ...] | None = None,
        d2: tuple[int, ...] | None = None,
        d3_socle: int | None = None,
        d2_socle: int | None = None,
        d3_presentation: str | None = None,
        d2_presentation: str | None = None,
    ) -> FilteredPoissonStratum:
        return FilteredPoissonStratum(
            name=name,
            role=role,
            normalized_dimension=dimension,
            parameters=parameters,
            nonzero_factors=nonzero,
            substitution_items=_coefficient_substitution(*bands, t),
            d3_hilbert_vector=d3,
            d2_hilbert_vector=d2,
            d3_socle_dimension=d3_socle,
            d2_socle_dimension=d2_socle,
            d3_presentation=d3_presentation,
            d2_presentation=d2_presentation,
        )

    ledger = (
        record(
            "T",
            "minimal",
            3,
            (a, e, r),
            (a, tangent_d),
            tangent,
            d2=(1,),
            d2_socle=1,
            d2_presentation="reduced_field",
        ),
        record(
            "C0",
            "minimal",
            3,
            (alpha, delta, k),
            (alpha, delta),
            component_c_zero,
        ),
        record(
            "A0",
            "minimal",
            3,
            (delta, k, lam),
            (delta, lam),
            component_a_zero,
        ),
        record(
            "S",
            "embedded_surface",
            2,
            (delta, k),
            (delta, k),
            core,
            d3=(1, 2),
            d2=(1, 3, 3, 1),
            d3_socle=2,
            d2_socle=2,
            d3_presentation="square_zero_two_generators",
            d2_presentation="three_quadrics_plus_one_cubic",
        ),
        record(
            "S_C",
            "embedded_surface",
            2,
            (alpha, b),
            (alpha,),
            constant_d,
            d3=(1, 2, 2),
            d2=(1, 2, 1),
            d3_socle=2,
            d2_socle=1,
            d3_presentation="two_generators_three_relations",
            d2_presentation="four_relations_on_three_displayed_generators",
        ),
        record(
            "S_A",
            "embedded_surface",
            2,
            (delta, c),
            (delta, c),
            constant_bc,
            d3=(1, 2, 3, 1),
            d2=(1, 1),
            d3_socle=2,
            d2_socle=1,
            d3_presentation="two_generators_four_groebner_relations",
            d2_presentation="dual_numbers",
        ),
        record(
            "K_C",
            "embedded_curve",
            1,
            (b,),
            (),
            core_constant_d,
            d3=(1, 3, 5, 7, 6, 3),
            d3_socle=4,
            d3_presentation="three_generators_nine_groebner_relations",
        ),
        record(
            "K_A",
            "embedded_curve",
            1,
            (delta,),
            (delta,),
            core_constant_bc,
            d3=(1, 4, 8, 10, 6, 1),
            d3_socle=6,
            d3_presentation="four_generators_eighteen_groebner_relations",
        ),
    )

    problem = standard_three_band_problem()
    for stratum in ledger:
        if stratum.substitution[sp.Symbol("q_1_0")] != 1:
            raise AssertionError(f"{stratum.name} is not normalized by d0=1")
        if any(
            sp.cancel(equation.subs(stratum.substitution)) != 0
            for equation in problem.coefficient_equations
        ):
            raise AssertionError(f"{stratum.name} left the Poisson-square locus")
    return ledger


def _localized_unit_ideal(
    restrictions: tuple[sp.Expr, ...],
    stratum: FilteredPoissonStratum,
) -> bool:
    nonzero_restrictions = tuple(
        sp.factor(sp.together(expression).as_numer_denom()[0])
        for expression in restrictions
        if expression != 0
    )
    if not nonzero_restrictions:
        return False
    z = sp.Symbol(f"z_{stratum.name}")
    nonzero_product = sp.prod(stratum.nonzero_factors)
    generators = list(nonzero_restrictions)
    if nonzero_product != 1:
        generators.append(sp.expand(z * nonzero_product - 1))
    variables = (z, *stratum.parameters)
    basis = sp.groebner(generators, *variables, order="grevlex")
    return any(polynomial.as_expr() == 1 for polynomial in basis.polys)


def lower_band_survival_audit(
    equations: Iterable[sp.Expr],
) -> tuple[LowerBandStratumAudit, ...]:
    """Reduce proposed lower equations on all eight associated strata.

    ``preserved`` means every equation vanishes identically on the dense
    chart. ``eliminated`` means the localized restriction ideal is the unit
    ideal. ``cut`` means a proper closed sublocus may survive.
    """

    equations = tuple(sp.expand(equation) for equation in equations)
    reports: list[LowerBandStratumAudit] = []
    for stratum in associated_stratum_ledger():
        restrictions = tuple(
            sp.factor(sp.cancel(equation.subs(stratum.substitution)))
            for equation in equations
        )
        if all(expression == 0 for expression in restrictions):
            status = "preserved"
        elif _localized_unit_ideal(restrictions, stratum):
            status = "eliminated"
        else:
            status = "cut"
        reports.append(
            LowerBandStratumAudit(
                stratum=stratum.name,
                status=status,
                restrictions=restrictions,
                d3_hilbert_vector=stratum.d3_hilbert_vector,
                d2_hilbert_vector=stratum.d2_hilbert_vector,
            )
        )
    return tuple(reports)


def main() -> int:
    ledger = associated_stratum_ledger()
    print(
        "STRATUM|ROLE|DIM|D3_HILBERT|D3_SOCLE|D3_PRESENTATION|"
        "D2_HILBERT|D2_SOCLE|D2_PRESENTATION"
    )
    for stratum in ledger:
        print(
            f"{stratum.name}|{stratum.role}|{stratum.normalized_dimension}|"
            f"{stratum.d3_hilbert_vector}|{stratum.d3_socle_dimension}|"
            f"{stratum.d3_presentation}|{stratum.d2_hilbert_vector}|"
            f"{stratum.d2_socle_dimension}|{stratum.d2_presentation}"
        )
    print("PASS: eight normalized associated strata compile for lower-band use")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

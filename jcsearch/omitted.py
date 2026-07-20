"""Exact omitted-value classification for inverse pencils H(W)-sW+t."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

import sympy as sp


@dataclass(frozen=True)
class OmittedValue:
    """A certified parameter value whose inverse polynomial has no simple root."""

    slope: sp.Expr
    intercept: sp.Expr
    partition: tuple[int, ...]
    roots: tuple[sp.Expr, ...]
    factorization: sp.Expr


def multiplicity_partitions(total: int, maximum: int | None = None) -> Iterator[tuple[int, ...]]:
    """Yield unordered partitions of ``total`` with every part at least two."""
    if total == 0:
        yield ()
        return
    maximum = total if maximum is None else min(maximum, total)
    for first in range(maximum, 1, -1):
        for tail in multiplicity_partitions(total - first, first):
            yield (first,) + tail


def classify_omitted_values(primitive, variable) -> tuple[OmittedValue, ...]:
    """Classify all ``(s,t)`` for which ``H-sW+t`` has no simple root.

    The input must be an exact characteristic-zero polynomial of degree at
    least three.  Completeness follows by enumerating every multiplicity
    partition of its degree with parts at least two.  The coefficients of
    degrees two and higher determine a zero-dimensional factor-matching
    system; the remaining linear and constant coefficients recover ``s,t``.
    """
    H = sp.Poly(sp.sympify(primitive), variable)
    degree = H.degree()
    if degree < 3:
        raise ValueError("omitted values are not finite for inverse degree below three")
    if H.domain.characteristic() != 0:
        raise ValueError("the classifier requires characteristic zero")

    leading = H.LC()
    results: list[OmittedValue] = []
    seen: set[tuple[str, str]] = set()

    for partition in multiplicity_partitions(degree):
        roots = sp.symbols(f"_omitted_r0:{len(partition)}")
        monic = sp.prod(
            (variable - root) ** multiplicity
            for root, multiplicity in zip(roots, partition)
        )
        difference = sp.Poly(sp.expand(H.as_expr() - leading * monic), variable)
        equations = [
            difference.coeff_monomial(variable**power)
            for power in range(2, degree)
        ]
        basis = sp.groebner(equations, *roots, order="lex")
        if len(basis.polys) == 1 and basis.polys[0].as_expr() == 1:
            continue
        if not basis.is_zero_dimensional:
            raise ValueError(
                f"factor-matching system for partition {partition} is not zero-dimensional"
            )
        solutions = sp.solve_poly_system(list(basis.polys), *roots) or []

        for solution in solutions:
            substitutions = dict(zip(roots, solution))
            factorization = sp.factor(leading * monic.subs(substitutions), extension=True)
            residual = sp.Poly(sp.expand(H.as_expr() - factorization), variable, extension=True)
            assert residual.degree() <= 1
            slope = sp.factor(residual.coeff_monomial(variable), extension=True)
            intercept = sp.factor(-residual.coeff_monomial(1), extension=True)
            assert sp.expand(
                H.as_expr() - slope * variable + intercept - factorization
            ) == 0

            key = (sp.srepr(sp.simplify(slope)), sp.srepr(sp.simplify(intercept)))
            if key in seen:
                continue
            seen.add(key)
            results.append(
                OmittedValue(
                    slope=slope,
                    intercept=intercept,
                    partition=partition,
                    roots=tuple(solution),
                    factorization=factorization,
                )
            )

    return tuple(results)

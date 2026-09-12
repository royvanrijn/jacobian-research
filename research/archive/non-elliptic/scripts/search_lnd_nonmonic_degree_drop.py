#!/usr/bin/env python3
"""Exact bounded search on the first unresolved nonmonic slice-LND charts.

The derivation is D=d/ds on Q[x,s].  We use the degree-drop carrier

    p = x*s - 1

and test I=x^c p^d J for (c,d)=(0,1),(0,2),(1,1) and three small primary
residual ideals J at (x,s)=(0,0).  Membership is exact: the unique primitive
vanishing at the rational root s=1/x is divided by x^c p^d and then reduced
modulo J.

The d=1 run is a regression for the support-weight theorem in the
accompanying note when c=0.  The repeated-carrier and invariant-content
runs are candidate generation only; a bounded prefix is never reported as
a proof or counterexample.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import sympy as sp


x, s = sp.symbols("x s")
p = x * s - 1


@dataclass(frozen=True)
class ResidualChart:
    name: str
    generators: tuple[sp.Expr, ...]

    @property
    def basis(self) -> sp.GroebnerBasis:
        return sp.groebner(self.generators, s, x, domain=sp.QQ)


CHARTS = (
    ResidualChart("x-double", (s, x**2)),
    ResidualChart("s-double", (x, s**2)),
    ResidualChart("double-point", (x**2, x * s, s**2)),
)


def polynomial_or_none(expr: sp.Expr) -> sp.Expr | None:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    if denominator.has(x, s):
        return None
    return sp.expand(numerator / denominator)


def carrier_quotient(
    h: sp.Expr, content_power: int, carrier_power: int
) -> sp.Expr | None:
    primitive = sp.integrate(sp.expand(h), s)
    root_value = sp.cancel(primitive.subs(s, 1 / x))
    return polynomial_or_none(
        (primitive - root_value) / (x**content_power * p**carrier_power)
    )


def in_residual(expr: sp.Expr, chart: ResidualChart) -> bool:
    remainder = chart.basis.reduce(sp.Poly(expr, s, x, domain=sp.QQ).as_expr())[1]
    return sp.expand(remainder) == 0


def in_image(
    h: sp.Expr,
    content_power: int,
    carrier_power: int,
    chart: ResidualChart,
) -> bool:
    quotient = carrier_quotient(h, content_power, carrier_power)
    return quotient is not None and in_residual(quotient, chart)


def support_weight(f: sp.Expr) -> int:
    poly = sp.Poly(sp.expand(f), x, s, domain=sp.QQ)
    return min(i - j for (i, j), coefficient in poly.terms() if coefficient)


def seed_family() -> list[sp.Expr]:
    monomials = [
        x**i * s**j
        for i in range(5)
        for j in range(4)
        if i + j > 0 and i + j <= 5
    ]
    seeds = list(monomials)
    for left, right in combinations(monomials, 2):
        seeds.append(left + right)
        seeds.append(left - right)
    return seeds


def multiplier_family() -> tuple[sp.Expr, ...]:
    return (sp.Integer(1), x, s, x**2, x * s, s**2)


def main() -> None:
    pure_bound = 6
    mixed_start = 4
    seeds = seed_family()
    multipliers = multiplier_family()

    print(f"SEARCH: {len(seeds)} monomial/binomial seeds; exact powers 1..{pure_bound}")
    total_survivors = 0
    total_tail_obstructions = 0

    carrier_profiles = ((0, 1), (0, 2), (1, 1))
    for content_power, carrier_power in carrier_profiles:
        for chart in CHARTS:
            survivors: list[sp.Expr] = []
            tail_obstructions: list[tuple[sp.Expr, sp.Expr]] = []

            for f in seeds:
                if all(
                    in_image(
                        f**m,
                        content_power,
                        carrier_power,
                        chart,
                    )
                    for m in range(1, pure_bound + 1)
                ):
                    survivors.append(f)
                    if content_power == 0 and carrier_power == 1:
                        assert support_weight(f) >= 1

                    for g in multipliers:
                        failures = [
                            not in_image(
                                g * f**m,
                                content_power,
                                carrier_power,
                                chart,
                            )
                            for m in range(mixed_start, pure_bound + 1)
                        ]
                        if all(failures):
                            tail_obstructions.append((f, g))

            total_survivors += len(survivors)
            total_tail_obstructions += len(tail_obstructions)
            print(
                "CHART:",
                f"x^{content_power}",
                f"p^{carrier_power}",
                chart.name,
                f"pure-prefix survivors={len(survivors)}",
                f"bounded-tail obstructions={len(tail_obstructions)}",
            )

            if carrier_power == 2 or content_power > 0:
                for f in survivors:
                    print("  SURVIVOR:", f"f={f}", f"weight={support_weight(f)}")

            for f, g in tail_obstructions[:3]:
                print("  CANDIDATE-ONLY:", f"f={f}", f"g={g}")

    print(f"SUMMARY: pure-prefix survivors={total_survivors}")
    print(f"SUMMARY: bounded-tail obstructions={total_tail_obstructions}")
    print("NOTE: bounded survival or failure is not an all-order certificate")


if __name__ == "__main__":
    main()

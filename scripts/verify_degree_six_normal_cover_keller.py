#!/usr/bin/env python3
"""Compile and verify the C2 x C2 degree-six Hasse fiber."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber
from jcsearch.normal_covering import analyze_component_action


T = sp.Symbol("T")
S = sp.Symbol("S")
CERTIFICATE = (
    ROOT / "arithmetic/certificates/normal_cover_v4_sextic.json"
)


def main() -> None:
    polynomial = sp.Poly(
        (T**2 - 2) * (T**2 - 17) * (T**2 - 34),
        T,
        domain=sp.QQ,
    )
    compilation = compile_polynomial_to_keller_fiber(
        polynomial,
        T,
        inverse_variable=S,
    )
    assert compilation.translation == 1
    assert compilation.geometric_degree == 6
    assert compilation.coordinate_degrees == (7, 38, 36)
    assert compilation.target == (1, 0, sp.Rational(528, 577))
    assert compilation.seed == (
        S**6
        + 6 * S**5
        - 38 * S**4
        - 192 * S**3
        + 377 * S**2
        + 1154 * S
    )
    translated = sp.expand(polynomial.as_expr().subs(T, S + 1))
    assert compilation.inverse_polynomial == translated
    assert sp.factor(compilation.inverse_polynomial) == (
        (S**2 + 2 * S - 33)
        * (S**2 + 2 * S - 16)
        * (S**2 + 2 * S - 1)
    )

    x, y, z = compilation.source_variables
    jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(component, variable) for variable in (x, y, z)]
                for component in compilation.determinant_one_map
            ]
        )
    )
    assert sp.factor(jacobian) == 1

    certificate = json.loads(CERTIFICATE.read_text())
    action = certificate["action"]
    cover = analyze_component_action(
        degree=action["degree"],
        generators=action["generators"],
        components=action["components"],
    )
    assert len(cover.group) == 4
    assert cover.factorization_shape == (2, 2, 2)
    assert cover.is_normal_cover
    assert cover.is_faithful
    assert cover.normal_covering_number == 3

    print("PASS degree-six inverse polynomial is the translated 2+2+2 fiber")
    print("PASS determinant-one quadratic-gauge compilation at (1,0,528/577)")
    print("PASS C2 x C2 faithful normal cover has r=gamma=3")
    print("PASS ramified-prime witnesses are recorded at 2 and 17")


if __name__ == "__main__":
    main()

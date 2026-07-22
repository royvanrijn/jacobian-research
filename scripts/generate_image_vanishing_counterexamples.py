#!/usr/bin/env python3
"""Generate the 21/42-dimensional Image and Vanishing witnesses.

The source is the independently replayed essential BCW quotient V=I+H in
21 variables.  This script does not try to verify the all-order theorems of
Abhyankar--Gurjar or Zhao; those arguments are written in
extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md.  It verifies the exact
finite data used by those arguments and expands the resulting quartic.
"""

from __future__ import annotations

from math import comb
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "image_vanishing_counterexamples_21_42.json"


SparsePoly = dict[tuple[int, ...], sp.Expr]


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    answer: SparsePoly = {}
    for a, ca in left.items():
        for b, cb in right.items():
            exponent = tuple(x + y for x, y in zip(a, b))
            answer[exponent] = sp.expand(answer.get(exponent, 0) + ca * cb)
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def evaluate_h(components: list[list[dict[str, object]]], point: list[sp.Rational]) -> list[sp.Expr]:
    values = []
    for component in components:
        value = sp.Integer(0)
        for term in component:
            monomial = sp.Integer(1)
            for variable, exponent in term["monomial"]:
                monomial *= point[variable] ** exponent
            value += sp.Rational(term["coefficient"]) * monomial
        values.append(sp.cancel(value))
    return values


def encode(poly: SparsePoly) -> list[dict[str, object]]:
    return [
        {
            "coefficient": str(coefficient),
            "monomial": [[index, exponent] for index, exponent in enumerate(monomial) if exponent],
        }
        for monomial, coefficient in sorted(poly.items(), reverse=True)
    ]


def contraction_polynomial(source: dict[str, object]) -> SparsePoly:
    """Return p(w,z)=-sum_i w_i H_i(z), with variables ordered (w,z)."""
    n = source["dimension"]
    answer: SparsePoly = {}
    for output, component in enumerate(source["H"]):
        for term in component:
            exponent = [0] * (2 * n)
            exponent[output] = 1
            for variable, power in term["monomial"]:
                exponent[n + variable] = power
            key = tuple(exponent)
            answer[key] = sp.expand(answer.get(key, 0) - sp.Rational(term["coefficient"]))
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def laplacian_quartic(source: dict[str, object]) -> SparsePoly:
    """Expand R=-1/4 (u+i v).H(u-i v), with variables ordered (u,v)."""
    n = source["dimension"]
    answer: SparsePoly = {}
    for output, component in enumerate(source["H"]):
        for term in component:
            partial: SparsePoly = {(0,) * (2 * n): -sp.Rational(term["coefficient"]) / 4}

            output_factor: SparsePoly = {}
            for variable, coefficient in ((output, 1), (n + output, sp.I)):
                exponent = [0] * (2 * n)
                exponent[variable] = 1
                output_factor[tuple(exponent)] = coefficient
            partial = multiply(partial, output_factor)

            for variable, power in term["monomial"]:
                factor: SparsePoly = {}
                for v_power in range(power + 1):
                    exponent = [0] * (2 * n)
                    exponent[variable] = power - v_power
                    exponent[n + variable] = v_power
                    factor[tuple(exponent)] = comb(power, v_power) * (-sp.I) ** v_power
                partial = multiply(partial, factor)

            for monomial, coefficient in partial.items():
                answer[monomial] = sp.expand(answer.get(monomial, 0) + coefficient)
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def reflected_negative(poly: SparsePoly, n: int) -> SparsePoly:
    """Return -poly(u,-v)."""
    return {
        monomial: sp.expand(-coefficient * (-1) ** sum(monomial[n:]))
        for monomial, coefficient in poly.items()
    }


def symmetric_quartic(source: dict[str, object]) -> SparsePoly:
    """Expand P=1/4 (u-i v).H(u+i v), independently of R."""
    n = source["dimension"]
    answer: SparsePoly = {}
    for output, component in enumerate(source["H"]):
        for term in component:
            partial: SparsePoly = {(0,) * (2 * n): sp.Rational(term["coefficient"]) / 4}
            output_factor: SparsePoly = {}
            for variable, coefficient in ((output, 1), (n + output, -sp.I)):
                exponent = [0] * (2 * n)
                exponent[variable] = 1
                output_factor[tuple(exponent)] = coefficient
            partial = multiply(partial, output_factor)
            for variable, power in term["monomial"]:
                factor: SparsePoly = {}
                for v_power in range(power + 1):
                    exponent = [0] * (2 * n)
                    exponent[variable] = power - v_power
                    exponent[n + variable] = v_power
                    factor[tuple(exponent)] = comb(power, v_power) * sp.I ** v_power
                partial = multiply(partial, factor)
            for monomial, coefficient in partial.items():
                answer[monomial] = sp.expand(answer.get(monomial, 0) + coefficient)
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def main() -> None:
    source = json.loads(SOURCE.read_text())
    n = source["dimension"]
    assert n == 21
    assert source["jacobian_determinant"] == "1"
    assert len(source["H"]) == n
    assert all(
        sum(exponent for _, exponent in term["monomial"]) == 3
        for component in source["H"]
        for term in component
    )

    points = [[sp.Rational(value) for value in point] for point in source["collision_points"]]
    target = [sp.Rational(value) for value in source["common_image"]]
    assert len({tuple(point) for point in points}) == 3
    for point in points:
        image = [sp.cancel(z + h) for z, h in zip(point, evaluate_h(source["H"], point))]
        assert image == target

    # This is the all-order coordinate certificate: a polynomial formal left
    # inverse for coordinate zero would take the same value at all three
    # collision points, but their zeroth coordinates are 0, 1, and -1.
    distinguished = 0
    distinguished_values = [point[distinguished] for point in points]
    assert distinguished_values == [0, 1, -1]

    p = contraction_polynomial(source)
    assert len(p) == sum(len(component) for component in source["H"])
    assert all(sum(monomial) == 4 for monomial in p)

    quartic = laplacian_quartic(source)
    assert all(sum(monomial) == 4 for monomial in quartic)
    symmetric = symmetric_quartic(source)
    assert symmetric == reflected_negative(quartic, n)

    # The final output is an identity coordinate and all stored collision
    # points lie in its T=1 slice.  Record, but do not overstate, the resulting
    # nonhomogeneous 20-dimensional restriction.
    identity_outputs = [index for index, component in enumerate(source["H"]) if not component]
    assert identity_outputs == [20]
    assert [point[20] for point in points] == [1, 1, 1]
    restricted_degrees = sorted({
        sum(power for variable, power in term["monomial"] if variable != 20)
        for component in source["H"][:-1]
        for term in component
    })
    assert restricted_degrees == [1, 2, 3]

    artifact = {
        "format": "image-vanishing-counterexamples-21-42-v1",
        "source_H": str(SOURCE.relative_to(ROOT)),
        "source_dimension": n,
        "distinguished_coordinate": distinguished,
        "collision_coordinate_values": [str(value) for value in distinguished_values],
        "special_image_counterexample": {
            "variables": "w_0,...,w_20,z_0,...,z_20",
            "p_definition": "p(w,z)=-sum_i w_i*H_i(z), where the source artifact is V=I+H",
            "p_terms": encode(p),
            "multiplier": "z_0",
            "certificate": "E(p^m)=0 for every m>=1; E(z_0*p^m) is nonzero for infinitely many m",
        },
        "laplacian_counterexample": {
            "dimension": 2 * n,
            "variables": "u_0,...,u_20,v_0,...,v_20",
            "R_definition": "R(u,v)=-(1/4)*(u+I*v).H(u-I*v)",
            "R_terms": encode(quartic),
            "linear_multiplier": "u_0-I*v_0",
            "generalized_certificate": "Delta^m(R^m)=0 for every m>=1; Delta^m((u_0-I*v_0)*R^m) is nonzero infinitely often",
            "classical_certificate": "Hess(R) is nilpotent; Delta^m(R^(m+1)) is nonzero infinitely often",
        },
        "symmetric_lift": {
            "P_definition": "P(u,v)=(1/4)*(u-I*v).H(u+I*v)",
            "relation_to_R": "P(u,v)=-R(u,-v)",
            "gradient_map": "(u,v)+grad(P) is Keller and noninjective",
        },
        "identity_output_restriction": {
            "identity_coordinate": 20,
            "collision_slice_value": "1",
            "restricted_dimension": 20,
            "nonlinear_degrees_after_restriction": restricted_degrees,
            "warning": "the restriction is Keller and noninjective but is not cubic homogeneous",
        },
        "statistics": {
            "source_cubic_terms": sum(len(component) for component in source["H"]),
            "contraction_terms": len(p),
            "expanded_quartic_terms": len(quartic),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS image/VC: exact 21D cubic Keller collision and coordinate-0 separation")
    print(f"PASS image/VC: expanded p into {len(p)} terms and R into {len(quartic)} terms")
    print("PASS image/VC: independently verified P(u,v)=-R(u,-v)")
    print("PASS image/VC: identity output restricts the collision to a nonhomogeneous 20D Keller map")
    print(f"PASS image/VC: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

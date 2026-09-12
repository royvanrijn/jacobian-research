#!/usr/bin/env python3
"""Generate the 20/40-dimensional witnesses from the identity-output slice."""

from __future__ import annotations

from math import comb
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "image_vanishing_counterexamples_20_40.json"

SparsePoly = dict[tuple[int, ...], sp.Expr]


def add_term(poly: SparsePoly, exponent: tuple[int, ...], coefficient: sp.Expr) -> None:
    value = sp.expand(poly.get(exponent, 0) + coefficient)
    if value:
        poly[exponent] = value
    elif exponent in poly:
        del poly[exponent]


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    answer: SparsePoly = {}
    for a, ca in left.items():
        for b, cb in right.items():
            add_term(answer, tuple(x + y for x, y in zip(a, b)), ca * cb)
    return answer


def restrict_components(source: dict[str, object]) -> list[SparsePoly]:
    """Set the final source variable s to one and drop its identity output."""
    n = source["dimension"] - 1
    components: list[SparsePoly] = []
    for stored_component in source["H"][:n]:
        component: SparsePoly = {}
        for term in stored_component:
            exponent = [0] * n
            for variable, power in term["monomial"]:
                if variable < n:
                    exponent[variable] = power
                else:
                    assert variable == n
            add_term(component, tuple(exponent), sp.Rational(term["coefficient"]))
        components.append(component)
    return components


def evaluate(poly: SparsePoly, point: list[sp.Rational]) -> sp.Expr:
    return sp.expand(sum(
        coefficient * sp.prod(value**power for value, power in zip(point, exponent))
        for exponent, coefficient in poly.items()
    ))


def contraction_polynomial(components: list[SparsePoly]) -> SparsePoly:
    n = len(components)
    answer: SparsePoly = {}
    for output, component in enumerate(components):
        for old, coefficient in component.items():
            exponent = [0] * (2 * n)
            exponent[output] = 1
            exponent[n:] = old
            add_term(answer, tuple(exponent), -coefficient)
    return answer


def transformed_polynomial(components: list[SparsePoly], symmetric: bool) -> SparsePoly:
    """Expand -(u+iv).h(u-iv), or independently +(u-iv).h(u+iv)."""
    n = len(components)
    answer: SparsePoly = {}
    for output, component in enumerate(components):
        for old, coefficient in component.items():
            initial = coefficient if symmetric else -coefficient
            partial: SparsePoly = {(0,) * (2 * n): initial}
            output_factor: SparsePoly = {}
            output_sign = -sp.I if symmetric else sp.I
            for variable, factor_coefficient in ((output, 1), (n + output, output_sign)):
                exponent = [0] * (2 * n)
                exponent[variable] = 1
                output_factor[tuple(exponent)] = factor_coefficient
            partial = multiply(partial, output_factor)

            input_sign = sp.I if symmetric else -sp.I
            for variable, power in enumerate(old):
                if not power:
                    continue
                factor: SparsePoly = {}
                for v_power in range(power + 1):
                    exponent = [0] * (2 * n)
                    exponent[variable] = power - v_power
                    exponent[n + variable] = v_power
                    factor[tuple(exponent)] = comb(power, v_power) * input_sign**v_power
                partial = multiply(partial, factor)
            for exponent, partial_coefficient in partial.items():
                add_term(answer, exponent, partial_coefficient)
    return answer


def reflected_negative(poly: SparsePoly, n: int) -> SparsePoly:
    return {
        exponent: sp.expand(-coefficient * (-1) ** sum(exponent[n:]))
        for exponent, coefficient in poly.items()
    }


def encode_rational(poly: SparsePoly) -> list[dict[str, object]]:
    return [
        {
            "coefficient": str(sp.Rational(coefficient)),
            "monomial": [[index, power] for index, power in enumerate(exponent) if power],
        }
        for exponent, coefficient in sorted(poly.items(), reverse=True)
    ]


def encode_complex(poly: SparsePoly) -> list[dict[str, object]]:
    encoded = []
    for exponent, coefficient in sorted(poly.items(), reverse=True):
        real, imaginary = sp.expand_complex(coefficient).as_real_imag()
        encoded.append({
            "coefficient": {"real": str(sp.Rational(real)), "imaginary": str(sp.Rational(imaginary))},
            "monomial": [[index, power] for index, power in enumerate(exponent) if power],
        })
    return encoded


def main() -> None:
    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == 21
    assert source["H"][20] == []
    points = [[sp.Rational(value) for value in point] for point in source["collision_points"]]
    assert [point[20] for point in points] == [1, 1, 1]

    components = restrict_components(source)
    n = len(components)
    assert n == 20 and all(components)
    degrees = sorted({sum(exponent) for component in components for exponent in component})
    assert degrees == [1, 2, 3]
    monomials = sorted({exponent for component in components for exponent in component})
    output_coefficient_matrix = sp.Matrix([
        [component.get(exponent, 0) for exponent in monomials]
        for component in components
    ])
    assert output_coefficient_matrix.rank() == n

    restricted_points = [point[:n] for point in points]
    restricted_images = [
        [sp.expand(point[i] + evaluate(components[i], point)) for i in range(n)]
        for point in restricted_points
    ]
    assert restricted_images[0] == restricted_images[1] == restricted_images[2]
    assert [point[0] for point in restricted_points] == [0, 1, -1]

    p = contraction_polynomial(components)
    assert all(2 <= sum(exponent) <= 4 for exponent in p)
    laplacian = transformed_polynomial(components, symmetric=False)
    symmetric = transformed_polynomial(components, symmetric=True)
    assert symmetric == reflected_negative(laplacian, n)
    assert sorted({sum(exponent) for exponent in laplacian}) == [2, 3, 4]

    artifact = {
        "format": "identity-slice-image-vanishing-counterexamples-20-40-v1",
        "source_H": str(SOURCE.relative_to(ROOT)),
        "slice": {"identity_coordinate": 20, "value": "1"},
        "source_dimension": n,
        "collision_points": [[str(value) for value in point] for point in restricted_points],
        "common_image": [str(value) for value in restricted_images[0]],
        "distinguished_coordinate": 0,
        "collision_coordinate_values": ["0", "1", "-1"],
        "restricted_nonlinear_part": {
            "definition": "h(q)=source_H(q,1), dropping the final zero output",
            "degrees": degrees,
            "output_span_rank": n,
            "linear_identity_output_functionals": 0,
            "components": [encode_rational(component) for component in components],
        },
        "special_image_counterexample": {
            "variables": "w_0,...,w_19,q_0,...,q_19",
            "p_definition": "p_20(w,q)=-sum_i w_i*h_i(q)",
            "p_terms": encode_rational(p),
            "multiplier": "q_0",
            "certificate": "specialization at s=1 of the homogeneous 21D contractions",
        },
        "laplacian_counterexample": {
            "dimension": 2 * n,
            "variables": "u_0,...,u_19,v_0,...,v_19",
            "R_definition": "R_20(u,v)=-(u+I*v).h(u-I*v)",
            "R_terms": encode_complex(laplacian),
            "degrees": [2, 3, 4],
            "linear_multiplier": "u_0-I*v_0",
            "operator_change": "sum d_w_i d_q_i = (1/4)*Delta_40",
            "generalized_certificate": "Delta^m(R_20^m)=0 always; Delta^m((u_0-I*v_0)*R_20^m) is nonzero infinitely often",
            "nonhomogeneous_hn_certificate": "Hess(R_20) is nilpotent and Delta^m(R_20^(m+1)) is nonzero infinitely often",
        },
        "symmetric_lift": {
            "P_definition": "P_20(u,v)=(u-I*v).h(u+I*v)",
            "relation_to_R": "P_20(u,v)=-R_20(u,-v)",
        },
        "statistics": {
            "restricted_h_terms": sum(len(component) for component in components),
            "contraction_terms": len(p),
            "expanded_laplacian_terms": len(laplacian),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS identity slice: exact noninjective 20D Keller restriction at s=1")
    print("PASS identity slice: no further linear identity-output functional")
    print(f"PASS identity slice: p_20 has {len(p)} terms and R_20 has {len(laplacian)} terms")
    print("PASS identity slice: independently expanded P_20(u,v)=-R_20(u,-v)")
    print(f"PASS identity slice: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

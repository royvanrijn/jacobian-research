#!/usr/bin/env sage
"""Audit the compact CM24 fourth-q12 moving cubic over GF(73).

The input artifact is emitted by
``derive_q80_fourth_q12_local_gates_gf73.sage`` after Smith saturation and
the exact component-valuation solve.  This small checker does not rebuild the
local arcs.  It reconstructs the bidegree-(14,3) moving equation and factors
every rational reducible fiber, retaining vertical factors in the old base
that disappear over ``GF(73)(V)``.

With ``--discriminant`` it also removes the persistent square factors from
the cubic-in-X discriminant.  The remaining degree-six branch factor assigns
the finite semistable orders; the exact ``A5+A3+4A2`` child lattice assigns
the final ``I3`` to infinity.  ``--write-artifact`` pins this compact derived
certificate without modifying the parent moving-cubic artifact.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-moving-cubic-gf73.json"
)
artifact_bytes = ARTIFACT.read_bytes()
artifact_hash = hashlib.sha256(artifact_bytes).hexdigest()
assert artifact_hash == (
    "c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89"
)
data = json.loads(artifact_bytes)
assert data["schema"] == "q80-fourth-q12-cm24-moving-cubic-gf73-v1"
assert data["prime"] == 73

finite = GF(73, impl="modn")
parser = argparse.ArgumentParser()
parser.add_argument("--sections", action="store_true")
parser.add_argument("--discriminant", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()
if arguments.write_artifact:
    arguments.discriminant = True
ring = PolynomialRing(finite, names=("T", "V", "X"))
T, V, X = ring.gens()
moving = sum(
    finite(coefficient)*T**t_degree*V**v_degree*X**x_degree
    for t_degree, v_degree, x_degree, coefficient
    in data["moving_terms_T_v_x_coefficient"]
)
assert (
    moving.degree(T), moving.degree(V), moving.degree(X)
) == tuple(data["moving_degrees_T_v_x"])

support = tuple(item[0] for item in data["rational_reducible_support"])
factor_rows = []
for parameter in support:
    specialized = ring(moving(T=finite(parameter)))
    factorization = tuple(specialized.factor())
    row = (
        int(parameter),
        tuple(
            (
                int(factor.degree(V)),
                int(factor.degree(X)),
                int(exponent),
                str(factor),
            )
            for factor, exponent in factorization
        ),
    )
    factor_rows.append(row)

infinity = ring(moving.coefficient({T: moving.degree(T)}))
infinity_factorization = tuple(infinity.factor())
infinity_row = tuple(
    (
        int(factor.degree(V)),
        int(factor.degree(X)),
        int(exponent),
        str(factor),
    )
    for factor, exponent in infinity_factorization
)

generic = ring(moving(T=finite(1)))
assert len(tuple(generic.factor())) == 1

print(
    "Q80FOURTHMOVINGGF73|"
    f"artifact={ARTIFACT}|sha256={artifact_hash}|degrees_T,V,X=2,14,3|"
    "generic_T=1_irreducible=1",
    flush=True,
)
for parameter, factors in factor_rows:
    print(
        f"Q80FOURTHMOVINGGF73|T={parameter}|factors={factors}|"
        "status=PASS_SPECIAL_FACTOR",
        flush=True,
    )
print(
    f"Q80FOURTHMOVINGGF73|T=infinity|factors={infinity_row}|"
    "status=PASS_INFINITY_FACTOR",
    flush=True,
)

if arguments.discriminant:
    parameter_field = FunctionField(finite, "tau")
    tau = parameter_field.gen()
    branch_ring = PolynomialRing(parameter_field, "v")
    branch_v = branch_ring.gen()
    cubic_ring = PolynomialRing(branch_ring, "x")
    branch_x = cubic_ring.gen()
    cubic_equation = cubic_ring(sum(
        finite(coefficient)
        *tau**t_degree*branch_v**v_degree*branch_x**x_degree
        for t_degree, v_degree, x_degree, coefficient
        in data["moving_terms_T_v_x_coefficient"]
    ))
    cubic_discriminant_factorization = tuple(
        cubic_equation.discriminant().factor()
    )
    branch_factors = tuple(
        factor for factor, exponent in cubic_discriminant_factorization
        if factor.degree() == 6 and exponent == 1
    )
    assert len(branch_factors) == 1
    branch_sextic = branch_factors[0]
    branch_discriminant = branch_sextic.discriminant()
    branch_numerator = branch_discriminant.numerator()
    branch_denominator = branch_discriminant.denominator()

    def order_at(polynomial, value):
        divisor = polynomial.parent().gen()-finite(value)
        order = 0
        while polynomial % divisor == 0:
            polynomial //= divisor
            order += 1
        return order

    finite_orders = {
        parameter: order_at(branch_numerator, parameter)
        for parameter in support
    }
    assert finite_orders == {14: 3, 25: 4, 47: 6, 58: 3, 67: 3}
    finite_fibers = {
        parameter: f"I{order}" for parameter, order in finite_orders.items()
    }
    print(
        "Q80FOURTHMOVINGGF73|"
        f"branch_sextic_discriminant_numerator={branch_numerator.factor()}|"
        f"denominator={branch_denominator.factor()}|"
        f"finite_orders={finite_orders}|finite_fibers={finite_fibers}|"
        "infinity=I3(lattice_remaining_A2)|"
        "signature=I6+I4+4I3|"
        "status=PASS_BRANCH_DISCRIMINANT_SIGNATURE",
        flush=True,
    )
    if arguments.write_artifact:
        output = {
            "schema": "q80-fourth-q12-cm24-discriminant-gf73-v1",
            "prime": int(73),
            "parent_artifact": str(ARTIFACT.relative_to(ROOT)),
            "parent_sha256": artifact_hash,
            "moving_degrees_T_v_x": list(map(int, (2, 14, 3))),
            "finite_fibers": [
                {
                    "T": int(parameter),
                    "branch_discriminant_order": int(finite_orders[parameter]),
                    "fiber": finite_fibers[parameter],
                    "factor_degrees_v_x_exponent": [
                        list(map(int, factor[:3]))
                        for row_parameter, factors in factor_rows
                        if row_parameter == parameter
                        for factor in factors
                    ],
                }
                for parameter in support
            ],
            "infinity": {
                "factor_degrees_v_x_exponent": [
                    list(map(int, factor[:3])) for factor in infinity_row
                ],
                "fiber": "I3",
                "assignment": (
                    "remaining A2 in exact specialized child lattice "
                    "A5+A3+4A2"
                ),
            },
            "semistable_signature": "I6+I4+4I3",
            "stripped_branch_sextic_coefficients_low_to_high": [
                {
                    "numerator_coefficients_low_to_high": list(
                        map(int, coefficient.numerator().list())
                    ),
                    "denominator_coefficients_low_to_high": list(
                        map(int, coefficient.denominator().list())
                    ),
                }
                for coefficient in branch_sextic.list()
            ],
            "branch_discriminant_numerator_factorization": str(
                branch_numerator.factor()
            ),
            "branch_discriminant_denominator_factorization": str(
                branch_denominator.factor()
            ),
            "reproduce": (
                "sage elkies-k3/scripts/"
                "analyze_q80_fourth_q12_moving_cubic_gf73.sage "
                "--discriminant --write-artifact"
            ),
        }
        output_path = (
            ROOT / "artifacts/generated-results/"
            "q80-fourth-q12-cm24-discriminant-gf73.json"
        )
        encoded = json.dumps(output, indent=2, sort_keys=True)+"\n"
        output_path.write_text(encoded)
        output_hash = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
        print(
            "Q80FOURTHMOVINGGF73|"
            f"artifact={output_path}|sha256={output_hash}|"
            "status=PASS_DISCRIMINANT_ARTIFACT_WRITE",
            flush=True,
        )

if arguments.sections:
    search_script = (
        ROOT
        / "elkies-k3/scripts/"
        "search_q80_third_child_polynomial_sections_gf73.sage"
    )
    saved_arguments = sys.argv
    sys.argv = [str(search_script), "--two-node"]
    load(str(search_script))
    sys.argv = saved_arguments
    assert len(two_node_candidates) == 30

    section_rows = []
    seen_x = set()
    section_ring = PolynomialRing(finite, names=("Ts", "Vs"))
    Ts, Vs = section_ring.gens()
    for section_x, section_y in two_node_candidates:
        x_key = tuple(map(int, section_x.list()))
        if x_key in seen_x:
            continue
        seen_x.add(x_key)
        restriction = section_ring(
            moving(T=Ts, V=Vs, X=section_x(Vs))
        )
        restriction_factorization = tuple(restriction.factor())
        degree_data = tuple(
            (
                int(factor.degree(Ts)),
                int(factor.degree(Vs)),
                int(exponent),
            )
            for factor, exponent in restriction_factorization
        )
        degree_one_maps = []
        for factor, exponent in restriction_factorization:
            if factor.degree(Ts) != 1 or factor.degree(Vs) > 1:
                continue
            coefficient_zero = section_ring(factor(Ts=0))
            coefficient_one = section_ring(factor(Ts=1))-coefficient_zero
            degree_one_maps.append(
                str(-coefficient_zero/coefficient_one)
            )
        if degree_one_maps:
            section_rows.append(
                (
                    x_key,
                    tuple(map(int, section_y.list())),
                    degree_data,
                    tuple(degree_one_maps),
                )
            )
    print(
        "Q80FOURTHMOVINGGF73|"
        f"old_polynomial_x_classes={len(seen_x)}|"
        f"section_or_vertical_hits={len(section_rows)}|"
        f"nonconstant_degree_one_sections={sum(any('/' in value for value in row[3]) for row in section_rows)}|"
        f"degree_one_section_hits={tuple(section_rows)}|"
        "status=PASS_OLD_SECTION_RESTRICTION_SCAN",
        flush=True,
    )
    assert len(section_rows) == 9
    assert sum(
        any("/" in value for value in row[3]) for row in section_rows
    ) == 5
print("Q80FOURTHMOVINGGF73|status=PASS", flush=True)

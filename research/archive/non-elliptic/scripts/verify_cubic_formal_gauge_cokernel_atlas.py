#!/usr/bin/env python3
"""Exact formal-gauge cokernel atlas for all ternary-cubic symbols.

The smooth universal quartic theorem uses the determinant-twisted gauge
differential G and the fact that ker(C)/im(G) is concentrated in collision
degree three.  This checker computes the same graded quotient for every
orbit representative in the existing cubic-symbol atlas.  It records the
precise boundary of that rigidity argument; it does not infer failure of
cotangent saturation from the existence of nongauge tensor corrections.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import research_universal_cubic_quartic_kernel_saturation as frontier  # noqa: E402
import verify_cubic_symbol_double_saturation as cubic  # noqa: E402
import verify_universal_cubic_cotangent_saturation as smooth  # noqa: E402
from verify_universal_cubic_filtered_syzygy_frontier import (  # noqa: E402
    singular_polynomial,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "cubic_formal_gauge_cokernel_atlas.json"
)

EXPECTED_NUMERATORS = {
    "smooth": [1, -3, 3, -1],
    "nodal": [1, -1],
    "cuspidal": [2, -2],
    "line-transverse-conic": [2, -2],
    "line-tangent-conic": [3, -3],
    "triangle": [3, -3],
    "concurrent-lines": [4, -4],
    "double-line": [5, -4],
    "triple-line": [7, -5],
    "zero": [10, -6],
}

EXPECTED_ANNIHILATORS = {
    "smooth": ("z,y,x", "(x,y,z)"),
    "nodal": ("x", "(x)"),
    "cuspidal": ("x2", "(x^2)"),
    "line-transverse-conic": ("yz", "(yz)"),
    "line-tangent-conic": ("y3", "(y^3)"),
    "triangle": ("xyz", "(xyz)"),
    "concurrent-lines": ("x3", "(x^3)"),
    "double-line": ("0", "(0)"),
    "triple-line": ("0", "(0)"),
    "zero": ("0", "(0)"),
}


def symbol_tensor(cubic_symbol: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [
            cubic.polarized_value(
                cubic_symbol,
                *(
                    cubic.RELATION.cross(cubic.STANDARD_BASIS[index])
                    for index in triple
                ),
            )
            for triple in smooth.TRIPLES
        ]
    )


def singular_module(matrix: sp.Matrix) -> str:
    variables = cubic.BASE_VARIABLES
    return ",".join(
        "["
        + ",".join(
            singular_polynomial(matrix[row, column], variables)
            for row in range(matrix.rows)
        )
        + "]"
        for column in range(matrix.cols)
    )


def parse_numerator(line: str) -> list[int]:
    values = [int(value) for value in line.split(",")]
    while values and values[-1] == 0:
        values.pop()
    return values


def singular_quotient_certificate(
    compatibility: sp.Matrix,
    gauge: sp.Matrix,
) -> dict[str, Any]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    program = f"""
ring gauge_atlas=0,(x,y,z),dp;
module C={singular_module(compatibility)};
module G={singular_module(gauge)};
module K=syz(C);
module Q=std(modulo(K,G));
ideal ANN=std(quotient(G,K));
int generator_index;
int kernel_min_degree=1000000;
int kernel_max_degree=-1;
for (generator_index=1;generator_index<=size(K);generator_index++)
{{
  if (deg(K[generator_index])<kernel_min_degree)
  {{
    kernel_min_degree=deg(K[generator_index]);
  }}
  if (deg(K[generator_index])>kernel_max_degree)
  {{
    kernel_max_degree=deg(K[generator_index]);
  }}
}}
print("@@KERNEL_GENERATORS="+string(size(K)));
print("@@KERNEL_MIN_DEGREE="+string(kernel_min_degree));
print("@@KERNEL_MAX_DEGREE="+string(kernel_max_degree));
print("@@VDIM="+string(vdim(Q)));
print("@@ANNIHILATOR="+string(ANN));
print("@@HILBERT_NUMERATOR");
hilb(Q,1);
print("@@COMPLETE=1");
quit;
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "? error occurred" not in completed.stdout
    lines = completed.stdout.splitlines()

    def integer(name: str) -> int:
        prefix = f"@@{name}="
        values = [
            int(line[len(prefix) :])
            for line in lines
            if line.startswith(prefix)
        ]
        assert len(values) == 1, (name, completed.stdout)
        return values[0]

    marker_index = lines.index("@@HILBERT_NUMERATOR")
    numerator = parse_numerator(lines[marker_index + 1])
    annihilator_prefix = "@@ANNIHILATOR="
    annihilators = [
        line[len(annihilator_prefix) :]
        for line in lines
        if line.startswith(annihilator_prefix)
    ]
    assert len(annihilators) == 1
    assert "@@COMPLETE=1" in lines
    return {
        "kernel_generator_count": integer("KERNEL_GENERATORS"),
        "kernel_generator_min_degree": integer("KERNEL_MIN_DEGREE"),
        "kernel_generator_max_degree": integer("KERNEL_MAX_DEGREE"),
        "vector_space_dimension": integer("VDIM"),
        "annihilator_singular_syntax": annihilators[0],
        "hilbert_numerator_over_one_minus_t_cubed": numerator,
    }


def vanishing_order_at_one(coefficients: list[int]) -> int:
    polynomial = sp.Poly(
        sum(
            coefficient * sp.Symbol("t") ** degree
            for degree, coefficient in enumerate(coefficients)
        ),
        sp.Symbol("t"),
    )
    order = 0
    while polynomial.eval(1) == 0:
        polynomial = sp.Poly(
            sp.div(polynomial, sp.Poly(1 - sp.Symbol("t")))[0],
            sp.Symbol("t"),
        )
        order += 1
    return order


def hilbert_value(coefficients: list[int], internal_degree: int) -> int:
    """Coefficient of t^internal_degree in p(t)/(1-t)^3."""

    return sum(
        coefficient
        * (
            (internal_degree - index + 2)
            * (internal_degree - index + 1)
            // 2
        )
        for index, coefficient in enumerate(coefficients)
        if internal_degree >= index
    )


def support_multiplicity(
    coefficients: list[int], cancellation_order: int
) -> int:
    t = sp.Symbol("t")
    polynomial = sp.Poly(
        sum(
            coefficient * t**degree
            for degree, coefficient in enumerate(coefficients)
        ),
        t,
    )
    for _ in range(cancellation_order):
        quotient, remainder = sp.div(polynomial, sp.Poly(1 - t, t))
        assert remainder.is_zero
        polynomial = quotient
    return int(polynomial.eval(1))


def matrix_record(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [sp.sstr(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def main() -> None:
    cubic.FACTOR_SINGULAR_EXPRESSIONS = False
    compatibility = smooth.compatibility_matrix()
    rows: dict[str, Any] = {}
    exact_gauges: dict[str, list[list[str]]] = {}
    for name, cubic_symbol in cubic.CUBIC_STRATA.items():
        tensor = symbol_tensor(cubic_symbol)
        gauge = smooth.gauge_matrix(tensor)
        assert (
            compatibility * tensor
        ).applyfunc(sp.expand) == sp.zeros(6, 1)
        assert (
            compatibility * gauge
        ).applyfunc(sp.expand) == sp.zeros(6, 9)
        smooth.verify_dual_number_gauge_action(tensor, gauge)
        certificate = singular_quotient_certificate(
            compatibility, gauge
        )
        expected_numerator = EXPECTED_NUMERATORS[name]
        assert certificate[
            "hilbert_numerator_over_one_minus_t_cubed"
        ] == expected_numerator
        expected_annihilator, display_annihilator = (
            EXPECTED_ANNIHILATORS[name]
        )
        assert (
            certificate["annihilator_singular_syntax"]
            == expected_annihilator
        )
        assert certificate["kernel_generator_count"] == 10
        assert certificate["kernel_generator_min_degree"] == 3
        assert certificate["kernel_generator_max_degree"] == 3
        cancellation_order = vanishing_order_at_one(expected_numerator)
        support_dimension = 3 - cancellation_order
        multiplicity = support_multiplicity(
            expected_numerator, cancellation_order
        )
        quartic_nongauge_dimension = hilbert_value(
            expected_numerator, internal_degree=1
        )
        rows[name] = {
            **certificate,
            "collision_degree_shift": 3,
            "annihilator": display_annihilator,
            "support_dimension": support_dimension,
            "support_multiplicity": multiplicity,
            "generic_rank_over_Q[x,y,z]": (
                multiplicity if support_dimension == 3 else 0
            ),
            "quartic_nongauge_dimension": quartic_nongauge_dimension,
            "formal_rigidity_above_degree_three": (
                support_dimension == 0
            ),
        }
        exact_gauges[name] = matrix_record(gauge)

    assert rows["smooth"]["quartic_nongauge_dimension"] == 0
    assert {
        name: rows[name]["quartic_nongauge_dimension"]
        for name in cubic.SQUAREFREE_STRATA
        if name != "smooth"
    } == {
        "nodal": 2,
        "cuspidal": 4,
        "line-transverse-conic": 4,
        "line-tangent-conic": 6,
        "triangle": 6,
        "concurrent-lines": 8,
    }

    parameters, universal = frontier.universal_tensor()
    _directions, _lift, quartic_data = (
        smooth.universal_quartic_gauge_lift(
            smooth.gauge_matrix(smooth.central_tensor()),
            parameters,
            universal,
        )
    )
    assert quartic_data["quartic_compatible_space_dimension"] == 24

    exact_data = {
        "compatibility_matrix": matrix_record(compatibility),
        "gauge_matrices": exact_gauges,
    }
    exact_sha256 = hashlib.sha256(
        json.dumps(
            exact_data,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    artifact = {
        "schema": "cubic-formal-gauge-cokernel-atlas.v2",
        "mathematical_status": "exact graded gauge-cokernel theorem",
        "basis_conventions": {
            "collision_variables": ["x", "y", "z"],
            "tensor_component_order": [
                list(triple) for triple in smooth.TRIPLES
            ],
            "gauge_matrix_unit_order": [
                [row, column]
                for row, column in itertools.product(range(3), repeat=2)
            ],
            "external_collision_degree_shift": 3,
            "singular_monomial_order": "dp",
        },
        "universal_compatible_quartic_space_dimension": quartic_data[
            "quartic_compatible_space_dimension"
        ],
        "rows": rows,
        "exact_matrices": exact_data,
        "exact_matrix_sha256": exact_sha256,
        "proved": [
            (
                "the smooth gauge cokernel is Q in collision degree "
                "three and vanishes in all higher degrees"
            ),
            (
                "every singular squarefree symbol has a dimension-two "
                "persistent gauge cokernel with the recorded principal "
                "annihilator and multiplicity"
            ),
            (
                "double-line, triple-line, and zero symbols have "
                "dimension-three faithful gauge cokernels of generic "
                "ranks one, two, and four"
            ),
            (
                "the universal 24 quartic directions span K_4, so the "
                "recorded quartic quotient dimensions are exact essential "
                "nongauge directions"
            ),
        ],
        "not_proved": [
            (
                "a nongauge tensor direction does not imply failure of "
                "cotangent saturation"
            ),
            "formal classification of the singular-symbol deformations",
            "normality or Keller-open compatibility",
        ],
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_cubic_formal_gauge_cokernel_atlas.py"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS: determinant-twisted gauge action derived on ten symbols")
    print("PASS: exact gauge-cokernel Hilbert numerators match the atlas")
    print("PASS: exact gauge-cokernel annihilators match the atlas")
    print("PASS: singular squarefree quartic nongauge dimensions are 2,4,4,6,6,8")
    print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

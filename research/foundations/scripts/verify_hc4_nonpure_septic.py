#!/usr/bin/env python3
"""Verify every non-pure septic weighted Schur face in HC4RSD33.

The root-valuation sieve leaves eleven root partitions and only transverse
degrees four and five.  This checker retains every term of the same weight,
including the binary-linear z^3 tail in degree five.  Cross-ratio strata are
localized exactly, and Singular verifies characteristic-zero radical
membership for every coefficient of the transverse face.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "artifacts" / "generated-results" / "hc4_nonpure_septic.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


def resonance(degree: int, multiplicity: int, transverse_degree: int, order: int) -> int:
    d = degree
    m = multiplicity
    e = transverse_degree
    n = order
    return (
        d * d * m
        + d * d * n * n
        - 2 * d * e * m * n
        - 2 * d * e * m
        - d * m * m
        - d * m
        + e * e * m * m
        + 2 * e * m * m
        + m * m
    )


def root_weight(degree: int, multiplicity: int, transverse_degree: int) -> int:
    for order in range((multiplicity + 1) // 2):
        if resonance(degree, multiplicity, transverse_degree, order) == 0:
            return order
    return (multiplicity + 1) // 2


partitions = {
    "61": (6, 1),
    "52": (5, 2),
    "511": (5, 1, 1),
    "43": (4, 3),
    "421": (4, 2, 1),
    "4111": (4, 1, 1, 1),
    "331": (3, 3, 1),
    "322": (3, 2, 2),
    "3211": (3, 2, 1, 1),
    "2221": (2, 2, 2, 1),
    "22111": (2, 2, 1, 1, 1),
}

expected_weight_sums = {
    "61": (4, 4, 4, 4, 4, 4),
    "52": (4, 4, 4, 4, 4, 4),
    "511": (5, 5, 5, 5, 5, 5),
    "43": (4, 4, 4, 4, 4, 4),
    "421": (4, 4, 4, 4, 4, 4),
    "4111": (5, 5, 5, 5, 5, 5),
    "331": (5, 5, 5, 5, 5, 5),
    "322": (4, 4, 4, 4, 4, 4),
    "3211": (5, 5, 5, 5, 5, 5),
    "2221": (4, 4, 4, 4, 4, 4),
    "22111": (5, 5, 5, 5, 5, 5),
}
for name, multiplicities in partitions.items():
    assert tuple(
        sum(root_weight(7, multiplicity, degree) for multiplicity in multiplicities)
        for degree in range(6)
    ) == expected_weight_sums[name]

# The other three non-pure partitions vanish directly under the same sieve.
assert tuple(
    sum(root_weight(7, multiplicity, degree) for multiplicity in (3, 1, 1, 1, 1))
    for degree in range(6)
) == (6, 6, 6, 6, 6, 6)
assert tuple(
    sum(root_weight(7, multiplicity, degree) for multiplicity in (2, 1, 1, 1, 1, 1))
    for degree in range(6)
) == (6, 6, 6, 6, 6, 6)
assert tuple(
    sum(root_weight(7, multiplicity, degree) for multiplicity in (1, 1, 1, 1, 1, 1, 1))
    for degree in range(6)
) == (7, 7, 7, 7, 7, 7)


SINGULAR = shutil.which("Singular")
assert SINGULAR is not None, "HC4RSD33 verification requires Singular"


def singular_string(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def verify_radical_origin(
    equations: list[sp.Expr],
    ring_variables: tuple[sp.Symbol, ...],
    targets: tuple[sp.Symbol, ...],
) -> None:
    """Prove target^8 belongs to the exact characteristic-zero ideal."""

    program = "ring R=0,(" + ",".join(map(str, ring_variables)) + "),dp;\n"
    for index, equation in enumerate(equations):
        numerator, _denominator = sp.fraction(sp.together(equation))
        program += f"poly p{index}={singular_string(numerator)};\n"
    program += "ideal I=" + ",".join(
        f"p{index}" for index in range(len(equations))
    ) + ";\n"
    program += "option(redSB); ideal G=std(I);\n"
    for target in targets:
        program += (
            f'if(reduce({target}^8,G)==0){{print("PASS_{target}");}}'
            f'else{{print("FAIL_{target}");}}\n'
        )
    program += "exit;\n"
    result = subprocess.run(
        [SINGULAR, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
        timeout=300,
    )
    assert result.returncode == 0, result.stderr
    assert "?" not in result.stdout, result.stdout
    for target in targets:
        assert f"PASS_{target}" in result.stdout, result.stdout
        assert f"FAIL_{target}" not in result.stdout, result.stdout


a, cross_b, inverse = sp.symbols("a cross_b inverse")
root_forms = (x, y, x - y, x - a * y, x - cross_b * y)
verified_layers: list[dict[str, object]] = []

for name, multiplicities in partitions.items():
    top = sp.prod(
        root_forms[index] ** multiplicity
        for index, multiplicity in enumerate(multiplicities)
    )
    root_factor = sp.prod(
        root_forms[index] ** ((multiplicity + 1) // 2)
        for index, multiplicity in enumerate(multiplicities)
    )
    root_factor_degree = sp.Poly(root_factor, x, y).total_degree()
    for transverse_degree in range(6):
        if expected_weight_sums[name][transverse_degree] > transverse_degree:
            continue

        residual_degree = transverse_degree - root_factor_degree
        g_coefficients = sp.symbols(f"g0:{residual_degree + 1}")
        transverse = root_factor * sum(
            g_coefficients[index]
            * x ** (residual_degree - index)
            * y**index
            for index in range(residual_degree + 1)
        )

        q_degree = 2 * transverse_degree - 7
        q_coefficients = sp.symbols(f"q0:{q_degree + 1}")
        passive_quadratic = sum(
            q_coefficients[index] * x ** (q_degree - index) * y**index
            for index in range(q_degree + 1)
        )

        r_degree = 3 * transverse_degree - 14
        if r_degree >= 0:
            r_coefficients = sp.symbols(f"r0:{r_degree + 1}")
            passive_cubic = sum(
                r_coefficients[index] * x ** (r_degree - index) * y**index
                for index in range(r_degree + 1)
            )
        else:
            r_coefficients = ()
            passive_cubic = 0

        initial_potential = (
            top
            + z * transverse
            + z**2 * passive_quadratic / 2
            + z**3 * passive_cubic
        )
        equations = sp.Poly(
            bordered_invariant(initial_potential), x, y, z
        ).coeffs()

        moduli_count = len(multiplicities) - 3
        if moduli_count <= 0:
            localization_variables: tuple[sp.Symbol, ...] = ()
        elif moduli_count == 1:
            equations.append(inverse * a * (a - 1) - 1)
            localization_variables = (inverse, a)
        else:
            equations.append(
                inverse
                * a
                * cross_b
                * (a - 1)
                * (cross_b - 1)
                * (a - cross_b)
                - 1
            )
            localization_variables = (inverse, a, cross_b)

        coefficient_variables = (
            *q_coefficients,
            *r_coefficients,
            *g_coefficients,
        )
        verify_radical_origin(
            equations,
            (*coefficient_variables, *localization_variables),
            coefficient_variables,
        )
        verified_layers.append(
            {
                "partition": name,
                "transverse_degree": transverse_degree,
                "root_weight": expected_weight_sums[name][transverse_degree],
                "g_dimension": len(g_coefficients),
                "q_dimension": len(q_coefficients),
                "r_dimension": len(r_coefficients),
                "cross_ratio_dimension": max(moduli_count, 0),
            }
        )

assert len(verified_layers) == 17


payload = {
    "format": "hc4-nonpure-septic-v1",
    "status": {
        "id": "HC4RSD33",
        "kind": "hybrid theorem",
        "scope": "all non-pure binary septic leading forms",
        "result": (
            "the valuation sieve leaves seventeen complete weighted faces; "
            "their exact discriminant-open coefficient ideals have the full "
            "transverse coefficient origin as radical"
        ),
    },
    "direct_weight_partitions": ["1111111", "211111", "31111"],
    "weighted_layers": verified_layers,
    "same_weight_tail_rule": "deg(r_k)=7-k*(7-e)",
    "residual": "the pure-seventh leading form",
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: the septic root-valuation sieve leaves exactly seventeen faces")
print("PASS: all zero-, one-, and two-cross-ratio saturations close")
print("PASS: every same-weight quadratic and cubic transverse tail is retained")
print("THEOREM: every non-pure septic leading form is a fixed cylinder")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")

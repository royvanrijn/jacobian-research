#!/usr/bin/env python3
"""Generate and verify the finite certificate ledger for one quintic map.

The ledger is a presentation artifact, not a new proof implementation.  This
script checks its elementary row data directly and runs the four canonical
exact checkers for the map, Galois certificates, signatures, modulo-seven
partitions, and the clean Q(sqrt(-31)) Hasse row.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import warnings
from functools import reduce
from pathlib import Path

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "fixed_quintic_certificate_ledger.json"
)
NOTE = ROOT / "verified" / "FIXED_QUINTIC_ARITHMETIC_ZOO.md"
BEGIN_MARKER = "<!-- BEGIN GENERATED FIXED QUINTIC LEDGER -->"
END_MARKER = "<!-- END GENERATED FIXED QUINTIC LEDGER -->"

T = sp.symbols("T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write",
        action="store_true",
        help="refresh the JSON artifact and generated Markdown table",
    )
    parser.add_argument(
        "--emit-json",
        action="store_true",
        help="print the exact JSON artifact",
    )
    return parser.parse_args()


def normalized_polynomial(target: tuple[sp.Rational, ...]) -> sp.Poly:
    pi, b, c = target
    return sp.Poly(
        T**5
        - 5 * T**3
        - 2 * pi * b * T**2
        + 4 * pi**3 * T
        - 2 * pi**5 * c,
        T,
        domain=sp.QQ,
    )


def integral_model(poly: sp.Poly) -> sp.Poly:
    _, cleared = poly.clear_denoms(convert=True)
    return cleared.primitive()[1]


def factor_degrees_mod_prime(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    model = integral_model(poly)
    assert int(model.LC()) % prime
    assert int(sp.discriminant(model.as_expr(), T)) % prime
    factors = sp.factor_list(model.as_expr(), modulus=prime)[1]
    return tuple(
        sorted(
            (
                int(sp.degree(factor, T))
                for factor, exponent in factors
                for _ in range(exponent)
            ),
            reverse=True,
        )
    )


def q(value: int, denominator: int = 1) -> sp.Rational:
    return sp.Rational(value, denominator)


def target_text(target: tuple[sp.Rational, ...]) -> str:
    return "(" + ",".join(str(value) for value in target) + ")"


def projective_target(target: tuple[sp.Rational, ...]) -> tuple[int, ...]:
    denominator = math.lcm(*(int(value.q) for value in target))
    coordinates = (denominator,) + tuple(
        int(value * denominator) for value in target
    )
    content = reduce(math.gcd, (abs(value) for value in coordinates))
    return tuple(value // content for value in coordinates)


def row(
    *,
    purpose: str,
    target: tuple[sp.Rational, ...],
    classification: str,
    witnesses: list[str],
    real_roots: int,
    local_certificate: list[str],
    normalized_display: str,
) -> dict[str, object]:
    poly = normalized_polynomial(target)
    assert poly.degree() == 5
    assert sp.discriminant(poly.as_expr(), T)
    assert poly.count_roots(-sp.oo, sp.oo) == real_roots
    projective = projective_target(target)
    return {
        "purpose": purpose,
        "target": [str(value) for value in target],
        "target_display": target_text(target),
        "projective_target": list(projective),
        "projective_height": max(abs(value) for value in projective),
        "classification": classification,
        "normalized_polynomial": str(poly.as_expr()),
        "normalized_display": normalized_display,
        "witnesses": witnesses,
        "real_roots": real_roots,
        "local_certificate": local_certificate,
    }


def build_artifact() -> dict[str, object]:
    split_target = (q(1), q(0), q(0))
    s5_one_target = (q(1), q(-1), q(-1))
    s5_three_target = (q(-1), q(-1), q(-1))
    s5_five_target = (q(1), q(0), q(-1, 2))
    a5_target = (q(1), q(0), q(-2, 5))
    c5_target = (q(1), q(0), q(-7, 10))
    d5_target = (q(2, 5), q(-21, 10), q(2))
    f20_target = (q(1, 2), q(3, 2), q(2, 5))
    product_target = (q(1), q(-3, 2), q(-9, 2))
    hasse_target = (q(5), q(-144, 5), q(-188, 3125))

    rows = [
        row(
            purpose="split",
            target=split_target,
            classification="Q^5",
            witnesses=["T(T-1)(T+1)(T-2)(T+2)"],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T(T-1)(T+1)(T-2)(T+2)",
        ),
        row(
            purpose="signature",
            target=s5_one_target,
            classification="S_5",
            witnesses=["11:(5)", "7:(4,1)", "3:(3,2)"],
            real_roots=1,
            local_certificate=[],
            normalized_display=r"T^5-5T^3+2T^2+4T+2",
        ),
        row(
            purpose="signature",
            target=s5_three_target,
            classification="S_5",
            witnesses=["5:(5)", "43:(2,1,1,1)", "nonsquare discriminant"],
            real_roots=3,
            local_certificate=[],
            normalized_display=r"T^5-5T^3-2T^2-4T-2",
        ),
        row(
            purpose="signature",
            target=s5_five_target,
            classification="S_5",
            witnesses=["2:(5)", "7:(4,1)", "19:(3,2)"],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T^5-5T^3+4T+1",
        ),
        row(
            purpose="alternating",
            target=a5_target,
            classification="A_5",
            witnesses=["3:(5)", "23:(3,1,1)", "discriminant=232^2"],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T^5-5T^3+4T+4/5",
        ),
        row(
            purpose="cyclic",
            target=c5_target,
            classification="C_5",
            witnesses=["2:(5)", "explicit order-five automorphism"],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T^5-5T^3+4T+7/5",
        ),
        row(
            purpose="dihedral",
            target=d5_target,
            classification="D_5",
            witnesses=[
                "square discriminant",
                "pair-sum resolvent split 5+5 (both 3:(5))",
                "11:(2,2,1)",
            ],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T^5-5T^3+(42/25)T^2+(32/125)T-128/3125",
        ),
        row(
            purpose="Frobenius",
            target=f20_target,
            classification="F_20",
            witnesses=[
                "29:(5)",
                "Dummit resolvent root -13/2",
                "nonsquare discriminant",
            ],
            real_roots=5,
            local_certificate=[],
            normalized_display=r"T^5-5T^3-(3/2)T^2+(1/2)T-1/40",
        ),
        row(
            purpose="product",
            target=product_target,
            classification="K_2 x K_3",
            witnesses=[
                "(T^2+T+1)(T^3-T^2-5T+9)",
                "cubic irreducible modulo 5",
            ],
            real_roots=1,
            local_certificate=[],
            normalized_display=r"(T^2+T+1)(T^3-T^2-5T+9)",
        ),
        row(
            purpose="Hasse failure",
            target=hasse_target,
            classification="irreducible 2+3",
            witnesses=[
                "common quadratic resolvent Q(sqrt(-31))",
                "cubic irreducible modulo 5",
            ],
            real_roots=1,
            local_certificate=[
                "2: quadratic splits",
                "31: cubic simple root 15",
                "all other finite primes: unramified common-resolvent argument",
            ],
            normalized_display=r"(T^2-8T+47)(T^3+8T^2+12T+8)",
        ),
    ]

    by_purpose = {entry["purpose"]: entry for entry in rows}
    assert normalized_polynomial(split_target).as_expr() == sp.expand(
        T * (T - 1) * (T + 1) * (T - 2) * (T + 2)
    )
    assert sp.discriminant(normalized_polynomial(a5_target).as_expr(), T) == 232**2
    assert factor_degrees_mod_prime(normalized_polynomial(a5_target), 3) == (5,)
    assert factor_degrees_mod_prime(normalized_polynomial(a5_target), 23) == (
        3,
        1,
        1,
    )

    signature_patterns = (
        (s5_one_target, {11: (5,), 7: (4, 1), 3: (3, 2)}),
        (s5_three_target, {5: (5,), 43: (2, 1, 1, 1)}),
        (s5_five_target, {2: (5,), 7: (4, 1), 19: (3, 2)}),
    )
    for target, patterns in signature_patterns:
        for prime, expected in patterns.items():
            assert factor_degrees_mod_prime(normalized_polynomial(target), prime) == expected

    assert factor_degrees_mod_prime(normalized_polynomial(c5_target), 2) == (5,)
    assert factor_degrees_mod_prime(normalized_polynomial(d5_target), 3) == (5,)
    assert factor_degrees_mod_prime(normalized_polynomial(d5_target), 11) == (
        2,
        2,
        1,
    )
    assert factor_degrees_mod_prime(normalized_polynomial(f20_target), 29) == (5,)
    assert [
        max(abs(value) for value in projective_target(target))
        for target in (
            s5_three_target,
            a5_target,
            c5_target,
            d5_target,
            f20_target,
        )
    ] == [1, 5, 10, 21, 15]

    product_poly = normalized_polynomial(product_target)
    product_factorization = (T**2 + T + 1) * (T**3 - T**2 - 5 * T + 9)
    assert sp.expand(product_poly.as_expr() - product_factorization) == 0
    assert factor_degrees_mod_prime(sp.Poly(T**3 - T**2 - 5 * T + 9, T), 5) == (
        3,
    )

    hasse_poly = normalized_polynomial(hasse_target)
    hasse_q = sp.Poly(T**2 - 8 * T + 47, T)
    hasse_h = sp.Poly(T**3 + 8 * T**2 + 12 * T + 8, T)
    assert sp.expand(hasse_poly.as_expr() - hasse_q.as_expr() * hasse_h.as_expr()) == 0
    assert sp.discriminant(hasse_q.as_expr(), T) == -31 * 2**2
    assert sp.discriminant(hasse_h.as_expr(), T) == -31 * 8**2
    assert factor_degrees_mod_prime(hasse_h, 5) == (3,)
    assert (-31) % 8 == 1
    assert int(hasse_h.eval(15)) % 31 == 0
    assert int(hasse_h.diff().eval(15)) % 31
    assert by_purpose["Hasse failure"]["real_roots"] == 1

    mod_seven = [
        {"partition": list(partition), "target": list(target)}
        for partition, target in (
            ((5,), (1, 0, 1)),
            ((4, 1), (1, 0, 3)),
            ((3, 2), (1, 0, 2)),
            ((3, 1, 1), (1, 1, 3)),
            ((2, 2, 1), (3, 2, 0)),
            ((2, 1, 1, 1), (1, 2, 6)),
            ((1, 1, 1, 1, 1), (1, 0, 0)),
        )
    ]
    for entry in mod_seven:
        target = tuple(q(value) for value in entry["target"])
        assert factor_degrees_mod_prime(normalized_polynomial(target), 7) == tuple(
            entry["partition"]
        )

    pi, b, c = sp.symbols("Pi B C")
    coefficient_map = (-2 * pi * b, 4 * pi**3, -2 * pi**5 * c)
    coefficient_jacobian = sp.factor(
        sp.Matrix(coefficient_map).jacobian((pi, b, c)).det()
    )
    assert coefficient_jacobian == -48 * pi**8

    return {
        "schema_version": 1,
        "map": {
            "seed": "S^5-5*S^3+4*S",
            "jacobian": "-2",
            "inverse_polynomial": "Pi^5*S^5-5*Pi*S^3-2*B*S^2+4*S-2*C",
            "normalized_coefficient_map": [
                "-2*Pi*B",
                "4*Pi^3",
                "-2*Pi^5*C",
            ],
            "normalized_coefficient_jacobian": "-48*Pi^8",
        },
        "rows": rows,
        "modulo_7_unramified_partitions": mod_seven,
        "source_checkers": [
            "scripts/verify_fixed_quintic_moduli_dominance.py",
            "scripts/verify_fixed_quintic_arithmetic_zoo.py",
            "scripts/verify_universal_quintic_calculator.py",
            "scripts/verify_fixed_quintic_hasse_minus_thirty_one.py",
        ],
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_fixed_quintic_certificate_ledger.py"
        ),
        "status": {
            "finite_examples": "exact",
            "universal_transitive_quintic_height_bound": 21,
            "height_minimality": "bounded exhaustive PARI computation",
            "fixed_map_hilbert_local_engineering": "externally referenced theorem",
            "infinitely_many_hasse_failures_in_this_pencil": "open",
        },
    }


def markdown_table(artifact: dict[str, object]) -> str:
    lines = [
        BEGIN_MARKER,
        "| purpose | target `(Pi,B,C)` | `H_proj` | type | witness primes / exact certificate | real | local certificate |",
        "|---|---|---:|---|---|---:|---|",
    ]
    for entry in artifact["rows"]:
        witnesses = "; ".join(entry["witnesses"])
        local = "; ".join(entry["local_certificate"]) or "—"
        lines.append(
            f"| {entry['purpose']} | `{entry['target_display']}` | "
            f"{entry['projective_height']} | `{entry['classification']}` | {witnesses} | "
            f"{entry['real_roots']} | {local} |"
        )
    lines.append(END_MARKER)
    return "\n".join(lines)


def run_source_checkers(artifact: dict[str, object]) -> None:
    for checker in artifact["source_checkers"]:
        completed = subprocess.run(
            [sys.executable, str(ROOT / checker)],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        if completed.returncode:
            raise SystemExit(
                f"{checker} failed\n{completed.stdout}{completed.stderr}"
            )


def replace_marked_block(text: str, block: str) -> str:
    assert text.count(BEGIN_MARKER) == 1
    assert text.count(END_MARKER) == 1
    before, tail = text.split(BEGIN_MARKER, 1)
    _, after = tail.split(END_MARKER, 1)
    return before + block + after


def main() -> None:
    args = parse_args()
    artifact = build_artifact()
    run_source_checkers(artifact)
    expected_json = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    expected_table = markdown_table(artifact)

    if args.emit_json:
        print(expected_json, end="")
        return

    if args.write:
        ARTIFACT.write_text(expected_json)
        NOTE.write_text(replace_marked_block(NOTE.read_text(), expected_table))
    else:
        assert ARTIFACT.read_text() == expected_json, (
            f"{ARTIFACT.relative_to(ROOT)} is stale; rerun with --write"
        )
        note = NOTE.read_text()
        assert expected_table in note, (
            f"generated ledger in {NOTE.relative_to(ROOT)} is stale; rerun with --write"
        )

    print("PASS: ten finite ledger rows are exact")
    print("PASS: every real-root count and witness-prime pattern is exact")
    print("PASS: every unramified partition of five occurs modulo 7")
    print("PASS: the clean Hasse row has exceptional primes 2 and 31")
    print("PASS: the JSON artifact and Markdown ledger agree")


if __name__ == "__main__":
    main()

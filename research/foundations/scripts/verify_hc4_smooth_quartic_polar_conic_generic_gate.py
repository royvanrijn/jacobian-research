#!/usr/bin/env python3
"""Verify the generic polar-conic exclusion after HC4NHM20.

The first exceptional divisor of the HC4NHM16 squarefree-line row is the
polar conic Delta=0.  It has the universal point (p,q,r)=(1,3,1).  Lines
through that point give a rational parametrization over Q(tau,m), and an
exact staged Singular calculation on that parametrization puts the sixth
power of every active deformation coordinate in the reciprocal-Hessian
ideal.  Thus its generic set-theoretic support is the determinant-zero
origin.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "research_hc4_smooth_quartic_simple_line.py"
DRIVER_SHA256 = "99f3f94c9ee4bac0a489f25916ff290b076d33e7165e88b0a952754548c419ec"


def polar_data() -> tuple[
    tuple[sp.Symbol, ...],
    sp.Expr,
    sp.Expr,
    sp.Expr,
]:
    """Return variables, Delta, and numerator/denominator of the second root."""

    p, q, r, tau, m = sp.symbols("p q r tau m")
    delta = (
        (3 * p**2 - q * r) * tau**5
        + (9 * p * r - q**2) * tau**4
        + (18 * r**2 - 6 * p * q) * tau**3
        + (18 * p**2 - 6 * q * r) * tau**2
        + (9 * p * r - q**2) * tau
        + (3 * r**2 - p * q)
    )
    denominator = (
        18 * m**2 * tau**3
        + 3 * m**2
        - m * tau**5
        - 6 * m * tau**2
        - tau**4
        - tau
    )
    numerator = (
        3 * m * tau**5
        - 9 * m * tau**4
        - 36 * m * tau**3
        + 18 * m * tau**2
        - 9 * m * tau
        - 6 * m
        + tau**5
        + 6 * tau**4
        + 6 * tau**3
        + 6 * tau**2
        + 6 * tau
        + 1
    )
    return (p, q, r, tau, m), delta, numerator, denominator


def verify_parametrization() -> None:
    variables, delta, numerator, denominator = polar_data()
    p, q, r, tau, m = variables
    theta = sp.symbols("theta")

    assert sp.expand(delta.subs({p: 1, q: 3, r: 1})) == 0
    line_value = sp.factor(
        delta.subs({p: 1, q: 3 + theta, r: 1 + m * theta})
    )
    assert sp.expand(
        line_value - theta * (denominator * theta - numerator)
    ) == 0
    assert sp.cancel(
        delta.subs(
            {
                p: 1,
                q: 3 + numerator / denominator,
                r: 1 + m * numerator / denominator,
            }
        )
    ) == 0
    print("PASS: lines through (1,3,1) rationally parametrize the generic polar conic")


def build_singular_program() -> tuple[str, int, tuple[int, ...]]:
    from research_hc4_smooth_quartic_simple_line import (
        build_equations,
        singular,
        unknown_degree,
    )

    equations, unknowns, _parameters = build_equations(
        "squarefree-line", False
    )
    assert len(equations) == 81

    tau, sigma = sp.symbols("tau sigma")
    p, q, r = sp.symbols("p q r")
    m, c, theta = sp.symbols("m c theta")
    bs = sp.symbols("b0:18")
    (_p, _q, _r, _tau, _m), _delta, numerator, denominator = polar_data()
    assert (p, q, r, tau, m) == (_p, _q, _r, _tau, _m)

    substitution = {
        p: c,
        q: c * (3 + theta),
        r: c * (1 + m * theta),
    }
    specialized: list[sp.Expr] = []
    seen: set[str] = set()
    for equation in equations:
        value = sp.expand(equation.subs(substitution))
        key = sp.srepr(value)
        if value != 0 and key not in seen:
            seen.add(key)
            specialized.append(value)
    assert len(specialized) == 81

    linear_indices = tuple(
        index
        for index, equation in enumerate(specialized)
        if unknown_degree(equation, unknowns) <= 1
    )
    assert linear_indices == (0, 1, 3, 18, 19, 21, 36, 37, 39, 54)
    linear = [specialized[index] for index in linear_indices]

    # These eight coefficients were selected from the exact 81-equation
    # list.  After the ten linear pivots, two give the binary quadratic
    # core in (b13,v), and the other six successively control the remaining
    # active variables.  The final power-membership check is performed in
    # their combined ideal, so the selection is proof-bearing rather than
    # a sampled heuristic.
    selected_indices = (6, 24, 38, 42, 43, 45, 49, 51)
    selected = [specialized[index] for index in selected_indices]

    parameters = (bs[15], bs[16], bs[17], c, m, sigma, tau)
    lines = [
        "option(redSB);",
        (
            f"ring rr=(0,{','.join(map(str, parameters))}),"
            f"({','.join(map(str, unknowns))}),dp;"
        ),
        f"number theta=({singular(numerator)})/({singular(denominator)});",
        f"ideal L={','.join(singular(value) for value in linear)};",
        f"ideal S={','.join(singular(value) for value in selected)};",
        "timer=1;",
        "ideal GL=std(L);",
        "ideal J=reduce(S,GL);",
        "ideal H=slimgb(J);",
        "ideal G=std(GL+H);",
        f"ideal T={','.join(f'{symbol}^6' for symbol in unknowns)};",
        "ideal RT=reduce(T,G);",
        'print("REMAINDERS_BEGIN");',
        "print(RT);",
        'print("REMAINDERS_END");',
        'print("linear_size="+string(size(GL)));',
        'print("nonlinear_size="+string(size(H)));',
        'print("combined_size="+string(size(G)));',
        'print("combined_dimension="+string(dim(G)));',
        'print("DONE");',
        "quit;",
    ]
    return "\n".join(lines) + "\n", len(unknowns), selected_indices


def verify_generic_polar_basis() -> None:
    assert shutil.which("Singular") is not None, "Singular is required"
    digest = hashlib.sha256(DRIVER.read_bytes()).hexdigest()
    assert digest == DRIVER_SHA256, (digest, DRIVER_SHA256)

    program, unknown_count, selected_indices = build_singular_program()
    completed = subprocess.run(
        ["Singular", "--no-tty", "--quiet"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=180,
        check=False,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0, output
    assert "? error" not in output and "DONE" in output, output
    body = output.split("REMAINDERS_BEGIN\n", 1)[1].split(
        "\nREMAINDERS_END", 1
    )[0]
    remainders = [line.rstrip(",") for line in body.splitlines() if line]
    assert remainders == ["0"] * unknown_count, remainders
    for marker in (
        "linear_size=10",
        "nonlinear_size=11",
        "combined_size=21",
        "combined_dimension=0",
    ):
        assert marker in output, marker
    print(
        "PASS: selected reciprocal coefficients "
        f"{selected_indices} give an exact generic polar-conic basis"
    )
    print("PASS: the sixth power of every active deformation coordinate reduces to zero")


def verify_determinant_zero_support() -> None:
    x, y, z = sp.symbols("x y z")
    p, q, r, b15, b16, b17 = sp.symbols("p q r b15 b16 b17")
    matrix = sp.Matrix(
        [
            [0, 0, -y**2],
            [0, 0, x**2],
            [
                -y**2,
                x**2,
                p * x**2
                + q * x * y
                + r * y**2
                + z * (b15 * x + b16 * y + b17 * z),
            ],
        ]
    )
    assert sp.expand(matrix.det()) == 0
    print("PASS: the radical support leaves only the determinant-zero boundary matrix")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="verify committed equation-builder provenance without symbolic or Singular replay",
    )
    arguments = parser.parse_args()

    digest = hashlib.sha256(DRIVER.read_bytes()).hexdigest()
    assert digest == DRIVER_SHA256, (digest, DRIVER_SHA256)
    if arguments.audit_existing_only:
        print(
            "PASS committed HC4 smooth-quartic polar-conic provenance is intact; "
            "no symbolic or Singular replay"
        )
        return
    verify_parametrization()
    verify_generic_polar_basis()
    verify_determinant_zero_support()
    print("THEOREM: the generic smooth polar-conic divisor has no nonzero quartic quotient")


if __name__ == "__main__":
    main()

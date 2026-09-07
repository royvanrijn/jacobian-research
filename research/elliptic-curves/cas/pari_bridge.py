#!/usr/bin/env python3
"""Small subprocess bridge to the installed PARI/GP executable."""

from __future__ import annotations

from research_runtime.supervisor import Limits, capture, capture_record, captured_run, run as supervised_run

from fractions import Fraction
import shutil
import subprocess
from typing import Any, Sequence

from ek_k3 import rational_to_string


def _gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(value)})"


def pari_version(timeout: float = 5.0) -> str | None:
    executable = shutil.which("gp")
    if executable is None:
        return None
    result = captured_run(
        [executable, "-q"],
        input="print(version());\nquit\n",
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    return result.stdout.strip()


def minimal_curve_data(
    coefficients: Sequence[Fraction],
    *,
    timeout: float = 30.0,
    rank_effort: int | None = None,
    known_points: Sequence[tuple[Fraction, Fraction]] | None = None,
    local_primes: Sequence[int] = (),
    stack_bytes: int = 256_000_000,
) -> dict[str, Any]:
    """Return exact minimal-model and conductor data from PARI/GP.

    ``ellrank`` is optional because even a small effort can be expensive.  Its
    output is recorded as a computational lower/upper bound, not promoted to
    a theorem by this bridge.  If ``known_points`` are supplied they are passed
    to ``ellrank`` and separately checked with ``ellisoncurve``.
    """

    if len(coefficients) != 5:
        raise ValueError("an extended Weierstrass vector has five coefficients")
    if stack_bytes < 8_000_000:
        raise ValueError("the PARI stack must be at least 8,000,000 bytes")
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    vector = ",".join(_gp_rational(value) for value in coefficients)
    points = tuple(known_points or ())
    for point in points:
        if len(point) != 2:
            raise ValueError("an affine point has two coordinates")
    primes = tuple(local_primes)
    if any(prime < 2 for prime in primes):
        raise ValueError("local reduction primes must be at least two")
    commands = [
        "default(realprecision,60);",
        f"E=ellinit([{vector}]);",
        "Em=ellminimalmodel(E);",
        "G=ellglobalred(Em);",
        'print("MODEL_BEGIN");',
        "print(Em.a1);print(Em.a2);print(Em.a3);print(Em.a4);print(Em.a6);",
        'print("MODEL_END");',
        'print("CONDUCTOR_BEGIN");print(G[1]);print("CONDUCTOR_END");',
        'print("LOG_CONDUCTOR_BEGIN");print(log(G[1]));print("LOG_CONDUCTOR_END");',
        'print("DISCRIMINANT_BEGIN");print(Em.disc);print("DISCRIMINANT_END");',
        'print("ROOT_NUMBER_BEGIN");print(ellrootno(Em));print("ROOT_NUMBER_END");',
    ]
    for prime in primes:
        commands.extend(
            [
                f"L=elllocalred(Em,{prime});",
                f'print("LOCAL_{prime}_BEGIN");',
                (
                    f"print(L[1]);print(L[2]);print(L[4]);"
                    f"print(valuation(Em.c4,{prime}));"
                    f"print(valuation(Em.disc,{prime}));print(ellap(Em,{prime}));"
                ),
                f'print("LOCAL_{prime}_END");',
            ]
        )
    if points:
        gp_points = ",".join(
            f"[{_gp_rational(x_value)},{_gp_rational(y_value)}]"
            for x_value, y_value in points
        )
        commands.extend(
            [
                f"P=[{gp_points}];",
                'print("POINTS_BEGIN");',
                "print(vecsum(vector(#P,i,ellisoncurve(E,P[i]))));",
                "print(matdet(ellheightmatrix(E,P)));",
                'print("POINTS_END");',
            ]
        )
    if rank_effort is not None:
        if rank_effort < 0:
            raise ValueError("PARI rank effort must be nonnegative")
        rank_call = (
            f"R=ellrank(E,{rank_effort},P);"
            if points
            else f"R=ellrank(E,{rank_effort});"
        )
        commands.extend(
            [
                rank_call,
                'print("RANK_BEGIN");print(R[1]);print(R[2]);print(#R[4]);print("RANK_END");',
            ]
        )
    commands.append("quit")
    result = captured_run(
        [executable, "-q", "-s", str(stack_bytes)],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP failed: {result.stderr.strip()}")

    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]

    def block(name: str) -> list[str]:
        start = lines.index(f"{name}_BEGIN") + 1
        end = lines.index(f"{name}_END")
        return lines[start:end]

    model = tuple(int(value) for value in block("MODEL"))
    answer: dict[str, Any] = {
        "minimal_model": model,
        "conductor": int(block("CONDUCTOR")[0]),
        "log_conductor": block("LOG_CONDUCTOR")[0],
        "minimal_discriminant": int(block("DISCRIMINANT")[0]),
        "root_number": int(block("ROOT_NUMBER")[0]),
    }
    if primes:
        answer["local_reduction"] = {
            str(prime): {
                "conductor_exponent": int(block(f"LOCAL_{prime}")[0]),
                "kodaira_code": int(block(f"LOCAL_{prime}")[1]),
                "tamagawa_number": int(block(f"LOCAL_{prime}")[2]),
                "minimal_c4_valuation": int(block(f"LOCAL_{prime}")[3]),
                "minimal_discriminant_valuation": int(block(f"LOCAL_{prime}")[4]),
                "ellap": int(block(f"LOCAL_{prime}")[5]),
            }
            for prime in primes
        }
    if points:
        point_data = block("POINTS")
        points_on_curve = int(point_data[0])
        answer["supplied_points"] = {
            "count": len(points),
            "on_curve_count": points_on_curve,
            "all_on_curve": points_on_curve == len(points),
            "height_pairing_determinant_approx": point_data[1],
        }
    if rank_effort is not None:
        rank = block("RANK")
        answer["pari_ellrank"] = {
            "lower_bound": int(rank[0]),
            "upper_bound": int(rank[1]),
            "returned_independent_points": int(rank[2]),
            "effort": rank_effort,
        }
    return answer

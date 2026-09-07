#!/usr/bin/env python3
"""Exact rank-17 certificate for the Nagao-family lead ``T=956/9``.

The script independently replays the original-quartic height-one-million
search, checks and maps all returned points over ``QQ``, and reproduces the
reported precision-stable 17-point subset.  That numerical selection is then
small-prime saturated.  Exact finite reductions of the returned basis certify
rank at least 17.  The conductor and root number are recomputed independently.

Every PARI process runs synchronously in a fresh process group with a strict
timeout and joined TERM/KILL cleanup.  Numerical height rank is used only for
selection and is not itself a rank certificate.
"""

from __future__ import annotations

from research_runtime.supervisor import Limits, capture, capture_record, captured_run, run as supervised_run

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time
from typing import Any, Iterable, Sequence

from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from nagao_linear_sections import point_on_short_curve


Q = Fraction
PARAMETER_T = Q(956, 9)
UNIFORM_HEIGHT = 1_000_000
EXPECTED_SIGNED_POINTS = 110
EXPECTED_NEW_IMAGES = 43
EXPECTED_POOL_SIZE = 55
EXPECTED_HEIGHT_SUBSET = (
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    13,
    14,
    15,
    16,
    17,
    23,
)
TARGET_LOG_CONDUCTOR = "182.72"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/certify_nagao_rank21_t956.py "
    "--output archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t956_rank17_certificate.json"
)


def gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(Q(value))})"


def gp_point(point: tuple[Fraction, Fraction]) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def gp_curve(coefficients: Sequence[Fraction]) -> str:
    return ",".join(gp_rational(Q(value)) for value in coefficients)


def gp_quartic(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"{gp_rational(Q(coefficient))}*x^{power}"
        for power, coefficient in enumerate(coefficients)
    )


def point_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    text = "\n".join(
        f"{rational_to_string(x_value)},{rational_to_string(y_value)}"
        for x_value, y_value in points
    )
    return hashlib.sha256(text.encode()).hexdigest()




def run_gp_capped(program, *, timeout, stack_bytes):
    if timeout<=0 or timeout>60:raise ValueError("PARI timeout must lie in (0,60]")
    executable=shutil.which('gp')
    if executable is None:raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    result=capture([executable,'-q','-s',str(stack_bytes)],input_text=program,
        limits=Limits(timeout,max(512_000_000,2*stack_bytes),pari_stack_bytes=stack_bytes),
        separate_stderr=True,check=False)
    if result.returncode or '***' in result.stderr:
        raise RuntimeError(f"PARI/GP failed: {' '.join(result.stderr.split())[:1000]}")
    return result.stdout,result.supervision['wall_seconds']


POINT_PATTERN = re.compile(r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]")


def parse_points(text: str) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((Q(x_value), Q(y_value)) for x_value, y_value in POINT_PATTERN.findall(text))


def signless_points(
    points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    answer: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in points:
        answer.setdefault(point[0], point)
    return tuple(answer.values())


def uniform_quartic_search(
    coefficients: Sequence[Fraction],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    program = "\n".join(
        (
            f"Q={gp_quartic(coefficients)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{UNIFORM_HEIGHT});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp_capped(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    marker = re.search(r"POINTS_BEGIN\n(.*?)\nPOINTS_END", output, re.DOTALL)
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if marker is None or milliseconds is None:
        raise AssertionError("PARI omitted uniform-search markers")
    return parse_points(marker.group(1)), {
        "status": "completed",
        "height_bound": UNIFORM_HEIGHT,
        "pari_milliseconds": int(milliseconds.group(1)),
        "wall_seconds": wall_seconds,
        "timeout_seconds": timeout,
        "pari_stack_bytes": stack_bytes,
    }


def _parse_index_vector(text: str) -> list[int]:
    match = re.search(r"\[(.*?)\]", text)
    if match is None:
        raise AssertionError("PARI omitted a subset vector")
    return [int(value.strip()) for value in match.group(1).split(",") if value.strip()]


def height_replay(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a height input point is off the exact curve")
    commands = [
        f"E=ellinit([{gp_curve(coefficients)}]);",
        f"P=[{','.join(gp_point(point) for point in points)}];",
    ]
    for precision in (72, 120):
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);IX=matindexrank(H);K=vecextract(P,IX[2]);HK=ellheightmatrix(E,K);",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(H));print(Vec(IX[2]));print(matdet(HK));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.append("quit")
    output, wall_seconds = run_gp_capped(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    records = []
    for precision in (72, 120):
        start = lines.index(f"HEIGHT_{precision}_BEGIN") + 1
        end = lines.index(f"HEIGHT_{precision}_END")
        values = lines[start:end]
        records.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": _parse_index_vector(values[1]),
                "subset_height_determinant": values[2],
                "gp_process_wall_seconds": wall_seconds,
            }
        )
    if len({record["numerical_rank"] for record in records}) != 1 or len(
        {tuple(record["subset_indices_one_based"]) for record in records}
    ) != 1:
        raise AssertionError("height selection changed with precision")
    return tuple(records)


def saturate_basis(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    program = "\n".join(
        (
            "default(realprecision,100);",
            f"E=ellinit([{gp_curve(coefficients)}]);",
            f"P=[{','.join(gp_point(point) for point in points)}];",
            "H0=ellheightmatrix(E,P);gettime();S=ellsaturation(E,P,20);",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("RETURNED_COUNT ",#S);',
            'print("ON_CURVE ",vecsum(vector(#S,i,ellisoncurve(E,S[i]))));',
            'print("DET_RATIO ",matdet(H0)/matdet(ellheightmatrix(E,S)));',
            'print("SATURATED_POINTS_BEGIN");print(S);print("SATURATED_POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp_capped(
        program, timeout=timeout, stack_bytes=stack_bytes
    )

    def value(label: str) -> str:
        match = re.search(rf"^{label} (.+)$", output, re.MULTILINE)
        if match is None:
            raise AssertionError(f"PARI omitted {label}")
        return match.group(1)

    marker = re.search(
        r"SATURATED_POINTS_BEGIN\n(.*?)\nSATURATED_POINTS_END", output, re.DOTALL
    )
    if marker is None:
        raise AssertionError("PARI omitted saturated basis")
    saturated = parse_points(marker.group(1))
    returned = int(value("RETURNED_COUNT"))
    on_curve = int(value("ON_CURVE"))
    if len(saturated) != returned or on_curve != returned or any(
        not point_on_short_curve(coefficients, point) for point in saturated
    ):
        raise AssertionError("PARI returned an invalid saturated basis")
    return saturated, {
        "status": "completed",
        "prime_bound_strict_upper_limit": 20,
        "input_point_count": len(points),
        "returned_point_count": returned,
        "exact_returned_points_on_curve": on_curve,
        "pari_milliseconds": int(value("PARI_MILLISECONDS")),
        "wall_seconds": wall_seconds,
        "timeout_seconds": timeout,
        "height_determinant_ratio": value("DET_RATIO"),
        "scope_warning": (
            "PARI documents ellsaturation under a finite-index hypothesis; exact finite reductions below certify independence of the returned basis"
        ),
    }


def conductor_probe(
    coefficients: Sequence[Fraction],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellminimalmodel(ellinit([{gp_curve(coefficients)}]));G=ellglobalred(E);",
            'print("MODEL_BEGIN");print(E.a1);print(E.a2);print(E.a3);print(E.a4);print(E.a6);print("MODEL_END");',
            'print("CONDUCTOR ",G[1]);print("LOG_CONDUCTOR ",log(G[1]));',
            'print("DISCRIMINANT ",E.disc);print("ROOT_NUMBER ",ellrootno(E));',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp_capped(
        program, timeout=timeout, stack_bytes=stack_bytes
    )

    def value(label: str) -> str:
        match = re.search(rf"^{label} (.+)$", output, re.MULTILINE)
        if match is None:
            raise AssertionError(f"PARI omitted {label}")
        return match.group(1)

    lines = [line.strip() for line in output.splitlines() if line.strip()]
    start = lines.index("MODEL_BEGIN") + 1
    end = lines.index("MODEL_END")
    return {
        "minimal_model": [int(item) for item in lines[start:end]],
        "conductor": value("CONDUCTOR"),
        "log_conductor": value("LOG_CONDUCTOR"),
        "minimal_discriminant": value("DISCRIMINANT"),
        "root_number": int(value("ROOT_NUMBER")),
        "wall_seconds": wall_seconds,
        "timeout_seconds": timeout,
    }


def signature_records(signatures: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": signature.prime,
            "group_order": signature.group_order,
            "doubled_subgroup_order": signature.doubled_subgroup_order,
            "quotient_dimension": signature.quotient_dimension,
            "rows": [list(row) for row in signature.rows],
        }
        for signature in signatures
    ]


def build_certificate(args: argparse.Namespace) -> dict[str, Any]:
    quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    short = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    raw_points, uniform_record = uniform_quartic_search(
        quartic, timeout=args.search_timeout, stack_bytes=args.stack_bytes
    )
    if len(raw_points) != EXPECTED_SIGNED_POINTS or any(
        y_value**2 != quartic_value(quartic, x_value)
        for x_value, y_value in raw_points
    ):
        raise AssertionError("the exact uniform quartic search changed")
    signless = signless_points(raw_points)
    visible_quartic = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
    visible_x = {point[0] for point in visible_quartic}
    visible_images = tuple(
        quartic_point_to_short_jacobian(RANK21_CONSTRUCTION, PARAMETER_T, point)
        for point in visible_quartic
    )
    seen_image_x = {point[0] for point in visible_images}
    new_images = []
    new_records = []
    for quartic_point in signless:
        if quartic_point[0] in visible_x or not quartic_point[1]:
            continue
        image = quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        if not point_on_short_curve(short, image):
            raise AssertionError("a mapped image missed the exact Jacobian")
        if image[0] in seen_image_x:
            continue
        seen_image_x.add(image[0])
        new_images.append(image)
        new_records.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "jacobian_x": rational_to_string(image[0]),
                "jacobian_y": rational_to_string(image[1]),
                "exact_quartic_and_jacobian_membership_checked": True,
            }
        )
    if len(new_images) != EXPECTED_NEW_IMAGES:
        raise AssertionError("the exact new-image count changed")
    pool = visible_images + tuple(new_images)
    if len(pool) != EXPECTED_POOL_SIZE:
        raise AssertionError("the exact point pool changed")
    height_runs = height_replay(
        short, pool, timeout=args.height_timeout, stack_bytes=args.stack_bytes
    )
    if height_runs[-1]["numerical_rank"] != 17 or tuple(
        height_runs[-1]["subset_indices_one_based"]
    ) != EXPECTED_HEIGHT_SUBSET:
        raise AssertionError("the reported stable height subset changed")
    selected = tuple(pool[index - 1] for index in EXPECTED_HEIGHT_SUBSET)
    saturated, saturation = saturate_basis(
        short,
        selected,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated) != 17:
        raise AssertionError("saturation changed the basis rank")
    signatures = find_mod2_reduction_certificate(short, saturated, prime_bound=500)
    exact_rank = combined_mod2_rank(signatures, 17)
    if exact_rank != 17:
        raise AssertionError("finite reductions did not certify all 17 directions")
    no_two_torsion_prime = find_two_torsion_certificate_prime(short)
    conductor = conductor_probe(
        short, timeout=args.conductor_timeout, stack_bytes=args.stack_bytes
    )
    if conductor["root_number"] != -1 or not conductor["log_conductor"].startswith(
        "140.6728210042633689"
    ):
        raise AssertionError("the pinned conductor or root number changed")
    exact_log_bound = exact_log_conductor_certificate(int(conductor["conductor"]))
    if not exact_log_bound["strict_target_proved_exactly"]:
        raise AssertionError("the exact conductor left the target range")
    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "status": "exact_rank17_certificate_complete",
        "theorem": (
            "the Nagao rank-21-family specialization T=956/9 has Mordell-Weil rank at least 17 and log conductor below 182.72"
        ),
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "short_weierstrass_coefficients": [
                rational_to_string(value) for value in short
            ],
            **conductor,
            "below_strict_log_conductor_target": True,
            "exact_log_conductor_bound": exact_log_bound,
            "target_rank": 21,
        },
        "primary_source": PRIMARY_SOURCE,
        "uniform_search": {
            **uniform_record,
            "signed_point_count": len(raw_points),
            "signless_point_count": len(signless),
            "visible_point_count": len(visible_images),
            "new_distinct_jacobian_images": len(new_images),
            "raw_point_sha256": point_digest(raw_points),
            "new_image_sha256": point_digest(new_images),
            "exact_pool_point_count": len(pool),
            "exact_pool_sha256": point_digest(pool),
            "all_memberships_checked_exactly": True,
            "new_points": new_records,
        },
        "height_selection": {
            "runs": list(height_runs),
            "stable_numerical_rank": 17,
            "selected_pool_indices_one_based": list(EXPECTED_HEIGHT_SUBSET),
            "selection_is_not_certification": True,
        },
        "exact_rank_certificate": {
            "small_prime_saturation": saturation,
            "saturated_basis_sha256": point_digest(saturated),
            "saturated_basis": [
                {
                    "jacobian_x": rational_to_string(point[0]),
                    "jacobian_y": rational_to_string(point[1]),
                    "exact_jacobian_membership_checked": True,
                }
                for point in saturated
            ],
            "two_torsion_certificate_prime": no_two_torsion_prime,
            "finite_reduction_signatures": signature_records(signatures),
            "combined_exact_rank_over_F2": exact_rank,
            "certified_algebraic_rank_lower_bound": exact_rank,
        },
        "process_safety": {
            "all_pari_processes_foreground": True,
            "fresh_process_group_per_call": True,
            "timeout_cleanup": "SIGTERM then SIGKILL after two seconds",
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": shutil.which("gp"),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t956_rank17_certificate.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "search_timeout",
        "height_timeout",
        "saturation_timeout",
        "conductor_timeout",
    ):
        if not 0 < getattr(args, name) <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must lie in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    certificate = build_certificate(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"rank>={certificate['exact_rank_certificate']['certified_algebraic_rank_lower_bound']} "
        f"logN={certificate['candidate']['log_conductor']}",
        flush=True,
    )


if __name__ == "__main__":
    main()

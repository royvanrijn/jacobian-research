#!/usr/bin/env python3
"""Direct skew-height search on the original ``T=956/9`` Nagao quartic.

This is deliberately disjoint from both the uniform height-one-million box
and all alternate-cover searches.  It runs the ten denominator slabs first
used for the ``u=42`` frontier, covering every denominator through 128000
with slab-dependent numerator bounds.  It then searches 36 determinant-one
charts centered at the twelve highest-projective-height nonvisible points in
the uniform checkpoint, using three orientations per center.

All PARI calls are one-shot foreground process groups with strict timeouts and
no retries.  Every returned point is mapped back to the original quartic and
then to its short Jacobian over ``QQ``.  Height computations are numerical
selection evidence only.  Any stable rank gain immediately triggers
small-prime saturation and exact finite-reduction certification.
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

from alternate_quartic_covers import point_on_short_curve, short_add
from certify_nagao_rank21_t956 import (
    PARAMETER_T,
    gp_curve,
    gp_point,
    gp_quartic,
    height_replay,
    parse_points,
    point_digest,
    saturate_basis,
    signless_points,
)
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
from search_nagao_u42_skew_height import (
    SEARCH_BOXES,
    SearchBox,
    centered_unimodular_matrix,
    map_chart_point,
    transform_binary_quartic,
)


Q = Fraction
CHECKPOINT_HEIGHT = 1_000_000
CHART_CENTER_COUNT = 12
CHART_SHIFTS = (0, -1, 1)
CHART_HEIGHT = 50_000
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_t956_skew.py "
    "--output archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t956_skew.json"
)




def run_gp_once(program, *, timeout, stack_bytes):
    if timeout<=0 or timeout>60:raise ValueError("PARI timeout must lie in (0,60]")
    executable=shutil.which('gp')
    if executable is None:return None,{'status':'unavailable','wall_seconds':0.0}
    record=capture_record([executable,'-q','-s',str(stack_bytes)],input_text=program,
        limits=Limits(timeout,max(512_000_000,2*stack_bytes),pari_stack_bytes=stack_bytes))
    status=record['outcome']
    fatal=[line for line in record['stderr'].splitlines() if '***' in line]
    if status=='strict_wall_timeout':return None,{**record,'status':'timeout'}
    if status!='completed':return None,{**record,'status':status}
    if fatal or record['returncode']:
        return None,{**record,'status':'pari_error','error':' '.join(fatal or record['stderr'].splitlines())[:2000]}
    return record['stdout'],{**record,'status':'completed'}


def search_original_quartic(
    coefficients: Sequence[Fraction],
    height_specification: str,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    program = "\n".join(
        (
            f"Q={gp_quartic(coefficients)};gettime();",
            f"R=hyperellratpoints(Q,{height_specification});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, process_record = run_gp_once(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record = dict(process_record)
    record.update(
        {
            "height_specification": height_specification,
            "timeout_seconds": timeout,
            "pari_stack_bytes": stack_bytes,
            "retried": False,
        }
    )
    if output is None:
        return (), record
    marker = re.search(r"POINTS_BEGIN\n(.*?)\nPOINTS_END", output, re.DOTALL)
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if marker is None or milliseconds is None:
        raise AssertionError("PARI omitted point-search markers")
    points = parse_points(marker.group(1))
    record.update(
        {
            "pari_milliseconds": int(milliseconds.group(1)),
            "signed_point_count": len(points),
            "distinct_abscissa_count": len(signless_points(points)),
        }
    )
    return points, record


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def load_checkpoint(
    path: Path,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    dict[str, Any],
]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if Q(data["candidate"]["parameter_t"]) != PARAMETER_T:
        raise AssertionError("the input certificate has the wrong parameter")
    certificate = data["exact_rank_certificate"]
    if certificate["certified_algebraic_rank_lower_bound"] != 17:
        raise AssertionError("the input does not certify rank at least 17")
    basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in certificate["saturated_basis"]
    )
    checkpoint_quartic = tuple(
        (Q(record["quartic_x"]), Q(record["quartic_z"]))
        for record in data["uniform_search"]["new_points"]
    )
    if (
        len(basis) != 17
        or len(checkpoint_quartic) != 43
        or data["uniform_search"]["height_bound"] != CHECKPOINT_HEIGHT
    ):
        raise AssertionError("the exact checkpoint dimensions changed")
    return basis, checkpoint_quartic, data


def chart_plan(
    checkpoint_points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[str, Fraction, tuple[int, int, int, int]], ...]:
    centers = sorted(
        (point[0] for point in checkpoint_points),
        key=lambda center: (-projective_height(center), center),
    )[:CHART_CENTER_COUNT]
    charts = []
    for center_index, center in enumerate(centers):
        for shift in CHART_SHIFTS:
            charts.append(
                (
                    f"high_{center_index:02d}_shift_{shift}",
                    center,
                    centered_unimodular_matrix(center, shift),
                )
            )
    return tuple(charts)


def short_multiply(
    coefficients: Sequence[Fraction],
    point: tuple[Fraction, Fraction],
    multiplier: int,
) -> tuple[Fraction, Fraction] | None:
    if multiplier < 0:
        point = point[0], -point[1]
        multiplier = -multiplier
    answer = None
    addend: tuple[Fraction, Fraction] | None = point
    while multiplier:
        if multiplier & 1:
            answer = short_add(coefficients, answer, addend)
        addend = short_add(coefficients, addend, addend)
        multiplier >>= 1
    return answer


def exact_linear_combination(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    relation: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    if len(basis) != len(relation):
        raise ValueError("relation length differs from the basis length")
    answer = None
    for point, coefficient in zip(basis, relation):
        answer = short_add(
            coefficients,
            answer,
            short_multiply(coefficients, point, int(coefficient)),
        )
    return answer


def discover_exact_relations(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[int, ...], ...], dict[str, Any]]:
    if not points:
        return (), {"status": "completed", "point_count": 0}
    commands = [
        "default(realprecision,120);",
        f"E=ellinit([{gp_curve(coefficients)}]);",
        f"B=[{','.join(gp_point(point) for point in basis)}];",
        "H=ellheightmatrix(E,B);",
    ]
    for index, point in enumerate(points):
        commands.extend(
            (
                f"Q={gp_point(point)};V=vector(#B,j,ellheight(E,B[j],Q))~;C=round(matsolve(H,V));",
                "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                f'print("RELATION_{index} ",Vec(C)," EXACT ",S==Q);',
            )
        )
    commands.append("quit")
    output, process_record = run_gp_once(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    if output is None:
        return (), process_record
    relations = []
    for index, point in enumerate(points):
        match = re.search(
            rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$", output, re.MULTILINE
        )
        if match is None or match.group(2) != "1":
            return (), {
                **process_record,
                "status": "no_exact_relation_for_every_point",
                "first_unresolved_index": index,
            }
        relation = tuple(int(value.strip()) for value in match.group(1).split(","))
        if exact_linear_combination(coefficients, basis, relation) != point:
            raise AssertionError("a PARI relation failed exact Fraction replay")
        relations.append(relation)
    return tuple(relations), {
        **process_record,
        "point_count": len(points),
        "all_relations_replayed_exactly": True,
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


def exact_gain_attempt(
    coefficients: Sequence[Fraction],
    pool: Sequence[tuple[Fraction, Fraction]],
    height_runs: Sequence[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    stable_rank = int(height_runs[-1]["numerical_rank"])
    if stable_rank <= 17:
        return {"status": "not_triggered", "stable_numerical_rank": stable_rank}
    indices = height_runs[-1]["subset_indices_one_based"]
    selected = tuple(pool[index - 1] for index in indices)
    saturated, saturation = saturate_basis(
        coefficients, selected, timeout=timeout, stack_bytes=stack_bytes
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, saturated, prime_bound=500
    )
    exact_rank = combined_mod2_rank(signatures, len(saturated))
    result: dict[str, Any] = {
        "status": (
            "certified" if exact_rank == len(saturated) else "finite_reductions_rank_deficient"
        ),
        "small_prime_saturation": saturation,
        "saturated_basis": [
            {
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_membership_checked": True,
            }
            for point in saturated
        ],
        "finite_reduction_signatures": signature_records(signatures),
        "combined_mod2_rank": exact_rank,
    }
    if exact_rank == len(saturated):
        result["no_rational_2_torsion_prime"] = find_two_torsion_certificate_prime(
            coefficients
        )
        result["certified_algebraic_rank_lower_bound"] = exact_rank
    return result


def build_search(args: argparse.Namespace) -> dict[str, Any]:
    basis, checkpoint_new, certificate = load_checkpoint(args.certificate_input)
    quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    short = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    if any(not point_on_short_curve(short, point) for point in basis):
        raise AssertionError("a certified basis point missed the exact curve")
    visible_quartic = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
    checkpoint_x = {point[0] for point in visible_quartic}
    checkpoint_x.update(point[0] for point in checkpoint_new)
    if len(checkpoint_x) != 55:
        raise AssertionError("the uniform checkpoint did not have 55 abscissae")

    visible_images = tuple(
        quartic_point_to_short_jacobian(RANK21_CONSTRUCTION, PARAMETER_T, point)
        for point in visible_quartic
    )
    checkpoint_images = visible_images + tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in certificate["uniform_search"]["new_points"]
    )
    seen_image_x = {point[0] for point in checkpoint_images}

    discovered: dict[Fraction, dict[str, Any]] = {}
    box_records = []
    for search_box in SEARCH_BOXES:
        raw_points, process_record = search_original_quartic(
            quartic,
            search_box.gp_height,
            timeout=args.box_timeout,
            stack_bytes=args.stack_bytes,
        )
        outside = []
        for point in signless_points(raw_points):
            if point[1] ** 2 != quartic_value(quartic, point[0]):
                raise AssertionError("a skew-box point missed the exact quartic")
            if point[0] in checkpoint_x:
                continue
            outside.append(point[0])
            entry = discovered.setdefault(
                point[0], {"quartic_point": point, "sources": []}
            )
            entry["sources"].append(f"skew:{search_box.identifier}")
        box_records.append(
            {
                "id": search_box.identifier,
                "numerator_absolute_bound": search_box.numerator_bound,
                "denominator_lower_bound": search_box.denominator_lower,
                "denominator_upper_bound": search_box.denominator_upper,
                **process_record,
                "outside_uniform_checkpoint_x_values": [
                    rational_to_string(value) for value in outside
                ],
            }
        )
        print(
            f"{search_box.identifier}: status={process_record['status']} outside={len(outside)}",
            flush=True,
        )

    charts = chart_plan(checkpoint_new)
    chart_records = []
    for identifier, center, matrix in charts:
        transformed = transform_binary_quartic(quartic, matrix)
        raw_points, process_record = search_original_quartic(
            transformed,
            str(CHART_HEIGHT),
            timeout=args.chart_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
        for transformed_point in signless_points(raw_points):
            mapped = map_chart_point(transformed_point, matrix)
            if mapped is None:
                continue
            if mapped[1] ** 2 != quartic_value(quartic, mapped[0]):
                raise AssertionError("a chart point missed the original quartic")
            mapped_by_x.setdefault(mapped[0], mapped)
        outside = []
        for x_value, point in mapped_by_x.items():
            if x_value in checkpoint_x:
                continue
            outside.append(x_value)
            entry = discovered.setdefault(
                x_value, {"quartic_point": point, "sources": []}
            )
            entry["sources"].append(f"chart:{identifier}")
        chart_records.append(
            {
                "id": identifier,
                "center": rational_to_string(center),
                "center_projective_height": projective_height(center),
                "matrix_a_b_c_d": list(matrix),
                "determinant": matrix[0] * matrix[3] - matrix[1] * matrix[2],
                **process_record,
                "distinct_finite_original_x_values": len(mapped_by_x),
                "outside_uniform_checkpoint_x_values": [
                    rational_to_string(value) for value in outside
                ],
            }
        )
        print(
            f"{identifier}: status={process_record['status']} outside={len(outside)}",
            flush=True,
        )

    image_records = []
    new_images = []
    for x_value in sorted(discovered):
        entry = discovered[x_value]
        quartic_point = entry["quartic_point"]
        image = quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        if not point_on_short_curve(short, image):
            raise AssertionError("a discovered point missed the exact Jacobian")
        duplicate_image = image[0] in seen_image_x
        record = {
            "quartic_x": rational_to_string(quartic_point[0]),
            "quartic_z": rational_to_string(quartic_point[1]),
            "jacobian_x": rational_to_string(image[0]),
            "jacobian_y": rational_to_string(image[1]),
            "sources": entry["sources"],
            "duplicate_checkpoint_or_prior_jacobian_sign_pair": duplicate_image,
            "exact_quartic_and_jacobian_membership_checked": True,
        }
        image_records.append(record)
        if not duplicate_image:
            seen_image_x.add(image[0])
            new_images.append(image)

    pool = basis + tuple(new_images)
    height_runs = height_replay(
        short, pool, timeout=args.height_timeout, stack_bytes=args.stack_bytes
    )
    stable_rank = int(height_runs[-1]["numerical_rank"])
    relations, relation_record = discover_exact_relations(
        short,
        basis,
        new_images,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if relations:
        for record, relation in zip(
            (record for record in image_records if not record["duplicate_checkpoint_or_prior_jacobian_sign_pair"]),
            relations,
        ):
            record["certified_basis_relation"] = list(relation)
            record["relation_replayed_exactly"] = True
    exact_gain = exact_gain_attempt(
        short,
        pool,
        height_runs,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "status": (
            "rank gain exact-certificate attempt completed"
            if stable_rank > 17
            else "bounded direct skew search complete; no numerical rank gain"
        ),
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "conductor": certificate["candidate"]["conductor"],
            "log_conductor": certificate["candidate"]["log_conductor"],
            "root_number": certificate["candidate"]["root_number"],
            "certified_input_rank_lower_bound": 17,
            "target_rank": 21,
        },
        "primary_source": PRIMARY_SOURCE,
        "certificate_input": {
            "path": str(args.certificate_input),
            "sha256": hashlib.sha256(args.certificate_input.read_bytes()).hexdigest(),
            "saturated_basis_sha256": point_digest(basis),
        },
        "uniform_checkpoint": {
            "height_bound": CHECKPOINT_HEIGHT,
            "distinct_abscissa_count": len(checkpoint_x),
            "excluded_from_all_new-yield_counts": True,
        },
        "skew_staircase": {
            "scope": (
                "ten disjoint denominator slabs covering every denominator 1..128000 with the recorded slab-dependent numerator bounds"
            ),
            "boxes": box_records,
        },
        "unimodular_charts": {
            "selection_rule": (
                "twelve nonvisible H=1000000 points with greatest max(abs(numerator(x)),denominator(x)); three shifts 0,-1,1 per center"
            ),
            "center_count": CHART_CENTER_COUNT,
            "shift_values": list(CHART_SHIFTS),
            "transformed_height_bound": CHART_HEIGHT,
            "chart_count": len(charts),
            "records": chart_records,
        },
        "new_point_analysis": {
            "outside_checkpoint_quartic_abscissa_count": len(discovered),
            "new_jacobian_sign_pair_count": len(new_images),
            "records": image_records,
            "relation_discovery": relation_record,
            "all_new_images_exactly_in_certified_rank17_span": (
                len(relations) == len(new_images)
            ),
        },
        "height_selection": {
            "pool_point_count": len(pool),
            "runs": list(height_runs),
            "stable_numerical_rank": stable_rank,
            "stable_numerical_rank_gain": stable_rank - 17,
            "selection_is_not_certification": True,
        },
        "exact_rank_gain_attempt": exact_gain,
        "bounded_scope": {
            "box_timeout_seconds_each": args.box_timeout,
            "chart_timeout_seconds_each": args.chart_timeout,
            "height_timeout_seconds": args.height_timeout,
            "relation_timeout_seconds": args.relation_timeout,
            "saturation_timeout_seconds_if_triggered": args.saturation_timeout,
            "no_retries": True,
            "all_pari_processes_foreground": True,
            "fresh_process_group_per_call": True,
            "timeout_cleanup": "SIGTERM then SIGKILL after two seconds",
        },
        "interpretation": {
            "exact": (
                "all stored quartic and Jacobian points were checked over QQ; stored basis relations, when present, were replayed with exact Fraction group arithmetic"
            ),
            "numerical": "height rank was replayed at 72 and 120 decimal digits",
            "not_claimed": (
                "a negative bounded search is not an upper rank bound and says nothing outside the explicit slabs and charts"
            ),
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
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=root
        / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t956_rank17_certificate.json",
    )
    parser.add_argument("--box-timeout", type=float, default=25.0)
    parser.add_argument("--chart-timeout", type=float, default=8.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--relation-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t956_skew.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "box_timeout",
        "chart_timeout",
        "height_timeout",
        "relation_timeout",
        "saturation_timeout",
    ):
        if not 0 < getattr(args, name) <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must lie in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    result = build_search(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"outside={result['new_point_analysis']['outside_checkpoint_quartic_abscissa_count']} "
        f"new_images={result['new_point_analysis']['new_jacobian_sign_pair_count']} "
        f"stable_rank={result['height_selection']['stable_numerical_rank']}",
        flush=True,
    )


if __name__ == "__main__":
    main()


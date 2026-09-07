#!/usr/bin/env python3
"""Exact auxiliary-Jacobian group orbits for two open Nagao section-7 slices.

The sibling accidental-slice search found rational points of height at most
200,000 on sixteen genus-one curves obtained by imposing ``x=m*T+n`` on
Nagao's section-7 quartic.  Three productive slices were already exhausted by
an independent full ``{-1,0,1}^r`` doubled-point search.  This script treats
the two remaining slices having more than sixteen known affine points:

``(m,n)=(1,-4471/339)`` and ``(-1,154687/447)``.

PARI's development ``ellfromeqn`` map realizes each quartic as a degree-four
cover of its auxiliary Jacobian.  Both signs of every pinned H=200,000 point
give a stable numerical rank-seven subgroup.  ``ellsaturation`` through 20
lowers each height determinant by 64^2.  For every coefficient vector in
``{-1,0,1}^7`` we form ``P0+2Q`` and pull it back by factoring the exact
quartic equation for the map's x-coordinate.  Thus this is a group-law orbit,
not a repeated rational-point box search.

Every returned parameter is replayed on the original quartic, canonicalized
using the even symmetry T -> -T, and decontaminated against all 21 generic
section abscissae and the sibling population.  Only projective height greater
than 200,000 survives.  Exact discriminant-radical proxies precede a bounded
exact-conductor pass.  Numerical auxiliary rank is explicitly heuristic.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from nagao_1994 import quartic_value, short_jacobian_coefficients
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    SECTION7_CONSTRUCTOR_PARAMETER,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import minimal_curve_data
from search_nagao_rank20_t5081_neighborhood import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from search_nagao_rank21_accidental_slices import (
    build_slices,
    select_minimum_intercept_priority_slices,
)


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

Q = Fraction
T0 = SECTION7_CONSTRUCTOR_PARAMETER
NAIVE_HEIGHT_EXCLUSION = 200_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")
PROXY_GATE = 190.0
PROXY_TRIAL_BOUND = 2_000
DEFAULT_CONDUCTOR_KEEP = 32
EXPECTED_INPUT_SHA256 = (
    "125a6b0df7941099547039302b6f1878b5009dcde774328527952699877b1670"
)
EXPECTED_INPUT_SCRIPT_SHA256 = (
    "edbc3c179e498cdf76426f4251ac0177cc7c3f1e59e405aba9a0a76c58da3258"
)
EXPECTED_DEV_GP_SHA256 = (
    "3bff0db14041b12b1af88ecd13b73ba09829c3abdde5ed9bd2b3112b368d7f88"
)
DEFAULT_DEV_GP = Path(
    "/private/tmp/pari-map-src.33iJSU/pari/Odarwin-aarch64/gp-dyn"
)
INPUT_ARTIFACT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_rank21_accidental_slices.json"
)
INPUT_SCRIPT_RELATIVE = Path(
    "elliptic-curves/cas/search_nagao_rank21_accidental_slices.py"
)
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_section7_auxiliary_group_orbits.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_auxiliary_group_orbits.py "
    "--gp /private/tmp/pari-map-src.33iJSU/pari/Odarwin-aarch64/gp-dyn"
)


@dataclass(frozen=True)
class OrbitSpecification:
    label: str
    priority_index_zero_based: int
    slope: Fraction
    intercept: Fraction
    expected_height_subset_indices: tuple[int, ...]
    expected_model: tuple[Fraction, ...]
    expected_saturated_basis: tuple[tuple[Fraction, Fraction], ...]


ORBIT_SPECIFICATIONS = (
    OrbitSpecification(
        "a04_sp01",
        3,
        Q(1),
        Q(-4471, 339),
        (1, 3, 5, 7, 9, 13, 17),
        (
            Q(0),
            Q(308236804344177057976),
            Q(0),
            Q(31283846573205158462859636953159685081792),
            Q(1044297476326303698701670888382203959163606498696263528114688),
        ),
        (
            (Q(-105124566564987195912), Q(-15014882519600127463709289600)),
            (Q(-67190859563465371272), Q(174750551360739203125733116800)),
            (Q(-109474269351777008392), Q(-40172350393279222905124867200)),
            (Q(-121431662267068836312), Q(3311608734892317642157435200)),
            (Q(-659495591773217649128, 9), Q(-3140772251486284527823078928000, 27)),
            (Q(-18924660459876049238728, 169), Q(-100644902505675467284433248137600, 2197)),
            (Q(1415934827874837544064952, 10609), Q(-3953144675474188410073986571363132800, 1092727)),
        ),
    ),
    OrbitSpecification(
        "a09_sm01",
        8,
        Q(-1),
        Q(154687, 447),
        (1, 3, 5, 7, 9, 11, 15),
        (
            Q(0),
            Q(674906638961445179536),
            Q(0),
            Q(150628140920708280631894264813465813098432),
            Q(11098749213681181177396216271409115760540451697617454196329728),
        ),
        (
            (Q(-245109617732461244712), Q(300466567130325345750858000)),
            (Q(-245089082889712195512), Q(559687633035660413294010000)),
            (Q(-244901591297101473112), Q(-36244550768908672179990000)),
            (Q(47869441627664759088), Q(4468273003491119701485911850000)),
            (Q(-245046004756840858312), Q(-784635144509410310964666000)),
            (Q(422610537838593770288), Q(-16455150890444258729225731242000)),
            (Q(-19837021240487970376172, 81), Q(-7033069897537959024212827000, 729)),
        ),
    ),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rational_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def gp_rational(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_vector(values: Sequence[Any]) -> str:
    parts = []
    for value in values:
        if isinstance(value, (tuple, list)):
            parts.append(gp_vector(value))
        else:
            parts.append(gp_rational(Q(value)))
    return "[" + ",".join(parts) + "]"


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def _eval_gp_ast(node: ast.AST) -> Any:
    if isinstance(node, ast.Expression):
        return _eval_gp_ast(node.body)
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return Q(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -Q(_eval_gp_ast(node.operand))
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return Q(_eval_gp_ast(node.left), _eval_gp_ast(node.right))
    if isinstance(node, (ast.List, ast.Tuple)):
        return tuple(_eval_gp_ast(item) for item in node.elts)
    raise ValueError(f"unsupported GP expression node: {ast.dump(node)}")


def parse_gp_exact(value: str) -> Any:
    return _eval_gp_ast(ast.parse(value.strip(), mode="eval"))


def polynomial_gp(coefficients: Sequence[Fraction], variable: str = "T") -> str:
    return "+".join(
        f"({gp_rational(Q(coefficient))})*{variable}^{power}"
        for power, coefficient in enumerate(coefficients)
    )


def load_inputs(root: Path) -> tuple[dict[str, Any], tuple[Any, ...]]:
    artifact_path = root / INPUT_ARTIFACT_RELATIVE
    input_script = root / INPUT_SCRIPT_RELATIVE
    if file_sha256(artifact_path) != EXPECTED_INPUT_SHA256:
        raise RuntimeError("the pinned accidental-slice artifact hash changed")
    if file_sha256(input_script) != EXPECTED_INPUT_SCRIPT_SHA256:
        raise RuntimeError("the pinned accidental-slice script hash changed")
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    accidental = tuple(
        (Q(record["x"]), Q(record["y"]))
        for record in artifact["decontamination_at_T0"]["accidental_points"]
    )
    priority = select_minimum_intercept_priority_slices(build_slices(accidental))
    if len(priority) != 16:
        raise AssertionError("the pinned priority-slice count changed")
    return artifact, priority


def h200000_parameters(
    artifact: dict[str, Any], priority_index: int
) -> tuple[Fraction, ...]:
    values = tuple(
        Q(record["T"])
        for record in artifact["auxiliary_search"]["pinned_identity_pilot"][
            "association_records"
        ]
        if priority_index
        in record["matching_priority_slice_indices_zero_based"]
    )
    return values + (T0,)


def build_gp_orbit_program(
    specification: OrbitSpecification,
    item: Any,
    known_parameters: Sequence[Fraction],
) -> str:
    if item.identifier != specification.label:
        raise AssertionError("the priority-slice identifier changed")
    if Q(item.slope) != specification.slope or Q(item.intercept) != specification.intercept:
        raise AssertionError("the priority-slice equation changed")
    coefficients = tuple(Q(value) for value in item.normalized.normalized_coefficients)
    point_expressions = []
    for parameter in known_parameters:
        ordinate = rational_square_root(item.normalized.normalized_value(parameter))
        if ordinate is None:
            raise AssertionError("a pinned H=200000 parameter lost its square")
        for signed_ordinate in (ordinate, -ordinate):
            point_expressions.append(
                "subst(subst(M,T,{T}),Y,{Y})".format(
                    T=gp_rational(parameter), Y=gp_rational(signed_ordinate)
                )
            )
    base_ordinate = rational_square_root(item.normalized.normalized_value(T0))
    if base_ordinate is None:
        raise AssertionError("the source point at T0 disappeared")
    expected_indices = gp_vector(specification.expected_height_subset_indices)
    expected_model = gp_vector(specification.expected_model)
    expected_basis = gp_vector(specification.expected_saturated_basis)
    return "\n".join(
        (
            "T='T;Y='Y;",
            f"F={polynomial_gp(coefficients)};",
            "Z=ellfromeqn(Y^2-F,[T,Y]);",
            "E=ellinit(Z[1]);M=Z[2];",
            f"if(Z[1]!={expected_model},error(\"MODEL_DRIFT\"));",
            f"V=[{','.join(point_expressions)}];",
            f"P0=subst(subst(M,T,{gp_rational(T0)}),Y,{gp_rational(base_ordinate)});",
            "if(!ellisoncurve(E,P0),error(\"BAD_P0\"));",
            'print("PARI_VERSION ",version());',
            'print("MODEL ",Z[1]);',
            'print("KNOWN_SIGNED_IMAGE_COUNT ",#V);',
            "default(realprecision,72)",
            "H72=ellheightmatrix(E,V);I72=matindexrank(H72)[2];",
            f"if(Vec(I72)!={expected_indices},error(\"HEIGHT72_SUBSET_DRIFT\"));",
            "K72=vecextract(V,I72);D72=matdet(ellheightmatrix(E,K72));",
            'print("HEIGHT72 ",matrank(H72),"|",Vec(I72),"|",D72);',
            "default(realprecision,120)",
            "H120=ellheightmatrix(E,V);I120=matindexrank(H120)[2];",
            f"if(Vec(I120)!={expected_indices},error(\"HEIGHT120_SUBSET_DRIFT\"));",
            "K=vecextract(V,I120);DK=matdet(ellheightmatrix(E,K));",
            'print("HEIGHT120 ",matrank(H120),"|",Vec(I120),"|",DK);',
            "W=ellsaturation(E,K,20);",
            f"if(W!={expected_basis},error(\"SATURATED_BASIS_DRIFT\"));",
            "DW=matdet(ellheightmatrix(E,W));SI=round(sqrt(DK/DW));",
            'print("SATURATION ",#W,"|",SI,"|",DW);',
            "A=(M[1]+Z[1][2]/4)*Y^2;",
            "if(poldegree(A,Y)!=0,error(\"MAP_X_REDUCTION_FAILED\"));",
            "L=List();C=0;",
            (
                "forvec(v=vector(#W,i,[-1,1]),{C++;Q0=[0];"
                "for(i=1,#W,if(v[i],Q0=elladd(E,Q0,ellmul(E,W[i],v[i]))));"
                "R=elladd(E,P0,ellmul(E,Q0,2));if(#R==2,"
                "G=A-(R[1]+Z[1][2]/4)*F;FA=factor(G);"
                "for(j=1,matsize(FA)[1],H0=FA[j,1];if(poldegree(H0)==1,"
                "listput(L,[-polcoeff(H0,0)/polcoeff(H0,1),Vec(v)]))))"
                "});"
            ),
            'print("COEFFICIENT_VECTOR_COUNT ",C);',
            'print("RAW_PULLBACK_COUNT ",#L);',
            'print("ORBIT ",Vec(L));',
            "quit",
            "",
        )
    )


def run_gp_orbit(
    gp_path: Path,
    program: str,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, tuple[int, ...]], ...]]:
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [str(gp_path), "-q", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"the one-shot GP orbit timed out after {timeout}s") from error
    wall_seconds = time.monotonic() - started
    if completed.returncode != 0 or "***" in completed.stderr:
        raise RuntimeError(
            "the one-shot GP orbit failed: " + completed.stderr.strip()[:1000]
        )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    by_prefix = {}
    for prefix in (
        "PARI_VERSION ",
        "MODEL ",
        "KNOWN_SIGNED_IMAGE_COUNT ",
        "HEIGHT72 ",
        "HEIGHT120 ",
        "SATURATION ",
        "COEFFICIENT_VECTOR_COUNT ",
        "RAW_PULLBACK_COUNT ",
        "ORBIT ",
    ):
        matches = [line[len(prefix) :] for line in lines if line.startswith(prefix)]
        if len(matches) != 1:
            raise AssertionError(f"GP omitted or duplicated marker {prefix!r}")
        by_prefix[prefix.strip()] = matches[0]
    orbit_raw = parse_gp_exact(by_prefix["ORBIT"])
    orbit = tuple(
        (Q(record[0]), tuple(int(value) for value in record[1]))
        for record in orbit_raw
    )
    for _, vector in orbit:
        if len(vector) != 7 or any(value not in (-1, 0, 1) for value in vector):
            raise AssertionError("GP emitted an out-of-budget coefficient vector")
    def parse_height(marker: str) -> dict[str, Any]:
        rank, indices, determinant = by_prefix[marker].split("|", 2)
        return {
            "decimal_precision": int(marker.removeprefix("HEIGHT")),
            "numerical_rank": int(rank),
            "subset_indices_one_based": [
                int(value) for value in parse_gp_exact(indices)
            ],
            "height_determinant": determinant,
        }
    saturation_rank, saturation_index, saturated_determinant = by_prefix[
        "SATURATION"
    ].split("|", 2)
    record = {
        "status": "completed",
        "wall_seconds": wall_seconds,
        "timeout_seconds": timeout,
        "stack_bytes": stack_bytes,
        "one_attempt_no_retry": True,
        "pari_version": by_prefix["PARI_VERSION"],
        "model": [rational_string(value) for value in parse_gp_exact(by_prefix["MODEL"])],
        "known_signed_image_count": int(by_prefix["KNOWN_SIGNED_IMAGE_COUNT"]),
        "height_replay": [parse_height("HEIGHT72"), parse_height("HEIGHT120")],
        "saturated_subgroup_rank": int(saturation_rank),
        "saturation_index_improvement": int(saturation_index),
        "saturated_height_determinant": saturated_determinant,
        "coefficient_vector_count": int(by_prefix["COEFFICIENT_VECTOR_COUNT"]),
        "raw_pullback_count": int(by_prefix["RAW_PULLBACK_COUNT"]),
    }
    return record, orbit


def generic_abscissa_labels(parameter: Fraction, x_value: Fraction) -> tuple[str, ...]:
    parameter = Q(parameter)
    x_value = Q(x_value)
    labels = []
    for index, root in enumerate(SECTION7_ROOTS):
        for sign, sign_label in ((1, "plus"), (-1, "minus")):
            if x_value == Q(root) + sign * parameter:
                labels.append(f"visible-{index:02d}-{sign_label}")
    for section in SECTION7_LINEAR_COMPANION_SECTIONS:
        if x_value == section.slope * parameter + section.intercept:
            labels.append(f"linear-{section.label}")
    for section in SECTION7_QUADRATIC_COMPANION_SECTIONS:
        candidate = (
            section.quadratic_coefficient * parameter**2
            + section.linear_coefficient * parameter
            + section.constant_coefficient
        )
        if x_value == candidate:
            labels.append(section.label)
    return tuple(labels)


def sibling_parameter_set(artifact: dict[str, Any]) -> set[Fraction]:
    values = {
        abs(Q(record["T"]))
        for record in artifact["auxiliary_search"]["pinned_identity_pilot"][
            "association_records"
        ]
    }
    values.add(abs(T0))
    return values


def validate_and_deduplicate_orbits(
    artifact: dict[str, Any],
    priority: Sequence[Any],
    raw_orbits: dict[str, tuple[tuple[Fraction, tuple[int, ...]], ...]],
) -> tuple[dict[Fraction, list[dict[str, Any]]], dict[str, int]]:
    sibling = sibling_parameter_set(artifact)
    candidates: dict[Fraction, list[dict[str, Any]]] = {}
    counts = {
        "raw_pullbacks": 0,
        "nonsquare_pullbacks": 0,
        "zero_parameters": 0,
        "sibling_H200000_or_T0_exclusions": 0,
        "projective_height_at_most_200000_exclusions": 0,
        "generic_section_intersection_exclusions": 0,
        "singular_parameter_exclusions": 0,
    }
    specifications = {item.label: item for item in ORBIT_SPECIFICATIONS}
    for label, orbit in raw_orbits.items():
        specification = specifications[label]
        item = priority[specification.priority_index_zero_based]
        for signed_parameter, vector in orbit:
            counts["raw_pullbacks"] += 1
            normalized_ordinate = rational_square_root(
                item.normalized.normalized_value(signed_parameter)
            )
            if normalized_ordinate is None:
                counts["nonsquare_pullbacks"] += 1
                continue
            parameter = abs(signed_parameter)
            if parameter == 0:
                counts["zero_parameters"] += 1
                continue
            if parameter in sibling:
                counts["sibling_H200000_or_T0_exclusions"] += 1
                continue
            if projective_height(parameter) <= NAIVE_HEIGHT_EXCLUSION:
                counts["projective_height_at_most_200000_exclusions"] += 1
                continue
            x_value = specification.slope * signed_parameter + specification.intercept
            removed = (
                item.normalized.removed_square_value(signed_parameter)
                * item.normalized.ordinate_constant_scale
            )
            original_ordinate = normalized_ordinate * removed
            primitive = section7_primitive_quartic_coefficients(parameter)
            if original_ordinate**2 != quartic_value(primitive, x_value):
                raise AssertionError("an auxiliary pullback missed the original quartic")
            generic_labels = generic_abscissa_labels(parameter, x_value)
            if generic_labels:
                counts["generic_section_intersection_exclusions"] += 1
                continue
            try:
                homogenized_discriminant(parameter)
            except ValueError:
                counts["singular_parameter_exclusions"] += 1
                continue
            source = {
                "slice": label,
                "signed_constructor_parameter": rational_string(signed_parameter),
                "coefficient_vector": list(vector),
                "forced_quartic_x_on_canonical_fiber": rational_string(x_value),
                "forced_quartic_y_one_sign": rational_string(original_ordinate),
                "exact_original_quartic_membership_checked": True,
                "all_21_generic_abscissas_checked": True,
            }
            if source not in candidates.setdefault(parameter, []):
                candidates[parameter].append(source)
    counts["unique_parameters_after_all_exact_exclusions"] = len(candidates)
    return candidates, counts


def stream_sha256(parameters: Iterable[Fraction]) -> str:
    text = "\n".join(
        rational_string(value)
        for value in sorted(
            parameters,
            key=lambda value: (projective_height(value), value.numerator, value.denominator),
        )
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def proxy_population(
    candidates: dict[Fraction, list[dict[str, Any]]],
) -> tuple[dict[str, Any], ...]:
    records = []
    for index, (parameter, sources) in enumerate(candidates.items(), start=1):
        proxy = conductor_radical_proxy(
            parameter, trial_prime_bound=PROXY_TRIAL_BOUND
        )
        records.append(
            {
                "constructor_parameter": rational_string(parameter),
                "projective_height": projective_height(parameter),
                "radical_proxy": proxy,
                "source_count": len(sources),
                "first_sources": sources[:4],
            }
        )
        if index % 500 == 0:
            print(f"proxy {index}/{len(candidates)}", flush=True)
    return tuple(
        sorted(
            records,
            key=lambda record: (
                record["radical_proxy"]["log_radical_upper_proxy"],
                record["projective_height"],
                Q(record["constructor_parameter"]),
            ),
        )
    )


def exact_conductors(
    records: Sequence[dict[str, Any]],
    *,
    keep: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    completed_records = []
    failures = []
    selected = records[: min(keep, len(records))]
    for index, record in enumerate(selected, start=1):
        parameter = Q(record["constructor_parameter"])
        coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
        try:
            data = minimal_curve_data(
                coefficients, timeout=timeout, stack_bytes=stack_bytes
            )
        except (subprocess.TimeoutExpired, RuntimeError) as error:
            failure = {
                "constructor_parameter": rational_string(parameter),
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "one_attempt_no_retry": True,
                "error": str(error)[:500],
            }
            failures.append(failure)
            print(
                f"conductor {index}/{len(selected)} T={parameter} "
                f"status={failure['status']}",
                flush=True,
            )
            continue
        exact = {
            **record,
            "status": "completed",
            "conductor": str(data["conductor"]),
            "log_conductor": data["log_conductor"],
            "root_number": data["root_number"],
            "minimal_discriminant": str(data["minimal_discriminant"]),
            "minimal_model": [str(value) for value in data["minimal_model"]],
            "below_strict_log_conductor_target": (
                Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
        }
        completed_records.append(exact)
        print(
            f"conductor {index}/{len(selected)} T={parameter} "
            f"lnN={data['log_conductor']} "
            f"subtarget={exact['below_strict_log_conductor_target']}",
            flush=True,
        )
    return tuple(completed_records), tuple(failures)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", type=Path, default=DEFAULT_DEV_GP)
    parser.add_argument("--orbit-timeout", type=float, default=60.0)
    parser.add_argument("--gp-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--conductor-keep", type=int, default=DEFAULT_CONDUCTOR_KEEP)
    parser.add_argument("--skip-conductors", action="store_true")
    parser.add_argument("--output", type=Path, default=root / OUTPUT_RELATIVE)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.orbit_timeout <= 60:
        raise SystemExit("--orbit-timeout must be in (0,60]")
    if not 0 < args.conductor_timeout <= 60:
        raise SystemExit("--conductor-timeout must be in (0,60]")
    if not 1 <= args.conductor_keep <= 128:
        raise SystemExit("--conductor-keep must be in [1,128]")
    if args.gp_stack_bytes > 1_000_000_000 or args.gp_stack_bytes < 256_000_000:
        raise SystemExit("--gp-stack-bytes must be in [256MB,1GB]")
    if not args.gp.is_file():
        raise SystemExit(f"development GP binary not found: {args.gp}")
    if file_sha256(args.gp) != EXPECTED_DEV_GP_SHA256:
        raise SystemExit("the development GP binary hash changed")

    root = Path(__file__).resolve().parents[2]
    artifact_input, priority = load_inputs(root)
    raw_orbits = {}
    orbit_records = []
    for specification in ORBIT_SPECIFICATIONS:
        item = priority[specification.priority_index_zero_based]
        known = h200000_parameters(
            artifact_input, specification.priority_index_zero_based
        )
        program = build_gp_orbit_program(specification, item, known)
        orbit_record, orbit = run_gp_orbit(
            args.gp,
            program,
            timeout=args.orbit_timeout,
            stack_bytes=args.gp_stack_bytes,
        )
        orbit_record.update(
            {
                "slice": specification.label,
                "slope": rational_string(specification.slope),
                "intercept": rational_string(specification.intercept),
                "known_signless_parameter_count_including_T0": len(known),
                "stable_numerical_auxiliary_rank": 7,
                "rank_status": (
                    "two-precision numerical height rank; not an exact rank upper bound"
                ),
            }
        )
        raw_orbits[specification.label] = orbit
        orbit_records.append(orbit_record)
        print(
            f"orbit {specification.label}: vectors={orbit_record['coefficient_vector_count']} "
            f"pullbacks={orbit_record['raw_pullback_count']} "
            f"rank={orbit_record['saturated_subgroup_rank']} "
            f"sat_index={orbit_record['saturation_index_improvement']}",
            flush=True,
        )

    candidates, exclusion_counts = validate_and_deduplicate_orbits(
        artifact_input, priority, raw_orbits
    )
    proxies = proxy_population(candidates)
    proxy_under_gate = sum(
        record["radical_proxy"]["log_radical_upper_proxy"] < PROXY_GATE
        for record in proxies
    )
    print(
        "GENERATION_CHECKPOINT "
        + json.dumps(
            {
                **exclusion_counts,
                "candidate_stream_sha256": stream_sha256(candidates),
                "proxy_population_count": len(proxies),
                "proxy_below_190_count": proxy_under_gate,
                "minimum_proxy": (
                    proxies[0]["radical_proxy"]["log_radical_upper_proxy"]
                    if proxies
                    else None
                ),
            },
            sort_keys=True,
        ),
        flush=True,
    )

    if args.skip_conductors:
        conductor_records: tuple[dict[str, Any], ...] = ()
        conductor_failures: tuple[dict[str, Any], ...] = ()
    else:
        conductor_records, conductor_failures = exact_conductors(
            proxies,
            keep=args.conductor_keep,
            timeout=args.conductor_timeout,
            stack_bytes=args.conductor_stack_bytes,
        )
    subtarget = tuple(
        record
        for record in conductor_records
        if record["below_strict_log_conductor_target"]
    )
    generation = {
        **exclusion_counts,
        "candidate_stream_sha256": stream_sha256(candidates),
        "proxy_population_count": len(proxies),
        "proxy_below_190_count": proxy_under_gate,
    }
    output = {
        "schema_version": 1,
        "status": "bounded_exact_auxiliary_group_orbit_complete",
        "claim_scope": {
            "exact": (
                "ellfromeqn models/maps, B=20 saturation replay, all 2*3^7 "
                "coefficient vectors, rational pullbacks, quartic membership, "
                "generic-section exclusions, radical proxies and completed conductors"
            ),
            "numerical": "stable auxiliary height ranks at 72 and 120 digits",
            "bounded": (
                "only the two declared open high-yield slices and coefficients "
                "in {-1,0,1}; a negative result is not a rank upper bound"
            ),
        },
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "inputs": {
            "accidental_slice_artifact": str(INPUT_ARTIFACT_RELATIVE),
            "accidental_slice_artifact_sha256": EXPECTED_INPUT_SHA256,
            "accidental_slice_script": str(INPUT_SCRIPT_RELATIVE),
            "accidental_slice_script_sha256": EXPECTED_INPUT_SCRIPT_SHA256,
            "excluded_closed_slice_equations": [
                ["1", "-12930/91"],
                ["-1", "130932/491"],
                ["1", "6380/47"],
            ],
            "coordination_scope": (
                "a04_sp01 and a09_sm01 only; the other eleven open priority "
                "slices were delegated to a separate nonoverlapping lane"
            ),
        },
        "method": {
            "map": (
                "PARI ellfromeqn quartic-to-Jacobian map; solve exact x-map "
                "equation for every P0+2Q"
            ),
            "coefficient_alphabet": [-1, 0, 1],
            "coefficient_dimension": 7,
            "coefficient_vectors_per_slice": 3**7,
            "saturation_prime_bound": 20,
            "naive_projective_height_exclusion": NAIVE_HEIGHT_EXCLUSION,
            "all_21_generic_sections_decontaminated": True,
            "parameter_even_symmetry_canonicalization": "T -> |T|",
            "proxy_trial_prime_bound": PROXY_TRIAL_BOUND,
            "proxy_gate": PROXY_GATE,
        },
        "orbit_records": orbit_records,
        "generation": generation,
        "proxy_selection": {
            "records_stored": min(128, len(proxies)),
            "top_records": list(proxies[:128]),
            "exact_conductor_keep": 0 if args.skip_conductors else min(args.conductor_keep, len(proxies)),
        },
        "exact_conductors": {
            "skipped": args.skip_conductors,
            "attempted": 0 if args.skip_conductors else min(args.conductor_keep, len(proxies)),
            "completed": len(conductor_records),
            "failures": list(conductor_failures),
            "records": list(conductor_records),
            "subtarget_count": len(subtarget),
            "subtarget_records": list(subtarget),
        },
        "outcome": {
            "subtarget_parameters_for_rank_triage": [
                record["constructor_parameter"] for record in subtarget
            ],
            "rank21_certified_in_this_lane": False,
        },
        "bounded_process_policy": {
            "orbit_timeout_seconds": args.orbit_timeout,
            "orbit_stack_bytes": args.gp_stack_bytes,
            "conductor_timeout_seconds": args.conductor_timeout,
            "conductor_stack_bytes": args.conductor_stack_bytes,
            "one_attempt_no_retry_per_declared_call": True,
            "foreground_only": True,
            "no_detached_processes": True,
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "script_sha256": file_sha256(Path(__file__).resolve()),
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "development_gp_path": str(args.gp),
            "development_gp_sha256": EXPECTED_DEV_GP_SHA256,
            "pari_version_from_orbits": sorted(
                {record["pari_version"] for record in orbit_records}
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {args.output}: candidates={len(candidates)} "
        f"proxy<190={proxy_under_gate} conductors={len(conductor_records)} "
        f"subtarget={len(subtarget)}",
        flush=True,
    )


if __name__ == "__main__":
    main()

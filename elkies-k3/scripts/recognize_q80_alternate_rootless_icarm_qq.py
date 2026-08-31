#!/usr/bin/env python3
"""Test a reconstructed Q80 rootless family against the selected ICARM curves.

The rank-29 curve is always tested first.  Curves 398--400 follow, then 273
and 302.  For each target this script factors the exact projective j-matching
equation over Q, checks every rational point of P^1(Q), and distinguishes a
mere equal-j quadratic twist from a Q-isomorphic specialization.

The default model path is the intended output of the generic Q80 equation
reconstruction and therefore need not exist before that reconstruction is
complete.  ``--model`` permits fail-closed regression tests on another exact
short Weierstrass family.  No target is declared recognized from a modular
residue or from equal j alone.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

import sympy


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = (
    ROOT / "artifacts/generated-results/q80-alternate-rootless-mw17-qq.json"
)
DEFAULT_TARGETS = (
    ROOT / "elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "q80-alternate-rootless-icarm-recognition.json"
)
CORE_PATH = (
    ROOT
    / "elkies-k3/scripts/recognize_h92_q12o5867_icarm_specializations_qq.py"
)
TARGET_ORDER = ("rank29", "curve398", "curve399", "curve400", "curve273", "curve302")


def load_core():
    specification = importlib.util.spec_from_file_location(
        "elkies_k3_exact_j_recognition_core", CORE_PATH
    )
    if specification is None or specification.loader is None:
        raise ImportError(f"cannot import exact recognition core from {CORE_PATH}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


core = load_core()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--require-rank29",
        action="store_true",
        help="exit unsuccessfully unless a Q-isomorphic rank-29 specialization exists",
    )
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_path(path: Path) -> str:
    """Render either a relative or absolute in-repository path canonically."""
    return str(path.resolve().relative_to(ROOT))


def get_path(payload: dict[str, object], path: tuple[str, ...]):
    value: object = payload
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def extract_short_model(payload: dict[str, object]):
    """Read one unambiguous short Weierstrass coefficient pair.

    The first layout is shared with the exact H3 rootless model and is the
    preferred contract for the future Q80 artifact.  The remaining layouts
    allow a compact standalone rootless-family certificate.
    """
    layouts = (
        (
            ("child", "minimal_A_coefficients_low_to_high"),
            ("child", "minimal_B_coefficients_low_to_high"),
        ),
        (
            ("minimal_A_coefficients_low_to_high",),
            ("minimal_B_coefficients_low_to_high",),
        ),
        (
            ("rootless_family", "A_coefficients_low_to_high"),
            ("rootless_family", "B_coefficients_low_to_high"),
        ),
        (
            ("A_coefficients_low_to_high",),
            ("B_coefficients_low_to_high",),
        ),
    )
    matches = []
    for a_path, b_path in layouts:
        a_values = get_path(payload, a_path)
        b_values = get_path(payload, b_path)
        if a_values is not None and b_values is not None:
            matches.append((a_path, b_path, a_values, b_values))
    if len(matches) != 1:
        raise ValueError(
            f"expected one supported short-Weierstrass layout, found {len(matches)}"
        )
    a_path, b_path, a_values, b_values = matches[0]
    if not isinstance(a_values, list) or not isinstance(b_values, list):
        raise TypeError("short-Weierstrass coefficients must be JSON lists")
    A = core.trim(list(map(Q, a_values)))
    B = core.trim(list(map(Q, b_values)))
    if len(A) > 9 or len(B) > 13:
        raise ValueError("a K3 short model must have deg(A)<=8 and deg(B)<=12")
    if A == [0] or B == [0]:
        raise ValueError("this recognizer requires nonzero A and B")
    return A, B, ".".join(a_path), ".".join(b_path)


def target_table(payload: dict[str, object]) -> list[dict[str, object]]:
    rows = payload.get("targets")
    if not isinstance(rows, list):
        raise ValueError("recognition target payload has no target list")
    by_label = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("label"), str):
            raise ValueError("malformed recognition target row")
        by_label[row["label"]] = row
    missing = [label for label in TARGET_ORDER if label not in by_label]
    if missing:
        raise ValueError(f"missing required recognition targets: {missing}")
    return [by_label[label] for label in TARGET_ORDER]


def main() -> None:
    arguments = parse_args()
    started = time.monotonic()
    if not arguments.model.exists():
        raise SystemExit(f"missing reconstructed Q80 model: {arguments.model}")
    model_raw = arguments.model.read_bytes()
    targets_raw = arguments.targets.read_bytes()
    model = json.loads(model_raw)
    target_payload = json.loads(targets_raw)
    A, B, a_path, b_path = extract_short_model(model)
    targets = target_table(target_payload)

    c4 = core.poly_scale(A, -48)
    c6 = core.poly_scale(B, -864)
    discriminant = core.poly_scale(
        core.poly_add(
            core.poly_scale(core.poly_pow(A, 3), 4),
            core.poly_scale(core.poly_pow(B, 2), 27),
        ),
        -16,
    )
    if discriminant == [0]:
        raise ArithmeticError("the supplied family is identically singular")
    variable = sympy.symbols("u")
    common = core.from_sympy(
        core.to_sympy(core.poly_pow(c4, 3), variable).gcd(
            core.to_sympy(discriminant, variable)
        )
    )
    j_numerator = core.poly_div_exact(core.poly_pow(c4, 3), common)
    j_denominator = core.poly_div_exact(discriminant, common)
    j_num_ZZ, j_den_ZZ = core.primitive_integer_pair(j_numerator, j_denominator)

    target_records = []
    for target_index, target in enumerate(targets):
        label = str(target["label"])
        ainvs = tuple(map(Q, target["ainvs"]))
        invariants = core.generalized_invariants(ainvs)
        if invariants["discriminant"] == 0:
            raise ArithmeticError(f"target {label} is singular")
        j_target = invariants["j"]
        matching = core.poly_add(
            core.poly_scale(core.poly_pow(c4, 3), j_target.denominator),
            core.poly_scale(discriminant, -j_target.numerator),
        )
        if matching == [0]:
            raise ArithmeticError(f"family is isotrivial with target j for {label}")
        primitive = core.primitive_integer_polynomial(matching)
        primitive_sympy = sympy.Poly.from_list(
            list(reversed(primitive)), gens=variable, domain="ZZ"
        )
        factor_unit, factorization = sympy.factor_list(primitive_sympy)
        product_check = sympy.Poly(factor_unit, variable, domain="ZZ")
        factor_records = []
        finite_roots: list[tuple[Fraction, int]] = []
        for factor, multiplicity in factorization:
            factor_poly = core.from_sympy(factor)
            product_check *= factor**multiplicity
            factor_record = core.polynomial_record(factor_poly)
            factor_record["multiplicity"] = int(multiplicity)
            if len(factor_poly) == 2:
                root = -factor_poly[0] / factor_poly[1]
                finite_roots.append((root, int(multiplicity)))
                factor_record["rational_root"] = core.rational_text(root)
            factor_records.append(factor_record)
        if product_check != primitive_sympy:
            raise ArithmeticError("factorization product check failed")
        infinity_multiplicity = 24 - (len(primitive) - 1)
        if infinity_multiplicity < 0:
            raise ArithmeticError("j-matching equation exceeds projective degree 24")
        roots = [
            (root.numerator, root.denominator, multiplicity)
            for root, multiplicity in finite_roots
        ]
        if infinity_multiplicity:
            roots.append((1, 0, infinity_multiplicity))

        solutions = []
        for numerator, denominator, multiplicity in roots:
            specialized_A = core.projective_value(A, 8, numerator, denominator)
            specialized_B = core.projective_value(B, 12, numerator, denominator)
            specialized_delta = -16 * (
                4 * specialized_A**3 + 27 * specialized_B**2
            )
            solution: dict[str, object] = {
                "parameter_projective": [str(numerator), str(denominator)],
                "parameter": (
                    "infinity"
                    if denominator == 0
                    else core.rational_text(Q(numerator, denominator))
                ),
                "root_multiplicity": multiplicity,
                "specialized_A": core.rational_text(specialized_A),
                "specialized_B": core.rational_text(specialized_B),
                "specialized_discriminant": core.rational_text(specialized_delta),
                "nonsingular": specialized_delta != 0,
                "equal_j": False,
                "q_isomorphic": False,
            }
            if specialized_delta:
                specialized_j = (-48 * specialized_A) ** 3 / specialized_delta
                if specialized_j != j_target:
                    raise ArithmeticError("rational root failed exact j substitution")
                solution["equal_j"] = True
                target_A = invariants["short_A"]
                target_B = invariants["short_B"]
                if not all((specialized_A, specialized_B, target_A, target_B)):
                    raise ArithmeticError("j=0 or 1728 twist case is not implemented")
                twist_parameter = (target_B / specialized_B) / (
                    target_A / specialized_A
                )
                if target_A != twist_parameter**2 * specialized_A:
                    raise ArithmeticError("A twist identity failed")
                if target_B != twist_parameter**3 * specialized_B:
                    raise ArithmeticError("B twist identity failed")
                scale = core.rational_square_root(twist_parameter)
                solution["twist_parameter_in_Qmod_squares"] = core.rational_text(
                    twist_parameter
                )
                solution["q_isomorphic"] = scale is not None
                if scale is not None:
                    solution["isomorphism_family_to_target_short"] = {
                        "x_target_over_x_family": core.rational_text(scale**2),
                        "y_target_over_y_family": core.rational_text(scale**3),
                    }
            solutions.append(solution)

        target_records.append(
            {
                "test_index": target_index,
                "label": label,
                "curve_id": int(target["icarm_id"]),
                "rank_lower_bound": int(target["certified_rank_lower_bound"]),
                "source": target["source"],
                "source_sha256": target["source_sha256"],
                "ainvs": list(map(core.rational_text, ainvs)),
                "j": core.rational_text(j_target),
                "primitive_affine_matching_polynomial": core.polynomial_record(
                    list(map(Q, primitive))
                ),
                "factorization_over_Q": factor_records,
                "irreducible_degree_24_over_Q": (
                    len(factorization) == 1
                    and factorization[0][0].degree() == 24
                    and factorization[0][1] == 1
                ),
                "infinity_root_multiplicity": infinity_multiplicity,
                "rational_projective_solution_count": len(solutions),
                "q_isomorphic_solution_count": sum(
                    bool(solution["q_isomorphic"]) for solution in solutions
                ),
                "rational_projective_solutions": solutions,
            }
        )

    if [row["label"] for row in target_records] != list(TARGET_ORDER):
        raise ArithmeticError("target order changed unexpectedly")
    rank29_recognized = target_records[0]["q_isomorphic_solution_count"] > 0
    recognized = [
        row["label"]
        for row in target_records
        if row["q_isomorphic_solution_count"] > 0
    ]
    status = (
        "PASS_EXACT_Q80_RANK29_POSITIVE_CONTROL"
        if rank29_recognized
        else "PASS_EXACT_Q80_TARGET_TESTS_WITHOUT_RANK29_CONTROL"
    )
    payload = {
        "schema": "elkies-k3.q80-alternate-rootless-icarm-recognition.v1",
        "status": status,
        "inputs": {
            repository_path(arguments.model): hashlib.sha256(model_raw).hexdigest(),
            repository_path(arguments.targets): hashlib.sha256(targets_raw).hexdigest(),
            str(CORE_PATH.relative_to(ROOT)): sha256(CORE_PATH),
        },
        "software": {"python": sys.version.split()[0], "sympy": sympy.__version__},
        "rootless_family": {
            "equation": "y^2=x^3+A(u)*x+B(u)",
            "A_json_path": a_path,
            "B_json_path": b_path,
            "A": core.polynomial_record(A),
            "B": core.polynomial_record(B),
            "reduced_j": {
                "formula": "c4(u)^3/Delta(u)",
                "cancelled_gcd_degree": len(common) - 1,
                "primitive_integer_numerator": core.polynomial_record(
                    list(map(Q, j_num_ZZ)), True
                ),
                "primitive_integer_denominator": core.polynomial_record(
                    list(map(Q, j_den_ZZ)), True
                ),
            },
        },
        "mandatory_test_order": list(TARGET_ORDER),
        "rank29_positive_control": {
            "tested_first": True,
            "recognized_over_Q": rank29_recognized,
            "mandatory_for_rank32_sieve_if_recognized": True,
        },
        "recognized_targets_over_Q": recognized,
        "targets": target_records,
        "proof_boundary": (
            "Exact projective j-recognition, nonsingularity, and Q-twist class "
            "for every rational solution. A target is recognized only when the "
            "specialization is Q-isomorphic, not merely equal-j. No Mordell--Weil "
            "rank or section independence is inferred by this script."
        ),
        "reproducing_command": (
            f".venv/bin/python {Path(__file__).relative_to(ROOT)} "
            f"--model {repository_path(arguments.model)}"
        ),
        "runtime_seconds": time.monotonic() - started,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.exists():
            raise SystemExit(f"missing recognition artifact: {arguments.output}")
        existing = json.loads(arguments.output.read_text())
        existing.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if existing != payload:
            raise SystemExit("stale Q80 ICARM recognition artifact")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)

    print(
        "Q80ICARM|"
        + "|".join(
            f"{row['label']}={row['q_isomorphic_solution_count']}"
            for row in target_records
        )
        + f"|rank29_control={int(rank29_recognized)}|status={status}",
        flush=True,
    )
    if arguments.require_rank29 and not rank29_recognized:
        raise SystemExit("rank-29 positive control was not recognized over Q")


if __name__ == "__main__":
    main()

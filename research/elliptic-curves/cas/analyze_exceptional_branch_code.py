#!/usr/bin/env python3
"""Build exact branch-divisor codes for available exceptional directions.

The Fermigier E22 and rank-20 anchors, plus ICARM 245 and 275, have pinned
exceptional quartic abscissas.  Holding one such abscissa fixed gives an exact
condition ``z^2=g(T)`` on the corresponding one-parameter Mestre family.
This script factors all of those divisors modulo two, enumerates every
elementary quadratic quotient, and tests code-level branch identifications.

The generic engine in :mod:`branch_divisor_code` is intentionally independent
of this data loader.  Future certified records (including ICARM 262 and K3
bisections) can be supplied as exact numerator/denominator coefficient lists
without changing its genus or cancellation calculations.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
from typing import Any, Sequence

import sympy as sp

from analyze_fermigier_exceptional_transport import (
    E22_INDEPENDENT_EXCEPTIONAL_LABELS,
    family_expression,
    primitive_poly,
    e22_exceptional_quotient,
    rank20_exceptional_quotient,
)
from branch_divisor_code import (
    SquareCondition,
    analyze_square_conditions,
    conditions_from_json_records,
)
from verify_icarm_curve275_mestre_rank20 import (
    EXPECTED_EXCEPTIONAL_X as ICARM275_EXCEPTIONAL_X,
    ROOTS as ICARM275_ROOTS,
)


Q = Fraction
ICARM245_ROOTS = tuple(map(Q, (0, 106, 344, 475, 594, 731)))
ICARM245_EXCEPTIONAL_X = tuple(
    map(
        Q,
        (
            "-1069/530", "1107/2", "-1217/10", "-2901/10",
            "3085/46", "-5773/2", "9517/170", "9933/10",
        ),
    )
)
ICARM243_ROOTS = tuple(
    map(Q, ("-1455/4", "2955/4", "1437/2", "-1149/4", "-1851/4", "-687/2"))
)
ICARM243_ANCHOR = Q(3895, 6)
ICARM243_BOUNDED_ACCIDENTAL_X = tuple(map(Q, (
    "-199739/204", "-160446/167", "138283/1788", "102053/1077",
    "93189/836", "17485/116", "19376/93", "195209/708",
    "108959/372", "165249/547", "106276/285", "145273/372",
    "30215/59", "115053/157", "2761/3", "49127/36",
)))
ICARM226_ROOTS = tuple(map(Q, (-138, -90, -60, -12, 138, 162)))
ICARM226_ANCHOR = Q(10167, 350)
ICARM226_BOUNDED_ACCIDENTAL_X = tuple(map(Q, (
    "-177663/350", "-42081/350", "-21169/210", "-4797/70",
    "-2663/50", "-20091/1750", "-78971/7000", "-1211/150",
    "-37671/30905", "46143/22750", "1581/700", "59937/21350",
    "14661/1300", "2487/175", "102279/3850", "14239/525",
    "33497/350", "70661/525", "69081/350",
)))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def result_digest(payload: dict[str, Any]) -> str:
    stable = dict(payload)
    stable.pop("result_sha256", None)
    stable.pop("generated_at_utc", None)
    return hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def mestre_quartic_expression(
    roots: Sequence[Fraction], parameter: sp.Symbol, abscissa: sp.Symbol
) -> sp.Expr:
    """Reconstruct ``(g^2-q(X-T)q(X+T))/T^2`` over ``QQ(T,X)``."""

    product = sp.prod(
        (abscissa - (sp.Rational(root.numerator, root.denominator) + sign * parameter))
        for root in roots
        for sign in (-1, 1)
    )
    product = sp.Poly(sp.expand(product), abscissa, domain=sp.QQ.frac_field(parameter))
    approximant = [sp.Integer(0)] * 7
    approximant[6] = sp.Integer(1)
    for index in range(5, -1, -1):
        square = sp.Poly(
            sum(approximant[degree] * abscissa**degree for degree in range(7)) ** 2,
            abscissa,
            domain=sp.QQ.frac_field(parameter),
        )
        approximant[index] = sp.cancel((product.nth(6 + index) - square.nth(6 + index)) / 2)
    remainder = sp.Poly(
        sum(approximant[degree] * abscissa**degree for degree in range(7)) ** 2
        - product.as_expr(),
        abscissa,
        domain=sp.QQ.frac_field(parameter),
    )
    if any(remainder.nth(index) != 0 for index in range(5, 13)):
        raise AssertionError("the supplied six-root tuple is not a Mestre quartic family")
    return sp.cancel(sum(remainder.nth(index) * abscissa**index for index in range(5)) / parameter**2)


def conditions_from_expression(
    labels_and_x: Sequence[tuple[str, Fraction]], expression: sp.Expr, parameter: sp.Symbol, x: sp.Symbol
) -> list[SquareCondition]:
    conditions = []
    for label, abscissa in labels_and_x:
        value = sp.cancel(expression.subs(x, sp.Rational(abscissa.numerator, abscissa.denominator)))
        numerator, denominator = sp.fraction(value)
        conditions.append(
            SquareCondition(
                label,
                primitive_poly(numerator, parameter),
                primitive_poly(denominator, parameter),
            )
        )
    return conditions


def _is_rational_square(value: sp.Rational) -> bool:
    """Test a nonnegative rational square without floating arithmetic."""

    return (
        value >= 0
        and isqrt(int(value.p)) ** 2 == int(value.p)
        and isqrt(int(value.q)) ** 2 == int(value.q)
    )


def bounded_observed_conditions(
    labels_and_x: Sequence[tuple[str, Fraction]],
    expression: sp.Expr,
    parameter: sp.Symbol,
    x: sp.Symbol,
    anchor: Fraction,
) -> list[SquareCondition]:
    """Decode a finite observed preimage list and recheck it at its anchor.

    The resulting cover conditions are exact.  The caller remains responsible
    for describing the source list as bounded whenever it came from a bounded
    rational-point search; this function makes no completeness or rank claim.
    """

    conditions = conditions_from_expression(labels_and_x, expression, parameter, x)
    anchor_value = sp.Rational(anchor.numerator, anchor.denominator)
    for condition in conditions:
        value = sp.cancel(
            condition.numerator.as_expr().subs(parameter, anchor_value)
            / condition.denominator.as_expr().subs(parameter, anchor_value)
        )
        if not _is_rational_square(sp.Rational(value)):
            raise AssertionError(
                f"{condition.label} is not a rational square at the declared anchor"
            )
    return conditions


def _external_families(paths: Sequence[Path]) -> tuple[dict[str, Any], dict[str, str]]:
    """Load exact future-condition records without weakening the pinned run."""

    families: dict[str, Any] = {}
    hashes: dict[str, str] = {}
    for path in paths:
        payload = json.loads(path.read_text())
        entries = payload.get("families")
        if not isinstance(entries, list):
            raise ValueError(f"{path}: expected a top-level families list")
        hashes[str(path)] = sha256_file(path)
        for entry in entries:
            if not isinstance(entry, dict):
                raise ValueError(f"{path}: family entry must be an object")
            name = entry.get("name")
            records = entry.get("conditions")
            if not isinstance(name, str) or not name:
                raise ValueError(f"{path}: family name must be a nonempty string")
            if name in families:
                raise ValueError(f"duplicate external branch-code family {name}")
            if not isinstance(records, list):
                raise ValueError(f"{path}: {name} has no conditions list")
            families[name] = analyze_square_conditions(
                conditions_from_json_records(
                    records, parameter_name=str(entry.get("parameter", "u"))
                ),
                maximum_records=0,
            )
    return families, hashes


def run(root: Path, *, extra_condition_files: Sequence[Path] = ()) -> dict[str, Any]:
    parameter, x = sp.symbols("T X")
    _, e22_all = e22_exceptional_quotient()
    _, rank20 = rank20_exceptional_quotient(root)
    e22 = [item for item in e22_all if item[0] in E22_INDEPENDENT_EXCEPTIONAL_LABELS]
    fermigier_expression = family_expression(parameter, x)
    icarm245_expression = mestre_quartic_expression(ICARM245_ROOTS, parameter, x)
    icarm275_expression = mestre_quartic_expression(ICARM275_ROOTS, parameter, x)
    icarm243_expression = mestre_quartic_expression(ICARM243_ROOTS, parameter, x)
    icarm226_expression = mestre_quartic_expression(ICARM226_ROOTS, parameter, x)

    e22_conditions = conditions_from_expression(
        e22, fermigier_expression, parameter, x
    )
    rank20_conditions = conditions_from_expression(
        rank20, fermigier_expression, parameter, x
    )
    families = {
        "fermigier_E22": analyze_square_conditions(
            e22_conditions,
            maximum_records=0,
        ),
        "fermigier_rank20_anchor": analyze_square_conditions(
            rank20_conditions,
            maximum_records=0,
        ),
        "fermigier_E22_rank20_combined": analyze_square_conditions(
            [*e22_conditions, *rank20_conditions],
            maximum_records=0,
        ),
        "icarm_245": analyze_square_conditions(
            conditions_from_expression(
                [(f"IC245E{index}", value) for index, value in enumerate(ICARM245_EXCEPTIONAL_X, 1)],
                icarm245_expression,
                parameter,
                x,
            ),
            maximum_records=0,
        ),
        "icarm_275": analyze_square_conditions(
            conditions_from_expression(
                [(f"IC275E{index}", value) for index, value in enumerate(ICARM275_EXCEPTIONAL_X, 1)],
                icarm275_expression,
                parameter,
                x,
            ),
            maximum_records=0,
        ),
    }
    supplemental_bounded_families = {
        "icarm_243_observed_preimages_h200000": analyze_square_conditions(
            bounded_observed_conditions(
                [(f"IC243O{index}", value) for index, value in enumerate(ICARM243_BOUNDED_ACCIDENTAL_X, 1)],
                icarm243_expression,
                parameter,
                x,
                ICARM243_ANCHOR,
            ),
            maximum_records=0,
        ),
        "icarm_226_observed_preimages_h200000": analyze_square_conditions(
            bounded_observed_conditions(
                [(f"IC226O{index}", value) for index, value in enumerate(ICARM226_BOUNDED_ACCIDENTAL_X, 1)],
                icarm226_expression,
                parameter,
                x,
                ICARM226_ANCHOR,
            ),
            maximum_records=0,
        ),
    }
    external_families, external_hashes = _external_families(extra_condition_files)
    overlap = set(families).intersection(external_families)
    if overlap:
        raise ValueError(f"external family names clash with pinned data: {sorted(overlap)}")
    families.update(external_families)
    pending = {
        "icarm_262_and_other_subcutoff_rank20_fibres": (
            "pending: this checkout has no pinned Mestre/root data and independently "
            "certified exceptional direction list from which z^2=g(T) can be formed"
        ),
        "eventual_rootless_K3_bisections": (
            "pending: lattice bisection exclusion data is present, but no explicit "
            "bisection double-cover equations z^2=g(u) are recorded yet"
        ),
    }
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.exceptional-branch-code.v1",
        "status": (
            "complete exact computation for five pinned codes plus two exact, "
            "bounded-source supplementary preimage lists; extensible inputs pending"
        ),
        "claim_level": "exact divisor computation; no rational-point search used",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "families": families,
        "supplemental_bounded_input_families": supplemental_bounded_families,
        "pending_inputs": pending,
        "external_condition_inputs": external_hashes,
        "outcome": {
            name: {
                "code_dimension": data["code_dimension"],
                "minimum_quadratic_quotient_genus": data["minimum_quadratic_quotient"]["genus"],
                "low_genus_quotient_count": len(data["low_genus_quadratic_quotients"]),
                "complete_cancellation_count": len(data["complete_branch_cancellations"]),
                "shared_branch_pair_count": len(data["shared_branch_pairs"]),
            }
            for name, data in families.items()
        },
        "supplemental_bounded_input_outcome": {
            name: {
                "code_dimension": data["code_dimension"],
                "minimum_quadratic_quotient_genus": data["minimum_quadratic_quotient"]["genus"],
                "low_genus_quotient_count": len(data["low_genus_quadratic_quotients"]),
                "complete_cancellation_count": len(data["complete_branch_cancellations"]),
                "shared_branch_pair_count": len(data["shared_branch_pairs"]),
            }
            for name, data in supplemental_bounded_families.items()
        },
        "sources": {
            "branch_code_engine_sha256": sha256_file(
                root / "elliptic-curves/cas/branch_divisor_code.py"
            ),
            "fermigier_transport_script_sha256": sha256_file(
                root / "elliptic-curves/cas/analyze_fermigier_exceptional_transport.py"
            ),
            "icarm275_rank20_verifier_sha256": sha256_file(
                root / "elliptic-curves/cas/verify_icarm_curve275_mestre_rank20.py"
            ),
            "icarm245_275_cross_shape_source_sha256": sha256_file(
                root / "archive/elliptic-curves/cas/analyze_icarm245_275_cross_shape_transport.sage"
            ),
            "icarm243_accidental_slice_source_sha256": sha256_file(
                root / "archive/elliptic-curves/cas/search_icarm_curve243_accidental_slices.py"
            ),
            "icarm226_accidental_slice_source_sha256": sha256_file(
                root / "archive/elliptic-curves/cas/search_icarm_curve226_accidental_slices.py"
            ),
            "script_sha256": sha256_file(Path(__file__)),
        },
        "software": {"python": platform.python_version(), "sympy": sp.__version__},
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/analyze_exceptional_branch_code.py"
        ),
    }
    artifact["result_sha256"] = result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/generated-results/elliptic-curves/elliptic_exceptional_branch_code.json",
    )
    parser.add_argument(
        "--conditions-file",
        action="append",
        type=Path,
        default=[],
        help=(
            "exact JSON input with families:[{name,parameter,conditions:[{label,"
            "numerator_coefficients_ascending,denominator_coefficients_ascending}]}]; "
            "may be supplied more than once"
        ),
    )
    args = parser.parse_args()
    artifact = run(root, extra_condition_files=args.conditions_file)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))


if __name__ == "__main__":
    main()

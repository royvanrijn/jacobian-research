#!/usr/bin/env python3
"""Test higher corrected moments on the exact cubic Q-residual survivor.

The general t0-open moment export becomes very large at order nine.  On the
slice s1=ell=Q=0 we instead specialize inside Singular *before* applying the
three linear chart pivots.  The resulting polynomial is reduced in the
checkpointed length-four Kummer quotient, then tested on the cubic factor of
the mu6 norm by exact characteristic-zero arithmetic.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import sympy as sp

from explore_two_pair_sic_bidegree33_full_anchor import PARAMETERS
from research_two_pair_sic_bidegree33_t0_Q_residual import (
    ExtensionArithmetic,
    determinant_over_extension,
    evaluate_fraction_field_element_at_algebraic_root,
    normal_form_bipolynomial,
    parse_extension_bipolynomial,
    specialize_extension_element_to_number_field,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    exact_moment_terms,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    ROOT,
    fixed_chart_expression_exact,
    localize_inverse_pair,
)


DEFAULT_CHECKPOINT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "slice_s1_0_ell_0_matrices678_exact.json"
    )
)
NORM_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "slice_s1_0_ell_0_exact.json"
    )
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--order",
        type=int,
        choices=(9, 10, 11, 12, 14),
        default=9,
    )
    parser.add_argument("--factor-degree", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--normal-form-only", action="store_true")
    parser.add_argument("--resume-normal-form", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def specialized_moment_export(
    singular: str,
    order: int,
    timeout: int,
) -> str:
    """Export one moment after early specialization on the Kummer slice."""

    expression = fixed_chart_expression_exact(
        exact_moment_terms(order),
        PARAMETERS.index("t0"),
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring source=0,(
  s6,s5,t4,t3,s4,s0,s1,s2,s3,t1,t2,uinv,a
),dp;
poly value={expression};
value=subst(value,s1,0);
value=subst(value,t1,0);
value=subst(value,s2,(-13/3)*uinv);
value=subst(value,t2,a*uinv^2);

poly t3Value=-(2*s0*s3-3*s0)*uinv^2;
poly s4Value=(
  -169*s0*uinv^2+(182/3)*uinv-2*a*s0*uinv^2
)*uinv^2*(1/3);
poly t4Value=(
  3*s0*s6+45*((-13/3)*uinv)*s4-30*s3^2
  -42*(a*uinv^2)^2-70
)*(1/14);

value=subst(value,t4,t4Value);
value=subst(value,t3,t3Value);
value=subst(value,s4,s4Value);
ideal inverseBasis=std(s0*uinv-1);
value=reduce(value,inverseBasis);
print("VALUE "+string(value));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        "\n".join(
            line
            for line in completed.stdout.splitlines()
            if "?" in line or "error" in line.lower()
        )[-4000:],
        completed.stderr[-2000:],
    )
    marker = re.search(r"(?m)^VALUE (.*)$", completed.stdout)
    assert marker is not None
    return localize_inverse_pair(marker.group(1))


def parse_matrix(
    serialized: list[list[str]],
    arithmetic: ExtensionArithmetic,
    locals_map: dict[str, sp.Symbol],
) -> list[list[sp.Poly]]:
    return [
        [
            arithmetic.make(sp.sympify(entry, locals=locals_map))
            for entry in row
        ]
        for row in serialized
    ]


def specialize_matrix(
    matrix: list[list[sp.Poly]],
    source: ExtensionArithmetic,
    target: ExtensionArithmetic,
) -> list[list[sp.Poly]]:
    return [
        [
            specialize_extension_element_to_number_field(
                entry,
                source,
                target,
            )
            for entry in row
        ]
        for row in matrix
    ]


def full_rank_minor_with_new_column(
    arithmetic: ExtensionArithmetic,
    matrices: list[list[list[sp.Poly]]],
    column: list[sp.Poly],
) -> list[int] | None:
    size = len(column)
    old_columns = [
        [matrix[row][matrix_column] for row in range(size)]
        for matrix in matrices
        for matrix_column in range(size)
    ]
    for selected in combinations(range(len(old_columns)), size - 1):
        minor_columns = [old_columns[index] for index in selected] + [column]
        minor = [
            [minor_columns[column_index][row] for column_index in range(size)]
            for row in range(size)
        ]
        if not determinant_over_extension(arithmetic, minor).is_zero:
            return [*selected, len(old_columns)]
    return None


def main() -> None:
    arguments = parse_arguments()
    assert 1 <= arguments.timeout <= 60
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")

    checkpoint_path = arguments.checkpoint
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    checkpoint = payload["matrix_checkpoint"]
    if "monic_fiber_basis" not in checkpoint or "pivot_s3" not in checkpoint:
        raise ValueError(
            "checkpoint predates the resumable fibre basis; regenerate it"
        )

    parameter = sp.symbols(checkpoint["base_parameter"])
    root = sp.symbols(checkpoint["extension_symbol"])
    base = sp.QQ.frac_field(parameter)
    locals_map = {str(parameter): parameter, str(root): root}
    arithmetic = ExtensionArithmetic(
        base,
        root,
        sp.sympify(checkpoint["modulus"], locals=locals_map),
    )
    pivot = arithmetic.make(
        sp.sympify(checkpoint["pivot_s3"], locals=locals_map)
    )
    basis = [
        {
            tuple(map(int, monomial.split(","))): arithmetic.make(
                sp.sympify(coefficient, locals=locals_map)
            )
            for monomial, coefficient in polynomial.items()
        }
        for polynomial in checkpoint["monic_fiber_basis"]
    ]
    monomials = [(0, 0), (0, 1), (1, 0), (0, 2)]

    if arguments.resume_normal_form is None:
        exported = specialized_moment_export(
            singular,
            arguments.order,
            arguments.timeout,
        )
        projected = re.sub(r"\ba\b", "T", exported)
        parsed = parse_extension_bipolynomial(
            projected,
            arithmetic,
            [parameter],
            pivot,
            root,
        )
        normal = normal_form_bipolynomial(arithmetic, parsed, basis)
        assert set(normal) <= set(monomials)
        vector = [
            normal.get(monomial, arithmetic.zero) for monomial in monomials
        ]
        export_length = len(exported)
        export_terms = len(re.findall(r"[+-]?[^+-]+", exported))
        normal_support = len(normal)
    else:
        normal_form_path = arguments.resume_normal_form
        if not normal_form_path.is_absolute():
            normal_form_path = ROOT / normal_form_path
        normal_payload = json.loads(
            normal_form_path.read_text(encoding="utf-8")
        )
        assert normal_payload["order"] == arguments.order
        vector = [
            arithmetic.make(sp.sympify(entry, locals=locals_map))
            for entry in normal_payload["normal_form_vector"]
        ]
        export_length = normal_payload["specialized_export_length"]
        export_terms = normal_payload["specialized_export_term_count"]
        normal_support = normal_payload["normal_form_support"]

    if arguments.normal_form_only:
        result = {
            "format": (
                "two-pair-sic-bidegree33-t0-Q-residual-"
                "higher-normal-form-v1"
            ),
            "status": (
                "exact higher-moment normal form in the generic "
                "length-four Kummer quotient; factor-field test not run"
            ),
            "slice": {"s1": 0, "ell": 0, "Q": 0},
            "order": arguments.order,
            "specialized_export_length": export_length,
            "specialized_export_term_count": export_terms,
            "normal_form_support": normal_support,
            "normal_form_vector": [
                sp.sstr(entry.as_expr()) for entry in vector
            ],
            "checkpoint": str(checkpoint_path),
            "reproduction_command": " ".join(sys.argv),
        }
        if arguments.output is not None:
            output = arguments.output
            if not output.is_absolute():
                output = ROOT / output
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(result, indent=2, sort_keys=True))
        return

    matrices = [
        parse_matrix(checkpoint[name], arithmetic, locals_map)
        for name in ("M6", "M7", "M8")
    ]
    norm_payload = json.loads(NORM_ARTIFACT.read_text(encoding="utf-8"))
    factors = []
    for record in norm_payload["mu6_norm"]["numerator_factors"]:
        factor = sp.Poly(
            sp.sympify(record["factor"], locals=locals_map),
            parameter,
            domain=sp.QQ,
        )
        if factor.degree() == arguments.factor_degree:
            factors.append((record, factor))
    assert len(factors) == 1
    factor_record, factor = factors[0]

    number_field = sp.QQ.alg_field_from_poly(
        factor,
        alias=f"theta_{factor.degree()}",
    )
    specialized_modulus = sp.Poly.from_dict(
        {
            (power,): evaluate_fraction_field_element_at_algebraic_root(
                arithmetic.base.convert(arithmetic.modulus.nth(power)),
                arithmetic.base,
                number_field,
            )
            for power in range(arithmetic.degree + 1)
            if arithmetic.modulus.nth(power) != 0
        },
        root,
        domain=number_field,
    )
    specialized_arithmetic = ExtensionArithmetic(
        number_field,
        root,
        specialized_modulus.as_expr(),
    )
    specialized_matrices = [
        specialize_matrix(matrix, arithmetic, specialized_arithmetic)
        for matrix in matrices
    ]
    specialized_vector = [
        specialize_extension_element_to_number_field(
            entry,
            arithmetic,
            specialized_arithmetic,
        )
        for entry in vector
    ]
    minor = full_rank_minor_with_new_column(
        specialized_arithmetic,
        specialized_matrices,
        specialized_vector,
    )

    result = {
        "format": (
            "two-pair-sic-bidegree33-t0-Q-residual-higher-moment-v1"
        ),
        "status": (
            "exact characteristic-zero higher-moment test on the cubic "
            "mu6-norm component of the one-parameter Q-residual slice"
        ),
        "slice": {"s1": 0, "ell": 0, "Q": 0},
        "factor": str(factor_record["factor"]),
        "factor_degree": factor.degree(),
        "order": arguments.order,
        "specialized_export_length": export_length,
        "specialized_export_term_count": export_terms,
        "normal_form_support": normal_support,
        "normal_form_vector": [
            sp.sstr(entry.as_expr()) for entry in vector
        ],
        "full_rank_after_adjoining_moment": minor is not None,
        "nonzero_minor_columns_zero_based": minor,
        "component_excluded_through_order": minor is not None,
        "checkpoint": str(checkpoint_path),
        "reproduction_command": " ".join(sys.argv),
    }
    if arguments.output is not None:
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

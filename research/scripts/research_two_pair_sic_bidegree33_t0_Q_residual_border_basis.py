#!/usr/bin/env python3
"""Specialize the exact Q-stratum fibre basis to the residual R20 branch.

The generic Q-stratum quotient by (mu4,mu5) has a three-element lifted
Groebner basis with polynomial transformation matrix.  Its only
nonconstant leading border Delta satisfies

    Res_s3(mu3,Delta) = c*u^20*J_Q^4*R20^2.

On the irreducible R20 component, the first subresultant gives the dense
linear pivot s3=-B/A.  This driver first inverts A modulo R20 over the
function field GF(p)(s1,ell,u), then specializes the already-computed
three-element basis.  The resulting standard-basis calculation has only
the ordinary variables (T,s6,s5).

This is a modular research calculation.  A successful generic basis or
unit-ideal test is not a characteristic-zero certificate, and every
cleared coefficient denominator must be treated separately.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import sympy as sp

from research_two_pair_sic_bidegree33_t0_Q_residual import (
    LEADING_ARTIFACT,
    ROOT,
    adapted_polynomial,
    residual_data,
    specialized_polynomial,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)
from research_two_pair_sic_bidegree33_t0_fitting_interpolation import (
    polynomial_evaluate,
    rational_reconstruct,
)


BASIS_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_basis_exact.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_residual_border_basis_mod1000003.json"
)
REDUCED_GENERATORS_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_residual_generators_mod1000003.json"
)
MOMENTS_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1000003)
    parser.add_argument(
        "--stage",
        choices=("specialize", "basis", "moments", "unit", "pencil"),
        default="basis",
    )
    parser.add_argument(
        "--pivot-mode",
        choices=("equation", "substitution", "algebraic", "quotient"),
        default="equation",
        help=(
            "keep A*s3+B as a linear base equation, or first invert A "
            "modulo R20 and substitute s3; algebraic first loads the "
            "small equation-mode generators and then treats T as a "
            "degree-five algebraic coefficient; quotient uses the "
            "same small generators in a block-ordered quotient ring"
        ),
    )
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--specialize",
        action="append",
        default=[],
        metavar="VARIABLE=VALUE",
        help="fix any of s1,ell,u modulo the selected prime",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--singular-output",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--include-specialized",
        action="store_true",
        help="print the three specialized generic-basis generators",
    )
    parser.add_argument(
        "--precondition",
        action="store_true",
        help=(
            "in quotient mode, adjoin the explicit first quadratic "
            "S-pair before the standard-basis calculation"
        ),
    )
    parser.add_argument(
        "--original-only",
        action="store_true",
        help="build the residual basis from the original mu4,mu5 only",
    )
    parser.add_argument(
        "--moments-artifact",
        type=Path,
        default=None,
        help="reuse an explicitly exported exact corrected-moment input",
    )
    parser.add_argument(
        "--write-moments-artifact",
        type=Path,
        default=None,
        help="write the exact corrected-moment input used by this run",
    )
    parser.add_argument(
        "--scan-variable",
        choices=("s1", "ell", "u"),
        default=None,
        help=(
            "batch the closed-fibre pencil over the listed values of "
            "one base variable; specialize the other two"
        ),
    )
    parser.add_argument(
        "--scan-values",
        default=None,
        help="comma-separated integers or START:STOP for --scan-variable",
    )
    parser.add_argument(
        "--reconstruct-training-count",
        type=int,
        default=None,
        help=(
            "rationally reconstruct every pencil coefficient along a "
            "scan line and validate on all remaining good points"
        ),
    )
    parser.add_argument(
        "--reduced-generators-artifact",
        type=Path,
        default=REDUCED_GENERATORS_ARTIFACT,
        help="equation-mode generator export used by algebraic pivot mode",
    )
    return parser.parse_args()


def invert_modulus(
    singular: str,
    prime: int,
    coefficient_parameters: list[str],
    modulus: str,
    value: str,
    timeout: int,
) -> tuple[str, float]:
    started = time.monotonic()
    coefficient_declaration = ",".join(
        [str(prime), *coefficient_parameters]
    )
    if coefficient_parameters:
        coefficient_declaration = f"({coefficient_declaration})"
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring inverseRing={coefficient_declaration},(T),dp;
option(redSB);
poly modulus={modulus};
poly value={value};
ideal inputIdeal=modulus,value;
matrix transformation;
ideal G=liftstd(inputIdeal,transformation);
poly inverse=transformation[2,1]/G[1];
ideal modulusBasis=std(ideal(modulus));
inverse=reduce(inverse,modulusBasis);
print(
  "INVERSE_META "+string(size(G))+" "
  +string(reduce(value*inverse-1,modulusBasis)==0)
);
print("INVERSE "+string(inverse));
""",
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )
    if completed.returncode != 0 or "\n   ? " in completed.stdout:
        raise RuntimeError(
            completed.stdout[-8000:] + completed.stderr[-4000:]
        )
    meta = re.search(
        r"(?m)^INVERSE_META (\d+) ([01])$",
        completed.stdout,
    )
    inverse = re.search(r"(?m)^INVERSE (.*)$", completed.stdout)
    if (
        meta is None
        or meta.groups() != ("1", "1")
        or inverse is None
    ):
        raise RuntimeError(completed.stdout[-8000:])
    return inverse.group(1), time.monotonic() - started


def adapted_expression(expression: str) -> str:
    return substitute(
        expression,
        (
            ("t1", "(s1*u-ell)"),
            ("t2", "(T*u^2)"),
        ),
    )


def exact_expand_expression(expression: str) -> str:
    """Expand a closed-fibre expression over QQ before modular parsing."""

    s6, s5, s3, T = sp.symbols("s6 s5 s3 T")
    parsed = sp.sympify(
        expression.replace("^", "**"),
        locals={"s6": s6, "s5": s5, "s3": s3, "T": T},
    )
    polynomial = sp.Poly(
        sp.expand(parsed),
        s6,
        s5,
        s3,
        T,
        domain=sp.QQ,
    )
    _, primitive = polynomial.primitive()
    return sp.sstr(primitive.as_expr()).replace("**", "^")


def interpolate_mod_prime(
    values: list[tuple[int, int]],
    prime: int,
) -> list[int]:
    """Return ascending coefficients of the polynomial through ``values``."""

    answer = [0] * len(values)
    for index, (node, value) in enumerate(values):
        basis = [1]
        denominator = 1
        for other_index, (other_node, _) in enumerate(values):
            if other_index == index:
                continue
            multiplied = [0] * (len(basis) + 1)
            for degree, coefficient in enumerate(basis):
                multiplied[degree] -= other_node * coefficient
                multiplied[degree + 1] += coefficient
            basis = [coefficient % prime for coefficient in multiplied]
            denominator = denominator * (node - other_node) % prime
        scale = value * pow(denominator, -1, prime) % prime
        for degree, coefficient in enumerate(basis):
            answer[degree] = (
                answer[degree] + scale * coefficient
            ) % prime
    return answer


def parse_scan_values(specification: str, prime: int) -> list[int]:
    """Parse a deterministic list of distinct finite-field values."""

    range_match = re.fullmatch(r"(-?\d+):(-?\d+)", specification)
    if range_match is not None:
        start, stop = map(int, range_match.groups())
        raw_values = list(range(start, stop))
    else:
        raw_values = [
            int(value)
            for value in specification.split(",")
            if value.strip()
        ]
    values = [value % prime for value in raw_values]
    if not values or len(set(values)) != len(values):
        raise ValueError("scan values must be nonempty and distinct modulo p")
    return values


def run_pencil_scan(
    arguments: argparse.Namespace,
    singular: str,
    assignments: dict[str, int],
    residual: str,
    pivot_a: str,
    pivot_b: str,
    pivot_scale: sp.Rational,
) -> None:
    """Evaluate the closed-fibre determinant pencil at a batch of points."""

    if (
        arguments.stage != "pencil"
        or arguments.pivot_mode != "equation"
        or not arguments.original_only
        or arguments.scan_values is None
    ):
        raise ValueError(
            "a pencil scan requires --stage pencil --pivot-mode equation "
            "--original-only and --scan-values"
        )
    scan_variable = arguments.scan_variable
    assert scan_variable is not None
    fixed_variables = {
        "s1",
        "ell",
        "u",
    } - {scan_variable}
    if set(assignments) != fixed_variables:
        raise ValueError(
            "specialize exactly the two base variables not being scanned"
        )
    values = parse_scan_values(
        arguments.scan_values,
        arguments.prime,
    )
    if scan_variable == "u" and 0 in values:
        raise ValueError("the residual chart requires u != 0")

    if arguments.moments_artifact is not None:
        moments_path = arguments.moments_artifact
        if not moments_path.is_absolute():
            moments_path = ROOT / moments_path
        moments_payload = json.loads(
            moments_path.read_text(encoding="utf-8")
        )
        if (
            moments_payload.get("format")
            != "two-pair-sic-bidegree33-t0-Q-corrected-moments-v1"
            or moments_payload.get("coordinate_ell") != "s1*u-t1"
            or int(moments_payload.get("through", 0)) < 7
        ):
            raise ValueError("invalid corrected-moment scan input")
        moments = {
            int(order): polynomial
            for order, polynomial in moments_payload[
                "moment_polynomials"
            ].items()
        }
        moments_source = str(moments_path)
    else:
        export = t0_open_localized_export(
            singular,
            tuple(range(2, 8)),
            0,
            min(arguments.timeout, 300),
        )
        moments = dict(
            zip(range(3, 8), export["polynomials"][:-1], strict=True)
        )
        moments_source = "fresh exact corrected-moment export"

    q_replacement = (("s2", "(s1^2*u-(13/3)*u)"),)
    source_moments = {
        order: adapted_expression(
            substitute(moments[order], q_replacement)
        )
        for order in range(3, 8)
    }
    leading_payload = json.loads(
        LEADING_ARTIFACT.read_text(encoding="utf-8")
    )
    border = adapted_expression(
        leading_payload["leading_coefficient_lcm"]
    )
    pivot_scale_text = str(pivot_scale).replace("**", "^")
    declarations = "\n".join(
        f"poly sourceP{order}={source_moments[order]};"
        for order in range(3, 8)
    )
    procedure = r"""
proc residualPencilSample(ideal baseData, ideal moments, int sampleId)
{
  option(redSB);
  poly modulus=baseData[1];
  poly pivotA=baseData[2];
  poly pivotB=baseData[3];
  poly border=baseData[4];
  poly p3=baseData[5];
  poly p4=moments[1];
  poly p5=moments[2];
  poly p6=moments[3];
  poly p7=moments[4];
  poly pivotRelation=pivotA*s3+pivotB;
  ideal Rbasis=std(ideal(modulus,pivotRelation));
  ideal pivotOpen=std(ideal(modulus,pivotA));
  border=reduce(border,Rbasis);
  p3=reduce(p3,Rbasis);
  int pivotOk=(reduce(1,pivotOpen)==0);
  int borderOk=(border==0);
  int baseOk=(p3==0);
  ideal G=std(ideal(modulus,pivotRelation,p4,p5));
  int quotientDimension=vdim(G);
  if(
    (pivotOk==0) || (borderOk==0) || (baseOk==0)
    || (dim(G)!=0) || (quotientDimension!=20)
  )
  {
    print(
      "SCAN_BAD "+string(sampleId)+" "+string(pivotOk)+" "
      +string(borderOk)+" "+string(baseOk)+" "
      +string(dim(G))+" "+string(quotientDimension)
    );
    return();
  }
  p6=reduce(p6,G);
  p7=reduce(p7,G);
  ideal standardMonomials=kbase(G);
  matrix M6[20][20];
  matrix M7[20][20];
  poly remainder;
  matrix coefficients;
  int columnIndex;
  int termIndex;
  int rowIndex;
  int found;
  for(columnIndex=1;columnIndex<=20;columnIndex++)
  {
    remainder=reduce(p6*standardMonomials[columnIndex],G);
    coefficients=coef(remainder,T*s3*s6*s5);
    for(termIndex=1;termIndex<=ncols(coefficients);termIndex++)
    {
      found=0;
      for(rowIndex=1;rowIndex<=20;rowIndex++)
      {
        if(coefficients[1,termIndex]==standardMonomials[rowIndex])
        {
          M6[rowIndex,columnIndex]=coefficients[2,termIndex];
          found=1;
        }
      }
      if(found==0) { print("SCAN_UNMATCHED_M6"); }
    }
    remainder=reduce(p7*standardMonomials[columnIndex],G);
    coefficients=coef(remainder,T*s3*s6*s5);
    for(termIndex=1;termIndex<=ncols(coefficients);termIndex++)
    {
      found=0;
      for(rowIndex=1;rowIndex<=20;rowIndex++)
      {
        if(coefficients[1,termIndex]==standardMonomials[rowIndex])
        {
          M7[rowIndex,columnIndex]=coefficients[2,termIndex];
          found=1;
        }
      }
      if(found==0) { print("SCAN_UNMATCHED_M7"); }
    }
  }
  matrix jointMatrix[20][40];
  for(columnIndex=1;columnIndex<=20;columnIndex++)
  {
    for(rowIndex=1;rowIndex<=20;rowIndex++)
    {
      jointMatrix[rowIndex,columnIndex]=M6[rowIndex,columnIndex];
      jointMatrix[rowIndex,columnIndex+20]=M7[rowIndex,columnIndex];
    }
  }
  print(
    "SCAN_META "+string(sampleId)+" "+string(size(G))+" "
    +string(rank(M6))+" "+string(rank(M7))+" "
    +string(rank(jointMatrix))
  );
  int zValue;
  matrix pencilMatrix[20][20];
  poly pencilValue;
  for(zValue=0;zValue<=20;zValue++)
  {
    pencilMatrix=M6+zValue*M7;
    pencilValue=det(pencilMatrix);
    print(
      "SCAN_VALUE "+string(sampleId)+" "+string(zValue)+" "
      +string(pencilValue)
    );
  }
}
"""
    source_program = f"""
ring source={arguments.prime},(
  T,s3,s6,s5,s1,ell,u
),dp;
poly sourceModulus={residual};
poly sourcePivotA={pivot_a};
poly sourcePivotB={pivot_b};
poly sourceBorder={border};
{declarations}
{procedure}
"""
    blocks = [source_program]
    points: list[dict[str, int]] = []
    for sample_id, value in enumerate(values):
        point = {**assignments, scan_variable: value}
        points.append(point)
        blocks.append(
            f"""
ring sample{sample_id}={arguments.prime},(
  T,s3,s6,s5
),dp;
map evaluate{sample_id}=source,T,s3,s6,s5,
  {point["s1"]},{point["ell"]},{point["u"]};
poly modulus=evaluate{sample_id}(sourceModulus);
poly pivotA=evaluate{sample_id}(sourcePivotA);
poly pivotB=({pivot_scale_text})*evaluate{sample_id}(sourcePivotB);
poly border=evaluate{sample_id}(sourceBorder);
poly p3=evaluate{sample_id}(sourceP3);
poly p4=evaluate{sample_id}(sourceP4);
poly p5=evaluate{sample_id}(sourceP5);
poly p6=evaluate{sample_id}(sourceP6);
poly p7=evaluate{sample_id}(sourceP7);
ideal baseData=modulus,pivotA,pivotB,border,p3;
ideal moments=p4,p5,p6,p7;
residualPencilSample(baseData,moments,{sample_id});
"""
        )
    program = "\n".join(blocks)
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
        timeout=arguments.timeout,
    )
    elapsed = time.monotonic() - started
    transcript = completed.stdout
    if completed.stderr:
        transcript += "\nSTDERR\n" + completed.stderr
    if arguments.singular_output is not None:
        transcript_path = arguments.singular_output
        if not transcript_path.is_absolute():
            transcript_path = ROOT / transcript_path
        transcript_path.parent.mkdir(parents=True, exist_ok=True)
        transcript_path.write_text(transcript, encoding="utf-8")
    else:
        transcript_path = None
    failure_markers = (
        "\n   ? ",
        "Flint exception",
        "Singular : signal",
        "Segment fault",
        "Bus error",
        "SCAN_UNMATCHED",
    )
    if completed.returncode != 0 or any(
        marker in transcript for marker in failure_markers
    ):
        raise RuntimeError(transcript[-12000:])

    bad = {
        int(sample): {
            "pivot_open": pivot == "1",
            "border_zero": border_zero == "1",
            "mu3_zero": base == "1",
            "dimension": int(dimension),
            "quotient_dimension": int(length),
        }
        for sample, pivot, border_zero, base, dimension, length in re.findall(
            r"(?m)^SCAN_BAD (\d+) ([01]) ([01]) ([01]) (-?\d+) (-?\d+)$",
            completed.stdout,
        )
    }
    meta = {
        int(sample): {
            "basis_size": int(size),
            "rank_mu6": int(rank6),
            "rank_mu7": int(rank7),
            "joint_block_rank": int(joint),
        }
        for sample, size, rank6, rank7, joint in re.findall(
            r"(?m)^SCAN_META (\d+) (\d+) (\d+) (\d+) (\d+)$",
            completed.stdout,
        )
    }
    raw_values: dict[int, list[tuple[int, int]]] = {}
    for sample, node, determinant in re.findall(
        r"(?m)^SCAN_VALUE (\d+) (\d+) (-?\d+)$",
        completed.stdout,
    ):
        raw_values.setdefault(int(sample), []).append(
            (int(node), int(determinant) % arguments.prime)
        )
    if set(meta) | set(bad) != set(range(len(points))):
        raise RuntimeError(transcript[-12000:])

    records = []
    for sample_id, point in enumerate(points):
        if sample_id in bad:
            records.append(
                {
                    "sample_id": sample_id,
                    "point": point,
                    "good": False,
                    "failure": bad[sample_id],
                }
            )
            continue
        evaluations = sorted(raw_values.get(sample_id, []))
        if len(evaluations) != 21:
            raise RuntimeError(transcript[-12000:])
        coefficients = interpolate_mod_prime(
            evaluations,
            arguments.prime,
        )
        for node, value in evaluations:
            replay = sum(
                coefficient * pow(node, degree, arguments.prime)
                for degree, coefficient in enumerate(coefficients)
            ) % arguments.prime
            if replay != value:
                raise AssertionError("scan pencil interpolation failed")
        records.append(
            {
                "sample_id": sample_id,
                "point": point,
                "good": True,
                **meta[sample_id],
                "coefficients_ascending": coefficients,
                "identically_zero": not any(coefficients),
            }
        )
    reconstruction: dict[str, object] | None = None
    if arguments.reconstruct_training_count is not None:
        good_records = [record for record in records if record["good"]]
        training_count = arguments.reconstruct_training_count
        if (
            training_count < 2
            or training_count >= len(good_records)
        ):
            raise ValueError(
                "the reconstruction training count must leave at least "
                "one good validation point"
            )
        line_arguments = [
            record["point"][scan_variable] for record in good_records
        ]
        reconstructed = []
        numerator_polynomials = []
        denominator_polynomials = []
        reconstruction_variable = sp.symbols("X")
        for coefficient_index in range(21):
            line_values = [
                record["coefficients_ascending"][coefficient_index]
                for record in good_records
            ]
            numerator, denominator = rational_reconstruct(
                line_arguments[:training_count],
                line_values[:training_count],
                arguments.prime,
            )
            validation = [
                (
                    polynomial_evaluate(
                        denominator,
                        argument,
                        arguments.prime,
                    )
                    != 0
                    and polynomial_evaluate(
                        numerator,
                        argument,
                        arguments.prime,
                    )
                    == value
                    * polynomial_evaluate(
                        denominator,
                        argument,
                        arguments.prime,
                    )
                    % arguments.prime
                )
                for argument, value in zip(
                    line_arguments[training_count:],
                    line_values[training_count:],
                    strict=True,
                )
            ]
            if not all(validation):
                raise AssertionError(
                    "line reconstruction failed held-out validation for "
                    f"pencil coefficient {coefficient_index}"
                )
            reconstructed.append(
                {
                    "coefficient_index": coefficient_index,
                    "numerator_degree": len(numerator) - 1,
                    "denominator_degree": len(denominator) - 1,
                    "numerator_coefficients_ascending": numerator,
                    "denominator_coefficients_ascending": denominator,
                }
            )
            numerator_polynomials.append(
                sp.Poly.from_list(
                    list(reversed(numerator)),
                    gens=reconstruction_variable,
                    modulus=arguments.prime,
                )
            )
            denominator_polynomials.append(
                sp.Poly.from_list(
                    list(reversed(denominator)),
                    gens=reconstruction_variable,
                    modulus=arguments.prime,
                )
            )
        numerator_gcd = numerator_polynomials[0]
        for polynomial in numerator_polynomials[1:]:
            numerator_gcd = sp.gcd(numerator_gcd, polynomial)
        denominator_gcd = denominator_polynomials[0]
        for polynomial in denominator_polynomials[1:]:
            denominator_gcd = sp.gcd(denominator_gcd, polynomial)
        reconstruction = {
            "training_points": training_count,
            "validation_points": len(good_records) - training_count,
            "coefficients": reconstructed,
            "all_numerator_gcd_degree": numerator_gcd.degree(),
            "c0_c1_numerator_gcd_degree": sp.gcd(
                numerator_polynomials[0],
                numerator_polynomials[1],
            ).degree(),
            "common_denominator_gcd_degree": denominator_gcd.degree(),
            "line_excluded_away_from_denominators": (
                numerator_gcd.degree() == 0
            ),
        }
    payload = {
        "format": (
            "two-pair-sic-bidegree33-t0-Q-residual-pencil-scan-v1"
        ),
        "status": (
            f"bounded exact finite-field pencil scan modulo "
            f"{arguments.prime}; not a characteristic-zero certificate"
        ),
        "prime": arguments.prime,
        "scan_variable": scan_variable,
        "fixed_specializations": assignments,
        "points": len(points),
        "good_points": sum(record["good"] for record in records),
        "bad_points": sum(not record["good"] for record in records),
        "common_pencil_zero_points": [
            record["point"]
            for record in records
            if record["good"] and record["identically_zero"]
        ],
        "records": records,
        "line_reconstruction": reconstruction,
        "moments_source": moments_source,
        "solver_seconds": round(elapsed, 6),
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        "singular_output": (
            None if transcript_path is None else str(transcript_path)
        ),
        "reproduction_command": " ".join(sys.argv),
    }
    output = arguments.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> None:
    arguments = parse_arguments()
    if arguments.prime in (0, 2, 3, 5, 7, 13):
        raise ValueError("choose a prime avoiding the displayed denominators")
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")

    assignments: dict[str, int] = {}
    for item in arguments.specialize:
        match = re.fullmatch(r"(s1|ell|u)=(-?\d+)", item)
        if match is None:
            raise ValueError(f"invalid finite-field specialization: {item}")
        name, value = match.groups()
        assignments[name] = int(value) % arguments.prime
    coefficient_parameter_names = [
        name for name in ("s1", "ell", "u") if name not in assignments
    ]
    parameter_substitutions = tuple(
        (name, str(value)) for name, value in assignments.items()
    )

    started = time.monotonic()
    residual_raw, pivot_a_raw, pivot_b_raw = residual_data()
    residual, residual_terms, _, _ = adapted_polynomial(residual_raw)
    pivot_a, pivot_a_terms, _, pivot_a_content = adapted_polynomial(
        pivot_a_raw
    )
    pivot_b, pivot_b_terms, _, pivot_b_content = adapted_polynomial(
        pivot_b_raw
    )
    pivot_scale = pivot_b_content / pivot_a_content
    if arguments.scan_variable is not None:
        run_pencil_scan(
            arguments,
            singular,
            assignments,
            residual,
            pivot_a,
            pivot_b,
            pivot_scale,
        )
        return
    if parameter_substitutions:
        residual, _, _ = specialized_polynomial(residual, assignments)
        pivot_a, _, pivot_a_special_content = specialized_polynomial(
            pivot_a,
            assignments,
        )
        pivot_b, _, pivot_b_special_content = specialized_polynomial(
            pivot_b,
            assignments,
        )
        pivot_scale *= (
            pivot_b_special_content / pivot_a_special_content
        )
    pivot_scale_text = str(pivot_scale).replace("**", "^")
    if arguments.pivot_mode in (
        "substitution",
        "algebraic",
        "quotient",
    ):
        pivot_inverse, inverse_seconds = invert_modulus(
            singular,
            arguments.prime,
            coefficient_parameter_names,
            residual,
            pivot_a,
            arguments.timeout,
        )
    else:
        pivot_inverse = ""
        inverse_seconds = 0.0

    basis_payload = json.loads(BASIS_ARTIFACT.read_text(encoding="utf-8"))
    if basis_payload.get("prime") != 0 or basis_payload.get("stratum") != "Q":
        raise ValueError("the exact Q basis artifact is required")
    if arguments.pivot_mode in ("algebraic", "quotient"):
        generator_path = arguments.reduced_generators_artifact
        if not generator_path.is_absolute():
            generator_path = ROOT / generator_path
        generator_payload = json.loads(
            generator_path.read_text(encoding="utf-8")
        )
        if (
            generator_payload.get("prime") != arguments.prime
            or generator_payload.get("pivot_mode") != "equation"
            or len(generator_payload.get("specialized_generators", []))
            != 3
        ):
            raise ValueError(
                "algebraic mode requires a matching three-generator "
                "equation-mode export"
            )
        generic_basis = [
            substitute(polynomial, (("s3", "(pivotS)"),))
            for polynomial in generator_payload["specialized_generators"]
        ]
    else:
        generic_basis = [
            (
                substitute(
                    adapted_expression(polynomial),
                    (("s3", "(pivotS)"),),
                )
                if arguments.pivot_mode == "substitution"
                else adapted_expression(polynomial)
            )
            for polynomial in basis_payload["basis_polynomials"]
        ]
    if parameter_substitutions:
        generic_basis = [
            substitute(polynomial, parameter_substitutions)
            for polynomial in generic_basis
        ]
    if arguments.original_only:
        generic_basis = ["0", "0", "0"]
    leading_payload = json.loads(
        LEADING_ARTIFACT.read_text(encoding="utf-8")
    )
    border = (
        substitute(
            adapted_expression(leading_payload["leading_coefficient_lcm"]),
            (("s3", "(pivotS)"),),
        )
        if arguments.pivot_mode == "substitution"
        else adapted_expression(leading_payload["leading_coefficient_lcm"])
    )
    if parameter_substitutions:
        border = substitute(border, parameter_substitutions)
        if not coefficient_parameter_names:
            border = exact_expand_expression(border)

    input_moment_declarations = ""
    moment_declarations = ""
    moment_program = ""
    if arguments.stage != "specialize":
        moment_orders = (
            tuple(range(4, 8))
            if arguments.stage in ("moments", "unit", "pencil")
            else (4, 5)
        )
        if arguments.moments_artifact is not None:
            moments_path = arguments.moments_artifact
            if not moments_path.is_absolute():
                moments_path = ROOT / moments_path
            moments_payload = json.loads(
                moments_path.read_text(encoding="utf-8")
            )
            if (
                moments_payload.get("format")
                != "two-pair-sic-bidegree33-t0-Q-corrected-moments-v1"
                or moments_payload.get("coordinate_ell") != "s1*u-t1"
                or int(moments_payload.get("through", 0))
                < max(moment_orders)
            ):
                raise ValueError(
                    "the moments artifact does not contain the required "
                    "corrected exact moments"
                )
            moments = {
                int(order): polynomial
                for order, polynomial in moments_payload[
                    "moment_polynomials"
                ].items()
            }
        else:
            export = t0_open_localized_export(
                singular,
                tuple(range(2, max(moment_orders) + 1)),
                0,
                min(arguments.timeout, 300),
            )
            moments = dict(
                zip(
                    range(3, max(moment_orders) + 1),
                    export["polynomials"][:-1],
                    strict=True,
                )
            )
        if arguments.write_moments_artifact is not None:
            moments_path = arguments.write_moments_artifact
            if not moments_path.is_absolute():
                moments_path = ROOT / moments_path
            moments_path.parent.mkdir(parents=True, exist_ok=True)
            cached_moments = {
                order: moments[order]
                for order in range(3, max(moment_orders) + 1)
            }
            moments_path.write_text(
                json.dumps(
                    {
                        "format": (
                            "two-pair-sic-bidegree33-t0-Q-"
                            "corrected-moments-v1"
                        ),
                        "status": (
                            "exact characteristic-zero corrected moment "
                            "input; not a component certificate"
                        ),
                        "coordinate_ell": "s1*u-t1",
                        "through": max(moment_orders),
                        "moment_polynomials": cached_moments,
                        "moment_sha256": {
                            order: hashlib.sha256(
                                polynomial.encode()
                            ).hexdigest()
                            for order, polynomial in cached_moments.items()
                        },
                        "reproduction_command": " ".join(sys.argv),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
        q_replacement = (("s2", "(s1^2*u-(13/3)*u)"),)
        specialized_moments = {}
        for order in (3, *moment_orders):
            value = adapted_expression(
                substitute(moments[order], q_replacement)
            )
            if arguments.pivot_mode in (
                "substitution",
                "algebraic",
                "quotient",
            ):
                value = substitute(value, (("s3", "(pivotS)"),))
            if parameter_substitutions:
                value = substitute(value, parameter_substitutions)
                if not coefficient_parameter_names:
                    value = exact_expand_expression(value)
            specialized_moments[order] = value
        input_moment_declarations = "\n".join(
            f"poly p{order}={specialized_moments[order]};"
            for order in (3, 4, 5)
        )
        if arguments.stage in ("moments", "unit", "pencil"):
            moment_declarations = "\n".join(
                f"poly p{order}={specialized_moments[order]};"
                for order in (6, 7)
            )
            moment_program = """
p6=reduce(p6,G);
p7=reduce(p7,G);
print(
  "MOMENT_META "+string(size(p6))+" "+string(size(p7))+" "
  +string(deg(p6))+" "+string(deg(p7))
);
"""
            if arguments.stage == "unit":
                moment_program += """
ideal U=std(G+ideal(p6,p7));
print(
  "UNIT_META "+string(size(U))+" "
  +string(reduce(1,U)==0)
);
"""
            elif arguments.stage == "pencil":
                if (
                    arguments.pivot_mode != "equation"
                    or set(assignments) != {"s1", "ell", "u"}
                    or not arguments.original_only
                ):
                    raise ValueError(
                        "the pencil stage requires equation mode, all "
                        "three base specializations, and --original-only"
                    )
                moment_program += """
ideal standardMonomials=kbase(G);
int quotientDimension=size(standardMonomials);
matrix M6[quotientDimension][quotientDimension];
matrix M7[quotientDimension][quotientDimension];
poly remainder;
matrix coefficients;
int columnIndex;
int termIndex;
int rowIndex;
int found;
for(columnIndex=1;columnIndex<=quotientDimension;columnIndex++)
{
  remainder=reduce(p6*standardMonomials[columnIndex],G);
  coefficients=coef(remainder,T*s3*s6*s5);
  for(termIndex=1;termIndex<=ncols(coefficients);termIndex++)
  {
    found=0;
    for(rowIndex=1;rowIndex<=quotientDimension;rowIndex++)
    {
      if(coefficients[1,termIndex]==standardMonomials[rowIndex])
      {
        M6[rowIndex,columnIndex]=coefficients[2,termIndex];
        found=1;
      }
    }
    if(found==0) { print("UNMATCHED_M6"); }
  }
  remainder=reduce(p7*standardMonomials[columnIndex],G);
  coefficients=coef(remainder,T*s3*s6*s5);
  for(termIndex=1;termIndex<=ncols(coefficients);termIndex++)
  {
    found=0;
    for(rowIndex=1;rowIndex<=quotientDimension;rowIndex++)
    {
      if(coefficients[1,termIndex]==standardMonomials[rowIndex])
      {
        M7[rowIndex,columnIndex]=coefficients[2,termIndex];
        found=1;
      }
    }
    if(found==0) { print("UNMATCHED_M7"); }
  }
}
matrix jointMatrix[quotientDimension][2*quotientDimension];
for(columnIndex=1;columnIndex<=quotientDimension;columnIndex++)
{
  for(rowIndex=1;rowIndex<=quotientDimension;rowIndex++)
  {
    jointMatrix[rowIndex,columnIndex]=M6[rowIndex,columnIndex];
    jointMatrix[rowIndex,columnIndex+quotientDimension]=
      M7[rowIndex,columnIndex];
  }
}
print(
  "PENCIL_META "+string(quotientDimension)+" "
  +string(rank(M6))+" "+string(rank(M7))+" "
  +string(rank(jointMatrix))
);
int zValue;
matrix pencilMatrix[quotientDimension][quotientDimension];
poly pencilValue;
for(zValue=0;zValue<=quotientDimension;zValue++)
{
  pencilMatrix=M6+zValue*M7;
  pencilValue=det(pencilMatrix);
  print("PENCIL_VALUE "+string(zValue)+" "+string(pencilValue));
}
"""

    basis_declarations = "\n".join(
        f"poly g{index}={polynomial};"
        for index, polynomial in enumerate(generic_basis, 1)
    )
    specialized_prints = (
        "\n".join(
            (
                f'print("SPECIALIZED_GENERATOR {index} "'
                f"+string(g{index}));"
            )
            for index in range(1, 4)
        )
        if arguments.include_specialized
        else ""
    )
    coefficient_setup = ""
    precondition_program = ""
    if arguments.pivot_mode == "substitution":
        ring_variables = "T,s6,s5"
        ring_order = "(dp(1),dp(2))"
        pivot_setup = f"""
ideal Rbasis=std(ideal(modulus));
poly pivotInverse={pivot_inverse};
poly pivotS=reduce(-({pivot_scale_text})*pivotB*pivotInverse,Rbasis);
print(
  "PIVOT_META "+string(reduce(pivotA*pivotInverse-1,Rbasis)==0)+" "
  +string(deg(pivotS,T))+" "+string(size(pivotS))
);
"""
        basis_generators = "modulus,g1,g2,g3"
    elif arguments.pivot_mode == "algebraic":
        if arguments.stage != "specialize":
            raise ValueError(
                "algebraic mode only exports post-pivot coefficient "
                "representatives; use quotient mode for basis stages"
            )
        ring_variables = "s6,s5"
        ring_order = "dp"
        coefficient_setup = ""
        pivot_setup = f"""
ideal Rbasis=0;
poly pivotS={pivot_inverse};
poly pivotA={pivot_a};
poly pivotB={pivot_b};
print(
  "PIVOT_META 1 4 "+string(size(pivotS))
);
"""
        basis_generators = "g1,g2,g3"
    elif arguments.pivot_mode == "quotient":
        ring_variables = "s6,s5,T"
        ring_order = "(dp(2),dp(1))"
        coefficient_setup = f"""
ideal modulusIdeal=std(ideal({residual}));
qring residualQuotient=modulusIdeal;
"""
        pivot_setup = f"""
ideal Rbasis=0;
poly pivotS={pivot_inverse};
poly pivotA={pivot_a};
poly pivotB={pivot_b};
print(
  "PIVOT_META "
  +string(pivotA*pivotS+({pivot_scale_text})*pivotB==0)+" "
  +"4 "+string(size(pivotS))
);
"""
        precondition_program = """
poly b2=subst(subst(diff(diff(g2,s6),s5),s6,0),s5,0);
poly b3=subst(subst(diff(diff(g3,s6),s5),s6,0),s5,0);
poly h23=b2*g3-b3*g2;
print(
  "PRECONDITION_META "+string(size(b2))+" "
  +string(size(b3))+" "+string(size(h23))
);
""" if arguments.precondition else ""
        basis_generators = (
            "h23,g1,g2,g3"
            if arguments.precondition
            else "g1,g2,g3"
        )
    else:
        ring_variables = "T,s3,s6,s5"
        ring_order = "(dp(2),dp(2))"
        pivot_setup = f"""
poly pivotRelation=pivotA*s3+({pivot_scale_text})*pivotB;
ideal Rbasis=std(ideal(modulus,pivotRelation));
ideal pivotOpen=std(ideal(modulus,pivotA));
print(
  "PIVOT_META "+string(reduce(pivotRelation,Rbasis)==0)+" "
  +string(deg(pivotRelation))+" "+string(size(pivotRelation))
);
print(
  "PIVOT_OPEN "+string(reduce(1,pivotOpen)==0)+" "
  +string(size(pivotOpen))
);
"""
        basis_generators = "modulus,pivotRelation,g1,g2,g3"
    if arguments.stage == "specialize":
        basis_computation = ""
    else:
        if arguments.original_only:
            if arguments.pivot_mode != "equation":
                raise ValueError(
                    "--original-only currently requires equation mode"
                )
            basis_generators = "modulus,pivotRelation"
        basis_generators = f"{basis_generators},p4,p5"
        basis_computation = f"""
print(
  "ORIGINAL_INPUT "+string(size(p4))+" "+string(size(p5))
);
p3=reduce(p3,Rbasis);
print("BASE_CHECK "+string(size(p3)));
ideal G=std(ideal({basis_generators}));
print(
  "BASIS_META "+string(size(G))+" "+string(dim(G))+" "
  +string(vdim(G))
);
int basisIndex;
for(basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print(
    "BASIS_LEAD "+string(basisIndex)+" "
    +string(leadexp(G[basisIndex]))+" "
    +string(size(G[basisIndex]))
  );
}}
{moment_declarations}
{moment_program}
"""
    coefficient_parameters = (
        ",".join(
            [
                str(arguments.prime),
                *coefficient_parameter_names,
                "T",
            ]
        )
        if arguments.pivot_mode == "algebraic"
        else ",".join(
            [str(arguments.prime), *coefficient_parameter_names]
        )
    )
    coefficient_declaration = coefficient_parameters
    if coefficient_parameter_names or arguments.pivot_mode == "algebraic":
        coefficient_declaration = f"({coefficient_parameters})"
    base_declarations = (
        ""
        if arguments.pivot_mode in ("algebraic", "quotient")
        else f"""
poly modulus={residual};
poly pivotA={pivot_a};
poly pivotB={pivot_b};
"""
    )
    program = f"""
ring residualRing={coefficient_declaration},(
  {ring_variables}
),{ring_order};
{coefficient_setup}
option(redSB);
{base_declarations}
{pivot_setup}
poly border={
    "0"
    if arguments.pivot_mode in ("algebraic", "quotient")
    else border
};
border=reduce(border,Rbasis);
print("BORDER_META "+string(border==0)+" "+string(size(border)));
{basis_declarations}
g1=reduce(g1,Rbasis);
g2=reduce(g2,Rbasis);
g3=reduce(g3,Rbasis);
print(
  "SPECIALIZED_INPUT "+string(size(g1))+" "
  +string(size(g2))+" "+string(size(g3))
);
{precondition_program}
{specialized_prints}
{input_moment_declarations}
{basis_computation}
"""
    print(
        "RESIDUAL_BORDER_START "
        f"prime={arguments.prime} stage={arguments.stage} "
        f"pivot_mode={arguments.pivot_mode} "
        f"program_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )
    solver_started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
        timeout=arguments.timeout,
    )
    solver_seconds = time.monotonic() - solver_started
    transcript = completed.stdout
    if completed.stderr:
        transcript += "\nSTDERR\n" + completed.stderr
    if arguments.singular_output is not None:
        transcript_path = arguments.singular_output
        if not transcript_path.is_absolute():
            transcript_path = ROOT / transcript_path
        transcript_path.parent.mkdir(parents=True, exist_ok=True)
        transcript_path.write_text(transcript, encoding="utf-8")
    else:
        transcript_path = None
    singular_failure_markers = (
        "\n   ? ",
        "Flint exception",
        "Singular : signal",
        "Segment fault",
        "Bus error",
        "int overflow",
        "UNMATCHED_M6",
        "UNMATCHED_M7",
    )
    if completed.returncode != 0 or any(
        marker in completed.stdout
        for marker in singular_failure_markers
    ):
        raise RuntimeError(transcript[-12000:])

    pivot_meta = re.search(
        r"(?m)^PIVOT_META ([01]) (-?\d+) (\d+)$",
        completed.stdout,
    )
    pivot_open = re.search(
        r"(?m)^PIVOT_OPEN ([01]) (\d+)$",
        completed.stdout,
    )
    border_meta = re.search(
        r"(?m)^BORDER_META ([01]) (\d+)$",
        completed.stdout,
    )
    basis_meta = re.search(
        r"(?m)^BASIS_META (\d+) (-?\d+) (-?\d+)$",
        completed.stdout,
    )
    specialized_input = re.search(
        r"(?m)^SPECIALIZED_INPUT (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    original_input = re.search(
        r"(?m)^ORIGINAL_INPUT (\d+) (\d+)$",
        completed.stdout,
    )
    base_check = re.search(
        r"(?m)^BASE_CHECK (\d+)$",
        completed.stdout,
    )
    if (
        pivot_meta is None
        or border_meta is None
        or specialized_input is None
        or (
            arguments.pivot_mode == "equation"
            and pivot_open is None
        )
        or (
            arguments.stage != "specialize"
            and (
                basis_meta is None
                or original_input is None
                or base_check is None
            )
        )
    ):
        raise RuntimeError(transcript[-12000:])
    if border_meta.group(1) != "1":
        raise RuntimeError(
            "the specialized leading border did not vanish on the "
            f"residual pivot:\n{transcript[-12000:]}"
        )
    if (
        arguments.pivot_mode == "equation"
        and pivot_open is not None
        and pivot_open.group(1) != "1"
    ):
        raise RuntimeError(
            "the dense pivot coefficient is not invertible modulo "
            f"the specialized R20:\n{transcript[-12000:]}"
        )
    if base_check is not None and base_check.group(1) != "0":
        raise RuntimeError(
            "mu3 did not vanish modulo the residual pivot:\n"
            f"{transcript[-12000:]}"
        )
    if (
        basis_meta is not None
        and (
            int(basis_meta.group(1)) == 0
            or int(basis_meta.group(2)) != 0
            or int(basis_meta.group(3)) <= 0
        )
    ):
        raise RuntimeError(
            "the computed basis is not a nonzero zero-dimensional "
            f"quotient:\n{transcript[-12000:]}"
        )
    leads = [
        {
            "index": int(index),
            "exponents": [
                int(value) for value in exponent.split(",")
            ],
            "terms": int(terms),
        }
        for index, exponent, terms in re.findall(
            r"(?m)^BASIS_LEAD (\d+) ([0-9,]+) (\d+)$",
            completed.stdout,
        )
    ]
    payload: dict[str, object] = {
        "format": "two-pair-sic-bidegree33-t0-Q-residual-border-basis-v1",
        "status": (
            (
                f"exact finite-field closed-fibre calculation modulo "
                f"{arguments.prime}; not a characteristic-zero "
                "certificate"
            )
            if set(assignments) == {"s1", "ell", "u"}
            else (
                f"exact finite-field function-field calculation modulo "
                f"{arguments.prime}; not a characteristic-zero "
                "certificate"
            )
        ),
        "prime": arguments.prime,
        "stage": arguments.stage,
        "pivot_mode": arguments.pivot_mode,
        "coefficient_field": (
            f"GF({arguments.prime})"
            + (
                f"({','.join(coefficient_parameter_names)})"
                if coefficient_parameter_names
                else ""
            )
        ),
        "specializations": assignments,
        "ordinary_variables": ring_variables.split(","),
        "input_profiles": {
            "R20_terms": residual_terms,
            "A_terms": pivot_a_terms,
            "B_terms": pivot_b_terms,
            "generic_basis_lengths": [
                len(polynomial)
                for polynomial in basis_payload["basis_polynomials"]
            ],
        },
        "pivot": {
            "inverse_verified": (
                pivot_meta.group(1) == "1"
                if arguments.pivot_mode != "equation"
                else None
            ),
            "relation_verified": (
                pivot_meta.group(1) == "1"
                if arguments.pivot_mode == "equation"
                else None
            ),
            "A_invertible_mod_R20": (
                pivot_open.group(1) == "1"
                if pivot_open is not None
                else None
            ),
            "degree_T": int(pivot_meta.group(2)),
            "terms": int(pivot_meta.group(3)),
            "inverse_seconds": round(inverse_seconds, 6),
            "content_ratio_B_over_A": str(pivot_scale),
        },
        "border_zero_mod_R20": border_meta.group(1) == "1",
        "border_remainder_terms": int(border_meta.group(2)),
        "specialized_input_terms": [
            int(value) for value in specialized_input.groups()
        ],
        "original_moment_input_terms": (
            None
            if original_input is None
            else [int(value) for value in original_input.groups()]
        ),
        "mu3_remainder_terms": (
            None if base_check is None else int(base_check.group(1))
        ),
        "specialized_basis": (
            None
            if basis_meta is None
            else {
                "size": int(basis_meta.group(1)),
                "dimension": int(basis_meta.group(2)),
                "vector_space_dimension": int(basis_meta.group(3)),
                "leads": leads,
            }
        ),
        "specialized_generators": [
            generator
            for _, generator in re.findall(
                r"(?m)^SPECIALIZED_GENERATOR (\d+) (.*)$",
                completed.stdout,
            )
        ],
        "solver_seconds": round(solver_seconds, 6),
        "program_bytes": len(program),
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        "singular_output": (
            None if transcript_path is None else str(transcript_path)
        ),
        "seconds": round(time.monotonic() - started, 6),
        "reproduction_command": " ".join(sys.argv),
    }
    moment_meta = re.search(
        r"(?m)^MOMENT_META (\d+) (\d+) (-?\d+) (-?\d+)$",
        completed.stdout,
    )
    if moment_meta is not None:
        payload["moment_normal_forms"] = {
            "mu6_terms": int(moment_meta.group(1)),
            "mu7_terms": int(moment_meta.group(2)),
            "mu6_total_degree": int(moment_meta.group(3)),
            "mu7_total_degree": int(moment_meta.group(4)),
        }
    unit_meta = re.search(
        r"(?m)^UNIT_META (\d+) ([01])$",
        completed.stdout,
    )
    if unit_meta is not None:
        payload["unit_test"] = {
            "basis_size": int(unit_meta.group(1)),
            "unit_ideal": unit_meta.group(2) == "1",
        }
    pencil_meta = re.search(
        r"(?m)^PENCIL_META (\d+) (\d+) (\d+) (\d+)$",
        completed.stdout,
    )
    pencil_values = [
        (int(node), int(value) % arguments.prime)
        for node, value in re.findall(
            r"(?m)^PENCIL_VALUE (\d+) (-?\d+)$",
            completed.stdout,
        )
    ]
    if arguments.stage == "pencil":
        if (
            pencil_meta is None
            or len(pencil_values) != int(pencil_meta.group(1)) + 1
        ):
            raise RuntimeError(transcript[-12000:])
        coefficients = interpolate_mod_prime(
            pencil_values,
            arguments.prime,
        )
        for node, value in pencil_values:
            reconstructed = sum(
                coefficient * pow(node, degree, arguments.prime)
                for degree, coefficient in enumerate(coefficients)
            ) % arguments.prime
            if reconstructed != value:
                raise AssertionError(
                    "determinant-pencil interpolation did not replay "
                    f"the value at z={node}"
                )
        quotient_dimension = int(pencil_meta.group(1))
        if (coefficients[0] != 0) != (
            int(pencil_meta.group(2)) == quotient_dimension
        ):
            raise AssertionError("constant pencil coefficient/rank mismatch")
        if (coefficients[-1] != 0) != (
            int(pencil_meta.group(3)) == quotient_dimension
        ):
            raise AssertionError("leading pencil coefficient/rank mismatch")
        payload["multiplication_pencil"] = {
            "quotient_dimension": quotient_dimension,
            "rank_mu6": int(pencil_meta.group(2)),
            "rank_mu7": int(pencil_meta.group(3)),
            "joint_block_rank": int(pencil_meta.group(4)),
            "evaluation_values": [
                {"z": node, "determinant": value}
                for node, value in pencil_values
            ],
            "coefficients_ascending": coefficients,
            "nonzero_coefficient_indices": [
                index
                for index, coefficient in enumerate(coefficients)
                if coefficient
            ],
            "identically_zero": not any(coefficients),
        }

    output = arguments.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

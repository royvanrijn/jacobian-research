#!/usr/bin/env python3
"""Exact bidegree-(2,1) exceptional-transport pilot for Fermigier.

The calculation fixes the independently exceptional pair P13 at the rank-22
anchor and R20E1 at the rank-20 anchor.  It studies the complete finite chart

    x(T) = (a(c) T + b(c) + k(T-T_22)(T-T_20)) / (cT+1),

where a(c), b(c) are forced by the two endpoint values.  Singular performs the
exact factorization of the two-parameter T-discriminant over QQ[c,k].  SymPy
then classifies every rational linear component of that discriminant by the
squareclass kernel of the specialized T-polynomial.

This is a bounded structural pilot, not a score or specialization sweep.  It
does not classify rational points on the remaining irreducible degree-32
discriminant component.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
from typing import Any

import sympy as sp

from analyze_fermigier_exceptional_transport import (
    E22_T,
    E22_U,
    RANK20_T,
    RANK20_U,
    e22_exceptional_quotient,
    factor_signature,
    family_expression,
    hyperelliptic_genus,
    rank20_exceptional_quotient,
    rational_text,
    squareclass_kernel_degree,
)


Q = Fraction
LEFT_LABEL = "P13"
RIGHT_LABEL = "R20E1"
AFFINE_CANCELLATION_SLOPE = Q(-5899690, 732683)
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic-curves/"
    "elliptic_fermigier_bidegree21_p13_r20e1_pilot.json"
)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def pencil_polynomial(
    left_x_q: Fraction, right_x_q: Fraction
) -> tuple[sp.Poly, dict[str, sp.Expr]]:
    """Return the primitive cleared polynomial in QQ[T,c,k]."""

    T, X, c, k = sp.symbols("T X c k")
    first = sp.Rational(E22_T.numerator, E22_T.denominator)
    second = sp.Rational(RANK20_T.numerator, RANK20_T.denominator)
    left_x = sp.Rational(left_x_q.numerator, left_x_q.denominator)
    right_x = sp.Rational(right_x_q.numerator, right_x_q.denominator)
    a = sp.cancel(
        (
            left_x * (c * first + 1)
            - right_x * (c * second + 1)
        )
        / (first - second)
    )
    b = sp.cancel(left_x * (c * first + 1) - a * first)
    denominator = c * T + 1
    numerator = sp.expand(a * T + b + k * (T - first) * (T - second))
    cleared = sp.together(
        family_expression(T, X).subs(X, numerator / denominator)
        * denominator**4
    ).as_numer_denom()[0]
    polynomial = sp.Poly(cleared, T, c, k, domain=sp.QQ).primitive()[1]
    if polynomial.LC() < 0:
        polynomial = -polynomial
    return polynomial, {
        "T": T,
        "c": c,
        "k": k,
        "a": a,
        "b": b,
        "numerator": numerator,
        "denominator": denominator,
    }


def _singular_number(value: Fraction) -> str:
    return rational_text(value)


def singular_factorization(
    left_x: Fraction, right_x: Fraction, *, timeout: int = 300
) -> dict[str, Any]:
    """Factor the exact discriminant with Singular and return stable metadata."""

    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the bivariate factorization")
    program = f"""
ring r=0,(c,k,T),dp;
number t1={_singular_number(E22_T)};
number t2={_singular_number(RANK20_T)};
number x1={_singular_number(left_x)};
number x2={_singular_number(right_x)};
poly a=(x1*(c*t1+1)-x2*(c*t2+1))/(t1-t2);
poly b=x1*(c*t1+1)-a*t1;
poly D=c*T+1;
poly N=a*T+b+k*(T-t1)*(T-t2);
poly f0=18103855887324900+102302344648*T^2-879500*T^4+T^6;
poly f1=-257843832010380-650709150*T^2+1860*T^4;
poly f2=1195214262641+1718550*T^2-2*T^4;
poly f3=-2051321790-1860*T^2;
poly f4=1149050+T^2;
poly P=f0*D^4+f1*N*D^3+f2*N^2*D^2+f3*N^3*D+f4*N^4;
poly rr=resultant(P,diff(P,T),T);
// rr=(-1)^45*LC_T(P)*Disc_T(P).  The omitted rational scalar is a unit.
poly dd=rr/((c-k)^2*(c+k)^2);
poly known=(k-c)^2*(k+c)^2*(28917c+10)^12*(39508c+39)^12
           *(5899690c+732683k)^12;
poly nonlinear=dd/known;
if (dd-known*nonlinear!=0) {{ "@@ERROR_DIVISION"; }}
list L=factorize(nonlinear);
int nonconstant=0;
int idx=0;
for (idx=1;idx<=size(L[1]);idx++) {{
  if (deg(L[1][idx])>0) {{ nonconstant=nonconstant+1; }}
}}
intvec wc=1,0,0;
intvec wk=0,1,0;
"@@RESULTANT_TERMS"; size(rr);
"@@RESULTANT_TOTAL_DEGREE"; deg(rr);
"@@DISCRIMINANT_TERMS"; size(dd);
"@@DISCRIMINANT_TOTAL_DEGREE"; deg(dd);
"@@NONLINEAR_NONCONSTANT_FACTOR_COUNT"; nonconstant;
"@@NONLINEAR_TOTAL_DEGREE"; deg(nonlinear);
"@@NONLINEAR_C_DEGREE"; deg(nonlinear,wc);
"@@NONLINEAR_K_DEGREE"; deg(nonlinear,wk);
"@@NONLINEAR_TERMS"; size(nonlinear);
for (idx=1;idx<=size(L[1]);idx++) {{
  if (deg(L[1][idx])>0) {{
    "@@NONLINEAR_EXPONENT"; L[2][idx];
    "@@NONLINEAR_BEGIN"; L[1][idx]; "@@NONLINEAR_END";
  }}
}}
"""
    completed = subprocess.run(
        [executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"Singular wrote to stderr: {completed.stderr.strip()}")
    output = completed.stdout
    if "@@ERROR_DIVISION" in output or "   ?" in output:
        raise AssertionError(f"Singular calculation failed:\n{output}")

    def integer_after(marker: str) -> int:
        match = re.search(rf"^{re.escape(marker)}\n(-?\d+)$", output, re.MULTILINE)
        if match is None:
            raise AssertionError(f"missing Singular marker {marker}")
        return int(match.group(1))

    factor_match = re.search(
        r"^@@NONLINEAR_BEGIN\n(.*?)\n@@NONLINEAR_END$",
        output,
        re.MULTILINE | re.DOTALL,
    )
    if factor_match is None:
        raise AssertionError("missing nonlinear discriminant factor")
    nonlinear_text = factor_match.group(1).strip()
    answer = {
        "resultant_term_count": integer_after("@@RESULTANT_TERMS"),
        "resultant_total_degree": integer_after("@@RESULTANT_TOTAL_DEGREE"),
        "discriminant_term_count": integer_after("@@DISCRIMINANT_TERMS"),
        "discriminant_total_degree": integer_after("@@DISCRIMINANT_TOTAL_DEGREE"),
        "nonlinear_nonconstant_factor_count": integer_after(
            "@@NONLINEAR_NONCONSTANT_FACTOR_COUNT"
        ),
        "nonlinear_total_degree": integer_after("@@NONLINEAR_TOTAL_DEGREE"),
        "nonlinear_degree_in_c": integer_after("@@NONLINEAR_C_DEGREE"),
        "nonlinear_degree_in_k": integer_after("@@NONLINEAR_K_DEGREE"),
        "nonlinear_term_count": integer_after("@@NONLINEAR_TERMS"),
        "nonlinear_exponent": integer_after("@@NONLINEAR_EXPONENT"),
        "nonlinear_primitive_factor_sha256": sha256_text(nonlinear_text + "\n"),
    }
    expected = {
        "resultant_term_count": 2493,
        "resultant_total_degree": 76,
        "discriminant_total_degree": 72,
        "nonlinear_nonconstant_factor_count": 1,
        "nonlinear_total_degree": 32,
        "nonlinear_exponent": 1,
    }
    for key, value in expected.items():
        if answer[key] != value:
            raise AssertionError(
                f"Singular {key} changed: {answer[key]} != {value}"
            )
    return answer


def component_record(
    polynomial: sp.Poly,
    symbols: dict[str, sp.Expr],
    *,
    label: str,
    substitutions: dict[sp.Symbol, sp.Expr],
    validity: str,
    interpretation: str,
) -> dict[str, Any]:
    T = symbols["T"]
    remaining = symbols["k"] if symbols["c"] in substitutions else symbols["c"]
    specialized = sp.Poly(
        sp.expand(polynomial.as_expr().subs(substitutions)),
        T,
        domain=sp.QQ.frac_field(remaining),
    ).primitive()[1]
    kernel_degree = squareclass_kernel_degree(specialized)
    return {
        "component": label,
        "validity": validity,
        "T_degree_after_specialization": specialized.degree(),
        "factor_signature_over_parameter_function_field": [
            list(item) for item in factor_signature(specialized)
        ],
        "squareclass_kernel_degree": kernel_degree,
        "squareclass_genus": hyperelliptic_genus(kernel_degree),
        "genus_at_most_one": kernel_degree <= 4,
        "interpretation": interpretation,
    }


def run(root: Path, *, timeout: int = 300) -> dict[str, Any]:
    e22_record, e22_points = e22_exceptional_quotient()
    rank20_record, rank20_points = rank20_exceptional_quotient(root)
    left_x = dict(e22_points)[LEFT_LABEL]
    right_x = dict(rank20_points)[RIGHT_LABEL]
    polynomial, symbols = pencil_polynomial(left_x, right_x)
    T, c, k = symbols["T"], symbols["c"], symbols["k"]
    generic = sp.Poly(
        polynomial.as_expr(), T, domain=sp.QQ.frac_field(c, k)
    )
    leading = sp.factor(generic.LC())
    expected_leading = sp.factor(
        sp.Rational(6666883836179368888567109693610000)
        * (c - k) ** 2
        * (c + k) ** 2
    )
    if leading != expected_leading:
        raise AssertionError("the generic leading coefficient changed")
    singular = singular_factorization(left_x, right_x, timeout=timeout)

    components = [
        component_record(
            polynomial,
            symbols,
            label="k-c=0",
            substitutions={k: c},
            validity="valid degree-drop locus",
            interpretation=(
                "The leading T^10 and T^9 terms vanish because x(T)~T at infinity; "
                "the remaining octic is squarefree over QQ(c)."
            ),
        ),
        component_record(
            polynomial,
            symbols,
            label="k+c=0",
            substitutions={k: -c},
            validity="valid degree-drop locus",
            interpretation=(
                "The leading T^10 and T^9 terms vanish because x(T)~-T at infinity; "
                "the remaining octic is squarefree over QQ(c)."
            ),
        ),
        component_record(
            polynomial,
            symbols,
            label="5899690*c+732683*k=0",
            substitutions={k: sp.Rational(AFFINE_CANCELLATION_SLOPE.numerator, AFFINE_CANCELLATION_SLOPE.denominator) * c},
            validity="valid but representation is non-reduced",
            interpretation=(
                "The numerator is divisible by c*T+1.  After cancellation this is "
                "the unique already-classified affine P13-to-R20E1 transport."
            ),
        ),
        component_record(
            polynomial,
            symbols,
            label="28917*c+10=0",
            substitutions={c: -sp.Rational(10, 28917)},
            validity="invalid at rank20 anchor",
            interpretation=(
                "Numerator and denominator both vanish at T=28917/10, so the "
                "uncancelled expression does not retain the declared endpoint value."
            ),
        ),
        component_record(
            polynomial,
            symbols,
            label="39508*c+39=0",
            substitutions={c: -sp.Rational(39, 39508)},
            validity="invalid at E22 anchor",
            interpretation=(
                "Numerator and denominator both vanish at T=39508/39, so the "
                "uncancelled expression does not retain the declared endpoint value."
            ),
        ),
    ]
    if any(item["genus_at_most_one"] for item in components):
        raise AssertionError("an unexpected low-genus rational component appeared")

    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-bidegree21-pilot.v1",
        "status": "complete exact bounded structural pilot",
        "claim_level": "exact computation; no new section or specialization found",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "family_id": "fermigier-mestre-v1",
        "anchors": {
            "E22": {
                "canonical_parameter": f"u={rational_text(E22_U)}",
                "alias": f"literal_shift_T={rational_text(E22_T)}",
                "certified_rank_lower_bound": 22,
                "exceptional_quotient_rank_lower_bound": e22_record[
                    "exceptional_quotient_rank_lower_bound"
                ],
            },
            "rank20": {
                "canonical_parameter": f"u={rational_text(RANK20_U)}",
                "alias": f"literal_shift_T={rational_text(RANK20_T)}",
                "certified_rank_lower_bound": 20,
                "exceptional_quotient_rank_lower_bound": rank20_record[
                    "exceptional_quotient_rank_lower_bound"
                ],
            },
        },
        "representative_pair": {
            "left_label": LEFT_LABEL,
            "left_x": rational_text(left_x),
            "right_label": RIGHT_LABEL,
            "right_x": rational_text(right_x),
            "both_independent_modulo_generic_sections": True,
        },
        "pencil": {
            "definition": (
                "x(T)=(a(c)*T+b(c)+k*(T-39508/39)*(T-28917/10))/(c*T+1), "
                "with a(c),b(c) fixed by x(39508/39)=256714/39 and "
                "x(28917/10)=-8545"
            ),
            "chart": "finite denominator chart d=1",
            "generic_cleared_T_degree": generic.degree(),
            "generic_squareclass_genus": hyperelliptic_genus(generic.degree()),
            "leading_coefficient_factorization": (
                "6666883836179368888567109693610000*(c-k)^2*(c+k)^2"
            ),
        },
        "discriminant_factorization_over_QQ_c_k": {
            **singular,
            "normalization": "up to a nonzero rational scalar",
            "factors": [
                {"factor": "irreducible nonlinear factor", "total_degree": 32, "exponent": 1},
                {"factor": "k-c", "total_degree": 1, "exponent": 2},
                {"factor": "k+c", "total_degree": 1, "exponent": 2},
                {"factor": "28917*c+10", "total_degree": 1, "exponent": 12},
                {"factor": "39508*c+39", "total_degree": 1, "exponent": 12},
                {"factor": "5899690*c+732683*k", "total_degree": 1, "exponent": 12},
            ],
        },
        "rational_linear_components": components,
        "nonlinear_component": {
            "irreducible_over_QQ": True,
            "discriminant_multiplicity": 1,
            "generic_behavior": (
                "one simple finite double root, hence squareclass kernel degree 8 "
                "and genus 3 at its generic point"
            ),
            "generic_genus_at_most_one": False,
            "rational_special_points_and_component_intersections_classified": False,
        },
        "scope": {
            "completed_pairs": [f"{LEFT_LABEL}x{RIGHT_LABEL}"],
            "completed_pair_count": 1,
            "possible_independent_pair_count": 80,
            "all_pairs_classified": False,
            "infinity_denominator_chart": (
                "already included in the complete Mobius classification at k=0; "
                "not separately extended to nonzero curvature in this pilot"
            ),
            "not_claimed": (
                "No claim is made about rational special points on the irreducible "
                "degree-32 component, its intersections, or the other 79 pairs."
            ),
        },
        "outcome": {
            "rational_linear_components": 5,
            "valid_new_rational_components": 2,
            "genus_at_most_one_components": 0,
            "new_sections": 0,
            "new_specializations": 0,
            "target_met": False,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "singular": subprocess.run(
                [shutil.which("Singular") or "Singular", "--version"],
                text=True,
                capture_output=True,
                check=True,
            ).stdout.splitlines()[0],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/analyze_fermigier_bidegree21_pilot.py"
        ),
    }
    stable = dict(artifact)
    stable.pop("generated_at_utc")
    artifact["result_sha256"] = sha256_text(
        json.dumps(stable, sort_keys=True, separators=(",", ":"))
    )
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=root / OUTPUT_RELATIVE)
    args = parser.parse_args()
    artifact = run(root, timeout=args.timeout)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Scout rational q-fibres of the descended Hurwitz base over GF(p).

This is a bounded modular routing calculation, not a characteristic-zero
certificate.  It constructs the already certified base equations

    H = K = J6 = J7 = 0

directly in the quotient module, specializes every q in GF(p), and tests
the open b0*P1*L4*S2*R1 != 0 by adjoining one inverse variable.  Any
surviving fibre is reported without attempting to interpret or lift it.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess

from flint import nmod_mpoly_ctx

import verify_two_pair_sic_bidegree33_rank_two_hurwitz_module_descent as descent


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_hurwitz_base_fibres_p29.json"
)


def specialize_to_two_variables(polynomial, value: int, prime: int, context):
    coefficients: dict[tuple[int, int], int] = {}
    for (b0_power, q_power, lambda_power), coefficient in polynomial.to_dict().items():
        exponent = (b0_power, lambda_power)
        updated = (
            coefficients.get(exponent, 0)
            + int(coefficient) * pow(value, q_power, prime)
        ) % prime
        if updated:
            coefficients[exponent] = updated
        else:
            coefficients.pop(exponent, None)
    return context.from_dict(coefficients)


def singular_fibre(
    value: int,
    prime: int,
    polynomials: dict[str, object],
    localizer,
    timeout: int,
    certificate: bool = False,
    saturate: bool = True,
) -> dict[str, object]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required")
    context = nmod_mpoly_ctx.get(["b0", "lambda"], prime)
    specialized = {
        name: specialize_to_two_variables(polynomial, value, prime, context)
        for name, polynomial in polynomials.items()
    }
    specialized_localizer = specialize_to_two_variables(
        localizer, value, prime, context
    )
    definitions = [
        f"poly {name}={polynomial};"
        for name, polynomial in specialized.items()
    ]
    equation_names = ",".join(specialized)
    ring_variables = "b0,lambda,s" if saturate else "b0,lambda"
    ideal_generators = equation_names + (",s*O-1" if saturate else "")
    source_lines = [
            f"ring R={prime},({ring_variables}),dp;",
            "option(redSB);",
            *definitions,
            f"poly O={specialized_localizer};",
            f"ideal I={ideal_generators};",
            "ideal G=std(I);",
            'print("FIBRE_BEGIN");',
            'print("basis_size="+string(size(G)));',
            'print("first_degree="+string(deg(G[1])));',
            'print("is_unit="+string(G[1]==1));',
            'print("quotient_length="+string(vdim(G)));',
            'print("FIBRE_END");',
            'print("BASIS_BEGIN");',
            "print(G);",
            'print("BASIS_END");',
    ]
    if certificate:
        if not saturate:
            source_lines.extend(
                [
                    "poly NP=1;",
                    "int nilpotency_index=0;",
                    "for (int n=1; n<=400; n++)",
                    "{",
                    "  NP=reduce(NP*O,G);",
                    "  if (NP==0)",
                    "  {",
                    "    nilpotency_index=n;",
                    "    break;",
                    "  }",
                    "}",
                    'print("NILPOTENCY_INDEX="+string(nilpotency_index));',
                    "if (nilpotency_index==1)",
                    "{",
                    "  matrix U=lift(I,ideal(O));",
                    '  print("LOCALIZER_LIFT_BEGIN");',
                    "  for (int i=1; i<=nrows(U); i++)",
                    "  {",
                    '    print("entry="+string(i)+",degree="'
                    '+string(deg(U[i,1]))+",terms="+string(size(U[i,1])));',
                    "  }",
                    '  print("LOCALIZER_LIFT_END");',
                    "}",
                ]
            )
        source_lines.extend(
            [
                "if (G[1]==1)",
                "{",
                "  matrix T=lift(I,G);",
                '  print("LIFT_BEGIN");',
                "  for (int i=1; i<=nrows(T); i++)",
                "  {",
                '    print("entry="+string(i)+",degree="'
                '+string(deg(T[i,1]))+",terms="+string(size(T[i,1])));',
                "  }",
                '  print("LIFT_END");',
                "}",
            ]
        )
    source_lines.append("quit;")
    source = "\n".join(source_lines)
    try:
        completed = subprocess.run(
            [executable, "-q"],
            input=source,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "q": value,
            "status": "timeout",
            "timeout_seconds": timeout,
            "source_sha256": sha256(source.encode()).hexdigest(),
        }
    if completed.returncode != 0:
        raise AssertionError(
            f"Singular failed on q={value}:\n"
            + completed.stdout[-2000:]
            + completed.stderr[-2000:]
        )
    match = re.search(
        r"FIBRE_BEGIN\s+basis_size=(\d+)\s+first_degree=(-?\d+)\s+"
        r"is_unit=(\d+)\s+quotient_length=(-?\d+)\s+FIBRE_END",
        completed.stdout,
    )
    if match is None:
        raise AssertionError(
            f"could not parse Singular output on q={value}:\n"
            + completed.stdout[-2000:]
        )
    basis_size, first_degree, is_unit, quotient_length = map(int, match.groups())
    basis_match = re.search(
        r"BASIS_BEGIN\s+(.*?)\s+BASIS_END",
        completed.stdout,
        re.DOTALL,
    )
    if basis_match is None:
        raise AssertionError("could not parse the specialized basis")
    answer = {
        "q": value,
        "status": "unit" if is_unit else "survivor",
        "basis_size": basis_size,
        "first_degree": first_degree,
        "quotient_length": quotient_length,
        "source_sha256": sha256(source.encode()).hexdigest(),
        "specialized_terms": {
            name: len(polynomial.to_dict())
            for name, polynomial in specialized.items()
        },
        "localizer_terms": len(specialized_localizer.to_dict()),
        "basis": basis_match.group(1).strip() if not is_unit else "_[1]=1",
    }
    lift_match = re.search(
        r"LIFT_BEGIN\s+(.*?)\s+LIFT_END", completed.stdout, re.DOTALL
    )
    if lift_match is not None:
        answer["lift_profile"] = [
            {
                "entry": int(entry),
                "degree": int(degree),
                "terms": int(terms),
            }
            for entry, degree, terms in re.findall(
                r"entry=(\d+),degree=(-?\d+),terms=(\d+)",
                lift_match.group(1),
            )
        ]
    nilpotency_match = re.search(
        r"NILPOTENCY_INDEX=(\d+)", completed.stdout
    )
    if nilpotency_match is not None:
        answer["localizer_nilpotency_index"] = int(nilpotency_match.group(1))
    localizer_lift_match = re.search(
        r"LOCALIZER_LIFT_BEGIN\s+(.*?)\s+LOCALIZER_LIFT_END",
        completed.stdout,
        re.DOTALL,
    )
    if localizer_lift_match is not None:
        answer["localizer_lift_profile"] = [
            {
                "entry": int(entry),
                "degree": int(degree),
                "terms": int(terms),
            }
            for entry, degree, terms in re.findall(
                r"entry=(\d+),degree=(-?\d+),terms=(\d+)",
                localizer_lift_match.group(1),
            )
        ]
    return answer


def construct_base(prime: int, through: int):
    later_orders = tuple(range(6, through + 1))
    modular = descent.modular_descent(prime, later_orders=later_orders)
    context = nmod_mpoly_ctx.get(["b0", "q", "lambda"], prime)
    h_polynomial = descent.flint_coefficient_polynomial(
        modular["norm"], context
    )
    linear_norm = descent.flint_linear_norm_descent(
        modular["descended_mu4"],
        modular["mu5_subresultant"],
        context,
    )
    later = {
        order: descent.flint_later_base_equation(
            modular[f"descended_mu{order}"],
            modular["descended_mu4"],
            modular["mu5_subresultant"],
            linear_norm,
            context,
            order - 5 if order <= 7 else None,
        )
        for order in later_orders
    }
    cubic_lead = descent.flint_coefficient_polynomial(
        descent.plain_coefficient_groups(
            modular["descended_mu4"], 2
        )[3],
        context,
    )
    quadratic_lead = descent.flint_coefficient_polynomial(
        descent.plain_coefficient_groups(
            modular["mu5_subresultant"], 2
        )[2],
        context,
    )
    birational_coordinate = descent.flint_coefficient_polynomial(
        modular["degree_drop_V"], context
    )
    b0, q, parameter = context.gens()
    first_pivot = 3 * parameter + 3 + 4 * q
    localizer = (
        b0
        * first_pivot
        * cubic_lead
        * quadratic_lead
        * birational_coordinate
        * linear_norm["remainder_linear"]
    )
    return (
        {
            "H": h_polynomial,
            "K": linear_norm["residual"],
            **{
                f"J{order}": later[order]["residual"]
                for order in later_orders
            },
        },
        localizer,
        {
            str(order): {
                "removed_P1_power": modular[
                    "later_removed_P1_powers"
                ][str(order)],
                "removed_determinant_factors": later[order][
                    "removed_factors"
                ],
            }
            for order in later_orders
        },
    )


def p29_small_lift_rejection() -> dict[str, object]:
    values = [
        Fraction(-12),   # lambda
        Fraction(-7),    # a1
        Fraction(5, 8),  # a2=z/b0
        Fraction(0),     # a3, solved below
        Fraction(-8),    # b0
        Fraction(56),    # b1=b0*q
        Fraction(-672),  # b2 from A=0
    ]
    moments = descent.exact_moment_polynomials((2, 3))

    def evaluate(polynomial) -> Fraction:
        return sum(
            Fraction(coefficient)
            * reduce(
                lambda product, item: product
                * values[item[0]] ** item[1],
                enumerate(exponent),
                Fraction(1),
            )
            for exponent, coefficient in polynomial.items()
        )

    pivot, rest = descent.split_linear_rational_polynomial(
        descent.rational_polynomial(moments[2]), 3
    )
    pivot_value = evaluate(pivot)
    values[3] = -evaluate(rest) / pivot_value
    mu2 = evaluate(descent.rational_polynomial(moments[2]))
    mu3 = evaluate(descent.rational_polynomial(moments[3]))
    assert mu2 == 0 and mu3 != 0
    return {
        "symmetric_base_lift": {
            "q": -7,
            "b0": -8,
            "lambda": -12,
            "z": -5,
        },
        "reconstructed_coefficients": {
            "a1": str(values[1]),
            "a2": str(values[2]),
            "a3": str(values[3]),
            "b0": str(values[4]),
            "b1": str(values[5]),
            "b2": str(values[6]),
        },
        "mu2": str(mu2),
        "mu3": str(mu3),
        "conclusion": "the symmetric small lift is not a QQ moment-zero point",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=29)
    parser.add_argument("--through", type=int, default=7)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--emit-msolve", type=Path)
    parser.add_argument("--emit-singular-parametric", type=Path)
    parser.add_argument("--emit-singular-membership", type=Path)
    parser.add_argument(
        "--membership-generators",
        type=int,
        help="use only the first N base equations in the membership source",
    )
    parser.add_argument("--emit-unsaturated", action="store_true")
    parser.add_argument("--skip-fibres", action="store_true")
    parser.add_argument("--certificate-fibres", action="store_true")
    parser.add_argument("--fibre-unsaturated", action="store_true")
    parser.add_argument(
        "--q-values",
        help="comma-separated rational q-fibres (default: every value)",
    )
    args = parser.parse_args()
    if args.prime <= 3 * args.through:
        raise ValueError("the prime must exceed three times --through")
    if args.through < 7:
        raise ValueError("--through must be at least 7")

    polynomials, localizer, factor_removal = construct_base(
        args.prime, args.through
    )
    if args.emit_msolve is not None:
        expressions = [str(polynomial) for polynomial in polynomials.values()]
        if not args.emit_unsaturated:
            expressions.append(str(localizer))
        source = "\n".join(
            [
                "b0,lambda,q",
                str(args.prime),
                ",\n".join(expressions),
            ]
        )
        args.emit_msolve.parent.mkdir(parents=True, exist_ok=True)
        args.emit_msolve.write_text(source + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "emitted": str(args.emit_msolve),
                    "source_sha256": sha256(source.encode()).hexdigest(),
                },
                sort_keys=True,
            )
        )
    if args.emit_singular_parametric is not None:
        definitions = [
            f"poly {name}={polynomial};"
            for name, polynomial in polynomials.items()
        ]
        equation_names = ",".join(polynomials)
        ring_variables = "b0,lambda" if args.emit_unsaturated else "b0,lambda,s"
        ideal_generators = equation_names
        if not args.emit_unsaturated:
            ideal_generators += ",s*O-1"
        source_lines = [
                'LIB "polylib.lib";',
                f"ring R=({args.prime},q),({ring_variables}),dp;",
                "option(redSB);",
                *definitions,
        ]
        if not args.emit_unsaturated:
            source_lines.append(f"poly O={localizer};")
        source = "\n".join(
            [
                *source_lines,
                f"ideal I={ideal_generators};",
                "ideal G=slimgb(I);",
                'print("BASIS_BEGIN");',
                "print(G);",
                'print("BASIS_END");',
                "if (G[1]==1)",
                "{",
                "  matrix T=lift(I,G);",
                "  number E=1;",
                "  number d;",
                "  for (int i=1; i<=nrows(T); i++)",
                "  {",
                "    for (int j=1; j<=ncols(T); j++)",
                "    {",
                "      if (T[i,j]!=0)",
                "      {",
                "        d=denominator(content(T[i,j]));",
                "        E=E*d/gcd(E,d);",
                "      }",
                "    }",
                "  }",
                '  print("DENOMINATOR_BEGIN");',
                "  print(E);",
                '  print("DENOMINATOR_END");',
                "}",
                "quit;",
            ]
        )
        args.emit_singular_parametric.parent.mkdir(parents=True, exist_ok=True)
        args.emit_singular_parametric.write_text(source + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "emitted": str(args.emit_singular_parametric),
                    "source_sha256": sha256(source.encode()).hexdigest(),
                },
                sort_keys=True,
            )
        )
    if args.emit_singular_membership is not None:
        membership_names = list(polynomials)
        if args.membership_generators is not None:
            if not 1 <= args.membership_generators <= len(membership_names):
                raise ValueError("--membership-generators is out of range")
            membership_names = membership_names[: args.membership_generators]
        definitions = [
            f"poly {name}={polynomial};"
            for name, polynomial in polynomials.items()
            if name in membership_names
        ]
        equation_names = ",".join(membership_names)
        source = "\n".join(
            [
                'LIB "polylib.lib";',
                f"ring R=({args.prime},q),(b0,lambda),dp;",
                *definitions,
                f"poly O={localizer};",
                f"ideal I={equation_names};",
                "matrix T=lift(I,ideal(O));",
                "poly check=O;",
                "for (int i=1; i<=nrows(T); i++)",
                "{",
                "  check=check-T[i,1]*I[i];",
                "}",
                'print("CHECK_BEGIN");',
                "print(check);",
                'print("CHECK_END");',
                "number E=1;",
                "number d;",
                "for (int i=1; i<=nrows(T); i++)",
                "{",
                "  if (T[i,1]!=0)",
                "  {",
                "    d=denominator(content(T[i,1]));",
                "    E=E*d/gcd(E,d);",
                "  }",
                "}",
                'print("DENOMINATOR_BEGIN");',
                "print(E);",
                'print("DENOMINATOR_END");',
                'print("LIFT_PROFILE_BEGIN");',
                "for (int i=1; i<=nrows(T); i++)",
                "{",
                '  print("entry="+string(i)+",degree="'
                '+string(deg(T[i,1]))+",terms="+string(size(T[i,1])));',
                "}",
                'print("LIFT_PROFILE_END");',
                "quit;",
            ]
        )
        args.emit_singular_membership.parent.mkdir(parents=True, exist_ok=True)
        args.emit_singular_membership.write_text(source + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "emitted": str(args.emit_singular_membership),
                    "source_sha256": sha256(source.encode()).hexdigest(),
                },
                sort_keys=True,
            )
        )
    if args.skip_fibres:
        return
    q_values = (
        range(args.prime)
        if args.q_values is None
        else [int(value) % args.prime for value in args.q_values.split(",")]
    )
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        fibres = list(
            executor.map(
                lambda value: singular_fibre(
                    value,
                    args.prime,
                    polynomials,
                    localizer,
                    args.timeout,
                    args.certificate_fibres,
                    not args.fibre_unsaturated,
                ),
                q_values,
            )
        )
    counts = {
        status: sum(fibre["status"] == status for fibre in fibres)
        for status in ("unit", "survivor", "timeout")
    }
    artifact = {
        "format": "two-pair-sic-bidegree33-rank-two-hurwitz-base-fibres-v1",
        "field": f"GF({args.prime})",
        "equations": list(polynomials),
        "open_localizer": "b0*P1*lc_z(N4)*lc_z(S5)*V*R1",
        "prime_guard": f"{args.prime}>3*{args.through}",
        "polynomial_profiles": {
            name: descent.flint_profile(polynomial)
            for name, polynomial in polynomials.items()
        },
        "descent_factor_removal": factor_removal,
        "counts": counts,
        "fibres": fibres,
        "conclusion": (
            "bounded exact modular routing only; unit fibres exclude "
            "geometric points with that GF(p)-rational q-coordinate on "
            "the displayed open, but do not exclude q in extension "
            "fields or imply a characteristic-zero unit ideal"
        ),
    }
    if args.prime == 29:
        artifact["small_lift_rejection"] = p29_small_lift_rejection()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(counts, sort_keys=True))


if __name__ == "__main__":
    main()

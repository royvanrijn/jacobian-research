#!/usr/bin/env python3
"""Research the remaining squarefree quartic Rabinowitsch membership.

This script reconstructs the ideal

  I=(f3,f4,f5,f6,z*p*(8*c-3*d^2)-1)

over Q and asks Singular's modular standard-basis routine for a candidate
basis.  It then checks both ideal containments exactly:

* every generator of I reduces to zero by the candidate basis; and
* ``lift(I,G)`` expresses every candidate generator as an exact
  Q[z,c,d,lambda]-linear combination of the original generators.

Only the second check turns modular reconstruction into a characteristic-zero
membership certificate.
"""

from __future__ import annotations

from math import factorial, gcd, isqrt
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e, lam = sp.symbols("a b c d e lam")
P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4
A = u * v * (u - v) * (u - lam * v)


def singular_support(expression: str) -> set[tuple[int, int, int]]:
    """Parse only the monomial support of a Singular polynomial."""
    compact = re.sub(r"\s+", "", expression)
    if compact == "0":
        return set()
    terms = compact.replace("-", "+-").split("+")
    support: set[tuple[int, int, int]] = set()
    variable_index = {"c": 0, "d": 1, "lam": 2}
    for term in terms:
        if not term:
            continue
        exponents = [0, 0, 0]
        for factor in term.lstrip("-").split("*"):
            if factor.lstrip("-").isdigit():
                continue
            if "^" in factor:
                variable, exponent_text = factor.split("^", 1)
                exponent = int(exponent_text)
            else:
                variable = factor
                exponent = 1
            exponents[variable_index[variable]] += exponent
        support.add(tuple(exponents))
    return support


def singular_polynomial(
    expression: str,
    prime: int,
) -> dict[tuple[int, int, int], int]:
    """Parse a printed Singular F_p[c,d,lam] polynomial."""
    compact = re.sub(r"\s+", "", expression)
    if compact == "0":
        return {}
    terms = compact.replace("-", "+-").split("+")
    result: dict[tuple[int, int, int], int] = {}
    variable_index = {"c": 0, "d": 1, "lam": 2}
    for term in terms:
        if not term:
            continue
        coefficient = -1 if term.startswith("-") else 1
        exponents = [0, 0, 0]
        for factor in term.lstrip("-").split("*"):
            if factor.isdigit():
                coefficient *= int(factor)
                continue
            if "^" in factor:
                variable, exponent_text = factor.split("^", 1)
                exponent = int(exponent_text)
            else:
                variable = factor
                exponent = 1
            exponents[variable_index[variable]] += exponent
        monomial = tuple(exponents)
        result[monomial] = (result.get(monomial, 0) + coefficient) % prime
    return {monomial: value for monomial, value in result.items() if value}


def rational_reconstruct(
    residue: int,
    modulus: int,
) -> tuple[int, int] | None:
    """Return the balanced rational reconstruction, when it exists."""
    bound = isqrt(modulus // 2)
    old_remainder, remainder = modulus, residue
    old_coefficient, coefficient = 0, 1
    while remainder > bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_coefficient, coefficient = (
            coefficient,
            old_coefficient - quotient * coefficient,
        )
    numerator, denominator = remainder, coefficient
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        denominator == 0
        or denominator > bound
        or abs(numerator) > bound
        or gcd(numerator, denominator) != 1
        or (numerator - residue * denominator) % modulus
    ):
        return None
    return numerator, denominator


def raw_lift_record(
    singular: str,
    prime: int,
    exponent: int,
    eliminated_polynomials: list[sp.Expr],
    pivot: sp.Expr,
    h: sp.Expr,
    target: sp.Expr,
    tracked_basis: bool = False,
) -> dict:
    """Compute and parse one deterministic Singular lift over F_p."""
    program = f"""
option(redSB);
ring r={prime},(c,d,lam),dp;
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        program += f"poly f{order}={singular_expression(polynomial)};\n"
    program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
ideal F=f3,f4,f5,f6;
ideal T=target;
"""
    if tracked_basis:
        program += """
matrix transformation;
ideal G=liftstd(F,transformation,"slimgb");
matrix targetCoordinates=lift(G,T);
matrix L=transformation*targetCoordinates;
"""
    else:
        program += """
matrix L=lift(F,T);
"""
    program += """
matrix reconstructed=matrix(F)*L;
int row;
if (reconstructed[1,1]-target!=0)
{
  print("BAD_FINITE_FIELD_TARGET_LIFT");
  exit(1);
}
for (row=1;row<=nrows(L);row++)
{
  print("ROW_BEGIN");
  print(L[row,1]);
  print("ROW_END");
}
"""
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    errors = [
        line
        for line in (completed.stdout + completed.stderr).splitlines()
        if line.lstrip().startswith("?")
    ]
    if errors:
        raise RuntimeError("\n".join(errors[:20]))
    blocks = completed.stdout.split("ROW_BEGIN\n")[1:]
    if len(blocks) != 4:
        raise RuntimeError(completed.stdout[:1000])
    rows = [
        singular_polynomial(
            block.split("\nROW_END", 1)[0].strip(),
            prime,
        )
        for block in blocks
    ]
    support = [
        [row_index, *monomial]
        for row_index, row in enumerate(rows)
        for monomial in sorted(row)
    ]
    coefficients = [
        rows[row_index][tuple(item[1:])]
        for item in support
        for row_index in [item[0]]
    ]
    support_bytes = json.dumps(
        support,
        separators=(",", ":"),
    ).encode()
    return {
        "prime": prime,
        "seconds": round(time.monotonic() - started, 3),
        "support": support,
        "support_sha256": hashlib.sha256(support_bytes).hexdigest(),
        "term_counts": [len(row) for row in rows],
        "coefficients": coefficients,
    }


def canonical_syzygy_record(
    singular: str,
    prime: int,
    exponent: int,
    eliminated_polynomials: list[sp.Expr],
    pivot: sp.Expr,
    h: sp.Expr,
    target: sp.Expr,
    algorithm: str,
    component_order: bool,
) -> dict:
    """Reduce a lift by a standard basis of the generator syzygies."""
    ordering = "(dp,C)" if component_order else "dp"
    program = f"""
option(redSB);
option(returnSB);
ring r={prime},(c,d,lam),{ordering};
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        program += f"poly f{order}={singular_expression(polynomial)};\n"
    program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
ideal F=f3,f4,f5,f6;
ideal T=target;
matrix L=lift(F,T);
module S=syz(F,"{algorithm}");
module M=[L[1,1],L[2,1],L[3,1],L[4,1]];
module N=reduce(M,S);
matrix C=N;
matrix reconstructed=matrix(F)*C;
if (reconstructed[1,1]-target!=0)
{{
  print("BAD_CANONICAL_SYZYGY_LIFT");
  exit(1);
}}
int row;
for (row=1;row<=nrows(C);row++)
{{
  print("ROW_BEGIN");
  print(C[row,1]);
  print("ROW_END");
}}
"""
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    errors = [
        line
        for line in (completed.stdout + completed.stderr).splitlines()
        if line.lstrip().startswith("?")
    ]
    if errors:
        raise RuntimeError("\n".join(errors[:20]))
    blocks = completed.stdout.split("ROW_BEGIN\n")[1:]
    if len(blocks) != 4:
        raise RuntimeError(completed.stdout[:1000])
    rows = [
        singular_polynomial(
            block.split("\nROW_END", 1)[0].strip(),
            prime,
        )
        for block in blocks
    ]
    support = [
        [row_index, *monomial]
        for row_index, row in enumerate(rows)
        for monomial in sorted(row)
    ]
    support_bytes = json.dumps(
        support,
        separators=(",", ":"),
    ).encode()
    return {
        "prime": prime,
        "seconds": round(time.monotonic() - started, 3),
        "support": support,
        "support_sha256": hashlib.sha256(support_bytes).hexdigest(),
        "term_counts": [len(row) for row in rows],
    }


def merge_crt(
    residues: list[int],
    modulus: int,
    coefficients: list[int],
    prime: int,
) -> int:
    """Merge one coefficient vector into a componentwise CRT state."""
    inverse = pow(modulus, -1, prime)
    for index, coefficient in enumerate(coefficients):
        correction = ((coefficient - residues[index]) * inverse) % prime
        residues[index] += modulus * correction
    return modulus * prime


def crt_diagnostic(state: dict) -> tuple[dict, list[tuple[int, int] | None]]:
    """Reconstruct the CRT state and test candidates at the holdout prime."""
    reconstructed = [
        rational_reconstruct(value, state["modulus"])
        for value in state["residues"]
    ]
    holdout = state["holdout"]
    holdout_prime = holdout["prime"]
    validated = 0
    maximum_numerator_bits = 0
    maximum_denominator_bits = 0
    for candidate, expected in zip(
        reconstructed,
        holdout["coefficients"],
        strict=True,
    ):
        if candidate is None:
            continue
        numerator, denominator = candidate
        if denominator % holdout_prime == 0:
            continue
        observed = (
            numerator * pow(denominator, -1, holdout_prime)
        ) % holdout_prime
        if observed == expected:
            validated += 1
            maximum_numerator_bits = max(
                maximum_numerator_bits,
                abs(numerator).bit_length(),
            )
            maximum_denominator_bits = max(
                maximum_denominator_bits,
                denominator.bit_length(),
            )
    return (
        {
            "build_primes": len(state["crt_primes"]),
            "build_modulus_bits": state["modulus"].bit_length(),
            "holdout_prime": holdout_prime,
            "candidates": sum(value is not None for value in reconstructed),
            "validated": validated,
            "total": len(reconstructed),
            "maximum_numerator_bits": maximum_numerator_bits,
            "maximum_denominator_bits": maximum_denominator_bits,
        },
        reconstructed,
    )


def singular_rational_polynomial(
    terms: list[tuple[list[int], tuple[int, int]]],
) -> str:
    """Format sparse rational terms for Singular."""
    rendered = []
    for support, coefficient in terms:
        numerator, denominator = coefficient
        coefficient_text = (
            str(numerator)
            if denominator == 1
            else f"({numerator}/{denominator})"
        )
        factors = [coefficient_text]
        for variable, exponent in zip(
            ("c", "d", "lam"),
            support,
            strict=True,
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        rendered.append("*".join(factors))
    return "+".join(rendered).replace("+-", "-") or "0"


def verify_rational_lift(
    singular: str,
    support: list[list[int]],
    reconstructed: list[tuple[int, int]],
    exponent: int,
    eliminated_polynomials: list[sp.Expr],
    pivot: sp.Expr,
    h: sp.Expr,
    target: sp.Expr,
) -> float:
    """Verify the reconstructed sparse certificate exactly over Q."""
    multipliers: list[list[tuple[list[int], tuple[int, int]]]] = [
        [] for _ in eliminated_polynomials
    ]
    for item, coefficient in zip(support, reconstructed, strict=True):
        multipliers[item[0]].append((item[1:], coefficient))
    program = """
ring r=0,(c,d,lam),dp;
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        program += f"poly f{order}={singular_expression(polynomial)};\n"
    for order, terms in enumerate(multipliers, start=3):
        program += (
            f"poly A{order}={singular_rational_polynomial(terms)};\n"
        )
    program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
if (A3*f3+A4*f4+A5*f5+A6*f6-target!=0)
{{
  print("BAD_RECONSTRUCTED_RATIONAL_LIFT");
  exit(1);
}}
print("EXACT_RECONSTRUCTED_RATIONAL_LIFT");
"""
    started = time.monotonic()
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=1200,
    )
    if "EXACT_RECONSTRUCTED_RATIONAL_LIFT" not in completed.stdout:
        raise RuntimeError(
            (completed.stdout + completed.stderr)[:2000]
        )
    return time.monotonic() - started


def apolar_moment(order: int) -> sp.Expr:
    symbol_power = sp.Poly(sp.expand(A**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    return sp.expand(
        sum(
            coefficient
            * polynomial_power.coeff_monomial(x**x_order * y**y_order)
            * factorial(x_order)
            * factorial(y_order)
            for (x_order, y_order), coefficient in symbol_power.terms()
        )
    )


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def main() -> None:
    moments = [apolar_moment(order) for order in range(1, 7)]
    b_solution = sp.Rational(2, 3) * c * (lam + 1) - d * lam
    chart_moments = [
        sp.cancel(moment.subs(b, b_solution).subs(e, 1))
        for moment in moments[1:]
    ]
    second = sp.Poly(chart_moments[0], a)
    pivot = sp.factor(second.coeff_monomial(a) / 576)
    constant = second.coeff_monomial(1)
    a_solution = sp.cancel(-constant / second.coeff_monomial(a))
    eliminated_polynomials = [
        sp.Poly(
            sp.fraction(sp.cancel(moment.subs(a, a_solution)))[0],
            c,
            d,
            domain=sp.QQ[lam],
        ).primitive()[1].as_expr()
        for moment in chart_moments[1:]
    ]
    h = 8 * c - 3 * d**2
    target = lam**4 * (lam - 1) ** 4
    extra_quadratic = (
        (25 * lam**4 - 50 * lam**3 + 75 * lam**2 - 50 * lam + 25)
        * d**2
        + (
            100 * lam**5
            - 350 * lam**4
            + 100 * lam**3
            + 100 * lam**2
            - 350 * lam
            + 100
        )
        * d
        - 44 * lam**6
        - 68 * lam**5
        + 386 * lam**4
        + 208 * lam**3
        + 386 * lam**2
        - 68 * lam
        - 44
    )

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    if len(sys.argv) >= 6 and sys.argv[1] == "--compare-syzygy-lifts":
        exponent = int(sys.argv[2])
        algorithm = sys.argv[3]
        component_order = sys.argv[4] == "component"
        if algorithm not in {"std", "slimgb"}:
            raise ValueError("algorithm must be std or slimgb")
        if sys.argv[4] not in {"dp", "component"}:
            raise ValueError("order must be dp or component")
        primes = [int(value) for value in sys.argv[5:]]
        first_support: list[list[int]] | None = None
        for prime in primes:
            record = canonical_syzygy_record(
                singular,
                prime,
                exponent,
                eliminated_polynomials,
                pivot,
                h,
                target,
                algorithm,
                component_order,
            )
            if first_support is None:
                first_support = record["support"]
            print(
                "CANONICAL_SYZYGY_LIFT",
                algorithm,
                sys.argv[4],
                prime,
                record["seconds"],
                record["term_counts"],
                record["support_sha256"],
                record["support"] == first_support,
                flush=True,
            )
        return

    if len(sys.argv) >= 4 and sys.argv[1] == "--compare-tracked-lifts":
        exponent = int(sys.argv[2])
        primes = [int(value) for value in sys.argv[3:]]
        first_support: list[list[int]] | None = None
        for prime in primes:
            record = raw_lift_record(
                singular,
                prime,
                exponent,
                eliminated_polynomials,
                pivot,
                h,
                target,
                tracked_basis=True,
            )
            if first_support is None:
                first_support = record["support"]
            print(
                "TRACKED_LIFT",
                prime,
                record["seconds"],
                record["term_counts"],
                record["support_sha256"],
                record["support"] == first_support,
                flush=True,
            )
        return

    if len(sys.argv) >= 5 and sys.argv[1] == "--crt-lift":
        exponent = int(sys.argv[2])
        checkpoint_path = Path(sys.argv[3])
        primes = [int(value) for value in sys.argv[4:]]
        if checkpoint_path.exists():
            payload = json.loads(checkpoint_path.read_text())
            if payload["format"] != "rank-two-squarefree-crt-v1":
                raise ValueError("unrecognized CRT checkpoint format")
            if payload["exponent"] != exponent:
                raise ValueError("checkpoint exponent mismatch")
            state = {
                "support_sha256": payload["support_sha256"],
                "support": payload["support"],
                "modulus": int(payload["crt"]["modulus"], 16),
                "residues": [
                    int(value, 16)
                    for value in payload["crt"]["residues"]
                ],
                "crt_primes": payload["crt"]["primes"],
                "holdout": payload["holdout"],
                "records": payload["records"],
            }
        else:
            state = {
                "support_sha256": None,
                "support": [],
                "modulus": 1,
                "residues": [],
                "crt_primes": [],
                "holdout": None,
                "records": [],
            }

        def save_checkpoint(extra: dict | None = None) -> None:
            checkpoint = {
                "format": "rank-two-squarefree-crt-v1",
                "exponent": exponent,
                "identity": "lambda^4*(lambda-1)^4*(p*h)^5=sum(A_i*f_i)",
                "support_sha256": state["support_sha256"],
                "support": state["support"],
                "crt": {
                    "modulus": hex(state["modulus"]),
                    "primes": state["crt_primes"],
                    "residues": [
                        hex(value) for value in state["residues"]
                    ],
                },
                "holdout": state["holdout"],
                "records": state["records"],
            }
            if extra is not None:
                checkpoint.update(extra)
            temporary_path = checkpoint_path.with_suffix(
                checkpoint_path.suffix + ".tmp"
            )
            temporary_path.write_text(
                json.dumps(checkpoint, indent=2) + "\n"
            )
            temporary_path.replace(checkpoint_path)

        seen_primes = {record["prime"] for record in state["records"]}
        for prime in primes:
            if prime in seen_primes:
                print("SKIP_EXISTING_PRIME", prime, flush=True)
                continue
            record = raw_lift_record(
                singular,
                prime,
                exponent,
                eliminated_polynomials,
                pivot,
                h,
                target,
            )
            status = "good"
            if state["support_sha256"] is None:
                state["support_sha256"] = record["support_sha256"]
                state["support"] = record["support"]
                state["residues"] = [0] * len(record["support"])
            elif (
                record["support_sha256"] != state["support_sha256"]
                or record["support"] != state["support"]
            ):
                status = "unlucky"
            if status == "good":
                if state["holdout"] is not None:
                    previous = state["holdout"]
                    state["modulus"] = merge_crt(
                        state["residues"],
                        state["modulus"],
                        previous["coefficients"],
                        previous["prime"],
                    )
                    state["crt_primes"].append(previous["prime"])
                state["holdout"] = {
                    "prime": prime,
                    "coefficients": record["coefficients"],
                }
            state["records"].append(
                {
                    "prime": prime,
                    "status": status,
                    "seconds": record["seconds"],
                    "term_counts": record["term_counts"],
                    "support_sha256": record["support_sha256"],
                }
            )
            seen_primes.add(prime)
            save_checkpoint()
            print(
                "CRT_IMAGE",
                prime,
                status,
                record["seconds"],
                record["term_counts"],
                flush=True,
            )
            if state["crt_primes"]:
                diagnostic, _ = crt_diagnostic(state)
                print(
                    "CRT_DIAGNOSTIC",
                    json.dumps(diagnostic, sort_keys=True),
                    flush=True,
                )

        if not state["crt_primes"] or state["holdout"] is None:
            print("NEED_AT_LEAST_TWO_GOOD_PRIMES")
            return
        diagnostic, reconstructed_optional = crt_diagnostic(state)
        print(
            "CRT_FINAL_DIAGNOSTIC",
            json.dumps(diagnostic, sort_keys=True),
        )
        if diagnostic["validated"] != diagnostic["total"]:
            print("MORE_PRIMES_REQUIRED")
            return
        reconstructed = [
            value
            for value in reconstructed_optional
            if value is not None
        ]
        exact_seconds = verify_rational_lift(
            singular,
            state["support"],
            reconstructed,
            exponent,
            eliminated_polynomials,
            pivot,
            h,
            target,
        )
        exact_certificate = {
            "exact_verification": {
                "result": "pass",
                "seconds": round(exact_seconds, 3),
                "coefficient_count": len(reconstructed),
                "maximum_numerator_bits": max(
                    abs(value[0]).bit_length()
                    for value in reconstructed
                ),
                "maximum_denominator_bits": max(
                    value[1].bit_length()
                    for value in reconstructed
                ),
            },
            "coefficients": reconstructed,
        }
        save_checkpoint(exact_certificate)
        print(
            "EXACT_RECONSTRUCTED_RATIONAL_LIFT",
            json.dumps(
                exact_certificate["exact_verification"],
                sort_keys=True,
            ),
        )
        return

    if len(sys.argv) == 2 and sys.argv[1] == "--generic-q-lift":
        program = """
ring r=(0,lam),(c,d),dp;
option(redSB);
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly q={singular_expression(extra_quadratic)};
ideal J=q,f3,f4,f5,f6;
ideal T=p^3;
matrix L=lift(J,T);
poly reconstructed=0;
int row;
for (row=1;row<=size(J);row++)
{{
  reconstructed=reconstructed+J[row]*L[row,1];
}}
if (reconstructed-p^3!=0)
{{
  print("BAD_GENERIC_Q_LIFT");
  exit(1);
}}
for (row=1;row<=size(J);row++)
{{
  poly work=L[row,1];
  while (work!=0)
  {{
    print("DENOMINATOR");
    print(denominator(leadcoef(work)));
    work=work-lead(work);
  }}
}}
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        denominator_texts = {
            line.strip()
            for marker, line in zip(
                completed.stdout.splitlines(),
                completed.stdout.splitlines()[1:],
            )
            if marker.strip() == "DENOMINATOR"
        }
        factors: set[str] = set()
        for denominator_text in denominator_texts:
            denominator_expression = sp.sympify(
                denominator_text.replace("^", "**"),
                locals={"lam": lam},
            )
            _, factor_list = sp.factor_list(denominator_expression, lam)
            factors.update(str(factor) for factor, _ in factor_list)
        print("GENERIC_Q_LIFT_UNIQUE_DENOMINATORS", len(denominator_texts))
        print("GENERIC_Q_LIFT_DENOMINATOR_FACTORS", sorted(factors))
        return

    if len(sys.argv) >= 4 and sys.argv[1] in {
        "--compare-lifts",
        "--compare-raw-lifts",
    }:
        exponent = int(sys.argv[2])
        primes = [int(value) for value in sys.argv[3:]]
        canonical = sys.argv[1] == "--compare-lifts"
        supports: list[list[set[tuple[int, int, int]]]] = []
        for prime in primes:
            program = f"""
option(redSB);
ring r={prime},(c,d,lam),dp;
"""
            for order, polynomial in enumerate(eliminated_polynomials, start=3):
                program += (
                    f"poly f{order}={singular_expression(polynomial)};\n"
                )
            program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
ideal F=f3,f4,f5,f6;
ideal T=target;
matrix L=lift(F,T);
"""
            if canonical:
                program += """
module S=std(syz(F));
module M=[L[1,1],L[2,1],L[3,1],L[4,1]];
module N=reduce(M,S);
matrix C=N;
"""
                lift_matrix = "C"
            else:
                lift_matrix = "L"
            program += f"""
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{{
  reconstructed=reconstructed+F[row]*{lift_matrix}[row,1];
}}
if (reconstructed-target!=0)
{{
  print("BAD_CANONICAL_LIFT");
  exit(1);
}}
for (row=1;row<=size(F);row++)
{{
  print("ROW_BEGIN");
  print({lift_matrix}[row,1]);
  print("ROW_END");
}}
"""
            completed = subprocess.run(
                [singular, "-q"],
                input=program,
                text=True,
                capture_output=True,
                check=True,
                timeout=300,
            )
            blocks = completed.stdout.split("ROW_BEGIN\n")[1:]
            assert len(blocks) == 4, completed.stdout[:1000]
            prime_supports: list[set[tuple[int, int, int]]] = []
            for block in blocks:
                expression_text = block.split("\nROW_END", 1)[0].strip()
                prime_supports.append(singular_support(expression_text))
            supports.append(prime_supports)
            print(
                "prime",
                prime,
                "term_counts",
                [len(value) for value in prime_supports],
            )
        for left in range(1, len(supports)):
            print(
                "support_matches_first",
                primes[left],
                [
                    supports[0][row] == supports[left][row]
                    for row in range(4)
                ],
            )
        return

    if len(sys.argv) == 4 and sys.argv[1] == "--prime-lift":
        prime = int(sys.argv[2])
        exponent = int(sys.argv[3])
        assert prime > 2 and exponent >= 0
        program = f"""
option(redSB);
ring r={prime},(c,d,lam),dp;
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
ideal F=f3,f4,f5,f6;
ideal T=target;
matrix L=lift(F,T);
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{{
  reconstructed=reconstructed+F[row]*L[row,1];
}}
if (reconstructed-target!=0)
{{
  print("BAD_FINITE_FIELD_TARGET_LIFT");
  exit(1);
}}
print("FINITE_FIELD_TARGET_LIFT");
for (row=1;row<=size(F);row++)
{{
  print(deg(L[row,1]));
  print(size(L[row,1]));
}}
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        return

    if len(sys.argv) == 3 and sys.argv[1] == "--direct":
        exponent = int(sys.argv[2])
        assert exponent >= 0
        program = """
option(redSB);
ring r=0,(c,d,lam),dp;
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)}*(p*h)^{exponent};
ideal F=f3,f4,f5,f6;
ideal T=target;
timer=1;
matrix L=lift(F,T);
int lift_seconds=timer;
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{{
  reconstructed=reconstructed+F[row]*L[row,1];
}}
if (reconstructed-target!=0)
{{
  print("BAD_DIRECT_TARGET_LIFT");
  exit(1);
}}
print("EXACT_DIRECT_TARGET_CERTIFICATE");
print("exponent");
print({exponent});
print("lift_seconds");
print(lift_seconds);
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=1200,
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        return

    if len(sys.argv) == 3 and sys.argv[1] == "--prime":
        prime = int(sys.argv[2])
        assert prime > 2
        program = f"""
option(redSB);
ring r={prime},(c,d,lam),dp;
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)};
poly multiplier=p*h;
ideal F=f3,f4,f5,f6;
ideal G=std(F);
int exponent=-1;
poly powered=target;
int candidate;
for (candidate=0;candidate<=32;candidate++)
{{
  if (reduce(powered,G)==0)
  {{
    exponent=candidate;
    break;
  }}
  powered=powered*multiplier;
}}
print("FINITE_FIELD_SATURATION_EXPONENT");
print(exponent);
print("FINITE_FIELD_BASIS_SIZE");
print(size(G));
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        return

    if len(sys.argv) == 2 and sys.argv[1] == "--saturation":
        program = """
LIB "modstd.lib";
option(redSB);
ring r=0,(c,d,lam),dp;
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            program += f"poly f{order}={singular_expression(polynomial)};\n"
        program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)};
poly multiplier=p*h;
ideal F=f3,f4,f5,f6;
timer=1;
ideal G=modStd(F);
int modular_seconds=timer;
ideal V=std(G);
int exponent=-1;
poly powered=target;
int candidate;
for (candidate=0;candidate<=32;candidate++)
{{
  if (reduce(powered,V)==0)
  {{
    exponent=candidate;
    break;
  }}
  powered=powered*multiplier;
}}
if (exponent<0)
{{
  print("NO_SATURATION_EXPONENT_THROUGH_32");
  exit(1);
}}
ideal T=powered;
timer=1;
matrix L=lift(F,T);
int lift_seconds=timer;
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{{
  reconstructed=reconstructed+F[row]*L[row,1];
}}
if (reconstructed-powered!=0)
{{
  print("BAD_TARGET_ONLY_LIFT");
  exit(1);
}}
print("EXACT_TARGET_SATURATION_CERTIFICATE");
print("exponent");
print(exponent);
print("modular_seconds");
print(modular_seconds);
print("lift_seconds");
print(lift_seconds);
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=1200,
        )
        print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())
        return

    program = """
LIB "modstd.lib";
option(redSB);
ring r=0,(z,c,d,lam),lp;
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        program += f"poly f{order}={singular_expression(polynomial)};\n"
    program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(h)};
poly target={singular_expression(target)};
ideal I=f3,f4,f5,f6,z*p*h-1;
timer=1;
ideal G=modStd(I);
int modular_seconds=timer;
ideal V=std(G);
if (size(reduce(I,V))!=0)
{{
  print("BAD_FORWARD_CONTAINMENT");
  exit(1);
}}
if (reduce(target,V)!=0)
{{
  print("BAD_TARGET_REDUCTION");
  exit(1);
}}
timer=1;
matrix T=lift(I,G);
int lift_seconds=timer;
int row;
int column;
poly reconstructed;
for (column=1;column<=size(G);column++)
{{
  reconstructed=0;
  for (row=1;row<=size(I);row++)
  {{
    reconstructed=reconstructed+I[row]*T[row,column];
  }}
  if (reconstructed-G[column]!=0)
  {{
    print("BAD_REVERSE_CONTAINMENT");
    exit(1);
  }}
}}
print("EXACT_BIDIRECTIONAL_CERTIFICATE");
print("basis_size");
print(size(G));
print("modular_seconds");
print(modular_seconds);
print("lift_seconds");
print(lift_seconds);
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=1200,
    )
    print(completed.stdout.strip())
    if completed.stderr.strip():
        print(completed.stderr.strip())


if __name__ == "__main__":
    main()

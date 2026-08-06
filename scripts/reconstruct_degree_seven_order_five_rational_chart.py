#!/usr/bin/env python3
"""Reconstruct and certify the degree-seven order-five chart over Q.

The input files are independently interpolated saturated Groebner bases over
finite fields.  Coefficientwise CRT and balanced rational reconstruction give
a candidate basis over Q.  We require stability after dropping the final
build prime, exact reduction to every build image, and (when supplied) exact
agreement with a prime that was not used in reconstruction.  Singular then
certifies the characteristic-zero dimension, length, and lexicographic shape.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt, prod
import json
from pathlib import Path
import shutil
import subprocess

import sympy as sp


DEFAULT_PRIMES = (
    1013,
    1019,
    1021,
    1031,
    1033,
    1039,
    1049,
    1051,
    1061,
    1063,
    1069,
    1087,
    1091,
    1093,
    1097,
)
PROFILE_PRIMES = (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)


def load_basis(path: Path) -> list[dict[tuple[int, int], int]]:
    data = json.loads(path.read_text())
    return [
        {
            (term["sigma_degree"], term["tau_degree"]): term["coefficient"]
            for term in polynomial
        }
        for polynomial in data["zero_scheme"]["groebner_basis_terms"]
    ]


def crt_merge(residue: int, modulus: int, image: int, prime: int):
    merged = residue + modulus * (
        ((image - residue) % prime) * pow(modulus, -1, prime) % prime
    )
    return merged, modulus * prime


def reconstruct(residue: int, modulus: int) -> Fraction | None:
    """Balanced rational reconstruction with the standard sqrt(M/2) bound."""

    residue %= modulus
    bound = isqrt(modulus // 2)
    old_remainder, remainder = modulus, residue
    old_coefficient, coefficient = 0, 1
    while abs(remainder) >= bound:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_coefficient, coefficient = (
            coefficient,
            old_coefficient - quotient * coefficient,
        )
    if coefficient < 0:
        remainder, coefficient = -remainder, -coefficient
    if (
        coefficient == 0
        or abs(remainder) >= bound
        or coefficient >= bound
        or gcd(remainder, coefficient) != 1
        or (remainder - residue * coefficient) % modulus
    ):
        return None
    return Fraction(remainder, coefficient)


def reconstruct_basis(
    primes: list[int],
    images: list[list[dict[tuple[int, int], int]]],
) -> list[dict[tuple[int, int], Fraction]]:
    result = []
    for polynomial_index in range(len(images[0])):
        support = sorted(
            set().union(
                *(set(image[polynomial_index]) for image in images)
            )
        )
        polynomial = {}
        for monomial in support:
            residue = images[0][polynomial_index].get(monomial, 0)
            modulus = primes[0]
            for prime, image in zip(primes[1:], images[1:], strict=True):
                residue, modulus = crt_merge(
                    residue,
                    modulus,
                    image[polynomial_index].get(monomial, 0),
                    prime,
                )
            coefficient = reconstruct(residue, modulus)
            if coefficient is None:
                raise AssertionError(
                    f"failed rational reconstruction for G{polynomial_index + 1}"
                    f" monomial {monomial} modulo {modulus}"
                )
            if coefficient:
                polynomial[monomial] = coefficient
        result.append(polynomial)
    return result


def reduction(
    basis: list[dict[tuple[int, int], Fraction]], prime: int
) -> list[dict[tuple[int, int], int]]:
    return [
        {
            monomial: (
                coefficient.numerator
                * pow(coefficient.denominator, -1, prime)
            )
            % prime
            for monomial, coefficient in polynomial.items()
            if coefficient.numerator % prime
        }
        for polynomial in basis
    ]


def rational_points(basis, prime: int):
    modular = reduction(basis, prime)
    return [
        (sigma, tau)
        for sigma in range(prime)
        for tau in range(prime)
        if all(
            sum(
                coefficient
                * pow(sigma, sigma_degree, prime)
                * pow(tau, tau_degree, prime)
                for (sigma_degree, tau_degree), coefficient in polynomial.items()
            )
            % prime
            == 0
            for polynomial in modular
        )
    ]


def render_polynomial(polynomial) -> str:
    terms = []
    for (sigma_degree, tau_degree), coefficient in sorted(polynomial.items()):
        factors = [f"({coefficient.numerator}/{coefficient.denominator})"]
        if sigma_degree:
            factors.append("sigma" if sigma_degree == 1 else f"sigma^{sigma_degree}")
        if tau_degree:
            factors.append("tau" if tau_degree == 1 else f"tau^{tau_degree}")
        terms.append("*".join(factors))
    return "+".join(terms).replace("+-", "-")


def serialize(polynomial):
    return [
        {
            "sigma_degree": sigma_degree,
            "tau_degree": tau_degree,
            "numerator": coefficient.numerator,
            "denominator": coefficient.denominator,
        }
        for (sigma_degree, tau_degree), coefficient in sorted(polynomial.items())
    ]


def parse_polynomial(text: str, variables: tuple[sp.Symbol, sp.Symbol]):
    expression = sp.sympify(
        text.replace("^", "**"),
        locals={str(variable): variable for variable in variables},
    )
    return sp.Poly(expression, *variables, domain=sp.QQ)


def singular_exact(basis):
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for exact certification")
    source = f"""ring r=0,(sigma,tau),dp;
option(redSB);
ideal I={','.join(render_polynomial(polynomial) for polynomial in basis)};
ideal G=std(I);
print(\"DP_DIMENSION=\"+string(dim(G)));
print(\"DP_VDIM=\"+string(vdim(G)));
ring l=0,(tau,sigma),lp;
option(redSB);
ideal J=imap(r,G);
ideal L=std(J);
print(\"LEX_SIZE=\"+string(size(L)));
print(\"LEX_DIMENSION=\"+string(dim(L)));
print(\"LEX_VDIM=\"+string(vdim(L)));
for (int i=1; i<=size(L); i++) {{ print(\"LEX_\"+string(i)+\"=\"+string(L[i])); }}
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
    )
    if completed.stderr.strip() or "   ?" in completed.stdout:
        raise AssertionError(completed.stdout + completed.stderr)
    values = {}
    lex = []
    for line in completed.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.startswith("LEX_") and key not in {
            "LEX_SIZE",
            "LEX_DIMENSION",
            "LEX_VDIM",
        }:
            lex.append(value)
        else:
            values[key] = int(value)
    assert values == {
        "DP_DIMENSION": 0,
        "DP_VDIM": 8,
        "LEX_SIZE": 2,
        "LEX_DIMENSION": 0,
        "LEX_VDIM": 8,
    }
    tau, sigma = sp.symbols("tau sigma")
    lex_polynomials = [parse_polynomial(item, (tau, sigma)) for item in lex]
    univariate = next(
        polynomial
        for polynomial in lex_polynomials
        if polynomial.degree(tau) == 0 and polynomial.degree(sigma) == 8
    )
    linear_tau = next(
        polynomial
        for polynomial in lex_polynomials
        if polynomial.degree(tau) == 1
    )
    factorization = sp.factor_list(univariate.as_expr())
    assert len(factorization[1]) == 1
    assert factorization[1][0][1] == 1
    assert sp.Poly(factorization[1][0][0], sigma).degree() == 8
    return values, lex_polynomials, univariate, linear_tau


def serialize_sympy(polynomial: sp.Poly):
    return [
        {
            "tau_degree": monomial[0],
            "sigma_degree": monomial[1],
            "coefficient": str(coefficient),
        }
        for monomial, coefficient in polynomial.terms()
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primes",
        type=int,
        nargs="+",
        default=list(DEFAULT_PRIMES),
    )
    parser.add_argument("--holdout-prime", type=int)
    parser.add_argument(
        "--artifact-directory",
        type=Path,
        default=Path("artifacts/generated-results"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if len(args.primes) < 3:
        raise SystemExit("at least three build primes are required")

    def path_for(prime: int):
        return args.artifact_directory / (
            f"degree_seven_order_five_chart_gf{prime}.json"
        )

    images = [load_basis(path_for(prime)) for prime in args.primes]
    basis = reconstruct_basis(args.primes, images)
    previous = reconstruct_basis(args.primes[:-1], images[:-1])
    assert previous == basis, "rational reconstruction is not yet stable"
    for prime, image in zip(args.primes, images, strict=True):
        assert reduction(basis, prime) == image

    holdout = None
    if args.holdout_prime is not None:
        holdout_image = load_basis(path_for(args.holdout_prime))
        assert reduction(basis, args.holdout_prime) == holdout_image
        holdout = {
            "prime": args.holdout_prime,
            "status": "exact agreement; prime excluded from CRT",
        }

    exact, lex, univariate, linear_tau = singular_exact(basis)
    prime_profile = {}
    for prime in PROFILE_PRIMES:
        scan_path = args.artifact_directory / (
            f"degree_seven_order_five_scan_gf{prime}.json"
        )
        scan = json.loads(scan_path.read_text())
        expected = [
            (point["sigma"], point["tau"])
            for point in scan["order_five_consistent_points"]
        ]
        points = rational_points(basis, prime)
        assert points == expected
        prime_profile[str(prime)] = [list(point) for point in points]
    certificate = {
        "scope": "exact characteristic-zero order-five pivot-chart certificate",
        "build_primes": args.primes,
        "crt_modulus": str(prod(args.primes)),
        "crt_modulus_bits": prod(args.primes).bit_length(),
        "stability": "identical after dropping final build prime",
        "holdout": holdout,
        "full_plane_prime_profile": prime_profile,
        "full_plane_points_matched": sum(map(len, prime_profile.values())),
        "degree_reverse_lexicographic_basis": [serialize(item) for item in basis],
        "exact_zero_scheme": {
            "dimension": exact["DP_DIMENSION"],
            "vector_space_dimension": exact["DP_VDIM"],
            "lexicographic_basis": [serialize_sympy(item) for item in lex],
            "primitive_sigma_polynomial": str(
                sp.primitive(univariate.as_expr(), expand=True)[1]
            ),
            "linear_tau_relation": str(linear_tau.as_expr()),
            "primitive_sigma_polynomial_irreducible_over_Q": True,
            "component": "one reduced closed point of degree eight over Q",
        },
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        "PASS: stable rational reconstruction from {} primes ({}-bit CRT)".format(
            len(args.primes), certificate["crt_modulus_bits"]
        )
    )
    if holdout:
        print(f"PASS: exact agreement at independent GF({args.holdout_prime})")
    print("PASS: exact Q zero scheme has dimension 0 and length 8")
    print("PASS: primitive degree-eight sigma polynomial is irreducible over Q")


if __name__ == "__main__":
    main()

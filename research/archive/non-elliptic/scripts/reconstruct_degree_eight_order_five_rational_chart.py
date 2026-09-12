#!/usr/bin/env python3
"""Reconstruct the degree-eight order-five chart over Q from modular bases."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import json
from math import gcd, prod
from pathlib import Path
import shutil
import subprocess

import sympy as sp

from reconstruct_degree_seven_order_five_rational_chart import (
    crt_merge,
    load_basis,
    reduction,
    render_polynomial,
    serialize,
)


DEFAULT_PRIMES = (
    2147483647,
    2099999999,
    2049999979,
    1019,
    1021,
    1013,
)


def simultaneous_reconstruct_basis(primes, images):
    """LLL-reconstruct each normalized basis vector with one denominator.

    Coefficientwise rational reconstruction needs a modulus roughly quadratic
    in the coefficient height.  A normalized Groebner polynomial instead has
    one common denominator, so its entire coefficient vector can be recovered
    as one short vector in a congruence lattice.  Exact reduction to every
    modular image and stability under dropping a prime remain mandatory.
    """

    result = []
    for polynomial_index in range(len(images[0])):
        support = sorted(
            set().union(*(set(image[polynomial_index]) for image in images))
        )
        leading_candidates = [
            monomial
            for monomial in support
            if all(
                image[polynomial_index].get(monomial) == 1
                for image in images
            )
        ]
        assert len(leading_candidates) == 1
        leading = leading_candidates[0]
        monomials = [monomial for monomial in support if monomial != leading]
        residues = []
        modulus = prod(primes)
        for monomial in monomials:
            residue = images[0][polynomial_index].get(monomial, 0)
            current_modulus = primes[0]
            for prime, image in zip(primes[1:], images[1:], strict=True):
                residue, current_modulus = crt_merge(
                    residue,
                    current_modulus,
                    image[polynomial_index].get(monomial, 0),
                    prime,
                )
            assert current_modulus == modulus
            residues.append(residue)
        lattice = [[1, *residues]] + [
            [0]
            + [modulus if index == diagonal else 0 for index in range(len(monomials))]
            for diagonal in range(len(monomials))
        ]
        reduced_lattice = sp.Matrix(lattice).lll().tolist()
        candidates = [row for row in reduced_lattice if row[0]]
        assert candidates
        row = min(candidates, key=lambda item: sum(value * value for value in item))
        common = reduce(gcd, (abs(value) for value in row if value), 0)
        row = [value // common for value in row]
        if row[0] < 0:
            row = [-value for value in row]
        denominator = row[0]
        polynomial = {leading: Fraction(1)}
        polynomial.update(
            {
                monomial: Fraction(numerator, denominator)
                for monomial, numerator in zip(monomials, row[1:], strict=True)
                if numerator
            }
        )
        result.append(polynomial)
    return result


def complete_exact_basis(seed_basis, leading_monomials):
    """Complete stable characteristic-zero generators to the reduced basis."""

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for exact reconstruction")
    source = f"""ring r=0,(sigma,tau),dp;
option(redSB);
ideal I={','.join(render_polynomial(polynomial) for polynomial in seed_basis)};
ideal G=std(I);
print("SIZE="+string(size(G)));
for (int i=1; i<=size(G); i++) {{ print("G_"+string(i)+"="+string(G[i])); }}
"""
    completed = subprocess.run(
        [singular, "-q"], input=source, text=True, capture_output=True, check=True
    )
    if completed.stderr.strip() or "   ?" in completed.stdout:
        raise AssertionError(completed.stdout + completed.stderr)
    sigma, tau = sp.symbols("sigma tau")
    raw_basis = []
    size = None
    for line in completed.stdout.splitlines():
        if line.startswith("SIZE="):
            size = int(line.split("=", 1)[1])
        elif line.startswith("G_"):
            polynomial = parse_polynomial(line.split("=", 1)[1], (sigma, tau))
            raw_basis.append(
                {
                    monomial: Fraction(coefficient.p, coefficient.q)
                    for monomial, coefficient in polynomial.terms()
                }
            )
    assert size == len(raw_basis) == len(leading_monomials) == 5
    basis = []
    for polynomial, leading in zip(raw_basis, leading_monomials, strict=True):
        leading_coefficient = polynomial[leading]
        basis.append(
            {
                monomial: coefficient / leading_coefficient
                for monomial, coefficient in polynomial.items()
            }
        )
    return basis


def parse_polynomial(text: str, variables):
    expression = sp.sympify(
        text.replace("^", "**"),
        locals={str(variable): variable for variable in variables},
    )
    return sp.Poly(expression, *variables, domain=sp.QQ)


def exact_shape(basis):
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for exact reconstruction")
    source = f"""ring r=0,(sigma,tau),dp;
option(redSB);
ideal I={','.join(render_polynomial(polynomial) for polynomial in basis)};
ideal G=std(I);
print("DP_DIMENSION="+string(dim(G)));
print("DP_VDIM="+string(vdim(G)));
ring l=0,(tau,sigma),lp;
option(redSB);
ideal J=imap(r,G);
ideal L=std(J);
print("LEX_SIZE="+string(size(L)));
print("LEX_DIMENSION="+string(dim(L)));
print("LEX_VDIM="+string(vdim(L)));
for (int i=1; i<=size(L); i++) {{ print("LEX_"+string(i)+"="+string(L[i])); }}
"""
    completed = subprocess.run(
        [singular, "-q"], input=source, text=True, capture_output=True, check=True
    )
    if completed.stderr.strip() or "   ?" in completed.stdout:
        raise AssertionError(completed.stdout + completed.stderr)
    values = {}
    lex_text = []
    for line in completed.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.startswith("LEX_") and key not in {
            "LEX_SIZE",
            "LEX_DIMENSION",
            "LEX_VDIM",
        }:
            lex_text.append(value)
        else:
            values[key] = int(value)
    assert values["DP_DIMENSION"] == values["LEX_DIMENSION"] == 0
    assert values["DP_VDIM"] == values["LEX_VDIM"] == 12
    tau, sigma = sp.symbols("tau sigma")
    lex = [parse_polynomial(text, (tau, sigma)) for text in lex_text]
    univariate = next(
        polynomial
        for polynomial in lex
        if polynomial.degree(tau) == 0 and polynomial.degree(sigma) == 12
    )
    factors = sp.factor_list(univariate.as_expr())[1]
    factor_degrees = [
        (sp.Poly(factor, sigma).degree(), multiplicity)
        for factor, multiplicity in factors
    ]
    return values, lex, univariate, factor_degrees


def serialize_sympy(polynomial):
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
    parser.add_argument("--primes", type=int, nargs="+", default=list(DEFAULT_PRIMES))
    parser.add_argument("--holdout-prime", type=int)
    parser.add_argument(
        "--artifact-directory",
        type=Path,
        default=Path("artifacts/generated-results"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    def path_for(prime):
        return args.artifact_directory / (
            f"degree_eight_order_five_chart_gf{prime}.json"
        )

    images = [load_basis(path_for(prime)) for prime in args.primes]
    lattice_basis = simultaneous_reconstruct_basis(args.primes, images)
    previous = simultaneous_reconstruct_basis(args.primes[:-1], images[:-1])
    stable_generators = (0, 1, 2)
    assert all(
        lattice_basis[index] == previous[index] for index in stable_generators
    ), "rational reconstruction of the generating subset is not stable"
    leading_monomials = []
    for polynomial_index in range(len(images[0])):
        candidates = [
            monomial
            for monomial in images[0][polynomial_index]
            if all(
                image[polynomial_index].get(monomial) == 1 for image in images
            )
        ]
        assert len(candidates) == 1
        leading_monomials.append(candidates[0])
    basis = complete_exact_basis(
        [lattice_basis[index] for index in stable_generators],
        leading_monomials,
    )
    for prime, image in zip(args.primes, images, strict=True):
        assert reduction(basis, prime) == image
    holdout = None
    if args.holdout_prime:
        image = load_basis(path_for(args.holdout_prime))
        assert reduction(basis, args.holdout_prime) == image
        holdout = {
            "prime": args.holdout_prime,
            "status": "exact agreement; excluded from CRT",
        }

    exact, lex, univariate, factor_degrees = exact_shape(basis)
    certificate = {
        "scope": "exact characteristic-zero degree-eight pivot-chart certificate",
        "build_primes": args.primes,
        "crt_modulus": str(prod(args.primes)),
        "crt_modulus_bits": prod(args.primes).bit_length(),
        "reconstruction_method": "simultaneous LLL common-denominator vectors",
        "stability": (
            "the first three normalized generators are identical after "
            "dropping the final build prime; exact Buchberger completion "
            "recovers the full five-polynomial basis"
        ),
        "stable_generator_indices_zero_based": list(stable_generators),
        "holdout": holdout,
        "degree_reverse_lexicographic_basis": [serialize(item) for item in basis],
        "exact_zero_scheme": {
            "dimension": exact["DP_DIMENSION"],
            "vector_space_dimension": exact["DP_VDIM"],
            "lexicographic_basis": [serialize_sympy(item) for item in lex],
            "primitive_sigma_polynomial": str(
                sp.primitive(univariate.as_expr(), expand=True)[1]
            ),
            "factor_degrees_over_Q": [list(item) for item in factor_degrees],
        },
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS: stable degree-eight reconstruction from {len(args.primes)} primes"
    )
    if holdout:
        print(f"PASS: independent holdout GF({args.holdout_prime})")
    print("EXACT: length=12 factor_degrees=", factor_degrees)


if __name__ == "__main__":
    main()

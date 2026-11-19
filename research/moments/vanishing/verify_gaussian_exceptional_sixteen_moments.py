#!/usr/bin/env python3
"""Replay the weighted slice and the two exact finite-field chart checks."""

from __future__ import annotations

import argparse
from functools import cache
import gzip
from hashlib import sha256
import json
from math import comb, factorial, isqrt, prod
from pathlib import Path
import platform
import subprocess
import tempfile

import flint


HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "gaussian_exceptional_sixteen_moment_certificate.json.gz"
PROOF = HERE / "GAUSSIAN_EXCEPTIONAL_SIXTEEN_MOMENTS.md"
NAMES = ("a", "b", "c", "d", "e", "f")
WEIGHTS = (3, 2, 1, 3, 2, 4)
PRIME = 32003
ORDERS = range(11, 17)
SPATIAL_ZERO = (0, 0, 0)
PARAMETER_ZERO = (0,) * len(NAMES)


@cache
def period(n, i, j, k):
    """Exact E[Z^i W^j Y^k (W+ZY-Z^3/2)^n]."""
    result = flint.fmpq(0)
    for quadratic in range(n + 1):
        for cubic in range(n - quadratic + 1):
            z_degree = i + quadratic + 3 * cubic
            w_degree = j + n - quadratic - cubic
            y_degree = k + quadratic
            if z_degree != w_degree or y_degree % 2:
                continue
            numerator = ((-1)**cubic * comb(n, cubic) * comb(n - cubic, quadratic)
                         * factorial(z_degree) * factorial(y_degree))
            denominator = 2**(cubic + y_degree // 2) * factorial(y_degree // 2)
            result += flint.fmpq(numerator, denominator)
    return result


def weighted_initials():
    parameters = flint.fmpq_mpoly_ctx.get(NAMES)
    space = flint.fmpq_mpoly_ctx.get(("z", "w", "y", *NAMES))
    z, w, y, *controls = space.gens()
    zero = parameters.constant(0)
    one = space.constant(1)
    pivots = (one, y, y**2, y**3, w*y, w*y**2, w**2*y, w**3, z, z*w**2)
    tangents = (
        w**2 + (11*y-y**3)/4,
        2*w*y-3*z+z*y**2,
        z*w+(y**2-3)/2,
        z*w*y+(y**3-7*y)/2,
        z**2*w-w*y-z,
        z**3+flint.fmpq(2, 9)*w*y**2+flint.fmpq(4, 9)*z*w**2,
    )

    def wick(n, polynomial):
        terms = {}
        for monomial, coefficient in polynomial.terms():
            value = period(n, *monomial[:3])
            if value:
                power = monomial[3:]
                terms[power] = terms.get(power, 0) + value * coefficient
        return parameters.from_dict(terms)

    def embed(polynomial):
        return space.from_dict({SPATIAL_ZERO + mon: coef
                                for mon, coef in polynomial.terms()})

    base = w+z*y-z**3/2
    # Compare the period formula with literal spatial powers at all used orders.
    def direct_wick(polynomial):
        value = flint.fmpq(0)
        for monomial, coefficient in polynomial.terms():
            i, j, k = monomial[:3]
            if i == j and k % 2 == 0:
                assert monomial[3:] == PARAMETER_ZERO
                value += coefficient * flint.fmpq(factorial(i) * factorial(k),
                                                   2**(k//2) * factorial(k//2))
        return value

    base_power = one
    for n in range(max(ORDERS) + 1):
        for marker in (*pivots, *tangents):
            assert direct_wick(base_power * marker) == wick(n, marker).to_dict().get(
                PARAMETER_ZERO, 0
            )
        if n:
            assert direct_wick(base_power) == 0
        base_power *= base

    assert all(not wick(n, tangent) for n in range(10) for tangent in tangents)
    assert all(not wick(n, marker) for n in range(10, 16)
               for marker in (*pivots, *tangents))
    jacobian = flint.fmpq_mat([
        [n * wick(n - 1, pivot).to_dict().get(PARAMETER_ZERO, 0)
         for pivot in pivots] for n in range(1, 11)
    ])
    determinant = jacobian.det()
    assert determinant and int(determinant.numerator) % PRIME
    inverse = jacobian.inv()

    # parts[r] is the weight-r perturbation. powers[k,r] is the
    # weight-r coefficient of its kth power. For k>=2 computing it
    # never uses parts[r], so the nonlinear moment terms are known
    # before the weight-r pivot correction.
    parts = {}
    powers = {}
    leading = {}
    for weight in range(1, max(ORDERS) + 1):
        for exponent in range(2, weight + 1):
            powers[exponent, weight] = sum(
                parts[i] * powers[exponent - 1, weight - i]
                for i in range(1, weight - exponent + 2)
            )
        nonlinear = {
            n: sum((comb(n, k) * wick(n - k, powers[k, weight])
                    for k in range(2, min(n, weight) + 1)), zero)
            for n in (*range(1, 11), *(n for n in ORDERS if weight <= n))
        }
        corrections = [
            -sum((inverse[i, n-1] * nonlinear[n] for n in range(1, 11)), zero)
            for i in range(len(pivots))
        ]
        parts[weight] = (
            sum(tangent * u for tangent, u, degree in zip(tangents, controls, WEIGHTS)
                if degree == weight)
            + sum(pivot * embed(value) for pivot, value in zip(pivots, corrections))
        )
        powers[1, weight] = parts[weight]
        for n in range(1, 11):
            assert n * wick(n - 1, parts[weight]) + nonlinear[n] == 0
        for n in ORDERS:
            if weight < n:
                assert not nonlinear[n], (weight, n, "unexpected lower term")
            elif weight == n:
                assert nonlinear[n]
                assert all(sum(a*b for a, b in zip(mon, WEIGHTS)) == n
                           for mon in nonlinear[n].to_dict())
                leading[n] = nonlinear[n]
    return determinant, leading


def verify_charts(leading, cross_check_singular):
    assert not any(PRIME % d == 0 for d in range(2, isqrt(PRIME) + 1))
    modular = flint.nmod_mpoly_ctx.get(NAMES, PRIME)

    def specialize(polynomial, value):
        terms = {}
        for mon, coefficient in polynomial.terms():
            if value == 0 and mon[2]:
                continue
            reduced_power = mon[:2] + (0,) + mon[3:]
            residue = int(coefficient.numerator) * pow(int(coefficient.denominator), -1, PRIME)
            terms[reduced_power] = (terms.get(reduced_power, 0) + residue) % PRIME
        return modular.from_dict(terms)

    def read_polynomial(encoded):
        terms = {}
        for monomial, coefficient in encoded:
            monomial = tuple(monomial)
            assert len(monomial) == len(NAMES) and monomial[2] == 0
            assert all(type(power) is int and power >= 0 for power in monomial)
            assert type(coefficient) is int and 1 <= coefficient < PRIME
            assert monomial not in terms
            terms[monomial] = coefficient
        return modular.from_dict(terms)

    def combination(encoded, generators):
        assert len(encoded) == len(generators)
        return sum(read_polynomial(row) * g for row, g in zip(encoded, generators))

    with gzip.open(CERTIFICATE, "rt", encoding="utf-8") as handle:
        certificate = json.load(handle)
    assert certificate["prime"] == PRIME
    zero_generators = [specialize(leading[n], 0) for n in ORDERS]
    one_generators = [specialize(leading[n], 1) for n in ORDERS]
    pure_powers = {}
    for target, multipliers in certificate["c_zero"]:
        assert len(target) == len(NAMES) and target[2] == 0
        assert all(type(exponent) is int and exponent >= 0 for exponent in target)
        indices = [i for i, exponent in enumerate(target) if exponent]
        assert len(indices) == 1 and indices[0] not in pure_powers
        pure_powers[indices[0]] = target[indices[0]]
        assert combination(multipliers, zero_generators) == modular.from_dict({
            tuple(target): 1
        })
    assert set(pure_powers) == {0, 1, 3, 4, 5}
    chart_input = "a,b,d,e,f\n" + str(PRIME) + "\n" + ",\n".join(
        str(polynomial).replace(" ", "") for polynomial in one_generators
    ) + "\n"
    with tempfile.TemporaryDirectory() as directory:
        input_path = Path(directory) / "affine.in"
        output_path = Path(directory) / "affine.out"
        input_path.write_text(chart_input, encoding="ascii")
        command = ["msolve", "-g", "2", "-l", "2",
                   "-f", str(input_path), "-o", str(output_path)]
        subprocess.run(command, check=True, capture_output=True, timeout=120)
        basis = [line.strip() for line in output_path.read_text().splitlines()
                 if line.strip() and not line.startswith("#")]
        assert basis == ["[1]:"], basis

    if cross_check_singular:
        lines = [f"ring R={PRIME},(a,b,c,d,e,f),dp;"]
        lines.extend(f"poly h{n}={leading[n]};" for n in ORDERS)
        lines.extend([
            "ideal I=h11,h12,h13,h14,h15,h16,c-1;",
            "ideal B=slimgb(I);",
            'if(size(B)==1 && B[1]==1){ print("UNIT"); }',
            "quit;",
        ])
        run = subprocess.run(["Singular", "-q"], input="\n".join(lines) + "\n",
                             text=True, capture_output=True, check=True, timeout=180)
        assert run.stdout.strip() == "UNIT" and not run.stderr, (run.stdout, run.stderr)
    return {
        "c_zero_powers": {NAMES[i]: degree for i, degree in sorted(pure_powers.items())},
        "c_one_reduced_groebner_basis": ["1"],
        "c_one_msolve_input_sha256": sha256(chart_input.encode("ascii")).hexdigest(),
        "c_one_msolve_options": ["-g", "2", "-l", "2"],
        "c_one_singular_cross_check": cross_check_singular,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    parser.add_argument("--cross-check-singular", action="store_true")
    args = parser.parse_args()
    determinant, leading = weighted_initials()
    chart_checks = verify_charts(leading, args.cross_check_singular)
    weighted_length, remainder = divmod(prod(ORDERS), prod(WEIGHTS))
    assert remainder == 0
    result = {
        "status": "passed",
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "sources": {path.name: sha256(path.read_bytes()).hexdigest()
                    for path in (Path(__file__), PROOF, CERTIFICATE)},
        "weights": dict(zip(NAMES, WEIGHTS)),
        "first_ten_moment_pivot_determinant": str(determinant),
        "initial_weights": list(leading),
        "initial_term_counts": {n: len(polynomial.to_dict()) for n, polynomial in leading.items()},
        "prime": PRIME,
        "modular_chart_checks": chart_checks,
        "weighted_complete_intersection_length": weighted_length,
        "scope": "Local transverse slice at the exceptional orbit only.",
    }
    output = json.dumps(result, indent=2) + "\n"
    if args.write_receipt:
        args.write_receipt.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()

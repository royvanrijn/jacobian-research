#!/usr/bin/env python3
"""Modular probe for the maximal derivative-gcd stratum in Case 2.

This script is intentionally a probe, not a characteristic-zero certificate.
If ``deg(gcd(C',G'))=7``, then the monic gcd is ``C'/8`` and the stratum is
cut out by the seven coefficients of ``remainder(G',C')`` in the three
high-layer parameters.  The exact algebraic-number-field standard basis is
currently expensive.  We reduce the degree-35 coefficient field at simple
linear primes, test small subsets over finite fields, and report stable unit
subsets that can be lifted by a later exact membership calculation.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import sympy as sp

from audit_case2_residue_strata import (
    default_replay_root,
    load_exact_core,
    singular_minpoly,
    solve_case2_through_j1,
)


def maximal_gcd_constraints(ec, C, G):
    divisor = ec.tder(C)
    remainder = list(ec.tder(G))
    divisor_degree = len(divisor) - 1
    leading = divisor[divisor_degree].d[(0,) * ec.NP]
    leading_inverse = leading.inv()
    for degree in range(len(remainder) - 1, divisor_degree - 1, -1):
        coefficient = remainder[degree] * leading_inverse
        shift = degree - divisor_degree
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
    if any(remainder[divisor_degree:]):
        raise AssertionError("polynomial division failed")
    constraints = tuple(remainder[:divisor_degree])
    if any(not polynomial for polynomial in constraints):
        raise AssertionError("expected seven nonzero remainder coefficients")
    return constraints


def rational_mod(value, prime: int) -> int:
    numerator = int(value.p) % prime
    denominator = int(value.q) % prime
    if denominator == 0:
        raise ZeroDivisionError
    return numerator * pow(denominator, -1, prime) % prime


def field_value_mod(coefficient, prime: int, root: int) -> int:
    value = 0
    for exponent in range(len(coefficient.p) - 1, -1, -1):
        value = (
            value * root + rational_mod(coefficient.p[exponent], prime)
        ) % prime
    return value


def simple_roots(ec, prime: int) -> tuple[int, ...]:
    coefficients = tuple(
        rational_mod(ec.MOD[index], prime)
        for index in range(len(ec.MOD))
    )

    def evaluate(root: int) -> int:
        value = 0
        for coefficient in reversed(coefficients):
            value = (value * root + coefficient) % prime
        return value

    derivative = tuple(
        index * coefficients[index] % prime
        for index in range(1, len(coefficients))
    )

    def evaluate_derivative(root: int) -> int:
        value = 0
        for coefficient in reversed(derivative):
            value = (value * root + coefficient) % prime
        return value

    return tuple(
        root
        for root in range(prime)
        if evaluate(root) == 0 and evaluate_derivative(root) != 0
    )


def polynomial_mod_source(polynomial, prime: int, root: int) -> str:
    terms = []
    names = ("r", "s", "h")
    for monomial, coefficient in polynomial.d.items():
        scalar = field_value_mod(coefficient, prime, root)
        if scalar == 0:
            continue
        factors = []
        for name, exponent in zip(names, monomial):
            if exponent:
                factors.append(name if exponent == 1 else f"{name}^{exponent}")
        if not factors:
            term = str(scalar)
        elif scalar == 1:
            term = "*".join(factors)
        else:
            term = f"{scalar}*" + "*".join(factors)
        terms.append(term)
    return "+".join(terms) or "0"


def finite_field_unit(
    constraints,
    subset: tuple[int, ...],
    prime: int,
    root: int,
    singular: str,
    timeout: int,
) -> bool:
    source = "\n".join(
        (
            f"ring R={prime},(r,s,h),dp;",
            "option(redSB);",
            "ideal I="
            + ",".join(
                polynomial_mod_source(constraints[index], prime, root)
                for index in subset
            )
            + ";",
            "ideal S=std(I);",
            'if(size(S)>0 && S[1]==1){print("UNIT");}',
            "quit;",
            "",
        )
    )
    with tempfile.TemporaryDirectory(prefix="jc2-gcd-mod-") as directory:
        path = Path(directory) / "probe.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    return "UNIT" in completed.stdout


def finite_field_profile(
    constraints,
    prime: int,
    root: int,
    singular: str,
    timeout: int,
    subset: tuple[int, ...] | None = None,
):
    selected = (
        constraints
        if subset is None
        else tuple(constraints[index] for index in subset)
    )
    source = "\n".join(
        (
            f"ring R={prime},(r,s,h),dp;",
            "option(redSB);",
            "ideal I="
            + ",".join(
                polynomial_mod_source(polynomial, prime, root)
                for polynomial in selected
            )
            + ";",
            "ideal S=std(I);",
            'print("BASIS_SIZE="+string(size(S)));',
            'print("DIMENSION="+string(dim(S)));',
            'if(size(S)>0 && S[1]==1){print("UNIT");}',
            'print("BASIS_BEGIN");',
            "S;",
            'print("BASIS_END");',
            "quit;",
            "",
        )
    )
    with tempfile.TemporaryDirectory(prefix="jc2-gcd-profile-") as directory:
        path = Path(directory) / "profile.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    values = {}
    basis_lines = []
    in_basis = False
    for line in completed.stdout.splitlines():
        if line == "BASIS_BEGIN":
            in_basis = True
            continue
        if line == "BASIS_END":
            in_basis = False
            continue
        if in_basis:
            basis_lines.append(line)
        if line.startswith(("BASIS_SIZE=", "DIMENSION=")):
            key, value = line.split("=", 1)
            values[key] = int(value)
    return {
        "prime": prime,
        "root": root,
        "unit_ideal": "UNIT" in completed.stdout,
        "basis_size": values["BASIS_SIZE"],
        "dimension": values["DIMENSION"],
        "basis": basis_lines,
    }


def exact_subset_profile(
    ec,
    constraints,
    subset: tuple[int, ...],
    singular: str,
    timeout: int,
):
    source = "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            'ring R=(0,u),(r,s,h),dp;',
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I="
            + ",".join(constraints[index].sing(3) for index in subset)
            + ";",
            "ideal S=nfmodStd(I);",
            'print("BASIS_SIZE="+string(size(S)));',
            'print("DIMENSION="+string(dim(S)));',
            'print("BASIS_BEGIN");',
            "S;",
            'print("BASIS_END");',
            "quit;",
            "",
        )
    )
    with tempfile.TemporaryDirectory(prefix="jc2-gcd-exact-") as directory:
        path = Path(directory) / "exact.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    return completed.stdout.strip().splitlines()


def finite_field_lift(
    constraints,
    subset: tuple[int, ...],
    prime: int,
    root: int,
    singular: str,
    timeout: int,
):
    source = "\n".join(
        (
            f"ring R={prime},(r,s,h),dp;",
            "ideal I="
            + ",".join(
                polynomial_mod_source(constraints[index], prime, root)
                for index in subset
            )
            + ";",
            "ideal T=r,s;",
            "matrix L=lift(I,T);",
            'print("LIFT_BEGIN");',
            "L;",
            'print("LIFT_END");',
            "quit;",
            "",
        )
    )
    with tempfile.TemporaryDirectory(prefix="jc2-gcd-lift-") as directory:
        path = Path(directory) / "lift.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    lines = completed.stdout.strip().splitlines()
    begin = lines.index("LIFT_BEGIN")
    end = lines.index("LIFT_END")
    return lines[begin + 1 : end]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path)
    parser.add_argument("--singular", default=shutil.which("Singular"))
    parser.add_argument("--prime-start", type=int, default=43)
    parser.add_argument("--prime-stop", type=int, default=300)
    parser.add_argument("--max-subset-size", type=int, default=7)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--exact-timeout", type=int, default=300)
    args = parser.parse_args()
    if not args.singular:
        raise RuntimeError("Singular is required")

    replay_root = (args.replay_root or default_replay_root()).resolve()
    ec = load_exact_core(replay_root)
    B, _, C, F, G, compatibility = solve_case2_through_j1(ec)
    if not compatibility:
        raise AssertionError("unused J1 compatibility equations are missing")
    constraints = maximal_gcd_constraints(ec, C, G)
    j0 = ec.tadd(
        ec.tmul(B, ec.tder(G)),
        ec.tscale(ec.tmul(ec.tder(C), F), ec.K(-1)),
    )
    j0_constraints = tuple(polynomial for polynomial in j0 if polynomial)
    if len(j0_constraints) != 18:
        raise AssertionError("unexpected J0 coefficient count")
    terminal_j0 = j0_constraints[-1]
    augmented_constraints = constraints + (terminal_j0,)
    degrees = tuple(max(map(sum, polynomial.d)) for polynomial in constraints)

    reductions = []
    for prime in sp.primerange(args.prime_start, args.prime_stop):
        try:
            roots = simple_roots(ec, prime)
        except ZeroDivisionError:
            continue
        for root in roots:
            reductions.append((prime, root))
    if not reductions:
        raise RuntimeError("no simple linear coefficient-field reductions")

    stable_subset = None
    tested_reductions = reductions[:3]
    profiles = tuple(
        finite_field_profile(
            constraints,
            prime,
            root,
            args.singular,
            args.timeout,
        )
        for prime, root in tested_reductions
    )
    for size in range(1, args.max_subset_size + 1):
        for subset in itertools.combinations(range(7), size):
            if all(
                finite_field_unit(
                    constraints,
                    subset,
                    prime,
                    root,
                    args.singular,
                    args.timeout,
                )
                for prime, root in tested_reductions
            ):
                stable_subset = subset
                break
        if stable_subset is not None:
            break

    stable_generator_subset = None
    expected_basis = ["S[1]=s", "S[2]=r"]
    for size in range(1, 8):
        for subset in itertools.combinations(range(7), size):
            subset_profiles = tuple(
                finite_field_profile(
                    constraints,
                    prime,
                    root,
                    args.singular,
                    args.timeout,
                    subset=subset,
                )
                for prime, root in tested_reductions
            )
            if all(
                profile["basis"] == expected_basis
                for profile in subset_profiles
            ):
                stable_generator_subset = subset
                break
        if stable_generator_subset is not None:
            break

    stable_augmented_unit_subset = None
    for size in range(1, 9):
        for subset in itertools.combinations(range(8), size):
            if 7 not in subset:
                continue
            subset_profiles = tuple(
                finite_field_profile(
                    augmented_constraints,
                    prime,
                    root,
                    args.singular,
                    args.timeout,
                    subset=subset,
                )
                for prime, root in tested_reductions
            )
            if all(profile["unit_ideal"] for profile in subset_profiles):
                stable_augmented_unit_subset = subset
                break
        if stable_augmented_unit_subset is not None:
            break

    print(
        json.dumps(
            {
                "claim_scope": "modular probe only",
                "constraint_count": len(constraints),
                "constraint_degrees": degrees,
                "field_reductions": tested_reductions,
                "full_ideal_profiles": profiles,
                "stable_unit_subset": stable_subset,
                "stable_unit_subset_size": (
                    None if stable_subset is None else len(stable_subset)
                ),
                "stable_generator_subset": stable_generator_subset,
                "stable_generator_subset_size": (
                    None
                    if stable_generator_subset is None
                    else len(stable_generator_subset)
                ),
                "terminal_j0_term_count": len(terminal_j0.d),
                "terminal_j0_degree": max(map(sum, terminal_j0.d)),
                "stable_augmented_unit_subset": (
                    stable_augmented_unit_subset
                ),
                "stable_augmented_unit_subset_size": (
                    None
                    if stable_augmented_unit_subset is None
                    else len(stable_augmented_unit_subset)
                ),
                "first_reduction_lift": (
                    None
                    if stable_generator_subset is None
                    else finite_field_lift(
                        constraints,
                        stable_generator_subset,
                        tested_reductions[0][0],
                        tested_reductions[0][1],
                        args.singular,
                        args.timeout,
                    )
                ),
                "uses_j1_compatibility": False,
                "uses_j0": False,
                "exact_subset_profile": (
                    None
                    if (
                        not args.exact
                        or stable_augmented_unit_subset is None
                    )
                    else exact_subset_profile(
                        ec,
                        augmented_constraints,
                        stable_augmented_unit_subset,
                        args.singular,
                        args.exact_timeout,
                    )
                ),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

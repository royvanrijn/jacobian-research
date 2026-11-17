#!/usr/bin/env python3
"""Discover small certificates on the Case-2 derivative-gcd degree-six row.

This is a finite-field discovery probe, not a characteristic-zero
certificate.  Since ``t|H`` and ``deg(C')=7``, a degree-six common divisor
can be reconstructed by writing

    C' = H(t) (t+v),  deg(H)=6,  H(0)=0.

Synthetic division expresses ``H`` in the original three parameters and
the one new parameter ``v``.  The probe imposes ``C'(0)=H(0)=0``, divides
``G'`` by this reconstructed ``H``, and searches for one additional ``J0``
coefficient that makes a stable finite-field unit ideal.  Any subset
found here must be lifted by a separate exact checker before it is used as a
mathematical result.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from audit_case2_residue_strata import (
    default_replay_root,
    load_exact_core,
    singular_minpoly,
    solve_case2_through_j1,
)
from probe_case2_maximal_gcd import (
    field_value_mod,
    simple_roots,
)


@dataclass(frozen=True)
class VPolynomial:
    """A polynomial in ``v`` with exact-core parameter-polynomial values."""

    coefficients: dict[int, object]


def vp_clean(coefficients):
    return VPolynomial(
        {degree: value for degree, value in coefficients.items() if value}
    )


def vp_constant(value):
    return vp_clean({0: value})


def vp_add(left, right):
    result = dict(left.coefficients)
    for degree, value in right.coefficients.items():
        result[degree] = result.get(degree, value * 0) + value
    return vp_clean(result)


def vp_scale(polynomial, scalar):
    return vp_clean(
        {
            degree: value * scalar
            for degree, value in polynomial.coefficients.items()
        }
    )


def vp_sub(left, right):
    return vp_add(left, vp_scale(right, -1))


def vp_mul(left, right):
    result = {}
    for left_degree, left_value in left.coefficients.items():
        for right_degree, right_value in right.coefficients.items():
            degree = left_degree + right_degree
            value = left_value * right_value
            result[degree] = result.get(degree, value * 0) + value
    return vp_clean(result)


def vp_times_v(polynomial):
    return VPolynomial(
        {
            degree + 1: value
            for degree, value in polynomial.coefficients.items()
        }
    )


def gcd6_constraints(ec, C, G):
    """Return ``H(0)`` and the six coefficients of ``G' mod H``."""

    cprime = ec.tder(C)
    gprime = ec.tder(G)
    if len(cprime) != 8 or len(gprime) != 12:
        raise AssertionError("unexpected derivative degrees")

    # Synthetic division by t+v.  The two first constraints impose the
    # triangular consequence C'(0)=0 and the required factor H(0)=0.
    quotient = [None] * 7
    quotient[6] = vp_constant(cprime[7])
    for degree in range(6, 0, -1):
        quotient[degree - 1] = vp_sub(
            vp_constant(cprime[degree]),
            vp_times_v(quotient[degree]),
        )

    leading = cprime[7].d[(0,) * ec.NP]
    leading_inverse = leading.inv()
    gcofactor = [None] * 6
    for total_degree in range(11, 5, -1):
        cofactor_degree = total_degree - 6
        known = vp_constant(ec.PP())
        for h_degree in range(6):
            other_degree = total_degree - h_degree
            if cofactor_degree < other_degree < 6:
                known = vp_add(
                    known,
                    vp_mul(
                        quotient[h_degree],
                        gcofactor[other_degree],
                    ),
                )
        gcofactor[cofactor_degree] = vp_scale(
            vp_sub(vp_constant(gprime[total_degree]), known),
            leading_inverse,
        )

    remainders = []
    for total_degree in range(6):
        product = vp_constant(ec.PP())
        for h_degree in range(7):
            other_degree = total_degree - h_degree
            if 0 <= other_degree < 6:
                product = vp_add(
                    product,
                    vp_mul(
                        quotient[h_degree],
                        gcofactor[other_degree],
                    ),
                )
        remainders.append(
            vp_sub(vp_constant(gprime[total_degree]), product)
        )
    constraints = (vp_constant(cprime[0]), quotient[0], *remainders)
    if any(not constraint.coefficients for constraint in constraints):
        raise AssertionError("unexpected zero gcd-six constraint")
    return constraints


def j0_constraints(ec, B, C, F, G):
    j0 = ec.tadd(
        ec.tmul(B, ec.tder(G)),
        ec.tscale(ec.tmul(ec.tder(C), F), ec.K(-1)),
    )
    return tuple(
        (degree, vp_constant(polynomial))
        for degree, polynomial in enumerate(j0)
        if polynomial
    )


def polynomial_mod_source(polynomial, prime: int, root: int) -> str:
    terms = []
    names = ("r", "s", "h")
    for v_degree, parameter_polynomial in polynomial.coefficients.items():
        for monomial, coefficient in parameter_polynomial.d.items():
            scalar = field_value_mod(coefficient, prime, root)
            if scalar == 0:
                continue
            factors = []
            for name, exponent in zip(names, monomial):
                if exponent:
                    factors.append(
                        name if exponent == 1 else f"{name}^{exponent}"
                    )
            if v_degree:
                factors.append("v" if v_degree == 1 else f"v^{v_degree}")
            monomial_source = "*".join(factors)
            if not monomial_source:
                term = str(scalar)
            elif scalar == 1:
                term = monomial_source
            else:
                term = f"{scalar}*{monomial_source}"
            terms.append(term)
    return "+".join(terms) or "0"


def finite_profile(
    constraints,
    prime: int,
    root: int,
    singular: str,
    timeout: int,
    show_basis: bool = False,
):
    source = "\n".join(
        (
            f"ring R={prime},(r,s,h,v),dp;",
            "option(redSB);",
            "ideal I="
            + ",".join(
                polynomial_mod_source(polynomial, prime, root)
                for polynomial in constraints
            )
            + ";",
            "ideal S=std(I);",
            'print("BASIS_SIZE="+string(size(S)));',
            'print("DIMENSION="+string(dim(S)));',
            'if(size(S)>0 && S[1]==1){print("UNIT");}',
            *(
                ('print("BASIS_BEGIN");', "S;", 'print("BASIS_END");')
                if show_basis
                else ()
            ),
            "quit;",
            "",
        )
    )
    try:
        with tempfile.TemporaryDirectory(
            prefix="jc2-gcd6-mod-"
        ) as directory:
            path = Path(directory) / "probe.sing"
            path.write_text(source)
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        return {
            "prime": prime,
            "root": root,
            "unit_ideal": False,
            "basis_size": None,
            "dimension": None,
            "basis": None,
            "timed_out": True,
        }
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    values = {}
    basis = []
    in_basis = False
    for line in completed.stdout.splitlines():
        if line == "BASIS_BEGIN":
            in_basis = True
            continue
        if line == "BASIS_END":
            in_basis = False
            continue
        if in_basis:
            basis.append(line)
        if line.startswith(("BASIS_SIZE=", "DIMENSION=")):
            key, value = line.split("=", 1)
            values[key] = int(value)
    return {
        "prime": prime,
        "root": root,
        "unit_ideal": "UNIT" in completed.stdout,
        "basis_size": values["BASIS_SIZE"],
        "dimension": values["DIMENSION"],
        "basis": basis if show_basis else None,
        "timed_out": False,
    }


def stable_unit(constraints, reductions, singular, timeout):
    return all(
        finite_profile(
            constraints,
            prime,
            root,
            singular,
            timeout,
        )["unit_ideal"]
        for prime, root in reductions
    )


def constraint_term_count(polynomial) -> int:
    return sum(
        len(parameter_polynomial.d)
        for parameter_polynomial in polynomial.coefficients.values()
    )


def polynomial_singular(polynomial) -> str:
    terms = []
    for v_degree, parameter_polynomial in polynomial.coefficients.items():
        coefficient = parameter_polynomial.sing(3)
        if v_degree == 0:
            terms.append(f"({coefficient})")
        elif v_degree == 1:
            terms.append(f"({coefficient})*v")
        else:
            terms.append(f"({coefficient})*v^{v_degree}")
    return "+".join(terms) or "0"


def exact_profile(
    ec,
    constraints,
    singular: str,
    timeout: int,
    backend: str,
):
    source = "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            "ring R=(0,u),(r,s,h,v),dp;",
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I="
            + ",".join(
                polynomial_singular(polynomial)
                for polynomial in constraints
            )
            + ";",
            f"ideal S={backend}(I);",
            'print("BASIS_SIZE="+string(size(S)));',
            'print("DIMENSION="+string(dim(S)));',
            'if(size(S)>0 && S[1]==1){print("UNIT");}',
            "quit;",
            "",
        )
    )
    metadata = {
        "singular_input_sha256": hashlib.sha256(
            source.encode()
        ).hexdigest(),
        "constraint_term_counts": [
            sum(
                len(parameter_polynomial.d)
                for parameter_polynomial in polynomial.coefficients.values()
            )
            for polynomial in constraints
        ],
        "constraint_parameter_degrees": [
            max(
                v_degree + sum(monomial)
                for v_degree, parameter_polynomial
                in polynomial.coefficients.items()
                for monomial in parameter_polynomial.d
            )
            for polynomial in constraints
        ],
    }
    try:
        with tempfile.TemporaryDirectory(
            prefix="jc2-gcd6-exact-"
        ) as directory:
            path = Path(directory) / "probe.sing"
            path.write_text(source)
            completed = subprocess.run(
                [singular, "-q", str(path)],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        return {
            **metadata,
            "unit_ideal": None,
            "basis_size": None,
            "dimension": None,
            "timed_out": True,
        }
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr)
    values = {}
    for line in completed.stdout.splitlines():
        if line.startswith(("BASIS_SIZE=", "DIMENSION=")):
            key, value = line.split("=", 1)
            values[key] = int(value)
    return {
        "unit_ideal": "UNIT" in completed.stdout,
        "basis_size": values["BASIS_SIZE"],
        "dimension": values["DIMENSION"],
        "timed_out": False,
        **metadata,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path)
    parser.add_argument("--singular", default=shutil.which("Singular"))
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument(
        "--exact-backend",
        choices=("std", "nfmodStd"),
        default="nfmodStd",
    )
    args = parser.parse_args()
    if not args.singular:
        raise RuntimeError("Singular is required")

    replay_root = (args.replay_root or default_replay_root()).resolve()
    ec = load_exact_core(replay_root)
    B, _, C, F, G, compatibility = solve_case2_through_j1(ec)
    if not compatibility:
        raise AssertionError("the unused J1 compatibility system is missing")
    gcd_constraints = gcd6_constraints(ec, C, G)
    j0 = j0_constraints(ec, B, C, F, G)

    requested_primes = (47, 61, 73)
    reductions = tuple(
        (prime, simple_roots(ec, prime)[0])
        for prime in requested_primes
    )
    base_profiles = tuple(
        finite_profile(
            gcd_constraints,
            prime,
            root,
            args.singular,
            args.timeout,
            show_basis=index == 0,
        )
        for index, (prime, root) in enumerate(reductions)
    )

    stable_j0_degrees = []
    for degree, constraint in j0:
        if stable_unit(
            (*gcd_constraints, constraint),
            reductions,
            args.singular,
            args.timeout,
        ):
            stable_j0_degrees.append(degree)

    minimal_subset = None
    minimum_cardinality_subset = None
    exact = None
    if stable_j0_degrees:
        degree = stable_j0_degrees[-1]
        j0_constraint = dict(j0)[degree]
        selected = list(range(len(gcd_constraints)))
        for index in tuple(selected):
            candidate = [
                gcd_constraints[position]
                for position in selected
                if position != index
            ]
            if stable_unit(
                (*candidate, j0_constraint),
                reductions,
                args.singular,
                args.timeout,
            ):
                selected.remove(index)
        minimal_subset = {
            "gcd_constraint_indices": selected,
            "j0_degree": degree,
        }
        ordered_subsets = []
        for size in range(4, len(gcd_constraints) + 1):
            candidates = sorted(
                itertools.combinations(range(len(gcd_constraints)), size),
                key=lambda subset: sum(
                    constraint_term_count(gcd_constraints[index])
                    for index in subset
                ),
            )
            for subset in candidates:
                candidate_constraints = tuple(
                    gcd_constraints[index] for index in subset
                ) + (j0_constraint,)
                first_prime, first_root = reductions[0]
                if not finite_profile(
                    candidate_constraints,
                    first_prime,
                    first_root,
                    args.singular,
                    min(args.timeout, 10),
                )["unit_ideal"]:
                    continue
                if stable_unit(
                    candidate_constraints,
                    reductions[1:],
                    args.singular,
                    min(args.timeout, 10),
                ):
                    ordered_subsets.append(subset)
            if ordered_subsets:
                best = ordered_subsets[0]
                minimum_cardinality_subset = {
                    "gcd_constraint_indices": list(best),
                    "gcd_constraint_term_counts": [
                        constraint_term_count(gcd_constraints[index])
                        for index in best
                    ],
                    "j0_degree": degree,
                }
                break
        if args.exact:
            exact_indices = (
                minimum_cardinality_subset["gcd_constraint_indices"]
                if minimum_cardinality_subset is not None
                else selected
            )
            exact = exact_profile(
                ec,
                tuple(
                    gcd_constraints[index] for index in exact_indices
                )
                + (j0_constraint,),
                args.singular,
                args.timeout,
                args.exact_backend,
            )

    print(
        json.dumps(
            {
                "reconstruction": "C'=H*(t+v), deg(H)=6, H(0)=0",
                "gcd_constraint_count": len(gcd_constraints),
                "j0_degrees_tested": [degree for degree, _ in j0],
                "reductions": [list(reduction) for reduction in reductions],
                "gcd_only_profiles": base_profiles,
                "stable_single_j0_degrees": stable_j0_degrees,
                "greedy_stable_subset": minimal_subset,
                "minimum_cardinality_stable_subset": (
                    minimum_cardinality_subset
                ),
                "characteristic_zero_certificate": exact,
            },
            sort_keys=True,
        )
    )
    print("CASE2_GCD6_MODULAR_PROBE_PASS")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify the primitive common-power carrier Wronskian theorem.

Suppose a source carrier has coordinates ``q=y``, ``v=x*y^k`` and common
leading edge

    P = a*q^(-k*m)*c(v)^m + lower terms,
    R = b*q^(-k*n)*c(v)^n + lower terms,

where ``a*b`` is nonzero, ``gcd(m,n)=1``, ``deg(c)=k`` and the nonzero
constant Jacobian is measured against

    dx wedge dy = -q^(-k) dq wedge dv.

If the root multiplicities of ``c`` have gcd one, a unimodular target
monomial and ordinary target shears force the first nonshear coefficient at

    delta = k*(m+n-1)+1.

Writing ``rad(c)`` for the squarefree radical, that coefficient is

    H = rad(c)*N/c^(m+n),       deg(N)=k-r-1,

where ``r`` is the number of distinct roots.  Its coefficients satisfy the
universal, ``m,n``-independent equation

    k*(rad(c)*N)'-(k-1)*(c'/c)*rad(c)*N = nonzero constant.

Consequently the carrier residue is a three-point map

    g = rad(c)^k*N^k/c^(k-1).

This checker verifies the symbolic identities, the local divisor argument,
the resulting universal passport, and a bounded census of multiplicity
partitions.  The proof is algebraic; the bounded census is a regression test,
not the source of the theorem.
"""

from __future__ import annotations

import argparse
from functools import reduce
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Iterator

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_common_power_carrier_wronskian.json"
)


def multiplicity_gcd(partition: tuple[int, ...]) -> int:
    """Return the gcd of a nonempty root-multiplicity partition."""

    return reduce(gcd, partition)


def integer_partitions(
    total: int, maximum: int | None = None
) -> Iterator[tuple[int, ...]]:
    """Yield the partitions of ``total`` in nonincreasing order."""

    if total == 0:
        yield ()
        return
    upper = min(total, maximum if maximum is not None else total)
    for first in range(upper, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


def bezout_target_monomial(m: int, n: int) -> dict[str, object]:
    """Find nonnegative ``A,B`` with ``B*n-A*m=1``."""

    if m < 1 or n < 1 or gcd(m, n) != 1:
        raise ValueError("m and n must be positive and coprime")
    if m == 1:
        b = 1
    else:
        b = pow(n, -1, m)
    a = (b * n - 1) // m
    if min(a, b) < 0 or b * n - a * m != 1:
        raise AssertionError("Bezout target monomial construction failed")

    matrix = ((a, -b), (n, -m))
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if determinant != 1:
        raise AssertionError("target exponent matrix is not unimodular")
    return {
        "m": m,
        "n": n,
        "A": a,
        "B": b,
        "identity": f"{b}*{n}-{a}*{m}=1",
        "target_exponent_matrix": [list(row) for row in matrix],
        "determinant": determinant,
        "pi": f"P^{a}/R^{b}",
        "h": f"P^{n}/R^{m}",
    }


def symbolic_reduction_audit() -> dict[str, object]:
    """Check the universal Wronskian and logarithmic-derivative identities."""

    v = sp.symbols("v")
    k, power_sum = sp.symbols("k power_sum", integer=True, positive=True)
    c = sp.Function("c")(v)
    carrier_numerator = sp.Function("carrier_numerator")(v)
    delta = k * (power_sum - 1) + 1
    h_coefficient = carrier_numerator * c ** (-power_sum)

    reduced = sp.simplify(
        c**power_sum
        * (
            k * sp.diff(h_coefficient, v)
            + delta * sp.diff(c, v) * h_coefficient / c
        )
    )
    universal_operator = (
        k * sp.diff(carrier_numerator, v)
        - (k - 1) * sp.diff(c, v) * carrier_numerator / c
    )
    if sp.simplify(reduced - universal_operator) != 0:
        raise AssertionError("the Wronskian did not reduce universally")

    residue = carrier_numerator**k * c ** (1 - k)
    logarithmic = sp.simplify(
        carrier_numerator * sp.diff(residue, v) / residue
    )
    if sp.simplify(logarithmic - universal_operator) != 0:
        raise AssertionError("the residue derivative identity changed")

    mu = sp.symbols("mu", integer=True, positive=True)
    local_h_order = 1 - mu * power_sum
    local_leading_coefficient = sp.expand(
        k * local_h_order + delta * mu
    )
    expected_local_coefficient = k - mu * (k - 1)
    if sp.simplify(
        local_leading_coefficient - expected_local_coefficient
    ) != 0:
        raise AssertionError("the local indicial coefficient changed")

    return {
        "forced_descent": "delta=k*(m+n-1)+1",
        "finite_root_order": "ord_alpha(H)=1-mu_alpha*(m+n)",
        "finite_indicial_coefficient": "k-mu_alpha*(k-1)",
        "exceptional_resonance": "k=2 and mu_alpha=2",
        "coefficient_shape": "H=rad(c)*N/c^(m+n)",
        "universal_operator": (
            "k*(rad(c)*N)'-(k-1)*(c'/c)*rad(c)*N=constant_nonzero"
        ),
        "residue_map": "g=rad(c)^k*N^k/c^(k-1)",
        "derivative_identity": (
            "rad(c)*N*g'/g="
            "k*(rad(c)*N)'-(k-1)*(c'/c)*rad(c)*N"
        ),
    }


def fixed_carrier_linear_system(
    c: sp.Expr, k: int
) -> tuple[sp.Matrix, tuple[sp.Symbol, ...], sp.Symbol]:
    """Build the coefficient matrix for the fixed-``c`` carrier equation."""

    v = sp.symbols("v")
    coefficients = sp.symbols(f"a0:{k}")
    kappa = sp.symbols("kappa")
    carrier_numerator = sum(
        coefficients[index] * v**index for index in range(k)
    )
    equation = sp.Poly(
        sp.expand(
            k * c * sp.diff(carrier_numerator, v)
            - (k - 1) * sp.diff(c, v) * carrier_numerator
            - kappa * c
        ),
        v,
    )
    scalar_equations = [
        equation.coeff_monomial(v**degree) for degree in range(2 * k - 1)
    ]
    matrix, right_side = sp.linear_eq_to_matrix(
        scalar_equations, (*coefficients, kappa)
    )
    if any(entry != 0 for entry in right_side):
        raise AssertionError("the fixed-carrier system ceased to be homogeneous")
    return matrix, coefficients, kappa


def fixed_carrier_kernel_audit() -> dict[str, object]:
    """Exhibit the linear reconstruction interface and its uniqueness gate."""

    v = sp.symbols("v")
    k = 5
    forced_c = sp.expand(v * (v - 1) ** 2 * (v**2 - 3 * v + 3))
    generic_c = sp.expand(v * (v - 1) ** 2 * (v - 2) * (v - 3))

    forced_matrix, coefficients, _ = fixed_carrier_linear_system(forced_c, k)
    forced_kernel = forced_matrix.nullspace()
    if len(forced_kernel) != 1:
        raise AssertionError("the forced F2 carrier kernel is no longer a line")
    kernel_vector = forced_kernel[0]
    forced_a = sp.factor(
        sum(kernel_vector[index] * v**index for index in range(k))
    )
    expected_radical = v * (v - 1) * (v**2 - 3 * v + 3)
    if sp.simplify(forced_a / expected_radical).diff(v) != 0:
        raise AssertionError("the forced F2 kernel lost its radical numerator")

    generic_matrix, _, _ = fixed_carrier_linear_system(generic_c, k)
    if generic_matrix.nullspace():
        raise AssertionError("the generic comparison carrier became admissible")

    if forced_matrix.cols != len(coefficients) + 1:
        raise AssertionError("the fixed-carrier unknown census changed")

    return {
        "linear_system": (
            "Phi_c(A,kappa)=k*c*A'-(k-1)*c'*A-kappa*c=0, "
            "deg(A)<=k-1"
        ),
        "unknown_count": "k coefficients of A plus kappa",
        "uniqueness": (
            "for primitive c, ker(k*c*d/dv-(k-1)*c')=0; hence any "
            "nonzero projective solution (A,kappa) is unique"
        ),
        "forced_k5_example": {
            "c": "v*(v-1)^2*(v^2-3*v+3)",
            "matrix_shape": list(forced_matrix.shape),
            "rank": forced_matrix.rank(),
            "kernel_dimension": len(forced_kernel),
            "A_up_to_scale": "v*(v-1)*(v^2-3*v+3)",
        },
        "generic_k5_comparison": {
            "c": "v*(v-1)^2*(v-2)*(v-3)",
            "matrix_shape": list(generic_matrix.shape),
            "rank": generic_matrix.rank(),
            "kernel_dimension": len(generic_matrix.nullspace()),
        },
    }


def passport_for_partition(
    k: int, partition: tuple[int, ...]
) -> dict[str, object]:
    """Return the universal carrier passport for a primitive partition."""

    if k < 3 or sum(partition) != k:
        raise ValueError("partition must have sum k>=3")
    if multiplicity_gcd(partition) != 1:
        raise ValueError("the clean carrier theorem requires a primitive partition")

    distinct_roots = len(partition)
    if distinct_roots == k:
        raise ValueError("the squarefree degree-k carrier is obstructed")
    numerator_degree = k - distinct_roots - 1
    if numerator_degree < 0:
        raise AssertionError("the infinity degree gate failed")

    simple_roots = sum(mu == 1 for mu in partition)
    zero_profile = [k] * numerator_degree + [1] * simple_roots
    pole_profile = [
        (k - 1) * mu - k for mu in partition if mu >= 2
    ]
    degree_from_zeros = sum(zero_profile)
    degree_from_poles = sum(pole_profile)
    if degree_from_zeros != degree_from_poles:
        raise AssertionError("zero and pole degrees disagree")
    degree = degree_from_zeros

    infinity_index = k - 2
    third_unramified = degree - infinity_index
    if third_unramified < 0:
        raise AssertionError("the third fiber is larger than the cover")
    third_profile = [infinity_index] + [1] * third_unramified

    ramification_zero = sum(index - 1 for index in zero_profile)
    ramification_pole = sum(index - 1 for index in pole_profile)
    ramification_third = sum(index - 1 for index in third_profile)
    total_ramification = (
        ramification_zero + ramification_pole + ramification_third
    )
    if total_ramification != 2 * degree - 2:
        raise AssertionError("Riemann-Hurwitz failed")

    parity = {
        "over_0": "even" if ramification_zero % 2 == 0 else "odd",
        "over_infinity": (
            "even" if ramification_pole % 2 == 0 else "odd"
        ),
        "over_third_value": (
            "even" if ramification_third % 2 == 0 else "odd"
        ),
    }
    if (k % 2 == 1) != all(value == "even" for value in parity.values()):
        raise AssertionError("the universal alternating-group parity changed")

    return {
        "k": k,
        "root_multiplicity_partition": list(partition),
        "distinct_roots": distinct_roots,
        "simple_roots": simple_roots,
        "numerator_degree": numerator_degree,
        "residue_degree": degree,
        "passport": {
            "over_0": zero_profile,
            "over_infinity": pole_profile,
            "over_third_value": third_profile,
        },
        "local_orders_at_c_roots": [
            k - (k - 1) * mu for mu in partition
        ],
        "ramification_check": {
            "over_0": ramification_zero,
            "over_infinity": ramification_pole,
            "over_third_value": ramification_third,
            "total": total_ramification,
            "two_degree_minus_two": 2 * degree - 2,
        },
        "branch_cycle_parity": parity,
        "geometric_monodromy_constraint": (
            "contained_in_alternating_group"
            if k % 2 == 1
            else "not_contained_in_alternating_group"
        ),
    }


def common_power_audit(max_power: int = 14) -> dict[str, object]:
    """Stress the Bezout target change and descent formulas."""

    checked = 0
    for m in range(1, max_power + 1):
        for n in range(1, max_power + 1):
            if gcd(m, n) != 1:
                continue
            row = bezout_target_monomial(m, n)
            a = int(row["A"])
            b = int(row["B"])
            if b * n - a * m != 1:
                raise AssertionError("bounded Bezout audit failed")
            for k in range(3, 13):
                delta = k * (m + n - 1) + 1
                if gcd(k, delta) != 1:
                    raise AssertionError("the carrier target ray is not primitive")
                removable = [k * index for index in range(1, m + n)]
                if removable and removable[-1] >= delta:
                    raise AssertionError("a target shear reached the Jacobian row")
            checked += 1

    samples = [
        bezout_target_monomial(m, n)
        for m, n in ((2, 3), (3, 5), (4, 7), (5, 8))
    ]
    return {
        "bounded_power_pairs_checked": checked,
        "maximum_m_or_n": max_power,
        "sample_target_monomials": samples,
        "leading_target_orders": {"pi": "k", "h": "0"},
        "target_two_form_order": "k*(m+n)",
        "forced_descent": "k*(m+n-1)+1",
        "removable_shears": "k,2k,...,k*(m+n-1)",
        "transverse_index": 1,
    }


def partition_census(max_k: int = 24) -> dict[str, object]:
    """Enumerate a bounded partition census as a regression audit."""

    rows: list[dict[str, int]] = []
    admitted_total = 0
    for k in range(3, max_k + 1):
        total = primitive = admitted = imprimitive = squarefree_obstructed = 0
        for partition in integer_partitions(k):
            total += 1
            if multiplicity_gcd(partition) != 1:
                imprimitive += 1
                continue
            primitive += 1
            if len(partition) == k:
                squarefree_obstructed += 1
                continue
            passport_for_partition(k, partition)
            admitted += 1
        admitted_total += admitted
        rows.append(
            {
                "k": k,
                "all_partitions": total,
                "primitive_partitions": primitive,
                "admitted_passports": admitted,
                "squarefree_obstructed": squarefree_obstructed,
                "imprimitive_deferred": imprimitive,
            }
        )

    return {
        "maximum_k": max_k,
        "admitted_profiles_checked": admitted_total,
        "rows": rows,
        "logical_role": (
            "bounded regression census; the theorem follows from the symbolic "
            "divisor and derivative identities"
        ),
    }


def f2_specializations() -> list[dict[str, object]]:
    """Recover the two F2 carrier passports at k=5."""

    squarefree_cofactor = passport_for_partition(5, (2, 1, 1, 1))
    double_cofactor = passport_for_partition(5, (2, 2, 1))
    if squarefree_cofactor["passport"] != {
        "over_0": [1, 1, 1],
        "over_infinity": [3],
        "over_third_value": [3],
    }:
        raise AssertionError("the F2 cyclic cubic passport changed")
    if double_cofactor["passport"] != {
        "over_0": [5, 1],
        "over_infinity": [3, 3],
        "over_third_value": [3, 1, 1, 1],
    }:
        raise AssertionError("the F2 degree-six passport changed")
    return [squarefree_cofactor, double_cofactor]


def build_payload() -> dict[str, object]:
    return {
        "schema": "plane-jc.common-power-carrier-wronskian.v1",
        "status": "proved-general-reduction-to-three-point-hurwitz-data",
        "hypotheses": {
            "carrier_chart": "q=y, v=x*y^k, x=v*q^-k",
            "two_form": "dx wedge dy=-q^-k*dq wedge dv",
            "leading_edge": [
                "P=a*q^(-k*m)*c(v)^m+lower",
                "R=b*q^(-k*n)*c(v)^n+lower",
            ],
            "nonvanishing": "a*b*Jacobian!=0",
            "coprime_common_powers": "gcd(m,n)=1",
            "carrier_degree": "deg(c)=k>=3",
            "primitive_root_partition": "gcd(mu_1,...,mu_r)=1",
        },
        "target_monomial_and_descent": common_power_audit(),
        "symbolic_reduction": symbolic_reduction_audit(),
        "universal_passport": {
            "notation": {
                "partition": "mu_1+...+mu_r=k",
                "s": "number of mu_i equal to one",
                "L": "k-r-1",
                "f": "k*L+s",
                "pole_parts": "p_i=(k-1)*mu_i-k for mu_i>=2",
            },
            "squarefree_gate": "r=k would force deg(N)=-1 and is impossible",
            "coefficient_polynomial": (
                "N is squarefree, coprime to c, and has degree L=k-r-1"
            ),
            "passport": {
                "over_0": "(k^L,1^s)",
                "over_infinity": "(p_i : mu_i>=2)",
                "over_third_value": "(k-2,1^(f-k+2))",
            },
            "third_source_point": "v=infinity",
            "monodromy_parity": {
                "odd_k": (
                    "all three branch permutations are even, so geometric "
                    "monodromy is contained in A_f"
                ),
                "even_k": (
                    "the branch permutation at the third value is odd, so "
                    "geometric monodromy is not contained in A_f"
                ),
                "arithmetic_boundary": (
                    "this constrains geometric monodromy only; Galois action "
                    "on constants may enlarge the arithmetic group"
                ),
            },
            "simple_double_corollary": {
                "hypothesis": (
                    "t double roots and s=k-2*t>=1 simple roots"
                ),
                "numerator_degree": "t-1",
                "residue_degree": "t*(k-2)",
                "passport": (
                    "(k^(t-1),1^(k-2*t)) | ((k-2)^t) | "
                    "(k-2,1^((t-1)*(k-2)))"
                ),
            },
        },
        "fixed_carrier_linear_reconstruction": fixed_carrier_kernel_audit(),
        "partition_regression_census": partition_census(),
        "f2_k5_specializations": f2_specializations(),
        "excluded_or_deferred_loci": [
            {
                "locus": "k=2",
                "reason": "a double root makes the local and infinity indicial coefficient resonant",
            },
            {
                "locus": "gcd(mu_1,...,mu_r)>1",
                "reason": (
                    "rational homogeneous coefficients can occur at descents "
                    "smaller than k and are not ordinary target shears"
                ),
            },
        ],
        "claim_boundary": (
            "This is a necessary carrier reduction.  It assigns an explicit "
            "Hurwitz passport to every solution of the universal operator, "
            "but it neither proves that a coefficient polynomial c,N exists "
            "for every partition nor realizes the lower Laurent system."
        ),
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/"
            "verify_common_power_carrier_wronskian.py"
        ),
        "software": {"python": "sympy exact symbolic arithmetic"},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display = artifact.relative_to(ROOT)
        except ValueError:
            display = artifact
        print(f"WROTE {display}")
    else:
        expected = json.loads(artifact.read_text())
        if expected != payload:
            raise AssertionError(
                "the pinned common-power carrier artifact is stale; inspect "
                "the change before using --refresh"
            )

    print("COMMON_POWER_TARGET_MONOMIAL_PASS")
    print("PRIMITIVE_PRETARGET_SHEAR_GATE_PASS")
    print("UNIVERSAL_CARRIER_WRONSKIAN_PASS")
    print("UNIVERSAL_CARRIER_BELYI_PASSPORT_PASS")
    print("F2_K5_PASSPORT_SPECIALIZATIONS_PASS")
    print(f"COMMON_POWER_CARRIER_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()

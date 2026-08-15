#!/usr/bin/env python3
"""Exact Brumer--Kramer local/class-group lower-bound audit for R20.

This script does not compute the full class group. It builds the cubic subfield
of the 2-division field, computes its maximal order and the exact local term in
the Brumer--Kramer 2-Selmer bound, and combines that with the pinned
rank-at-least-20 certificate. The resulting lower bound on the cubic field's
class-group 2-rank is unconditional.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.polys.numberfields import prime_decomp, round_two

from structural_search import two_division_cubic


ROOT = Path(__file__).resolve().parents[2]
R20_NEAR_MISS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "fermigier_rank20_near_miss_v1.json"
)
R20_NEAR_MISS_SHA256 = (
    "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1"
)
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic_r20_brumer_kramer.json"
PRIMARY_SOURCE = "https://arxiv.org/abs/1606.07178"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/audit_r20_brumer_kramer.py"
)
X = sp.symbols("x")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(record: dict[str, Any]) -> str:
    encoded = json.dumps(
        record, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _factor_record(value: int) -> list[dict[str, int]]:
    return [
        {"prime": int(prime), "exponent": int(exponent)}
        for prime, exponent in sorted(sp.factorint(abs(int(value))).items())
    ]


def _factor_map(value: int) -> dict[int, int]:
    return {
        int(prime): int(exponent)
        for prime, exponent in sp.factorint(abs(int(value))).items()
    }


def weierstrass_discriminant(ainvs: Sequence[int]) -> int:
    if len(ainvs) != 5:
        raise ValueError("five generalized Weierstrass coefficients are required")
    a1, a2, a3, a4, a6 = (int(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    return -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def monic_cubic_for_leading_scaled_root(
    primitive_coefficients_ascending: Sequence[int],
) -> tuple[int, ...]:
    """Return the monic polynomial for z=a*x from a*x^3+b*x^2+c*x+d."""

    if len(primitive_coefficients_ascending) != 4:
        raise ValueError("a cubic requires four ascending coefficients")
    d, c, b, a = (int(value) for value in primitive_coefficients_ascending)
    if a == 0:
        raise ValueError("the leading cubic coefficient must be nonzero")
    # a^2 f(z/a) = z^3 + b z^2 + a c z + a^2 d.
    coefficients = (a * a * d, a * c, b, 1)
    content = 0
    for value in coefficients:
        content = sp.igcd(content, abs(value))
    if content != 1:
        raise AssertionError("the leading-scaled polynomial unexpectedly lost primitivity")
    return coefficients


def _poly_from_ascending(coefficients: Sequence[int]) -> sp.Poly:
    expression = sum(
        int(value) * X**degree for degree, value in enumerate(coefficients)
    )
    return sp.Poly(expression, X, domain=sp.ZZ)


def _modular_irreducibility_witness(poly: sp.Poly) -> dict[str, Any]:
    for prime in sp.primerange(2, 500):
        reduced = sp.Poly(poly.as_expr(), X, modulus=prime)
        if reduced.is_irreducible:
            coefficients = [int(value) % prime for value in reduced.all_coeffs()]
            return {
                "prime": int(prime),
                "coefficients_descending_mod_prime": coefficients,
                "root_residues": [
                    residue
                    for residue in range(prime)
                    if int(reduced.eval(residue)) % prime == 0
                ],
                "proof": (
                    "a cubic over a field is irreducible iff it has no root; "
                    "the listed exhaustive residue test has no root"
                ),
            }
    raise RuntimeError("no small modular irreducibility witness was found")


def _load_input() -> dict[str, Any]:
    observed = sha256_file(R20_NEAR_MISS)
    if observed != R20_NEAR_MISS_SHA256:
        raise RuntimeError(
            f"pinned R20 input changed: {observed} != {R20_NEAR_MISS_SHA256}"
        )
    record = json.loads(R20_NEAR_MISS.read_text(encoding="utf-8"))
    if record["family"]["adapter_parameter"] != "28917/20":
        raise RuntimeError("the pinned R20 parameter changed")
    if "at least 20 independent rational points" not in record["conclusion"]:
        raise RuntimeError("the pinned rank-at-least-20 conclusion changed")
    return record


def _serialize_prime_decomposition(primes: Iterable[Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": int(prime.p),
            "ramification_index": int(prime.e),
            "residue_degree": int(prime.f),
            "alpha": str(prime.alpha),
        }
        for prime in primes
    ]


def build_record() -> dict[str, Any]:
    source = _load_input()
    ainvs = tuple(int(value) for value in source["global_curve"]["minimal_model"])
    conductor = int(source["global_curve"]["conductor"])
    stored_discriminant = int(source["global_curve"]["minimal_discriminant"])
    discriminant = weierstrass_discriminant(ainvs)
    if discriminant != stored_discriminant:
        raise AssertionError("the minimal discriminant did not replay from the model")

    primitive_two_division = two_division_cubic(ainvs)
    monic_coefficients = monic_cubic_for_leading_scaled_root(primitive_two_division)
    cubic = _poly_from_ascending(monic_coefficients)
    irreducibility = _modular_irreducibility_witness(cubic)
    if not cubic.is_irreducible:
        raise AssertionError("the 2-division cubic became reducible over Q")

    polynomial_discriminant = int(sp.discriminant(cubic.as_expr(), X))
    maximal_order, field_discriminant = round_two(cubic)
    field_discriminant = int(field_discriminant)
    quotient = polynomial_discriminant // field_discriminant
    field_index = sp.integer_nthroot(abs(quotient), 2)[0]
    if quotient != field_index * field_index:
        raise AssertionError("polynomial discriminant / field discriminant is not a square")
    if field_index != int(maximal_order.denom):
        raise AssertionError("the maximal-order index and denominator disagree")

    conductor_factor = _factor_map(conductor)
    discriminant_factor = _factor_map(discriminant)
    if set(conductor_factor) != set(discriminant_factor):
        raise AssertionError("conductor and minimal discriminant bad-prime supports differ")

    multiplicative_primes = sorted(
        prime for prime, exponent in conductor_factor.items() if exponent == 1
    )
    additive_primes = sorted(
        prime for prime, exponent in conductor_factor.items() if exponent > 1
    )
    phi_m = sorted(
        prime
        for prime in multiplicative_primes
        if discriminant_factor[prime] % 2 == 0
    )

    additive_records = []
    additive_term = 0
    for prime in additive_primes:
        decomposition = prime_decomp(
            prime, cubic, ZK=maximal_order, dK=field_discriminant
        )
        count = len(decomposition)
        additive_term += count - 1
        additive_records.append(
            {
                "prime": prime,
                "minimal_discriminant_valuation": discriminant_factor[prime],
                "conductor_exponent": conductor_factor[prime],
                "number_of_cubic_primes": count,
                "contribution_np_minus_one": count - 1,
                "prime_decomposition": _serialize_prime_decomposition(decomposition),
            }
        )

    u_term = 2 if discriminant > 0 else 1
    n_term = len(phi_m) + additive_term
    known_rank_lower_bound = 20
    rational_two_torsion_dimension = 0
    class_group_two_rank_lower_bound = (
        known_rank_lower_bound
        + rational_two_torsion_dimension
        - u_term
        - n_term
    )
    if class_group_two_rank_lower_bound < 0:
        raise AssertionError("the derived class-group lower bound is nonsensical")

    result: dict[str, Any] = {
        "schema_version": 1,
        "claim_level": (
            "unconditional exact lower bound on the 2-rank of the ideal class "
            "group of the R20 2-division cubic field; no class-group upper bound "
            "and no new elliptic-curve rank claim"
        ),
        "source": {
            "candidate": "fermigier-mestre-v1:u=28917/20",
            "path": str(R20_NEAR_MISS.relative_to(ROOT)),
            "sha256": R20_NEAR_MISS_SHA256,
            "rank_lower_bound": known_rank_lower_bound,
            "minimal_model": list(ainvs),
            "minimal_discriminant": str(discriminant),
            "conductor": str(conductor),
        },
        "two_division_field": {
            "primitive_two_division_cubic_coefficients_ascending": list(
                primitive_two_division
            ),
            "coordinate_change": (
                f"z={primitive_two_division[-1]}*x, where x is a nonzero "
                "2-torsion abscissa"
            ),
            "monic_defining_polynomial_coefficients_ascending": list(
                monic_coefficients
            ),
            "irreducibility_witness": irreducibility,
            "rational_two_torsion_dimension": rational_two_torsion_dimension,
            "polynomial_discriminant": str(polynomial_discriminant),
            "polynomial_discriminant_factorization": _factor_record(
                polynomial_discriminant
            ),
            "field_discriminant": str(field_discriminant),
            "field_discriminant_factorization": _factor_record(field_discriminant),
            "power_order_index": str(field_index),
            "power_order_index_factorization": _factor_record(field_index),
            "signature": {
                "real_places": int(cubic.count_roots(-sp.oo, sp.oo)),
                "complex_place_pairs": 0,
            },
            "maximal_order_basis": {
                "denominator": int(maximal_order.denom),
                "basis_matrix_columns": [
                    [int(value) for value in row]
                    for row in maximal_order.matrix.to_Matrix().tolist()
                ],
                "basis_interpretation": (
                    "columns are numerator coefficient vectors in the power "
                    "basis (1,z,z^2), divided by the common denominator"
                ),
            },
        },
        "brumer_kramer": {
            "primary_source": PRIMARY_SOURCE,
            "hypothesis_rational_two_torsion_trivial": True,
            "multiplicative_primes": multiplicative_primes,
            "additive_primes": additive_primes,
            "phi_m_even_discriminant_valuation": phi_m,
            "phi_m_count": len(phi_m),
            "additive_data": additive_records,
            "additive_term": additive_term,
            "u_term": u_term,
            "n_term": n_term,
            "inequality": "dim Sel_2(E/Q) <= dim Cl(K)[2] + u(E) + n(E)",
            "lower_chain": (
                "20 <= rank(E(Q)) <= dim Sel_2(E/Q) "
                f"<= dim Cl(K)[2] + {u_term} + {n_term}"
            ),
            "class_group_two_rank_lower_bound": class_group_two_rank_lower_bound,
            "conclusion": (
                "dim_F2 Cl(K)[2] >= "
                f"{class_group_two_rank_lower_bound} unconditionally"
            ),
        },
        "next_gate": {
            "target_class_group_two_rank": class_group_two_rank_lower_bound,
            "exact_closure_condition": (
                "an unconditional upper certificate dim_F2 Cl(K)[2] <= 13 "
                "would force dim Sel_2(E/Q)=rank(E(Q))=20 and close this fixed fibre"
            ),
            "conditional_closure_condition": (
                "a GRH-conditional class-group upper bound 13 would close the "
                "fixed fibre subject to that GRH assumption"
            ),
            "residual_cover_condition": (
                "if the class-group 2-rank exceeds 13, compute the relative "
                "2-Selmer quotient and explicit covers not explained by the "
                "twenty pinned generators"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": REPRODUCING_COMMAND,
    }
    result["result_sha256"] = canonical_digest(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = build_record()
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing pinned output: {args.output}")
        pinned = json.loads(args.output.read_text(encoding="utf-8"))
        if pinned != record:
            raise SystemExit("R20_BRUMER_KRAMER_PINNED_OUTPUT_MISMATCH")
        print("R20_BRUMER_KRAMER_CHECK_PASS")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print("R20_BRUMER_KRAMER_WRITE_PASS")
    print(record["brumer_kramer"]["conclusion"])


if __name__ == "__main__":
    main()

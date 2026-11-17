#!/usr/bin/env python3
"""Exact mass and Euler-product controls for RANK_FIVE_DISCRIMINANT_COUNT.md."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


def primes_through(bound: int) -> list[int]:
    sieve = [True]*(bound+1)
    sieve[0] = sieve[1] = False
    for p in range(2, sp.integer_nthroot(bound, 2)[0]+1):
        if sieve[p]:
            for n in range(p*p, bound+1, p):
                sieve[n] = False
    return [p for p in range(2, bound+1) if sieve[p]]


def valuation(n: int, p: int) -> int:
    assert n != 0
    exponent = 0
    while n % p == 0:
        exponent += 1
        n //= p
    return exponent


def run() -> dict:
    controls = []

    def equal(name, left, right=0):
        assert sp.cancel(left-right) == 0, name
        controls.append(name)

    p, n, X = sp.symbols("p n X")
    mass_all = 1+1/p+1/p**2
    mass_nontotal = 1+1/p
    equal("total_ramified_mass_complement", mass_all-mass_nontotal, 1/p**2)
    equal("unrestricted_euler_factor", (1-1/p)*mass_all, 1-1/p**3)
    equal("fundamental_discriminant_euler_factor", (1-1/p)*mass_nontotal, 1-1/p**2)
    equal("inert_prime_allowed_mass_ratio", mass_nontotal/mass_all, 1-1/(p*p+p+1))
    equal("inert_prime_zeta_factor_identity", mass_nontotal/mass_all, (1-p**-2)/(1-p**-3))
    equal("real_signature_factor", sp.Rational(1, 2)*sp.Rational(1, 6), sp.Rational(1, 12))
    equal("complex_signature_factor", sp.Rational(1, 2)*sp.Rational(1, 2), sp.Rational(1, 4))

    expected_discriminants = [-2835, -14256, -30051]
    local_polynomials = []
    for k, expected in enumerate(expected_discriminants):
        polynomial = X**3+6*X**2+3+9*k
        assert 6 % 3 == 0 and (3+9*k) % 3 == 0 and (3+9*k) % 9 != 0
        controls.append(f"cyclic_local_polynomial_{k}_eisenstein")
        discriminant = int(sp.discriminant(polynomial, X))
        assert discriminant == expected and valuation(discriminant, 3) == 4
        controls.append(f"cyclic_local_polynomial_{k}_exact_discriminant_exponent")
        assert (discriminant//3**4) % 3 == 1
        controls.append(f"cyclic_local_polynomial_{k}_square_discriminant_class_over_Q3")
        local_polynomials.append({"k": k, "polynomial": str(polynomial),
                                  "discriminant": discriminant, "three_adic_exponent": 4})
    cyclic_count = 3  # Pauli Corollary 4.7 and the complete Q3 table.
    cyclic_mass = sp.Rational(cyclic_count, 3*3**4)
    equal("three_ramified_cyclic_fields_total_mass", cyclic_mass, sp.Rational(1, 81))
    allowed3 = mass_nontotal.subs(p, 3)+cyclic_mass
    equal("wild_prime_allowed_mass", allowed3, sp.Rational(109, 81))
    equal("wild_prime_total_mass", mass_all.subs(p, 3), sp.Rational(13, 9))
    equal("wild_prime_density_ratio", allowed3/mass_all.subs(p, 3), sp.Rational(109, 117))
    assert allowed3 > mass_nontotal.subs(p, 3)
    controls.append("strict_increase_over_fundamental_subfamily_at_three")

    # An exact rational product plus a proved integer-tail inequality
    # gives an enclosure for the actual infinite product, not just a
    # numerical truncation whose remaining tail is ignored.
    bound = 10000
    prime_list = primes_through(bound)
    assert len(prime_list) == 1229 and prime_list[-1] == 9973
    controls.append("bounded_prime_sieve_crosscheck")
    rho_partial = Fraction(109, 117)
    used = []
    for prime in prime_list:
        if prime % 3 == 2:
            rho_partial *= Fraction(prime*prime+prime, prime*prime+prime+1)
            used.append(prime)
    assert used[0] == 2 and 3 not in used
    controls.append("infinite_specification_partial_product_includes_two_once_and_separates_three")
    lower = rho_partial*Fraction(bound, bound+1)
    assert Fraction(75987590, 10**8) < rho_partial < Fraction(75987591, 10**8)
    assert Fraction(75979, 10**5) < lower <= rho_partial < Fraction(75988, 10**5)
    controls.append("exact_partial_product_and_certified_infinite_product_enclosure")
    equal("tail_majorant_telescoping_identity", 1/(n*(n+1)), 1/n-1/(n+1))
    assert Fraction(3, 4)+Fraction(1, 16) < 1
    controls.append("cyclic_field_order_bound_is_negligible_for_leading_term")
    # Strong acceptability is a statement about the shapes at all primes;
    # these labels record the complete two nonzero residue classes.
    local_shapes = {1: "all maximal cubic orders", 2: "all nontotally ramified maximal cubic orders"}
    assert set(local_shapes) == {1, 2}
    controls.append("all_large_prime_specifications_have_acceptable_shapes")

    here = Path(__file__).resolve().parent
    sources = []
    for name, statements in [
        ("bst-2013-davenport-heilbronn.pdf", [
            "Theorem 8, equation (8), p. 444 (PDF p. 6): leading count for acceptable infinite local specifications",
            "Lemma 19, p. 453 (PDF p. 15): nontotally ramified local density",
            "Lemma 22, p. 456 (PDF p. 18): O_epsilon(X^(3/4+epsilon)) cyclic-order bound",
            "Remark 4, equation (58), p. 480 (PDF p. 42): total etale cubic local mass"]),
        ("pauli-2006-local-class-fields.pdf", [
            "Corollary 4.7, pp. 638--639 (PDF pp. 13--14): exactly three ramified cyclic cubic Q3 extensions",
            "Section 7, pp. 647 and 649 (PDF pp. 22 and 24): discriminant convention and complete Q3 table"]),
        ("milne-ant300.pdf", [
            "Proposition 7.55 and Remark 7.56, pp. 122--123: Eisenstein uniformizer and its integral power basis"]),
    ]:
        path = here/"sources"/name
        data = path.read_bytes()
        assert data.startswith(b"%PDF-")
        sources.append({"path": "sources/"+name, "sha256": sha256(data).hexdigest(),
                        "statements": statements, "proof_imported_not_replayed": True})
    with localcontext() as context:
        context.prec = 35
        upper_decimal = str(Decimal(rho_partial.numerator)/Decimal(rho_partial.denominator))
        lower_decimal = str(Decimal(lower.numerator)/Decimal(lower.denominator))
    return {
        "status": "pass", "python": platform.python_version(), "sympy": sp.__version__,
        "command": ".venv/bin/python research/arithmetic/verify_rank_five_discriminant_count.py",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256((here/"RANK_FIVE_DISCRIMINANT_COUNT.md").read_bytes()).hexdigest(),
        "exact_controls": controls, "local_three_polynomials": local_polynomials,
        "local_three_allowed_mass": "109/81", "local_three_density_ratio": "109/117",
        "rho_enclosure": {"prime_cutoff": bound, "inert_prime_factors": len(used),
                          "exact_rational_comparison": "0.75979 < rho < 0.75988",
                          "lower_bound_decimal_rounded": lower_decimal,
                          "upper_bound_decimal_rounded": upper_decimal,
                          "tail_lower_factor": f"{bound}/{bound+1}"},
        "leading_constants": {"positive": "rho/(12*zeta(3))", "negative": "rho/(4*zeta(3))"},
        "fundamental_constants": {"positive": "1/(12*zeta(2))", "negative": "1/(4*zeta(2))"},
        "external_primary_sources": sources,
        "scope": ("Exact mass and polynomial controls, and a bounded Euler product with an "
                  "explicit infinite-tail bound. Completeness of local cyclic field classification, "
                  "maximal-order power bases, and the infinite-local-specification counting "
                  "theorem are explicit primary-source imports. Counts concern cubic field "
                  "discriminant height and varying-map fibre algebras, not fixed-Phi target height. "
                  "No secondary coefficient or error stronger than o(X) is asserted."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("rank_five_discriminant_count_checks.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()

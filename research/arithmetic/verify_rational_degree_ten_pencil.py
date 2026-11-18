#!/usr/bin/env python3
"""Exact degree-ten rational paired-critical-value pencil certificate.

Run from the repository root:
    .venv/bin/python research/arithmetic/verify_rational_degree_ten_pencil.py
Optional --write-receipt PATH. Requires pinned sympy==1.14.0.
No numerical root continuation or finite parameter sampling is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform

import sympy as sp


def modular_power(base: sp.Poly, exponent: int, modulus: sp.Poly) -> sp.Poly:
    result = sp.Poly(1, base.gens[0], modulus=base.get_modulus())
    while exponent:
        if exponent & 1:
            result = (result * base).rem(modulus)
        base = (base * base).rem(modulus)
        exponent >>= 1
    return result


def rabin_irreducible(polynomial: sp.Poly) -> bool:
    """Rabin's exact criterion; independent of Poly.is_irreducible."""
    polynomial = polynomial.monic()
    prime = polynomial.get_modulus()
    degree = polynomial.degree()
    variable = polynomial.gens[0]
    x = sp.Poly(variable, variable, modulus=prime)
    for divisor in sp.factorint(degree):
        residual = modular_power(x, prime ** (degree // divisor), polynomial) - x
        if sp.gcd(polynomial, residual).degree() != 0:
            return False
    return (modular_power(x, prime ** degree, polynomial) - x).rem(polynomial).is_zero


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    if sp.__version__ != "1.14.0":
        raise RuntimeError("Replay requires pinned sympy==1.14.0")

    z, s, u = sp.symbols("z s u")
    h = (
        9*z**10 - 864*z**8 + 28161*z**6
        - 350150*z**4 + 1127500*z**2
    )
    a = sp.Integer(250000)
    branch_quartic = (
        59049*s**4 - 126938174857200000*s**2
        + 68130645548641000000000000000
    )
    square_root = 2187*s*branch_quartic
    f = h-s*z+a
    delta = sp.discriminant(f, z)
    assert sp.expand(delta-square_root**2) == 0
    assert sp.Poly(square_root, s).degree() == 5
    assert sp.Poly(square_root, s).LC() == 3**17
    assert sp.Poly(sp.gcd(square_root, sp.diff(square_root, s)),s).degree() == 0
    assert sp.discriminant(square_root, s) != 0

    critical_quartic = 81*u**4 - 4023*u**3 + 40230*u**2 - 44700*u + 10000
    critical_u = (u-25)*critical_quartic
    critical_z = sp.expand(z*sp.diff(h,z)-h-a)
    assert sp.expand(critical_z-critical_u.subs(u,z**2)) == 0
    assert sp.discriminant(critical_quartic,u) == 1801825579550588085000000
    assert critical_quartic.subs(u,25) == -7182500
    assert critical_quartic.subs(u,0) == 10000
    critical_test_points = [0, sp.Rational(1,2), 1, 25, 50]
    assert [sp.sign(critical_quartic.subs(u,x)) for x in critical_test_points] == [1,-1,1,-1,1]
    assert sp.discriminant(critical_u,u) == 92953129806939475496581031250000000000
    assert sp.discriminant(critical_z,z) != 0
    assert sp.Poly(critical_z,z).degree() == 10
    assert critical_z.subs(z,0) == -a
    assert sp.expand(h+a-(z**2-25)**2*(9*z**6-414*z**4+1836*z**2+400)) == 0

    branch_a = sp.Integer(59049)
    branch_b = -sp.Integer(126938174857200000)
    branch_c = sp.Integer(68130645548641000000000000000)
    pair_discriminant = branch_b**2-4*branch_a*branch_c
    assert pair_discriminant == 21114280070272604391840000000000
    assert branch_c != 0
    assert branch_a > 0 and branch_b < 0 and branch_c > 0 and pair_discriminant > 0

    # A second presentation, used to find the model by reciprocal pairing
    # of the four remaining critical u-coordinates after fixing u=1.
    h_monic = (
        z**10-sp.Rational(96,25)*z**8+sp.Rational(3129,625)*z**6
        -sp.Rational(14006,5625)*z**4+sp.Rational(1804,5625)*z**2
    )
    a_monic = sp.Rational(16,5625)
    assert sp.expand(h-9*5**10*h_monic.subs(z,z/5)) == 0
    assert a == 9*5**10*a_monic

    # Independent finite-field group certificate.
    # The leading coefficient is 9, so normalize to a monic polynomial.
    prime = 13
    cubic = sp.Poly(z**3+2*z**2+3*z+1,z,modulus=prime)
    septic = sp.Poly(
        z**7-2*z**6-4*z**5-3*z**3-3*z**2+4*z+4,
        z,modulus=prime,
    )
    specialized = sp.Poly(f.subs(s,12),z,modulus=prime).monic()
    assert specialized == cubic*septic
    assert rabin_irreducible(cubic)
    assert rabin_irreducible(septic)
    assert sp.gcd(specialized,specialized.diff()).degree() == 0

    # Wrong branch polynomial must fail the exact identity.
    assert sp.expand(delta-(square_root+1)**2) != 0
    # A square discriminant with index-three ramification is an invalid
    # answer to this paired-simple-critical-point problem.
    wrong_h = z**10+sp.Rational(3,2)*z**5
    wrong_a = sp.Integer(-1)
    wrong_r = z*sp.diff(wrong_h,z)-wrong_h-wrong_a
    assert sp.expand(wrong_r-(3*z**5+1)**2) == 0
    assert sp.discriminant(wrong_r,z) == 0
    assert not rabin_irreducible(sp.Poly(z**7-z,z,modulus=prime))

    report = {
        "schema_version": 1,
        "kind": "internal_exact_replay",
        "command": ".venv/bin/python research/arithmetic/verify_rational_degree_ten_pencil.py",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "arithmetic": "Z[z,s], Q[z,s] and F_13",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input": {
            "H": str(h), "a": str(a),
            "discriminant_square_root": str(sp.expand(square_root)),
            "factored_square_root": str(square_root),
            "critical_u": str(sp.expand(critical_u)),
        },
        "checks": {
            "full_discriminant_identity": True,
            "degree_discriminant_square_root": 5,
            "square_root_squarefree": True,
            "critical_degree": 10,
            "ten_distinct_nonzero_critical_points": True,
            "ten_real_critical_points": True,
            "five_real_branch_values": True,
            "critical_u_discriminant": str(sp.discriminant(critical_u,u)),
            "critical_quartic_discriminant": str(sp.discriminant(critical_quartic,u)),
            "critical_quartic_at_25": "-7182500",
            "branch_pair_discriminant": str(pair_discriminant),
            "rational_normalization_identity": True,
            "mod_13_s": 12,
            "mod_13_factor_degrees": [3,7],
            "mod_13_cubic": str(cubic.as_expr()),
            "mod_13_septic": str(septic.as_expr()),
            "mod_13_irreducibility_method": "Rabin modular-power/gcd criterion",
            "mod_13_squarefree": True,
            "negative_controls": 3,
        },
        "conclusions_with_written_proof": {
            "finite_branch_values": 5,
            "finite_branch_cycle_type": "two disjoint transpositions",
            "infinity_profile": [9,1],
            "geometric_monodromy": "A10",
            "arithmetic_monodromy": "A10",
            "algebraic_constant_field": "Q",
            "collinear_ordinary_nodes": 5,
        },
        "scope": (
            "The exact identities prove the literal rational square discriminant, "
            "simple critical points and distinct paired values. Full monodromy "
            "also uses irreducibility, the nine-cycle at infinity, the Frobenius "
            "(3,7) certificate, Jordan's prime-cycle theorem and simplicity of A10, "
            "as detailed in the accompanying note. No external human review, "
            "formal verification or historical priority is asserted."
        ),
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
        "result": "PASS",
    }
    rendered = json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(rendered,encoding="utf-8")
    print(rendered,end="")


if __name__ == "__main__":
    main()

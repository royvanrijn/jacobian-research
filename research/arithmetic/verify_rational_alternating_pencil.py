#!/usr/bin/env python3
"""Exact rational degree-eight alternating-pencil family certificate.

Run with the repository's pinned environment:
    .venv/bin/python research/arithmetic/verify_rational_alternating_pencil.py

Requires sympy==1.14.0. Identities are polynomial identities over Q[c,t,s].
The finite-field certificate is exact over F_11. The companion note supplies
the monodromy, descent, and infinitely-many-isomorphism-classes arguments.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import sympy as sp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    if sp.__version__ != "1.14.0":
        raise RuntimeError("Replay requires the pinned sympy==1.14.0")
    z, u, s, c, t = sp.symbols("z u s c t")
    d = (47*c**2 + 7*t**2) / 168
    e = (379*c**3 + 357*c*t**2 + 196*t**3) / 21168
    a = -(7*t**2 - 3*c**2)**2 / 87808
    h = z**8 + c*z**6 + d*z**4 + e*z**2
    b = (
        6912*c**7 - 14112*c**5*t**2 + 5488*c**4*t**3
        + 8575*c**3*t**4 - 7203*c**2*t**5
        + 31213*c*t**6 + 16807*t**7
    ) / 533655864
    q_constant = (
        (7*t**2-3*c**2)**2 * (4*c**2-7*c*t+7*t**2)**2
        * (48*c**2+56*c*t+35*t**2)**3 / 379718108242115328
    )
    q = s**4 + b*s**2 + q_constant
    critical = 7*u**4 + 5*c*u**3 + 3*d*u**2 + e*u - a
    expected_critical_disc = (
        t**2*(t-2*c)**6*(c+t)**2*(48*c**2+56*c*t+35*t**2)
        / 1259712
    )
    expected_branch_pair_disc = (
        -t**6*(t-2*c)**5*(c+t)**3 / 10806531246
    )
    delta = sp.discriminant(h-s*z+a, z)
    assert sp.expand(delta + 7**7*q**2) == 0
    assert sp.factor(sp.discriminant(critical, u)-expected_critical_disc) == 0
    assert sp.factor(b**2-4*q_constant-expected_branch_pair_disc) == 0

    # The positive-definite quadratic factors cannot vanish at a nonzero real
    # pair (c,t). These discriminants check the elementary real reduction used
    # in the note; the rational equation 7t^2=3c^2 has only the zero solution.
    assert sp.discriminant(4*c**2-7*c*t+7*t**2, t) == -63*c**2
    assert sp.discriminant(48*c**2+56*c*t+35*t**2, t) == -3584*c**2
    assert sp.sqrt(sp.Rational(3, 7)).is_rational is False

    seed = {c: 1, t: 1}
    seed_h = sp.expand(h.subs(seed))
    seed_a = a.subs(seed)
    seed_q = q.subs(seed)
    assert seed_h == z**8+z**6+sp.Rational(9,28)*z**4+sp.Rational(233,5292)*z**2
    assert seed_a == -sp.Rational(1,5488)
    assert sp.expand(seed_q - (
        s**4+sp.Rational(5960,66706983)*s**2
        +sp.Rational(2685619,1483273860320763)
    )) == 0
    assert expected_critical_disc.subs(seed) == sp.Rational(139,314928)
    assert sp.discriminant(seed_q,s) != 0

    prime = 11

    def modular_polynomial(expression: sp.Expr) -> sp.Poly:
        coefficients = {}
        for monomial, coefficient in sp.Poly(expression, z, domain=sp.QQ).terms():
            numerator, denominator = coefficient.as_numer_denom()
            coefficients[monomial] = (
                int(numerator) * pow(int(denominator), -1, prime) % prime
            )
        return sp.Poly.from_dict(coefficients, z, modulus=prime)

    cubic = sp.Poly(z**3+z**2-3*z-4, z, modulus=prime)
    quintic = sp.Poly(z**5-z**4+5*z**3-4*z**2-3, z, modulus=prime)
    specialized = modular_polynomial(seed_h-2*z+seed_a)
    assert specialized == sp.Poly(
        z**8+z**6+7*z**4+2*z**2-2*z+1, z, modulus=prime
    )
    assert cubic*quintic == specialized
    assert cubic.is_irreducible and quintic.is_irreducible
    assert sp.gcd(specialized, specialized.diff()).degree() == 0

    # Reducing t=1+11k (and c=1) fixes the entire reduction for every integer k.
    k = sp.symbols("k", integer=True)
    one_parameter = sp.Poly(
        sp.expand((h-s*z+a).subs({c:1,t:1+11*k,s:2})), z,k,domain=sp.QQ
    )
    for monomial, coefficient in one_parameter.terms():
        numerator, denominator = coefficient.as_numer_denom()
        assert int(denominator) % prime != 0
        if monomial[1] > 0:
            assert int(numerator) % prime == 0

    # Invalid inputs: changing the critical value certificate is detected;
    # t=2c degenerates the branch and critical polynomials; a reducible
    # polynomial cannot serve as the irreducible cubic Frobenius factor.
    assert sp.expand(delta+7**7*(q+1)**2) != 0
    assert expected_critical_disc.subs(t,2*c) == 0
    assert expected_branch_pair_disc.subs(t,2*c) == 0
    assert not sp.Poly(z**3-z,z,modulus=prime).is_irreducible

    report = {
        "schema_version": 1,
        "kind": "internal_exact_replay",
        "command": ".venv/bin/python research/arithmetic/verify_rational_alternating_pencil.py",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "arithmetic": "Q[c,t,s], Q and F_11",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameter_family": {
            "H": str(h), "a": str(a), "Q": str(q),
            "rational_nondegenerate_domain": "t != 0, t != 2c, t != -c",
            "infinite_A8_subfamily": "c=1, t=1+11k, integer k>=0",
        },
        "seed": {"c": 1, "t": 1, "H": str(seed_h), "a": str(seed_a), "Q": str(seed_q)},
        "checks": {
            "family_discriminant_identity": True,
            "critical_discriminant_factorization": str(expected_critical_disc),
            "branch_pair_discriminant_factorization": str(expected_branch_pair_disc),
            "seed_critical_discriminant": "139/314928",
            "seed_branch_squarefree": True,
            "mod_11_s": 2,
            "mod_11_cubic": str(cubic.as_expr()),
            "mod_11_quintic": str(quintic.as_expr()),
            "mod_11_factor_degrees": [3,5],
            "mod_11_factors_irreducible": True,
            "all_t_equal_1_mod_11_have_same_reduction": True,
            "negative_controls": 3,
        },
        "scope": (
            "Exact all-parameter polynomial identities, nondegeneracy factors "
            "and a finite-field Frobenius certificate. The companion written "
            "proof uses Jordan's theorem and simplicity of A8 for monodromy, "
            "and a coefficient invariant for infinitely many cover classes. "
            "No human external review or formal verification is asserted."
        ),
        "result": "PASS",
    }
    output = json.dumps(report, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(output, encoding="utf-8")
    print(output,end="")


if __name__ == "__main__":
    main()

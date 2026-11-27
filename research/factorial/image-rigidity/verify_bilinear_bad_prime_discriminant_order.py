#!/usr/bin/env python3
"""Exact single bad-lift obstruction; not a search for factorial counterexamples."""
import hashlib
import json
import platform
import time
from pathlib import Path

import sympy as sp


def main():
    started = time.monotonic()
    u, a = sp.symbols("u A")
    m, p, n = 5, 7, 34
    # Construct E_n=U_n(u/2) by its monic integer recurrence.
    previous, current = sp.Integer(1), u
    for _ in range(1, n):
        previous, current = current, sp.expand(u * current - previous)
    f = sum(sp.diff(current, u, j) * a**j for j in range(n + 1))
    print("Computing the exact degree-34 factorial-transform discriminant.", flush=True)
    delta = sp.Poly(sp.discriminant(f, a), u, domain=sp.ZZ)
    content, primitive = delta.primitive()
    integer, valuation = abs(int(content)), 0
    while integer % p == 0:
        integer //= p
        valuation += 1
    reduced = sp.Poly(primitive, u, modulus=p)
    factor = sp.Poly(u**4 + 2*u*u - 2, u, modulus=p)
    remaining, order = reduced, 0
    while remaining.rem(factor).is_zero:
        remaining = remaining.exquo(factor)
        order += 1
    assert delta.degree() == n * (n - 1) == 1122
    assert valuation == 136 and order == 49
    assert reduced.eval(0) != 0
    assert not remaining.rem(factor).is_zero
    # A separate F_7 polynomial identity verifies the reported exact order.
    assert factor**order * remaining == reduced
    assert sp.Poly(u**5-5*u**3+5*u, u, modulus=p) == sp.Poly(u, u, modulus=p)*factor
    receipt = {
        "status": "passed",
        "evidence_type": "internal_exact_single_polynomial_computation",
        "scope": "Only m=5,p=7,n=34; refutes the low-order proof bound, not a factorial window.",
        "command": ".venv/bin/python research/factorial/image-rigidity/verify_bilinear_bad_prime_discriminant_order.py",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "parameters": {"m": m, "p": p, "n": n},
        "input_E34_ascending": list(reversed([int(c) for c in sp.Poly(current, u).all_coeffs()])),
        "discriminant_degree": delta.degree(),
        "discriminant_integer_content_v7": valuation,
        "primitive_discriminant_mod7_ascending": list(reversed([int(c) % p for c in reduced.all_coeffs()])),
        "noncentral_D5_factor_mod7": "u^4+2*u^2-2",
        "factor_order": order,
        "quotient_mod7_ascending": list(reversed([int(c) % p for c in remaining.all_coeffs()])),
        "next_factor_remainder_mod7": [int(c) % p for c in remaining.rem(factor).all_coeffs()],
        "order_less_than_p_minus_one": order < p - 1,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic()-started, 3),
        "independent_replay": False, "external_review": False, "formal_verification": False,
    }
    output = Path(__file__).with_name("bilinear_bad_prime_discriminant_order_verification.json")
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "passed", "content_v7": valuation,
                      "factor_order": order, "receipt": str(output)}), flush=True)


if __name__ == "__main__":
    main()

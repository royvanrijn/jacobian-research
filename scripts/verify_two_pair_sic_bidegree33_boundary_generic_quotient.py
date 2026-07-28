#!/usr/bin/env python3
"""Modular generic rank-six quotient on the reduced s0 common boundary.

This checker works over two rational-function fields of positive
characteristic.  It is a structural reconstruction aid, not by itself a
characteristic-zero unit certificate.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp

from explore_two_pair_sic_bidegree33_full_anchor import (
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_boundary_generic_quotient.json"
)
PRIMES = (47, 101)


def certificate(prime: int) -> dict[str, object]:
    singular = shutil.which("Singular")
    assert singular is not None
    expressions = [
        chart_expression(moment_terms(order, prime), 0, prime)
        for order in range(2, 8)
    ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        120,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring r=({prime},s1,s2,s3,t0,t1,t2),(s5,t4),dp;
poly p4={polynomials[1]};
poly p5={polynomials[2]};
poly p6={polynomials[3]};
poly p7={polynomials[4]};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);
print("GBSIZE "+string(size(G)));
print("VDIM "+string(vdim(G)));
print("R6SIZE "+string(size(r6)));
print("R7SIZE "+string(size(r7)));
print("R6 "+string(r6));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    summaries = [
        re.search(rf"(?m)^{label} (\d+)$", completed.stdout)
        for label in ("GBSIZE", "VDIM", "R6SIZE", "R7SIZE")
    ]
    assert all(marker is not None for marker in summaries), completed.stdout[:1000]
    assert tuple(int(marker.group(1)) for marker in summaries if marker) == (3, 6, 6, 6)

    denominators: list[str] = []
    for denominator in re.findall(r"/\(([^()]*)\)\*", completed.stdout):
        if denominator not in denominators:
            denominators.append(denominator)
    assert len(denominators) == 2

    s1, s2, s3, t0, t1, t2 = sp.symbols("s1 s2 s3 t0 t1 t2")
    environment = {
        str(symbol): symbol for symbol in (s1, s2, s3, t0, t1, t2)
    }
    parsed = [
        sp.Poly(
            sp.sympify(value.replace("^", "**"), locals=environment),
            s1,
            s2,
            s3,
            t0,
            t1,
            t2,
            modulus=prime,
        )
        for value in denominators
    ]
    linear = s1 * t0 - t1
    quadratic = s1**2 - s2 - (13 * pow(3, -1, prime) % prime) * t0**2
    expected = [
        sp.Poly(linear * quadratic, *parsed[0].gens, modulus=prime),
        sp.Poly(linear * quadratic**2, *parsed[0].gens, modulus=prime),
    ]
    assert parsed == expected
    return {
        "prime": prime,
        "groebner_basis_size": 3,
        "quotient_length": 6,
        "mu6_remainder_terms": 6,
        "mu7_remainder_terms": 6,
        "denominator_factors": [
            "s1*t0-t1",
            "s1^2-s2-(13/3)*t0^2",
        ],
        "denominator_products": ["L*Q", "L*Q^2"],
    }


def main() -> None:
    payload = {
        "certificates": [certificate(prime) for prime in PRIMES],
        "scope": (
            "exact rational-function-field calculations in two finite "
            "characteristics; reconstruction evidence, not a QQ unit certificate"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("PASS common boundary: generic (mu_4,mu_5) quotient has length six")
    print("PASS mu_6 and mu_7 reduce to six fiber coordinates")
    print("PASS denominator divisor reconstructs as L*Q with Q coefficient 13/3")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify one explicit Hasse failure using only standard-library arithmetic.

General Keller reconstruction and the S3 local-cover criterion retain their
written proof scopes. This verifies every explicit coefficient and witness.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path,
                        default=Path(__file__).with_name("elliptic_exception_hasse_certificate.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cert = json.loads(args.certificate.read_text(encoding="utf-8"))
    checks: list[dict[str, object]] = []

    def require(name: str, condition: bool, **details) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append({"name": name, "passed": True, **details})

    def multiply(first, second):
        out = [F(0)] * (len(first)+len(second)-1)
        for i, aa in enumerate(first):
            for j, bb in enumerate(second):
                out[i+j] += aa*bb
        return out

    def evaluate(coefficients, value):
        result = 0
        for coefficient in reversed(coefficients):
            result = result*value + coefficient
        return result

    def discriminant_cubic(coefficients):
        dd, cc, bb, aa = coefficients
        return (bb*bb*cc*cc - 4*aa*cc**3 - 4*bb**3*dd
                - 27*aa*aa*dd*dd + 18*aa*bb*cc*dd)

    def valuation(value: int, prime: int) -> int:
        require("nonzero_witness_value", value != 0)
        value = abs(value)
        exponent = 0
        while value % prime == 0:
            exponent += 1
            value //= prime
        return exponent

    def is_prime(value: int) -> bool:
        if value < 2:
            return False
        if value % 2 == 0:
            return value == 2
        return all(value % divisor for divisor in range(3, math.isqrt(value)+1, 2))

    pi, ell, vv = (F(cert[key]) for key in ("pi", "ell", "v"))
    target = list(map(F, cert["target"]))
    d = int(cert["quadratic_radicand"])
    rq, rh = F(cert["quadratic_discriminant_multiplier"]), F(cert["cubic_discriminant_multiplier"])
    rr = -F(3,2)-ell
    ss = F(3,2)/pi**2-ell*rr-vv
    tt = (F(9,16)/pi-vv*ss)/ell
    quadratic = [vv,ell,F(1)]
    cubic = [tt,ss,rr,F(1)]
    require("original_cubic_coefficients", cubic == list(map(F,cert["cubic_coefficients_low_to_high"])))
    expected_target = [-pi,-F(32,9)*pi**2*(tt+ell*ss+vv*rr+F(5,4)/pi),-F(32,9)*vv*tt]
    require("exact_target_coordinates", target == expected_target)
    bb, cc = target[1:]
    inverse = [-F(9,32)*cc,F(9,16)/pi,-F(5,4)/pi-F(9,32)*bb/pi**2,
               F(3,2)/pi**2,-F(3,2),F(1)]
    require("whole_inverse_factorization", multiply(quadratic,cubic) == inverse)
    dq = ell**2-4*vv
    dh = discriminant_cubic(cubic)
    require("quadratic_discriminant_squareclass", dq == d*rq**2)
    require("cubic_discriminant_squareclass", dh == d*rh**2)
    require("quadratic_irreducible_over_Q", dq < 0)
    require("nonzero_factor_discriminants", dq != 0 and dh != 0)
    alpha = ell**2-vv-rr*ell+ss
    beta = (ell-rr)*vv+tt
    resultant = beta**2-ell*alpha*beta+vv*alpha**2
    require("factors_coprime", resultant != 0)
    different_class = F(d,-3)
    require("quadratic_field_differs_from_minus_three",
            math.isqrt(different_class.numerator)**2 != different_class.numerator
            or math.isqrt(different_class.denominator)**2 != different_class.denominator)

    scale = math.lcm(*(coefficient.denominator for coefficient in cubic))
    require("integer_scale_is_exact_LCM", scale == int(cert["integer_scale"]))
    integral = [coefficient*scale**(3-index) for index,coefficient in enumerate(cubic)]
    require("scaled_cubic_is_monic_integral", integral[-1] == 1 and all(value.denominator == 1 for value in integral))
    integral = [value.numerator for value in integral]
    require("exact_integral_cubic", integral == list(map(int,cert["integral_cubic_coefficients_low_to_high"])))
    dint = discriminant_cubic(integral)
    require("discriminant_scaling", dint == scale**6*dh)
    factorization = cert["integral_cubic_discriminant_factorization"]
    primes = {prime for prime,exponent in factorization}
    require("distinct_discriminant_factors", len(primes) == len(factorization))
    require("all_discriminant_factors_prime", all(is_prime(prime) for prime in primes))
    claimed_discriminant = cert["integral_cubic_discriminant_sign"]*math.prod(prime**exponent for prime,exponent in factorization)
    require("entire_discriminant_factorization", dint == claimed_discriminant)
    p_irred = int(cert["irreducible_mod_prime"])
    require("irreducibility_prime_is_prime", is_prime(p_irred))
    require("irreducible_cubic_mod_prime", integral[-1] % p_irred != 0
            and all(evaluate(integral,residue) % p_irred for residue in range(p_irred)),
            prime=p_irred, residues_checked=p_irred)

    q_primes = set()
    for witness in cert["quadratic_local_witnesses"]:
        prime, root = witness["prime"], int(witness["root"])
        require(f"quadratic_Hensel_prime_{prime}", prime in primes and (root**2-d) % prime == 0
                and (2*root) % prime != 0)
        q_primes.add(prime)
    derivative = [index*coefficient for index,coefficient in enumerate(integral)][1:]
    h_primes = set()
    for witness in cert["cubic_local_witnesses"]:
        prime, root = witness["prime"], int(witness["root"])
        val_poly = valuation(evaluate(integral,root),prime)
        val_derivative = valuation(evaluate(derivative,root),prime)
        require(f"cubic_strict_Hensel_prime_{prime}", prime in primes
                and val_poly == witness["polynomial_valuation"]
                and val_derivative == witness["derivative_valuation"]
                and val_poly > 2*val_derivative,
                polynomial_valuation=val_poly,derivative_valuation=val_derivative)
        h_primes.add(prime)
    require("all_discriminant_primes_have_local_witnesses", q_primes | h_primes == primes)
    require("witness_categories_disjoint", not (q_primes & h_primes))

    # The finite S3 element-cover property used for unramified places.
    for permutation in itertools.permutations(range(3)):
        odd = sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3)) % 2
        require(f"S3_unramified_cover_{permutation}", not odd or any(permutation[i] == i for i in range(3)))

    cover_y = d*rq*rh
    cover_z = cover_y/(819200*vv-243483)
    j_value = (400*vv-81)*(20480000*vv**2-110592000*vv+143467857)
    require("elliptic_normalization_exact_point", cover_z**2 == j_value/41231686041600000000)
    xx = F(6400,3)*vv-3984
    yy = F(20971520000,3)*cover_z
    require("short_Weierstrass_exact_point", yy**2 == xx**3-10758382*xx+6600759744)
    require("outside_two_counted_lines", target[0] not in (F(-1),F(2)))
    require("not_an_integral_target", target[0].denominator != 1)

    source = Path(__file__).with_name("ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md")
    report = {
        "status":"passed",
        "scope":"Standard-library exact verification of one explicit target, integral cubic, irreducibility, discriminant factors and all required local witnesses; general reconstruction and local-cover theorems remain written inputs.",
        "python":platform.python_version(),
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "certificate_sha256":hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
        "check_count":len(checks),"checks":checks,
        "discriminant_prime_count":len(primes),
        "elliptic_point":{"X":str(xx),"T":str(yy)},
    }
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"check_count":len(checks),
                      "discriminant_primes":len(primes),"output":str(args.output) if args.output else None}))


if __name__ == "__main__":
    main()

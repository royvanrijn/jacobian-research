#!/usr/bin/env python3
"""Exact rational certificates for elementary constants in the Euler analogue.

These calculations verify the finite ingredients of elementary proofs.
For the divisor constants, p>=32 contributes a factor <=1 because e+1<=2^e.
For each smaller relevant prime, both the corrected tau factor and the M
factor decrease after the first exponent E with (E+2)^5<=p(E+1)^5.
The corrected tau factor has an additional next-term ratio strictly below 1.
Consequently the enumerated exponent maxima are global maxima, not guesses.
"""

from fractions import Fraction
import hashlib
import json
from pathlib import Path


def certify(primes, kind, bound):
    product = Fraction(1)
    rows = []
    for p in primes:
        end = 0
        while (end + 2)**5 > p * (end + 1)**5:
            end += 1
        values = []
        for e in range(end + 1):
            if kind == "tau":
                factor = Fraction((e + 1) * (p - 1) * p**e, p**(e + 1) - 1)
                fifth = factor**5 / p**e
            else:
                fifth = Fraction((e + 1)**5, p**e)
            values.append(fifth)
        maximum_e = max(range(len(values)), key=values.__getitem__)
        maximum_fifth = values[maximum_e]
        product *= maximum_fifth
        rows.append({"prime": p, "last_exponent_required": end,
                     "maximizing_exponent": maximum_e,
                     "maximum_fifth_power": str(maximum_fifth)})
    assert product < bound**5
    return {"local_maxima": rows, "product_fifth_power": str(product),
            "upper_bound": str(bound), "upper_bound_fifth_power": str(bound**5),
            "strict_inequality_verified": True}


def main():
    tau = certify([2, 5, 11, 17, 23, 29], "tau", Fraction(5, 2))
    maximum = certify([7, 13, 19, 31], "M", Fraction(15, 8))
    assert 100 > Fraction(5, 2)**5
    union = Fraction(2, 5) + Fraction(17, 12) * Fraction(2, 5) + Fraction(22, 3) * Fraction(1, 250)
    assert union == Fraction(249, 250) < 1
    receipt = {
        "tau_scope": "n has only prime factors congruent to 2 mod3",
        "tau_bound": "tau(n)*n/sigma(n) < (5/2)*n^(1/5)",
        "tau_constant": tau,
        "M_scope": "M(n)=max_{1<=u<=n} product_{p^e||u,p=1mod3}(e+1)",
        "M_bound": "M(n) < (15/8)*n^(1/5)",
        "M_constant": maximum,
        "large_prime_control": "For p>=32, (e+1)^5<=32^e<=p^e, by e+1<=2^e for e>=0",
        "exponent_tail_control": "Next-term ratios for (e+1)/p^(e/5) decrease with e. The corrected tau ratio has the additional factor p*(p^(e+1)-1)/(p^(e+2)-1)<1.",
        "threshold": 100000000,
        "threshold_variable": "X=N^(1/4)>=100",
        "fifth_root_control": "100>(5/2)^5 implies X^(-1/5)<2/5; X^(-6/5)<1/250",
        "assumed_union_bound": "2/5+(17/12)X^(-1/5)+(22/3)X^(-6/5)",
        "resulting_strict_upper_bound": str(union),
        "lattice_bound_proof": "Represent each six-unit orbit by x>=1,0<=y<x. Distinct labels L=x(x-1)/2+y+1 satisfy Q-L=(y-(x+1)/2)^2+(x^2-5)/4>=0 for x>=3, with x=1,2 checked. Therefore N_E^*(t)<=6floor(t)<=6t.",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    path = Path(__file__).with_name("elementary_constants_receipt.json")
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

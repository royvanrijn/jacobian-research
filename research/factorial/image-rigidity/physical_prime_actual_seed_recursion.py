"""Exact controls of the finite Artin–Hasse recursion."""
from fractions import Fraction
from math import factorial
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ORDER = 30


def artin_hasse_coefficients(prime, cutoff):
    powers = []
    power = 1
    while power <= cutoff:
        powers.append(power)
        power *= prime
    coefficients = [Fraction(1)]
    for degree in range(1, cutoff + 1):
        coefficients.append(sum(coefficients[degree-power]
                                for power in powers if power <= degree) / degree)
    return coefficients


def residue(value, modulus):
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def multiply(left, right, modulus):
    return [
        sum(left[j] * right[i-j] for j in range(i+1)) % modulus
        for i in range(ORDER+1)
    ]


def direct_coefficient(degree, prime, value):
    quotient, row = divmod(degree, prime-1)
    return row, Fraction(factorial(row)*(prime*value)**quotient, factorial(degree))



count = 0
one = [1] + [0]*ORDER
for prime in [2, 3, 5, 7, 11]:
    q = prime - 1
    ah = artin_hasse_coefficients(prime, ORDER)
    for precision in range(1, 5):
        modulus = prime**precision
        for value in (0, 1, 2, -1):
            expected = [[0]*q for _ in range(ORDER+1)]
            for degree in range(ORDER+1):
                row, coefficient = direct_coefficient(degree, prime, value)
                expected[degree][row] = residue(coefficient, modulus)

            chi = [0]*(ORDER+1)
            power, i = prime, 1
            while power <= ORDER:
                E = (power-1)//q
                chi[power] = -prime**(E-i)*value**E % modulus
                power *= prime
                i += 1

            nested = [[0]*q for _ in range(ORDER+1)]
            chi_power = one
            for degree in range(ORDER+1):
                row, coefficient = direct_coefficient(degree, prime, value)
                coefficient = residue(coefficient, modulus)
                for k, value_at_k in enumerate(chi_power):
                    nested[k][row] = (nested[k][row] + coefficient*value_at_k) % modulus
                chi_power = multiply(chi_power, chi, modulus)

            actual = [[0]*q for _ in range(ORDER+1)]
            for degree in range(ORDER+1):
                coefficient = residue(ah[degree], modulus)
                for row in range(q):
                    wraps, target = divmod(degree+row, q)
                    entry = Fraction((prime*value)**wraps * factorial(target), factorial(row))
                    scale = coefficient * residue(entry, modulus) % modulus
                    if degree >= q*precision:
                        assert scale == 0
                    for k in range(ORDER+1-degree):
                        actual[k+degree][target] += scale*nested[k][row]
            for i, (direct, computed) in enumerate(zip(expected, actual)):
                assert direct == [c % modulus for c in computed], (prime, precision, value, i)
                count += q
source = HERE/"PHYSICAL_PRIME_ACTUAL_SEED_RECURSION.tex"
files = [
    Path(__file__), source, HERE/"PHYSICAL_PRIME_RAISING_CALCULUS.tex",
    HERE/"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
]
receipt = {
    "status": "PASS",
    "controls": count,
    "source_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
    "scope": "Direct cyclic coefficients checked against the finite Artin–Hasse recursion.",
    "independent_whole_proof_replay": False,
    "formal_verification": False,
    "external_review": False
}
(HERE/"physical_prime_actual_seed_recursion_verification.json").write_text(
    json.dumps(receipt,indent=2) + "\n"
)
print(json.dumps({"status": "PASS", "controls": count}))

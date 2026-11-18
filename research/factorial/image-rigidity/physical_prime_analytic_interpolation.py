"""Exact controls of the residue formulas and the nonpositive continuation."""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform

from flint import nmod_poly

HERE = Path(__file__).resolve().parent
COUNTS = {}


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    COUNTS[name] = COUNTS.get(name, 0) + 1


def residue_pair(prime, n):
    h = [0]
    w = []
    for k in range(n + 1):
        top = n + (prime - 1) * k - 1
        w.append((-1)**k * comb(top, prime * k) % prime)
        if k:
            h.append((-1)**k * comb(top, prime * k - 1) % prime)
    return nmod_poly(h, prime), nmod_poly(w, prime)


def substitute_reciprocal(poly, base, exponent, prime):
    """Evaluate poly(base^-1) and clear base^exponent."""
    result = nmod_poly([], prime)
    for k, coefficient in enumerate(poly):
        if coefficient:
            result += int(coefficient) * base**(exponent - k)
    return result


for prime in [2, 3, 5, 7, 11]:
    S = nmod_poly([0, 1], prime)
    B = S**prime - S
    quotients = [B // (S + offset) for offset in range(prime)]
    base = B**(prime - 1)
    for n in range(1, 18):
        h, w = residue_pair(prime, n)
        check("actual_W_reciprocal_trace",
              substitute_reciprocal(w, base, n - 1, prime)
              == -sum((qi**(prime*n-prime+1) for qi in quotients), nmod_poly([], prime)))
        check("actual_H_reciprocal_trace",
              substitute_reciprocal(h, base, n, prime)
              == sum((qi**n for qi in quotients), nmod_poly([], prime))**prime)
        check("W_residue_multiplicity", all(power <= prime - 1 for _, power in w.factor()[1]))
        check("H_residue_multiplicity",
              all(power % prime == 0 for factor, power in h.factor()[1] if factor[0]))
        paired = h.gcd(w)
        check("common_residue_multiplicity",
              all(power <= prime - 2 for _, power in paired.factor()[1]))

_, w145 = residue_pair(7, 145)
h145, _ = residue_pair(7, 145)
gcd145 = h145.gcd(w145)
check("nontrivial_common_residue_control", gcd145.degree() == 32)
check("common_residue_bound_at_145", all(power == 4 for _, power in gcd145.factor()[1]))


def weights(prime, limit, modulus):
    result = [1]
    for k in range(limit):
        value = result[-1]
        for i in range(1, prime):
            value = value * (prime*k + i) % modulus
        result.append(value)
    return result


for prime, precision in [(2, 1), (2, 4), (3, 3), (5, 2), (7, 2), (11, 1)]:
    modulus = prime**precision
    degree = modulus
    length = prime**(precision - 1)
    multiplier = 1 if prime == 2 and precision >= 3 else -1
    upper = 25
    factorial_units = weights(prime, max(degree, upper), modulus)
    for Y in [1, -1, 2, 4]:
        denominator = [(-1)**j * comb(degree, j) % modulus for j in range(degree + 1)]
        denominator[length] = (denominator[length] - multiplier * pow(Y, length, modulus)) % modulus
        inverse = pow(denominator[degree], -1, modulus)
        sequence = {
            n: sum(factorial_units[k] * comb(n+(prime-1)*k-1, prime*k)
                   * pow(Y, k, modulus) for k in range(n)) % modulus
            for n in range(1, degree + 1)
        }
        for n in range(0, -upper-1, -1):
            actual = -inverse * sum(denominator[j] * sequence[n+degree-j]
                                    for j in range(degree)) % modulus
            sequence[n] = actual
            u = -n
            expected = sum((-1)**(prime*k) * factorial_units[k] * comb(u+k, prime*k)
                           * pow(Y, k, modulus) for k in range(u//(prime-1) + 1)) % modulus
            check("two_sided_recurrence_continuation", actual == expected)


def valuation(n, prime):
    result = 0
    while n and n % prime == 0:
        n //= prime
        result += 1
    return result


for prime in [2, 3, 7]:
    for s in [0, 1, 2]:
        for b in [-4, 0, 3, 19]:
            product = [1]
            denominator = 1
            for k in range(1, 30):
                updated = [0] * (len(product) + 1)
                for i, coefficient in enumerate(product):
                    updated[i] += coefficient*(b-k+1)
                    updated[i+1] += coefficient*prime**s
                product = updated
                denominator *= k
                tail = 0
                j = prime**(s+1)
                while j <= k:
                    tail += k//j
                    j *= prime
                check("falling_factorial_Gauss_bound",
                      all(coefficient == 0
                          or valuation(Fraction(coefficient, denominator).numerator, prime)
                          - valuation(Fraction(coefficient, denominator).denominator, prime) >= -tail
                          for coefficient in product))


def inverse_in_quadratic(value, prime, polynomial):
    residue_polynomial = nmod_poly([int(c) for c in polynomial], prime)
    reduced = nmod_poly([int(c) for c in value], prime)
    gcd, inverse, _ = reduced.xgcd(residue_polynomial)
    assert gcd == nmod_poly([1], prime)
    lifted = nmod_poly([int(c) for c in inverse], prime**2)
    return (lifted * (2 - value*lifted)) % polynomial


def modular_evaluation(coefficients, argument, polynomial, modulus):
    value = nmod_poly([], modulus)
    for coefficient in reversed(coefficients):
        value = (value*argument + int(coefficient)) % polynomial
    return value


for prime in [3, 5, 7, 11]:
    modulus = prime**2
    w = factorial(prime-1)
    nonsquare = next(d for d in range(1, prime) if pow(d, (prime-1)//2, prime) == prime-1)
    polynomial = nmod_poly([-nonsquare, 0, 1], modulus)
    generator = nmod_poly([0, 1], modulus)
    for lift in [generator, generator + 1, generator + prime]:
        B = (pow(lift, prime, polynomial) - lift) % polynomial
        inverse = inverse_in_quadratic(B, prime, polynomial)
        Y = pow(inverse, prime-1, polynomial)
        omega = [pow(lift+i, prime**2, polynomial) for i in range(prime)]
        exponents = {(-ell*pow(prime, -1, prime**2-1)) % (prime**2-1): ell
                     for ell in range(1, prime-1)}
        for n in sorted(set(range(1,14)) | exponents.keys()):
            unit_coefficients = weights(prime, n, modulus)
            coefficients = [prime]
            for k in range(1, n+1):
                top = n+(prime-1)*k-1
                coefficients.append(unit_coefficients[k] *
                                    (prime*comb(top,prime*k)+comb(top,prime*k-1)) % modulus)
            actual = modular_evaluation(coefficients, Y, polynomial, modulus)
            tame_exponent = -prime*n % (prime**2-1)
            predicted = pow(-w*B, n, polynomial)
            predicted *= sum((pow(root,tame_exponent,polynomial) for root in omega),
                             nmod_poly([], modulus))
            check("unramified_trace_mod_p2", actual == predicted % polynomial)
            if n in exponents:
                check("common_residue_cluster_trace_valuation_one",
                      actual and all(int(c) % prime == 0 for c in actual))
                reduced_Y = nmod_poly([int(c) for c in Y], prime)
                reduced_polynomial = nmod_poly([int(c) for c in polynomial], prime)
                h, lower = residue_pair(prime, n)
                check("prescribed_H_residue", not modular_evaluation(h, reduced_Y, reduced_polynomial, prime))
                ell = exponents[n]
                for j in range(ell+1):
                    row = [int(lower[k])*comb(k,j) for k in range(j, len(lower))]
                    check("prescribed_exact_W_multiplicity",
                          bool(modular_evaluation(row, reduced_Y, reduced_polynomial, prime)) == (j == ell))

    for ell in range(1, prime-1):
        degree = prime*ell
        left = nmod_poly([comb(degree,j)*sum(i**(degree-j) for i in range(prime))
                          for j in range(degree+1)], modulus)
        expected = prime * (generator**prime-generator)**ell
        check("symbolic_affine_line_power_sum", left == expected)


def prime_divisors(value):
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            yield divisor
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        yield value


def certify_gap_six_initial_interval():
    prime, upper = 7, 529
    Y = nmod_poly([0, 1], prime)
    precisions = {n: prime**(3 + valuation(n, prime)) for n in range(2, upper+1)}
    unit_coefficients = weights(prime, upper, max(precisions.values()))
    zero_remainders = []
    for n, modulus in precisions.items():
        w_coefficients, h_coefficients = [], []
        for k in range(n+1):
            top = n + (prime-1)*k - 1
            unit = unit_coefficients[k]
            lower = unit * comb(top, prime*k)
            upper_term = unit * comb(top, prime*k-1) if k else 0
            w_coefficients.append(lower % modulus)
            h_coefficients.append((prime*lower + upper_term) % modulus)
        w = nmod_poly(w_coefficients, prime)
        h = nmod_poly(h_coefficients, prime)
        gcd = w.gcd(h)
        lifted_h = nmod_poly(h_coefficients, modulus)
        product = nmod_poly([1], prime)
        for factor, exponent in gcd.factor()[1]:
            degree = factor.degree()
            check("gap_six_irreducible_residue_factor",
                  all((pow(Y, prime**(degree//d), factor)-Y).gcd(factor).degree() == 0
                      for d in prime_divisors(degree))
                  and pow(Y, prime**degree, factor) == Y % factor)
            product *= factor**exponent
            lift = nmod_poly([int(c) for c in factor], modulus)
            remainder = lifted_h % lift
            if not remainder:
                zero_remainders.append((n, factor))
                check("gap_six_eisenstein_residue",
                      n == 529 and factor == Y-2 and n % prime != 0
                      and exponent == 3
                      and not w % factor**3 and bool(w % factor**4)
                      and not h % factor**14 and bool(h % factor**15)
                      and nmod_poly(w_coefficients, prime**2)(2) == prime)
            else:
                check("gap_six_unramified_point_value_exclusion", bool(remainder))
        check("gap_six_complete_residue_factorization", product == gcd)
    check("gap_six_unique_zero_point_value_remainder",
          zero_remainders == [(529, Y-2)])
    return {"terminal_min": 2, "terminal_max": upper,
            "eisenstein_exception_to_point_value_test": {"n": 529, "residue": 2}}


gap_six_interval = certify_gap_six_initial_interval()

files = [
    Path(__file__),
    HERE / "PHYSICAL_PRIME_ANALYTIC_INTERPOLATION.tex",
    HERE / "PHYSICAL_PRIME_MOMENT_PAIR_RECURRENCES.tex",
    HERE / "PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex",
    HERE / "GAP_SIX_TWO_DIGIT_LEADING_TWO_SFC.tex",
]
record = {
    "status": "PASS",
    "controls": sum(COUNTS.values()),
    "control_counts": COUNTS,
    "scope": "Exact controls of the written identities and complete replay of the finite gap-six interval. The analytic interpolation and resultant theorem have written proofs.",
    "gap_six_interval": gap_six_interval,
    "source_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
    "python": platform.python_version(),
    "independent_whole_proof_replay": False,
    "formal_verification": False,
    "external_review": False,
}
(HERE / "physical_prime_analytic_interpolation_verification.json").write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps({"status": "PASS", "controls": record["controls"], "counts": COUNTS}))

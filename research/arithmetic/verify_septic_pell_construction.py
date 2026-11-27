"""Exact finite inputs for SEPTIC_PELL_FIVE_DIRECTION_CONSTRUCTION.md.

Python 3.12, SymPy 1.14, PARI/GP 2.15.4. Default: complete F_5 sieve
and a one-direction characteristic-zero control. --system K emits the
five-row coefficient system with K quadratic and 5-K quintic Z rows.
No five-direction candidate or specialization theorem is certified.
"""

import argparse
from collections import Counter
import hashlib
import itertools
import json
import platform
import subprocess

import sympy as s

from verify_lyw_polynomial_pell import gp_frobenius


t = s.Symbol("t")
EXPECTED_COUNTS = {"bad": 5704, "0": 5689, "1": 2670, "2": 1019,
                   "3": 396, "4": 110, "5": 33, "6": 4}
EXPECTED_TRACE = "49837a583da277362346f09e77b0388112b813f55557254ec03c489d4cfabceb"


def valuation_three(integer):
    assert integer > 0
    exponent = 0
    while integer % 3 == 0:
        integer //= 3
        exponent += 1
    return exponent


def root_one_multiplicity(coefficients):
    polynomial = s.Poly.from_list(list(reversed(coefficients)), t, modulus=3)
    factor = s.Poly(t-1, t, modulus=3)
    multiplicity = 0
    while polynomial.eval(1) == 0:
        polynomial = polynomial.exquo(factor)
        multiplicity += 1
    return multiplicity


def validate_frobenius(coefficients, f, prime):
    """Necessary controls; these alone do not derive the characteristic polynomial."""
    assert len(coefficients) == 13 and coefficients[-1] == 1
    assert coefficients[0] == prime**6
    for i in range(6):
        assert coefficients[i] == prime**(6-i)*coefficients[12-i]
    points = 2
    for value in range(prime):
        residue = int(f.eval(value)) % prime
        points += 1 if residue == 0 else (
            2 if pow(residue, (prime-1)//2, prime) == 1 else 0)
    assert points == prime+1+coefficients[-2]


def finite_field_sieve():
    """All monic degree-seven F_5 polynomials with zero t^6 coefficient.

    Reflection identifies P(t) with -P(-t). Multiplicities restore the
    entire 5^6-element normalized domain. See the note for affine coverage.
    """
    counts = Counter()
    large_multiplicities = Counter()
    trace = hashlib.sha256()
    representatives = 0
    good_representatives = 0
    maximum_bound = 0
    for coefficients in itertools.product(range(5), repeat=6):
        reflected = tuple((-c if i % 2 == 0 else c) % 5
                          for i, c in enumerate(coefficients))
        if coefficients > reflected:
            continue
        representatives += 1
        weight = 1 if coefficients == reflected else 2
        p = s.Poly(t**7+sum(c*t**i for i, c in enumerate(coefficients)),
                   t, modulus=5)
        f = p*p+1
        if s.gcd(f, f.diff()).degree():
            counts["bad"] += weight
            continue
        good_representatives += 1
        frobenius = gp_frobenius([int(f.nth(i)) for i in range(15)], 5)
        validate_frobenius(frobenius, f, 5)
        order = sum(frobenius)
        assert order % 7 == 0
        exponent = valuation_three(order)
        counts[str(exponent)] += weight
        # Either integer is an upper bound for dim ker(Frob-1) on J[3].
        multiplicity = root_one_multiplicity(frobenius) if exponent >= 5 else 12
        upper_bound = min(exponent, multiplicity)
        maximum_bound = max(maximum_bound, upper_bound)
        assert upper_bound <= 4, (coefficients, frobenius)
        if exponent >= 5:
            large_multiplicities[str(multiplicity)] += weight
        trace.update(json.dumps([coefficients, frobenius],
                                separators=(",", ":")).encode())
    assert representatives == 7875
    assert sum(counts.values()) == 5**6
    assert dict(counts) == EXPECTED_COUNTS
    assert dict(large_multiplicities) == {"1": 20, "2": 8, "3": 9}
    assert maximum_bound == 4
    assert trace.hexdigest() == EXPECTED_TRACE
    return {"normalized_polynomials": 5**6, "reflection_representatives": representatives,
            "good_representatives": good_representatives,
            "weighted_counts_by_v3_order": dict(sorted(counts.items())),
            "multiplicity_at_one_when_v3_at_least_five": dict(sorted(large_multiplicities.items())),
            "maximum_certified_upper_bound": maximum_bound,
            "trace_sha256": trace.hexdigest()}


def one_direction_control():
    z = 2*t*t+2*t
    linear = 2*t+8
    h = s.expand(z*linear-1)
    p = s.expand((z*z*linear**3-3*z*linear**2+3*linear-z)/2)
    assert s.expand(2*p*z+z*z-1-h**3) == 0
    assert s.expand(2*p*z+z*z-1-(h+1)**3) != 0
    assert p == (16*t**7+224*t**6+1168*t**5+2740*t**4
                 +2708*t**3+735*t*t-190*t+12)
    f = s.Poly(p*p+1, t)
    assert s.gcd(f, f.diff()).degree() == 0
    reduced_five = s.Poly(f, t, modulus=5)
    assert s.gcd(reduced_five, reduced_five.diff()).degree() > 0
    reduced_seven = s.Poly(f, t, modulus=7)
    assert reduced_seven.degree() == 14
    assert s.gcd(reduced_seven, reduced_seven.diff()).degree() == 0
    frobenius = gp_frobenius([int(f.nth(i)) for i in range(15)], 7)
    assert frobenius == [117649, 0, 0, 1372, -833, 476, 392, 68, -17, 4, 0, 0, 1]
    validate_frobenius(frobenius, reduced_seven, 7)
    assert sum(frobenius) == 119112 and valuation_three(sum(frobenius)) == 1
    for value in (0, 9, 18):
        assert int(p.subs(t, value)) % 9 == 3
        assert int(z.subs(t, value)) % 9 == 0
        assert int(h.subs(t, value)) % 9 == 8
    # The cyclotomic quadratic branch also gives a polynomial P.
    z_other = t*t-t+1
    h_other = (t+1)*z_other+t
    quotient, remainder = s.div(h_other**3+1, z_other, t)
    assert remainder == 0
    p_other = s.expand((quotient-z_other)/2)
    assert s.degree(p_other, t) == 7
    assert s.expand(2*p_other*z_other+z_other**2-1-h_other**3) == 0
    # A corrupt middle coefficient must fail the reciprocity control.
    corrupted = frobenius.copy()
    corrupted[4] += 1
    try:
        validate_frobenius(corrupted, reduced_seven, 7)
    except AssertionError:
        pass
    else:
        raise AssertionError("corrupt Frobenius passed")
    return {"P_ascending": [int(s.Poly(p, t).nth(i)) for i in range(8)],
            "Z_ascending": [0, 2, 2], "H_ascending": [-1, 16, 20, 4],
            "bad_at_five": True, "good_seven_order": sum(frobenius),
            "rational_three_rank": 1,
            "unit_and_local_cube_slice": "n=0 mod9; virtuality and specialized independence not asserted"}


def coefficient_system(quadratic_rows):
    assert 0 <= quadratic_rows <= 5
    variables = [f"p{j}" for j in range(8)]
    p = sum(s.Symbol(f"p{j}")*t**j for j in range(8))
    equations = []
    nonzero = ["p7", "disc_t(P^2+1)"]
    for row in range(5):
        z_degree, h_degree = (2, 3) if row < quadratic_rows else (5, 4)
        z_names = [f"z{row}_{j}" for j in range(z_degree+1)]
        h_names = [f"h{row}_{j}" for j in range(h_degree+1)]
        variables.extend(z_names+h_names)
        z = sum(s.Symbol(name)*t**j for j, name in enumerate(z_names))
        h = sum(s.Symbol(name)*t**j for j, name in enumerate(h_names))
        residual = s.Poly(2*p*z+z*z-1-h**3, t)
        equations.extend(str(residual.nth(j)) for j in range(8+z_degree))
        nonzero.extend([z_names[-1], h_names[-1]])
    assert len(variables) == 63-4*quadratic_rows
    assert len(equations) == 65-3*quadratic_rows
    return {"quadratic_rows": quadratic_rows, "quintic_rows": 5-quadratic_rows,
            "coefficient_domain": "Q", "variables": variables,
            "equations_equal_zero": equations, "required_nonzero": nonzero,
            "separate_obligations": ["five independent divided divisor classes",
                                     "integral unit and virtuality on a common slice",
                                     "local separation at three", "distinct-field infinitude"],
            "five_integral_degree_seven_chart": "must have disc_t(P^2+1)=0 mod5"}


def verify():
    control = one_direction_control()
    for rows in range(6):
        coefficient_system(rows)
    sieve = finite_field_sieve()
    version = subprocess.run(["gp", "--version-short"], capture_output=True,
                             text=True, check=True).stdout.strip()
    return {"python": platform.python_version(), "sympy": s.__version__,
            "pari_gp": version, "control": control, "complete_mod_five_sieve": sieve,
            "scope": "degree profiles, complete good-reduction obstruction at five, one-direction control; rank-five target open"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", type=int, choices=range(6), metavar="K")
    arguments = parser.parse_args()
    result = verify() if arguments.system is None else coefficient_system(arguments.system)
    print(json.dumps(result, indent=2))

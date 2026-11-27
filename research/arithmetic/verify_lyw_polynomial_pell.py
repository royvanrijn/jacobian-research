"""Exact inputs for LY8--LY9. Python 3.12, SymPy 1.14, PARI/GP 2.15.4.

This is a finite torsion obstruction, not an integer Pell prefix search.
The written reduction, continued-fraction completeness and pushforward
arguments are in LEVIN_YAN_WILJANEN_DEGREE_TEN_MIRROR_SLICE.md.
"""

import hashlib
import json
import math
import platform
import subprocess

import sympy as s

from verify_lyw_degree_ten_mirror_slice import F_COEFFS, poly, t


# Ascending coefficients of det(X - Frobenius) on H^1.
FROBENIUS = {
    59: [12117361, 0, 163607, 35046, 528, 594, 47, 0, 1],
    71: [25411681, -4294932, 660371, -41535, 5847, -585, 131, -12, 1],
}


def gp_frobenius(coefficients, prime):
    expression = "+".join(f"({c})*x^{i}" for i, c in enumerate(coefficients))
    program = (
        f"F=Mod(1,{prime})*({expression});"
        'P=hyperellcharpoly(F);print(Vec(Vecrev(Vec(P))));\n'
    )
    run = subprocess.run(
        ["gp", "-q", "-f"], input=program, text=True,
        capture_output=True, check=True, timeout=30,
    )
    assert not run.stderr.strip(), run.stderr
    return json.loads(run.stdout.strip())


def pell_prefix(coefficients, prime, bound):
    """All polynomial CF convergents of numerator degree <= bound.

    A constant next denominator detects a polynomial unit of any constant
    nonzero norm, including nonsquare norms. No normalization to norm one
    over the finite field is assumed.
    """
    f = s.Poly(poly(coefficients), t, modulus=prime)
    assert f.degree() % 2 == 0 and s.gcd(f, f.diff()).degree() == 0
    half_degree = f.degree() // 2
    lead = next(k for k in range(1, prime) if k*k % prime == int(f.LC()) % prime)
    a0 = s.Poly(lead*t**half_degree, t, modulus=prime)
    for j in range(half_degree - 1, -1, -1):
        coefficient = (f - a0*a0).nth(half_degree + j)
        coefficient = int(coefficient)*pow(2*lead, -1, prime) % prime
        a0 += s.Poly(coefficient*t**j, t, modulus=prime)
    assert (f - a0*a0).degree() < half_degree
    p = s.Poly(0, t, modulus=prime)
    q = s.Poly(1, t, modulus=prime)
    a = a0
    degree = 0
    trace = hashlib.sha256()
    for step in range(bound + 1):
        assert a.degree() > 0
        degree += a.degree()
        p = a*q - p
        q = (f - p*p).exquo(q)
        assert not q.is_zero
        trace.update(json.dumps([
            degree,
            [int(c) % prime for c in p.all_coeffs()],
            [int(c) % prime for c in q.all_coeffs()],
        ], separators=(",", ":")).encode())
        if degree > bound:
            return {"unit_degree": None, "steps": step + 1,
                    "next_degree": degree, "next_denominator_degree": q.degree(),
                    "trace_sha256": trace.hexdigest()}
        if q.degree() == 0:
            return {"unit_degree": degree, "steps": step + 1}
        a = (a0 + p).quo(q)
    raise AssertionError("positive partial-quotient degrees did not exceed bound")


def verify():
    f = poly(F_COEFFS)
    orders = []
    for prime, expected in FROBENIUS.items():
        reduced = s.Poly(f, t, modulus=prime)
        assert reduced.degree() == 10
        assert s.gcd(reduced, reduced.diff()).degree() == 0
        actual = gp_frobenius(F_COEFFS, prime)
        assert actual == expected, (prime, actual, expected)
        # Reciprocal coordinates give the same smooth projective curve.
        assert gp_frobenius(list(reversed(F_COEFFS)), prime) == expected
        # Direct F_p point count independently checks the first trace.
        points = 2  # Both rational infinities: the leading coefficient is a square.
        for value in range(prime):
            residue = int(reduced.eval(value)) % prime
            points += 1 if residue == 0 else (2 if pow(residue, (prime - 1)//2, prime) == 1 else 0)
        assert points == prime + 1 + expected[-2]
        orders.append(sum(actual))
    bound = math.gcd(*orders)
    assert orders == [12317184, 21740967] and bound == 2187
    assert orders[0] % 71 and orders[1] % 59
    obstruction = pell_prefix(F_COEFFS, 59, bound)
    assert obstruction["unit_degree"] is None
    assert (obstruction["steps"], obstruction["next_degree"],
            obstruction["next_denominator_degree"]) == (2147, 2188, 4)

    # Positive controls: primitive degree two and a degree divisible by three.
    for power in (2, 3):
        control = s.Poly(t**(2*power) + 1, t)
        coefficients = [int(control.nth(i)) for i in range(control.degree() + 1)]
        assert pell_prefix(coefficients, 59, power)["unit_degree"] == power
        assert pell_prefix(coefficients, 59, power - 1)["unit_degree"] is None

    # LY9: universal norm identities; every relevant mod-nine residue.
    p, z = s.symbols("p z")
    assert s.expand((p + z)**2 - (p*p + 1) - (2*p*z + z*z - 1)) == 0
    residues = []
    for value in range(9):
        roots = [r for r in range(9) if (r*r - value*value - 1) % 9 == 0]
        if value % 3 == 0:
            assert len(roots) == 2
            cubes = [(value + r) % 9 in (1, 8) for r in roots]
            assert cubes == [value == 0]*2
            residues.append([value, roots, cubes])
        else:
            assert roots == []
    # If Z=0 mod9, each ratio (y-P-Z)/(y-P) is 1 mod9,
    # simultaneously with the noncube norm-minus-one unit at P=3,6.
    for value in (3, 6):
        for root in (1, 8):
            denominator = (root - value) % 9
            assert math.gcd(denominator, 9) == 1
            assert (root - value)*pow(denominator, -1, 9) % 9 == 1
    candidate = s.Poly((t**7 + 3)**2 + 1, t)
    assert s.gcd(candidate, candidate.diff()).degree() == 0
    candidate_coeffs = [int(candidate.nth(i)) for i in range(15)]
    candidate_frob = gp_frobenius(candidate_coeffs, 11)
    assert candidate_frob == [1771561, 0, 0, 0, 0, 0, -394, 0, 0, 0, 0, 0, 1]
    assert sum(candidate_frob) == 1771168 and sum(candidate_frob) % 3 != 0

    # This altered claimed certificate must fail exact comparison.
    corrupted = FROBENIUS[59].copy()
    corrupted[4] += 1
    assert gp_frobenius(F_COEFFS, 59) != corrupted
    version = subprocess.run(["gp", "--version-short"], text=True,
                             capture_output=True, check=True).stdout.strip()
    return {"python": platform.python_version(), "sympy": s.__version__,
            "pari_gp": version, "jacobian_orders": orders, "torsion_bound": bound,
            "mod59_polynomial_cf": obstruction, "norm_minus_one_residues": residues,
            "candidate_t7_plus_3_order_mod11": sum(candidate_frob),
            "scope": "exact inputs for LY8--LY9; no real rank-five family"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))

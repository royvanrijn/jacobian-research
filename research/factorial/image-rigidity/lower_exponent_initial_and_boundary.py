"""Exact internal replay for lower_exponent_initial_and_boundary.tex.

The universal statements are proved in that note. This bounded symbolic replay
checks the stated identities and boundary-space rank; it does not extrapolate.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import platform
import sympy as sp

HERE = Path(__file__).resolve().parent
s, z, t = sp.symbols("s z t")

def rising(x, n):
    out = 1
    for j in range(n):
        out *= x + j
    return out

def U(n, d, shape_minus_one, zz):
    return sp.expand(sum((-1)**k * comb(n, k)
                         * rising(shape_minus_one + 1, d*k) * zz**(n-k)
                         for k in range(n + 1)))

def equation(v, n, d, zz):
    return sp.expand(zz*(v - v.subs(s, s+1))
        - rising(s+2, d-1)*((s+d*n+1)*v.subs(s, s+d)
                           - (s+d+1)*v.subs(s, s+d+1)))

def replay():
    cases = 0
    for n, d in product(range(2, 7), range(1, 5)):
        v = U(n-1, d, s, z)
        m = d*(n-1)
        assert equation(v, n, d, z) == 0
        assert sp.Poly(v, s).degree() == m
        gap = sp.Poly(v - (-1)**(n-1)*rising(s+1, m), s)
        assert gap.is_zero or gap.degree() <= m-d
        # Independent coefficient-matrix nullspace check at nonzero and zero z.
        # The proof's degree argument, not these samples, gives all complex z.
        for zz in [0, 1, -2, sp.I]:
            columns = [sp.Poly(equation(s**j, n, d, zz), s)
                       for j in range(m+1)]
            max_degree = max(p.degree() for p in columns if not p.is_zero)
            mat = sp.Matrix([[p.nth(i) for p in columns]
                             for i in range(max_degree+1)])
            assert mat.rank() == m
            residual = mat * sp.Matrix([sp.Poly(v.subs(z, zz), s).nth(j)
                                        for j in range(m+1)])
            assert all(sp.expand(entry) == 0 for entry in residual)
        cases += 1

    inequality_cases = 0
    for d, x, delta in product(range(1, 9),
                               [Fraction(1, 7), Fraction(1, 2), 1, 3, 20],
                               [0, Fraction(1, 3), 1, 4, 17]):
        y = x + delta
        assert rising(y, 3*d) > rising(x, 2*d)*(3*rising(y, d)-2*rising(x, d))
        inequality_cases += 1

    moment_cases = 0
    for a, d, beta in product(range(0, 6), range(1, 7),
                              [sp.Rational(1, 3), sp.Rational(1, 2), 1, 3]):
        f2 = sum(comb(2, k)*rising(beta+2*a, d*k)*t**k for k in range(3))
        f3 = sum(comb(3, k)*rising(beta+3*a, d*k)*t**k for k in range(4))
        assert sp.Poly(sp.gcd(f2, f3), t).degree() == 0
        # The finite obstruction equations are direct consequences of a
        # hypothetical consecutive common zero; replay the underlying identity.
        n = 3
        ll = beta + a*(n-1)-1
        mm = ll + a
        assert sp.expand(U(n,d,mm,z) - z*U(n-1,d,mm,z)
                         + rising(mm+1,d)*U(n-1,d,mm+d,z)) == 0
        moment_cases += 1

    # Corrupt a selected polynomial without changing its leading coefficient.
    # The recurrence must reject it rather than accepting a degree-only control.
    assert equation(U(3, 3, s, z)+1, 4, 3, z) != 0
    report = {
        "kind": "internal bounded exact replay; universal proof is the TeX source",
        "claims": {
            "two_three_separation": "all integer d>=1 and real y>=x>0; elementary proof",
            "lower_exponent_initial_block": "all a>=0,d>=1 integers and beta>0, terminal<=3",
            "parameter_polynomial_selection": "all N>=2,d>=1,z complex; elementary degree proof",
            "all_index_binomial_sfc": False,
        },
        "coverage": {
            "symbolic_recurrence_cases": cases,
            "polynomial_kernel_rank_specializations_per_case": 4,
            "rational_inequality_cases": inequality_cases,
            "exact_moment_gcd_cases": moment_cases,
        },
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "source_sha256": {
            name: sha256((HERE/name).read_bytes()).hexdigest()
            for name in ["lower_exponent_initial_and_boundary.tex",
                         "lower_exponent_initial_and_boundary.py"]
        },
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
        "passed": True,
    }
    return report

if __name__ == "__main__":
    report = replay()
    (HERE/"lower_exponent_initial_and_boundary_verification.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))

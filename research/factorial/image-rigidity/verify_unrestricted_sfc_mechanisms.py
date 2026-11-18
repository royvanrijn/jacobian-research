#!/usr/bin/env python3
"""Exact controls for UNRESTRICTED_SFC_MECHANISMS.tex."""
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path
import hashlib, json

def moments(nmax):
    return [sum(comb(n,k)*factorial(2*k) for k in range(n+1))
            for n in range(nmax+1)]


def rising(q, k):
    return prod(q+j for j in range(k))


def gamma_block_moment(columns, k, n):
    """Mixed moment from independent Gamma/Beta moment formulas."""
    b = [sum(a[i]*kj for a, kj in zip(columns, k))
         for i in range(len(columns[0]))]
    return (Fraction(n)**(sum(b)-sum(k))
            * prod(rising(n+1, kj) for kj in k)
            * prod(Fraction(factorial(z), rising(n+1, z)) for z in b))


def euler_block_multiplier(columns, k, n):
    """The finite product of diagonal Euler factors, separately evaluated."""
    b = [sum(a[i]*kj for a, kj in zip(columns, k))
         for i in range(len(columns[0]))]
    value = Fraction(n)**(sum(b)-sum(k))
    for ell in range(1, n+1):
        value *= prod(Fraction(ell+kj, ell) for kj in k)
        value /= prod(Fraction(ell+z, ell) for z in b)
    return value


def check_gamma_controls():
    f = lambda k: (k+1)**2*(k+2)
    coefficients = [[f(k)*f(m-k) for k in range(m+1)]
                    for m in range(5)]
    assert coefficients[0] == [4]
    assert coefficients[1] == [24, 24]
    assert coefficients[2] == [72, 144, 72]
    values = [sum((-1)**k*term for k, term in enumerate(row))
              for row in coefficients]
    assert values == [4, 0, 0, 0, -24]
    scalar = [Fraction(factorial(k)*f(k), f(0)) for k in range(4)]
    assert scalar == [1, 6, 36, 240]
    assert scalar[2] == scalar[1]**2
    assert scalar[3] != scalar[1]**3
    bad_f = lambda k: f(k)-(1 if k == 2 else 0)
    assert sum((-1)**k*bad_f(k)*bad_f(2-k) for k in range(3)) != 0

    multiplier_controls = 0
    samples = [
        (((0,), (1,), (2,)), [(0,0,0), (1,2,1), (0,4,0)]),
        (((1,0), (0,1), (1,1)), [(1,1,1), (2,0,2), (0,0,3)]),
        (((0,0), (2,1)), [(0,0), (2,3), (1,1)]),
    ]
    for columns, ks in samples:
        for k in ks:
            for n in (1, 2, 5):
                assert (euler_block_multiplier(columns, k, n)
                        == gamma_block_moment(columns, k, n)
                           / prod(factorial(kj) for kj in k))
                multiplier_controls += 1

    columns = ((0,), (1,), (2,))
    for n in range(1, 33):
        a = gamma_block_moment(columns, (2,0,2), n)
        b = gamma_block_moment(columns, (1,2,1), n)
        c = gamma_block_moment(columns, (0,4,0), n)
        target = Fraction(144*(n+2), (n+3)*(n+4))
        assert a-2*b+c == target > 0
        assert a-b+c != target
    assert target < Fraction(144, 32)
    assert (gamma_block_moment(columns, (2,0,2), 1)
            - 2*gamma_block_moment(columns, (1,2,1), 1)
            + gamma_block_moment(columns, (0,4,0), 1)) == Fraction(108,5)
    return {
        "euler_values_at_minus_one_one": values,
        "finite_block_multiplier_controls": multiplier_controls,
        "toric_degeneration_controls": 32,
        "gamma_corrupt_coefficients_rejected": True,
    }



def valuation(a, p):
    count = 0
    while a % p == 0:
        a //= p
        count += 1
    return count


def prime_factors(a):
    factors = set()
    p = 2
    while p*p <= a:
        while a % p == 0:
            factors.add(p)
            a //= p
        p += 1
    if a > 1:
        factors.add(a)
    return factors


def determinant(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    value = Fraction(1)
    for j in range(len(a)):
        pivot = next(i for i in range(j, len(a)) if a[i][j])
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        factor = a[j][j]
        value *= factor
        a[j] = [x/factor for x in a[j]]
        for i in range(j+1, len(a)):
            factor = a[i][j]
            a[i] = [x-factor*y for x, y in zip(a[i], a[j])]
    return value


def check_diagonal_and_saddle_controls():
    a, b = [1,4,24], [1,6,72,720]
    assert [(-2*a[k] if k < 3 else 0)
            + (30*a[k-1] if k else 0)
            + ([3,-16][k] if k < 2 else 0) for k in range(4)] == b
    assert sum(Fraction(c)*Fraction(3,16)**k for k,c in enumerate(a)) == Fraction(83,32)
    sylvester = []
    for coefficients, shifts in ((a[::-1],3), (b[::-1],2)):
        for shift in range(shifts):
            sylvester.append([0]*shift+coefficients
                             +[0]*(5-shift-len(coefficients)))
    resultant = determinant(sylvester)
    assert resultant == 382464 and resultant % 5 != 0

    table = {
        2: (Fraction(-3,2),Fraction(-3,2),Fraction(-1),Fraction(-3,2)),
        3: (Fraction(0),Fraction(-1),Fraction(-2,3),Fraction(-2,3)),
        5: (Fraction(0),Fraction(0),Fraction(0),Fraction(-1)),
        7: (Fraction(0),Fraction(0),Fraction(0),Fraction(0)),
    }
    exceptional = set().union(*(prime_factors(c) for c in a+b))
    assert exceptional == {2,3,5}
    assert a[1]**2-a[0]*a[2] == -8
    assert b[1]**2-b[0]*b[2] == -36
    assert prime_factors(8) | prime_factors(36) <= exceptional
    for p, expected in table.items():
        row = []
        for coefficients in (a,b):
            v = [valuation(c,p) for c in coefficients]
            m = len(v)-1
            row += [max(-Fraction(v[k],k) for k in range(1,m+1)),
                    min(Fraction(v[k]-v[m],m-k) for k in range(m))]
        assert tuple(row) == expected
        l2,u2,l3,u3 = row
        assert l2 >= u3 and l3 >= u2
    for coefficients, primes in ((a,(2,3)), (b,(2,3,5))):
        for p in primes:
            v = [valuation(c,p) for c in coefficients]
            assert len({v[k+1]-v[k] for k in range(len(v)-1)}) > 1

    g = lambda r: -7*r**6+12*r**5-19*r+14
    gp = lambda r: -42*r**5+60*r**4-19
    assert g(Fraction(5,4)) == Fraction(689,4096) > 0
    assert g(Fraction(4,3)) == -Fraction(70,729) < 0
    assert gp(Fraction(5,4)) == -Fraction(353,512) < 0
    assert 8-7*Fraction(5,4) < 0
    u = Fraction(1,3)
    power_sum = 3*u**5
    diagonal = -7/u**2+20*u**3/power_sum-(5*u**4/power_sum)**2
    off_diagonal = -(5*u**4/power_sum)**2
    assert diagonal == -28 and off_diagonal == -25
    for v in ((1,-1,0),(1,1,-2)):
        quadratic = sum(v[i]*v[j]*(diagonal if i==j else off_diagonal)
                        for i in range(3) for j in range(3))
        assert quadratic == -3*sum(x*x for x in v) < 0
    assert -7+Fraction(20,3) == -Fraction(1,3)
    # An intentional sign corruption loses the isolating-interval test.
    assert not (-g(Fraction(5,4)) > 0 and -g(Fraction(4,3)) < 0)
    return {
        "diagonal_normal_form_resultant": int(resultant),
        "diagonal_threshold_prime_classes": ["2","3","5","all other primes"],
        "four_maxima_rational_controls": True,
        "saddle_sign_corruption_rejected": True,
    }



def gamma_resultant_probe():
    """One exact diagnostic; run under Sage, not the default stdlib mode."""
    from sage.all import PolynomialRing, QQ, matrix, factorial, binomial, prod
    from sage.env import SAGE_VERSION

    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polys = PolynomialRing(ring, names=("z", "b", "c"))
    z,b,c = polys.gens()
    variables = (z,b,c)
    degrees = (2,3,4)

    def mu(k):
        return prod(1+h*x for h in range(k))

    def moment(m):
        value = polys.zero()
        for k0 in range(m+1):
            for k1 in range(m-k0+1):
                k3 = m-k0-k1
                coefficient = factorial(m)//(
                    factorial(k0)*factorial(k1)*factorial(k3))
                value += coefficient*z**k0*b**k1*c**k3*mu(k1+3*k3)
        return value

    forms = tuple(moment(d) for d in degrees)
    total_degree = sum(degrees)-len(degrees)+1
    monomials = [z**a*b**bb*c**(total_degree-a-bb)
                 for a in range(total_degree,-1,-1)
                 for bb in range(total_degree-a,-1,-1)]
    rows, nonreduced = [], []
    for row_index, monomial in enumerate(monomials):
        exponents = monomial.exponents()[0]
        choices = [i for i in range(3) if exponents[i]>=degrees[i]]
        assert choices
        i = choices[0]
        row_form = polys(monomial//variables[i]**degrees[i])*forms[i]
        rows.append([row_form.monomial_coefficient(column)
                     for column in monomials])
        if len(choices)>=2:
            nonreduced.append(row_index)
    macaulay = matrix(ring,rows)
    extraneous = macaulay.matrix_from_rows_and_columns(nonreduced,nonreduced)
    numerator, denominator = macaulay.det(), extraneous.det()
    assert denominator != 0
    full = ring(numerator/denominator)
    assert numerator == denominator*full
    assert (macaulay.nrows(),extraneous.nrows()) == (36,10)

    # Independent dehomogenized binary quartic invariants.
    binary = PolynomialRing(ring,"B")
    B = binary.gen()
    raw = [binary(1)] + [
        sum(binomial(r,j)*B**(r-j)*mu(r+2*j) for j in range(r+1))
        for r in range(1,5)]
    m1,m2,m3,m4 = raw[1:]
    k2 = m2-m1*m1
    k3 = m3-3*m2*m1+2*m1**3
    k4 = m4-4*m3*m1-3*m2*m2+12*m2*m1*m1-6*m1**4
    invariant_i = k4+6*k2*k2
    invariant_j = k2*k4+2*k2**3-k3*k3
    assert invariant_i == m4-4*m1*m3+3*m2*m2
    assert (invariant_i.degree(),invariant_j.degree()) == (4,6)
    independent = ring(invariant_i.resultant(invariant_j))
    assert independent == full

    fixed_factor = (x**36*(1+x)**10*(1+2*x)**6*(1+3*x)**3
                    *(1+4*x)**2*(1+5*x)**2)
    content = 60183678025728
    residual = ring(full/(content*fixed_factor))
    assert full == content*fixed_factor*residual
    assert (full.degree(),full.valuation(),len(full.dict())) == (82,36,47)
    assert (residual.degree(),len(residual.dict())) == (23,24)
    assert all(co>0 and co.denominator()==1 for co in full.coefficients())
    assert all(co>0 and co.denominator()==1 for co in residual.coefficients())
    encoded = json.dumps([str(full[i]) for i in range(83)],
                         separators=(",",":")).encode("ascii")
    coefficient_hash = hashlib.sha256(encoded).hexdigest()
    assert coefficient_hash == (
        "2068938a73802314aea04764fb8fc5d62ae5d198d3d5e660bcf380bfd6ee84a0")
    print(json.dumps({
        "status":"PASS_exact_single_diagnostic_only",
        "input":{"A":[[0,1,3]],"s":3,"n":2,"window":[2,3,4]},
        "sage_version":SAGE_VERSION,
        "macaulay_matrix_size":36,
        "extraneous_matrix_size":10,
        "independent_binary_resultant_equal":True,
        "normalization_scalar":1,
        "degree":82,
        "positive_coefficients":47,
        "negative_coefficients":0,
        "coefficient_sha256":coefficient_hash,
        "ordinary_value":str(full(1)),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope":("One falsification diagnostic only. No all-support coefficient "
                 "positivity or unrestricted SFC theorem is established.")
    },indent=2))



def check_row_coalescence_controls():
    def mu(k, x):
        return prod(1+h*x for h in range(k))

    count = 0
    for q in (2,3):
        for u in (0,1,4):
            for x in (Fraction(1,3),Fraction(2)):
                grouped = prod((1+r*x)**u*mu(u,q*x/(1+r*x))
                               for r in range(q))
                assert mu(q*u,x) == grouped
                count += 1
    x = Fraction(1,3)
    assert mu(4,x) != prod((1+r*x)**2*mu(2,x/(1+r*x))
                           for r in range(2))
    # Nonnegative coefficients are needed for nonzero preservation.
    for x in (Fraction(1,3),Fraction(2)):
        y0,y1 = 2*x,2*x/(1+x)
        assert y0-y1-y0*y1/2 == 0
    # A nonprimitive row really can split; this does not test the all-index proof.
    for k in range(9):
        assert rising(1,2*k) == 4**k*rising(Fraction(1,2),k)*rising(1,k)
    assert rising(1,2) != rising(1,2)**2
    return {
        "row_grouping_normalization_controls":count,
        "nonprimitive_factorization_controls":9,
        "corrupt_row_grouping_factor_rejected":True,
        "mixed_sign_specialization_can_vanish":True,
    }


def main():
    row_controls = check_row_coalescence_controls()
    route_controls = check_diagonal_and_saddle_controls()
    gamma_controls = check_gamma_controls()
    ms=moments(32)
    for n in range(2,33):
        assert ms[n]-2*n*(2*n-1)*ms[n-1]+4*n*(n-1)*ms[n-2]==1
        assert ms[n]>=factorial(2*n)
    for n in range(3,33):
        assert (ms[n]-(4*n*n-2*n+1)*ms[n-1]
                +(n-1)*(8*n-6)*ms[n-2]
                -4*(n-1)*(n-2)*ms[n-3])==0
    assert ms[2]-(2*2*(2*2-1)+1)*ms[1]+4*2*(2-1)*ms[0]!=1
    print(json.dumps({
        "status":"PASS_exact_controls_only",
        "inhomogeneous_controls":31,
        "homogeneous_controls":30,
        "first_moments":ms[:5],
        "corrupt_coefficient_rejected":True,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        **gamma_controls,
        **route_controls,
        **row_controls,
        "scope":("Exact finite recurrence, operator, arithmetic and saddle controls only; "
                 "all-index proofs are written. Neither finite-block nor "
                 "ordinary unrestricted SFC is established.")
    },indent=2))

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma-resultant-probe",action="store_true",
                        help="Exact single-case Sage diagnostic; requires Sage Python")
    args = parser.parse_args()
    if args.gamma_resultant_probe:
        gamma_resultant_probe()
    else:
        main()

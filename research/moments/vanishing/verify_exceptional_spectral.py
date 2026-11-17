#!/usr/bin/env python3
"""Exact spectral tests for exceptional finite averages of the normalizer circuit.

The polynomial arithmetic below uses Fraction only. Orthogonality integrals
independently check the triangular associated-Legendre expansions. The
all-parameter orbit-span and all-k spectral formulas use the adjacent proof.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import hashlib
import json
from verify_exceptional_sl2 import group, average_circuit_on_diagonal

ROOT = Path(__file__).resolve().parent


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for n, c in poly.items():
            out[n] += c
    return {n: c for n, c in out.items() if c}


def scale(poly, c):
    c = F(c)
    return {n: c*x for n, x in poly.items() if c*x}


def mul(a, b):
    out = defaultdict(F)
    for n, c in a.items():
        for m, e in b.items():
            out[n+m] += c*e
    return {n: c for n, c in out.items() if c}


def derivative(poly, q):
    return {n-q: c*F(factorial(n), factorial(n-q))
            for n, c in poly.items() if n >= q}


def divide(poly, divisor):
    quotient = {}
    rest = dict(poly)
    top = max(divisor)
    while rest and max(rest) >= top:
        degree = max(rest)-top
        c = rest[max(rest)]/divisor[top]
        quotient[degree] = c
        rest = add(rest, scale({n+degree: x for n, x in divisor.items()}, -c))
    return quotient, rest


@lru_cache(None)
def legendre(l):
    return {l-2*r: F((-1)**r*comb(l, r)*comb(2*l-2*r, l), 2**l)
            for r in range(l//2+1)}


def one_minus_square(n):
    return {2*r: F((-1)**r*comb(n, r)) for r in range(n+1)}


def half_integral(poly):
    return sum((c/F(n+1) for n, c in poly.items() if n % 2 == 0), F(0))


@lru_cache(None)
def associated_derivative(l, q):
    return derivative(legendre(l), q)


@lru_cache(None)
def beta(k, j):
    terms = {
        2*r: F((-1)**r*comb(2*k, r)*comb(k+2*r, j+k))
        for r in range(2*k+1) if 0 <= j+k <= k+2*r
    }
    terms = {n: c for n, c in terms.items() if c}
    if j < 0:
        terms, rest = divide(terms, one_minus_square(-2*j))
        assert not rest, (k, j, rest)
    return terms


@lru_cache(None)
def expansion(k, j):
    """Expand beta in derivatives of ordinary Legendre polynomials."""
    rest = dict(beta(k, j))
    q = 2*abs(j)
    coefficients = {}
    while rest:
        l = max(rest)+q
        basis = associated_derivative(l, q)
        c = rest[max(rest)]/basis[max(basis)]
        coefficients[l] = c
        rest = add(rest, scale(basis, -c))
    return coefficients


def coefficient_by_integral(k, j, l):
    q = 2*abs(j)
    if l < q:
        return F(0)
    norm = F(factorial(l+q), (2*l+1)*factorial(l-q))
    integrand = mul(mul(beta(k, j), one_minus_square(q)),
                    associated_derivative(l, q))
    return half_integral(integrand)/norm


def harmonic_square(k, check_integrals=True, pair_factor=2):
    coeffs = {j: expansion(k, j) for j in range(-k, k+1)}
    checked = 0
    if check_integrals:
        for j, values in coeffs.items():
            for l in range(2*abs(j), 4*k+2*j+1, 2):
                assert values.get(l, F(0)) == coefficient_by_integral(k, j, l), (k, j, l)
                checked += 1
    out = {}
    for l in range(0, 4*k+1, 2):
        value = coeffs[0].get(l, F(0))**2
        for j in range(1, min(k, l//2)+1):
            value += pair_factor*coeffs[j].get(l, F(0))*coeffs[-j].get(l, F(0))*F(
                factorial(l+2*j), factorial(l-2*j))
        out[l] = value/F(2*l+1)
    assert all(isinstance(c, F) for c in out.values())
    if pair_factor == 2:
        # Trivial projection gives h(P_N^(2k))=0 exactly.
        assert sum(out.values(), F(0)) == 0
    return out, checked


def zonal_projection_weights(values, d, highest_l):
    def add_pair(a, b):
        return (a[0]+b[0], a[1]+b[1])

    def mul_pair(a, b):
        return (a[0]*b[0]+d*a[1]*b[1], a[0]*b[1]+a[1]*b[0])

    def scale_pair(a, c):
        return (a[0]*c, a[1]*c)

    tau = Counter(
        add_pair(add_pair(mul_pair(a, a), mul_pair(b, b)),
                 scale_pair(add_pair(mul_pair(c, c), mul_pair(e, e)), -1))
        for a, b, c, e in values)
    out = [(F(0), F(0)) for _ in range(highest_l+1)]
    for t, multiplicity in tau.items():
        previous = (F(0), F(0))
        current = (F(1), F(0))
        for l in range(highest_l+1):
            out[l] = add_pair(out[l], scale_pair(current, F(multiplicity, len(values))))
            # Independent explicit-coefficient evaluation checks recurrence.
            direct = (F(0), F(0))
            for n, c in legendre(l).items():
                power = (F(1), F(0))
                for _ in range(n):
                    power = mul_pair(power, t)
                direct = add_pair(direct, scale_pair(power, c))
            assert direct == current
            nxt = scale_pair(
                add_pair(scale_pair(mul_pair(t, current), 2*l+1),
                         scale_pair(previous, -l)), F(1, l+1))
            previous, current = current, nxt
    assert out[0] == (F(1), F(0))
    return out


def invariant_dimension(n, degrees):
    a, b, c, relation = degrees

    def count(target):
        if target < 0:
            return 0
        return sum(1 for i in range(target//a+1) for j in range(target//b+1)
                   if target-a*i-b*j >= 0 and (target-a*i-b*j) % c == 0)
    return count(n)-count(n-relation)


def unipotent_certificate(values, d):
    """Direct finite-group polynomial calculation for a_s=[[1,s],[0,1]]."""
    rz=(F(0),F(0))
    def ra(a,b): return (a[0]+b[0],a[1]+b[1])
    def rs(a,c): return (a[0]*c,a[1]*c)
    def rm(a,b): return (a[0]*b[0]+d*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    z=(rz,rz); o=((F(1),F(0)),rz)
    def ca(a,b): return (ra(a[0],b[0]),ra(a[1],b[1]))
    def cs(a,c): return (rs(a[0],c),rs(a[1],c))
    def cm(a,b): return (ra(rm(a[0],b[0]),rs(rm(a[1],b[1]),-1)),
                          ra(rm(a[0],b[1]),rm(a[1],b[0])))
    def pa(*polys):
        out={}
        for poly in polys:
            for n,c in poly.items(): out[n]=ca(out.get(n,z),c)
        return {n:c for n,c in out.items() if c!=z}
    def ps(poly,c): return {n:cs(x,c) for n,x in poly.items() if cs(x,c)!=z}
    def pm(a,b):
        out={}
        for n,c in a.items():
            for m,e in b.items(): out[n+m]=ca(out.get(n+m,z),cm(c,e))
        return {n:c for n,c in out.items() if c!=z}
    def pp(poly,n):
        out={0:o}
        for _ in range(n): out=pm(out,poly)
        return out
    def peval(l,poly):
        return pa(*(ps(pp(poly,n),c) for n,c in legendre(l).items()))
    neutral={}; moment={}
    for a,b,c,e in values:
        tau=ra(ra(rm(a,a),rm(b,b)),rs(ra(rm(c,c),rm(e,e)),-1))
        X=(a,b); Z=(rs(c,-1),e); delta=(rz,rs(b,2))
        right={0:(tau,rz),1:cs(cm(X,Z),2)}
        conjugated={0:(tau,rz),1:cs(cm(Z,delta),2),2:cs(cm(Z,Z),-2)}
        neutral=pa(neutral,{0:o},ps(pp(right,2),-6),ps(pp(right,4),5))
        moment=pa(moment,ps(pa(peval(4,conjugated),ps(peval(2,conjugated),-1)),F(64,441)))
    neutral=ps(neutral,F(3,2*len(values)))
    moment=ps(moment,F(1,len(values)))
    def rational(poly):
        assert all(c[0][1]==0 and c[1]==rz for c in poly.values())
        return {n:c[0][0] for n,c in poly.items()}
    return rational(neutral),rational(moment)


def run():
    cases = {
        '2T': (2, [], (6, 8, 12, 24)),
        '2O': (2, [[(0, F(1, 2)), (0, F(1, 2)), (0, 0), (0, 0)]], (8, 12, 18, 36)),
        '2I': (5, [[(0, 0), (F(1, 2), 0), (F(1, 4), -F(1, 4)),
                   (-F(1, 4), -F(1, 4))]], (12, 20, 30, 60)),
    }
    support = sorted({l for j in range(-1, 5) for l in expansion(1, j)})
    assert support == [2, 4, 6, 8, 10, 12]
    dimensions = {}
    weights = {}
    unipotent = {}
    for name, (d, extra, degrees) in cases.items():
        values = group(d, extra)[0]
        lam, moment = unipotent_certificate(values, d)
        if name == '2I':
            assert lam == moment == {}
        else:
            assert lam == {0: F(1), 4: F(5)}, (name, lam)
            assert moment == {0: F(16, 189), 4: F(160, 189), 8: F(400, 189)}, (name, moment)
            coefficients, axis = average_circuit_on_diagonal(values, d)
            rz = (F(0), F(0))
            expected = {tuple((F(sign), F(0)) if j == i else rz for j in range(3))
                        for i in range(3) for sign in [-1, 1]}
            assert axis == expected
            assert coefficients[0] == ((F(2, 3), F(0)), rz)
        unipotent[name] = {
            'neutral_coefficient_polynomial': {str(n): str(c) for n, c in lam.items()},
            'second_moment_polynomial': {str(n): str(c) for n, c in moment.items()}}
        if name != '2T':
            weights[name] = zonal_projection_weights(values, d, 24)
            for l, value in enumerate(weights[name]):
                if invariant_dimension(2*l, degrees) == 0:
                    assert value == (0, 0), (name, l, value)
        multiplicities = {2*l: invariant_dimension(2*l, degrees) for l in support}
        dimensions[name] = {
            'right_types_and_fixed_multiplicities': multiplicities,
            'orbit_span_dimension': sum(multiplicities.values()),
            'negative_type_V4_fixed_dimension': invariant_dimension(4, degrees),
            'neutral_type_V8_fixed_dimension': invariant_dimension(8, degrees),
        }
    assert [dimensions[n]['orbit_span_dimension'] for n in cases] == [9, 6, 3]
    assert all(x['negative_type_V4_fixed_dimension'] == 0 for x in dimensions.values())
    assert [dimensions[n]['neutral_type_V8_fixed_dimension'] for n in cases] == [1, 1, 0]
    rows = []
    integrals = 0
    for k in range(1, 7):
        b, n = harmonic_square(k)
        integrals += n
        moments = {}
        for name, w in weights.items():
            value = tuple(sum((c*w[l][i] for l, c in b.items()), F(0)) for i in (0, 1))
            assert value[1] == 0
            moments[name] = str(value[0])
        if k == 1:
            assert moments == {'2O': '16/189', '2I': '0'}
        if k == 2:
            assert moments == {'2O': '-15616/153153', '2I': '4304/459459'}
        if k >= 2:
            assert all(F(value) != 0 for value in moments.values())
        rows.append({'power_before_averaging': k, 'second_pure_moments': moments,
                     'unprojected_harmonic_square_contributions':
                         {str(l): str(c) for l, c in b.items()}})
    wrong, _ = harmonic_square(1, check_integrals=False, pair_factor=1)
    assert sum(wrong.values(), F(0)) != 0
    dependency = ROOT/'verify_exceptional_sl2.py'
    old_receipt = json.loads((ROOT/'exceptional_sl2_verification.json').read_text())
    assert old_receipt['checker_sha256'] == hashlib.sha256(dependency.read_bytes()).hexdigest()
    for name, (_, _, degrees) in cases.items():
        retained = old_receipt['groups'][name]['symmetric_invariant_dimensions_n0_to32']
        assert retained == [invariant_dimension(n, degrees) for n in range(33)]
    return {
        'completed': True,
        'normalizer_right_harmonic_degrees': support,
        'orbit_span_dimension_certificates': dimensions,
        'second_moment_rows': rows,
        'unipotent_all_parameter_polynomial_certificates': unipotent,
        'independent_associated_legendre_integral_checks': integrals,
        'zonal_weights': {name: [[str(v[0]), str(v[1])] for v in w]
                         for name, w in weights.items()},
        'normalizer_pure_identity_checks': 6,
        'missing_pair_factor_mutation_detected': True,
        'exact_fraction_types_enforced': True,
        'orbit_span_all_parameter_proof_verified_by_checker': False,
        'all_higher_powers_excluded': False,
        'exceptional_model_solved': False,
        'independent_replay': False, 'formal_verification': False, 'external_review': False,
        'dependency_checker_sha256': hashlib.sha256(dependency.read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == '__main__':
    receipt = run()
    (ROOT/'exceptional_spectral_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: v for k, v in receipt.items()
                      if k not in {'zonal_weights', 'second_moment_rows'}}, indent=2))
    for row in receipt['second_moment_rows']:
        print(row['power_before_averaging'], row['second_pure_moments'])

#!/usr/bin/env python3
"""Exact finite certificates for Haar-preserving shear witnesses.

Only the Python standard library is used. The all-group and all-order
arguments are written in HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent


class Cyclotomic4:
    """Q[z]/(z^4+c3*z^3+c2*z^2+c1*z+c0), for Phi5 or Phi8."""
    def __init__(self, modulus):
        self.modulus = tuple(F(x) for x in modulus)
        self.zero = (F(0),)*4
        self.one = (F(1), F(0), F(0), F(0))
        self.z = (F(0), F(1), F(0), F(0))

    def add(self, a, b):
        return tuple(x+y for x, y in zip(a, b))

    def scale(self, a, c):
        return tuple(x*F(c) for x in a)

    def mul(self, a, b):
        out = [F(0)]*7
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i+j] += x*y
        for i in range(6, 3, -1):
            for j, c in enumerate(self.modulus):
                out[i-4+j] -= out[i]*c
        return tuple(out[:4])

    def power(self, a, n):
        out = self.one
        for _ in range(n):
            out = self.mul(out, a)
        return out

    def conjugate(self, a):
        invz = self.power(self.z, 4) if self.modulus == (1, 1, 1, 1) else self.scale(self.power(self.z, 3), -1)
        out = self.zero
        for n, c in enumerate(a):
            out = self.add(out, self.scale(self.power(invz, n), c))
        return out


def matrix_group(field, gens, expected_order):
    add, mul, sc = field.add, field.mul, field.scale
    o, z = field.one, field.zero
    identity = (o, z, z, o)

    def mm(a, b):
        return (add(mul(a[0], b[0]), mul(a[1], b[2])),
                add(mul(a[0], b[1]), mul(a[1], b[3])),
                add(mul(a[2], b[0]), mul(a[3], b[2])),
                add(mul(a[2], b[1]), mul(a[3], b[3])))

    for a, b, c, d in gens:
        assert add(mul(a, d), sc(mul(b, c), -1)) == o
        assert d == field.conjugate(a) and c == sc(field.conjugate(b), -1)
    values = {identity}
    queue = deque([identity])
    while queue:
        a = queue.popleft()
        for g in gens:
            b = mm(a, g)
            if b not in values:
                values.add(b)
                queue.append(b)
                assert len(values) <= expected_order
    assert len(values) == expected_order
    for a, b, c, d in values:
        assert (d, sc(b, -1), sc(c, -1), a) in values

    def canonical(a):
        return min(a, tuple(sc(x, -1) for x in a))
    projective = {canonical(a) for a in values}
    counts = Counter()
    for a in projective:
        power = identity
        for k in range(1, 61):
            power = mm(power, a)
            if canonical(power) == canonical(identity):
                counts[k] += 1
                break
        else:
            raise AssertionError('Unexpected projective order')
    return values, dict(sorted(counts.items()))


def invariant_form(field, generators, coeffs, n):
    for a, b, c, d in generators:
        powers = [[field.power(x, j) for j in range(n+1)] for x in (a, b, c, d)]
        result = [field.zero]*(n+1)
        for i, coefficient in coeffs.items():
            for r in range(n-i+1):
                for s in range(i+1):
                    term = field.mul(field.mul(powers[0][n-i-r], powers[2][r]),
                                     field.mul(powers[1][i-s], powers[3][s]))
                    term = field.scale(term, coefficient*comb(n-i, r)*comb(i, s))
                    result[r+s] = field.add(result[r+s], term)
        assert result == [field.scale(field.one, coeffs.get(j, 0)) for j in range(n+1)]


def invariant_dimensions(field, group, n):
    sums = [field.zero]*(n+1)
    for a, b, c, d in group:
        trace = field.add(a, d)
        old, now = field.zero, field.one
        for j in range(n+1):
            sums[j] = field.add(sums[j], now)
            old, now = now, field.add(field.mul(trace, now), field.scale(old, -1))
    dims = []
    for total in sums:
        mean = field.scale(total, F(1, len(group)))
        assert mean[1:] == (0, 0, 0) and mean[0].denominator == 1
        dims.append(int(mean[0]))
    return dims


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for mon, c in poly.items():
            out[mon] += c
    return {mon: c for mon, c in out.items() if c}


def scale(poly, c):
    return {mon: F(c)*x for mon, x in poly.items() if c*x}


def mul(a, b):
    out = defaultdict(F)
    for ma, ca in a.items():
        for mb, cb in b.items():
            out[tuple(x+y for x, y in zip(ma, mb))] += ca*cb
    return {mon: c for mon, c in out.items() if c}


def power(poly, k):
    out = {(0, 0, 0, 0): F(1)}
    for _ in range(k):
        out = mul(out, poly)
    return out


def diff(poly, i):
    out = {}
    for mon, c in poly.items():
        if mon[i]:
            new = list(mon)
            new[i] -= 1
            out[tuple(new)] = c*mon[i]
    return out


VARS = [{tuple(int(i == j) for j in range(4)): F(1)} for i in range(4)]
X, Y, Z, W = VARS
DET = add(mul(X, W), scale(mul(Y, Z), -1))


def lower(poly):
    return add(mul(Z, diff(poly, 0)), mul(W, diff(poly, 1)))


def raise_(poly):
    return add(mul(X, diff(poly, 2)), mul(Y, diff(poly, 3)))


def omega(poly):
    return add(diff(diff(poly, 0), 3), scale(diff(diff(poly, 1), 2), -1))


def casimir(poly):
    return add({m: (m[0]+m[1]-m[2]-m[3])*(m[0]+m[1]-m[2]-m[3]+2)*c
                for m, c in poly.items()},
               scale(lower(raise_(poly)), 4))


def projection(poly, target, degree):
    out = poly
    eigen = target*(target+2)
    for d in range(degree % 2, degree+1, 2):
        if d != target:
            other = d*(d+2)
            out = scale(add(casimir(out), scale(out, -other)), F(1, eigen-other))
    return out


def haar(poly):
    return sum((c*((-1)**b)*F(factorial(a)*factorial(b), factorial(a+b+1))
                for (a, b, cc, d), c in poly.items() if a == d and b == cc), F(0))


def compact_star(poly):
    # On SU2: conjugate(x)=w, conjugate(y)=-z,
    # conjugate(z)=-y, conjugate(w)=x. Coefficients here are rational.
    return {(d, c, b, a): coefficient*((-1)**(b+c))
            for (a, b, c, d), coefficient in poly.items()}


def binary_power(coeffs, k):
    out = {0: F(1)}
    for _ in range(k):
        new = defaultdict(F)
        for i, c in out.items():
            for j, e in coeffs.items():
                new[i+j] += c*e
        out = {i: c for i, c in new.items() if c}
    return out


def norm_moment(coeffs, n, s):
    return sum((F(c*c, (n*s+1)*comb(n*s, i))
                for i, c in binary_power(coeffs, s).items()), F(0))


def transvectant(coeffs, n, r):
    # Binary homogeneous forms, indexed by the exponent of Y.
    def derivative(a, degree, dx, dy):
        return {i-dy: c*F(factorial(degree-i), factorial(degree-i-dx))*F(factorial(i), factorial(i-dy))
                for i, c in a.items() if degree-i >= dx and i >= dy}
    result = defaultdict(F)
    for j in range(r+1):
        left = derivative(coeffs, n, r-j, j)
        right = derivative(coeffs, n, j, r-j)
        for a, c in left.items():
            for b, e in right.items():
                result[a+b] += ((-1)**j)*comb(r, j)*c*e
    return {i: c for i, c in result.items() if c}


def verify_case(name):
    if name == '2O':
        field = Cyclotomic4((1, 0, 0, 0))
        z, o, zz = field.z, field.one, field.zero
        i = field.power(z, 2)
        half = lambda a: field.scale(a, F(1, 2))
        generators = [
            (z, zz, zz, field.scale(field.power(z, 3), -1)),
            (zz, o, field.scale(o, -1), zz),
            (half(field.add(o, i)), half(field.add(o, i)),
             half(field.add(field.scale(o, -1), i)), half(field.add(o, field.scale(i, -1)))),
        ]
        n, coeffs, order = 8, {0: F(1), 4: F(14), 8: F(1)}, 48
        expected_orders = {1: 1, 2: 9, 3: 8, 4: 6}
        expected_omega, expected_kappa, expected_binary = F(672), F(7, 4290), F(1568, 143)
        expected_norm, expected_transvectant = F(8, 15), F(1128960)
    else:
        field = Cyclotomic4((1, 1, 1, 1))
        z, o, zz = field.z, field.one, field.zero
        z2, z3, z4 = (field.power(z, k) for k in (2, 3, 4))
        root5 = field.add(o, field.scale(field.add(z, z4), 2))
        a = field.scale(field.mul(field.add(z4, field.scale(z, -1)), root5), F(1, 5))
        b = field.scale(field.mul(field.add(z2, field.scale(z3, -1)), root5), F(1, 5))
        generators = [
            (field.scale(z3, -1), zz, zz, field.scale(z2, -1)),
            (zz, o, field.scale(o, -1), zz),
            (a, b, b, field.scale(a, -1)),
        ]
        n, coeffs, order = 12, {1: F(1), 6: F(11), 11: F(-1)}, 120
        expected_orders = {1: 1, 2: 15, 3: 20, 5: 24}
        expected_omega, expected_kappa, expected_binary = F(-39600), F(-11, 3907008), F(-3630, 323)
        expected_norm, expected_transvectant = F(25, 1092), F(-26345088000)
    values, orders = matrix_group(field, generators, order)
    assert orders == expected_orders
    invariant_form(field, generators, coeffs, n)
    dims = invariant_dimensions(field, values, n)
    assert dims == [1]+[0]*(n-1)+[1]
    assert all(coeffs.get(i, 0) == (-1)**i*coeffs.get(n-i, 0) for i in range(n+1))
    p = {(n-i, i, 0, 0): c for i, c in coeffs.items()}
    q = {(0, 0, n-i, i): c for i, c in coeffs.items()}
    assert not lower(q) and not raise_(p)
    assert compact_star(p) == q
    norm_lowered = haar(mul(lower(p), compact_star(lower(p))))
    invalid_generator_value = haar(mul(p, lower(compact_star(lower(p)))))
    assert norm_lowered > 0 and invalid_generator_value == -norm_lowered
    lowered = p
    for _ in range(n):
        lowered = lower(lowered)
    assert lowered == scale(q, factorial(n))
    assert not lower(lowered)
    assert haar(mul(p, q)) == expected_norm == norm_moment(coeffs, n, 1)

    r = n//2
    product = mul(q, lower(p))
    contracted = product
    for _ in range(r):
        contracted = omega(contracted)
    descendant = p
    for _ in range(r+1):
        descendant = lower(descendant)
    assert contracted == scale(descendant, expected_omega)
    assert not omega(contracted)
    constant = F(factorial(r)*factorial(n+r+1), factorial(n+1))
    assert expected_omega/constant == expected_kappa
    assert expected_kappa*F(factorial(n), factorial(r-1)) == expected_binary
    direct = projection(product, n, 2*n)
    expected = mul(power(DET, r), scale(descendant, expected_kappa))
    assert direct == expected  # independent Casimir versus Omega projection
    tv = transvectant(coeffs, n, r)
    assert tv == {i: expected_transvectant*c for i, c in coeffs.items()}

    first_m = r+1
    relevant = power(p, first_m)
    for _ in range(r):
        relevant = lower(relevant)
    relevant = scale(mul(power(q, r), relevant), F(1, factorial(r)))
    contributing = sum(1 for a, b, c, d in relevant if a == d and b == c)
    assert contributing >= 2
    assert haar(relevant) == 0  # first weight-allowed pure moment, directly
    norms = {s: norm_moment(coeffs, n, s) for s in (1, 2, 3, n+1, 2*(n+1))}
    assert all(c > 0 for c in norms.values())
    sharp = []
    for s in (1, 2):
        e, m = n*s, (n+1)*s
        assert n*m == (n+1)*e
        extremal = power(p, s)
        for _ in range(n*s):
            extremal = lower(extremal)
        assert extremal == scale(power(q, s), factorial(n*s))
        assert haar(mul(power(p, m), power(q, m))) == norms[m]
        sharp.append({'marker': 'p^'+str(s), 'matrix_entry_degree': e,
                      'last_nonzero_index': m, 'moment': str(norms[m])})
    return {
        'n': n, 'binary_form_coefficients_by_Y_exponent': {str(i): str(c) for i, c in coeffs.items()},
        'cyclotomic_modulus_low_coefficients': [str(c) for c in field.modulus],
        'group_order': len(values), 'projective_element_orders': orders,
        'group_generator_form_invariance': True, 'compact_reality_condition': True,
        'invalid_coefficient_mutation_detected': True,
        'invalid_pF_generator_integral': str(invalid_generator_value),
        'fixed_dimensions_n0_through_n': dims,
        'shear_witness': 'P=p(x+q*z,y+q*w), q=p(z,w)',
        'exact_matrix_entry_degree': n*(n+1),
        'principal_substitution_weights': [n-j*(n+2) for j in range(n+1)],
        'separating_invariant_degree_in_P': n+2,
        'separating_invariant_value': str(norms[n+1]),
        'marked_norm_C1': str(norms[1]),
        'mixed_cutoff_slope': str(F(n+1, n)),
        'sharp_marker_controls': sharp,
        'first_resonant_pure_index': first_m,
        'first_resonant_direct_haar': '0',
        'first_resonant_balanced_monomials': contributing,
        'omega_projection_scalar': str(expected_omega),
        'determinant_projection_normalizer': str(constant),
        'Vn_descendant_projection_coefficient': str(expected_kappa),
        'Vn_binary_form_second_coefficient': str(expected_binary),
        'maximum_projected_root_multiplicity': r-1,
        'nullcone_root_threshold': r,
        'casimir_and_omega_projection_agree': True,
        'self_transvectant_scalar': str(expected_transvectant),
    }


def principal_sl3_control():
    def mm(a, b):
        return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    def sc(a, c):
        return [[c*x for x in row] for row in a]
    def sub(a, b):
        return [[x-y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]
    def bracket(a, b):
        return sub(mm(a, b), mm(b, a))
    E = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
    D = [[0, 0, 0], [2, 0, 0], [0, 2, 0]]
    H = [[2, 0, 0], [0, 0, 0], [0, 0, -2]]
    assert bracket(E, D) == H and bracket(H, E) == sc(E, 2) and bracket(H, D) == sc(D, -2)
    top = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
    v = top
    for _ in range(4):
        v = bracket(D, v)
    assert v == [[0, 0, 0], [0, 0, 0], [96, 0, 0]]
    assert bracket(D, v) == [[0]*3 for _ in range(3)]
    return {'principal_triple_checked': True, 'unique_top_adjoint_SL2_weight': 4,
            'fourth_lowering_endpoint_coefficient': 96}


def main():
    rows = {}
    for name in ('2O', '2I'):
        rows[name] = verify_case(name)
        print(name, 'degree', rows[name]['exact_matrix_entry_degree'],
              'stable projection coefficient', rows[name]['Vn_binary_form_second_coefficient'],
              'first resonant moment', rows[name]['first_resonant_direct_haar'], flush=True)
    receipt = {
        'completed': True, 'cases': rows,
        'principal_sl3_control': principal_sl3_control(),
        'all_group_classification_verified_by_checker': False,
        'all_order_pure_and_mixed_proofs_verified_by_checker': False,
        'exceptional_nullcone_counterexamples_constructed': True,
        'Mathieu_Zhao_failure_witness_constructed': False,
        'independent_replay': False, 'formal_verification': False, 'external_review': False,
        'python_version': sys.version.split()[0],
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (ROOT/'haar_shear_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'cases'}, indent=2))


if __name__ == '__main__':
    main()

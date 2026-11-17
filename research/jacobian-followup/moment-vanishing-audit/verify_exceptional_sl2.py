#!/usr/bin/env python3
"""Exact finite-group and all-degree character-series certificates.

Finite group generation uses rational pairs in Q(sqrt(d)), with no floating
point arithmetic. The Molien identity is a rational-function equality, not
an extrapolation from a finite power-series prefix.
"""
from collections import Counter, deque
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sympy as S

ROOT = Path(__file__).resolve().parent


def group(d, extra):
    zero = (F(0), F(0))
    one = (F(1), F(0))

    def scalar(a, b=0):
        return (F(a), F(b))

    def add(x, y):
        return (x[0]+y[0], x[1]+y[1])

    def neg(x):
        return (-x[0], -x[1])

    def mul(x, y):
        a, b = x
        c, e = y
        return (a*c+d*b*e, a*e+b*c)

    def qneg(q):
        return tuple(neg(x) for x in q)

    def conj(q):
        return (q[0],)+tuple(neg(x) for x in q[1:])

    def qmul(q, r):
        # Quaternion formula a0*b0-dot(a,b), a0*b+b0*a+a cross b.
        a, b, c, e = q
        f, g, h, j = r
        return (
            add(add(mul(a, f), neg(mul(b, g))), neg(add(mul(c, h), mul(e, j)))),
            add(add(mul(a, g), mul(b, f)), add(mul(c, j), neg(mul(e, h)))),
            add(add(mul(a, h), mul(c, f)), add(mul(e, g), neg(mul(b, j)))),
            add(add(mul(a, j), mul(e, f)), add(mul(b, h), neg(mul(c, g)))),
        )

    identity = (one, zero, zero, zero)
    i = (zero, one, zero, zero)
    j = (zero, zero, one, zero)
    hurwitz = (scalar(F(1, 2)),)*4
    generators = [i, j, hurwitz] + [tuple(scalar(*x) for x in e) for e in extra]
    values = {identity}
    queue = deque([identity])
    while queue:
        x = queue.popleft()
        for gen in generators:
            y = qmul(x, gen)
            assert qmul(y, conj(y)) == identity
            if y not in values:
                values.add(y)
                queue.append(y)
                assert len(values) <= 120, 'Unexpected infinite/oversized generated group'
    assert all(conj(q) in values for q in values)
    assert qmul(qmul(qmul(i, j), conj(i)), conj(j)) == qneg(identity)
    # A common invariant line for i and j would make their commutator act
    # trivially there, whereas -1 acts as -Id in the standard module.

    def canonical(q):
        for x in q:
            if x != zero:
                # This chooses between q and -q lexicographically; it does
                # not purport to order real quadratic numbers.
                return min(q, qneg(q))
        raise AssertionError('zero quaternion')

    projective = {canonical(q) for q in values}
    order_counts = Counter()
    for q in projective:
        power = identity
        for n in range(1, 61):
            power = qmul(power, q)
            if power in (identity, qneg(identity)):
                order_counts[n] += 1
                break
        else:
            raise AssertionError('Projective order exceeded 60')
    traces = Counter(q[0] for q in values)
    return values, traces, order_counts


def molien(d, traces, order):
    t = S.symbols('t')
    rational = lambda x: S.Rational(x.numerator, x.denominator)
    result = S.Integer(0)
    used = set()
    for (a, b), multiplicity in traces.items():
        if (a, b) in used:
            continue
        A = 1-2*rational(a)*t+t**2
        if not b:
            result += multiplicity/A
        else:
            assert traces[(a, -b)] == multiplicity
            result += 2*multiplicity*A/(A**2-4*d*rational(b)**2*t**2)
            used.add((a, -b))
        used.add((a, b))
    return S.cancel(result/order)


def average_circuit_on_diagonal(values, d):
    """Evaluate each left-weight component of the averaged circuit at 1."""
    zero = (F(0), F(0))
    def add(x, y):
        return (x[0]+y[0], x[1]+y[1])
    def neg(x):
        return (-x[0], -x[1])
    def mul(x, y):
        return (x[0]*y[0]+d*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    def scale(x, n):
        return (n*x[0], n*x[1])
    def cmul(x, y):
        return (add(mul(x[0], y[0]), neg(mul(x[1], y[1]))),
                add(mul(x[0], y[1]), mul(x[1], y[0])))
    totals = {k: (zero, zero) for k in [-4, 0, 4, 8, 12, 16]}
    axis = set()
    for a, b, c, e in values:
        tau = add(add(mul(a,a), mul(b,b)), neg(add(mul(c,c), mul(e,e))))
        axis.add((tau, scale(add(mul(b,c), mul(a,e)),2),
                  scale(add(mul(b,e), neg(mul(a,c))),2)))
        ss = mul(tau,tau)
        s2 = mul(ss,ss)
        ur = scale(add(mul(a,c),neg(mul(b,e))),-2)
        ui = scale(add(mul(a,e),mul(b,c)),-2)
        U = cmul((ur,ui),(ur,ui))
        U2=cmul(U,U); U3=cmul(U2,U); U4=cmul(U3,U)
        def times(z, x):
            return (mul(z[0],x),mul(z[1],x))
        coefficients = {
            -4:(U[0],neg(U[1])),
            0:(add((F(1),F(0)), add(scale(ss,-6),scale(s2,5))),zero),
            4:times(U,add(scale(ss,-6),scale(s2,10))),
            8:times(U2,add(scale(ss,-2),scale(s2,10))),
            12:times(U3,scale(s2,5)),
            16:times(U4,s2),
        }
        for k, z in coefficients.items():
            totals[k]=(add(totals[k][0],z[0]),add(totals[k][1],z[1]))
    totals={k:(scale(z[0],F(1,len(values))),scale(z[1],F(1,len(values))))
            for k,z in totals.items()}
    return totals, axis


def run():
    half = F(1, 2)
    quarter = F(1, 4)
    specifications = {
        '2T': (2, [], 24, (6, 8, 12, 24), {1: 1, 2: 3, 3: 8}),
        '2O': (2, [[(0, half), (0, half), (0, 0), (0, 0)]],
               48, (8, 12, 18, 36), {1: 1, 2: 9, 3: 8, 4: 6}),
        # Flicker's printed p.215 theta, using
        # (a,b,c,d) -> [[a+bi,c+di],[-c+di,a-bi]].
        '2I': (5, [[(0, 0), (half, 0), (quarter, -quarter), (-quarter, -quarter)]],
               120, (12, 20, 30, 60), {1: 1, 2: 15, 3: 20, 5: 24}),
    }
    t = S.symbols('t')
    rows = {}
    groups = {}
    for name, (d, extra, order, degrees, distribution) in specifications.items():
        values, traces, projective = group(d, extra)
        assert len(values) == order
        assert dict(projective) == distribution
        actual = molien(d, traces, order)
        a, b, c, relation = degrees
        expected = (1-t**relation)/((1-t**a)*(1-t**b)*(1-t**c))
        assert S.cancel(actual-expected) == 0
        # Independent symmetric-character recurrence through n=32.
        zero, one = (F(0), F(0)), (F(1), F(0))
        counts = []
        previous = {x: zero for x in traces}
        current = {x: one for x in traces}
        for n in range(33):
            total = [F(0), F(0)]
            for (ra, rb), count in traces.items():
                ca, cb = current[(ra, rb)]
                total[0] += count*ca
                total[1] += count*cb
            assert total[1] == 0 and total[0] % order == 0
            counts.append(int(total[0]/order))
            new = {}
            for x in traces:
                ra, rb = x
                ca, cb = current[x]
                pa, pb = previous[x]
                new[x] = (2*ra*ca+2*d*rb*cb-pa, 2*ra*cb+2*rb*ca-pb)
            previous, current = current, new
        assert counts[0] == 1 and not any(counts[1:a]) and counts[a] == 1
        groups[name] = values
        rows[name] = {
            'quadratic_radicand': d, 'group_order': len(values),
            'projective_order': sum(projective.values()),
            'projective_element_order_counts': dict(sorted(projective.items())),
            'real_part_multiplicities': [
                {'rational': str(x[0]), 'sqrt_coefficient': str(x[1]), 'count': count}
                for x, count in sorted(traces.items())],
            'molien_rational_function': str(S.factor(actual)),
            'hypersurface_series_degrees': list(degrees),
            'all_degree_molien_identity': True,
            'symmetric_invariant_dimensions_n0_to32': counts,
            'least_nontrivial_left_highest_weight': a,
            'standard_invariant_line_obstructed_by_Q8_commutator': True,
        }
    assert groups['2T'] < groups['2O']
    # 2T has rational coordinates, so field labels do not affect inclusion.
    assert groups['2T'] < groups['2I']
    averaged = {}
    for name, d in [('2O',2),('2I',5)]:
        coefficients, axis = average_circuit_on_diagonal(groups[name],d)
        assert coefficients[-4] == ((0,0),(0,0))
        if name == '2I':
            assert coefficients[0] == ((0,0),(0,0))
            assert any(z != ((0,0),(0,0)) for k,z in coefficients.items() if k>0)
        else:
            z=(F(0),F(0)); o=(F(1),F(0))
            expected_axis={tuple((sign*o[0],o[1]) if j==i else z
                                 for j in range(3))
                           for i in range(3) for sign in [-1,1]}
            assert axis == expected_axis
            assert coefficients[0] == ((F(2,3),0),(0,0))
        averaged[name]={str(k):[[str(x) for x in part] for part in z]
                        for k,z in coefficients.items()}
    # Independent sphere monomial calculation for the O neutral polynomial
    # C=(5/3)(x^4+y^4+z^4)-1. E(x^8)=1/9, E(x^4*y^4)=1/105.
    def odd_factorial(n):
        result=1
        for k in range(1,n+1,2):
            result*=k
        return result
    def sphere_even(exponents):
        numerator=1
        for e in exponents:
            numerator*=odd_factorial(2*e-1)
        return F(numerator,odd_factorial(2*sum(exponents)+1))
    fourth=sphere_even((2,0,0))
    eighth=sphere_even((4,0,0))
    mixed=sphere_even((2,2,0))
    assert (fourth,eighth,mixed)==(F(1,5),F(1,9),F(1,105))
    second = F(25,9)*(3*eighth+6*mixed)-F(10,3)*3*fourth+1
    assert second == F(16,189)
    # Removing the added generator fails to produce the exceptional groups.
    assert len(group(2, [])[0]) == 24 != 48
    # Replacing the group by its abstract projective quotient would falsely
    # remove the central -Id and change the odd-degree invariant conclusion.
    assert all(not row['symmetric_invariant_dimensions_n0_to32'][n]
               for row in rows.values() for n in range(1, 33, 2))
    return {
        'completed': True, 'groups': rows,
        'binary_tetrahedral_in_both_larger_groups': True,
        'exact_quaternion_generation': True,
        'averaged_normalizer_circuit_left_diagonal_coefficients': averaged,
        'octahedral_average_second_moment': str(second),
        'icosahedral_average_nonzero_positive_weight_confirmed': True,
        'all_degree_series_proved_by_rational_identity': True,
        'missing_generator_mutation_detected': True,
        'no_exceptional_moment_counterexample_constructed': True,
        'general_reductive_reduction_verified_by_checker': False,
        'independent_replay': False, 'formal_verification': False, 'external_review': False,
        'sympy_version': S.__version__,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == '__main__':
    receipt = run()
    (ROOT/'exceptional_sl2_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'groups'}, indent=2))
    for name, data in receipt['groups'].items():
        print(name, data['group_order'], data['molien_rational_function'],
              'first nonconstant weight', data['least_nontrivial_left_highest_weight'])

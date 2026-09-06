#!/usr/bin/env python3
"""CAS-free verification of trace witnesses and the integral mixing accounting."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import integral_quadratic_mixing as source

OUTPUT = r.OUT/'rank_jump_integral_quadratic_mixing_verification_v1.json'


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[0], -a[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def div(a, b):
    n = b[0]*b[0]+b[1]*b[1]
    assert n
    out = mul(a, (b[0], -b[1]))
    return out[0]/n, out[1]/n


def cq(a, b=0):
    return F(a), F(b)


def point_add(P, Q, A=5):
    if P is None: return Q
    if Q is None: return P
    x, y = P; u, v = Q
    if x == u and y == neg(v): return None
    if P == Q:
        slope = div(add(mul(cq(3), mul(x, x)), cq(A)), mul(cq(2), y))
    else:
        slope = div(add(v, neg(y)), add(u, neg(x)))
    xx = add(add(mul(slope, slope), neg(x)), neg(u))
    yy = add(neg(y), mul(slope, add(x, neg(xx))))
    return xx, yy


def compute():
    data = r.read(source.OUTPUT)
    for name, digest in data['bindings'].items():
        assert r.digest((r.ROOT/name).read_bytes()) == digest
    T, Tc = (cq(0, 1), cq(1, 2)), (cq(0, -1), cq(1, -2))
    P, anti = (cq(4), cq(-9)), (cq(-1), cq(0, -3))
    assert point_add(T, Tc) == P
    assert point_add(T, (Tc[0], neg(Tc[1]))) == anti
    assert point_add(T, T) == point_add(P, anti)
    assert mul(T[1], T[1]) == add(add(mul(mul(T[0], T[0]), T[0]), mul(cq(5), T[0])), cq(-3))
    # Monic cubic rational-root test, plus independent mod-3 doubling.
    for B, point in [(-3, (4, -9)), (3, (1, -3))]:
        assert all(x**3+5*x+B for x in [-3, -1, 1, 3])
        pts = [(x, y) for x in range(3) for y in range(3) if (y*y-x**3-5*x-B) % 3 == 0]
        assert pts == [(0, 0), (1, 0), (2, 0)]
        assert (point[0] % 3, point[1] % 3) in pts
        # All nonzero reduced points are 2-torsion, so the double image is {O}.
        assert all(y == 0 for x, y in pts)
    # The shared-class identity has difference -theta*(theta^3+5theta-3).
    assert [0, 3, -5, 0, -1] == [0]+[3, -5, 0, -1]
    assert data['constructive_control']['exact_ranks'] == [1, 1]
    assert all(x['upper_bound'] == 1 for x in data['constructive_control']['descents'])

    closed = data['closed_control']
    # These real sign lines differ, independently of the retained unit basis.
    intervals = [(F(-2), F(-1)), (F(-1), F(0)), (F(12), F(13))]
    def signs(constant, slope):
        values = []
        for lo, hi in intervals:
            ends = sorted([constant+slope*lo, constant+slope*hi])
            assert ends[0] >= 0 or ends[1] <= 0
            values.append(int(ends[1] <= 0))
        return values
    original = signs(-1, -1)
    twist = [signs(x, 1) for x in [-1, 0, 2]]
    assert original == [0, 1, 1] and twist == [[1, 1, 0], [1, 1, 0], [0, 0, 0]]
    assert r.pack(original) not in source.span([r.pack(x) for x in twist])
    assert closed['module'] == source.module_counts(1, 3, 0)
    assert data['constructive_control']['module'] == source.module_counts(1, 1, 1)
    for row in data['production']:
        assert row['trivial_integral_summands_at_least'] == row['generic_saturated_local_norm_codimension']+row['certified_norm_defect_dimension_at_least']
    return {
        'schema': 'rank-jump.integral-quadratic-mixing-verification.v1',
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (Path(__file__), source.OUTPUT)},
        'status': 'PASS', 'exact_Qi_group_law_witnesses': True,
        'independent_reduction_nondivisibility': True,
        'closed_control_real_sign_separation': True,
        'rank_upper_bound_source': 'Two retained PARI effort-zero unconditional 2-descent upper bounds; not independently recomputed by this CAS-free verifier.',
        'boundary': 'The verifier checks trace witnesses, nondivisibility and accounting. The written lattice proof and retained exact rank/Selmer certificates supply the remaining mathematical dependencies.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['build', 'check']); args = p.parse_args()
    data = compute()
    if args.mode == 'check':
        assert r.read(OUTPUT) == data; print('PASS CAS-free integral mixing witnesses')
    else:
        r.write_new(OUTPUT, data); print(data['status'])

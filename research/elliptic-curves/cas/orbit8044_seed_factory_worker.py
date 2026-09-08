#!/usr/bin/env sage-python
"""Exact orbit8044 chart and two-branch seed producer. No point-search backend."""
import argparse
from fractions import Fraction as F
import json
from math import gcd, isqrt
from pathlib import Path

from sage.all import QQ, ZZ, EllipticCurve, PolynomialRing, matrix, pari
import det1092_funnel as f
import orbit8044_seed_factory as factory
from det1092_funnel_worker import checked_protocol, prepare, first_seed, seal_seed
from det1092_v3_contract import evaluate
from v3_warm_support import atomic, bindings, read, require, sha

ANCHOR = f.ART/'det1092_small_conic_seed_v1/conic-solve.json'


def build_chart(output):
    cover, chart = read(f.COVER), read(f.CHART)
    bindings(f.ROOT, cover['inputs'])
    a, b, c, d = map(int, chart['parameter_matrix'])
    q0, q1, q2 = map(int, cover['curve_over_Q']['q_coefficients'])
    transformed = [q0*d*d+q1*b*d+q2*b*b,
        2*q0*c*d+q1*(a*d+b*c)+2*q2*a*b, q0*c*c+q1*a*c+q2*a*a]
    content = gcd(*transformed); root = isqrt(content)
    require(root*root == content, 'content must be an exact square')
    q0, q1, q2 = [v//content for v in transformed]
    G = matrix(ZZ, [[2*q2, q1, 0], [q1, 2*q0, 0], [0, 0, -2]])
    v = list(map(ZZ, read(ANCHOR)['vector']))
    require((matrix(ZZ, 1, 3, v)*G*matrix(ZZ, 3, 1, v))[0, 0] == 0, 'conic anchor is not isotropic')
    M = matrix(ZZ, pari.qfparam(G, v, 3))
    R = PolynomialRing(QQ, 'u')
    X, Y, Z = [R(list(row)) for row in M.rows()]
    require(Z**2 == q0*Y**2+q1*X*Y+q2*X**2, 'parametrization identity failed')
    require(X.gcd(Y).degree() == 0 and max(X.degree(), Y.degree()) == 2, 'degree-two basepoint-free map required')
    result = dict(status='PASS_EXACT_ORBIT8044_PARAMETER_CHART',
        method='PARI qfparam flag3 on the primitive transformed conic; polynomial identity independently checked',
        coordinate='Factory u is the reduced conic parameter, distinct from the old slope coordinate; s=X(u)/Y(u)',
        s_numerator=list(map(str, M.row(0))), s_denominator=list(map(str, M.row(1))),
        z_numerator=list(map(str, M.row(2))), primitive_conic=list(map(str, [q0, q1, q2])),
        square_content=str(content), original_parameter_matrix=list(map(str, [a, b, c, d])),
        degree_to_reduced_parameter=2,
        inputs={str(p.relative_to(f.ROOT)):sha(p) for p in [f.COVER, f.CHART, ANCHOR, Path(__file__)]},
        boundary='Reparametrization of the existing conic, no new generic-rank theorem. Every specialized extra direction still requires a certificate.')
    atomic(output, result, immutable=True)
    print('PASS_EXACT_ORBIT8044_PARAMETER_CHART', flush=True)


def both_points(row):
    """Solve the residual quadratic itself, retaining both finite branches."""
    cover, chart = read(f.COVER), read(f.CHART)
    s = F(row['parameter']); tau = F(row['original_parameter'])
    R = PolynomialRing(QQ, 'x')
    rc, rb, ra = [QQ(evaluate(v, tau)) for v in cover['lift']['residual_coefficients']]
    f0, f1, f2 = [R(v)(QQ(tau)) for v in cover['lift']['line_coefficients']]
    if not ra or not f2:
        return None, 'EXCEPTIONAL_RESIDUAL_OR_LINE_MAP'
    disc = rb*rb-4*ra*rc
    if not disc or not disc.is_square():
        return None, 'UNRESOLVED_TWO_DISTINCT_RATIONAL_BRANCHES'
    roots = sorted([(-rb-disc.sqrt())/(2*ra), (-rb+disc.sqrt())/(2*ra)])
    parent = read(f.ART/'curve302_recovered_mw17_parent_v1.json')
    old = EllipticCurve(QQ, [QQ(evaluate(v, tau)) for v in parent['a_invariants']])
    _, _, c, d = map(int, chart['parameter_matrix'])
    h = QQ(c*s.numerator+d*s.denominator); w = QQ(chart['weierstrass_u'])
    points = []
    for x in roots:
        y = -(f0+f1*x)/f2
        old([x, y])
        X = (x+old.b2()/12)*h**4/w**2
        Y = (y+(old.a1()*x+old.a3())/2)*h**6/w**3
        P = F(str(X)), F(str(Y))
        require(P[1]**2 == P[0]**3+F(row['model'][3])*P[0]+F(row['model'][4]), 'transported point is off curve')
        points.append(P)
    return sorted(points), None


def confirm(folder, case):
    p = checked_protocol(folder)
    require(p['factory']['chart_sha256'] == sha(folder/'parameter-chart.json'), 'parameter chart changed')
    chart = read(folder/'parameter-chart.json'); bindings(f.ROOT, chart['inputs'])
    address = next(r for r in factory.parameters(p['factory']['max_height'], p['factory']['maximum_parameters']) if r['id'] == case)
    seed = folder/'seeds'/case
    atomic(seed/'address.json', address, immutable=True)
    s = factory.image(chart, address['a'], address['b'])
    if s is None:
        atomic(seed/'result.json', dict(status='PARAMETER_POLE', address=address), immutable=True)
        return
    arithmetic = f.Arithmetic(p)
    record = arithmetic.record(address['index'], s.numerator, s.denominator)
    if record[3] != 'SMOOTH':
        atomic(seed/'result.json', dict(status='SINGULAR_FIBRE', address=address, parameter=str(s)), immutable=True)
        return
    if record[9] == 'CHART_POLE':
        atomic(seed/'result.json', dict(status='ORIGINAL_PARAMETER_POLE', address=address, parameter=str(s)), immutable=True)
        return
    require(record[9] == 'SPLIT', 'parametrized conic failed splitting identity')
    row = arithmetic.candidate(record, 'orbit8044_parametrized_factory')
    row.update(id=case, factory_u=address['u'], factory_u_height=address['height'])
    prepared = prepare(seed, row)
    if prepared is None:
        return
    model, base, proof = prepared
    points, exception = both_points(row)
    if exception:
        atomic(seed/'result.json', dict(status=exception, parameter=str(s), rank_lower_bound=17), immutable=True)
        return
    atomic(seed/'conic-points.json', [list(map(str, P)) for P in points], immutable=True)
    trace_word = list(map(int, read(f.COVER)['lift']['trace_word']))
    dependencies = factory.obvious_dependencies(model, base, points, trace_word)
    atomic(seed/'generic-dependence.json', dependencies, immutable=True)
    candidates = [P for P, d in zip(points, dependencies['points']) if d['status'] == 'UNKNOWN']
    found = first_seed(model, base, proof, candidates)
    if found:
        seal_seed(seed, row, model, base, *found,
            evidence=dict(kind='frozen_conic', sha256=sha(seed/'conic-points.json'),
                          factory_address_sha256=sha(seed/'address.json')))
    else:
        atomic(seed/'result.json', dict(status='EXACT_INHERITED_CONIC_POINTS' if not candidates else 'EXTRA_DIRECTION_UNRESOLVED',
            parameter=str(s), rank_lower_bound=17, quartic_seed_charts=0), immutable=True)
    print('ORBIT8044_SEED', address['u'], read(seed/'result.json')['status'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['chart', 'confirm'])
    p.add_argument('--directory', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--case')
    a = p.parse_args()
    if a.action == 'chart':
        require(a.output is not None, 'chart output required')
        build_chart(a.output.resolve())
    else:
        require(a.directory is not None and a.case is not None, 'directory and case required')
        confirm(a.directory.resolve(), a.case)

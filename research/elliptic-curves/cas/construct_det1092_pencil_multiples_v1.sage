#!/usr/bin/env sage-python
"""Frozen generic pencil operations n=2,3; no specialized point inputs.

Run each n separately under timeout 25s. Checkpoints are immutable; failures
or killed jobs are not negative mathematical results. No search is performed.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
GENERIC = ART / 'det1092_norm8_seed_cover_v2/generic.json'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
OUT = ART / 'det1092_pencil_multiples_v1'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')


def rec(value):
    return {key: list(map(str, getattr(value, key)().list()))
            for key in ['numerator', 'denominator']}


def run(n):
    started = time.monotonic()
    OUT.mkdir(exist_ok=True)
    protocol = {
        'classification': 'prospective generic-only construction',
        'multiples': [2, 3], 'pencil_orbit': 20124,
        'origin': 'generic inherited +P14',
        'B': 'conjugate inherited -P15 minus origin',
        'limits': {'seconds_per_process': 25, 'new_parameters': 0,
                   'point_searches': 0, 'exceptional_point_inputs': 0},
        'inputs': {str(p.relative_to(ROOT)): digest(p)
                   for p in [GENERIC, PARENT, Path(__file__)]},
    }
    protocol_path = OUT / 'protocol.json'
    if protocol_path.exists():
        assert json.loads(protocol_path.read_text()) == protocol
    else:
        save(protocol_path, protocol)
    destination = OUT / ('multiple%d.json' % n)
    assert not destination.exists()
    data = json.loads(GENERIC.read_text())
    parent = json.loads(PARENT.read_text())
    R = PolynomialRing(QQ, 'z'); z = R.gen(); K = R.fraction_field()
    def value(record):
        return K(R(record['numerator'])) / R(record['denominator'])
    Q = PolynomialRing(K, 't'); t = Q.gen()
    f = Q([R(row) for row in data['quartic_t_coefficients_in_z']])
    origin = data['sections'][0]
    t0 = value(origin['t_of_z']); s = value(origin['W_of_z'])
    q0, q1, q2, q3, q4 = f(t + t0).list()
    assert q0 == s*s and s
    J = EllipticCurve(K, [0, q2, 0, q1*q3-4*s*s*q4,
                          s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4])
    xb = q1*q1/(4*s*s)-q2
    yb = -(q1*xb+2*s*s*q3)/(2*s)
    B = J([xb, yb]); N = n*B
    print('n', n, 'group operation done', time.monotonic()-started, flush=True)
    assert not N.is_zero()
    xn, yn = N[:2]
    inverse_den = xn*xn-4*s*s*q4
    assert inverse_den
    u = (2*s*yn+q1*xn+2*s*s*q3)/inverse_den
    T = t0+u
    W = xn*u*u/(2*s)-s-q1*u/(2*s)
    assert W*W == f(T)
    degree = max(T.numerator().degree(), T.denominator().degree())
    print('n', n, 'quartic identity, map degree', degree,
          time.monotonic()-started, flush=True)
    # Evaluate the original family and the source-only inverse map.
    aa = [value(row) for row in parent['a_invariants']]
    def at(v):
        return K(v.numerator()(T))/v.denominator()(T)
    E = EllipticCurve(K, list(map(at, aa)))
    h = R(data['pole_h'])(T)
    cx = R(data['nx'])(T)/(h*h)
    cy = R(data['ny'])(T)/(h*h*h)
    shift = R(data['shift'])(T)
    m = h*z-shift/h
    X = (h*W-cx+m*m)/2
    Y = m*(X-cx)-cy
    x = X-E.b2()/12; y = Y-(E.a1()*x+E.a3())/2
    assert E([x, y])
    assert (h*(Y+cy)/(X-cx)+shift)/(h*h) == z
    print('n', n, 'original equation and birationality checked',
          time.monotonic()-started, flush=True)
    numerator, denominator = T.numerator(), T.denominator()
    ramification = numerator.derivative()*denominator-numerator*denominator.derivative()
    assert ramification
    # Clearing Delta(T)'s denominator introduces no zeros at affine finite T.
    delta = E.discriminant()
    bad = denominator*delta.numerator()
    ram_support = ramification.squarefree_part()
    smooth_support = ram_support // ram_support.gcd(bad)
    result = {
        'classification': 'new verified construction' if smooth_support.degree()>0
                          else 'bounded construction; ramification test inconclusive',
        'status': 'PASS_RANK18_RATIONAL_BASE_CHANGE' if smooth_support.degree()>0
                  else 'NO_SMOOTH_RAMIFICATION_CERTIFIED',
        'n': n, 'parameter_map_degree': int(degree),
        't_of_z': rec(T), 'W_of_z': rec(W),
        'original_point': [rec(x), rec(y)],
        'ramification_polynomial': list(map(str, ramification.list())),
        'smooth_ramification_support': list(map(str, smooth_support.list())),
        'ramification_support_degree': int(ram_support.degree()),
        'smooth_ramification_support_degree': int(smooth_support.degree()),
        'identities': ['W^2=F_z(T)', 'original elliptic equation',
                       'forward pencil parameter recovers z'],
        'boundary': 'A generic rank18 lower bound over Q(z), not an exact rank or a specialized rank certificate. No302 witness or later search artifacts read.',
        'seconds': time.monotonic()-started,
        'inputs': protocol['inputs'],
    }
    save(destination, result)
    print(result['status'], 'degree', degree, 'smooth ramification support',
          smooth_support.degree(), 'seconds', result['seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--multiple', type=int, choices=[2, 3], required=True)
    run(parser.parse_args().multiple)

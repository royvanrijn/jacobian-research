#!/usr/bin/env sage-python
"""Two frozen generic elliptic-pencil multiples, factored-map certificates.

Run each n separately under timeout25s. No point/parameter search. A bounded
prime list is used only to certify ramification of the already fixed maps.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
from sage.all import QQ, GF, PolynomialRing, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
GENERIC = ART / 'det1092_norm8_seed_cover_v2/generic.json'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
OUT = ART / 'det1092_pencil_multiples_v2'
PRIMES = [47, 53, 101, 107, 109, 127, 131, 137, 139, 149, 157, 163]


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
        'ramification_certificate_primes': PRIMES,
        'limits': {'seconds_per_process': 25, 'new_parameters': 0,
                   'point_searches': 0, 'exceptional_point_inputs': 0},
        'inputs': {str(p.relative_to(ROOT)): digest(p)
                   for p in [GENERIC, PARENT, Path(__file__)]},
    }
    path = OUT / 'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        save(path, protocol)
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
    assert not N.is_zero()
    xn, yn = N[:2]
    inverse_den = xn*xn-4*s*s*q4
    assert inverse_den
    u = (2*s*yn+q1*xn+2*s*s*q3)/inverse_den
    T = t0+u
    W = xn*u*u/(2*s)-s-q1*u/(2*s)
    assert W*W == f(T)
    degree = max(T.numerator().degree(), T.denominator().degree())
    print('n', n, 'quartic identity, degree', degree,
          time.monotonic()-started, flush=True)
    # The original rational map is retained as a composition, not expanded.
    # Its birational identities are checked in the independent verifier.
    E = EllipticCurve(K, [value(row) for row in parent['a_invariants']])
    delta = E.discriminant(); assert delta.denominator().degree()==0
    delta = delta.numerator()
    numerator, denominator = T.numerator(), T.denominator()
    ram = numerator.derivative()*denominator-numerator*denominator.derivative()
    assert ram
    sf = ram.squarefree_part().monic()
    attempts = []; certificate = None
    for p in PRIMES:
        try:
            Rp = PolynomialRing(GF(p), 'z')
            nn, dd, rr, ds = map(Rp, [numerator, denominator, sf, delta])
        except (ZeroDivisionError, ValueError):
            attempts.append({'prime': p, 'status': 'denominator_bad'}); continue
        if any(a.degree()!=b.degree() for a,b in
               [(nn,numerator),(dd,denominator),(rr,sf),(ds,delta)]):
            attempts.append({'prime': p, 'status': 'degree_drop'}); continue
        if nn.gcd(dd).degree() or rr.gcd(rr.derivative()).degree():
            attempts.append({'prime': p, 'status': 'map_or_ramification_bad'}); continue
        dd_power = ds.degree()
        bad = dd * sum((ds[i]*nn**i*dd**(dd_power-i)
                        for i in range(dd_power+1)), Rp.zero())
        common = rr.gcd(bad)
        lower = rr.degree()-common.degree()
        attempts.append({'prime': p, 'status': 'checked',
                         'smooth_ramification_support_lower_bound': int(lower)})
        if lower>0:
            certificate = {'prime': p, 'ramification_support_mod_p': list(map(int, rr.list())),
                'bad_support_gcd_mod_p': list(map(int, common.list())),
                'smooth_ramification_support_lower_bound': int(lower)}
            break
    result = {
        'classification': 'new verified construction' if certificate else 'bounded inconclusive check',
        'status': 'PASS_RANK18_RATIONAL_BASE_CHANGE' if certificate else 'UNKNOWN_RAMIFICATION',
        'n': n, 'parameter_map_degree': int(degree),
        't_of_z': rec(T), 'W_of_z': rec(W),
        'map_to_original': data['inverse_map'],
        'ramification_support': list(map(str, sf.list())),
        'ramification_certificate': certificate, 'prime_attempts': attempts,
        'ramification_support_degree': int(sf.degree()),
        'boundary': 'Generic rank18 lower bound over Q(z), not exact rank or a specialized certificate. No302 witness or later search artifact read. Independent replay required.',
        'seconds': time.monotonic()-started, 'inputs': protocol['inputs'],
    }
    save(destination, result)
    print(result['status'], 'degree', degree, 'ramification', certificate,
          'seconds', result['seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--multiple', type=int, choices=[2, 3], required=True)
    run(parser.parse_args().multiple)

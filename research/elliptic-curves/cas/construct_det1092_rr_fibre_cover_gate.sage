#!/usr/bin/env sage-python
"""One RR branch-divisor gate, with a frozen 11-fibre local comparison.

No point search, Selmer computation, class group, or pilot read/write.
The control roster is the already frozen, completed V1 roster, not V3.
"""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
PREFIX = 'det1092_rr_fibre_cover_gate'
PROTOCOL = ART / (PREFIX + '_protocol_v1.json')
OUT = ART / (PREFIX + '_v1.json')
HALVES = ART / 'det1092_rr_net_halving_gate_v1.json'
ROSTER = ROOT / 'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v1/roster.json'
TRANSPORT = ART / 'det1092_reduced_parameter_chart_v1/generic-proof.json'
PRIMES = [17, 47, 53, 61, 67, 71, 79, 83, 89, 101, 107, 113,
          127, 137, 149, 179, 191, 197]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write_new(path, data):
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def algebra():
    # Parametrize a generic point Q and its tangent; C=2Q.
    R = PolynomialRing(QQ, ['m', 'qx', 'qy', 'n', 'x'])
    m, qx, qy, n, x = R.gens()
    K = R.fraction_field()
    a = 2*m*qy - 3*qx*qx
    b = qy*qy - qx**3 - a*qx
    cx = m*m - 2*qx
    cy = m*(qx-cx) - qy
    assert cy*cy == cx**3 + a*cx + b
    H = lambda z: z**4 - 6*cx*z*z - 8*cy*z - 3*cx*cx - 4*a
    assert H(m) == 0
    e = K(qx) + 2*qy/(n-m)
    fe = e**3 + a*e + b
    assert (n-m)**4 * fe == qy*qy*H(n)
    # Exact addition, reduced modulo the equation of R=(x,y).
    S = PolynomialRing(K, 'y')
    y = S.gen()
    curve = y*y - (x**3+a*x+b)
    slope = (y-qy)/K(x-qx)
    px = slope*slope-x-qx
    py = -qy+slope*(qx-px)
    assert (py+cy-(m+K(2*qy)/(x-qx))*(px-cx)) % curve == 0
    # The two elliptic quotients of the resulting even sextic.
    Z = PolynomialRing(K, 'z')
    z = Z.gen()
    sextic = (z*z+e)**3+a*(z*z+e)+b
    assert sextic == z**6+3*e*z**4+(3*e*e+a)*z*z+fe
    # Second quotient: U=fe/z^2, V=fe*y/z^3.
    assert fe**2*sextic == fe**3+(3*e*e+a)*fe**2*z*z+3*e*fe**2*z**4+fe**2*z**6
    return {
        'half_coordinates': ['qx=(m^2-cx)/2', 'qy=m*(qx-cx)-cy'],
        'branch_parameter': 'e=qx+2*qy/(n-m)',
        'sextic': 'y^2=z^6+3*e*z^4+(3*e^2+a)*z^2+e^3+a*e+b',
        'map_to_marked_fibre': 'P=(z^2+e,y)+Q, where 2Q=C',
        'branch_identity': '(n-m)^4*(e^3+a*e+b)=qy^2*H(n)',
        'second_elliptic_factor': 'V^2=U^3+(3*e^2+a)*U^2+3*e*(e^3+a*e+b)*U+(e^3+a*e+b)^2',
        'second_quotient_map': 'U=(e^3+a*e+b)/z^2, V=(e^3+a*e+b)*y/z^3',
        'excluded': ['Delta_E=0', 'qy=0', 'n=m', 'H(n)=0', 'poles of displayed charts'],
        'built_in_points': 'Two rational infinity points over the half field, both mapping to Q; no independent direction is certified.'}


def freeze():
    if PROTOCOL.exists() or OUT.exists():
        raise FileExistsError('Preserve the frozen experiment')
    roster = read(ROSTER)
    assert len(roster) == 11
    a, b, c, d = map(QQ, read(TRANSPORT)['parameter_matrix'])
    cases = []
    for row in roster:
        s = QQ(row['parameter'])
        t = s if row['presentation'] == 'normalized' else (a*s+b)/(c*s+d)
        cases.append({'case': 'case-' + hashlib.sha256(str(t).encode()).hexdigest()[:12],
                      'parameter': str(t)})
    assert len({r['case'] for r in cases}) == 11
    data = {'classification': 'calibration-informed design; outcome-blind arithmetic execution',
            'cases': sorted(cases, key=lambda r: r['case']), 'primes': PRIMES,
            'rule': 'One fixed historical generic RR centre. For each of the unchanged V1 roster fibres, test the monic half quartic at every frozen prime. Record denominator exclusions and all root counts. Any integral no-root reduction certifies no local half and hence no genus2 double cover branched on one RR pair over Q. Never replace a case or extend the prime list.',
            'boundary': 'Not a genus2 Selmer panel. The branch-divisor construction must exist over Q before such descent is meaningful. No control rank or point-search result is an execution input.',
            'limits': {'cases': 11, 'primes_per_case': 18, 'point_searches': 0,
                       'Selmer_runs': 0, 'class_group_runs': 0, 'pilot_changes': 0},
            'inputs': {str(p.relative_to(ROOT)): sha(p) for p in
                       [HALVES, ROSTER, TRANSPORT, Path(__file__)]}}
    write_new(PROTOCOL, data)
    print('FROZEN', len(cases), 'cases;', len(PRIMES), 'primes per case', flush=True)


def run():
    p = read(PROTOCOL)
    assert p['primes'] == PRIMES
    for path, expected in p['inputs'].items():
        assert sha(ROOT/path) == expected
    formula = algebra()
    h = read(HALVES)
    R = PolynomialRing(QQ, 't')
    def ev(rec, t):
        return R(rec['numerator'])(t)/R(rec['denominator'])(t)
    rows = []
    for row in p['cases']:
        t = QQ(row['parameter'])
        aa, bb = ev(h['short_A'], t), ev(h['short_B'], t)
        cx, cy = [ev(v, t) for v in h['centre_short']]
        assert cy*cy == cx**3+aa*cx+bb
        delta = -16*(4*aa**3+27*bb**2)
        assert delta and cy
        coeff = [-3*cx*cx-4*aa, -8*cy, -6*cx, QQ(0), QQ(1)]
        assert coeff == [ev(v, t) for v in h['halving_quartic_coefficients']]
        trials = []
        for prime in PRIMES:
            if any(v.denominator() % prime == 0 for v in coeff):
                trials.append({'prime': prime, 'status': 'EXCLUDED_NONINTEGRAL_CHART'})
                continue
            S = PolynomialRing(GF(prime), 'm')
            fp = S(coeff)
            roots = [int(z) for z in GF(prime) if fp(z) == 0]
            trials.append({'prime': prime, 'coefficients': list(map(int, fp.list())),
                           'roots': roots, 'squarefree': bool(fp.discriminant()),
                           'status': 'NO_LOCAL_HALF' if not roots else 'ROOT_MOD_P_ONLY'})
        witnesses = [v['prime'] for v in trials if v['status'] == 'NO_LOCAL_HALF']
        rows.append({**row, 'short_curve': [str(aa), str(bb)], 'centre': [str(cx), str(cy)],
                     'halving_coefficients': list(map(str, coeff)), 'trials': trials,
                     'status': 'Q_BRANCH_COVER_OBSTRUCTED' if witnesses else 'UNKNOWN',
                     'obstruction_primes': witnesses})
        print(row['case'], rows[-1]['status'], witnesses, flush=True)
    # This irreducibility statement is checked only after all masked rows finish.
    at_zero = next(r for r in rows if QQ(r['parameter']) == 0)
    F17 = PolynomialRing(GF(17), 'm')
    f17 = F17(list(map(QQ, at_zero['halving_coefficients'])))
    assert f17.is_irreducible()
    data = {'classification': 'verified application and new deduction',
            'status': 'PASS_RR_BRANCH_DIVISIBILITY_GATE_AND_FROZEN_LOCAL_PANEL',
            'protocol_sha256': sha(PROTOCOL), 'construction_over_half_field': formula,
            'criterion': 'For a reduced RR pair D=P+(C-P), a smooth degree2 genus2 cover of E branched exactly at D exists over a field k iff C belongs to 2E(k).',
            'rows': rows,
            'zero_half_field': {'degree': 4, 'irreducibility_prime': 17,
                                'polynomial_mod_prime': list(map(int, f17.list())),
                                'no_extension_of_degree_below_four': True},
            'boundary': 'Only covers branched on one RR pair are excluded. Other branch divisors, higher-degree covers, or correspondences remain possible. No Selmer data or rank-incidence conclusion is asserted.',
            'limits': p['limits']}
    write_new(OUT, data)
    print(data['status'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['algebra', 'freeze', 'run'])
    mode = parser.parse_args().mode
    if mode == 'algebra':
        algebra()
        print('PASS_UNIVERSAL_BIELLIPTIC_IDENTITIES', flush=True)
    elif mode == 'freeze':
        freeze()
    else:
        run()

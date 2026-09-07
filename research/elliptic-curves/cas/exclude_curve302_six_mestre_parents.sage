#!/usr/bin/env sage-python
"""All-rational-parameter inverse test on six certified Mestre MW11 models.

The positive controls are their saved T=1 fibres. The negative certificate
is a primitive integral degree24 polynomial without a projective root at
one small prime. No height bound on T, new fibration or search sweep.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd, lcm, prime_range

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/mestre_468_replay_bundle_v1.json'
PUBLIC = ROOT/'elliptic-curves/cas/icarm_curve302.py'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/curve302_six_mestre_parent_inverse_v1.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def primitive(poly):
    poly *= lcm([c.denominator() for c in poly])
    poly /= gcd([ZZ(c) for c in poly])
    return -poly if poly.leading_coefficient() < 0 else poly


def compute():
    bundle = json.loads(INPUT.read_text())
    public = runpy.run_path(str(PUBLIC))
    E = EllipticCurve(QQ, [QQ(str(c)) for c in public['GENERAL_WEIERSTRASS_COEFFICIENTS']])
    target = E.j_invariant()
    ring = PolynomialRing(QQ, 'T')
    records = []
    for row in bundle['rows']:
        model = row['generic_heights']
        A = ring(model['curve_A_coefficients'])
        B = ring(model['curve_B_coefficients'])
        j = 6912*A**3/(4*A**3+27*B**2)
        numerator, denominator = j.numerator(), j.denominator()
        assert numerator.gcd(denominator) == 1
        degree = max(numerator.degree(), denominator.degree())
        assert degree == 24
        seed = row['specialized_seed']
        assert QQ(str(seed['fibre_T'])) == 1
        control_curve = EllipticCurve(QQ, list(map(QQ, seed['curve'])))
        control_j = control_curve.j_invariant()
        assert j(1) == j(-1) == control_j
        assert EllipticCurve(QQ, [A(1), B(1)]).is_isomorphic(control_curve)
        control = primitive(numerator*control_j.denominator()-denominator*control_j.numerator())
        assert control(1) == control(-1) == 0
        equation = primitive(numerator*target.denominator()-denominator*target.numerator())
        assert equation.degree() == degree
        witness = None
        for p in prime_range(5, 2000):
            fp = equation.change_ring(GF(p))
            if fp.degree() != degree or any(fp(t) == 0 for t in GF(p)):
                continue
            assert control.change_ring(GF(p))(1) == 0
            assert control.change_ring(GF(p))(-1) == 0
            witness = {'prime': int(p), 'coefficients_low_to_high': list(map(int, fp.list())),
                       'finite_root_count': 0, 'infinity_root': False,
                       'control_roots_plus_minus_one_retained': True}
            break
        assert witness is not None, 'No obstruction in declared prime list; leave unresolved'
        records.append({'outer_u': row['outer_u'], 'j_map_degree': int(degree),
                        'j_numerator_low_to_high': list(map(str, numerator.list())),
                        'j_denominator_low_to_high': list(map(str, denominator.list())),
                        'target_equation_low_to_high': list(map(str, equation.list())),
                        'control_j': str(control_j), 'control_parameters': ['1', '-1'],
                        'control_Q_isomorphism_at_T1': True, 'obstruction': witness})
        print('MESTRE_INVERSE', row['outer_u'], witness['prime'], flush=True)
    assert list(map(lambda r: int(r['outer_u']), records)) == [11, 13, 17, 19, 23, 29]
    return {'schema': 'curve302.six-mestre-parent-inverse.v1',
            'status': 'SIX_FIXED_MW11_FIBRATIONS_EXCLUDE302_AT_ALL_RATIONAL_PARAMETERS',
            'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), INPUT, PUBLIC]},
            'target_j': str(target), 'records': records,
            'proof': 'A rational parameter giving302 must annihilate the primitive integral homogeneous degree24 j-comparison polynomial. Reduction of a primitive parameter pair gives a projective root at every prime. Each saved reduction has no finite root and nonzero leading coefficient, hence no infinity root.',
            'boundary': 'Only these six displayed MW11 fibrations. No exclusion of their other elliptic fibrations, the determinant468 lattice type, other outer parameters, or all302 parents. No generic-rank14–16 fibration is constructed by this test.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    text = json.dumps(compute(), indent=2, sort_keys=True)+'\n'
    if args.check:
        assert OUT.read_text() == text
    else:
        assert not OUT.exists(), 'Use --check; preserve the certificate'
        OUT.write_text(text)
    print('PASS_SIX_MESTRE_INVERSE', flush=True)

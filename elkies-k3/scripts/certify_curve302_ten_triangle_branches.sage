#!/usr/bin/env sage-python
"""Ten pairwise different branch divisors on the fixed target j-line."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import GF,PolynomialRing

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
SCREEN = ART/'elkies-k3-curve302-shortword-triangles-v1.json'
HELPER = ROOT/'elkies-k3/scripts/certify_curve302_triangle_branch_separation.sage'
OUT = ART/'elkies-k3-curve302-ten-triangle-branches-v1.json'


def build():
    screen = json.loads(SCREEN.read_text())
    helper = runpy.run_path(str(HELPER)); R = PolynomialRing(GF(1013),'s'); s = R.gen()
    records = []; paths = [Path(__file__),SCREEN,HELPER]
    for row in screen['records']:
        if 'control_source' in row:
            path = ROOT/row['control_source']; paths.append(path)
            assert sha256(path.read_bytes()).hexdigest() == row['control_sha256']
            d = json.loads(path.read_text())
            n = R(d['j_numerator_coefficients_low_to_high']); den = R(d['j_denominator_coefficients_low_to_high'])
        else:
            d = next(d for d in row['results'] if d['prime'] == 1013)
            n = R(d['j_numerator']); den = R(d['j_denominator'])
        assert n.gcd(den) == 1 and max(n.degree(),den.degree()) == 24
        h,samples = helper['branch_polynomial'](n,den)
        residual,remainder = h.quo_rem(s**16*(s-R.base_ring()(1728))**12)
        assert not remainder and h.degree() == 43 and residual.degree() == 15
        records.append({'index':row['index'],'normalized_branch_polynomial':list(map(int,h.list())),
                        'residual_branch_polynomial':list(map(int,residual.list())),
                        'binary_discriminant_samples':samples})
    assert len(records) == len({tuple(r['normalized_branch_polynomial']) for r in records}) == 10
    return {'schema':'curve302.ten-triangle-branches.v1','status':'TEN_PAIRWISE_INEQUIVALENT_TRIANGLE_J_MAPS',
            'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
            'prime':1013,'records':records,
            'boundary':'Pairwise inequivalence of the ten elliptic fibrations follows from different branch divisors of their degree-preserving j-map reductions. All are on the same determinant948 K3; this is not ten surface classes or a302 parent.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--check',action='store_true'); args = parser.parse_args()
    signal.alarm(120)
    result = build()
    if args.check: assert result == json.loads(OUT.read_text())
    else:
        assert not OUT.exists(); OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'],flush=True)

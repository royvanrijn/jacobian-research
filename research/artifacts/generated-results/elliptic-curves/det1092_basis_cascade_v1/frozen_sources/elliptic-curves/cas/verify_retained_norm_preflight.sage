#!/usr/bin/env python3
"""Independent resultant norms and product-gcd replay, without gate imports."""
from sage.all import QQ, ZZ, PolynomialRing, prod, gcd, lcm
from pathlib import Path
from hashlib import sha256
import json
import argparse

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = ART/'retained_norm_preflight_inputs_v1.json'
OUTPUT = ART/'retained_norm_preflight_v1.json'
REPLAY = ART/'retained_norm_preflight_sage_v1.json'


def calculate():
    inputs = json.loads(INPUT.read_text())['cases']
    results = json.loads(OUTPUT.read_text())['rows']
    assert len(inputs) == len(results) == 12
    ring = PolynomialRing(QQ, 'z')
    checked = rejected = 0
    for case, result in zip(inputs, results):
        assert case['id'] == result['id']
        f = ring(case['cubic_ascending'])
        disc = f.discriminant()
        assert f.degree() == 3 and f.is_monic() and disc
        elements = [ring(v) for v in case['elements']]
        norms = [f.resultant(a) for a in elements]
        assert all(norms) and list(map(str, norms)) == result['norms']
        forbidden = 2*abs(disc.numerator())*disc.denominator()
        forbidden *= prod(x.denominator() for x in f.list())
        for alpha, n in zip(elements, norms):
            den = lcm(x.denominator() for x in alpha.list())
            content = gcd([ZZ(x*den) for x in alpha.list()])
            forbidden *= den*abs(content)*n.denominator()
        active = list(range(len(elements)))
        for wave in result['rounds']:
            assert wave and len({w['index'] for w in wave}) == len(wave)
            for w in wave:
                i = w['index']
                assert i in active
                avoid = forbidden*prod(abs(norms[j].numerator()) for j in active if j != i)
                remaining = abs(norms[i].numerator())
                while gcd(remaining, avoid) > 1:
                    remaining //= gcd(remaining, avoid)
                assert remaining == ZZ(w['remainder']) and gcd(remaining, avoid) == 1
                root = ZZ(w['floor_square_root'])
                assert root >= 0 and root**2 < remaining < (root+1)**2
                rejected += 1
            active = [i for i in active if i not in {w['index'] for w in wave}]
        assert active == result['unresolved_indices']
        assert result['input_generator_count'] == len(elements)
        assert result['forced_zero_count'] == len(elements)-len(active)
        assert result['unramified_coefficient_dimension_upper_bound'] == len(active)
        assert result['additional_strict_classes'] == ('UNKNOWN' if active else 0)
        assert result['whole_curve_rank_decision'] == 'UNKNOWN'
        checked += len(elements)
    assert checked == 428
    paths = [INPUT, OUTPUT, Path(__file__).resolve()]
    return {'status': 'PASS', 'schema': 'elliptic-curves.retained-norm-preflight-sage.v1',
            'bindings': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
            'independent_resultant_norms': checked, 'independent_isolation_witnesses': rejected,
            'boundary': 'Necessary coefficient restrictions; no class admissibility or elliptic rank exclusion.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.check:
        assert result == json.loads(REPLAY.read_text())
    else:
        with REPLAY.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write('\n')
    print(result['status'], result['independent_resultant_norms'], result['independent_isolation_witnesses'])

#!/usr/bin/env python3
"""Independent Sage arithmetic check of compression and common-field certificates."""
import argparse
import gzip
import json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, matrix, GF
import retrospective as r

HERE = Path(__file__).resolve().parent
GEOMETRY = r.OUT / 'rank_jump_solubility_first_geometry_v1.json.gz'
INPUT = r.OUT / 'rank_jump_soluble_quartet_compression_inputs_v1.json'
RESULT = r.OUT / 'rank_jump_soluble_quartet_compression_v1.json'
ATLAS = r.OUT / 'rank_jump_atlas_common_cover_v1.json'
OUTPUT = r.OUT / 'rank_jump_quartet_compression_verification_v1.json'


def verify():
    inp = r.read(INPUT); result = r.read(RESULT); atlas = r.read(ATLAS)
    for data in (inp, result, atlas):
        for path, sha in data['bindings'].items(): assert r.digest((r.ROOT / path).read_bytes()) == sha
    raw = json.loads(gzip.decompress(GEOMETRY.read_bytes()))
    # Divide by the leading coefficient rather than by integer content.
    monic = set()
    for c in raw:
        a, b, d = map(ZZ, c['integer_quadratic'])
        assert d and b*b != 4*a*d
        monic.add((QQ(a)/d, QQ(b)/d))
    assert len(monic) == len(raw) == atlas['distinct_branch_polynomials'] == atlas['distinct_quadratic_function_fields'] == 39119
    assert atlas['same_branch_groups'] == atlas['same_field_groups'] == []
    R = PolynomialRing(QQ, 't'); checks = 0
    for case, row in zip(inp['cases'], result['rows'], strict=True):
        forms = [R(c['form']) for c in case['covers']]
        assert row['id'] == case['id']
        assert list(map(str, (f.discriminant() for f in forms))) == row['quadratic_discriminants']
        assert all(f.gcd(g).degree() == 0 for i, f in enumerate(forms) for g in forms[i+1:])
        for q in row['character_quotients']:
            f = R(1)
            for i in range(4):
                if q['mask'] >> i & 1: f *= forms[i]
            assert f.is_squarefree() and list(map(str, f.list())) == q['coefficients']
            assert (f.degree()-2)//2 == q['genus']
        assert sum(q['genus'] for q in row['character_quotients']) == row['carrier_genus'] == 17
        for f, w in zip(forms, row['individual_conic_witnesses'], strict=True):
            assert f(QQ(w['parameter'])) == QQ(w['root'])**2
        for fibre in row['fibres']:
            t = QQ(fibre['published_parameter']); a, b = t.numerator(), t.denominator()
            numerators = [ZZ(f[0]*b*b+f[1]*a*b+f[2]*a*a) for f in forms]
            assert all(numerators)
            kernel = [0]
            for mask in range(1, 16):
                n = ZZ(1)
                for i in range(4):
                    if mask >> i & 1: n *= numerators[i]
                if n.is_square(): kernel.append(mask)
                checks += 1
            assert kernel == fibre['square_character_masks']
            rank = matrix(GF(2), [[(m >> i) & 1 for i in range(4)] for m in kernel]).rank()
            assert 4-rank == fibre['radicand_squareclass_rank']
            assert fibre['full_lift'] == all(n.is_square() for n in numerators)
            assert fibre['all_four_product_square'] == (15 in kernel)
    assert checks == result['product_square_tests'] == 1440
    return {'schema': 'rank-jump.quartet-compression-verification.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (INPUT, RESULT, ATLAS, GEOMETRY, Path(__file__), HERE/'retrospective.py')},
            'distinct_monic_branch_polynomials': len(monic), 'independent_integer_square_tests': checks,
            'boundary': 'Independent exact polynomial and arithmetic checks. Riemann-Hurwitz, Castelnuovo-Severi and the finiteness theorem are mathematical dependencies explained in the note, not reimplemented algorithms.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['build', 'check']); mode = p.parse_args().mode
    result = verify()
    if mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS: 39119 distinct quadratic fields and 1440 independent square tests')

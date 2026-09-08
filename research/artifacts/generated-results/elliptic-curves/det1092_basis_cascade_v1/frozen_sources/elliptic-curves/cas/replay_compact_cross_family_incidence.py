#!/usr/bin/env python3
"""Sage-free verification of projective j-map exclusions and all rational preimages."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from fixed_cubic_geometry import poly, add, mul, power, scale, value
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint
ROOT = Path(__file__).resolve().parents[2]

def mod_trim(a, p):
    a = [x % p for x in a]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def coprime_mod(a, b, p):
    a, b = mod_trim(a, p), mod_trim(b, p)
    while b != [0]:
        while a != [0] and len(a) >= len(b):
            shift = len(a)-len(b); coefficient = a[-1]*pow(b[-1], -1, p) % p
            for i, x in enumerate(b):
                a[shift+i] = (a[shift+i]-coefficient*x) % p
            a = mod_trim(a, p)
        a, b = b, a
    return len(a) == 1 and a != [0]

def mod_value(a, x, p):
    answer = 0
    for y in reversed(a):
        answer = (answer*x+y) % p
    return answer

def replay(path, output):
    if output.exists():
        raise FileExistsError('preserve incidence replay')
    data = cert.read(path)
    for name, h in data['sources'].items():
        if cert.hashed(ROOT/name) != h:
            raise ArithmeticError('frozen input/source changed')
    if data['status'] != 'COMPLETE_DECLARED_INCIDENCE_AUDIT':
        raise ArithmeticError('incidence audit incomplete')
    original = cert.read(ROOT/'artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v1.json')
    if data['targets'] != [{k:r[k] for k in ('id','curve','j_invariant','family','parameter','rank_lower_bound')} for r in original['curves']]:
        raise ArithmeticError('new-curve roster changed')
    maps = {row['family']:row for row in data['maps']}
    targets = {row['id']:row for row in data['targets']}
    if len(maps) != 12 or len(targets) != 32 or len(data['pairs']) != 384 or {(r['target'],r['family']) for r in data['pairs']} != {(i,f) for i in targets for f in maps}:
        raise ArithmeticError('incomplete incidence product')
    for target in targets.values():
        inv = cert.weierstrass_invariants(target['curve'])
        if inv['c4']**3/inv['discriminant'] != cert.F(target['j_invariant']):
            raise ArithmeticError('target j differs from the curve')
    for family in maps.values():
        source = cert.read(ROOT/family['source'])
        row = next(r for r in source['families'] if r.get('family',r.get('fibration_id')) == family['family']) if 'families' in source else source
        if family['A'] != row['A_coefficients_low_to_high'] or family['B'] != row['B_coefficients_low_to_high']:
            raise ArithmeticError('family coefficients differ')
        a, b = poly(family['A']), poly(family['B'])
        n, d = family['numerator'], family['denominator']
        if len(n) != family['degree']+1 or len(d) != len(n):
            raise ArithmeticError('homogeneous degree changed')
        if mul(poly(n), add(scale(power(a,3),4),scale(power(b,2),27))) != mul(poly(d),scale(power(a,3),6912)):
            raise ArithmeticError('rational j-map identity failed')
        if not family['good_map_reductions']:
            raise ArithmeticError('no coprimality witness')
        for row in family['good_map_reductions']:
            p = row['prime']
            if not _is_prime(p) or not coprime_mod(n,d,p) or (n[-1] % p == 0 and d[-1] % p == 0):
                raise ArithmeticError('map is not a morphism of unchanged projective degree')
            image = set()
            for x in range(p):
                u,v = mod_value(n,x,p),mod_value(d,x,p)
                image.add(None if v == 0 else u*pow(v,-1,p) % p)
            image.add(None if d[-1] % p == 0 else n[-1]*pow(d[-1],-1,p) % p)
            if sorted(x for x in image if x is not None) != row['finite_images'] or (None in image) != row['infinity_in_image']:
                raise ArithmeticError('projective image differs')
    for pair in data['pairs']:
        family, target = maps[pair['family']], targets[pair['target']]
        j = cert.F(target['j_invariant'])
        if pair['status'] == 'NO_RATIONAL_J_PREIMAGE':
            witness = pair['modular_exclusion']; p = witness['prime']
            row = next(r for r in family['good_map_reductions'] if r['prime'] == p)
            v = None if j.denominator % p == 0 else j.numerator*pow(j.denominator,-1,p) % p
            if v != witness['target_reduction'] or (v is None and row['infinity_in_image']) or (v is not None and v in row['finite_images']):
                raise ArithmeticError('target is not excluded')
        elif pair['status'] == 'RATIONAL_J_PREIMAGES_CERTIFIED':
            equation = add(scale(poly(family['numerator']),j.denominator),scale(poly(family['denominator']),-j.numerator))
            product = poly([1])
            for root in pair['rational_roots']:
                t = cert.F(root['parameter'])
                if value(equation,t) or not value(poly(family['denominator']),t):
                    raise ArithmeticError('false rational preimage')
                product = mul(product,power(poly([-t,1]),root['multiplicity']))
            for factor in pair['residual_factors']:
                f = factor['coefficients']; p = factor['no_projective_root_prime']
                if not _is_prime(p) or not f[-1] % p or any(mod_value(f,x,p) == 0 for x in range(p)):
                    raise ArithmeticError('residual rational-root exclusion failed')
                product = mul(product,power(poly(f),factor['multiplicity']))
            if scale(product,1/product[-1]) != scale(equation,1/equation[-1]) or pair['infinity'] != (len(equation)-1 < family['degree']):
                raise ArithmeticError('rational preimage factorization incomplete')
        else:
            raise ArithmeticError('unresolved preimage claim')
    counts = dict(Counter(r['status'] for r in data['pairs']))
    paths = [Path(__file__).resolve(), path, ROOT/'elliptic-curves/cas/fixed_cubic_geometry.py']
    checkpoint(output, {'schema':'elliptic-curves.compact-cross-family-j-incidence-replay.v1',
        'status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'pairs_checked':384,'status_counts':counts,
        'claim_boundary':'Exact j-incidence only in the recorded twelve presentations, including infinity; no rational isomorphism or section independence claim.'})
    print('REPLAYED ALL CROSS-FAMILY J INCIDENCE', counts, flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args(); replay(a.input.resolve(),a.output.resolve())

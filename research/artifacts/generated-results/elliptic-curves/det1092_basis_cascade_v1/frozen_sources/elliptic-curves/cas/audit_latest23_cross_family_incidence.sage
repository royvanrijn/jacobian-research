#!/usr/bin/env sage-python
"""Exact incidence extension for the twenty-three curves added after the47-curve incidence audit."""
import argparse
from pathlib import Path
import sys
from sage.all import QQ, ZZ, GF, PolynomialRing, prime_range, gcd, lcm
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART = ROOT/'artifacts/generated-results/elliptic-curves'
R = PolynomialRing(QQ, 't')

def primitive_pair(n, d):
    common = n.gcd(d)
    n, d = n//common, d//common
    scale = lcm([q.denominator() for q in list(n)+list(d)])
    nn, dd = n*scale, d*scale
    content = gcd([ZZ(q) for q in list(nn)+list(dd)])
    nn, dd = nn/content, dd/content
    degree = max(nn.degree(), dd.degree())
    return [int(nn[i]) for i in range(degree+1)], [int(dd[i]) for i in range(degree+1)]

def finite_images(n, d, p):
    F = GF(p); P = PolynomialRing(F, 't'); a, b = P(n), P(d)
    if a.gcd(b).degree() > 0 or (n[-1] % p == 0 and d[-1] % p == 0):
        return None
    images = set()
    for x in F:
        y, z = a(x), b(x)
        images.add(None if not z else int(y/z))
    images.add(None if d[-1] % p == 0 else n[-1]*pow(d[-1], -1, p) % p)
    return sorted(x for x in images if x is not None), None in images

def no_root_prime(poly):
    for p in prime_range(3, 1000):
        p = int(p)
        if ZZ(poly.leading_coefficient()) % p == 0:
            continue
        q = PolynomialRing(GF(p), 't')(poly)
        if all(q(x) for x in GF(p)):
            return p
    return None

def run(output):
    if output.exists():
        raise FileExistsError('preserve the incidence audit')
    source_paths = [ART/'compact_six_r17_atlas_v1.json', ART/'compact_five_mw16_atlas_v1.json',
                    ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json',
                    ART/'new_high_rank_curve_index_v6.json']
    families = []
    for path in source_paths[:2]:
        for row in cert.read(path)['families']:
            families.append((row.get('family', row.get('fibration_id')), row, path))
    families.append(('published-R17', cert.read(source_paths[2]), source_paths[2]))
    targets = [{k:r[k] for k in ('id', 'curve', 'j_invariant', 'family', 'parameter', 'rank_lower_bound')}
               for r in cert.read(source_paths[3])['curves'] if int(r['id'].split('-')[-1])>47]
    if len(targets)!=23:raise ArithmeticError('latest-twenty-three roster changed')
    maps = []
    for name, row, path in families:
        a = R([QQ(q) for q in row['A_coefficients_low_to_high']])
        b = R([QQ(q) for q in row['B_coefficients_low_to_high']])
        n, d = primitive_pair(6912*a**3, 4*a**3+27*b**2)
        witnesses = []
        for p in prime_range(3, 252):
            p = int(p); image = finite_images(n, d, p)
            if image is not None:
                witnesses.append({'prime':p, 'finite_images':image[0], 'infinity_in_image':image[1]})
        maps.append({'family':name, 'source':str(path.relative_to(ROOT)),
            'A':row['A_coefficients_low_to_high'], 'B':row['B_coefficients_low_to_high'],
            'numerator':n, 'denominator':d, 'degree':len(n)-1, 'good_map_reductions':witnesses})
    result = {'schema':'elliptic-curves.compact-cross-family-j-incidence.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in source_paths+[Path(__file__).resolve()]},
        'targets':targets, 'maps':maps, 'pairs':[], 'status':'RUNNING',
        'claim_boundary':'Exact rational j-incidence within twelve recorded family presentations. Modular exclusions require a morphism of unchanged degree on projective space. Equal j does not assert rational isomorphism or independent sections. No exclusion of other families or additional points.'}
    checkpoint(output, result)
    for target in targets:
        j = QQ(target['j_invariant'])
        for family in maps:
            witness = None
            for row in family['good_map_reductions']:
                p = row['prime']
                value = None if j.denominator() % p == 0 else int(j.numerator())*pow(int(j.denominator()), -1, p) % p
                if (value is None and not row['infinity_in_image']) or (value is not None and value not in row['finite_images']):
                    witness = {'prime':p, 'target_reduction':value}; break
            pair = {'target':target['id'], 'family':family['family']}
            if witness is not None:
                pair.update(status='NO_RATIONAL_J_PREIMAGE', modular_exclusion=witness)
            else:
                poly = j.denominator()*R(family['numerator'])-j.numerator()*R(family['denominator'])
                if not poly:
                    raise ArithmeticError('unexpected constant j map')
                roots, residual = [], []
                factors = list(poly.factor())
                for factor, multiplicity in factors:
                    if factor.degree() == 1:
                        roots.append({'parameter':str(-factor[0]/factor[1]), 'multiplicity':int(multiplicity)})
                    else:
                        integer = factor*factor.denominator()
                        residual.append({'coefficients':[int(q) for q in integer], 'multiplicity':int(multiplicity),
                                         'no_projective_root_prime':no_root_prime(integer)})
                infinity = poly.degree() < family['degree']
                complete = all(r['no_projective_root_prime'] is not None for r in residual)
                pair.update(status='RATIONAL_J_PREIMAGES_CERTIFIED' if complete else 'RESIDUAL_ROOT_COMPLETENESS_UNKNOWN',
                    rational_roots=roots, infinity=infinity, residual_factors=residual)
            result['pairs'].append(pair)
        checkpoint(output, result)
        print('CROSS-FAMILY INCIDENCE', target['id'], [(r['family'],r.get('rational_roots'),r.get('infinity')) for r in result['pairs'][-12:] if 'rational_roots' in r], flush=True)
    result['status'] = 'COMPLETE_DECLARED_INCIDENCE_AUDIT'
    checkpoint(output, result)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())

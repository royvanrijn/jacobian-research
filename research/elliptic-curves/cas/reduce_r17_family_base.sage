#!/usr/bin/env sage-python
"""Bounded Cremona--Stoll reduction of the base coordinate of compiled R17 models.

This reduces the auxiliary genus-three curve y^2=A(t) to propose an SL2(Z)
base change. Its rank is irrelevant. Both elliptic coefficient identities
are exact; no search outcomes or known high-rank parameters select the map.
"""
import argparse
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ, ZZ, pari

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import audit_r17_constant_scaling as scaling
from research_runtime.store import checkpoint

INPUT = ROOT / 'artifacts/generated-results/elliptic-curves/r17_constant_scaling_audit_v1.json'
DIRECTORY = ROOT / 'artifacts/local/elliptic-curves/r17-base-reduction-v1'


def sources():
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in
            (Path(__file__).resolve(), INPUT, Path(scaling.__file__).resolve())}


def prepare(directory):
    if (directory / 'protocol.json').exists():
        raise FileExistsError('base-reduction protocol already frozen')
    checkpoint(directory / 'protocol.json', {'schema': 'elliptic-curves.r17-base-reduction.v1',
        'sources': sources(), 'families': ['103b2', '11952', '074d9', '07ca9', '08234', '08f72'],
        'method': 'PARI Cremona-Stoll hyperellred of primitive A polynomial; exact SL2 base substitution on weights 8 and 12; bounded constant-content scaling',
        'precision_bits': 8192, 'pari_stack_bytes': 256000000,
        'worker_wall_seconds': 120, 'worker_rss_bytes': 1073741824, 'maximum_concurrent_workers': 2,
        'gate': 'Start with 103b2, the smallest literal coefficient-size model. Continue to other families only if its checked reduction lowers the largest coefficient bit size by at least 25 percent.',
        'scope': 'At most six auxiliary polynomial reductions. No parameter sweep, point search, generic-rank increase, minimality or novelty claim.'})


def homogeneous(poly, weight, a, b, c, d):
    t = poly.parent().gen()
    return sum(poly[i] * (a*t+b)**i * (c*t+d)**(weight-i) for i in range(poly.degree()+1))


def run(directory, family):
    protocol = cert.read(directory / 'protocol.json')
    if protocol['sources'] != sources() or family not in protocol['families']:
        raise ArithmeticError('frozen source or family changed')
    output = directory / (family + '.json')
    if output.exists():
        raise FileExistsError('preserve the existing reduction attempt')
    row = next(r for r in cert.read(INPUT)['rows'] if r['family'] == family)
    R = PolynomialRing(QQ, 't'); t = R.gen()
    A, B = (R(list(map(QQ, row[k]))) for k in ('A_coefficients_low_to_high', 'B_coefficients_low_to_high'))
    primitive = A / QQ(scaling.content(list(map(cert.F, row['A_coefficients_low_to_high']))))
    if primitive.gcd(primitive.derivative()).degree() != 0 or any(q not in ZZ for q in primitive):
        raise ArithmeticError('auxiliary polynomial is not integral squarefree')
    pari.allocatemem(protocol['pari_stack_bytes'], silent=True)
    pari.set_real_precision_bits(protocol['precision_bits'])
    checkpoint(output, {'status': 'RUNNING', 'family': family, 'protocol_sha256': cert.hashed(directory / 'protocol.json')})
    reducer = pari('(P)->{my(m);my(Q=hyperellred(P,&m));[Q,m]}')
    reduced, transform = reducer(pari(primitive))
    matrix = transform[1]
    a, b, c, d = (QQ(matrix[i,j]) for i,j in ((0,0),(0,1),(1,0),(1,1)))
    if a*d-b*c != 1:
        raise ArithmeticError('base matrix is not SL2')
    AA, BB = homogeneous(A,8,a,b,c,d), homogeneous(B,12,a,b,c,d)
    aa, bb = [cert.F(str(x)) for x in AA.list()], [cert.F(str(x)) for x in BB.list()]
    u, diag = scaling.scale_for(aa,bb); new_A, new_B = AA/QQ(str(u))**4, BB/QQ(str(u))**6
    # Reverse substitution verifies the entire coefficient identity, not samples.
    if homogeneous(new_A,8,d,-b,-c,a)*QQ(str(u))**4 != A or homogeneous(new_B,12,d,-b,-c,a)*QQ(str(u))**6 != B:
        raise ArithmeticError('inverse weighted base identity failed')
    new_bits = scaling.bits([cert.F(str(x)) for x in new_A.list()+new_B.list()])
    result = {'status': 'PASS_EXACT_BASE_CHANGE', 'family': family, 'protocol_sha256': cert.hashed(directory/'protocol.json'),
        'constant_scaled_input_bits': row['after_bits'], 'after_bits': new_bits,
        'base_matrix_a_b_c_d': list(map(str,(a,b,c,d))), 'additional_scale_u': str(u),
        'total_scale_from_literal_source': str(u*cert.F(row['scale_u'])),
        'A_coefficients_low_to_high': list(map(str,new_A.list())), 'B_coefficients_low_to_high': list(map(str,new_B.list())),
        'auxiliary_reduction': str(reduced), 'auxiliary_transform': str(transform), 'constant_scaling': diag,
        'source': row['source'], 'source_sha256': row['source_sha256'],
        'claim_boundary': 'Exact base automorphism and constant Weierstrass isomorphism only. Section transport must be checked before prospective use.'}
    checkpoint(output, result)
    print('BASE REDUCTION',family,row['after_bits'],'->',new_bits,flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default='103b2')
    args=p.parse_args()
    prepare(args.directory) if args.stage=='prepare' else run(args.directory,args.family)

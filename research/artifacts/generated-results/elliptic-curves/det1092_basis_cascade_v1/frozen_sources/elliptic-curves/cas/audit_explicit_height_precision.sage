#!/usr/bin/env sage-python
"""Measure the effect of explicitly requesting PARI height precision on fixed curves."""
import argparse
import inspect
from pathlib import Path
import sys
from sage.all import QQ, EllipticCurve, pari, RealField
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
INDEX = ROOT/'artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v1.json'

def run(output):
    if output.exists():
        raise FileExistsError('preserve height precision audit')
    selected = [r for r in cert.read(INDEX)['curves'] if r['rank_lower_bound'] == 25]
    if len(selected) != 3:
        raise ArithmeticError('fixed rank25 roster changed')
    RF = RealField(512); rows = []
    pari.default('realprecision',110)
    for row in selected:
        model = tuple(map(cert.F,row['curve']))
        points = [tuple(map(cert.F,p)) for p in row['points']]
        proof = row['rank_certificate']
        cert.checked_rank(model,points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
        E = pari(EllipticCurve(QQ,[QQ(str(q)) for q in model]))
        pairs = [[QQ(str(q)) for q in point] for point in points]
        default = E.ellheightmatrix(pairs)
        explicit64 = E.ellheightmatrix(pairs,precision=64)
        high = E.ellheightmatrix(pairs,precision=384)
        old = [[str(default[i,j]) for j in range(25)] for i in range(25)]
        new = [[str(high[i,j]) for j in range(25)] for i in range(25)]
        if any(default[i,j] != explicit64[i,j] for i in range(25) for j in range(25)):
            raise ArithmeticError('default64 API observation changed')
        changes = [[i,j] for i in range(25) for j in range(25) if round(cert.F(old[i][j])*1000000) != round(cert.F(new[i][j])*1000000)]
        error = max(abs(RF(old[i][j])-RF(new[i][j])) for i in range(25) for j in range(25))
        rows.append({'id':row['id'],'family':row['family'],'parameter':row['parameter'],
            'curve':row['curve'],'points':row['points'],'default_gram':old,'explicit384_gram':new,
            'maximum_observed_entry_difference':str(error),'rounded_million_scale_changes':changes,
            'method_signature':str(inspect.signature(E.ellheightmatrix))})
        print('HEIGHT PRECISION',row['parameter'],'rounded entries changed',len(changes),'max difference',error,flush=True)
    paths = [Path(__file__).resolve(),INDEX,ROOT/'elliptic-curves/cas/prospective_half_lattice.sage',ROOT/'elliptic-curves/cas/prospective_half_lattice_v2.sage']
    checkpoint(output,{'schema':'elliptic-curves.explicit-height-precision-audit.v1',
        'status':'COMPLETE_FIXED_NUMERICAL_COMPARISON','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'curves':rows,'claim_boundary':'Measured default64 versus explicit384-bit PARI Gram entries on three exactly certified point sets. These are numerical comparisons, not rigorous height intervals or rank proofs. Fixed rounded-Gram equality applies only to these inputs.'})

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())

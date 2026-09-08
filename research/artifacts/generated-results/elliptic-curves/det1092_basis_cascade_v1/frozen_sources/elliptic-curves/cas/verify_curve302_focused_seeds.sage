#!/usr/bin/env sage-python
"""Copied-input seed membership, generic specialization and finite independence."""
import json,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,EllipticCurve,PolynomialRing

def main(d):
    read=lambda n:json.loads((d/n).read_bytes())
    rank=SourceFileLoader('standalone_rank',str(d/'verify_factor_free_rank.sage')).load_module()
    intake=read('intake.json');rows=[]
    for row in intake['rows']:
        ident=row['id'];seed=read(ident+'-seed.json');cloud=read(ident+'-cloud.json')
        assert seed['curve']==cloud['curve'] and seed['points']==cloud['points']
        assert len(seed['points'])==row['initial_rank']
        rows.append(rank.check_one(d/(ident+'-cloud.json')))
    parent=read('generic-parent.json');R=PolynomialRing(QQ,'t')
    value=lambda v:R(v['numerator'])(0)/R(v['denominator'])(0)
    a=[value(v) for v in parent['a_invariants']];E=EllipticCurve(QQ,a)
    raw=[E([value(v) for v in P]) for P in parent['basis_weierstrass_coordinates']]
    seed=read('curve302-generic17-seed.json');S=EllipticCurve(QQ,seed['curve'])
    assert S.a_invariants()==(0,0,0,-27*E.c4(),-54*E.c6())
    assert [[str(36*P[0]+3*E.b2()),str(108*(2*P[1]+E.a1()*P[0]+E.a3()))] for P in raw]==seed['points']
    assert read('curve302-full31-seed.json')['curve']==seed['curve']
    result=dict(status='PASS',rows=rows,scope='Exact seed points, generic17 specialization and short transport; independent complete finite-group rank. No point search.')
    out=d/'result.json';assert not out.exists();out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print('PASS four focused seeds17/31/19/19',flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);a=ap.parse_args();main(a.directory)

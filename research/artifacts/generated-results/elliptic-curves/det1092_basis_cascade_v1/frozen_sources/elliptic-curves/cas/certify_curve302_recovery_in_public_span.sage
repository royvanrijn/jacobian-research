#!/usr/bin/env sage-python
"""Post-search numerical proposals accepted only by exact group identities."""
import json,argparse,sys,math
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from sage.all import QQ,EllipticCurve,RealField,matrix,vector
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))

def main(d):
    data=json.loads((d/'input.json').read_bytes());E=EllipticCurve(QQ,data['curve'])
    public=[E(P) for P in data['public_points']];found=[E(P) for P in data['recovered_points']]
    assert len(public)==31 and len(found)==24
    geometry=SourceFileLoader('recovery_height',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
    model=tuple(F(str(v)) for v in E.a_invariants())
    pts=[tuple(F(str(v)) for v in P.xy()) for P in public+found]
    gram,asym=geometry.canonical_height_gram(model,pts);RF=RealField(384)
    G=matrix(RF,[[str(gram[i][j]) for j in range(31)] for i in range(31)])
    words=[]
    for j,P in enumerate(found):
        approx=G.solve_right(vector(RF,[str(gram[i][31+j]) for i in range(31)]))
        q=[F(str(v)).limit_denominator(64) for v in approx];den=math.lcm(*(v.denominator for v in q))
        word=[int(v*den) for v in q]
        assert den<=64 and max(map(abs,word))<=64
        assert den*P==sum((word[i]*public[i] for i in range(31)),E(0))
        words.append(dict(denominator=den,word=word))
    W=matrix(QQ,[[QQ(v)/r['denominator'] for v in r['word']] for r in words])
    assert W.rank()==24 and W[:17,:].rank()==17
    result=dict(status='PASS',relations=words,public_span_dimension=31,recovered_dimension=24,
        generic_dimension=17,recovered_exceptional_directions=7,
        boundary='Public points used only after terminal target-free point waves. All24 recovered basis points lie in the rational span of the independent public31; their quotient over the original generic17 has dimension7. This is known-curve calibration, not a new rank24 or rank32 discovery.')
    out=d/'result.json';assert not out.exists();out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('EXACT24 public-span identities;7 of14 exceptional directions',flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);a=ap.parse_args();main(a.directory)

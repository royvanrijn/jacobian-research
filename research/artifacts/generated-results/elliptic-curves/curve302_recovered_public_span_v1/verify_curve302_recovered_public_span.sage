#!/usr/bin/env sage-python
"""Standalone exact rational group identities; no numerical height inputs."""
import json,argparse
from pathlib import Path
from sage.all import QQ,EllipticCurve,matrix
def main(d):
    data=json.loads((d/'input.json').read_bytes());proof=json.loads((d/'result.json').read_bytes())
    E=EllipticCurve(QQ,data['curve']);public=[E(P) for P in data['public_points']];found=[E(P) for P in data['recovered_points']]
    assert len(public)==31 and len(found)==24 and len(proof['relations'])==24
    rows=[]
    for P,r in zip(found,proof['relations']):
        den=r['denominator'];word=r['word'];assert 1<=den<=64 and len(word)==31 and max(map(abs,word))<=64
        assert den*P==sum((word[i]*public[i] for i in range(31)),E(0))
        rows.append([QQ(v)/den for v in word])
    W=matrix(QQ,rows);assert W.rank()==24 and W[:17,:].rank()==17
    assert proof['recovered_exceptional_directions']==7 and proof['status']=='PASS'
    print('PASS independent exact public-span replay:7/14 recovered',flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);a=ap.parse_args();main(a.directory)

"""Two explicitly selected nodal witnesses; run each under a 60-second cap."""
import os
os.environ['SAGE_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import hashlib
import json
import time
from pathlib import Path
from sage.all import PolynomialRing,QQ
from sage.env import SAGE_VERSION
from sage.schemes.curves.zariski_vankampen import braid_monodromy

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('family',choices=['L','K'])
    family=parser.parse_args().family;start=time.monotonic()
    row=next(r for r in json.loads((HERE/'witnesses.json').read_text()) if r['family']==family)
    R=PolynomialRing(QQ,names=('P','Q','t'));P,Q,t=R.gens()
    p=R(row['p'].replace('**','^'));q=R(row['q'].replace('**','^'))
    implicit=(P-p).resultant(Q-q,t)
    A=PolynomialRing(QQ,names=('P','Q'));equation=A(implicit)
    assert equation.total_degree()==10 and equation.degree(A.gen(1))==6
    print('IMPLICIT_READY',family,'degree10, six strands',flush=True)
    braids,components,vertical,degree=braid_monodromy(equation,arrangement=(equation,))
    assert degree==6 and vertical=={} and set(components.values())=={0}
    payload=dict(family=family,parameter=row['parameter'],degree=int(degree),
                 words=[list(map(int,b.Tietze())) for b in braids],
                 components={int(k):int(v) for k,v in components.items()},vertical={},
                 sage_version=SAGE_VERSION,elapsed_seconds=round(time.monotonic()-start,3),
                 witness_input_sha256=hashlib.sha256((HERE/'witnesses.json').read_bytes()).hexdigest(),
                 assurance='Sage/SIROCCO certified continuation for one explicitly selected nodal witness. No claim about a union with unspecified additional images.')
    (HERE/f'braids-nodal-{family}.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('CERTIFIED_BRAIDS',family,len(braids),'seconds',payload['elapsed_seconds'],flush=True)

if __name__=='__main__':main()

"""One explicitly selected six-strand SIROCCO braid calculation.

Run under Sage's Python, with an external 60-second cap. Save braids before
any monodromy enumeration. No claim about the whole parameter locus follows
from one witness.
"""
import argparse
import json
from pathlib import Path
from sage.all import PolynomialRing,QQ
from sage.schemes.curves.zariski_vankampen import braid_monodromy

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('case',choices=['generic','cusp0','cusp1'])
    args=parser.parse_args()
    ring=PolynomialRing(QQ,names=('P','Q','t'));P,Q,t=ring.gens()
    if args.case=='generic':
        p=t**6+t**3
        q=t**10+QQ(5)/3*t**7+QQ(7)/9*t**4+QQ(7)/81*t
    else:
        data=json.loads((ROOT/'critical_classes_discovery.json').read_text())[int(args.case[-1])]
        p=ring(data['p'].replace('z','t'));q=ring(data['q'].replace('z','t'))
    implicit=(P-p).resultant(Q-q,t)
    target=PolynomialRing(QQ,names=('P','Q'));equation=target(implicit)
    assert equation.degree(target.gen(1))==6
    print('implicit ready',args.case,'total degree',equation.total_degree(),flush=True)
    braids,components,vertical,degree=braid_monodromy(equation,arrangement=(equation,))
    assert degree==6 and vertical=={} and set(components.values())=={0}
    words=[list(map(int,braid.Tietze())) for braid in braids]
    payload={'case':args.case,'words':words,'degree':int(degree),'components':dict(components),'vertical':{},
             'assurance':'Sage/SIROCCO certified continuation; one fixed rational witness only'}
    (ROOT/('braids-'+args.case+'.json')).write_text(json.dumps(payload,indent=2)+'\n')
    print('CERTIFIED BRAIDS',args.case,len(words),flush=True)


if __name__=='__main__':main()

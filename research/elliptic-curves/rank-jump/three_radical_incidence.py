#!/usr/bin/env python3
"""Point-free additive norm constructors on the frozen generic panel."""
import argparse
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as panel

PROTOCOL=Path(__file__).with_name('THREE_RADICAL_INCIDENCE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_three_radical_incidence_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-three-radical-incidence-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,panel.INPUT,Path(panel.__file__),Path(r.__file__))}


def primitive(poly):
    from sage.all import ZZ,lcm,gcd
    d=lcm([c.denominator() for c in poly.list()]);cs=[ZZ(c*d) for c in poly.list()]
    content=gcd(cs);assert content
    if cs[-1]<0:content=-content
    return poly.parent()([c/content for c in cs])


def gate(f,q):
    from sage.all import ZZ
    q=primitive(q);norm=ZZ(f.resultant(q));assert norm
    n=abs(norm);delta=abs(ZZ(f.discriminant()));steps=[]
    while True:
        common=n.gcd(delta)
        if common==1:break
        steps.append(str(common));n//=common
    square=bool(n.is_square())
    return {'primitive_coefficients':list(map(str,q.list())),
        'norm':str(norm),'discriminant_gcd_steps':steps,'outside_discriminant_norm':str(n),
        'outside_norm_is_square':square,'outside_norm_square_root':str(n.sqrt()) if square else None,
        'status':'UNKNOWN_UNRAMIFIEDNESS' if square else 'PROVED_RAMIFIED_OUTSIDE_DISCRIMINANT'}


def worker(token):
    f,pts,scale=panel.model_data(token);z=f.parent().gen();rows=[]
    for triple in combinations(range(len(pts)),3):
        xs=[pts[i][0] for i in triple];s=sum(xs);pairs=sum(a*b for a,b in combinations(xs,2))
        q=-3*z*z+2*s*z+s*s-4*pairs
        # The symmetric Heron polynomial is the norm of a three-radical sum,
        # up to its forced square from simultaneous sign reversal.
        a,b,c=[x-z for x in xs];assert q==(a+b-c)**2-4*a*b
        rows.append({'generic_indices':list(triple),**gate(f,q)})
    controls=[{'generic_index':i,**gate(f,x-z)} for i,(x,y) in enumerate(pts)]
    assert all(x['outside_norm_is_square'] for x in controls)
    derivative=gate(f,f.derivative());assert derivative['outside_norm_is_square']
    return {'token':token,'status':'PASS','bindings':bindings(),'cubic_ascending':list(map(str,f.list())),
        'discriminant':str(f.discriminant()),'generic_dimension':len(pts),'rows':rows,
        'generic_controls':controls,'derivative_false_positive_control':derivative,
        'summary':{'triples':len(rows),'ramified_exclusions':sum(not x['outside_norm_is_square'] for x in rows),
                   'unresolved_survivors':sum(x['outside_norm_is_square'] for x in rows)}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'{args.token}.json',worker(args.token))
    else:
        WORK.mkdir(parents=True,exist_ok=True);rows=[]
        for case in r.read(panel.INPUT)['cases']:
            token=case['token'];path=WORK/f'{token}.json'
            if not path.exists():
                with (WORK/f'{token}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,__file__,'worker','--token',token],stdout=log,stderr=log,timeout=30)
                        reason=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason:r.write_new(path,{'token':token,'status':'UNKNOWN','reason':reason,'bindings':bindings()})
            row=r.read(path);assert row['bindings']==bindings();rows.append(row)
            print(token,row['status'],row.get('summary'),flush=True)
        r.write_new(OUTPUT,{'schema':'rank-jump.three-radical-incidence.v1','bindings':bindings(),
            'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'PARTIAL','rows':rows,
            'boundary':'Ramification exclusions for this explicit additive constructor, not a class-group bound or a rank predictor. No exceptional-point inputs.'})

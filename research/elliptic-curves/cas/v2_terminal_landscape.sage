#!/usr/bin/env sage-python
"""Bounded oracle-only full-domain CVP and exact quartic vetting after V2."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import matrix, vector, ZZ, QQ, pari
from research_runtime.store import checkpoint
from terminal_affine_cvp import AffineCVP
from visibility_lattice_v2 import ExactParity

CAS=Path(__file__).resolve().parent
b=SourceFileLoader('terminal_base',str(CAS/'diagnose_v2_terminal.sage')).load_module()
D=b.D;read=b.read;load=b.load
group=load('half_lattice_pointed_sieve.py')
geo=load('prospective_half_lattice_v3.sage')
mapper=load('factor_free_pari_mapping.sage')
helper=load('audit_curve302_m25_remaining_strict_cvp.sage')
mapper.pari.allocatemem(256000000,silent=True)


def setup(points,target,model):
    n=len(points)
    h=matrix(ZZ,geo.rounded_gram(geo.canonical_height_gram(model,(*points,target))[0],1000000))
    g=h[:n,:n];u=matrix(ZZ,pari(g).qflllgram()).transpose();inv=u.inverse()
    t=vector(QQ,g.solve_right(2*vector(ZZ,list(h.column(n))[:n])))*inv
    reduced=u*g*u.transpose()
    constant=4*h[n,n]-(t*reduced*t)
    return g,u,inv,t,reduced,constant


def evaluate(model,points,target,cword,qword,mapping=None):
    translation=[(c-q)//2 for c,q in zip(cword,qword)]
    assert all((c-q)%2==0 for c,q in zip(cword,qword))
    moved=group.short_add(model,target,group.linear_combination(model,points,translation))
    if moved is None:raise ArithmeticError('missing quotient became contained')
    class Fixed:
        @staticmethod
        def mapping(*args):return mapping
    exact=helper.coordinate(Fixed if mapping else mapper,model,points,cword,moved)
    return {'centre_word':cword,'closest_double_target_word':qword,'translation_word':translation,
            'target':list(map(str,target)),'translated_target':list(map(str,moved)),
            'height':exact['reduced_coordinate_height'],'coordinate':exact['reduced_coordinate'],
            'square_root':exact['reduced_quartic_square_root_absolute'],'mapping':exact['mapping'],
            'inside_125000':int(exact['reduced_coordinate_height'])<=125000 and int(exact['reduced_coordinate'][1])!=0}


def vet(points,target,model,out,k):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'result.json').exists(): return read(out/'result.json')
    g,u,inv,t,reduced,constant=setup(points,target,model)
    if (out/'nearest.json').exists(): nearest=read(out/'nearest.json')
    else:
        nearest=AffineCVP(reduced.rows()).nearest(t,count=k,node_limit=30000000)
        checkpoint(out/'nearest.json',nearest)
    centre_oracle=ExactParity(reduced.rows())
    distinct={}
    for row in nearest['rows']:
        qword=list(map(int,vector(ZZ,row['word'])*u));p=tuple(x%2 for x in qword)
        if p not in distinct:distinct[p]=(qword,row)
    trials=[]
    for index,(p,(qword,row)) in enumerate(distinct.items()):
        path=out/f'trial-{index:04d}.json'
        if path.exists(): trials.append(read(path));continue
        if not any(p):
            record={'parity':list(p),'status':'ZERO_PARITY_MINIMUM_CENTRE_IS_INFINITY','distance':row['distance']}
        else:
            rp=list(map(lambda x:int(x)%2,vector(ZZ,p)*inv))
            seed,_=centre_oracle.babai([rp]);proof=centre_oracle.solve(rp,seed[0],2000000)
            choices={}
            for minimum in proof['minima']:
                word=list(map(int,vector(ZZ,minimum)*u))
                point=group.linear_combination(model,points,word)
                if point[1]<0:word=[-x for x in word];point=(point[0],-point[1])
                choices[point]=word
            candidates=[]
            for point,word in sorted(choices.items()):
                candidates.append(evaluate(model,points,target,word,qword))
            record={'parity':list(p),'status':'EXACT_QUARTIC_WITNESSES',
                'half_distance_numerator':str(QQ(row['distance'])+constant),'centre_minimum_norm':proof['norm'],
                'centre_CVP':proof,'policy_centre':0,'witnesses':candidates,
                'best':min(candidates,key=lambda x:int(x['height']))}
        checkpoint(path,record);trials.append(record)
        print('VETTED',out.name,index+1,'/',len(distinct),'height',record.get('best',{}).get('height'),flush=True)
    winners=[r['best'] for r in trials if 'best' in r]
    result={'rank':len(points),'parity_universe':1<<len(points),'nearest':nearest,
        'orthogonal_constant':str(constant),'global_minimum_half_distance_numerator':str(QQ(nearest['rows'][0]['distance'])+constant),
        'distinct_vetted_parities':len(distinct),'best':min(winners,key=lambda x:int(x['height'])) if winners else None,
        'cheap_parities':sum(r.get('best',{}).get('inside_125000',False) for r in trials),
        'boundary':'Full-domain exact k-nearest half-lattice words, with all boundary ties. Quartic minimum only over these vetted representatives; not a global chart-height minimum over every parity and translate.'}
    checkpoint(out/'result.json',result);return result


def terminal(k):
    data=read(D/'coordinates.json');fixed=read(D/'fixed-M30.json')
    return vet(tuple(tuple(map(F,p)) for p in fixed['points']),tuple(map(F,data['missing_point'])),tuple(map(F,fixed['curve'])),D/'terminal-landscape',k)


def policy_charts():
    data=read(D/'coordinates.json');fixed=read(D/'fixed-M30.json')
    points=tuple(tuple(map(F,p)) for p in fixed['points']);target=tuple(map(F,data['missing_point']));model=tuple(map(F,fixed['curve']))
    g,u,inv,t,reduced,constant=setup(points,target,model)
    oracle=AffineCVP(reduced.rows());out=D/'terminal-searched-charts';out.mkdir(exist_ok=True)
    rows=[]
    for j,path in enumerate(sorted((b.V2/'calibration302/epoch-13').glob('chart-*.json'))):
        output=out/path.name
        if output.exists(): rows.append(read(output));continue
        chart=read(path);word=chart['centre']['representative'];p=vector(ZZ,[int(x)%2 for x in vector(ZZ,word)*inv])
        nearest=oracle.nearest((t-p)/2,count=1,node_limit=30000000)
        witnesses=[]
        for row in nearest['rows']:
            qword=list(map(int,(p+2*vector(ZZ,row['word']))*u))
            witnesses.append(evaluate(model,points,target,word,qword,chart['mapping']))
        result={'chart':b.rel(path),'chart_sha256':b.sha(path),'target_CVP':nearest,
            'best':min(witnesses,key=lambda x:int(x['height'])),'witnesses':witnesses}
        checkpoint(output,result);rows.append(result)
        print('SEARCHED CHART AUDIT',j+1,'height',result['best']['height'],flush=True)
    checkpoint(out/'result.json',{'count':len(rows),'cheap':sum(r['best']['inside_125000'] for r in rows),
        'best':min(rows,key=lambda r:int(r['best']['height']))})


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['terminal','searched']);ap.add_argument('--k',type=int,default=256);a=ap.parse_args()
    terminal(a.k) if a.action=='terminal' else policy_charts()

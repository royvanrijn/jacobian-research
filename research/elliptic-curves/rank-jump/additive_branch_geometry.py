#!/usr/bin/env python3
"""Branch divisors of generic-only additive norms in five fixed families."""
import argparse
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_symbolic_discriminant as atlases
import branch_divisibility_capacity as published

PROTOCOL=Path(__file__).with_name('ADDITIVE_BRANCH_GEOMETRY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_additive_branch_geometry_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_additive_branch_geometry_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-additive-branch-geometry-v1'


def export():
    specs=r.read(PROTOCOL)['families'];families={}
    for path in atlases.ATLASES:
        for f in r.read(path)['families']:
            name=f.get('family',f.get('fibration_id'))
            if name in specs:
                sections=sorted(f['sections'],key=lambda x:x['basis_index'])
                families[name]={'family':name,'A':f['A_coefficients_low_to_high'],'B':f['B_coefficients_low_to_high'],
                    'sections':[{'x':x['X'],'y':x['Y']} for x in sections]}
    pub=r.read(published.INPUT)
    families['published-r17']={'family':'published-r17',**{k:pub[k] for k in ('A','B','sections')}}
    r.write_new(INPUT,{'schema':'rank-jump.additive-branch-geometry-inputs.v1','families':[families[s] for s in specs],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
                    [Path(__file__),PROTOCOL,*atlases.ATLASES,published.INPUT]},
        'boundary':'Generic polynomial data only. No specialized points, parameters or outcome labels.'})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,Path(r.__file__))}


def primitive(q):
    from sage.all import ZZ,gcd
    cs=[ZZ(c*q.denominator()) for c in q.list()];g=gcd(cs)
    if cs[-1]<0:g=-g
    return q.parent()([c/g for c in cs])


def worker(family):
    from sage.all import QQ,GF,PolynomialRing,matrix
    source=next(x for x in r.read(INPUT)['families'] if x['family']==family)
    R=PolynomialRing(QQ,'t');A=R(source['A']);B=R(source['B']);D=primitive(-4*A**3-27*B**2)
    sections=[(R(x['x']),R(x['y'])) for x in source['sections']]
    assert A.degree()<=8 and B.degree()<=12 and D.degree()==24
    assert all(y*y==x**3+A*x+B and x.degree()<=4 and y.degree()<=6 for x,y in sections)
    norms=[];records=[]
    for triple in combinations(range(len(sections)),3):
        x,y,z=[sections[i][0] for i in triple]
        c=x*x+y*y+z*z-2*(x*y+x*z+y*z);b=2*(x+y+z);a=R(-3)
        M=matrix(R,[[c,-B*a,-B*b],[b,c-A*a,-A*b-B*a],[a,b,c-A*a]])
        q=primitive(M.det());assert q.degree()<=24
        norms.append(q);records.append({'generic_indices':list(triple),'degree':int(q.degree()),
            'norm_ascending':list(map(str,q.list())),'squarefree_witness_prime':None,'discriminant_coprime_witness_prime':None})
    primes=r.read(PROTOCOL)['limits']['primes'];reductions={}
    for p in primes:
        F=GF(p);P=PolynomialRing(F,'t');dp=P(D);values=[P(n) for n in norms];reductions[p]=values
        for row,v in zip(records,values):
            if v.degree()!=row['degree']:continue
            if row['squarefree_witness_prime'] is None and v.gcd(v.derivative())==1:row['squarefree_witness_prime']=p
            if row['discriminant_coprime_witness_prime'] is None and dp.degree()==24 and v.gcd(dp)==1:row['discriminant_coprime_witness_prime']=p
    exceptions=[];unresolved=[];pair_count=0
    for i in range(len(norms)):
        for j in range(i):
            pair_count+=1;witness=None
            for p in primes:
                a,b=reductions[p][i],reductions[p][j]
                if a.degree()==records[i]['degree'] and b.degree()==records[j]['degree'] and a.gcd(b)==1:witness=p;break
            if witness is None:unresolved.append([j,i])
            elif witness!=primes[0]:exceptions.append({'pair':[j,i],'prime':witness})
    all_pass=not unresolved and all(x['degree']==24 and x['squarefree_witness_prime'] and x['discriminant_coprime_witness_prime'] for x in records)
    return {'family':family,'status':'PASS' if all_pass else 'PARTIAL','bindings':bindings(),
        'generic_dimension':len(sections),'norms':records,'discriminant_ascending':list(map(str,D.list())),
        'pair_count':pair_count,'default_pair_witness_prime':primes[0],
        'pair_witness_exceptions':exceptions,'unresolved_pairs':unresolved,
        'summary':{'norm_count':len(norms),'degree24_count':sum(x['degree']==24 for x in records),
            'squarefree_count':sum(x['squarefree_witness_prime'] is not None for x in records),
            'discriminant_coprime_count':sum(x['discriminant_coprime_witness_prime'] is not None for x in records),
            'proved_coprime_pairs':pair_count-len(unresolved),
            'single_norm_cover_genus':11 if all_pass else 'UNKNOWN',
            'simultaneous_norm_cover_genera_first_four':[11,45,137,369] if all_pass else 'UNKNOWN'},
        'boundary':'Geometric branch and necessary norm-cover data, not actual Selmer incidence or rational elliptic solubility.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for family in r.read(PROTOCOL)['families']:
        path=WORK/f'{family}.json'
        if not path.exists():
            with (WORK/f'{family}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,__file__,'worker','--family',family],stdout=log,stderr=log,timeout=60)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='60-second timeout'
            if reason:r.write_new(path,{'family':family,'status':'UNKNOWN','reason':reason,'bindings':bindings()})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(family,row['status'],row.get('summary'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.additive-branch-geometry.v1','status':'PASS' if all(x['status']=='PASS' for x in rows) else 'PARTIAL',
        'bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'{args.family}.json',worker(args.family))
    else:globals()[args.mode]()

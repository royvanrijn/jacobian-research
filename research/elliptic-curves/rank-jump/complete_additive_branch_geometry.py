#!/usr/bin/env python3
"""Extend the frozen branch test to rational generic sections."""
import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r
import additive_branch_geometry as first

PROTOCOL=Path(__file__).with_name('ADDITIVE_BRANCH_GEOMETRY_COMPLETION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_additive_branch_geometry_v2.json'
WORK=r.ROOT/'artifacts/local/rank-jump-additive-branch-geometry-v2'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,first.INPUT,first.OUTPUT,Path(first.__file__),Path(r.__file__))}


def functions(row,R):
    F=R.fraction_field()
    def get(x):
        if isinstance(x,list):return F(R(x))
        return F(R(x['numerator_coefficients_low_to_high']))/R(x['denominator_coefficients_low_to_high'])
    return [(get(s['x']),get(s['y'])) for s in row['sections']]


def worker(family):
    from sage.all import QQ,GF,PolynomialRing,matrix,gcd,lcm
    source=next(x for x in r.read(first.INPUT)['families'] if x['family']==family)
    R=PolynomialRing(QQ,'t');A=R(source['A']);B=R(source['B']);D=first.primitive(-4*A**3-27*B**2)
    sections=functions(source,R)
    assert all(y*y==x**3+A*x+B for x,y in sections)
    norms=[];records=[]
    for triple in combinations(range(len(sections)),3):
        xs=[sections[i][0] for i in triple];e=lcm([x.denominator() for x in xs])
        x,y,z=[R(v*e) for v in xs]
        c=x*x+y*y+z*z-2*(x*y+x*z+y*z);b=2*e*(x+y+z);a=-3*e*e
        content=gcd([a,b,c]);a,b,c=[v//content for v in (a,b,c)]
        assert gcd([a,b,c])==1
        M=matrix(R,[[c,-B*a,-B*b],[b,c-A*a,-A*b-B*a],[a,b,c-A*a]])
        q=first.primitive(M.det());norms.append(q)
        records.append({'generic_indices':list(triple),'degree':int(q.degree()),
            'common_x_denominator_degree':int(e.degree()),'removed_quadratic_content_degree':int(content.degree()),
            'norm_ascending':list(map(str,q.list())),'squarefree_witness_prime':None,'discriminant_coprime_witness_prime':None})
    primes=r.read(PROTOCOL)['limits']['primes'];reductions={}
    for p in primes:
        P=PolynomialRing(GF(p),'t');dp=P(D);values=[P(n) for n in norms];reductions[p]=values
        for row,v in zip(records,values):
            if v.degree()!=row['degree']:continue
            if row['squarefree_witness_prime'] is None and v.gcd(v.derivative())==1:row['squarefree_witness_prime']=p
            if row['discriminant_coprime_witness_prime'] is None and dp.degree()==D.degree() and v.gcd(dp)==1:row['discriminant_coprime_witness_prime']=p
    exceptions=[];unresolved=[];pair_count=0
    for i in range(len(norms)):
        for j in range(i):
            pair_count+=1;witness=None
            for p in primes:
                a,b=reductions[p][i],reductions[p][j]
                if a.degree()==records[i]['degree'] and b.degree()==records[j]['degree'] and a.gcd(b)==1:witness=p;break
            if witness is None:unresolved.append([j,i])
            elif witness!=primes[0]:exceptions.append({'pair':[j,i],'prime':witness})
    all_pass=not unresolved and all(x['degree']%2==0 and x['squarefree_witness_prime'] and x['discriminant_coprime_witness_prime'] for x in records)
    degrees=sorted(x['degree'] for x in records)
    genera=[]
    if all_pass:
        for k in range(1,5):
            values=[1+(2**k)*(sum(ds)-4)//4 for ds in (degrees[:k],degrees[-k:])]
            genera.append({'classes':k,'genus_range':values})
    return {'family':family,'status':'PASS' if all_pass else 'PARTIAL','bindings':bindings(),
        'generic_dimension':len(sections),'norms':records,'discriminant_ascending':list(map(str,D.list())),
        'pair_count':pair_count,'default_pair_witness_prime':primes[0],
        'pair_witness_exceptions':exceptions,'unresolved_pairs':unresolved,
        'summary':{'norm_count':len(norms),'degree_distribution':{str(k):v for k,v in sorted(Counter(degrees).items())},
            'squarefree_count':sum(x['squarefree_witness_prime'] is not None for x in records),
            'discriminant_coprime_count':sum(x['discriminant_coprime_witness_prime'] is not None for x in records),
            'proved_coprime_pairs':pair_count-len(unresolved),'simultaneous_norm_cover_genus_ranges_first_four':genera},
        'boundary':'Primitive rational-section norm branch geometry. No exceptional inputs, parameter sweep, class group, or sufficient solubility claim.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for family in r.read(PROTOCOL)['limits']['families_to_complete']:
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
    published=next(x for x in r.read(first.OUTPUT)['rows'] if x['family']=='published-r17')
    assert published['status']=='PASS';rows.append(published)
    r.write_new(OUTPUT,{'schema':'rank-jump.additive-branch-geometry.v2','status':'PASS' if all(x['status']=='PASS' for x in rows) else 'PARTIAL',
        'bindings':bindings(),'rows':rows,'reuse':'The published-R17 row is preserved from v1, including its original bindings.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'{args.family}.json',worker(args.family))
    else:capture()

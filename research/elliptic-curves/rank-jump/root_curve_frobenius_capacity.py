#!/usr/bin/env python3
"""Bounded equation-only counts on normalized cubic root curves."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import subprocess
import sys
import retrospective as r
import generic_selmer_capacity as prior
import complete_generic_selmer_geometry as completion

PROTOCOL=Path(__file__).with_name('ROOT_CURVE_FROBENIUS_CAPACITY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_root_curve_frobenius_capacity_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_root_curve_frobenius_capacity_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-root-curve-frobenius-capacity-v1'


def trim(a):
    while len(a)>1 and not a[-1]:a.pop()
    return a


def rem(a,b,p):
    a=trim(list(a));b=trim(list(b));assert b!=[0]
    while len(a)>=len(b) and a!=[0]:
        k=len(a)-len(b);c=a[-1]*pow(b[-1],-1,p)%p
        for j,x in enumerate(b):a[j+k]=(a[j+k]-c*x)%p
        trim(a)
    return a


def gcd(a,b,p):
    while b!=[0]:a,b=b,rem(a,b,p)
    return [(x*pow(a[-1],-1,p))%p for x in a]


def value(a,t,p):
    z=0
    for x in reversed(a):z=(z*t+x)%p
    return z


def deriv(a,p):return trim([i*a[i]%p for i in range(1,len(a))] or [0])
def mod(a,p):
    return [Q(x).numerator*pow(Q(x).denominator,-1,p)%p for x in a]
def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def export():
    inputs=r.read(prior.INPUT);old=r.read(prior.OUTPUT)
    replacement=r.read(completion.INPUT)['families'][0];geometry=r.read(completion.OUTPUT)['rows'][0]
    eqs={x['family']:x for x in inputs['families']};eqs[replacement['family']]=replacement
    rows={x['family']:x for x in old['rows']};rows[geometry['family']]=geometry
    used=sorted({x['family'] for x in inputs['cases']});families=[]
    for name in used:
        e=eqs[name];g=rows[name];assert g['status']=='PASS'
        families.append({'family':name,'A':e['A'],'B':e['B'],
            'discriminant':g['discriminant_ascending'],'squarefree_factors':g['squarefree_factors'],
            'root_curve_genus':g['root_cover_genus']})
    r.write_new(INPUT,{'schema':'rank-jump.root-curve-frobenius-capacity-inputs.v1','families':families,
        'bindings':prior.bind([Path(__file__),PROTOCOL,prior.INPUT,prior.OUTPUT,completion.INPUT,completion.OUTPUT,prior.VERIFICATION])})


def eligible(f,p):
    if p<=3:return None
    try:
        A=mod(f['A'],p);B=mod(f['B'],p);D=mod(f['discriminant'],p)
        factors=[(mod(q['coefficients_ascending'],p),q['multiplicity']) for q in f['squarefree_factors']]
    except ValueError:return None
    if not D[-1] or gcd(trim(A[:]),trim(D[:]),p)!=[1]:return None
    seen=[];nodes=[]
    for q,e in factors:
        if not q[-1] or gcd(q,deriv(q,p),p)!=[1] or any(gcd(q,x,p)!=[1] for x in seen):return None
        seen.append(q)
        if e==1:continue
        assert e==2 and len(q)==2
        t=-q[0]*pow(q[1],-1,p)%p;a=value(A,t,p);b=value(B,t,p)
        x=-3*b*pow(2*a,-1,p)%p
        at=value(deriv(A,p),t,p);bt=value(deriv(B,p),t,p)
        att=value(deriv(deriv(A,p),p),t,p);btt=value(deriv(deriv(B,p),p),t,p)
        assert (x**3+a*x+b)%p==0 and (3*x*x+a)%p==0 and (at*x+bt)%p==0
        delta=(at*at-6*x*(att*x+btt))%p
        if not x or not delta:return None
        leg=pow(delta,(p-1)//2,p);correction=1 if leg==1 else -1;assert leg in (1,p-1)
        nodes.append({'base':t,'double_root':x,'tangent_discriminant':delta,'normalization_correction':correction})
    return A,B,D,nodes


def worker(family):
    f=next(x for x in r.read(INPUT)['families'] if x['family']==family);pr=r.read(PROTOCOL);counts=[]
    for p in range(5,pr['limits']['prime_bound']+1,2):
        if not prime(p):continue
        data=eligible(f,p)
        if data is None:continue
        A,B,D,nodes=data;hist=[0,0,0,0];fibres=[]
        cubes=[x*x*x%p for x in range(p)]
        for t in range(p):
            a,b=value(A,t,p),value(B,t,p)
            roots=[x for x in range(p) if (cubes[x]+a*x+b)%p==0]
            hist[len(roots)]+=1;fibres.append(len(roots))
        inf=[x for x in range(p) if (cubes[x]+A[8]*x+B[12])%p==0]
        raw=sum(fibres)+len(inf);total=raw+sum(n['normalization_correction'] for n in nodes)
        g=f['root_curve_genus'];assert abs(total-p-1)<=isqrt(4*g*g*p)
        counts.append({'prime':p,'A_mod_p':A,'B_mod_p':B,'cubic_root_histogram':hist,
            'finite_fibre_counts':fibres,'infinity_roots':inf,'nodes':nodes,
            'raw_root_curve_count':raw,'normalized_root_curve_count':total,
            'odd_count':bool(total%2),'frobenius_trace':p+1-total})
        if len(counts)==pr['limits']['good_primes_per_family']:break
    status='PASS' if len(counts)==pr['limits']['good_primes_per_family'] else 'UNKNOWN'
    result={'family':family,'root_curve_genus':f['root_curve_genus'],'status':status,'counts':counts,
        'odd_count_primes':[x['prime'] for x in counts if x['odd_count']],
        'bindings':prior.bind([Path(__file__),PROTOCOL,INPUT,Path(r.__file__)])}
    r.write_new(WORK/(family+'.json'),result)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for f in r.read(INPUT)['families']:
        family=f['family'];path=WORK/(family+'.json')
        if not path.exists():
            error=None
            with (WORK/(family+'.log')).open('x') as log:
                try:
                    p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--family',family],stdout=log,stderr=log,timeout=30)
                    if p.returncode:error='worker failure'
                except subprocess.TimeoutExpired:error='bounded timeout'
            if error:r.write_new(path,{'family':family,'status':'UNKNOWN','reason':error})
        row=r.read(path);rows.append(row);print(family,row['status'],[(x['prime'],x['normalized_root_curve_count']) for x in row.get('counts',[])],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.root-curve-frobenius-capacity.v1','rows':rows,
        'bindings':prior.bind([Path(__file__),PROTOCOL,INPUT]),
        'boundary':'Normalized finite-field root-curve counts only. No elliptic rational points, rank labels or exceptional classes entered workers.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':worker(args.family)
    else:globals()[args.mode]()

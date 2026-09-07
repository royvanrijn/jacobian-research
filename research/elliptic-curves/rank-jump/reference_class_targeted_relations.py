#!/usr/bin/env python3
"""One equation-only targeted norm wave for positive strict-class construction."""
import argparse
import hashlib
from pathlib import Path
from math import gcd,prod
import runpy
import subprocess
import sys
import time
import retrospective as r
import seeded_reference_class as seed

PROTOCOL=Path(__file__).with_name('REFERENCE_CLASS_TARGETED_RELATIONS_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_reference_class_targeted_relations_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_reference_class_targeted_relations_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-class-targeted-relations-v1'
LEGACY=r.ROOT/'elliptic-curves/cas/pursue_small_conductor_class_target_strips.sage'

def prepare():
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    old=runpy.run_path(str(LEGACY));p=r.read(PROTOCOL);ref=r.read(seed.REFERENCE)
    matrix=old['prior_state'](6)
    selected=old['selection'](matrix,6,p['bounds']['vmax'],p['bounds']['targets_max'])
    data,nf,w,lookup=old['context']()
    assert str(nf.disc())==ref['field_discriminant']
    assert [str(nf.nf_get_pol().polcoef(i)) for i in range(4)]==list(map(str,ref['cubic_ascending']))
    out={'schema':'rank-jump.reference-class-targeted-relations-inputs.v1',
         'cubic_ascending':ref['cubic_ascending'],'field_discriminant':ref['field_discriminant'],
         'S_finite':ref['S_finite'],'selection':selected,
         'fixed_a':data['integral_norm_generator']['fixed_a'],
         'w_power_basis':data['integral_norm_generator']['w_power_basis'],
         'sl2_matrix':data['sl2_matrix'],'binary_cubic_descending':data['reduced_binary_cubic_descending'],
         'slope_scale':str(old['SLOPE_SCALE']),
         'bindings':seed.bindings([Path(__file__),PROTOCOL,LEGACY,seed.REFERENCE,seed.POOL])}
    r.write_new(INPUT,out);print('Prepared',len(selected['targets']),'equation-only targets',flush=True)

def residues(values,modulus):
    levels=[values]
    while len(levels[-1])>1:
        old=levels[-1];levels.append([prod(old[i:i+2]) for i in range(0,len(old),2)])
    now=[modulus%levels[-1][0]]
    for level in reversed(levels[:-1]):now=[now[i//2]%v for i,v in enumerate(level)]
    return now

def strip(n,support):
    while True:
        d=gcd(n,support)
        if d==1:return n
        n//=d

def worker():
    from sage.all import QQ,PolynomialRing,pari,prime_range
    p=r.read(PROTOCOL);inp=r.read(INPUT)
    pari.allocatemem(64000000,p['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(inp['cubic_ascending']);nf=pari.nfinit([pari(f),inp['S_finite']])
    assert str(nf.disc())==inp['field_discriminant']
    w=pari.Mod(pari(R(inp['w_power_basis'])),pari(f));a=int(inp['fixed_a']);M=inp['sl2_matrix']
    c0,c1,c2,c3=map(int,inp['binary_cubic_descending'])
    primorial=prod(map(int,prime_range(p['bounds']['smooth_bound']+1)));scale=int(inp['slope_scale'])
    completed=[];total=0;hits=0
    for k,target in enumerate(inp['selection']['targets']):
        start=time.monotonic();v1,v2=target['lattice_basis'];pairs=[];values=[];seen=set()
        for v in range(target['old_vmax']+1,p['bounds']['vmax']+1):
            for slope in target['root_slopes_scaled']:
                center=slope*v//scale
                for u in range(center-1,center+2):
                    if (u,v) in seen or gcd(u,v)!=1:continue
                    seen.add((u,v));m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                    if gcd(m,n)!=1:continue
                    value=((c0*m+c1*n)*m+c2*n*n)*m+c3*n*n*n
                    assert value and value%target['p']==0 and (m-target['norm_slope']*n)%target['p']==0
                    pairs.append((u,v,m,n));values.append(value)
        digest=hashlib.sha256();relations=[]
        if values:
            for pair,value,remainder in zip(pairs,values,residues(list(map(abs,values)),primorial)):
                rem=strip(abs(value),gcd(abs(value),remainder))
                digest.update(('%s,%s,%s,%s\n'%(pair[0],pair[1],value,rem)).encode())
                if rem!=1:continue
                u,v,m,n=pair;alpha=a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
                assert pari.nfeltnorm(nf,alpha)==a*a*value
                fac=pari.idealfactor(nf,alpha)
                assert pari.idealhnf(nf,pari.idealfactorback(nf,fac))==pari.idealhnf(nf,alpha)
                factors=[]
                for j in range(fac.nrows()):
                    P,e=fac[j,0],int(fac[j,1]);q=int(P[0]);assert e>0 and q<=p['bounds']['smooth_bound']
                    factors.append({'p':q,'hnf':str(pari.idealhnf(nf,P)),'e':int(P[2]),'f':int(P[3]),'valuation':e})
                assert any(c['p']==target['p'] and c['hnf']==target['hnf'] for c in factors)
                relations.append({'u':u,'v':v,'m':m,'n':n,'norm':str(a*a*value),
                    'alpha_ascending':[str(pari.lift(alpha).polcoef(i)) for i in range(3)],'ideal_factorization':factors})
        chunk={'target':target,'candidate_count':len(values),'population_digest_sha256':digest.hexdigest(),
               'relations':relations,'elapsed_seconds':time.monotonic()-start}
        path=WORK/('target_%04d.json'%k);r.write_new(path,chunk)
        completed.append({'path':str(path.relative_to(r.ROOT)),'sha256':r.digest(path.read_bytes())})
        total+=len(values);hits+=len(relations)
        print('TARGET',k+1,'CANDIDATES',total,'PRINCIPAL_RELATIONS',hits,flush=True)
    r.write_new(WORK/'worker.json',{'status':'COMPLETED','chunks':completed,'candidate_count':total,'relation_occurrences':hits})

def capture():
    WORK.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    with (WORK/'worker.log').open('x') as log:
        try:
            proc=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,
                timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            reason=None if proc.returncode==0 else 'worker failure'
        except subprocess.TimeoutExpired:reason='bounded timeout'
    chunks=[];count=0;hits=0
    for path in sorted(WORK.glob('target_*.json')):
        row=r.read(path);chunks.append({'path':str(path.relative_to(r.ROOT)),'sha256':r.digest(path.read_bytes())})
        count+=row['candidate_count'];hits+=len(row['relations'])
    out={'schema':'rank-jump.reference-class-targeted-relations.v1',
        'status':'COMPLETED' if reason is None else 'PARTIAL','reason':reason,
        'elapsed_seconds':time.monotonic()-start,'completed_targets':len(chunks),
        'candidate_count':count,'principal_relation_occurrences':hits,'chunks':chunks,
        'bindings':seed.bindings([Path(__file__),PROTOCOL,INPUT]),
        'positive_endpoint':'PENDING_STRICT_CLASS_EXTRACTION','boundary':'No new strict class or rank asserted from principal-relation count.'}
    r.write_new(OUTPUT,out);print(len(chunks),'targets',hits,'principal relations',reason,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['prepare','capture','worker']);args=a.parse_args();globals()[args.mode]()

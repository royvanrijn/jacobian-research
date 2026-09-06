#!/usr/bin/env python3
"""Bounded point-independent norm-relation pilot on the frozen103b2 pair."""
import argparse
from fractions import Fraction as Q
from math import gcd,prod
from pathlib import Path
import subprocess
import sys
import retrospective as r
import matched103b2_class_boundary as prior
import fresh_governing_panel as base
from verify_unpointed_governing_norm import Algebra

CAS=r.ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
import prepare_small_conductor_norm_form as forms

PROTOCOL=Path(__file__).with_name('MATCHED103B2_NORM_RELATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_matched103b2_norm_relations_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-matched103b2-norm-relations-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,prior.OUTPUT,base.INPUT,base.OUTPUT,Path(prior.__file__),Path(base.__file__),Path(forms.__file__),Path(r.__file__))}


def form_worker(token):
    row,f,pts,primes,nf=prior.setup(token)
    K=Algebra(list(map(str,f.list())));one=K.elt([1])
    bs=[K.elt([str(v.polcoef(i)) for i in range(3)]) for v in nf.nf_get_zk()]
    coords=lambda v,basis:[sum(x*y for x,y in zip(line,v)) for line in forms.inverse(list(map(list,zip(*basis))))]
    pair=coords(K.mul(bs[1],bs[2]),bs)
    w=K.add(bs[1],K.elt([-pair[2]]));t=K.add(bs[2],K.elt([-pair[1]]));normal=[one,w,t]
    table=[[coords(K.mul(u,v),normal) for v in normal] for u in normal]
    assert all(v.denominator==1 for line in table for a in line for v in a)
    a,b,c,d=-table[1][1][2],table[1][1][1],-table[2][2][2],table[2][2][1]
    assert table[1][2]==[-a*d,0,0] and table[1][1]==[-a*c,b,-a] and table[2][2]==[-b*d,d,-c]
    initial=list(map(int,[a,b,c,d]));reduced,M,history=forms.reduce_form(initial)
    assert forms.discriminant(reduced)==int(nf.disc())
    assert M[0]*M[3]-M[1]*M[2]==1
    for m,n in [(1,0),(0,1),(1,1),(-1,1)]:
        u,v=M[0]*m+M[1]*n,M[2]*m+M[3]*n
        beta=K.add(K.elt([a*u]),tuple(v*x for x in w));value=sum(x*m**(3-i)*n**i for i,x in enumerate(reduced))
        assert K.norm(beta)==a*a*value
    mono=next(x['reduction']['reduced_cubic_ascending'] for x in r.read(prior.OUTPUT)['rows'] if x['token']==token)
    return {'status':'PASS','cubic_ascending':list(map(str,f.list())),'maximal_order_basis':[list(map(str,v)) for v in bs],
       'field_discriminant':str(nf.disc()),'binary_coefficients_descending':list(map(str,reduced)),
       'monic_coefficients_descending':list(reversed(mono)),'initial_binary_coefficients':initial,'SL2_matrix':M,'reduction_history':history,
       'fixed_a':str(a),'w_ascending':list(map(str,w)),
       'coefficient_bits':{'binary':max(abs(x).bit_length() for x in reduced),'monic':max(abs(int(x)).bit_length() for x in mono)}}


def primes_to(n):
    sieve=bytearray(b'\1')*(n+1);sieve[:2]=b'\0\0'
    for p in range(2,int(n**0.5)+1):
        if sieve[p]:sieve[p*p:n+1:p]=b'\0'*(((n-p*p)//p)+1)
    return [p for p in range(2,n+1) if sieve[p]]


def norm_worker(token):
    import hashlib
    source=r.read(WORK/f'{token}-form.json');assert source['status']=='PASS'
    p=r.read(PROTOCOL);primes=primes_to(p['smooth_bound']);primorial=prod(primes);rows=[]
    for name,key in [('reduced_monic','monic_coefficients_descending'),('maximal_order_binary','binary_coefficients_descending')]:
        coefficients=list(map(int,source[key]));digest=hashlib.sha256();count=0;hits=[];small=[]
        for m in range(-p['box'],p['box']+1):
            for n in range(1,p['box']+1):
                if gcd(m,n)!=1:continue
                value=sum(v*m**(3-i)*n**i for i,v in enumerate(coefficients));assert value
                remain=abs(value)
                while True:
                    factor=gcd(remain,primorial)
                    if factor==1:break
                    remain//=factor
                digest.update(f'{m},{n},{value},{remain}\n'.encode());count+=1
                small.append((remain.bit_length(),m,n,str(remain)));small=sorted(small)[:10]
                if remain==1:
                    cofactor=abs(value);factors=[]
                    for ell in primes:
                        exponent=0
                        while cofactor%ell==0:cofactor//=ell;exponent+=1
                        if exponent:factors.append([ell,exponent])
                        if cofactor==1:break
                    assert cofactor==1;hits.append({'m':m,'n':n,'value':str(value),'factors':factors})
        rows.append({'arm':name,'count':count,'digest':digest.hexdigest(),'smooth_count':len(hits),'smooth_values':hits,'smallest_remainders':small})
        print(token,name,'smooth',len(hits),'of',count,flush=True)
    return {'status':'PASS','rows':rows}


def ideal_worker(token):
    from sage.all import pari,QQ,PolynomialRing,matrix,GF,AA,ZZ
    from research_runtime.local_kummer import LocalSquareclasses
    norm=r.read(WORK/f'{token}-norm.json')
    if norm['status']!='PASS':return {'status':'UNKNOWN','reason':'norm gate incomplete'}
    records=norm['rows'][1]['smooth_values'][:r.read(PROTOCOL)['limits']['maximum_relations']]
    if not records:return {'status':'NO_SMOOTH_RELATIONS','verified_principal_relations':0,'independent_strict_classes':0,'scope':'Empty frozen input span, not absence of field classes.'}
    old,f,pts,primes,nf=prior.setup(token);source=r.read(WORK/f'{token}-form.json');S=old['local']['S_finite']
    a=QQ(source['fixed_a']);theta=pari.Mod('z',pari(f));w=pari(f.parent()(list(map(QQ,source['w_ascending']))))(theta);M=source['SL2_matrix']
    # Full norm support: fixed a^2 cannot be omitted from ideal verification.
    fixed_factors=list(ZZ(a).abs().factor(proof=True));assert all(p.is_prime(proof=True) for p,e in fixed_factors)
    betas=[];columns=set();rows=[]
    for rec in records:
        m,n=rec['m'],rec['n'];u,v=M[0]*m+M[1]*n,M[2]*m+M[3]*n;beta=pari(a*u)+v*w
        assert pari.nfeltnorm(nf,beta)==a*a*int(rec['value'])
        support=sorted({int(p) for p,e in fixed_factors}|{p for p,e in rec['factors']});I=pari.idealhnf(nf,1);vals={};ideal_rows=[]
        for p in support:
            for j,P in enumerate(pari.idealprimedec(nf,p)):
                e=int(pari.idealval(nf,beta,P));assert e>=0
                if e:
                    I=pari.idealmul(nf,I,pari.idealpow(nf,P,e));ideal_rows.append({'prime':p,'prime_index':j,'exponent':e})
                    if p not in S and e%2:vals[p,j]=1;columns.add((p,j))
        assert pari.idealhnf(nf,I)==pari.idealhnf(nf,beta)
        betas.append(beta);rows.append({'m':m,'n':n,'ideal_factorization':ideal_rows,'outside_parity':vals})
    columns=sorted(columns);constraints=[[int(col in row['outside_parity']) for row in rows] for col in columns]
    for p in S:
        chars=LocalSquareclasses(nf,p);sigs=[chars.signature(b) for b in betas]
        constraints.extend([[s[j] for s in sigs] for j in range(len(sigs[0]))])
    roots=f.roots(AA,multiplicities=False)
    for root in roots:
        constraints.append([int(f.parent()([QQ(pari.lift(b).polcoef(i)) for i in range(3)])(root)<0) for b in betas])
    kernel=matrix(GF(2),constraints).right_kernel().basis_matrix()
    # A nonzero kernel is only a candidate squareclass source. Do not assert
    # independence until separate characters certify nontriviality.
    for row in rows:row['outside_parity']=[list(x) for x in sorted(row['outside_parity'])]
    return {'status':'PASS_RELATION_AUDIT','fixed_a_factors':[[int(p),int(e)] for p,e in fixed_factors],
       'verified_principal_relations':len(rows),'relations':rows,'outside_columns':[list(c) for c in columns],
       'strict_constraint_matrix':constraints,'strict_product_kernel':[[int(x) for x in v] for v in kernel],
       'candidate_strict_products':int(kernel.nrows()),'independent_strict_classes':'UNKNOWN' if kernel.nrows() else 0,
       'boundary':'Dependency products satisfy strict conditions; nontriviality, quotient independence and elliptic rationality are separate gates.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(PROTOCOL)['cases']:
        row={'token':token}
        for stage in ['form','norm','ideal']:
            path=WORK/f'{token}-{stage}.json'
            if not path.exists():
                with (WORK/f'{token}-{stage}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token,'--stage',stage],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits'][stage+'_seconds'])
                        error=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:error='bounded timeout'
                if error:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':error})
            value=r.read(path);assert value['bindings']==bindings();row[stage]=value;print(token,stage,value['status'],flush=True)
        rows.append(row)
    r.write_new(OUTPUT,{'schema':'rank-jump.matched103b2-norm-relations.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');p.add_argument('--stage');args=p.parse_args()
    if args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        result={'form':form_worker,'norm':norm_worker,'ideal':ideal_worker}[args.stage](args.token)
        r.write_new(WORK/f'{args.token}-{args.stage}.json',{'bindings':bindings(),**result})

#!/usr/bin/env python3
"""Bounded real-root relation pilot and independent strict-class circuit gates."""
import argparse
from pathlib import Path
from math import gcd,prod
import hashlib
import subprocess
import sys
import retrospective as r
import matched103b2_norm_relations as old
import matched103b2_class_boundary as prior

PROTOCOL=Path(__file__).with_name('MATCHED103B2_ROOT_CIRCUIT_PROTOCOL.json')
INPUT=old.OUTPUT
SEED=r.OUT/'rank_jump_matched103b2_retained_norm_audit_v1.json'
OUTPUT=r.OUT/'rank_jump_matched103b2_root_circuits_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-matched103b2-root-circuits-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT,SEED,Path(old.__file__),Path(prior.__file__),Path(r.__file__))}


def norm_worker(token):
    from sage.all import QQ,AA,PolynomialRing,RealIntervalField
    row=next(x for x in r.read(INPUT)['rows'] if x['token']==token);form=row['form'];spec=r.read(PROTOCOL)['limits']
    c=list(map(int,form['binary_coefficients_descending']));F=PolynomialRing(QQ,'z')(list(reversed(c)));roots=F.roots(AA,multiplicities=False);assert len(roots)==3
    primes=old.primes_to(spec['smooth_bound']);P=prod(primes);digest=hashlib.sha256();seen=set();accepted=[];root_records=[];count=0;best=[]
    for j,root in enumerate(roots):
        interval=RealIntervalField(256)(root);lo=interval.lower().exact_rational();hi=interval.upper().exact_rational();assert AA(lo)<root<AA(hi)
        floors=[]
        for n in range(1,spec['denominators']+1):
            m0=int((root*n).floor());assert QQ(m0)<root*n<QQ(m0+1);floors.append(m0)
            for shift in [-1,0,1]:
                m=m0+shift
                if gcd(m,n)!=1 or (m,n) in seen:continue
                seen.add((m,n));count+=1
                value=((c[0]*m+c[1]*n)*m+c[2]*n*n)*m+c[3]*n**3;assert value
                rem=abs(value)
                while gcd(rem,P)>1:rem//=gcd(rem,P)
                digest.update(f'{j},{m},{n},{value},{rem}\n'.encode())
                best=sorted(best+[(rem.bit_length(),j,m,n,str(rem))])[:10]
                if rem.bit_length()<=spec['residual_bits']:
                    accepted.append({'root_index':j,'m':m,'n':n,'value':str(value),'remainder':str(rem)})
        root_records.append({'lower':str(lo),'upper':str(hi),'floors':floors})
    return {'status':'PASS','count':count,'roots':root_records,'digest':digest.hexdigest(),'accepted_count':len(accepted),
             'accepted':accepted[:spec['accepted_relations_per_case']],'smallest_remainders':best}


def audit_worker(token):
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,matrix,pari
    from research_runtime.local_kummer import LocalSquareclasses
    spec=r.read(PROTOCOL)['limits'];norm=r.read(WORK/f'{token}-norm.json')
    if norm['status']!='PASS':return {'status':'UNKNOWN','reason':'norm stage incomplete'}
    raw=next(x for x in r.read(INPUT)['rows'] if x['token']==token)['form']
    oldrow,f,pts,primes,nf=prior.setup(token);S=oldrow['local']['S_finite'];assert oldrow['local']['strict_generic_dimension']==0
    a=QQ(raw['fixed_a']);theta=pari.Mod('z',pari(f));w=pari(f.parent()(list(map(QQ,raw['w_ascending']))))(theta);M=raw['SL2_matrix']
    records=[dict(x) for x in norm['accepted']]
    seed=next(x for x in r.read(SEED)['rows'] if x['token']==token)
    for x in seed['retained_principal_relations']:
        if not any((q['m'],q['n'])==(x['m'],x['n']) for q in records):records.append({'m':x['m'],'n':x['n'],'value':x['polynomial_value'],'remainder':x['residual'],'source':'retained_box'})
    if not records:return {'status':'NO_RELATIONS','relation_count':0,'valuation_rank_Q':0,'valuation_rank_F2':0,'strict_product_kernel_dimension':0,'certified_independent_strict_classes':0}
    smallprimes=old.primes_to(spec['smooth_bound']);fixed=list(ZZ(a).abs().factor(proof=True));betas=[];divisors=[];columns=set();proofs=[]
    for rec in records:
        m,n=rec['m'],rec['n'];beta=pari(a*(M[0]*m+M[1]*n))+(M[2]*m+M[3]*n)*w
        assert pari.nfeltnorm(nf,beta)==pari(a*a*ZZ(rec['value']))
        factors={};q=abs(int(rec['value']))//int(rec['remainder'])
        for p in smallprimes:
            e=0
            while q%p==0:q//=p;e+=1
            if e:factors[p]=e
            if q==1:break
        assert q==1
        for p,e in ZZ(rec['remainder']).factor(proof=True):factors[int(p)]=factors.get(int(p),0)+int(e)
        for p,e in fixed:factors[int(p)]=factors.get(int(p),0)+2*int(e)
        assert prod(p**e for p,e in factors.items())==abs(int(pari.nfeltnorm(nf,beta)))
        I=pari.idealhnf(nf,1);vals={}
        for p in sorted(factors):
            assert ZZ(p).is_prime(proof=True)
            for j,P in enumerate(pari.idealprimedec(nf,p)):
                v=int(pari.idealval(nf,beta,P));assert v>=0
                if v:I=pari.idealmul(nf,I,pari.idealpow(nf,P,v));vals[p,j]=v;columns.add((p,j))
        assert pari.idealhnf(nf,I)==pari.idealhnf(nf,beta)
        betas.append(beta);divisors.append(vals);proofs.append({**rec,'beta_ascending':[str(pari.lift(beta).polcoef(i)) for i in range(3)],'norm_factors':[[p,e] for p,e in sorted(factors.items())],
                                                            'ideal_valuations':[[p,j,e] for (p,j),e in sorted(vals.items())]})
    columns=sorted(columns);A=[[v.get(c,0) for v in divisors] for c in columns]
    rQ=int(matrix(QQ,A).rank());r2=int(matrix(GF(2),A).rank());constraints=[[x%2 for x in row] for row in A];local=[]
    for p in S:
        chars=LocalSquareclasses(nf,p);sigs=[list(chars.signature(b)) for b in betas]
        constraints.extend([[s[j] for s in sigs] for j in range(len(sigs[0]))]);local.append({'place':p,'signatures':sigs})
    roots=f.roots(AA,multiplicities=False)
    signs=[[int(f.parent()([QQ(pari.lift(b).polcoef(i)) for i in range(3)])(root)<0) for root in roots] for b in betas]
    constraints.extend([[s[j] for s in signs] for j in range(3)]);local.append({'place':'infinity','signatures':signs})
    kernel=matrix(GF(2),constraints).right_kernel().basis_matrix();products=[]
    for v in kernel.rows()[:spec['maximum_candidate_classes']]:
        b=pari.Mod(1,pari(f))
        for bit,element in zip(v,betas):
            if bit:b*=element
        assert ZZ(pari.nfeltnorm(nf,b)).is_square()
        products.append((list(map(int,v)),b))
    signatures=[0]*len(products);character_places=[];offset=0
    for p in r.primes(spec['character_prime_bound']):
        if p in S or int(f.discriminant())%p==0:continue
        for j,P in enumerate(pari.idealprimedec(nf,p)):
            bits=[int(pari.nfislocalpower(nf,P,b,2)==0) for v,b in products]
            for i,bit in enumerate(bits):signatures[i]|=bit<<offset
            character_places.append([p,j]);offset+=1
        if r.rank(signatures)==len(products):break
    return {'status':'PASS','relation_count':len(betas),'relations':proofs,'valuation_columns':[list(c) for c in columns],
       'integer_valuation_matrix':A,'valuation_rank_Q':rQ,'valuation_rank_F2':r2,'ideal_two_torsion_capacity':rQ-r2,
       'strict_squareclass_dimension_upper_bound_from_units_and_relations':2+rQ-r2,
       'local':local,'strict_constraint_matrix':constraints,'strict_product_kernel':[[int(x) for x in v] for v in kernel],
       'strict_product_kernel_dimension':int(kernel.nrows()),'tested_product_masks':[v for v,b in products],
       'character_places':character_places,'product_character_signatures':signatures,
       'certified_independent_strict_classes':r.rank(signatures),
       'boundary':'Independence is certified only by the retained characters. Strict class solubility on the elliptic curve is UNKNOWN. Matrix capacity bounds only products of these supplied principal generators.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(PROTOCOL)['cases']:
        row={'token':token}
        for stage in ['norm','audit']:
            path=WORK/f'{token}-{stage}.json'
            if not path.exists():
                with (WORK/f'{token}-{stage}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token,'--stage',stage],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits'][stage+'_seconds'])
                        error=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:error='bounded timeout'
                if error:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':error})
            value=r.read(path);assert value['bindings']==bindings();row[stage]=value;print(token,stage,value['status'],value.get('accepted_count',value.get('certified_independent_strict_classes')),flush=True)
        rows.append(row)
    r.write_new(OUTPUT,{'schema':'rank-jump.matched103b2-root-circuits.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');p.add_argument('--stage');args=p.parse_args()
    if args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        result={'norm':norm_worker,'audit':audit_worker}[args.stage](args.token)
        r.write_new(WORK/f'{args.token}-{args.stage}.json',{'bindings':bindings(),**result})

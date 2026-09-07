#!/usr/bin/env python3
"""Replay bounded norm forms and audit retained small-cofactor principal ideals."""
import argparse
from fractions import Fraction as Q
from math import prod,gcd
from pathlib import Path
import hashlib
import sys
import retrospective as r
import matched103b2_norm_relations as source
import matched103b2_class_boundary as prior
from verify_unpointed_governing_norm import Algebra

PROTOCOL=Path(__file__).with_name('MATCHED103B2_RETAINED_NORM_AUDIT_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_matched103b2_retained_norm_audit_v1.json'


def compute():
    from sage.all import ZZ,QQ,pari,prime_range,RealBallField
    pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    data=r.read(source.OUTPUT)
    for name,sha in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    spec=r.read(source.PROTOCOL);prime_list=list(map(int,prime_range(spec['smooth_bound']+1)));P=prod(prime_list)
    rows=[]
    for row in data['rows']:
        token=row['token'];form=row['form'];K=Algebra(form['cubic_ascending']);w=K.elt(form['w_ascending']);a=Q(form['fixed_a']);M=form['SL2_matrix']
        coefficients=list(map(int,form['binary_coefficients_descending']))
        # Verify the norm cubic at four distinct projective points, so all four coefficients agree.
        def element(m,n):
            u,v=M[0]*m+M[1]*n,M[2]*m+M[3]*n
            return K.add(K.elt([a*u]),tuple(v*x for x in w))
        for m,n in [(1,0),(0,1),(1,1),(-1,1)]:
            assert K.norm(element(m,n))==a*a*sum(c*m**(3-i)*n**i for i,c in enumerate(coefficients))
        counts=[]
        for arm in row['norm']['rows']:
            c=list(map(int,form['binary_coefficients_descending' if arm['arm']=='maximal_order_binary' else 'monic_coefficients_descending']))
            h=hashlib.sha256();count=0;hits=0
            for m in range(-spec['box'],spec['box']+1):
                for n in range(1,spec['box']+1):
                    if gcd(m,n)!=1:continue
                    value=((c[0]*m+c[1]*n)*m+c[2]*n*n)*m+c[3]*n**3
                    rem=abs(value)
                    while gcd(rem,P)>1:rem//=gcd(rem,P)
                    h.update(f'{m},{n},{value},{rem}\n'.encode());count+=1;hits+=int(rem==1)
            assert h.hexdigest()==arm['digest'] and count==arm['count'] and hits==arm['smooth_count']
            counts.append({'arm':arm['arm'],'replayed_values':count,'smooth_values':hits})
        old,f,pts,primes,nf=prior.setup(token);S=old['local']['S_finite'];theta=pari.Mod('z',pari(f))
        retained=[v for v in row['norm']['rows'][1]['smallest_remainders'] if v[0]<=r.read(PROTOCOL)['limits']['maximum_remainder_bits']]
        audits=[]
        for bits,m,n,remainder in retained:
            beta=pari(f.parent()([QQ(x) for x in element(m,n)]))(theta)
            value=sum(c*m**(3-i)*n**i for i,c in enumerate(coefficients));N=ZZ(str(K.norm(element(m,n))))
            assert str(N)==str(a*a*value)==str(pari.nfeltnorm(nf,beta))
            small=abs(value)//int(remainder);factorization={}
            for p in prime_list:
                e=0
                while small%p==0:small//=p;e+=1
                if e:factorization[p]=e
                if small==1:break
            assert small==1
            for p,e in ZZ(remainder).factor(proof=True):factorization[int(p)]=factorization.get(int(p),0)+int(e)
            assert prod(p**e for p,e in factorization.items())==abs(value)
            for p,e in ZZ(str(a)).abs().factor(proof=True):factorization[int(p)]=factorization.get(int(p),0)+2*int(e)
            assert prod(p**e for p,e in factorization.items())==abs(N)
            ideal=pari.idealhnf(nf,1);entries=[];localized=[]
            for p in sorted(factorization):
                assert ZZ(p).is_prime(proof=True)
                for j,prime in enumerate(pari.idealprimedec(nf,p)):
                    v=int(pari.idealval(nf,beta,prime));assert v>=0
                    if not v:continue
                    ideal=pari.idealmul(nf,ideal,pari.idealpow(nf,prime,v))
                    entry={'prime':p,'prime_index':j,'residue_degree':int(prime[3]),'exponent':v,'norm':str(ZZ(p)**int(prime[3]))}
                    entries.append(entry)
                    if p not in S and v%2:localized.append(entry)
            assert pari.idealhnf(nf,ideal)==pari.idealhnf(nf,beta)
            audits.append({'m':m,'n':n,'residual':remainder,'polynomial_value':str(value),'full_norm':str(N),
              'beta_ascending':list(map(str,element(m,n))),'norm_factors':[[p,e] for p,e in sorted(factorization.items())],
              'principal_ideal_factors':entries,'localized_odd_prime_ideal_support':localized,
              'principal_relation_verified':True,'is_strict_class_candidate':False if localized else 'NEEDS_LOCAL_CHECKS'})
        R=RealBallField(96);D=ZZ(nf.disc()).abs();bound=12*R(D).log()**2
        ceil_lower=ZZ(bound.lower().ceil());ceil_upper=ZZ(bound.upper().ceil());assert ceil_lower==ceil_upper
        B=int(ceil_upper)
        rows.append({'token':token,'coefficient_bits':form['coefficient_bits'],'norm_digest_replays':counts,
          'retained_principal_relations':audits,'localized_relation_rank':int(bool(audits) and bool(audits[0]['localized_odd_prime_ideal_support'])) if len(audits)<=1 else 'UNKNOWN',
          'Bach_GRH_generating_cutoff':B,'Bach_interval':str(bound),
          'audited_ideal_factors_within_GRH_cutoff':all(int(c['norm'])<=B for a in audits for c in a['principal_ideal_factors']),
          'independent_additional_strict_classes':0,
          'boundary':'No strict incidence claim from norm smoothness. Prime-ideal generation at the displayed cutoff is conditional on GRH; the finite ideal identities are unconditional.'})
        print(token,form['coefficient_bits'],'retained relations',len(audits),'Bach',B,flush=True)
    files=(Path(__file__),PROTOCOL,source.OUTPUT,Path(source.__file__),prior.OUTPUT,Path(prior.__file__),Path(r.__file__))
    return {'schema':'rank-jump.matched103b2-retained-norm-audit.v1','status':'PASS',
      'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},'rows':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();out=compute()
    if a.mode=='build':r.write_new(OUTPUT,out)
    else:assert r.read(OUTPUT)==out

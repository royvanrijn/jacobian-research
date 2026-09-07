#!/usr/bin/env python3
"""Replay targeted principal witnesses, then reuse the strict-class extractor."""
from pathlib import Path
from math import gcd,prod
import retrospective as r
import reference_partial_class_extraction as extract
import reference_class_targeted_relations as wave

WORK=r.ROOT/'artifacts/local/rank-jump-reference-targeted-class-extraction-v1'
OUTPUT=r.OUT/'rank_jump_reference_targeted_class_extraction_v1.json'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    WORK.mkdir(parents=True,exist_ok=True);inp=r.read(wave.INPUT);out=r.read(wave.OUTPUT)
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    R=PolynomialRing(QQ,'z');f=R(inp['cubic_ascending']);nf=pari.nfinit([pari(f),inp['S_finite']])
    assert str(nf.disc())==inp['field_discriminant']
    inherited=r.read(extract.ACCEPTED);elements=list(inherited['elements_GP']);a=int(inp['fixed_a'])
    w=pari.Mod(pari(R(inp['w_power_basis'])),pari(f));M=inp['sl2_matrix'];c=list(map(int,inp['binary_cubic_descending']))
    bounds=r.read(wave.PROTOCOL)['bounds'];scale=int(inp['slope_scale']);checks=0
    for k,item in enumerate(out['chunks']):
        path=r.ROOT/item['path'];assert r.digest(path.read_bytes())==item['sha256'];chunk=r.read(path)
        target=inp['selection']['targets'][k];assert target==chunk['target'];v1,v2=target['lattice_basis']
        for row in chunk['relations']:
            u,v,m,n=[row[key] for key in ['u','v','m','n']]
            assert target['old_vmax']<v<=bounds['vmax'] and gcd(u,v)==gcd(m,n)==1
            assert (m,n)==(u*v1[0]+v*v2[0],u*v1[1]+v*v2[1])
            assert any(abs(u-s*v//scale)<=1 for s in target['root_slopes_scaled'])
            alpha=a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
            assert [str(pari.lift(alpha).polcoef(i)) for i in range(3)]==row['alpha_ascending']
            N=ZZ(pari.nfeltnorm(nf,alpha));assert str(N)==row['norm']
            assert N==a*a*sum(c[i]*m**(3-i)*n**i for i in range(4))
            factors=row['ideal_factorization'];assert prod(q['p']**(q['f']*q['valuation']) for q in factors)==abs(N)
            assert all(QQ(x).denominator()==1 for x in pari.nfalgtobasis(nf,alpha))
            for q in factors:
                p=q['p'];assert ZZ(p).is_prime(proof=True) and p<=bounds['smooth_bound']
                P=next(P for P in pari.idealprimedec(nf,p) if str(pari.idealhnf(nf,P))==q['hnf'])
                assert int(P[2])==q['e'] and int(P[3])==q['f']
                assert int(pari.idealval(nf,alpha,P))==q['valuation']>0;checks+=1
            assert any(q['p']==target['p'] and q['hnf']==target['hnf'] for q in factors)
            elements.append(str(pari.nfalgtobasis(nf,alpha)))
    combined=WORK/'accepted-elements.json'
    r.write_new(combined,{'basis_GP':list(map(str,nf.nf_get_zk())),'elements_GP':elements,
        'sources':{'seeded_occurrences':len(inherited['elements_GP']),'targeted_occurrences':len(elements)-len(inherited['elements_GP'])}})
    extract.ACCEPTED=combined
    result=extract.compute();result['schema']='rank-jump.reference-targeted-class-extraction.v1'
    result['targeted_ideal_valuation_checks']=checks
    result['bindings'].update(wave.seed.bindings([Path(__file__),wave.PROTOCOL,wave.INPUT,wave.OUTPUT]))
    r.write_new(OUTPUT,result)
    print('PASS targeted witnesses;',result['new_distinct_principal_elements'],'new elements;',
        result['strict_character_rank'],'detected strict directions;',result['positive_endpoint'],flush=True)

if __name__=='__main__':compute()

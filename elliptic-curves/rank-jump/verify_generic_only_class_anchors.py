#!/usr/bin/env python3
"""Independent local-power verification of generic unramified characters."""
import argparse
from pathlib import Path
import retrospective as r
SOURCE=r.OUT/'rank_jump_generic_only_class_anchors_v1.json'
OUTPUT=r.OUT/'rank_jump_generic_only_class_anchors_verification_v1.json'

def compute():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,matrix,pari
    data=r.read(SOURCE);cases=[]
    for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    for case in data['cases']:
        d=r.read(r.ROOT/case['source']);R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending'])
        nf=pari.nfinit([pari(f),d['S_finite']]);roots=f.roots(AA,multiplicities=False)
        gammas=[pari.Mod(pari(R(g['beta_ascending'])),pari(f)) for g in d['generic_classes']]
        checks=0;chars=[]
        for encoded in case['generic_coefficient_masks']:
            mask=int(encoded);beta=pari.Mod(1,pari(f))
            for i,g in enumerate(gammas):
                if mask>>i&1:beta*=g
            assert ZZ(pari.nfeltnorm(nf,beta)).is_square()
            for p in d['S_finite']:
                for P in pari.idealprimedec(nf,p):
                    assert int(pari.idealval(nf,beta,P))%2==0;checks+=1
                    if p==2:
                        assert int(P[3])==1
                        # Q_2(sqrt(5)) remains the unramified quadratic extension
                        # after a totally ramified base change with residue F_2.
                        assert pari.nfislocalpower(nf,P,beta,2)==1 or pari.nfislocalpower(nf,P,5*beta,2)==1
                        checks+=1
            b=R([QQ(pari.lift(beta).polcoef(i)) for i in range(3)])
            assert all(b(x)>0 for x in roots)
            bits=[]
            for a in case['anchors']:
                p=a['p'];v=int(b.change_ring(GF(p))(a['root']));assert v
                bit=int(pow(v,(p-1)//2,p)==p-1)
                assert bit==(int(a['generic_character_mask'])&mask).bit_count()%2
                bits.append(bit)
            chars.append(bits)
        assert matrix(GF(2),chars).rank()==case['ordinary_unramified_dimension']
        cases.append({'source':case['source'],'independent_unramified_characters':len(chars),'local_checks':checks})
    return {'schema':'rank-jump.generic-only-class-anchors-verification.v1','status':'PASS','cases':cases,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),SOURCE]},
        'boundary':'Certifies inherited unramified characters. Does not construct an additional strict class.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS independently verified generic unramified characters')

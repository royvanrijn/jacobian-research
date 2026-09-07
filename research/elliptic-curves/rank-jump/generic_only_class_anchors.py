#!/usr/bin/env python3
"""Ordinary unramified class characters constructed only from generic sections."""
import argparse
from math import prod
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('GENERIC_ONLY_CLASS_ANCHORS_PROTOCOL.json')
INPUTS=[r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json',
        r.OUT/'rank_jump_curve302_strict_constructor_arithmetic_v1.json']
OUTPUT=r.OUT/'rank_jump_generic_only_class_anchors_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-generic-only-class-anchors-v1'

def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}

def compute():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,matrix,pari,prime_range
    policy=r.read(PROTOCOL);pari.allocatemem(64000000,policy['bounds']['pari_stack_bytes'],silent=True)
    cases=[]
    for source in INPUTS:
        d=r.read(source);R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);S=d['S_finite']
        assert f.is_irreducible() and all(ZZ(p).is_prime(proof=True) for p in S)
        remaining=abs(ZZ(f.discriminant()))
        for p in S:
            while remaining%p==0:remaining//=p
        assert remaining==1
        nf=pari.nfinit([pari(f),S]);polys=[R(c['beta_ascending']) for c in d['generic_classes']]
        gammas=[pari.Mod(pari(b),pari(f)) for b in polys];n=len(gammas)
        for b,g,c in zip(polys,gammas,d['generic_classes']):
            assert b.degree()==1 and b[1]<0 and ZZ(-b[1]).is_square()
            assert all(x in ZZ for x in b.list()) and ZZ(b[0]).gcd(ZZ(-b[1]).sqrt())==1
            norm=ZZ(pari.nfeltnorm(nf,g));assert norm==ZZ(c['norm']) and norm.is_square() and norm>0
        constraints=[];valuation_records=[];dyadic=[]
        for p in S:
            for P in pari.idealprimedec(nf,p):
                vals=[int(pari.idealval(nf,g,P)) for g in gammas]
                constraints.append([v%2 for v in vals])
                valuation_records.append({'p':p,'hnf':str(pari.idealhnf(nf,P)),'valuations':vals})
                if p!=2:continue
                e=int(P[2]);residue_degree=int(P[3]);assert residue_degree==1
                pi=pari(2) if e==1 else pari.nfbasistoalg(nf,pari.idealappr(nf,P))
                assert int(pari.idealval(nf,pi,P))==1
                units=[]
                for k in range(1,2*e+1):
                    u=1+pi**k;assert int(pari.idealval(nf,u,P))==0
                    bits=[int(pari.nfhilbert(nf,g,u,P)==-1) for g in gammas]
                    constraints.append(bits);units.append({'k':k,'element_GP':str(u),'hilbert_bits':bits})
                dyadic.append({'hnf':str(pari.idealhnf(nf,P)),'e':e,'f':residue_degree,
                               'uniformizer_GP':str(pi),'unit_filtration':units})
        roots=f.roots(AA,multiplicities=False);assert len(roots)==3
        signs=[[int(b(x)<0) for b in polys] for x in roots];constraints.extend(signs)
        kernel=matrix(GF(2),constraints).right_kernel().basis_matrix()
        masks=[r.pack(list(map(int,row))) for row in kernel];dim=len(masks)
        anchors=[];labels=[]
        for p0 in prime_range(policy['bounds']['anchor_prime_start'],policy['bounds']['anchor_prime_bound']+1):
            p=int(p0)
            if p in S or f.discriminant()%p==0:continue
            for root0 in sorted(f.change_ring(GF(p)).roots(multiplicities=False)):
                root=int(root0);vals=[int(b.change_ring(GF(p))(root0)) for b in polys]
                if not all(vals):continue
                raw=r.pack([int(pow(v,(p-1)//2,p)==p-1) for v in vals])
                label=r.pack([(raw&m).bit_count()%2 for m in masks])
                if r.rank(labels+[label])==len(labels):continue
                Ps=[P for P in pari.idealprimedec(nf,p) if int(P[3])==1
                    and int(pari.idealval(nf,pari.Mod(pari(R.gen()-root),pari(f)),P))>0]
                assert len(Ps)==1
                labels.append(label);anchors.append({'p':p,'root':root,'hnf':str(pari.idealhnf(nf,Ps[0])),
                    'generic_residues':vals,'generic_character_mask':str(raw),'unramified_character_mask':str(label)})
                if len(labels)==dim:break
            if len(labels)==dim:break
        assert len(labels)==dim and r.rank(labels)==dim
        record={'source':str(source.relative_to(r.ROOT)),'generic_dimension':n,
            'ordinary_unramified_dimension':dim,'generic_coefficient_masks':list(map(str,masks)),
            'constraints':constraints,'bad_prime_valuations':valuation_records,'dyadic':dyadic,
            'real_sign_constraints':signs,'anchors':anchors,
            'proof':'Norm-square primitive linear generic Kummer representatives have even valuations outside discriminant support. At odd primes even valuation implies unramified quadratic extension. At dyadic residue degree one, 1+pi^k for 1<=k<=2e generate units modulo squares, since U^(2e+1) consists of squares. Hilbert orthogonality to these units is precisely unramifiedness. Positivity splits all real places. Independent good-prime residue characters certify independence as ordinary class-group characters.'}
        cases.append(record);r.write_new(WORK/('case_%d.json'%len(cases)),record)
        print('CASE',len(cases),'GENERIC',n,'UNRAMIFIED',dim,'ANCHORS',[(a['p'],a['root']) for a in anchors],flush=True)
    return {'schema':'rank-jump.generic-only-class-anchors.v1','status':'PASS',
        'bindings':bindings([Path(__file__),PROTOCOL,*INPUTS]),'cases':cases,
        'additional_strict_class':'NOT_CONSTRUCTED_BY_THIS_CONTROL',
        'boundary':'Inherited ordinary unramified characters and independent ideal anchors only. No exceptional point, historical protected anchor, class-group completeness assertion, or new rank direction.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error,'bindings':bindings([Path(__file__),PROTOCOL,*INPUTS,WORK/'worker.log'])})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()

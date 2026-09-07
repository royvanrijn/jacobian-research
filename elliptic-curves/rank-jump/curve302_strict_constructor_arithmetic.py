#!/usr/bin/env python3
"""Equation-only maximal-order and strict-local input for302 class construction."""
import argparse
from pathlib import Path
from math import prod
import subprocess
import sys
import retrospective as r
from early_relation_pool import elimination

PROTOCOL=Path(__file__).with_name('CURVE302_STRICT_CONSTRUCTOR_ARITHMETIC_PROTOCOL.json')
PARENT=r.OUT/'rank_jump_curve302_new_parent_inputs_v1.json'
FACTORS=r.OUT/'record_prime_factor_proofs_20260904.json'
INPUT=r.OUT/'rank_jump_curve302_strict_constructor_arithmetic_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-curve302-strict-constructor-arithmetic-v1'

def binding(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}

def prepare():
    p=r.read(PARENT);f=r.read(FACTORS)['records']['302']['factorizations']
    r.write_new(INPUT,{'schema':'rank-jump.curve302-strict-constructor-arithmetic-inputs.v1',
        'model':p['specialized_model'],'generic_points':p['specialized_generic_points'],
        'discriminant_factors':f['DISCRIMINANT_FACTORIZATION'],'conductor_factors':f['CONDUCTOR_FACTORIZATION'],
        'bindings':binding([Path(__file__),PROTOCOL,PARENT,FACTORS])})

def compute():
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari,EllipticCurve,GF,matrix
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    pari.allocatemem(64000000,268435456,silent=True)
    d=r.read(INPUT);aa=list(map(QQ,d['model']));assert aa[:3]==[1,1,1]
    E=EllipticCurve(aa);Df=abs(ZZ(E.discriminant()));fac=[(ZZ(p),int(e)) for p,e in d['discriminant_factors']]
    assert prod(p**e for p,e in fac)==Df and all(p.is_prime(proof=True) for p,e in fac)
    S=sorted(int(p) for p,e in d['conductor_factors']);assert S==[int(p) for p,e in fac]
    R=PolynomialRing(QQ,'z');f=R([64*aa[4]+16,16*aa[3]+8,5,1]);assert f.discriminant()==256*E.discriminant()
    nf=pari.nfinit([pari(f),S]);index=ZZ(nf[3]);field_disc=ZZ(nf.disc())
    assert f.discriminant()==index**2*field_disc
    decomposition=[];fieldfac=[];local=[];gammas=[];polys=[]
    for x,y in d['generic_points']:
        x,y=QQ(x),QQ(y);assert E.is_on_curve(x,y)
        xx,yy=4*x,8*y+4*x+4;assert yy**2==f(xx)
        den=ZZ(xx.denominator()).sqrt();assert den in ZZ
        a,b=ZZ(xx*den**2),ZZ(yy*den**3);beta=R([a,-den**2])
        gamma=pari.Mod(pari(beta),pari(f));assert pari.nfeltnorm(nf,gamma)==b*b
        gammas.append(gamma);polys.append(beta)
    for p in S:
        ev=int(abs(field_disc).valuation(p));iv=int(index.valuation(p))
        if ev:fieldfac.append([p,ev])
        ps=pari.idealprimedec(nf,p)
        decomposition.append({'p':p,'field_discriminant_valuation':ev,'order_index_valuation':iv,
            'primes':[{'hnf':str(pari.idealhnf(nf,P)),'ramification_index':int(P[2]),'residue_degree':int(P[3])} for P in ps]})
        L=LocalSquareclasses(nf,p);sig=[list(map(int,L.signature(g))) for g in gammas]
        local.append({'place':p,'signatures':sig,'generic_image_rank':r.rank(list(map(r.pack,sig))),
                      'point_kummer_dimension':L.point_kummer_dimension})
    assert prod(ZZ(p)**e for p,e in fieldfac)==abs(field_disc)
    roots=f.roots(AA,multiplicities=False)
    signs=[[int(beta(x)<0) for x in roots] for beta in polys]
    local.append({'place':'infinity','signatures':signs,'generic_image_rank':r.rank(list(map(r.pack,signs))),
                  'point_kummer_dimension':int(len(roots)==3)})
    locrows=[r.pack([b for c in local for b in c['signatures'][i]]) for i in range(17)]
    localrank,kernel=elimination(locrows)
    assert matrix(GF(2),[[int(v>>i&1) for i in range(max(locrows).bit_length())] for v in locrows]).rank()==localrank
    chars=[0]*17;blocks=[];width=0
    for p in r.primes(r.read(PROTOCOL)['bounds']['proof_prime_bound']):
        if p in S:continue
        ff=f.change_ring(GF(p));rr=ff.roots(multiplicities=False)
        if len(rr)!=3:continue
        vals=[[int(beta.change_ring(GF(p))(x)) for x in rr] for beta in polys]
        if any(v==0 for row in vals for v in row):continue
        for i,row in enumerate(vals):
            for j,v in enumerate(row):chars[i]|=int(pow(v,(p-1)//2,p)==p-1)<<(width+j)
        blocks.append({'prime':p,'roots':list(map(int,rr))});width+=3
    assert r.rank(chars)==17
    return {'schema':'rank-jump.curve302-strict-constructor-arithmetic.v1','status':'PASS',
        'bindings':binding([Path(__file__),PROTOCOL,INPUT,r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py']),
        'cubic_ascending':list(map(str,f.list())),'S_finite':S,'field_discriminant':str(field_disc),
        'field_discriminant_factors':fieldfac,'defining_order_index':str(index),'prime_decomposition':decomposition,
        'maximal_order_basis':list(map(str,nf.nf_get_zk())),
        'generic_classes':[{'beta_ascending':list(map(str,b.list())),'norm':str(pari.nfeltnorm(nf,g))} for b,g in zip(polys,gammas)],
        'local':local,'generic_local_rows':[str(v) for v in locrows],'generic_local_rank':localrank,
        'generic_strict_dimension':len(kernel),'generic_strict_masks':[str(k) for k in kernel],
        'independence_blocks':blocks,'generic_characters':[str(v) for v in chars],'generic_mod2_rank':17,
        'additional_strict_class':'NOT_YET_CONSTRUCTED','boundary':'Complete equation-only constructor arithmetic and inherited strict controls. No new class, Selmer upper bound or rank prediction.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    out=r.read(OUTPUT);print(out['status'],'strict G',out.get('generic_strict_dimension'),'field bits',len(bin(abs(int(out.get('field_discriminant',0)))))-2,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:globals()[a.mode]()

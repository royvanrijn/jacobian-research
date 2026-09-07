#!/usr/bin/env python3
"""Complete v1 unchanged: correct PARI basis input from t_VEC to t_COL."""
import argparse
from pathlib import Path
from itertools import product
from math import gcd
import subprocess
import sys
import retrospective as r
import bounded_gain_reference as ref
import early_relation_pool as pool

PROTOCOL=Path(__file__).with_name('PRIME_SQUARE_CLASS_CONSTRUCTOR_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_prime_square_class_constructor_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_prime_square_class_constructor_v2.json'
WORK=r.ROOT/'artifacts/local/rank-jump-prime-square-class-constructor-v2'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def directions():
    return [v for v in product(range(-2,3),repeat=3) if gcd(*v)==1 and next(x for x in v if x)!=0 and next(x for x in v if x)>0]


def worker():
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    pari.allocatemem(64000000,268435456,silent=True)
    data=r.read(INPUT);R=PolynomialRing(QQ,'z');f=R(data['cubic_ascending']);assert f.discriminant()>0
    nf=pari.nfinit([pari(f),data['S_finite']]);assert str(nf.disc())==data['field_discriminant']
    th=pari.Mod('z',pari(f));zk=list(nf.nf_get_zk());T=matrix(ZZ,3,3,lambda i,j:ZZ(pari.nfelttrace(nf,zk[i]*zk[j])))
    assert T.is_positive_definite();ds=directions();assert len(ds)==49
    enc=lambda v:[str(pari.lift(v).polcoef(i)) for i in range(3)]
    mat=lambda M:[[str(M[i,j]) for j in range(3)] for i in range(3)]
    targets=[]
    for c in data['prime_ideals']:
        dec=pari.idealprimedec(nf,c['p']);P=next(q for q in dec if str(pari.idealhnf(nf,q))==c['hnf'])
        assert int(P[2])==c['e'] and int(P[3])==1
        I=pari.idealpow(nf,P,2);targets.append((f"prime-column-{c['retained_column']}",'candidate',I,{'prime':c['p'],'prime_hnf':c['hnf']}))
    for i,g in enumerate(data['generic_classes'][:3]):
        beta=pari(R(g['beta_ascending']))(th);assert pari.nfeltnorm(nf,beta)==pari(QQ(g['norm']))
        targets.append((f'generic-{i}','positive_control',pari.idealhnf(nf,beta),{'generic_index':i}))
    rows=[]
    for label,kind,I,extra in targets:
        H=matrix(ZZ,3,3,lambda i,j:ZZ(I[i,j]));G=H.transpose()*T*H
        up=pari.qflllgram(pari(G));U=matrix(ZZ,3,3,lambda i,j:ZZ(up[i,j]));assert abs(U.det())==1
        reduced=H*U;target=ZZ(pari.idealnorm(nf,I));assert target.is_square()
        norms=[];hits=[]
        for v in ds:
            coords=reduced*vector(ZZ,v);alpha=pari.nfbasistoalg(nf,pari(list(coords)).Col())
            N=ZZ(pari.nfeltnorm(nf,alpha));norms.append(str(N))
            if abs(N)!=target:continue
            alpha=alpha if N>0 else -alpha
            assert pari.idealhnf(nf,alpha)==I and pari.nfeltnorm(nf,alpha)==target
            hits.append({'direction':list(v),'sign_correction':1 if N>0 else -1,'alpha_ascending':enc(alpha),'norm':str(target)})
        rows.append({'id':label,'kind':kind,**extra,'ideal_hnf':mat(H),'lll_unimodular_matrix':mat(U),
            'ideal_norm':str(target),'tested_norms':norms,'hits':hits})
        print(label,len(hits),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.prime-square-class-constructor.v2','status':'PASS',
        'maximal_order_basis':[enc(x) for x in zk],'trace_gram':mat(T),'directions':[list(v) for v in ds],'rows':rows,
        'bindings':bindings([Path(__file__),Path(__file__).with_name('prime_square_class_constructor.py'),PROTOCOL,INPUT,Path(r.__file__)]),
        'completion_reason':'v1 stopped before enumeration: nfbasistoalg requires t_COL, not t_VEC. Selection, basis reduction, coefficient box and time limit unchanged.',
        'boundary':'Exact principal-square witnesses or bounded generator misses. Additional independence, strict local solubility and CT have not yet been asserted.'})


def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error:r.write_new(OUTPUT,{'schema':'rank-jump.prime-square-class-constructor.v2','status':'UNKNOWN','reason':error})
    result=r.read(OUTPUT);print(result['status'],[(x['id'],len(x['hits'])) for x in result.get('rows',[])],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['worker','capture']);args=p.parse_args();globals()[args.mode]()

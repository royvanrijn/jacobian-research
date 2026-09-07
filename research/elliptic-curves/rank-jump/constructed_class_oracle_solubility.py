#!/usr/bin/env python3
"""Retrospective rational-point witnesses for already frozen constructed covers."""
import argparse
from functools import reduce
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('CONSTRUCTED_CLASS_ORACLE_SOLUBILITY_PROTOCOL.json')
CLASSES=r.OUT/'rank_jump_constructed_class_compaction_v1.json'
REFERENCE=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
ORACLE=r.OUT/'small_conductor_rank22_proof_v1.json'
OUTPUT=r.OUT/'rank_jump_constructed_class_oracle_solubility_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-constructed-class-oracle-solubility-v1'

def compute():
    from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,pari,lcm,gcd
    data=r.read(CLASSES);ref=r.read(REFERENCE);oracle=r.read(ORACLE);policy=r.read(PROTOCOL)
    for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,policy['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    E=EllipticCurve(list(map(QQ,oracle['integral_model'])));assert list(E.a_invariants()[:3])==[1,0,0]
    assert f==R([64*E.a6(),16*E.a4(),1,1])
    points=[E(list(map(QQ,p))) for p in oracle['integral_points']];assert len(points)==22
    polys=[4*QQ(P[0])-R.gen() for P in points]+[R(c['beta_ascending']) for c in data['cases']]
    chars=[0]*len(polys);blocks=[];width=0
    for p in r.primes(policy['bounds']['character_prime_bound']):
        if p==2 or f.discriminant()%p==0 or any(a.denominator()%p==0 for b in polys for a in b):continue
        rr=f.change_ring(GF(p)).roots(multiplicities=False)
        if len(rr)!=3:continue
        vals=[[int(b.change_ring(GF(p))(x)) for x in rr] for b in polys]
        if any(v==0 for row in vals for v in row):continue
        for i,row in enumerate(vals):
            for j,v in enumerate(row):chars[i]|=int(pow(v,(p-1)//2,p)==p-1)<<(width+j)
        blocks.append({'p':p,'roots':list(map(int,rr))});width+=3
    A=matrix(GF(2),[[(c>>j)&1 for j in range(width)] for c in chars[:22]])
    assert A.rank()==22
    sqrt=pari('(nf,b)->{my(y); if(nfeltissquare(nf,b,&y),y,0)}')
    V=PolynomialRing(QQ,['u','v','w','s']);results=[]
    for k,case in enumerate(data['cases']):
        try:mask=A.transpose().solve_right(vector(GF(2),[(chars[22+k]>>j)&1 for j in range(width)]))
        except ValueError:
            results.append({'column':case['column'],'status':'UNKNOWN','reason':'no finite-character match in the oracle span'});continue
        indices=[i for i,b in enumerate(mask) if b];P=E(0)
        for i in indices:P+=points[i]
        assert not P.is_zero();xx=4*QQ(P[0]);yy=8*QQ(P[1])+4*QQ(P[0]);assert yy**2==f(xx)
        beta=pari.Mod(pari(R(case['beta_ascending'])),pari(f));ratio=pari.Mod(pari(xx-R.gen()),pari(f))/beta
        root=sqrt(nf,ratio)
        if root==0:
            results.append({'column':case['column'],'status':'UNKNOWN','reason':'finite-character match is not an exact square'});continue
        xi=pari.nfbasistoalg(nf,root);assert beta*xi**2==pari.Mod(pari(xx-R.gen()),pari(f))
        Nbeta=QQ(pari.nfeltnorm(nf,beta));q=Nbeta.sqrt();assert q in QQ and q*q==Nbeta
        if q*QQ(pari.nfeltnorm(nf,xi))==-yy:xi=-xi
        assert q*QQ(pari.nfeltnorm(nf,xi))==yy
        coords=[QQ(pari.lift(xi).polcoef(i)) for i in range(3)]+[QQ(1)]
        scale=lcm([x.denominator() for x in coords]);integral=[ZZ(scale*x) for x in coords]
        content=gcd(integral);integral=[x//content for x in integral]
        quadrics=[V(q) for q in case['cover_quadrics']];assert all(q(*integral)==0 for q in quadrics)
        assert V(case['x_numerator'])(*integral)/integral[3]**2==xx
        result={'column':case['column'],'status':'RATIONAL_CERTIFIED_RETROSPECTIVELY',
            'oracle_point_indices':indices,'point_on_original_model':[str(P[0]),str(P[1])],
            'point_on_cubic_model':[str(xx),str(yy)],'cubic_field_square_root_ascending':[str(pari.lift(xi).polcoef(i)) for i in range(3)],
            'primitive_cover_point':list(map(str,integral)),
            'proof':'beta*xi^2=4*x(P)-theta exactly; both cover quadrics vanish at the displayed primitive rational point.'}
        results.append(result);r.write_new(WORK/('column_%02d.json'%case['column']),result)
        print('COLUMN',case['column'],'RATIONAL','ORACLE_MASK',indices,flush=True)
    return {'schema':'rank-jump.constructed-class-oracle-solubility.v1','status':'PASS','cases':results,
        'oracle_character_rank':22,'character_blocks':blocks,'rational_classes_certified':sum(x['status']=='RATIONAL_CERTIFIED_RETROSPECTIVELY' for x in results),
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,CLASSES,REFERENCE,ORACLE]},
        'boundary':'Oracle-assisted retrospective solubility evaluation only. The classes were constructed and independently certified before this oracle was admitted. These witnesses are forbidden as prospective selector inputs. No new point search, curve rank or global Selmer upper bound.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()

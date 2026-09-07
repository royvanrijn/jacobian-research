#!/usr/bin/env python3
"""Identically selected inherited pair governing polynomials across a masked panel."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base

PROTOCOL=Path(__file__).with_name('FRESH_GOVERNING_OCTIC_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fresh_governing_octics_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fresh-governing-octics-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,base.INPUT,Path(base.__file__),Path(r.__file__))}


def worker(token):
    from sage.all import QQ,ZZ,PolynomialRing,GF,lcm
    f,points,scale=base.model_data(token);P,Q=points[:2];x,y=P;u,v=Q;c=u-x
    assert c and y and v and not f.discriminant().is_square()
    R=PolynomialRing(QQ,'T');T=R.gen()
    h=T**8-4*(y+v)*T**6+6*c*c*(x+u)*T**4-4*c**3*(v-y)*T**2+c**6
    d=lcm(z.denominator() for z in h.list());H=R([h[i]*d**(8-i) for i in range(9)])
    disc=ZZ(H.discriminant());assert disc
    model=['0','0','0',str(f[1]),str(f[0])];gal=r.galois(model);assert gal['galois_group']=='S3'
    blocks=[]
    for p in r.primes(1999):
        roots=r.roots_at(model[3],model[4],p)
        if roots:blocks.append((p,roots))
        sigs=[r.point_signature(model,list(map(str,point)),blocks) for point in [P,Q]]
        if r.rank(sigs)==2:break
    assert r.rank(sigs)==2
    table=[];excluded=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['prime_bound']):
        values=[d,disc,f.discriminant(),c.numerator(),c.denominator(),*[a.numerator() for a in [y,v]],*[a.denominator() for a in [x,y,u,v]]]
        if any(a%p==0 for a in values):excluded.append(p);continue
        Rp=PolynomialRing(GF(p),'z');fp=Rp(f.list())
        if not fp.is_irreducible():continue
        degrees=sorted(int(g.degree()) for g,m in Rp(H.list()).factor() for _ in range(m))
        assert degrees in ([1,1,3,3],[2,6])
        F=GF(p**3,name='theta',modulus=fp);theta=F.gen();e=p*p+p+1
        a=(F(x)-theta).sqrt();b=(F(u)-theta).sqrt()
        if a**e!=F(y):a=-a
        if b**e!=F(v):b=-b
        assert a**e==F(y) and b**e==F(v)
        norm=(a+b)**e;assert norm and norm**p==norm
        scalar=int(GF(p)(norm));bit=int(pow(scalar,(p-1)//2,p)==p-1)
        assert bit==int(degrees==[2,6])
        table.append({'prime':p,'factor_degrees':degrees,'radical_norm':scalar,'psi':bit})
    assert table
    return {'status':'PASS','token':token,'generic_pair_indices':[0,1],'relative_dimension_mod_generic':0,
       'model_scale':str(scale),'rational_octic_ascending':list(map(str,h.list())),
       'integral_octic_root_scale':str(d),'integral_octic_ascending':list(map(str,H.list())),
       'integral_octic_discriminant':str(disc),'galois':gal,'joint_class_field_degree':96,'governing_field_degree':192,
       'independence_blocks':blocks,'independence_signatures':sigs,'inert_prime_table':table,'excluded_primes':excluded,
       'boundary':'Explicit inherited-pair cochain; no exceptional quotient class. Inert psi is not a CT value without the relevant local twist hypotheses.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for c in r.read(base.INPUT)['cases']:
        token=c['token'];path=WORK/f'{token}.json'
        if not path.exists():
            with (WORK/f'{token}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['worker_seconds'])
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
            if reason:r.write_new(path,{'bindings':bindings(),'token':token,'status':'UNKNOWN','reason':reason})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(token,row['status'],len(row.get('inert_prime_table',[])),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fresh-governing-octics.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');a=p.parse_args()
    if a.mode=='capture':capture()
    else:r.write_new(WORK/f'{a.token}.json',{'bindings':bindings(),**worker(a.token)})

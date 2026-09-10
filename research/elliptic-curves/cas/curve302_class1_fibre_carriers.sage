#!/usr/bin/env sage-python
"""Exact fourteen class1 fibres as quadratic carriers of the original A-fibration.

One frozen prime <=2000 certifies all14 squarefree quartics and all91 pairwise
coprimalities. No carrier enumeration or factorization. Common-cover misses
concern only these fourteen specified fibres.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,prime_range
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'curve302_class1_point_transport_v1.json'
OUTPUT=ART/'curve302_class1_fibre_carriers_v1.json'
WORK=ROOT/'artifacts/local/elliptic-curves/curve302-class1-fibre-carriers-v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    if p.exists():assert read(p)==d,p
    else:
        p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x')as stream:json.dump(d,stream,indent=2,sort_keys=True);stream.write('\n')
def evaluate(c,t):
    ans=QQ(0)
    for v in reversed(c):ans=ans*t+QQ(v)
    return ans
def dec(c,t):return evaluate(c['numerator'],t)/evaluate(c['denominator'],t)
def enc(f):
    return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
def compute():
    paths=[Path(__file__),INPUT]+[ART/('x1092_class1_realization_'+s+'_v1.json')for s in ['trace','rr','equation']]
    bindings={str(p.relative_to(ROOT)):sha(p)for p in paths}
    save(WORK/'protocol.json',{'bindings':bindings,'limits':{'wall_seconds':120,'workers':1,'carriers':14,'pairs':91,'good_prime_bound':2000,'new_parameters':0,'factorizations':0}})
    data=read(INPUT);tr,rr,eq=map(read,paths[2:]);R=PolynomialRing(QQ,'t');F=R.fraction_field();t=R.gen()
    A,B,h,nx,ny=[R(tr[k])for k in ['A','B','h','Nx','Ny']]
    a0,b0,a1,b1=[R(rr[k])for k in ['a0','b0','a1','b1']]
    xp,yp=F(nx/h**2),F(ny/h**3);polys=[];carriers=[]
    for row in data['cases']:
        u=QQ(row['u']);q=R([dec(c,u)for c in eq['quartic_coefficients']])
        sf=R([dec(c,u)for c in eq['radical_square_factor']]);ds=R([dec(c,u)for c in eq['radical_denominator_sqrt']])
        m=F((a1-u*a0)/((u*b0-b1)*h))
        x0=(m*m-xp)/2;x1=F(sf/(2*ds));y0=m*(x0-xp)-yp;y1=m*x1
        assert q.degree()==4
        assert y0*y0+y1*y1*q==x0**3+3*x0*x1*x1*q+A*x0+B
        assert 2*y0*y1==3*x0*x0*x1+x1**3*q+A*x1
        v=QQ(row['quartic_point'][1]);assert v and q(0)==v*v
        at0=lambda f:f.numerator()(0)/f.denominator()(0)
        xs,ys=map(QQ,row['source_short_point'])
        assert at0(x0)+at0(x1)*v==xs and at0(y0)+at0(y1)*v==ys
        rec={'id':row['id'],'u':row['u'],'branch_quartic':list(map(str,q.list())),
             'rational_point_at_A_zero':['0',str(v)],'source_maps':{k:enc(vv)for k,vv in [('x0',x0),('x1',x1),('y0',y0),('y1',y1)]},
             'deck_trace':'P_w','deck_formula':'sigma(P)=P_w-P','native_anti_invariant_generator':'2P-P_w',
             'native_added_rank':1,'full_twist_rank':'UNKNOWN; at least1'}
        save(WORK/(row['id']+'.json'),rec);carriers.append(rec);polys.append(q)
        print('CARRIER',row['id'],'quartic',flush=True)
    witnesses=None
    for prime in prime_range(3,2001):
        p=int(prime)
        if any(c.denominator()%p==0 for q in polys for c in q.list()):continue
        qs=[q.change_ring(GF(p))for q in polys]
        if any(q.degree()!=4 or q.gcd(q.derivative()).degree()!=0 for q in qs):continue
        pairs=[{'i':i,'j':j,'resultant_mod_p':int(qs[i].resultant(qs[j])),'joint_normalization_genus':5}
               for i in range(14)for j in range(i+1,14)]
        if any(r['resultant_mod_p']==0 for r in pairs):continue
        witnesses={'prime':p,'squarefree_resultants':[int(q.resultant(q.derivative()))for q in qs],'pairs':pairs};break
    assert witnesses is not None
    result={'schema':'curve302.class1-fibre-carriers.v1','bindings':bindings,'carriers':carriers,
            'finite_certificate':witnesses,'status':'PASS_FOURTEEN_DISJOINT_GENUS_ONE_CARRIERS',
            'shared_quadratic_squareclasses':0,'genus_zero_or_one_joint_pairs':0,
            'boundary':'All91 pairs of the14 B-fibres through the frozen points have disjoint branch divisors over A, so joint genus5. This does not exclude other multisections, translated representatives, or additional sections on any one twist.'}
    save(OUTPUT,result);print('PASS14 carriers,91 genus5 pairs; witness prime',witnesses['prime'],flush=True)
if __name__=='__main__':compute()

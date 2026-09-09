#!/usr/bin/env sage-python
"""Resolve all final survivors, within the predeclared12-map/25s bound.

Exact generic word maps and rational quadratic equations. The old first-seed
carrier label is the only seed-derived input; no seed coordinate is read.
"""
import hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_low_degree_source_obstruction_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def rec(v):return dict(numerator=list(map(str,v.numerator().list())),denominator=list(map(str,v.denominator().list())))
signal.alarm(25);begun=time.monotonic()
protocol=read(OUT/'incidence-protocol-v2.json');sources=read(OUT/'sources.json')
for name,digest in protocol['inputs'].items():assert sha(ROOT/name)==digest
survivors=[i for i,r in enumerate(sources['sources']) if r['degree']==2]
for p in protocol['primes']:
    row=read(OUT/('prime-%d.json'%p));assert row['incoming']==survivors;survivors=row['survivors']
assert len(survivors)<=protocol['rational_lift_cap']
save('lift-protocol.json',dict(classification='exact old-carrier incidence on all modular survivors',
    sources=survivors,source_cap=12,seconds=25,point_searches=0,
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'incidence-protocol-v2.json',OUT/'sources.json',PARENT,PENCIL,
        *[OUT/('prime-%d.json'%p) for p in protocol['primes']],Path(__file__)]}))
R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
parent=read(PARENT);pencil=read(PENCIL);original=EllipticCurve(K,[dec(a) for a in parent['a_invariants']])
A,B=-original.c4()/48,-original.c6()/864;E=EllipticCurve(K,[A,B]);basis=[]
for row in parent['basis_weierstrass_coordinates']:
    x,y=map(dec,row);basis.append(E([x+original.b2()/12,y+(original.a1()*x+original.a3())/2]))
h=R(pencil['pole_h']);shift=R(pencil['shift']);cx=K(R(pencil['nx']))/h**2;cy=K(R(pencil['ny']))/h**3
zstar=QQ(protocol['carrier_label']);results=[]
for i in survivors:
    word=sources['sources'][i]['word'];P=sum((n*P for n,P in zip(word,basis)),E(0))
    z=(h*(P[1]+cy)/(P[0]-cx)+shift)/(h*h);assert max(z.numerator().degree(),z.denominator().degree())==2
    f=R(z.numerator()-zstar*z.denominator());f*=lcm([c.denominator() for c in f.list()]);f=R(f/gcd([ZZ(c) for c in f.list()]))
    if f.leading_coefficient()<0:f=-f
    roots=[];infinity=f.degree()<2;certificate={}
    if f.degree()==2:
        disc=f[1]**2-4*f[2]*f[0];certificate['discriminant']=str(disc)
        if disc<0:certificate['negative']=True
        else:
            n,d=ZZ(disc.numerator()),ZZ(disc.denominator());a,b=n.isqrt(),d.isqrt()
            certificate.update(numerator=str(n),denominator=str(d),numerator_floor_sqrt=str(a),denominator_floor_sqrt=str(b))
            if a*a==n and b*b==d:
                roots=sorted(set([(-f[1]+sgn*QQ(a)/b)/(2*f[2]) for sgn in [-1,1]]))
    elif f.degree()==1:roots=[-f[0]/f[1]]
    elif not f:raise ArithmeticError('constant map cannot equal this smooth carrier')
    row=dict(source_index=i,word=word,z_map=rec(z),section=[rec(v) for v in P.xy()],
        primitive_incidence_polynomial=list(map(str,f.list())),square_certificate=certificate,
        rational_preimages=list(map(str,roots)),infinity_preimage=infinity,
        status='RATIONAL_INCIDENCE_REQUIRES_CLASS_TEST' if roots or infinity else 'EXACT_NO_RATIONAL_SOURCE_POINT')
    save('lift-%d.json'%i,row);results.append(row)
    print(i,row['status'],flush=True)
save('lift-summary.json',dict(status='PASS_ALL_SURVIVORS_EXACTLY_DECIDED',sources=survivors,
    counts={s:sum(r['status']==s for r in results) for s in sorted({r['status'] for r in results})},
    inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'lift-protocol.json',*[OUT/('lift-%d.json'%i) for i in survivors]]}))
print('seconds',round(time.monotonic()-begun,3),flush=True)

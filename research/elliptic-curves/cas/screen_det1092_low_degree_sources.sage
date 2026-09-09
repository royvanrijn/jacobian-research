#!/usr/bin/env sage-python
"""Bounded modular incidence proof on the complete degree-two source set.

Retrospective first-seed carrier only. Reuses the old twelve-prime pool in
its old order. No point search, parameter scan, later seed or map expansion
over Q. Each prime is a separate immutable checkpoint with a25s cap.
"""
import argparse,hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_low_degree_source_obstruction_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
OLD=ART/'det1092_signed_source_orbits_v1/protocol.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,row):
    path=OUT/name
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as s:json.dump(row,s,indent=2,sort_keys=True);s.write('\n')
def freeze():
    old=read(OLD);sources=read(OUT/'sources.json')
    assert sources['status']=='PASS_COMPLETE_LOW_DEGREE_SOURCE_ELLIPSOID'
    save('incidence-protocol.json',dict(classification='retrospective exact source-incidence sieve, not a prospective seed selection',
        carrier_label=old['carrier_label'],primes=old['limits']['prime_pool'],
        seconds_per_prime=25,source_count=sources['counts']['2'],rational_lift_cap=12,
        rule='Test every degree-two source not already excluded by an earlier prime. Use the complete source ellipsoid and the unchanged old prime order; retain all degree-drop/chart failures as UNKNOWN. No Q map reconstruction in this stage.',
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,PENCIL,OLD,OUT/'geometry-protocol.json',OUT/'sources.json',Path(__file__)]}))
def run(prime):
    protocol=read(OUT/'incidence-protocol.json')
    for name,digest in protocol['inputs'].items():assert sha(ROOT/name)==digest
    assert prime in protocol['primes'];earlier=protocol['primes'][:protocol['primes'].index(prime)]
    source=read(OUT/'sources.json');active=[i for i,r in enumerate(source['sources']) if r['degree']==2]
    for p in earlier:
        prev=read(OUT/('prime-%d.json'%p));assert prev['incoming']==active;active=prev['survivors']
    parent=read(PARENT);pencil=read(PENCIL);base=PolynomialRing(QQ,'t');K=base.fraction_field()
    def dec(v):return K(base(v['numerator']))/base(v['denominator'])
    original=EllipticCurve(K,[dec(v) for v in parent['a_invariants']])
    Aq,Bq=-original.c4()/48,-original.c6()/864
    fp=GF(prime);R=PolynomialRing(fp,'t');F=R.fraction_field()
    def gauss(v):
        v=K(v)
        if not v:return 1,F(0)
        num,den=v.numerator(),v.denominator()
        a=min(c.valuation(prime) for c in num.list() if c)
        b=min(c.valuation(prime) for c in den.list() if c)
        val=int(a-b)
        if val<0:return val,None
        if val>0:return val,F(0)
        return val,F(R([fp(c/QQ(prime)**a) for c in num.list()]))/R([fp(c/QQ(prime)**b) for c in den.list()])
    record=dict(prime=prime,incoming=active,protocol_sha256=sha(OUT/'incidence-protocol.json'))
    va,A=gauss(Aq);vb,B=gauss(Bq)
    if va<0 or vb<0 or not (4*A**3+27*B**2):
        record.update(status='UNKNOWN_BAD_GENERIC_REDUCTION',survivors=active,trials=[]);save('prime-%d.json'%prime,record);return
    hq=base(pencil['pole_h']);shiftq=base(pencil['shift'])
    cxq=K(base(pencil['nx']))/hq**2;cyq=K(base(pencil['ny']))/hq**3
    converted=[gauss(v) for v in [hq,shiftq,cxq,cyq,QQ(protocol['carrier_label'])]]
    if any(n<0 for n,v in converted) or not converted[0][1]:
        record.update(status='UNKNOWN_PENCIL_GAUSS_CHART',survivors=active,trials=[]);save('prime-%d.json'%prime,record);return
    h,shift,cx,cy,zstar=[v for n,v in converted]
    E=EllipticCurve(F,[A,B]);basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,row);X=x+original.b2()/12;Y=y+(original.a1()*x+original.a3())/2
        vx,xx=gauss(X);vy,yy=gauss(Y)
        if vx<0 or vy<0:
            assert vx<0 and vy<0 and 3*vx==2*vy;basis.append(E(0))
        else:basis.append(E([xx,yy]))
    cache={(0,)*17:E(0)}
    def point(word):
        word=tuple(word)
        if word not in cache:
            i=next(i for i,n in enumerate(word) if n);sign=1 if word[i]>0 else -1
            shorter=list(word);shorter[i]-=sign;cache[word]=point(shorter)+sign*basis[i]
        return cache[word]
    trials=[];survivors=[]
    for i in active:
        v=source['sources'][i]['word'];P=point(v);row=dict(source_index=i)
        if P.is_zero():row['status']='UNKNOWN_REDUCED_SECTION_IS_O'
        elif P[0]==cx:row['status']='UNKNOWN_PENCIL_DENOMINATOR'
        else:
            z=(h*(P[1]+cy)/(P[0]-cx)+shift)/(h*h)
            num,den=z.numerator(),z.denominator();degree=max(num.degree(),den.degree())
            row.update(map_numerator=list(map(int,num.list())),map_denominator=list(map(int,den.list())),degree=int(degree))
            if degree!=2:row['status']='UNKNOWN_MAP_DEGREE_DROP'
            else:
                f=num-zstar*den;disc=f[1]**2-4*f[2]*f[0]
                row.update(incidence_coefficients=[int(f[j]) for j in range(3)],discriminant=int(disc))
                row['status']='EXCLUDED_NO_PROJECTIVE_ROOT' if f[2] and disc and not disc.is_square() else 'UNKNOWN_LOCAL_ROOT'
        if row['status']!='EXCLUDED_NO_PROJECTIVE_ROOT':survivors.append(i)
        trials.append(row)
    record.update(status='PASS_BOUNDED_MODULAR_SOURCE_INCIDENCE',survivors=survivors,trials=trials,
        cache_size=len(cache),status_counts={s:sum(r['status']==s for r in trials) for s in sorted({r['status'] for r in trials})})
    save('prime-%d.json'%prime,record)
    print('prime',prime,'incoming',len(active),'remaining',len(survivors),record['status_counts'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','run']);p.add_argument('--prime',type=int);args=p.parse_args()
    signal.alarm(25);begun=time.monotonic()
    if args.action=='freeze':freeze()
    else:run(args.prime)
    print('seconds',round(time.monotonic()-begun,3),flush=True)

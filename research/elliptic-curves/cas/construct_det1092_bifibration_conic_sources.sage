#!/usr/bin/env sage-python
"""Three existing rational conics as sources for both elliptic fibrations.

Separate generic construction, retrospective first-carrier incidence, and
one fixed equation-only specialization on the uniquely least alternate
degree source. No new curve selection or point search. Each action <=25s.
"""
import argparse,hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd,lcm
from rational_root_lattice_certificate import root_certificate
from split_seed_descent import build_frame,Classifier,point_record
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_bifibration_conic_sources_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
CONICS=ART/'det1092_rational_bisection_index_v1';OLD=ART/'det1092_signed_source_orbits_v1/protocol.json'
MASKS=[8044,47755,103186]
R=PolynomialRing(QQ,'u');K=R.fraction_field();u=R.gen()
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(v,T):return R(v['numerator'])(T)/R(v['denominator'])(T)
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
def rec(v):
    v=K(v);return dict(numerator=list(map(str,v.numerator().list())),denominator=list(map(str,v.denominator().list())))
def save(name,row):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(row,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    parent=read(PARENT);G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0]);rows=[]
    for mask in MASKS:
        a=vector(ZZ,read(CONICS/('orbit-%d.json'%mask))['word']);assert a*G*a==10
        rows.append(dict(mask=mask,alternate_degree=int(8-w*G*a)))
    selected=min(rows,key=lambda r:r['alternate_degree']);assert sum(r['alternate_degree']==selected['alternate_degree'] for r in rows)==1
    sources=[PARENT,PENCIL,*[CONICS/('orbit-%d.json'%m) for m in MASKS],Path(__file__),
             Path(__file__).with_name('split_seed_descent.py'),Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage'),
             Path(__file__).with_name('rational_root_lattice_certificate.py')]
    save('generic-protocol.json',dict(classification='generic conic maps and fixed one-address construction',
        sources=rows,selected=selected,parameter='0',seconds_per_action=25,prime_cap=1009,halving_steps=8,
        independent_nonhalving_proof_prime_cap=257,point_searches=0,parameter_replacements=0,
        rule='Use all three existing generic conics. For the one specialization choose the unique least alternate-degree source and its existing rational parameter u=0. No control or exceptional data in this choice.',
        inputs={str(p.relative_to(ROOT)):sha(p) for p in sources}))
def protocol():
    p=read(OUT/'generic-protocol.json')
    for name,h in p['inputs'].items():assert sha(ROOT/name)==h
    return p
def construct(mask):
    p=protocol();source=read(CONICS/('orbit-%d.json'%mask));parent=read(PARENT);pencil=read(PENCIL)
    T,W=dec(source['base_map']),dec(source['conic_ordinate'])
    assert max(T.numerator().degree(),T.denominator().degree())==2 and W*W==R(source['q'])(T)
    x0,x1,y0,y1=[R(v)(T) for v in source['elliptic_quadratic_maps']]
    x,y=x0+x1*W,y0+y1*W
    old=EllipticCurve(K,[val(a,T) for a in parent['a_invariants']]);old([x,y])
    X,Y=x+old.b2()/12,y+(old.a1()*x+old.a3())/2
    h=R(pencil['pole_h'])(T);shift=R(pencil['shift'])(T)
    cx=R(pencil['nx'])(T)/h**2;cy=R(pencil['ny'])(T)/h**3
    Z=(h*(Y+cy)/(X-cx)+shift)/(h*h)
    expected=next(r['alternate_degree'] for r in p['sources'] if r['mask']==mask)
    assert max(Z.numerator().degree(),Z.denominator().degree())==expected
    save('map-%d.json'%mask,dict(status='PASS_EXACT_BIFIBRATION_SOURCE_MAP',mask=mask,
        original_degree=2,alternate_degree=expected,T=rec(T),conic_ordinate=rec(W),
        original_point=[rec(x),rec(y)],short_point=[rec(X),rec(Y)],Z=rec(Z),
        source_sha256=sha(CONICS/('orbit-%d.json'%mask)),generic_protocol_sha256=sha(OUT/'generic-protocol.json')))
    print(mask,'bidegree',(2,expected),flush=True)
def incidence(mask):
    protocol();old=read(OLD)
    save('incidence-protocol.json',dict(classification='retrospective first-seed carrier incidence only',
        carrier_label=old['carrier_label'],primes=old['limits']['prime_pool'],seconds=25,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OLD,*[OUT/('map-%d.json'%m) for m in MASKS]]}))
    p=read(OUT/'incidence-protocol.json');mp=read(OUT/('map-%d.json'%mask));Z=dec(mp['Z'])
    f=R(Z.numerator()-QQ(p['carrier_label'])*Z.denominator());infinity=f.degree()<mp['alternate_degree']
    f*=lcm([c.denominator() for c in f.list()]);f=R(f/gcd([ZZ(c) for c in f.list()]))
    if f.leading_coefficient()<0:f=-f
    roots=[];g=f
    if not g[0]:
        roots.append('0')
        while g.degree()>0 and not g[0]:g=R(g/u)
    proof=root_certificate(g,p['primes']) if g.degree()>0 else None
    if proof:roots+=proof['rational_roots']
    status='RATIONAL_SOURCE_POINT_REQUIRES_CLASS_TEST' if roots or infinity else 'EXACT_NO_RATIONAL_SOURCE_POINT'
    save('incidence-%d.json'%mask,dict(mask=mask,status=status,polynomial=list(map(str,f.list())),
        rational_preimages=roots,infinity_preimage=infinity,root_proof=proof,
        map_sha256=sha(OUT/('map-%d.json'%mask)),incidence_protocol_sha256=sha(OUT/'incidence-protocol.json')))
    print(mask,status,flush=True)
def seed():
    p=protocol();mask=p['selected']['mask'];mp=read(OUT/('map-%d.json'%mask));parameter=QQ(p['parameter'])
    T=dec(mp['T']);tau=T(parameter);parent=read(PARENT)
    old=EllipticCurve(QQ,[val(a,tau) for a in parent['a_invariants']]);E=EllipticCurve(QQ,[-old.c4()/48,-old.c6()/864])
    basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=[val(v,tau) for v in row];basis.append(E([x+old.b2()/12,y+(old.a1()*x+old.a3())/2]))
    frame=build_frame(E,basis,p['prime_cap']);save('seed-generic-frame.json',frame)
    # Candidate evaluation follows the generic-only certificate checkpoint.
    P=E([dec(v)(parameter) for v in mp['short_point']]);result=Classifier(frame).classify(P,p['halving_steps'])
    result.update(mask=mask,parameter=str(parameter),original_parameter=str(tau),
        alternate_parameter=str(dec(mp['Z'])(parameter)),
        frame_sha256=sha(OUT/'seed-generic-frame.json'),map_sha256=sha(OUT/('map-%d.json'%mask)),
        generic_protocol_sha256=sha(OUT/'generic-protocol.json'),
        boundary='Original-fibre independence only. Generic alternate independence does not certify this particular alternate specialization.')
    save('seed-result.json',result);print('fixed seed',mask,result['status'],result.get('reason'),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['freeze','map','incidence','seed']);p.add_argument('--mask',type=int,choices=MASKS);a=p.parse_args()
    signal.alarm(25);begun=time.monotonic()
    if a.action=='freeze':freeze()
    elif a.action=='map':construct(a.mask)
    elif a.action=='incidence':incidence(a.mask)
    else:seed()
    print('seconds',round(time.monotonic()-begun,3),flush=True)

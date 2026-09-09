#!/usr/bin/env sage-python
"""Two frozen degree6 maps from a short alternate translating section.

Checkpoint the generic section, then each map; only afterward evaluate the
old302/control addresses. Every action has a25-second cap, no point search.
"""
import argparse,hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd,lcm
from rational_root_lattice_certificate import root_certificate
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_short_alternate_translations_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
SOURCE=ART/'det1092_bifibration_conic_sources_v1/map-47755.json'
R=PolynomialRing(QQ,'u');u=R.gen();K=R.fraction_field()
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
def rec(v):
    v=K(v);return dict(numerator=list(map(str,v.numerator().list())),denominator=list(map(str,v.denominator().list())))
def deg(v):return max(v.numerator().degree(),v.denominator().degree())
def save(name,row):
    path=OUT/name
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(row,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    g=read(OUT/'geometry-v2.json');assert g['status']=='PASS_GENERIC_MINIMUM_IMAGE_SET'
    assert len(g['selected_classes'])==2 and g['minimum_new_original_degree']==6
    save('map-protocol-v2.json',dict(classification='generic-only rational maps, no specialized inputs',
        cases=[dict(index=i,**row) for i,row in enumerate(g['selected_classes'])],
        limits=dict(seconds_per_action=25,maps=2,point_searches=0,target_inputs=0),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,PENCIL,SOURCE,OUT/'protocol-v2.json',OUT/'geometry-v2.json',Path(__file__)]}))
def protocol():
    p=read(OUT/'map-protocol-v2.json')
    for name,h in p['inputs'].items():assert sha(ROOT/name)==h
    return p
def section():
    p=protocol();g=read(OUT/'geometry-v2.json');words={tuple(r['word']) for r in g['short_sections']};assert len(words)==1
    ww=list(words.pop());parent=read(PARENT);G=matrix(ZZ,parent['generic_height_gram'])
    old=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);basis=[old([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    # Arithmetic order only: minimize the generic norm of each partial sum.
    todo=[(i,1 if n>0 else -1) for i,n in enumerate(ww) for j in range(abs(n))]
    word=vector(ZZ,17);point=old(0);order=[]
    while todo:
        def score(item):
            i,s=item;v=vector(ZZ,list(word));v[i]+=s;return (v*G*v,i,s)
        item=min(todo,key=score);todo.remove(item);i,s=item;word[i]+=s
        point+=s*basis[i];order.append(dict(index=i,sign=s,partial_norm=int(word*G*word)))
    assert list(word)==ww and point
    X=point[0]+old.b2()/12;Y=point[1]+(old.a1()*point[0]+old.a3())/2
    pencil=read(PENCIL);h=R(pencil['pole_h']);shift=R(pencil['shift'])
    cx=R(pencil['nx'])/h**2;cy=R(pencil['ny'])/h**3
    m=(Y+cy)/(X-cx);z=(m+shift/h)/h;W=(2*X+cx-m*m)/h
    assert deg(z)==1 and W*W==sum(R(v)(z)*u**i for i,v in enumerate(pencil['quartic_t_coefficients_in_z']))
    N,D=z.numerator(),z.denominator();T=(N[0]-u*D[0])/(u*D[1]-N[1]);assert z(T)==u
    save('short-section.json',dict(status='EXACT_GENERIC_SHORT_ALTERNATE_SECTION',word=ww,
        original_point=[rec(v) for v in point[:2]],short_point=[rec(X),rec(Y)],
        z_of_t=rec(z),W_of_t=rec(W),t_of_z=rec(T),W_of_z=rec(W(T)),
        arithmetic_order=order,map_protocol_sha256=sha(OUT/'map-protocol-v2.json')))
    print('short section, original norm',word*G*word,'degree z(t)',deg(z),flush=True)
def maps(index):
    p=protocol();case=p['cases'][index];assert len(case['translators'])==1
    sign=case['translators'][0]['sign'];src=read(SOURCE);section=read(OUT/'short-section.json');pencil=read(PENCIL)
    tsource=dec(src['T']);z=dec(src['Z']);X,Y=map(dec,src['short_point'])
    h=R(pencil['pole_h'])(tsource);cx=R(pencil['nx'])(tsource)/h**2;cy=R(pencil['ny'])(tsource)/h**3
    slope=(Y+cy)/(X-cx);Wsource=(2*X+cx-slope*slope)/h
    t0=dec(pencil['sections'][0]['t_of_z'])(z);s=dec(pencil['sections'][0]['W_of_z'])(z)
    S=PolynomialRing(K,'v');v=S.gen();f=S([R(row)(z) for row in pencil['quartic_t_coefficients_in_z']])
    q0,q1,q2,q3,q4=f(v+t0).list();assert q0==s*s and s
    J=EllipticCurve(K,[0,q2,0,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4])
    def forward(t,W):
        d=t-t0;assert d
        x=(2*s*(W+s)+q1*d)/(d*d)
        y=((x*x-4*s*s*q4)*d-q1*x-2*s*s*q3)/(2*s)
        return J([x,y])
    Q=forward(tsource,Wsource);Tq=dec(section['t_of_z'])(z);Wq=dec(section['W_of_z'])(z);Qshort=forward(Tq,Wq)
    target=Q+sign*Qshort;assert target
    x,y=target[:2];d=(2*s*y+q1*x+2*s*s*q3)/(x*x-4*s*s*q4)
    T=t0+d;W=x*d*d/(2*s)-s-q1*d/(2*s)
    assert W*W==f(T) and deg(T)==case['original_degree']==6
    save('map-%d.json'%index,dict(status='EXACT_GENERIC_SHORT_TRANSLATION_MAP',index=index,sign=sign,
        T=rec(T),W=rec(W),Z=rec(z),original_degree=int(deg(T)),alternate_degree=int(deg(z)),
        source_sha256=sha(SOURCE),section_sha256=sha(OUT/'short-section.json'),map_protocol_sha256=sha(OUT/'map-protocol-v2.json')))
    print('map',index,'sign',sign,'degree',deg(T),flush=True)
def incidence(index,label):
    protocol();old=read(ART/'det1092_rational_bisection_index_v1/controls.json')
    cases=old['rows'][0]['cases'];pool=read(ART/'det1092_pencil_multiples_v2/protocol.json')['ramification_certificate_primes']
    save('incidence-protocol.json',dict(classification='old nine-address retrospective incidence, after generic maps',
        cases=[dict(label=r['label'],parameter=r['parameter']) for r in cases],prime_pool=pool,seconds_per_case=25,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'map-0.json',OUT/'map-1.json',
            ART/'det1092_rational_bisection_index_v1/controls.json',ART/'det1092_pencil_multiples_v2/protocol.json',
            Path(__file__).with_name('rational_root_lattice_certificate.py')]}))
    ip=read(OUT/'incidence-protocol.json');case=next(r for r in ip['cases'] if r['label']==label)
    mp=read(OUT/('map-%d.json'%index));T=dec(mp['T']);f=T.numerator()-QQ(case['parameter'])*T.denominator()
    infinity=bool(f.degree()<mp['original_degree']);f*=lcm([v.denominator() for v in f]);f=R(f/gcd(list(f)))
    if f.leading_coefficient()<0:f=-f
    roots=[];g=f
    if g[0]==0:
        roots.append('0')
        while g.degree()>0 and not g[0]:g=R(g/u)
    proof=root_certificate(g,pool) if g.degree()>0 else None
    if proof:roots+=proof['rational_roots']
    status='RATIONAL_INCIDENCE_REQUIRES_ADMISSION' if roots or infinity else 'EXACT_NO_RATIONAL_INCIDENCE'
    save('incidence-%d-%s.json'%(index,label),dict(status=status,index=index,case=case,
        polynomial=list(map(str,f.list())),infinity_preimage=infinity,rational_preimages=roots,root_proof=proof,
        map_sha256=sha(OUT/('map-%d.json'%index)),incidence_protocol_sha256=sha(OUT/'incidence-protocol.json')))
    print(index,label,status,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','section','map','incidence']);parser.add_argument('--index',type=int,choices=[0,1]);parser.add_argument('--label');args=parser.parse_args()
    signal.alarm(25);start=time.monotonic()
    if args.action=='freeze':freeze()
    elif args.action=='section':section()
    elif args.action=='map':maps(args.index)
    else:incidence(args.index,args.label)
    print('seconds',round(time.monotonic()-start,3),flush=True)

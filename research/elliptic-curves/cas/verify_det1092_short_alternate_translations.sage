#!/usr/bin/env sage-python
"""Independent short-translation geometry, exact maps, and18 incidences.

Use --action geometry, map --index0/1, incidence, then assemble. Each25s.
No constructor imports or point searches; explicit polynomial chord checks
replace the constructor's elliptic-curve addition for the map verification.
"""
import argparse,hashlib,json,runpy,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,identity_matrix,block_diagonal_matrix,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_short_alternate_translations_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
SOURCE=ART/'det1092_bifibration_conic_sources_v1/map-47755.json'
R=PolynomialRing(QQ,'u');u=R.gen();K=R.fraction_field()
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
def deg(v):return max(v.numerator().degree(),v.denominator().degree())
def hashes(row):
    for name,h in row['inputs'].items():assert sha(ROOT/name)==h
def save(name,row):
    path=OUT/name
    if path.exists():assert read(path)==row,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(row,f,indent=2,sort_keys=True);f.write('\n')
def preflight():
    for name in ['protocol-v2.json','map-protocol-v2.json','incidence-protocol.json']:hashes(read(OUT/name))
    save('replay-protocol.json',dict(checker_sha256=sha(Path(__file__)),seconds_per_action=25,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol-v2.json',OUT/'geometry-v2.json',
            OUT/'map-protocol-v2.json',OUT/'short-section.json',SOURCE,PARENT,PENCIL,
            OUT/'map-0.json',OUT/'map-1.json',OUT/'incidence-protocol.json',
            Path(__file__).with_name('verify_det1092_bifibration_conic_sources.sage')]}))
def geometry():
    parent=read(PARENT);G=matrix(QQ,parent['generic_height_gram'])
    NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G);I=identity_matrix(QQ,19)
    frame=read(ART/'det1092_genus1_picard_image_v1/frame.json')
    D=vector(QQ,frame['D']);O=vector(QQ,frame['new_zero_section']);F=I.column(1)
    dot=lambda a,b:a*NS*b
    vertical=[vector(QQ,r['NS_coordinates']) for r in frame['vertical_old_sections']]
    roots=[v for v in vertical if dot(v,O)==0];assert len(roots)==5
    assert matrix(QQ,[[dot(a,b) for b in roots] for a in roots])==-2*identity_matrix(QQ,5)
    raw=read(ART/'det1092_low_degree_source_obstruction_v1/sources.json')
    prior=read(ART/'det1092_low_degree_source_obstruction_v1/independent-replay.json')
    assert prior['status'].startswith('PASS_')
    data=read(OUT/'geometry-v2.json');assert data['protocol_sha256']==sha(OUT/'protocol-v2.json')
    sources=[(i,r) for i,r in enumerate(raw['sources']) if r['degree']==1];assert len(sources)==130
    histogram={};records=[]
    for index,row in sources:
        v=vector(QQ,row['word']);S=vector(QQ,[1,v*G*v/2,*v]);assert dot(S,D)==1
        bits=[dot(S,r) for r in roots];assert set(bits)<={0,1}
        h=QQ(0) if S==O else 4+2*dot(S,O)-sum(bits)/2
        histogram[str(h)]=histogram.get(str(h),0)+1
        records.append((h,index,S,bits,row['word']))
    assert histogram==data['height_histogram']==read(OUT/'failure-v1.json')['height_histogram']
    assert '3/2' not in histogram
    minimum=min(h for h,index,S,bits,word in records if h>0);assert minimum==2
    selected=[r for r in records if r[0]==minimum];assert len(selected)==1
    h,index,S,bits,word=selected[0]
    assert word==read(OUT/'short-section.json')['word'] and index==126
    # Independent explicit Shioda formula, not a trivial-space projector.
    phi=S-O-(dot(S,O)+2)*D+sum((b*r/2 for b,r in zip(bits,roots)),vector(QQ,19))
    assert all(dot(phi,v)==0 for v in [D,O,*roots]) and -dot(phi,phi)==h
    perpendicular=matrix(QQ,[D,O,*roots])*NS
    perps=perpendicular.right_kernel().basis()
    C=vector(QQ,[2,4]+read(ART/'det1092_rational_bisection_index_v1/orbit-47755.json')['word'])
    actions={}
    for row in data['short_sections']:
        assert row['source_index']==index and row['word']==word and QQ(row['height'])==h
        sign=row['sign'];A=matrix(QQ,row['action']);actions[sign]=A
        assert A.transpose()*NS*A==NS and A*D==D and all(v in ZZ for v in A.list())
        if sign==1:
            assert A*O==S
            for r,b in zip(roots,bits):assert A*r==(r if b==0 else D-r)
            for v in perps:assert A*v==v-dot(v,phi)*D
        image=A*C
        assert list(image)==row['image_NS'] and dot(image,F)==row['original_degree']==6
    assert actions[-1]*actions[1]==I
    # All-integer degree formula via a nilpotent identity on F, not a fit.
    polynomials=[]
    oldA=matrix(QQ,read(ART/'det1092_two_fibration_action_v1/action.json')['alternate_B_translation'])
    for label,A,coefficient in [('short',actions[1],4),('old_B',oldA,19)]:
        delta=A.inverse()-I
        assert delta**3*F==0 and dot(C,F)==2
        assert dot(C,delta*F)==coefficient and dot(C,delta**2*F)==2*coefficient
        Q=A*O;qbits=[dot(Q,r) for r in roots]
        height=4+2*dot(Q,O)-sum(qbits)/2
        assert 2*height==coefficient
        polynomials.append(dict(translator=label,height=str(height),constant=2,
            linear=0,quadratic=coefficient,all_integer_n=True))
    save('replay-geometry.json',dict(status='PASS_SHORT_SECTION_ACTION_AND_ALL_INTEGER_DEGREE_LAW',
        common_sections=130,height_histogram=histogram,minimum_in_common_section_set='2',
        unique_section_word=word,component_bits=list(map(int,bits)),degree_laws=polynomials,
        no_height_three_halves_in_common_set=True,
        boundary='No global minimum assertion for the full alternate MW12 lattice, and no new search policy.',
        replay_protocol_sha256=sha(OUT/'replay-protocol.json')))
    print('PASS_SHORT_SECTION_ACTION_AND_ALL_INTEGER_DEGREE_LAW',polynomials,flush=True)
def check_section():
    parent=read(PARENT);pencil=read(PENCIL);data=read(OUT/'short-section.json')
    assert data['map_protocol_sha256']==sha(OUT/'map-protocol-v2.json')
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);A=-E.c4()/48;B=-E.c6()/864
    basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,row);basis.append((x+E.b2()/12,y+(E.a1()*x+E.a3())/2))
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;xx,yy=Q
        if x==xx and y==-yy:return None
        m=(yy-y)/(xx-x) if x!=xx else (3*x*x+A)/(2*y)
        xnew=m*m-x-xx;return (xnew,m*(x-xnew)-y)
    P=None
    for n,Q in reversed(list(zip(data['word'],basis))):
        for unused in range(abs(n)):P=add(P,(Q[0],Q[1] if n>0 else -Q[1]))
    X,Y=P;assert list(P)==[dec(v) for v in data['short_point']]
    x,y=map(dec,data['original_point']);assert X==x+E.b2()/12 and Y==y+(E.a1()*x+E.a3())/2
    assert Y*Y==X**3+A*X+B
    h=R(pencil['pole_h']);cx=R(pencil['nx'])/h**2;cy=R(pencil['ny'])/h**3
    slope=(Y+cy)/(X-cx);Z=(slope+R(pencil['shift'])/h)/h;W=(2*X+cx-slope*slope)/h
    assert Z==dec(data['z_of_t']) and W==dec(data['W_of_t']) and deg(Z)==1
    T=dec(data['t_of_z']);assert Z(T)==u and dec(data['W_of_z'])==W(T)
    return data
def maps(index):
    section=check_section();mp=read(OUT/('map-%d.json'%index));src=read(SOURCE);pencil=read(PENCIL)
    assert mp['source_sha256']==sha(SOURCE) and mp['section_sha256']==sha(OUT/'short-section.json')
    assert mp['map_protocol_sha256']==sha(OUT/'map-protocol-v2.json')
    assert mp['sign']==(-1 if index==0 else 1)
    T=dec(mp['T']);W=dec(mp['W']);z=dec(mp['Z'])
    assert z==dec(src['Z']) and deg(T)==mp['original_degree']==6 and deg(z)==2
    ts=dec(src['T']);Xs,Ys=map(dec,src['short_point'])
    h=R(pencil['pole_h'])(ts);cx=R(pencil['nx'])(ts)/h**2;cy=R(pencil['ny'])(ts)/h**3
    slope=(Ys+cy)/(Xs-cx);Ws=(2*Xs+cx-slope*slope)/h
    t0=dec(pencil['sections'][0]['t_of_z'])(z);s=dec(pencil['sections'][0]['W_of_z'])(z)
    V=PolynomialRing(K,'v');v=V.gen();f=V([R(row)(z) for row in pencil['quartic_t_coefficients_in_z']])
    q0,q1,q2,q3,q4=f(v+t0).list();assert q0==s*s and W*W==f(T)
    a4=q1*q3-4*s*s*q4;a6=s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4
    def forward(t,W):
        d=t-t0;assert d
        xx=(2*s*(W+s)+q1*d)/(d*d)
        yy=((xx*xx-4*s*s*q4)*d-q1*xx-2*s*s*q3)/(2*s)
        assert yy*yy==xx**3+q2*xx*xx+a4*xx+a6
        return xx,yy
    xs,ys=forward(ts,Ws)
    xq,yq=forward(dec(section['t_of_z'])(z),dec(section['W_of_z'])(z));yq*=mp['sign']
    xt,yt=forward(T,W);assert xq!=xs
    chord=(yq-ys)/(xq-xs)
    assert xt==chord*chord-q2-xs-xq and yt==-ys+chord*(xs-xt)
    # Translate back exactly; no exceptional target point or rank label.
    slopeback=(-yq-yt)/(xq-xt)
    assert xs==slopeback*slopeback-q2-xt-xq and ys==-yt+slopeback*(xt-xs)
    save('replay-map-%d.json'%index,dict(status='PASS_INDEPENDENT_DEGREE6_RATIONAL_TRANSLATION_MAP',
        index=index,sign=mp['sign'],original_degree=6,alternate_degree=2,
        generic_original_rank_lower_bound=18,
        proof='Verified rational source and inverse elliptic translation give a genuine genus0 original multisection; the existing all-prime theorem proves independence.',
        map_sha256=sha(OUT/('map-%d.json'%index)),replay_protocol_sha256=sha(OUT/'replay-protocol.json')))
    print('PASS_INDEPENDENT_DEGREE6_RATIONAL_TRANSLATION_MAP',index,flush=True)
def incidence():
    helpers=runpy.run_path(str(Path(__file__).with_name('verify_det1092_bifibration_conic_sources.sage')))
    rootcheck=helpers['check_root_certificate'];ip=read(OUT/'incidence-protocol.json');hashes(ip)
    old=read(ART/'det1092_rational_bisection_index_v1/controls.json')['rows'][0]['cases']
    assert ip['cases']==[dict(label=r['label'],parameter=r['parameter']) for r in old]
    assert len(ip['cases'])==9 and sum(QQ(r['parameter'])==0 for r in ip['cases'])==1
    rows=[]
    for index in [0,1]:
        mp=read(OUT/('map-%d.json'%index));T=dec(mp['T'])
        for case in ip['cases']:
            name='incidence-%d-%s.json'%(index,case['label']);record=read(OUT/name)
            assert record['map_sha256']==sha(OUT/('map-%d.json'%index))
            assert record['incidence_protocol_sha256']==sha(OUT/'incidence-protocol.json')
            f=T.numerator()-QQ(case['parameter'])*T.denominator()
            f*=lcm([v.denominator() for v in f]);f=R(f/gcd(list(f)))
            if f.leading_coefficient()<0:f=-f
            assert f.degree()==6 and f==R(record['polynomial']) and f[0]
            assert not record['infinity_preimage'] and record['rational_preimages']==[]
            assert record['status']=='EXACT_NO_RATIONAL_INCIDENCE'
            proof=rootcheck(f,record['root_proof'],ip['prime_pool'])
            rows.append(dict(index=index,label=case['label'],rational_preimages=0,proof=proof,
                certificate_sha256=sha(OUT/name)))
    save('replay-incidence.json',dict(status='PASS_INDEPENDENT_ALL18_EXACT_NONSPLITS',cases=rows,
        replay_protocol_sha256=sha(OUT/'replay-protocol.json'),
        boundary='These two fixed covers only. No exclusion of the fibres, other translates, or generic point existence.'))
    print('PASS_INDEPENDENT_ALL18_EXACT_NONSPLITS',flush=True)
def assemble():
    stages=['replay-geometry.json','replay-map-0.json','replay-map-1.json','replay-incidence.json']
    for name in stages:
        row=read(OUT/name);assert row['status'].startswith('PASS_')
        assert row['replay_protocol_sha256']==sha(OUT/'replay-protocol.json')
    save('independent-replay.json',dict(status='PASS_SHORT_TRANSLATION_DEGREE_LAW_TWO_COVERS_AND_PANEL',
        classification='new deduction and verified generic-only construction; no prospective302 seed',
        degree_laws=read(OUT/'replay-geometry.json')['degree_laws'],generic_cover_ranks=[18,18],
        rational_control_incidences=0,control_cases=18,new_specialized_rank_claims=0,
        inputs={str((OUT/name).relative_to(ROOT)):sha(OUT/name) for name in ['replay-protocol.json',*stages]},
        point_searches=0,source_enumerations=0,
        boundary='Minimum height2 among130 common sections, not a full MW12 shortest-vector theorem. No conclusion for n beyond the two incidence-tested signs.'))
    print('PASS_SHORT_TRANSLATION_DEGREE_LAW_TWO_COVERS_AND_PANEL',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['geometry','map','incidence','assemble']);parser.add_argument('--index',type=int,choices=[0,1]);args=parser.parse_args()
    signal.alarm(25);start=time.monotonic();preflight()
    if args.action=='geometry':geometry()
    elif args.action=='map':maps(args.index)
    elif args.action=='incidence':incidence()
    else:assemble()
    print('seconds',round(time.monotonic()-start,3),flush=True)

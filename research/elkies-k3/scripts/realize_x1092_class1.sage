#!/usr/bin/env sage-python
"""Bounded, checkpointed realization of X1092 J2 class 1; no fibre searches.

The requested embedding 33 is bound to a degree-two divisor on the known
rational surface by exact integral isometries. This need not extend the
particular common-core identification with embedding 25.
"""
import argparse
import csv
import hashlib
import json
import runpy
import time
import zipfile
from pathlib import Path
from sage.all import *
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
WORK = ROOT / 'artifacts/local/elkies-k3/x1092-class1-realization-v1'
PREFIX = 'x1092_class1_realization'
SOURCE = ART / 'curve302_recovered_mw17_parent_v1.json'
CENSUS = ART / 'det1092_pruned_rootless_j2_census_v1.json'
PRIORITY = ART / 'det1092_frame_realization_priority_v1.json'
PACKETS = ART / 'det1092_pruned_anchor_packets_v1.zip'
ORBITS = ART / 'curve302_parent_degree2_multisection_orbits_v1.tsv'
MASK = 109158
WORD = [-1,1,0,-1,-1,-1,0,1,1,0,0,0,0,0,-1,1,1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rows(m):
    return [list(map(int,r)) for r in m.rows()]

def poly(p):
    return list(map(str,p.list())) or ['0']

def rf(f):
    return {'numerator':poly(f.numerator()),'denominator':poly(f.denominator())}

def emit(stage, payload):
    payload.update(schema='x1092.class1.'+stage+'.v1', sage_version=SAGE_VERSION)
    p=ART/(PREFIX+'_'+stage+'_v1.json')
    tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(payload,indent=2)+'\n'); tmp.replace(p)
    print(stage, payload['status'], flush=True)

def negative_walls(D,G):
    # Same complete horizontal-wall test as certify_h92_marked_degree2_candidate.
    b=ZZ(D[1]); w=vector(ZZ,D[2:]); walls=[]; shells=[]
    for m in range(1,int(b)+1):
        cross=-b*m*G*w.column()
        Q=block_matrix(ZZ,[[b*b*G,cross],[cross.T,matrix(ZZ,[[m*m*(w*G*w)+1]])]])
        result=pari(Q).qfminim(2*b*b-1)
        shells.append({'degree':m,'full_short_vector_count':int(result[0])})
        for raw in matrix(ZZ,result[2]).T.rows():
            if abs(raw[-1])!=1: continue
            x=vector(ZZ,raw[:-1])*raw[-1]; n=x*G*x
            if (n-2)%(2*m): continue
            k=(n-2)//(2*m); pairing=D[0]*m+b*k-w*G*x
            if pairing<0: walls.append([int(k),m]+list(map(int,x)))
    return walls,shells

def marking():
    source=json.loads(SOURCE.read_text()); G=matrix(ZZ,source['generic_height_gram'])
    census=json.loads(CENSUS.read_text()); target=matrix(ZZ,census['rootless_classes'][0]['gram'])
    priority=json.loads(PRIORITY.read_text())['ranked_new_types'][0]
    assert priority['class_index']==1
    witness=priority['witness']
    assert [witness[k] for k in ('known_embedding_index','new_embedding_index','shared_sixth_index','common_core_determinant')]==[25,33,155,4100]
    with zipfile.ZipFile(PACKETS) as z: packet=json.loads(z.read('anchor-16.json'))
    left,right=[packet['embeddings'][i] for i in (25,33)]
    assert left['sixth_index']==right['sixth_index']==155
    C=matrix(ZZ,witness['common_core_basis_in_niemeier'])
    for e in (left,right):
        B=matrix(ZZ,e['complement_basis_in_ambient'])
        coords=matrix(ZZ,B.T.solve_right(C.T).T)
        assert coords*matrix(ZZ,e['gram'])*coords.T==matrix(ZZ,witness['common_core_gram'])
        assert coords.row_module().saturation()==coords.row_module()
    assert matrix(ZZ,witness['common_core_gram']).det()==4100
    orbit=next(r for r in csv.DictReader(ORBITS.open(),delimiter='\t') if int(r['orbit_mask'])==MASK)
    assert list(map(int,orbit['parent_MW17_w'].split()))==WORD
    w=vector(ZZ,WORD); assert w*G*w==12
    J=matrix(ZZ,[[0,1],[1,0]]); N=block_diagonal_matrix(J,-G)
    D=vector(ZZ,[3,2]+WORD); O=vector(ZZ,[-1,1]+[0]*17)
    U=matrix(ZZ,[D,D+O]); W=(U*N).right_kernel_matrix()
    H=-W*N*W.T; T=U.stack(W)
    assert abs(T.det())==1 and U*N*U.T==J
    assert T*N*T.T==block_diagonal_matrix(J,-H)
    iso=matrix(ZZ,pari(target).qfisom(pari(H)).sage())
    assert abs(iso.det())==1 and iso.T*H*iso==target
    requested=matrix(ZZ,right['gram'])
    reqiso=matrix(ZZ,pari(requested).qfisom(pari(H)).sage())
    assert abs(reqiso.det())==1 and reqiso.T*H*reqiso==requested
    requested_rows=reqiso.T*W
    assert requested_rows*N*requested_rows.T==-requested
    assert U*N*requested_rows.T==0
    assert pari(H).qfminim(2)[0]==0
    walls,shells=negative_walls(D,G); assert not walls
    payload={'status':'PASS_EXACT_RATIONAL_MARKED_NEF_U_AND_CLASS1_TRANSPORT',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,CENSUS,PRIORITY,PACKETS,ORBITS)},
      'requested_shared_core_witness':witness,
      'transport_scope':'Same requested J2 class; not asserted to extend the prescribed common-core identification.',
      'source_NS_gram':rows(N),'fibre_D':list(map(int,D)),'rational_zero_O':list(map(int,O)),
      'orbit_mask':MASK,'generic_trace_word':WORD,'transport_rows_D_D_plus_O_complement':rows(T),
      'transport_determinant':int(T.det()),'frame_gram':rows(H),
      'class1_representative_columns_in_complement':rows(iso),
      'embedding33_columns_in_complement':rows(reqiso),
      'embedding33_frame_rows_in_source_NS':rows(requested_rows),
      'full_gram_identity_verified':True,'complement_integral_isometry_verified':True,
      'nef_certificate':{'old_fibre_degree':2,'old_zero_intersection':1,
        'primitive':True,'effective':'RR and D.F_old=2>0 select D, since F_old is nef.',
        'old_vertical_root_rank':0,'complete_horizontal_wall_shells':shells,'negative_walls':walls,
        'completeness':'A negative-intersection irreducible curve is a (-2) component of effective D, hence has old fibre degree at most 2. Vertical roots do not exist.'},
      'equation':'UNKNOWN','rational_sections':'UNKNOWN','parameter_search':'BLOCKED'}
    emit('marking',payload)
    WORK.mkdir(parents=True,exist_ok=True)
    save((G,N,D,O,H,T),str(WORK/'marking.sobj'))

def trace():
    G,N,D,O,H,T=load(str(WORK/'marking.sobj'))
    loader=runpy.run_path(str(ROOT/'elliptic-curves/cas/load_curve302_recovered_parent.sage'))
    E,basis,unused=loader['load_curve302_recovered_parent'](SOURCE)
    K=E.base_ring(); R=K.ring(); t=R.gen()
    aa=E.a_invariants(); b2=aa[0]**2+4*aa[1]
    A=R(-E.c4()/48); B=R(-E.c6()/864); Es=EllipticCurve(K,[A,B])
    pts=[Es(p[0]+b2/12,p[1]+(aa[0]*p[0]+aa[2])/2) for p in basis]
    # Greedy exact norm descent keeps intermediate multiples small.
    todo=vector(ZZ,WORD); accum=vector(ZZ,[0]*17); P=Es(0); steps=[]
    while any(todo):
        choices=[]
        for i,c in enumerate(todo):
            if c:
                step=ZZ(1 if c>0 else -1); trial=vector(ZZ,list(accum)); trial[i]+=step
                choices.append((trial*G*trial,i,step,trial))
        _,i,step,trial=min(choices,key=lambda x:(x[0],x[1]))
        P+=step*pts[i]; accum=trial; todo[i]-=step
        steps.append([i,int(step),int(accum*G*accum)])
        print('trace partial',len(steps),'norm',steps[-1][-1],flush=True)
    assert list(accum)==WORD and accum*G*accum==12
    x,y=K(P[0]),K(P[1]); h=R(x.denominator()).sqrt()
    Nx=R(x*h*h); Ny=R(y*h**3)
    assert h.degree()==4 and Nx.degree()<=12 and Ny.degree()<=18
    assert Ny**2==Nx**3+A*Nx*h**4+B*h**6
    emit('trace',{'status':'PASS_EXACT_GENERIC_NORM12_TRACE','word':WORD,'steps':steps,
      'A':poly(A),'B':poly(B),'h':poly(h),'Nx':poly(Nx),'Ny':poly(Ny),
      'source_short_transform':{'x':'x_old+b2/12','y':'y_old+(a1*x_old+a3)/2'},
      'coefficient_max_bits':int(max(abs(c.numerator()).nbits()+c.denominator().nbits() for p in (h,Nx,Ny) for c in p))})
    save((Es,pts,A,B,h,Nx,Ny),str(WORK/'trace.sobj'))

def equation():
    Es,pts,Aold,Bold,h,Nx,Ny=load(str(WORK/'trace.sobj'))
    R=h.parent(); t=R.gen()
    columns=[(t**k*Nx)%(h*h) for k in range(8)]+[(-t**k*Ny)%(h*h) for k in range(2)]
    M=matrix(QQ,8,10,lambda i,j:columns[j][i]); ker=M.right_kernel_matrix()
    assert M.rank()==8 and ker.nrows()==2
    (a0,b0),(a1,b1)=[(R(list(row[:8])),R(list(row[8:]))) for row in ker]
    assert all((a*Nx-b*Ny)%(h*h)==0 for a,b in ((a0,b0),(a1,b1)))
    emit('rr',{'status':'PASS_EXACT_TWO_DIMENSIONAL_RR_KERNEL','a0':poly(a0),'b0':poly(b0),'a1':poly(a1),'b1':poly(b1),
      'kernel':[[str(c) for c in r] for r in ker],'constraint_rank':8})
    Ru=PolynomialRing(QQ,'u'); K=Ru.fraction_field(); u=Ru.gen()
    S=PolynomialRing(K,'t'); F=S.fraction_field(); tt=S.gen()
    lift=lambda p:S([K(c) for c in p])
    hh,nx,ny,AA=map(lift,(h,Nx,Ny,Aold))
    n=lift(a1)-u*lift(a0); d=u*lift(b0)-lift(b1)
    m=F(n/(d*hh)); xp=F(nx/hh**2); yp=F(ny/hh**3)
    print('RR complete; forming quartic',flush=True)
    rad=m**4-6*xp*m*m-8*yp*m-3*xp*xp-4*AA
    nr,dr=S(rad.numerator()),S(rad.denominator())
    sf=nr.gcd(nr.derivative()).monic(); q,rem=nr.quo_rem(sf**2)
    assert not rem and q.degree()==4 and q.gcd(q.derivative()).degree()==0
    ds=dr.sqrt(); assert rad==F(q*(sf/ds)**2)
    assert d.degree()==1
    t0=K(-d[0]/d[1]); norm=F(ds/d**2)
    assert norm.numerator().degree()==norm.denominator().degree()==0
    q0=K(norm.numerator()[0]/norm.denominator()[0]*n(t0)**2/(hh(t0)**2*sf(t0)))
    assert q0*q0==q(t0)
    # Pointed quartic to cubic, retaining genuinely birational formulas.
    shifted=q(tt+t0); e,dd,c,b,a=[K(shifted[i]) for i in range(5)]
    assert e==q0*q0
    ac2=c; ac4=b*dd-4*e*a; ac6=e*b*b+a*dd*dd-4*e*a*c
    E=EllipticCurve(K,[0,ac2,0,ac4,ac6])
    assert E.discriminant()!=0
    # Exact function-field checks, including inverse composition.
    L=F.extension(S(q).parent().fraction_field()['v'].gen()**2-F(q),'v')
    v=L.gen(); z=L(tt-t0)
    X=(2*q0*(v+q0)+dd*z)/z**2
    Y=((X*X-4*e*a)*z-dd*X-2*e*b)/(2*q0)
    assert Y*Y==X**3+c*X*X+ac4*X+ac6
    zi=(2*q0*Y+dd*X+2*e*b)/(X*X-4*e*a)
    vi=(X*zi*zi-dd*zi-2*e)/(2*q0)
    assert zi==z and vi==v
    oldx=(L(m*m-xp)+L(sf/ds)*v)/2
    oldy=L(m)*(oldx-L(xp))-L(yp)
    assert oldy*oldy==oldx**3+L(AA)*oldx+L(lift(Bold))
    L0=L(lift(a0))*(oldx*L(hh)**2-L(nx))+L(lift(b0))*(oldy*L(hh)**3+L(ny))
    L1=L(lift(a1))*(oldx*L(hh)**2-L(nx))+L(lift(b1))*(oldy*L(hh)**3+L(ny))
    assert L1/L0==u
    emit('equation',{'status':'PASS_EXACT_RATIONAL_EQUATION_AND_BIRATIONAL_MAPS',
      'a_invariants':[rf(K(x)) for x in E.a_invariants()],
      'quartic_coefficients':[rf(K(q[i])) for i in range(5)],'quartic_zero':[rf(t0),rf(q0)],
      'shifted_quartic_coefficients':[rf(K(shifted[i])) for i in range(5)],
      'radical_square_factor':[rf(K(x)) for x in sf], 'radical_denominator_sqrt':[rf(K(x)) for x in ds],
      'maps':{'source_to_pencil':'u=(a1*(x*h^2-Nx)+b1*(y*h^3+Ny))/(a0*(x*h^2-Nx)+b0*(y*h^3+Ny))',
       'source_to_quartic':'v=(2*x-(m^2-xP))*ds/sf; z=t-t0',
       'quartic_to_cubic':'X=(2*q0*(v+q0)+dd*z)/z^2; Y=((X^2-4*e*a)*z-dd*X-2*e*b)/(2*q0)',
       'cubic_to_quartic':'z=(2*q0*Y+dd*X+2*e*b)/(X^2-4*e*a); t=z+t0; v=(X*z^2-dd*z-2*e)/(2*q0)',
       'quartic_to_source':'x=(m^2-xP+v*sf/ds)/2; y=m*(x-xP)-yP',
       'verified_in_exact_function_field':True},
      'generic_sections':'UNKNOWN','carrier_strict_gate':'NOT_RUN','parameter_search':'BLOCKED'})
    save((E,q,t0,q0,sf,ds,m,xp,yp,a0,b0,a1,b1),str(WORK/'equation.sobj'))

def section_plan():
    G,N,D,O,H,T=load(str(WORK/'marking.sobj')); inv=T.inverse()
    short=matrix(ZZ,pari(G).qfminim(8)[2]); candidates=[]
    for col in short.columns()+[vector(ZZ,WORD)]:
        for v in (col,-col):
            height=v*G*v; div=vector(ZZ,[(height-2)//2,1]+list(v))
            if div*N*D!=1: continue
            child=div*inv
            assert child[1]==1 and all(x in ZZ for x in child)
            candidates.append((v,vector(ZZ,child[2:])))
    candidates.sort(key=lambda x:(x[0]*G*x[0],sum(abs(c) for c in x[0]),tuple(x[0])))
    helper=runpy.run_path(str(ROOT/'elkies-k3/scripts/plan_r17_norm12_direct_section_basis.sage'))
    result=helper['reduce_basis']([],candidates)
    glue_word=None
    if result is None:
        glue=[]
        for row in csv.DictReader(ORBITS.open(),delimiter='\t'):
            if int(row['minimum_norm'])!=10: continue
            v=vector(ZZ,list(map(int,row['parent_MW17_w'].split())))
            for sign in (1,-1):
                vv=sign*v; div=vector(ZZ,[2,2]+list(vv))
                if div*N*D==1:
                    child=div*inv
                    glue.append((sum(abs(c) for c in vv),int(row['orbit_mask']),sign,vv,vector(ZZ,child[2:])))
        glue.sort(key=lambda x:tuple(x[:3]))
        for _,mask,sign,v,child in glue:
            result=helper['reduce_basis']([child],candidates)
            if result is not None:
                glue_word=list(map(int,v));break
        if result is None:
            emit('section_plan',{'status':'UNKNOWN_BOUNDED_OLD_SECTION_AND_BISECTION_SPAN','candidate_count':len(candidates),
              'bisection_candidate_count':len(glue),'rank':int(matrix(ZZ,[v for _,v in candidates]).rank())}); return
    selected,B=result
    assert B.rank()==17
    emit('section_plan',{'status':'PASS_EXACT_RANK17_SECTION_CLASS_PLAN','height_bound':12,
      'candidate_count':len(candidates),'selected_source_words':rows(matrix(ZZ,selected)),
      'glue_bisection_word':glue_word,
      'child_frame_coordinates':rows(B),'subgroup_index':int(abs(B.det())),
      'height_gram':rows(B*H*B.T),'height_determinant':int((B*H*B.T).det()),
      'rational_coordinate_recovery':'PENDING'})

def sections():
    E,q,t0,q0,sf,ds,m,xp,yp,a0,b0,a1,b1=load(str(WORK/'equation.sobj'))
    Es,pts,Aold,Bold,h,Nx,Ny=load(str(WORK/'trace.sobj'))
    G,N,D,O,H,T=load(str(WORK/'marking.sobj'))
    plan=json.loads((ART/(PREFIX+'_section_plan_v1.json')).read_text())
    assert plan['status']=='PASS_EXACT_RANK17_SECTION_CLASS_PLAN'
    K=E.base_ring(); u=K.gen(); tt=q.parent().gen()
    shifted=q(tt+t0); e,dd,c,b,a=[K(shifted[i]) for i in range(5)]
    helpers=runpy.run_path(str(ROOT/'elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage'))
    ev=helpers['evaluate_rational']; invert=helpers['invert_mobius']
    records=[]
    sources=[]
    glue=plan.get('glue_bisection_word')
    if glue:
        # The residual line through P_{-v} has class (2,2,v).
        P=sum((-ZZ(c)*p for c,p in zip(glue,pts)),Es(0)); R=h.parent()
        hp=R(P[0].denominator()).sqrt(); nx=R(P[0]*hp**2); ny=R(P[1]*hp**3)
        assert hp.degree()==3 and vector(ZZ,glue)*G*vector(ZZ,glue)==10
        f1=(-ny*nx.inverse_mod(hp*hp))%(hp*hp)
        f0,rem=(-(f1*nx+ny)).quo_rem(hp*hp); assert not rem
        branch,rem=(f1**4-6*nx*f1*f1-8*ny*f1-3*nx*nx-4*Aold*hp**4).quo_rem(hp**6)
        assert not rem and branch.degree()<=2
        bb,rem=(f1*f1-nx).quo_rem(hp*hp); assert not rem
        x0,x1=bb/2,hp/2
        y0,rem=(-(2*f0+f1*bb)).quo_rem(2*hp); assert not rem
        y1=-f1/2
        assert y0*y0+y1*y1*branch==x0**3+3*x0*x1*x1*branch+Aold*x0+Bold
        assert 2*y0*y1==3*x0*x0*x1+x1**3*branch+Aold*x1
        S=q.parent(); F=S.fraction_field(); lift=lambda p:S([K(c) for c in p])
        lc0=a0*(x0*h*h-Nx)+b0*(y0*h**3+Ny)
        ls0=a0*x1*h*h+b0*y1*h**3
        lc1=a1*(x0*h*h-Nx)+b1*(y0*h**3+Ny)
        ls1=a1*x1*h*h+b1*y1*h**3
        ss=F(-(lift(lc1)-u*lift(lc0))/(lift(ls1)-u*lift(ls0)))
        incidence=S((ss*ss-lift(branch)).numerator())
        print('glue incidence degree',incidence.degree(),flush=True)
        candidates=[]
        for factor,exponent in incidence.factor():
            if factor.degree()!=1: continue
            tv=K(-factor[0]/factor[1])
            if not ss.denominator()(tv): continue
            sv=K(ev(ss,tv)); assert sv*sv==lift(branch)(tv)
            xx=K(lift(x0)(tv)+lift(x1)(tv)*sv); yy=K(lift(y0)(tv)+lift(y1)(tv)*sv)
            if not (lift(lc0)(tv)+lift(ls0)(tv)*sv): continue
            assert (lift(lc1)(tv)+lift(ls1)(tv)*sv)/(lift(lc0)(tv)+lift(ls0)(tv)*sv)==u
            candidates.append((tv,xx,yy))
        assert len(candidates)==1
        sources.append(('rational_bisection',glue,*candidates[0]))
        emit('glue',{'status':'PASS_EXACT_RATIONAL_BISECTION_SECTION','source_word':glue,
          'branch':poly(branch),'x0':poly(x0),'x1':poly(x1),'y0':poly(y0),'y1':poly(y1),
          'new_base_degree':1,'source_t':rf(candidates[0][0])})
    for word in plan['selected_source_words']:
        P=sum((ZZ(c)*p for c,p in zip(word,pts)),Es(0)); x,y=P[0],P[1]
        l0=a0*(x*h*h-Nx)+b0*(y*h**3+Ny)
        l1=a1*(x*h*h-Nx)+b1*(y*h**3+Ny)
        base=l1/l0; tnew=K(invert(base,u))
        sources.append(('old_section',word,tnew,K(ev(x,tnew)),K(ev(y,tnew))))
    for i,(kind,word,tnew,xx,yy) in enumerate(sources):
        print('transport source',i,kind,flush=True)
        assert yy*yy==xx**3+ev(Es.a4(),tnew)*xx+ev(Es.a6(),tnew)
        mm=K(ev(m,tnew)); px=K(ev(xp,tnew))
        assert yy+ev(yp,tnew)==mm*(xx-px)
        vv=(2*xx-(mm*mm-px))*K(ds(tnew))/K(sf(tnew))
        assert vv*vv==q(tnew)
        z=tnew-t0
        X=(2*q0*(vv+q0)+dd*z)/z**2
        Y=((X*X-4*e*a)*z-dd*X-2*e*b)/(2*q0)
        point=E(X,Y)
        assert not point.is_zero()
        records.append({'source_kind':kind,'source_word':word,'child_frame_coordinates':plan['child_frame_coordinates'][i],
          't':rf(tnew),'X':rf(K(X)),'Y':rf(K(Y)),'equation_verified':True})
        print('section',i+1,'verified',flush=True)
        emit('sections_progress',{'status':'PARTIAL_EXACT_SECTIONS','sections':records})
    B=matrix(ZZ,plan['child_frame_coordinates']); gram=B*H*B.T
    assert len(records)==17 and gram.is_positive_definite() and gram.det()==1092*B.det()**2
    emit('sections',{'status':'PASS_EXACT_RATIONAL_GENERIC_RANK17',
      'sections':records,'height_gram':rows(gram),'height_determinant':int(gram.det()),
      'index_in_geometric_frame':int(abs(B.det())),
      'height_proof':'Exact divisor transport; all fibres irreducible by rootless complement, so Shioda projections have height -NS pairing.',
      'rank_lower':17,'geometric_Picard_rank':19,'root_rank':0,'rank_upper':17,
      'carrier_strict_gate':'NOT_RUN','parameter_search':'BLOCKED'})

def normalize():
    E,*_=load(str(WORK/'equation.sobj')); K=E.base_ring(); R=K.ring()
    A=-E.c4()/48; B=-E.c6()/864
    # Only factors that can affect the 4/6 scaling, never factor Delta.
    factors=set(f for p in (A.numerator().gcd(B.numerator()),A.denominator(),B.denominator()) for f,e in p.factor())
    gauge=K(1)
    for f in factors:
        power=min(A.valuation(f)//4,B.valuation(f)//6); gauge*=f**(-power)
    aa,bb=R(A*gauge**4),R(B*gauge**6)
    assert aa.degree()<=8 and bb.degree()<=12
    section_data=json.loads((ART/(PREFIX+'_sections_v1.json')).read_text())
    decode=lambda r:K(R(r['numerator']))/R(r['denominator'])
    EE=EllipticCurve(K,[aa,bb]); points=[]
    for r in section_data['sections']:
        x,y=decode(r['X']),decode(r['Y'])
        P=EE(gauge**2*(x+E.b2()/12),gauge**3*y); points.append(P)
    emit('parent',{'status':'PASS_EXACT_POLYNOMIAL_MW17_PARENT','family':'x1092-j2-class1',
      'a_invariants':[rf(K(x)) for x in EE.a_invariants()],
      'basis_weierstrass_coordinates':[[rf(K(p[0])),rf(K(p[1]))] for p in points],
      'generic_height_gram':section_data['height_gram'],'generic_rank':17,
      'gauge_from_pointed_cubic':rf(gauge),'short_translation_b2_over12':rf(K(E.b2()/12)),
      'coefficient_bits':int(max(abs(c.numerator()).nbits()+c.denominator().nbits() for f in (aa,bb) for c in f)),
      'carrier_strict_gate':'NOT_RUN','parameter_search':'BLOCKED'})
    save((EE,points,matrix(ZZ,section_data['height_gram'])),str(WORK/'parent.sobj'))

def compact():
    E,points,G=load(str(WORK/'parent.sobj'))
    eq=load(str(WORK/'equation.sobj')); b0,b1=eq[-3],eq[-1]
    R=PolynomialRing(QQ,'s'); K=R.fraction_field(); s=R.gen()
    u_of_s=K(R(b1.list())/R(b0.list())); den=R(u_of_s.denominator())
    helper=runpy.run_path(str(ROOT/'elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage'))
    ev=helper['evaluate_rational']
    aa=R(ev(E.a4(),u_of_s)*den**8); bb=R(ev(E.a6(),u_of_s)*den**12)
    def content(p):
        return QQ(gcd([c.numerator() for c in p]))/lcm([c.denominator() for c in p])
    ca,cb=abs(content(aa)),abs(content(bb))
    n=gcd(ca.numerator()**3,cb.numerator()**2)
    d=lcm(ca.denominator()**3,cb.denominator()**2)
    def powerpart(value):
        scale=ZZ(1)
        for p in prime_range(1000):
            exponent=value.valuation(p)
            scale*=p**(exponent//12); value//=p**exponent
        root=value.nth_root(12,truncate_mode=True)[0]
        if root**12==value: scale*=root
        return scale
    scale=QQ(powerpart(n))/powerpart(d)
    aa=R(aa/scale**4);bb=R(bb/scale**6); EE=EllipticCurve(K,[aa,bb])
    transformed=[EE(ev(p[0],u_of_s)*den**4/scale**2,ev(p[1],u_of_s)*den**6/scale**3) for p in points]
    emit('compact_parent',{'status':'PASS_EXACT_POLYNOMIAL_MW17_PARENT','family':'x1092-j2-class1',
      'a_invariants':[rf(K(x)) for x in EE.a_invariants()],
      'basis_weierstrass_coordinates':[[rf(K(p[0])),rf(K(p[1]))] for p in transformed],
      'generic_height_gram':rows(G),'generic_rank':17,'u_of_s':rf(u_of_s),
      'coordinate_scale':'x_compact=x_parent(u(s))*den(u(s))^4/scale^2; y analog power6/scale^3',
      'scale':str(scale),'coefficient_bits':int(max(abs(c.numerator()).nbits()+c.denominator().nbits() for f in (aa,bb) for c in f)),
      'carrier_strict_gate':'NOT_RUN','parameter_search':'BLOCKED'})
    save((EE,transformed,G),str(WORK/'compact_parent.sobj'))

def carrier_plan():
    from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice
    from itertools import combinations,product
    E,pts,G=load(str(WORK/'compact_parent.sobj')); candidates=set()
    for indices in combinations(range(17),3):
        for signs in product((1,-1),repeat=2):
            v=vector(ZZ,[0]*17)
            for i,c in zip(indices,(1,)+signs):v[i]=c
            if v*G*v==10:candidates.add(tuple(map(int,v)))
    ordered=sorted(candidates,key=lambda w:hashlib.sha256(('x1092-class1-carrier-v1/'+str(w)).encode()).hexdigest())
    selected=[]; checked=[]; L=IntegralLattice(G)
    for word in ordered:
        v=vector(ZZ,word); q=vector(ZZ,next(L.enumerate_close_vectors(vector(QQ,v)/2)))
        minimum=(v-2*q)*G*(v-2*q)
        checked.append({'word':list(word),'coset_minimum':int(minimum)})
        if minimum==10:selected.append(list(word))
        if len(selected)==4:break
    emit('carrier_protocol_triples',{'status':'FROZEN_GENERIC_ONLY_CARRIER_PANEL',
      'parent_sha256':sha(ART/(PREFIX+'_compact_parent_v1.json')),
      'rule':'First four hash-ordered norm10 signed sums of exactly three distinct basis sections with exact coset minimum10. The earlier two-section candidate window had no survivors.',
      'words':selected,'candidate_count':len(candidates),'checked':checked,
      'limits':{'cold_trace_seconds_per_word':180,'words':4,'parameters':0,'point_searches':0},
      'strict_dimension':'UNKNOWN','boundary':'A rational NS bisection class is not a constructed arithmetic strict class.'})

def carrier():
    E,pts,G=load(str(WORK/'compact_parent.sobj')); K=E.base_ring(); R=K.ring()
    protocol=json.loads((ART/(PREFIX+'_carrier_protocol_triples_v1.json')).read_text()); word=protocol['words'][args.index]
    assert protocol['parent_sha256']==sha(ART/(PREFIX+'_compact_parent_v1.json'))
    P=sum((-ZZ(c)*p for c,p in zip(word,pts)),E(0))
    Q=sum((-ZZ(c)*p for c,p in reversed(list(zip(word,pts)))),E(0));assert P==Q
    h=R(P[0].denominator()).sqrt().monic(); nx=R(P[0]*h*h);ny=R(P[1]*h**3)
    assert h.degree()==3
    import sys
    sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
    from euclidean_seed_sieve import conic_from_trace
    conic=conic_from_trace(poly(R(E.a4())),poly(R(E.a6())),poly(h),poly(nx),poly(ny))
    emit('carrier_'+str(args.index),{'status':'PASS_EXACT_GENERIC_RATIONAL_BISECTION_CARRIER',
      'word':word,'parent_sha256':protocol['parent_sha256'],'conic':conic,
      'generic_trace_constructed_cold':True,'reversed_order_trace_verified':True,
      'cold_total_seconds':time.monotonic()-start,
      'constructed_strict_class':'UNKNOWN','cover_solubility_over_Q':'UNKNOWN',
      'new_rational_point':'UNKNOWN','independent_quotient_direction':'UNKNOWN',
      'parameters_tested':0})

def strict_preflight():
    dispatcher=runpy.run_path(str(ROOT/'elliptic-curves/rank-jump/prospective_constructed_strict_seed_gate.py'))
    carriers=[ART/(PREFIX+'_carrier_'+str(i)+'_v1.json') for i in range(4)]
    parent=ART/(PREFIX+'_compact_parent_v1.json')
    for p in carriers:
        r=json.loads(p.read_text())
        assert r['status']=='PASS_EXACT_GENERIC_RATIONAL_BISECTION_CARRIER' and r['parent_sha256']==sha(parent)
    signal={'schema':'rank-jump.prospective-constructed-strict-signals.v1',
      'provenance':'PROSPECTIVE','historical_point_inputs':0,'candidate_id':'x1092-j2-class1',
      'generic_strict_dimension':'UNKNOWN','principal_dependency':{'status':'UNKNOWN'},
      'constructed_strict_class':{'status':'UNKNOWN'},'cover_solubility':{'status':'UNKNOWN'},
      'new_rational_point':{'status':'UNKNOWN'},'independent_quotient_direction':{'status':'UNKNOWN'}}
    decision=dispatcher['evaluate'](signal)
    assert decision['next_action']=='WAIT_FOR_CONSTRUCTED_STRICT_CLASS'
    data=json.loads(parent.read_text())
    emit('strict_preflight',{'status':'UNKNOWN_NO_APPLICABLE_PROSPECTIVE_STRICT_CONSTRUCTOR',
      'parent_sha256':sha(parent),'generic_cubic_coefficients_in_x':[data['a_invariants'][4],data['a_invariants'][3],{'numerator':['0'],'denominator':['1']},{'numerator':['1'],'denominator':['1']}],
      'input_signal':signal,'dispatcher_receipt':decision,
      'complete_inherited_everywhere_even_half_ideal_image':'UNKNOWN',
      'strict_character_subspace':'UNKNOWN','ordinary_unramified_character_subspace':'UNKNOWN',
      'required_next_certificate':'An explicit principal-square identity and annihilating character proving a nonzero strict class modulo the full inherited image, bound to this cubic algebra.',
      'reason':'Available successful strict-class solvers have pinned arithmetic-field inputs. Neither the curve302 strict dimension nor its class certificates transfer to a different marked fibration. Exact norm10 NS carriers do not construct arithmetic strict classes.',
      'bindings':{str(p.relative_to(ROOT)):sha(p) for p in [parent,*carriers,ROOT/'elliptic-curves/rank-jump/prospective_constructed_strict_seed_gate.py']},
      'parameter_panel':'BLOCKED_NOT_COMMISSIONED','parameters_tested':0})

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('stage',choices=['marking','trace','equation','section_plan','sections','normalize','compact','carrier_plan','carrier','strict_preflight'])
    ap.add_argument('--index',type=int,default=0)
    args=ap.parse_args(); start=time.monotonic()
    globals()[args.stage]()
    print('elapsed_seconds',time.monotonic()-start,flush=True)

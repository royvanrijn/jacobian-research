#!/usr/bin/env sage-python
"""Independent genus-one carrier, rational maps, heights and finite groups.

No constructor import, CVP repetition, exceptional-point execution oracle,
point search, class group or later artifact.25-second total cap.
"""
import csv,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_norm8_seed_cover_v2'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance(d):
    for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest,(path,'HASH_MISMATCH')
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def residue(v,p):
    v=QQ(v);assert v.denominator()%p
    return int(v.numerator()*v.denominator().inverse_mod(p)%p)
def finite_group(p,A,B):
    points=[None]+[(x,y) for x in range(p) for y in range(p) if (y*y-x**3-A*x-B)%p==0]
    def add(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;xx,yy=Q
        if x==xx and (y+yy)%p==0:return None
        if x==xx:m=(3*x*x+A)*pow(2*y,-1,p)%p
        else:m=(yy-y)*pow(xx-x,-1,p)%p
        xxx=(m*m-x-xx)%p;return (xxx,(m*(x-xxx)-y)%p)
    def times(n,P):
        Q=None
        while n:
            if n&1:Q=add(Q,P)
            P=add(P,P);n//=2
        return Q
    return points,times
def verify():
    paths=[DIR/f for f in ['protocol.json','generic.json','calibration-protocol.json','calibration.json','covers.json']]
    protocol,generic,cp,cal,panel=map(read,paths)
    for d in [protocol,generic,cp,cal,panel]:provenance(d)
    assert cp['limits']['fixed_members']==2 and cp['limits']['historical_first_point_inputs']==1
    assert cp['limits']['later_point_inputs']==cp['limits']['point_searches']==cp['limits']['V3_inputs']==0
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    roster=read(ART/'det1092_rr_generic_point_controls_v2/protocol.json')
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    G=matrix(QQ,parent['generic_height_gram']);word=vector(ZZ,[0]*14+[1,-1,0])
    assert list(word)==generic['selection']['word'] and word*G*word==8 and G.det()==1092
    # Recheck the source-only selection, without re-enumerating any orbit.
    with (ART/'curve302_parent_degree2_multisection_orbits_v1.tsv').open() as stream:
        def key(row):
            w=list(map(int,row['parent_MW17_w'].split()))
            return sum(map(abs,w)),max(map(abs,w)),int(row['orbit_mask'])
        chosen=min((r for r in csv.DictReader(stream,delimiter='\t') if r['minimum_norm']=='8'),key=key)
    assert int(chosen['orbit_mask'])==generic['selection']['orbit']==20124
    assert list(map(int,chosen['parent_MW17_w'].split()))==list(word)
    C=basis[14]-basis[15]
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(generic['pole_h']);nx=R(generic['nx']);ny=R(generic['ny']);shift=R(generic['shift'])
    assert h*h==C[0].denominator() and h.degree()==2
    assert nx==cx*h*h and ny==cy*h**3 and (shift*nx-ny)%(h*h)==0
    assert shift.degree()<4
    Z=PolynomialRing(QQ,'z');z=Z.gen();F=Z.fraction_field()
    # The two section restrictions are checked directly on the old sections.
    for row,P in zip(generic['sections'],[basis[14],-basis[15]]):
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        slope=(yy+cy)/(xx-cx);param=(h*slope+shift)/(h*h)
        T=F(Z(row['t_of_z']['numerator']))/Z(row['t_of_z']['denominator'])
        W=F(Z(row['W_of_z']['numerator']))/Z(row['W_of_z']['denominator'])
        def at(f):return F(f.numerator()(T))/f.denominator()(T)
        assert at(param)==z
        assert W==(2*at(xx)+at(cx)-at(slope)**2)/h(T)
    coefficients=list(map(Z,generic['quartic_t_coefficients_in_z']))
    assert len(panel['covers'])==2
    checks=[];reconstructed_first=None
    for index,row in enumerate(panel['covers']):
        assert row==read(DIR/f'cover-{index:02d}.json') and row['index']==index
        zvalue=QQ(row['z']);assert zvalue==([QQ(0),QQ(cal['z_star'])][index])
        q=R(row['quartic_coefficients'])
        assert q==R([v(zvalue) for v in coefficients])
        assert q.degree()==4 and q.gcd(q.derivative()).degree()==0
        assert q.gcd(R(E.discriminant())).degree()==0 and E.j_invariant().derivative()
        M=h*h*zvalue-shift;m=K(M)/h
        raw=M**4-6*nx*M*M-8*ny*M-3*nx*nx-4*R(-E.c4()/48)*h**4
        assert raw==h**6*q
        x0,x1,y0,y1=[dec(row['maps'][k]) for k in ['x0','x1','y0','y1']]
        assert x0==(-cx+m*m)/2-E.b2()/12 and x1==h/2
        assert y0==m*(x0+E.b2()/12-cx)-cy-(E.a1()*x0+E.a3())/2
        assert y1==(m-E.a1()/2)*x1
        a1,a2,a3,a4,a6=ai
        assert y0*y0+y1*y1*q+a1*(x0*y0+x1*y1*q)+a3*y0==x0**3+3*x0*x1*x1*q+a2*(x0*x0+x1*x1*q)+a4*x0+a6
        assert 2*y0*y1+a1*(x0*y1+x1*y0)+a3*y1==3*x0*x0*x1+x1**3*q+2*a2*x0*x1+a4*x1
        line=list(map(R,row['RR_line']));f0,f1,f2=line
        assert gcd(line).degree()==0 and f2==h
        assert f0+f1*x0+f2*y0==0 and f1*x1+f2*y1==0
        assert all(f.degree()<=d for f,d in zip(line,[8,4,2])) and f2.degree()==2
        minus=-C;assert f0+f1*minus[0]+f2*minus[1]==0
        # Residual class D=(3O+8F)-(-C)=2O+4F+phi(word).
        # Smooth quartic normalization has genus1=p_a(D), hence D.O=0.
        cross=G*word
        H=(2*G).augment(matrix(QQ,17,1,list(cross))).stack(matrix(QQ,1,18,list(cross)+[8]))
        assert H==matrix(QQ,row['pullback_height_gram']) and H.is_positive_definite()
        assert H.det()==4*2**17*1092 and row['height_Schur_complement']=='4'
        assert row['anti_invariant_height']=='16'
        t0,s=map(QQ,row['inherited_basepoint']);assert s and s*s==q(t0)
        assert x0(t0)+x1(t0)*s==basis[14][0](t0)
        assert y0(t0)+y1(t0)*s==basis[14][1](t0)

        # Derive a birational Weierstrass model anchored at this inherited
        # point, independently of the constructor's quartic covariants.
        shifted=q(t+t0);q0,q1,q2,q3,q4=shifted.list();assert q0==s*s
        A4=q1*q3-4*q0*q4;A6=q0*q3*q3+q1*q1*q4-4*q0*q2*q4
        J0=EllipticCurve(QQ,[0,q2,0,A4,A6])
        I=12*q0*q4-3*q1*q3+q2*q2
        J=72*q0*q2*q4+9*q1*q2*q3-27*q4*q1*q1-27*q0*q3*q3-2*q2**3
        assert J0.c4()==16*I and J0.c6()==32*J
        proof=row['positive_rank_base'];ainv=list(map(QQ,proof['Jacobian_a_invariants']))
        assert ainv==[0,0,0,-27*I,-27*J]
        # The conjugate inherited point maps to this Jacobian point.
        Xbar=q1*q1/(4*s*s)-q2
        Ybar=-(q1*Xbar+2*s*s*q3)/(2*s)
        assert J0([Xbar,Ybar])
        Tpoint=[9*Xbar+3*q2,27*Ybar]
        candidate=list(map(QQ,proof['rational_nontorsion_point']))
        assert candidate==[Tpoint[0],-Tpoint[1]]
        # Exact inverse over the quadratic function field of J0.
        U=PolynomialRing(QQ,'X');XX=U.gen();UF=U.fraction_field()
        V=PolynomialRing(UF,'Y');YY=V.gen();cub=XX**3+q2*XX**2+A4*XX+A6
        algebra=V.quotient(YY*YY-cub,'y');yy=algebra.gen()
        du=(2*s*yy+q1*XX+2*s*s*q3)/(XX*XX-4*s*s*q4)
        ww=XX*du*du/(2*s)-s-q1*du/(2*s)
        assert ww*ww==sum(shifted[i]*du**i for i in range(5))
        assert 2*s*(ww+s)+q1*du==XX*du*du
        assert ((XX*XX-4*s*s*q4)*du-q1*XX-2*s*s*q3)/(2*s)==yy
        # Nonsquare leading coefficient removes rational denominator-zero
        # points in the inverse; the identity maps separately to(t0,s).
        assert q4>0
        sn,sd=q4.numerator().isqrt(),q4.denominator().isqrt()
        assert sn*sn!=q4.numerator() or sd*sd!=q4.denominator()
        counts=[];bound=ZZ(0);good={}
        assert [r['prime'] for r in proof['finite_exposures']]==cp['proof_primes']
        for exposure in proof['finite_exposures']:
            p=exposure['prime']
            if any(v.denominator()%p==0 for v in [*ainv,*candidate]):
                assert exposure['status']=='SKIP_DENOMINATOR';continue
            AA,BB=residue(ainv[3],p),residue(ainv[4],p)
            if (4*AA**3+27*BB**2)%p==0:
                assert exposure['status']=='SKIP_BAD_REDUCTION';continue
            pts,times=finite_group(p,AA,BB)
            assert exposure['status']=='GOOD' and len(pts)==exposure['order']
            pt=tuple(residue(v,p) for v in candidate);assert pt in pts
            bound=gcd(bound,len(pts));good[p]=(pts,times,pt)
            counts.append({'prime':p,'points':pts,'order':len(pts)})
        witness=proof['nontorsion_witness'];p=witness['prime']
        assert bound==witness['bound'] and p in good
        pts,times,pt=good[p];mul=times(int(bound),pt)
        assert mul is not None and list(mul)==witness['multiple'][:2]
        outcomes=[]
        assert len(row['control_outcomes'])==len(roster['cases'])==9
        for case,out in zip(roster['cases'],row['control_outcomes']):
            assert out['label']==case['label'] and out['parameter']==case['parameter']
            tau=QQ(case['parameter']);val=q(tau);assert str(val)==out['value']
            if val<0:assert not out['split']
            else:
                sn0,sd0=map(int,out['square_root_bounds']);n0,d0=val.numerator(),val.denominator()
                assert sn0*sn0<=n0<(sn0+1)**2 and sd0*sd0<=d0<(sd0+1)**2
                assert out['split']==(sn0*sn0==n0 and sd0*sd0==d0)
            if out['split']:
                root=QQ(sn0)/sd0;Ps=[x0(tau)+x1(tau)*root,y0(tau)+y1(tau)*root]
                assert list(map(str,Ps))==out['literal_point']
                if index==1 and tau==0:reconstructed_first=Ps
            outcomes.append({'label':case['label'],'split':out['split']})
        checks.append({'index':index,'genus':1,'pullback_rank_lower_bound':18,
          'height_gram_determinant':str(H.det()),'anti_invariant_height':16,
          'Jacobian_rank_lower_bound':1,'torsion_order_bound':int(bound),
          'Jacobian_generator_from_conjugate_inherited_point':list(map(str,Tpoint)),
          'explicit_base_generator':{'t0':str(t0),'s':str(s),'shifted_quartic':list(map(str,shifted.list())),
             'Jacobian_cubic_coefficients':list(map(str,[q2,A4,A6])),
             'infinite_rational_points':'Take nonzero multiples of (Xbar,Ybar) on Y^2=X^3+q2*X^2+A4*X+A6, then u=(2sY+q1X+2s^2q3)/(X^2-4s^2q4), t=t0+u, W=Xu^2/(2s)-s-q1u/(2s).',
             'generator_XY':list(map(str,[Xbar,Ybar])),
             'leading_nonsquare_bounds':list(map(str,[sn,sd]))},
          'finite_groups':counts,'control_outcomes':outcomes})
    # Only after witness-free reconstruction inspect the old first point
    # for diagnostic equality. Its completed rank18 proof is a dependency.
    first=read(ART/'det1092_first_centre_rr_net_replay_v1.json')
    E0=EllipticCurve(QQ,[v(0) for v in ai]);P=E0(first['reconstructed_point_literal302'])
    Q=E0(reconstructed_first);C0=E0([C[0](0),C[1](0)])
    assert Q==P or C0-Q==P
    assert first['independence']['rank']==18
    assert sum(r['split'] for r in checks[0]['control_outcomes'])==0
    assert [r['label'] for r in checks[1]['control_outcomes'] if r['split']]==['302-generic-section-control']
    result={'status':'PASS_INDEPENDENT_POSITIVE_RANK_GENUS_ONE_FIRST_SEED_COVER',
      'classification':'new constructive deduction and verified application',
      'covers':checks,'reconstructed_first_literal_point':list(map(str,reconstructed_first)),
      'first_or_centre_complement':True,'limits':cp['limits'],
      'boundary':'One generic pencil, a generic member and a retrospectively calibrated member through the first unlock. Both have explicit infinite rational base-point sources and elliptic function-field rank>=18. Only the calibrated member hits302; all eight controls are nonsplit. Calibration has not been removed from selection. No new302 seed, later point, V3 input, parameter sweep or rank upper bound.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'replay.json',result)
    print(result['status'],'both bases positive rank; first302 recovered from equations',flush=True)
if __name__=='__main__':
    signal.alarm(25);verify()

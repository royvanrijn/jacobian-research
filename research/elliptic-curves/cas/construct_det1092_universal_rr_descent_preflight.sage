#!/usr/bin/env sage-python
"""Universal RR sextic plus a separately labelled retrospective class audit.

The universal construction reads only generic RR data. Sealed first-point
data enter the diagnostic after that construction. No descent or point search.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,VectorSpace,matrix,vector,lcm,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_universal_rr_descent_preflight_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';loader=CAS/'load_curve302_recovered_parent.sage'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text());E,base,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    R=E.base_ring().ring();C=-sum((n*P for n,P in zip(net['trace_word'],base)),E(0))
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2;aa=-E.c4()/48;bb=-E.c6()/864
    h=R(C[0].denominator().sqrt());nx=R(cx*h*h);ny=R(cy*h**3)
    UV=PolynomialRing(QQ,['u','v']);u,v=UV.gens();S=PolynomialRing(UV,'T');T=S.gen()
    A=[S(R(row)) for row in net['A']];B=[S(R(row)) for row in net['B']]
    f1=B[1]+(u+v*T)*A[1];f2=B[2]+(u+v*T)*A[2];mnum=-f1+S(R(E.a1()))*f2/2
    L,rem=f2.quo_rem(S(h));assert not rem
    raw=mnum**4-6*S(nx)*mnum*mnum*L*L-8*S(ny)*mnum*L**3-(3*S(nx)**2+4*S(R(aa))*S(h)**4)*L**4
    q,rem=raw.quo_rem(S(h)**6);assert not rem;q=S(q);assert q.degree()==6
    vals=[QQ(c) for p in q.list() for c in p.coefficients()];den=lcm(c.denominator() for c in vals);content=gcd(ZZ(c*den) for c in vals)
    q=S(q*QQ(den)/content);scale=QQ(content)/den
    sparse=[]
    for i,p in enumerate(q.list()):
        for (j,k),c in p.dict().items():sparse.append({'T':i,'u':int(j),'v':int(k),'coefficient':str(c)})
    sparse.sort(key=lambda r:(r['T'],r['u'],r['v']))
    fixed=R([p(0,0) for p in q.list()]);assert fixed.degree()==6 and fixed.gcd(fixed.derivative()).degree()==0
    # The distinct, fibrewise pointed curve has genus1 and a built-in Q point.
    M=PolynomialRing(E.base_ring(),'lambda');lam=M.gen();F=lam**4-6*cx*lam*lam-8*cy*lam-3*cx*cx-4*aa
    c0,c1,c2,c3,c4=F.list();I=12*c4*c0-3*c3*c1+c2*c2;J=72*c4*c2*c0+9*c3*c2*c1-27*c4*c1*c1-27*c3*c3*c0-2*c2**3
    assert I==-48*aa and J==-1728*bb and F.discriminant()==256*E.discriminant()
    assert -27*I==6**4*aa and -27*J==6**6*bb
    # Retrospective diagnostic begins here, after the universal equation exists.
    ap=ART/'det1092_initial_unlock_construction_audit_v1.json';ip=ART/'curve302_recovered_public_span_v1/input.json';rp=ART/'curve302_recovered_public_span_v1/result.json';fp=ART/'curve302_recovered_quotient_local_filtration_v1.json'
    audit=json.loads(ap.read_text())['historical'];span=json.loads(ip.read_text());relations=json.loads(rp.read_text());fil=json.loads(fp.read_text())
    index=span['recovered_points'].index(audit['point_short_model']);word=vector(ZZ,relations['relations'][index]['word']);assert relations['relations'][index]['denominator']==1
    short=EllipticCurve(QQ,span['curve']);public=[short(P) for P in span['public_points']]
    point=sum((n*P for n,P in zip(word,public)),short(0));assert list(map(str,point[:2]))==audit['point_short_model']
    ambient=VectorSpace(GF(2),31);generic=matrix(GF(2),parent['basis_embedding_in_public_D'])
    M17=ambient.subspace(generic.columns());K=ambient.subspace(fil['local_filtration_mod_2']['strict_kernel_public_words']);Sspace=M17+K
    wf=vector(GF(2),word);assert wf not in Sspace and (M17.dimension(),K.dimension(),Sspace.dimension())==(17,10,27)
    char=next(r for r in Sspace.basis_matrix().right_kernel().basis() if r*wf)
    assert Sspace.basis_matrix()*char==0 and char*wf==1
    # Exact elliptic Kummer element in the retained literal302 cubic algebra.
    x,y=map(QQ,audit['point_literal302']);X=4*x;Y=8*y+4*x+4;dd=ZZ(X.denominator()).sqrt();assert dd*dd==X.denominator()
    n=ZZ(X*dd*dd);bnum=ZZ(Y*dd**3)
    E0=EllipticCurve(QQ,[a(0) for a in E.a_invariants()]);assert E0([x,y])
    cubic=R([64*E0.a6()+16,16*E0.a4()+8,5,1]);beta=R([n,-dd*dd])
    assert Y*Y==cubic(X) and cubic.resultant(beta)==bnum*bnum
    def rec(f):return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
    return {'classification':'verified application and retrospective diagnostic','status':'PASS_UNIVERSAL_GENUS2_FAMILY_AND_DESCENT_OBJECT_PREFLIGHT',
        'universal_genus2':{'base_parameters':['u','v'],'curve_coordinate':'T','equation':'s^2=scale*q(T;u,v)',
            'scale':str(scale),'sparse_q':sparse,'degree_T':6,'generic_genus':2,'fixed_0_0_squarefree':True,
            'branch_identity':'f2^4*F_T(-f1/f2+a1/2)=h^6*scale*q','h':list(map(str,h.list())),
            'selection_boundary':'The historical centre is calibrated; the universal equation uses no exceptional point or successful coordinate.'},
        'fibrewise_pointed_curve':{'equation':'w^2=lambda^4-6cx(t)*lambda^2-8cy(t)*lambda-3cx(t)^2-4a(t)',
            'coefficients':[rec(c) for c in F.list()],'genus':1,'Jacobian':'Q-isomorphic to E_t; scale6 gives coefficients -27I,-27J',
            'rational_infinity_points':[[1,0,1],[1,0,-1]],'everywhere_locally_soluble':'YES for every smooth specialization, since rationally pointed',
            'covering_class':'The degree2 divisor O+C gives the known Kummer lift delta(C); its torsor image is0 and its class modulo delta(M17) is0.'},
        'first_witness_diagnostic':{'public_word':list(map(int,word)),'strict_at_S':False,'strict_after_M17_translation':False,
            'separating_character':list(map(int,char)),'character_on_generic_and_strict_kernel':0,'character_on_first_word':1,
            'elliptic_Kummer_cubic':list(map(str,cubic.list())),'beta_coefficients':list(map(str,beta.list())),'norm_square_root':str(bnum)},
        'unresolved':'No target-blind map from elliptic fibre t to RR parameters(u,v) has been specified. No genus2 Selmer/control panel, full elliptic Selmer group, or identification with MW16 reference ideal classes is claimed.',
        'limits':{'universal_sextics':1,'point_searches':0,'Selmer_runs':0,'class_group_runs':0,'pilot_changes':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,np,loader,ap,ip,rp,fp,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'sextic terms',len(d['universal_genus2']['sparse_q']),flush=True)

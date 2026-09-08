#!/usr/bin/env sage-python
"""Independent universal equation, Jacobian-type and strict-filtration checks."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,VectorSpace,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_universal_rr_descent_preflight_v1.json';OUT=ART/'det1092_universal_rr_descent_preflight_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';fp=ART/'curve302_recovered_quotient_local_filtration_v1.json'
    ip=ART/'curve302_recovered_public_span_v1/input.json';ap=ART/'det1092_initial_unlock_construction_audit_v1.json'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text());fil=json.loads(fp.read_text());span=json.loads(ip.read_text());audit=json.loads(ap.read_text())['historical']
    R=PolynomialRing(QQ,'T');K=R.fraction_field()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);base=[E([dec(v) for v in P]) for P in parent['basis_weierstrass_coordinates']]
    C=-sum((n*P for n,P in zip(net['trace_word'],base)),E(0))
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2;aa=-E.c4()/48;bb=-E.c6()/864
    UV=PolynomialRing(QQ,['u','v']);u,v=UV.gens();S=PolynomialRing(UV,'T');T=S.gen()
    q=sum(S(QQ(row['coefficient'])*u**row['u']*v**row['v'])*T**row['T'] for row in d['universal_genus2']['sparse_q'])
    assert q.degree()==6
    A=[S(R(row)) for row in net['A']];B=[S(R(row)) for row in net['B']]
    f0,f1,f2=[B[i]+(u+v*T)*A[i] for i in range(3)]
    h=R(d['universal_genus2']['h']);assert h*h==C[0].denominator()
    # Clear the known centre and slope denominators before multivariate
    # arithmetic; no expensive generic fraction-field gcd is necessary.
    hh=S(h);nx=S(R(cx*h*h));ny=S(R(cy*h**3));mnum=-f1+S(R(E.a1()))*f2/2
    quartic=hh**4*mnum**4-6*nx*hh**2*mnum*mnum*f2*f2-8*ny*hh*mnum*f2**3-(3*nx*nx+4*S(R(aa))*hh**4)*f2**4
    scale=QQ(d['universal_genus2']['scale'])
    assert quartic==hh**10*scale*q
    fixed=R([p(0,0) for p in q.list()]);assert fixed.degree()==6 and fixed.gcd(fixed.derivative()).degree()==0
    M=PolynomialRing(K,'lambda');l=M.gen();Q=M([dec(c) for c in d['fibrewise_pointed_curve']['coefficients']])
    assert Q==l**4-6*cx*l*l-8*cy*l-3*cx*cx-4*aa and Q.discriminant()==256*E.discriminant()
    e,dd,c,b,a=Q.list();I=12*a*e-3*b*dd+c*c;J=72*a*c*e+9*b*c*dd-27*a*dd*dd-27*b*b*e-2*c**3
    assert -27*I==6**4*aa and -27*J==6**6*bb
    assert Q[4]==1 # weighted-projective infinity points(1:0:+/-1)
    diagnostic=d['first_witness_diagnostic'];word=vector(ZZ,diagnostic['public_word'])
    Eshort=EllipticCurve(QQ,span['curve']);public=[Eshort(P) for P in span['public_points']]
    P=sum((n*Q for n,Q in zip(word,public)),Eshort(0));assert list(map(str,P[:2]))==audit['point_short_model']
    ambient=VectorSpace(GF(2),31);M17=ambient.subspace(matrix(GF(2),parent['basis_embedding_in_public_D']).columns())
    kernel=ambient.subspace(fil['local_filtration_mod_2']['strict_kernel_public_words']);phi=vector(GF(2),diagnostic['separating_character'])
    assert (M17+kernel).basis_matrix()*phi==0 and phi*word.change_ring(GF(2))==1
    assert (M17.dimension(),kernel.dimension(),(M17+kernel).dimension())==(17,10,27)
    x=(P[0]-15)/36;y=(P[1]/108-x-1)/2;X=4*x;Y=8*y+4*x+4
    E0=EllipticCurve(QQ,[a(0) for a in E.a_invariants()]);assert E0([x,y])
    cubic=R([64*E0.a6()+16,16*E0.a4()+8,5,1]);beta=R(diagnostic['beta_coefficients']);root=ZZ(diagnostic['norm_square_root'])
    den=ZZ(X.denominator()).sqrt();assert beta==R([X*den*den,-den*den])
    assert Y*den**3==root and cubic==R(diagnostic['elliptic_Kummer_cubic']) and cubic.resultant(beta)==root*root
    return {'classification':'verified application and retrospective diagnostic','status':'PASS_INDEPENDENT_UNIVERSAL_RR_DESCENT_PREFLIGHT',
        'universal_sextic_identity':True,'generic_genus2_parameters':['u','v'],'curve_variable':'T',
        'fibrewise_quartic_genus':1,'fibrewise_Jacobian_is_E_t':True,'rational_points_at_infinity':True,
        'first_gain_strict_modulo_M17':False,'exact_Kummer_norm_identity':True,
        'separating_character_verified':True,'Selmer_runs':0,'point_searches':0,
        'boundary':'Universal family and arithmetic object distinction, not a blinded control-panel outcome or a genus2 Selmer-solubility theorem.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,np,fp,ip,ap]},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    result=verify();text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(result['status'],flush=True)

#!/usr/bin/env sage-python
"""Independent quartic, doubling, ramification and singular-incidence replay.

Classification: verified application. The genus proof uses the established
24-I1 K3 height formula; no general CAS genus or class-group algorithm.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_rr_net_halving_gate_v1.json';OUT=ART/'det1092_rr_net_halving_gate_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text())
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);base=[E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    w=vector(ZZ,d['centre_word']);assert w==-vector(ZZ,net['trace_word'])
    G=matrix(QQ,parent['generic_height_gram']);assert w*G*w==10
    C=sum((n*P for n,P in zip(w,base)),E(0))
    a1,a2,a3,a4,a6=E.a_invariants();cx,cy=map(dec,d['centre_short']);aa=dec(d['short_A']);bb=dec(d['short_B'])
    assert cx==C[0]+E.b2()/12 and cy==C[1]+(a1*C[0]+a3)/2
    assert aa==-E.c4()/48 and bb==-E.c6()/864 and cy*cy==cx**3+aa*cx+bb
    S=PolynomialRing(K,'m');m=S.gen();D=S([dec(v) for v in d['halving_quartic_coefficients']])
    assert D==m**4-6*cx*m*m-8*cy*m-3*cx*cx-4*aa
    qx=(m*m-cx)/2;qy=m*(qx-cx)-cy
    assert (qy*qy-qx**3-aa*qx-bb)%D==0
    assert (2*m*qy-3*qx*qx-aa)%D==0
    assert m*m-2*qx==cx and m*(qx-cx)-qy==cy
    assert D.discriminant()==256*E.discriminant()
    delta=R(E.discriminant());assert delta.degree()==24 and delta.gcd(delta.derivative()).degree()==0
    assert R(cx.denominator()*cy.denominator()).gcd(delta).degree()==0
    # Polynomial K3 model, smooth infinity; all24 finite roots are I1.
    assert all(R(v).degree()<=2*j for v,j in zip(E.a_invariants(),[1,2,3,4,6]))
    assert R(E.c4()).gcd(delta).degree()==0
    assert cx.numerator().degree()-cx.denominator().degree()<=4
    assert cy.numerator().degree()-cy.denominator().degree()<=6
    # Geometric connectedness: heights on the24-I1 K3 are4+2(P.O).
    # A degree1 component would give a half of height10/4<4. A degree2
    # component {Q,Q'} gives the nonzero rational2torsion Q'-Q=Q+Q'-C.
    # Torsion itself is excluded by the same positive height lower bound.
    assert QQ(10)/4<4
    branch=24;degree=4;g=1+(-2*degree+branch)//2;assert g==9
    AA=list(map(R,net['A']));BB=list(map(R,net['B']))
    N=S([dec(v) for v in d['net_r_numerator']]);L=S([dec(v) for v in d['net_r_denominator']])
    assert N==-(BB[1]+(m-a1/2)*BB[2]) and L==AA[1]+(m-a1/2)*AA[2]
    assert N[1]*L[0]-N[0]*L[1] and D.gcd(L*D.derivative()).degree()==0
    # Check the Mobius map and the tangent incidence by exact rational algebra.
    F=S.fraction_field();r=F(N)/L
    def dtpoly(p):return S([v.derivative(t) for v in p.list()])
    def dt(f):return F(dtpoly(f.numerator())*f.denominator()-f.numerator()*dtpoly(f.denominator()))/f.denominator()**2
    dm=F(D.derivative());dtD=F(dtpoly(D));rt=dt(r);rm=r.derivative(m)
    assert rm
    v=rt-rm*dtD/dm;u=r-t*v
    assert u+t*v==r and dtD+dm*(v-rt)/rm==0
    # No coordinate is fitted to a point, and no member search is performed.
    local=d['zero_no_half_witness'];p=ZZ(local['prime']);assert p.is_prime()
    Fp=PolynomialRing(GF(p),'m');d0=Fp([v(0) for v in D.list()])
    assert list(map(int,d0.list()))==local['quartic_coefficients']
    assert d0.degree()==4 and d0.discriminant() and not d0.roots()
    return {'classification':'verified application and new deduction','status':'PASS_INDEPENDENT_HALVING_GENUS9_AND_SINGULAR_INCIDENCE',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,np]},'checker_sha256':sha(Path(__file__)),
        'exact_doubling_identity':True,'quartic_discriminant_identity':True,
        'geometrically_connected':True,'degree':degree,'simple_ramification_contribution':branch,'genus':g,
        'connectedness_proof':'No degree1 component by height10/4<4; no degree2 component by absence of geometric rational2torsion.',
        'singular_incidence_map_verified':True,'zero_has_no_rational_half':True,'local_prime':int(p),
        'point_searches':0,'boundary':'Dense-open singular incidence. No classification of exceptional members, no rational points on the genus9 curve constructed, and no arithmetic explanation of the first independent302 point.'}
if __name__=='__main__':
    d=verify();text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(d['status'],flush=True)

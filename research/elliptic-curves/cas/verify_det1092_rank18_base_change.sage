#!/usr/bin/env sage-python
"""Independent coefficient/height verification; does not import the producer.

Classification: verified application of the established elliptic-surface
intersection height formula. No numerical height, point search or conic solver.
"""
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, gcd, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_orbit8044_rank18_base_change_v2.json'
OUT=ART/'det1092_orbit8044_rank18_base_change_replay_v2.json'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for name,expected in d['inputs'].items(): assert sha(ROOT/name)==expected
    parent=json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
    R=PolynomialRing(QQ,'t'); K=R.fraction_field()
    U=PolynomialRing(QQ,'u'); F=U.fraction_field()
    def decode(v,ring,field): return field(ring(v['numerator']))/ring(v['denominator'])
    a_invariants=[decode(v,R,K) for v in parent['a_invariants']]
    E=EllipticCurve(K,a_invariants)
    basis=[E([decode(v,R,K) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    q=R(d['curve_over_Q']['q_coefficients']); h=R(d['splitting']['discriminant_square_factor'])
    c,b,a=[decode(v,R,K) for v in d['lift']['residual_coefficients']]
    f0,f1,f2=[R(v) for v in d['lift']['line_coefficients']]
    assert q.degree()==2 and q.discriminant()!=0 and not q.is_square()
    assert b*b-4*a*c==q*h*h
    T=decode(d['parametrization']['t_of_u'],U,F); S=decode(d['parametrization']['s_of_u'],U,F)
    t0,s0=map(QQ,d['curve_over_Q']['rational_point'])
    assert s0*s0==q(t0) and S*S==q(T)
    assert (S-s0)/(T-t0)==U.gen()
    assert max(T.numerator().degree(),T.denominator().degree())==2
    delta=R(E.discriminant())
    assert delta.degree()==24 and gcd(delta,delta.derivative()).degree()==0
    assert gcd(q,delta).degree()==0
    def ev(v): return F(v.numerator()(T))/v.denominator()(T)
    old=[(ev(p[0]),ev(p[1])) for p in basis]
    x=decode(d['lift']['x_of_u'],U,F);y=decode(d['lift']['y_of_u'],U,F)
    assert x==(-ev(b)+h(T)*S)/(2*ev(a)) and f0(T)+f1(T)*x+f2(T)*y==0
    assert ev(a)*x*x+ev(b)*x+ev(c)==0
    w=vector(ZZ,d['lift']['trace_word'])
    trace=-sum((n*p for n,p in zip(w,basis)),E(0))
    # Recheck the RR identity and full elimination, including the removed root.
    assert f0+f1*trace[0]+f2*trace[1]==0
    V=PolynomialRing(K,'z');z=V.gen();a1,a2,a3,a4,a6=a_invariants
    line=f0+f1*z
    elimination=line**2-a1*z*line*f2-a3*line*f2-f2**2*(z**3+a2*z**2+a4*z+a6)
    assert elimination==(z-trace[0])*(a*z*z+b*z+c)
    pulled=EllipticCurve(F,[ev(v) for v in a_invariants]);Q=pulled([x,y])
    xb=(-ev(b)-h(T)*S)/(2*ev(a));yb=-(f0(T)+f1(T)*xb)/f2(T)
    Qb=pulled([xb,yb])
    assert Q+Qb==sum((n*pulled(p) for n,p in zip(w,old)),pulled(0))
    D=T.denominator();X=D**4*x;Y=D**6*y
    assert X==U(d['polynomial_model']['X']) and Y==U(d['polynomial_model']['Y'])
    assert U(X).degree()<=8 and U(Y).degree()<=12
    assert [D**(2*k)*ev(v) for k,v in zip([1,2,3,4,6],a_invariants)]==[U(v) for v in d['polynomial_model']['a_invariants']]
    G=matrix(QQ,parent['generic_height_gram']);cross=G*w
    H=(2*G).augment(matrix(QQ,17,1,list(cross))).stack(matrix(QQ,1,18,list(cross)+[8]))
    assert w*G*w==10 and G.det()==1092 and H==matrix(QQ,d['independence']['gram'])
    assert H.is_positive_definite() and H.det()==3*2**17*1092
    # Smooth branch fibres preserve semistability: 24I1 -> 48I1, chi 2 -> 4.
    # Polynomial X,Y of degrees <=8,12 show Q.O=0, including infinity.
    # Height(Q)=8, height(trace)=20, so height(Q-Qb)=4*8-20=12.
    assert d['independence']['anti_invariant_Q_minus_sigmaQ_height']==12
    assert not q(0).is_square()
    excluded=R(d['splitting']['excluded_t_polynomial'])
    for f in [delta,q,h,f2,R(a.numerator()),*[R(v.denominator()) for v in [a,b,c,*a_invariants]]]:
        radical=R(f / gcd(f,f.derivative()))
        assert (excluded % radical).is_zero()
    return {'classification':'verified application','status':'PASS_INDEPENDENT_COEFFICIENT_AND_INTERSECTION_HEIGHT_REPLAY',
        'inputs':{str(CERT.relative_to(ROOT)):sha(CERT)},'checker_sha256':sha(Path(__file__)),
        'rational_parametrization_identity':True,'RR_elimination_and_trace_identity':True,
        'all_17_inherited_section_identities':True,'new_section_identity':True,
        'semistable_fibres':48,'Q_dot_O':0,'anti_invariant_height':12,'positive_definite_gram_determinant':str(H.det()),
        'rank_lower_bound_over_Q_u':18,'point_searches_run':0,
        'literature_dependency':'Shioda intersection height formula and degree scaling; parent rank17 certificate remains a dependency.'}

if __name__=='__main__':
    result=verify();text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists(): assert OUT.read_text()==text
    else: OUT.write_text(text)
    print(result['status'],flush=True)

#!/usr/bin/env sage-python
"""Independent exact replay of the fixed two-cover rank/genus gate.

Classification: verified application. Reads coefficients, performs no conic
solver, rational point search, orbit selection, or population calculation.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,gcd,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
CERT=ART/'det1092_two_cover_rank19_genus_gate_v1.json'
OUT=ART/'det1092_two_cover_rank19_genus_gate_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify():
    d=json.loads(CERT.read_text())
    for name,digest in d['inputs'].items():assert sha(ROOT/name)==digest
    prior=runpy.run_path(str(CAS/'verify_det1092_rank18_base_change.sage'))['verify']()
    assert prior['rank_lower_bound_over_Q_u']==18
    parent=json.loads((ART/'curve302_recovered_mw17_parent_v1.json').read_text())
    rr=json.loads((ART/'det1092_initial_unlock_bisection_obstruction_v1.json').read_text())
    first=json.loads((ART/'det1092_orbit8044_rank18_base_change_v2.json').read_text())
    R=PolynomialRing(QQ,'t');K=R.fraction_field();U=PolynomialRing(QQ,'u');F=U.fraction_field()
    def rat(v,ring,field):return field(ring(v['numerator']))/ring(v['denominator'])
    E=EllipticCurve(K,[rat(v,R,K) for v in parent['a_invariants']])
    base=[E([rat(v,R,K) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    a2=d['second_cover'];q1=R(d['common_cover']['q1_coefficients']);q2=R(a2['q_coefficients'])
    assert q1==R(first['curve_over_Q']['q_coefficients']) and q2==R(d['common_cover']['q2_coefficients'])
    assert q1.degree()==q2.degree()==2 and q1.is_irreducible() and q2.is_irreducible()
    assert gcd(q1,q2).degree()==0 and gcd(q1*q2,(q1*q2).derivative()).degree()==0
    delta=R(E.discriminant());assert gcd(q1*q2,delta).degree()==0
    def olddec(v):return K(R(v['numerator_coefficients_low_to_high']))/R(v['denominator_coefficients_low_to_high'])
    c,b,a=map(olddec,rr['residual_coefficients']);h=R(a2['discriminant_square_factor'])
    assert b*b-4*a*c==q2*h*h
    f0,f1,f2=[R(rr['line_coefficients'][k]) for k in ['f0','f1','f2']]
    T=rat(a2['t_of_u'],U,F);S=rat(a2['s_of_u'],U,F);t0,s0=map(QQ,a2['rational_conic_point'])
    assert s0*s0==q2(t0) and S*S==q2(T) and (S-s0)/(T-t0)==U.gen()
    assert max(T.numerator().degree(),T.denominator().degree())==2
    def ev(v):return F(v.numerator()(T))/v.denominator()(T)
    pulled=EllipticCurve(F,[ev(v) for v in E.a_invariants()])
    x=rat(a2['x_of_u'],U,F);y=rat(a2['y_of_u'],U,F)
    assert x==(-ev(b)+h(T)*S)/(2*ev(a)) and f0(T)+f1(T)*x+f2(T)*y==0
    assert ev(a)*x*x+ev(b)*x+ev(c)==0
    Q=pulled([x,y]);xb=(-ev(b)-h(T)*S)/(2*ev(a));yb=-(f0(T)+f1(T)*xb)/f2(T);Qbar=pulled([xb,yb])
    w1=vector(ZZ,first['lift']['trace_word']);w2=vector(ZZ,a2['trace_word']);assert w2==-vector(ZZ,rr['trace_word'])
    assert Q+Qbar==sum((n*pulled([ev(P[0]),ev(P[1])]) for n,P in zip(w2,base)),pulled(0))
    X=U(T.denominator()**4*x);Y=U(T.denominator()**6*y)
    assert X==U(a2['X_polynomial']) and Y==U(a2['Y_polynomial']) and X.degree()<=8 and Y.degree()<=12
    G=matrix(QQ,parent['generic_height_gram']);assert G.det()==1092 and w1*G*w1==w2*G*w2==10
    assert w1*G*w2==0
    cross=matrix(QQ,17,2,lambda i,j:2*(G*[w1,w2][j])[i]);bottom=16*matrix.identity(QQ,2)
    H=(4*G).augment(cross).stack(cross.transpose().augment(bottom))
    assert H==matrix(QQ,d['independence']['gram']) and H.is_positive_definite()
    assert bottom-cross.transpose()*(4*G).inverse()*cross==6*matrix.identity(QQ,2)
    assert H.det()==36*4**17*1092
    assert d['independence']['anti_invariant_heights_on_common_cover']==[24,24]
    assert d['independence']['deck_characters']==[[-1,1],[1,-1]]
    # Four simple disjoint branch points; each has two ramification points of
    # index2 in the degree4 connected cover. 2g-2=-8+8=0.
    degree=4;ramification=4*2*(2-1);genus=1+(-2*degree+ramification)//2
    assert genus==d['common_cover']['genus']==1
    for row in d['bounded_local_checks']['rows']:
        p,k=row['prime'],row['exponent'];m=p**k;sq={x*x%m for x in range(m)}
        def condition(n,z):return all(int(q[0]*z*z+q[1]*n*z+q[2]*n*n)%m in sq for q in [q1,q2])
        count=sum(condition(n,1) for n in range(m))+sum(condition(1,z) for z in range(0,m,p))
        assert count==row['simultaneous_square_classes']>0
    return {'classification':'verified application','status':'PASS_INDEPENDENT_TWO_SECTION_HEIGHT_AND_GENUS_REPLAY',
        'inputs':{str(CERT.relative_to(ROOT)):sha(CERT)},'checker_sha256':sha(Path(__file__)),
        'both_section_and_trace_identities':True,'both_single_cover_anti_invariant_heights':12,
        'common_cover_anti_invariant_heights':[24,24],'independent_deck_characters':True,
        'rank_lower_bound_over_Q_common_cover':19,'gram_determinant':str(H.det()),
        'common_cover_genus':genus,'rational_parametrization':'IMPOSSIBLE','common_rational_point':'UNKNOWN',
        'point_searches_run':0,'claim_boundary':'Function-field rank19 and obstruction to parametrizing simultaneous splitting by P1. No rational-point existence claim for the genus1 base.'}

if __name__=='__main__':
    d=verify();text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(d['status'],flush=True)

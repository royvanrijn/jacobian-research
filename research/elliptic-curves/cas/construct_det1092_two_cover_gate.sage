#!/usr/bin/env sage-python
"""Two fixed covers only: independent sections and the rational-base obstruction.

Classification: verified application and new deduction. Uses completed orbit
8044 and retrospective initial-centre orbit127449. No new orbit selection,
point search, parameter population, or class-group computation.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,Conic,gcd,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
FIRST=ART/'det1092_orbit8044_rank18_base_change_v2.json'
SECOND=ART/'det1092_initial_unlock_bisection_obstruction_v1.json'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
REPLAY=ART/'det1092_orbit8044_rank18_base_change_replay_v2.json'
OUT=ART/'det1092_two_cover_rank19_genus_gate_v1.json'
PRIME_POWERS=((2,10),(3,6),(5,4),(7,3),(11,2),(13,2),(17,2),(19,2),(23,2),(29,2),(31,2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(v):return list(map(str,v.list())) or ['0']
def rec(v):return {'numerator':poly(v.numerator()),'denominator':poly(v.denominator())}

def build():
    first,second,parent,replay=[json.loads(p.read_text()) for p in (FIRST,SECOND,PARENT,REPLAY)]
    assert replay['status']=='PASS_INDEPENDENT_COEFFICIENT_AND_INTERSECTION_HEIGHT_REPLAY'
    assert replay['inputs'][str(FIRST.relative_to(ROOT))]==sha(FIRST)
    assert second['status']=='PASS_EXACT_FIRST_CENTRE_BISECTION_NONSPLIT_OBSTRUCTION'
    E,basis,_=runpy.run_path(str(CAS/'load_curve302_recovered_parent.sage'))['load_curve302_recovered_parent'](PARENT)
    K=E.base_ring();R=K.ring()
    def dec(v):return K(R(v['numerator_coefficients_low_to_high']))/R(v['denominator_coefficients_low_to_high'])
    c,b,a=map(dec,second['residual_coefficients'])
    f0,f1,f2=[R(second['line_coefficients'][key]) for key in ['f0','f1','f2']]
    disc=R(b*b-4*a*c);fac=disc.factor();q=R(fac.unit());h=R(1)
    for f,n in fac:q*=f**(n%2);h*=f**(n//2)
    content=gcd(q.list());assert content>0 and content.is_square()
    q=R(q/content);h*=content.sqrt()
    assert disc==q*h*h
    q1=R(first['curve_over_Q']['q_coefficients'])
    delta=R(E.discriminant())
    assert q.degree()==q1.degree()==2 and q.is_irreducible() and q1.is_irreducible()
    assert gcd(q,q1).degree()==0 and gcd(q,delta).degree()==0
    print('CHECKPOINT disjoint squarefree branch pairs; connected V4 cover of genus1',flush=True)
    conic=Conic(QQ,[-q[0],-q[1],0,-q[2],0,1]);ok,P=conic.has_rational_point(point=True)
    assert ok and P[0]!=0
    t0,s0=QQ(P[1]/P[0]),QQ(P[2]/P[0]);assert s0*s0==q(t0)
    U=PolynomialRing(QQ,'u');F=U.fraction_field();u=U.gen()
    T=F(t0)+(q.derivative()(t0)-2*s0*u)/(u*u-q[2]);S=F(s0)+u*(T-t0)
    assert S*S==q(T) and max(T.numerator().degree(),T.denominator().degree())==2
    def ev(v):return F(v.numerator()(T))/v.denominator()(T)
    pulled=EllipticCurve(F,[ev(v) for v in E.a_invariants()])
    x=(-ev(b)+h(T)*S)/(2*ev(a));y=-(f0(T)+f1(T)*x)/f2(T)
    xb=(-ev(b)-h(T)*S)/(2*ev(a));yb=-(f0(T)+f1(T)*xb)/f2(T)
    Q=pulled([x,y]);Qb=pulled([xb,yb]);w2=-vector(ZZ,second['trace_word'])
    old=[pulled([ev(P[0]),ev(P[1])]) for P in basis]
    assert Q+Qb==sum((n*P for n,P in zip(w2,old)),pulled(0))
    D=T.denominator();X=U(D**4*x);Y=U(D**6*y)
    assert X.degree()<=8 and Y.degree()<=12
    G=matrix(QQ,parent['generic_height_gram']);w1=vector(ZZ,first['lift']['trace_word'])
    assert w1*G*w1==w2*G*w2==10
    dot=w1*G*w2
    # Over degree4 common cover: height(Qi)=16, height(2Qi-Wi)=24.
    # Distinct nontrivial deck characters force anti-invariant cross pairing0.
    cross=matrix(QQ,17,2,lambda i,j:2*(G*[w1,w2][j])[i])
    bottom=matrix(QQ,[[16,dot],[dot,16]])
    H=(4*G).augment(cross).stack(cross.transpose().augment(bottom))
    schur=bottom-cross.transpose()*(4*G).inverse()*cross
    assert schur==6*matrix.identity(QQ,2) and H.is_positive_definite()
    assert H.det()==36*4**17*1092
    print('CHECKPOINT both positive-height anti-invariant directions; rank19 Gram verified',flush=True)
    local=[]
    for p,k in PRIME_POWERS:
        m=p**k;squares={x*x%m for x in range(m)}
        qs=[[int(v)%m for v in f.list()] for f in [q1,q]]
        def ok(n,d):return all((v[0]*d*d+v[1]*n*d+v[2]*n*n)%m in squares for v in qs)
        count=sum(ok(n,1) for n in range(m))+sum(ok(1,d) for d in range(0,m,p))
        local.append({'prime':p,'exponent':k,'modulus':m,'primitive_projective_base_classes':m+m//p,'simultaneous_square_classes':count})
    assert all(row['simultaneous_square_classes']>0 for row in local)
    return {'schema':'elliptic-curves.det1092-two-cover-rank19-genus-gate.v1','classification':'verified application',
        'status':'PASS_RANK19_FUNCTION_FIELD_AND_GENUS1_RATIONAL_PARAMETRIZATION_OBSTRUCTION',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [FIRST,SECOND,PARENT,REPLAY,CAS/'load_curve302_recovered_parent.sage',Path(__file__)]},
        'limits':{'cover_ids':[8044,127449],'new_cover_selections':0,'additional_quadratic_factorizations':1,'additional_conic_parametrizations':1,'rational_point_searches':0,'parameter_sweeps':0,'finite_congruence_moduli':[p**k for p,k in PRIME_POWERS]},
        'second_cover':{'q_coefficients':poly(q),'discriminant_square_factor':poly(h),'rational_conic_point':[str(t0),str(s0)],'t_of_u':rec(T),'s_of_u':rec(S),'x_of_u':rec(x),'y_of_u':rec(y),'X_polynomial':poly(X),'Y_polynomial':poly(Y),'trace_word':list(map(int,w2)),'Q_height':8,'anti_invariant_height':12,'constant_extension_required':False},
        'common_cover':{'equations':['s1^2=q1(t)','s2^2=q2(t)'],'q1_coefficients':poly(q1),'q2_coefficients':poly(q),'function_field_degree':4,'geometrically_connected':True,'branch_gcd_degree':0,'geometric_branch_points':4,'ramification_index':2,'ramification_points_per_branch':2,'genus':1,'field_of_definition':'Q','rational_point_status':'UNKNOWN','rational_parametrization':'IMPOSSIBLE even after constant-field extension, by Riemann-Hurwitz','splitting_condition':'Both q1(t) and q2(t) are rational squares, with chart denominators nonzero.'},
        'independence':{'inherited_rank':17,'trace_pairing_on_parent':str(dot),'trace_words':[list(map(int,w1)),list(map(int,w2))],'anti_invariant_heights_on_common_cover':[24,24],'deck_characters':[[-1,1],[1,-1]],'anti_invariant_cross_pairing':0,'gram':[list(map(str,row)) for row in H.rows()],'schur_complement':[list(map(str,row)) for row in schur.rows()],'gram_determinant':str(H.det()),'rank_lower_bound_over_Q_common_cover':19,'argument':'Positive heights prove both anti-invariant points non-torsion. Their distinct deck characters make them independent of one another and of the invariant old17. Distinct squareclasses alone are not used as a rank certificate.'},
        'bounded_local_checks':{'classification':'verified application','rows':local,'conclusion':'No obstruction at these finite moduli. This does not prove Q_p solubility, everywhere local solubility, a global rational point, or a positive-rank Jacobian.'},
        'general_branch_gate':{'classification':'new deduction','hypotheses':'k independent quadratic characters, each with two simple geometric branch points; the k branch pairs are disjoint; characteristic0.','genus_formula':'g=1+2^(k-1)*(k-2)','examples':{'1':0,'2':1,'3':5,'4':17},'meaning':'Combining positive-height character sections raises function-field rank but typically removes a rational parameter line. To retain a rational base, seek controlled branch sharing or several proved independent sections on one character.'},
        'literature':[{'classification':'established literature','title':'Riemann-Hurwitz','url':'https://stacks.math.columbia.edu/tag/0C1B'},{'classification':'established literature','title':'Schuett--Shioda, Elliptic Surfaces','url':'https://arxiv.org/abs/0907.0298'}],
        'boundary':'No claim of a rank19 rationally parametrized family, infinitely many rational points on the common genus1 base, specialized rank19, or new high-rank curve. The original eight-fibre pilot is untouched.'}

if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    payload=build();OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(payload['status'],flush=True)

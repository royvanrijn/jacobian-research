#!/usr/bin/env sage-python
"""One bounded verified application: lift the already constructed orbit 8044.

No orbit census, parameter sweep, point search, or V3 mutation. At most one
quadratic discriminant factorization and one conic parametrization are used.
Run from research/: sage -python elliptic-curves/cas/construct_det1092_rank18_base_change.sage
"""
import hashlib
import json
import runpy
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, Conic, gcd, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
SOURCE = ART / 'curve302_parent_cheapest_lattice_bisection_v1.json'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
LOADER = ROOT / 'elliptic-curves/cas/load_curve302_recovered_parent.sage'
OUTPUT = ART / 'det1092_orbit8044_rank18_base_change_v2.json'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p): return [str(v) for v in p.list()] or ['0']
def rat(f): return {'numerator': poly(f.numerator()), 'denominator': poly(f.denominator())}

def build():
    source = json.loads(SOURCE.read_text())
    assert source['status'] == 'PASS_EXACT_NONSPLIT_CHEAPEST_LATTICE_BISECTION'
    assert source['selection']['orbit_mask'] == 8044
    parent = json.loads(PARENT.read_text())
    assert source['inputs'][str(PARENT.relative_to(ROOT))] == sha(PARENT)
    E, basis, _ = runpy.run_path(str(LOADER))['load_curve302_recovered_parent'](PARENT)
    K = E.base_ring(); R = K.ring(); t = R.gen()
    def decode(v):
        return K(R(v['numerator_coefficients_low_to_high'])) / R(v['denominator_coefficients_low_to_high'])
    c,b,a = map(decode, source['RR_construction']['generic_residual_x_coefficients_low_to_high'])
    f0,f1,f2 = [R(source['RR_construction'][key+'_coefficients_low_to_high']) for key in ('f0','f1','f2')]
    D = R(b*b-4*a*c)
    factors = D.factor()
    q = R(factors.unit()); h = R(1)
    for f,e in factors:
        q *= f**(e%2); h *= f**(e//2)
    content = gcd(q.list())
    assert content > 0 and content.is_square()
    q /= content; h *= content.sqrt()
    q = R(q); h = R(h)
    assert D == q*h*h and q.degree() == 2 and q.discriminant() != 0
    assert gcd(q,R(E.discriminant())).degree() == 0
    assert R(E.discriminant()).degree() == 24
    print('CHECKPOINT quadratic cover and smooth branch fibres verified', flush=True)
    # Coordinates on this projective conic are [Z:T:S], with S^2=q(T/Z) Z^2.
    conic = Conic(QQ,[-q[0],-q[1],0,-q[2],0,1])
    ok,point = conic.has_rational_point(point=True)
    assert ok and point[0] != 0
    t0,s0 = QQ(point[1]/point[0]),QQ(point[2]/point[0])
    assert s0*s0 == q(t0)
    U = PolynomialRing(QQ,'u'); F = U.fraction_field(); u = U.gen()
    T = F(t0)+(q.derivative()(t0)-2*s0*u)/(u*u-q[2])
    S = F(s0)+u*(T-t0)
    assert S*S == q(T)
    assert max(T.numerator().degree(),T.denominator().degree()) == 2
    def subst(v): return F(v.numerator()(T))/v.denominator()(T)
    x = (-subst(b)+h(T)*S)/(2*subst(a))
    y = -(f0(T)+f1(T)*x)/f2(T)
    A = [subst(v) for v in E.a_invariants()]
    pulled = EllipticCurve(F,A)
    Q = pulled([x,y])
    xbar = (-subst(b)-h(T)*S)/(2*subst(a))
    ybar = -(f0(T)+f1(T)*xbar)/f2(T)
    Qbar = pulled([xbar,ybar])
    w = vector(ZZ,source['selection']['bisection_phi_w_parent_MW17_coordinates'])
    old = [pulled([subst(p[0]),subst(p[1])]) for p in basis]
    W = sum((n*p for n,p in zip(w,old)),pulled(0))
    assert Q+Qbar == W
    print('CHECKPOINT lifted section, deck conjugate and exact trace identity verified', flush=True)
    # Clear only the known base-map denominator. This gives the chi=4 model.
    denominator = T.denominator()
    X,Y = denominator**4*x,denominator**6*y
    assert X.denominator().degree() == 0 and Y.denominator().degree() == 0
    X,Y = U(X),U(Y)
    assert X.degree() <= 8 and Y.degree() <= 12
    # At infinity use v^8 X(1/v), v^12 Y(1/v). Finite polynomial coordinates
    # and these degree bounds prove Q.O=0 on the global minimal 48-I1 model.
    weights = [1,2,3,4,6]
    integral_A = [U(denominator**(2*k)*v) for k,v in zip(weights,A)]
    G = matrix(QQ,parent['generic_height_gram'])
    assert G.det() == 1092 and w*G*w == 10
    cross = G*w
    H = (2*G).augment(matrix(QQ,17,1,list(cross))).stack(matrix(QQ,1,18,list(cross)+[8]))
    assert H.is_positive_definite() and H.det() == 3*2**17*1092
    excluded = R(E.discriminant()) * R(a.numerator()) * R(f2) * h * q
    for v in [a,b,c,*E.a_invariants()]: excluded *= R(v.denominator())
    excluded = R(excluded / gcd(excluded,excluded.derivative())).monic()
    print('CHECKPOINT disjoint zero section; positive-definite rank18 Gram verified', flush=True)
    return {
        'schema':'elliptic-curves.det1092-orbit8044-rank18-base-change.v1',
        'status':'PASS_EXPLICIT_Q_RATIONAL_BASE_CHANGE_RANK_AT_LEAST_18',
        'classification':'verified application',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,PARENT,LOADER,Path(__file__))},
        'limits':{'orbits':1,'orbit_ids':[8044],'quadratic_factorizations':1,'conic_parametrizations':1,'point_searches':0,'parameter_sweeps':0},
        'curve_over_Q':{'equation':'s^2=q(t)','q_coefficients':poly(q),'genus':0,'rational_point':[str(t0),str(s0)],'constant_field_extension_required':False},
        'parametrization':{'t_of_u':rat(T),'s_of_u':rat(S),'inverse':'u=(s-s0)/(t-t0)','degree_to_parameter_line':2},
        'splitting':{'criterion':'For rational t off the excluded polynomial and infinity, the two fibre points are rational iff q(t) is a rational square. Ramified fibres q(t)=0 have one rational point when t is rational.','discriminant_square_factor':poly(h),'excluded_t_polynomial':poly(excluded),'also_exclude':'t=infinity in the affine formula; poles of t(u); basepoint t=t0 in the inverse chart','zero_fibre_split':False},
        'lift':{'x_formula':'(-b(t)+h(t)*s)/(2*a(t))','y_formula':'-(f0(t)+f1(t)*x)/f2(t)','residual_coefficients':[rat(v) for v in (c,b,a)],'line_coefficients':[poly(v) for v in (f0,f1,f2)],'x_of_u':rat(x),'y_of_u':rat(y),'deck_action':'s -> -s','trace_word':list(map(int,w)),'trace_identity':'Q+sigma(Q)=sum(w_i P_i)'},
        'polynomial_model':{'base_denominator':poly(denominator),'a_invariants':[poly(v) for v in integral_A],'X':poly(X),'Y':poly(Y),'coordinate_rule':'X=den(t(u))^4*x, Y=den(t(u))^6*y'},
        'independence':{'method':'Shioda intersection height and degree-two base-change scaling','chi':4,'fibres':'48 I1','branch_discriminant_gcd_degree':0,'Q_dot_O':0,'Q_height':8,'trace_height_on_cover':20,'Q_sigmaQ_pairing':2,'anti_invariant_Q_minus_sigmaQ_height':12,'inherited_rank':17,'gram':[list(map(str,row)) for row in H.rows()],'gram_determinant':str(H.det()),'schur_complement':'3','rank_lower_bound':18},
        'rational_points_source':'The displayed rational conic point and line parametrization give infinitely many Q-points. Excluding finitely many parameters preserves infinitude. This is a function-field rank statement; rank18 is not asserted at every rational specialization.',
        'literature':[{'classification':'established literature','source':'Schuett--Shioda, Elliptic Surfaces, sections 8 and 11','url':'https://arxiv.org/abs/0907.0298'}],
        'boundary':'One already constructed rational bisection. No new multisection census or search exposure. Its zero fibre remains nonsplit over Q, so this cover does not explain the initial rational302 discovery. No combined-cover independence claim.'}

if __name__ == '__main__':
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    payload=build()
    OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print('PASS',OUTPUT,flush=True)

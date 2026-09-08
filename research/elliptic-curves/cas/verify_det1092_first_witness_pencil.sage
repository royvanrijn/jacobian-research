#!/usr/bin/env sage-python
"""Independent exact branch identity and modular discriminant reconstruction.

Degree28 is proved by the four-root degeneration, not guessed modulo p.
At191 reconstruct that polynomial from29 direct Sylvester determinants.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_first_witness_pencil_genus_gate_v1.json';OUT=ART/'det1092_first_witness_pencil_genus_gate_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';bp=ART/'det1092_first_centre_rr_net_replay_v1.json'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text());bridge=json.loads(bp.read_text())
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def dec(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[dec(row) for row in parent['a_invariants']]);base=[E([dec(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    C=-sum((n*P for n,P in zip(net['trace_word'],base)),E(0))
    z0=QQ(bridge['witness_coordinate']);uf=bridge['chart_to_net_u'];u0=R(uf['numerator'])(z0)/R(uf['denominator'])(z0)
    assert u0==QQ(d['u0']) and z0==QQ(d['witness_coordinate'])
    V=PolynomialRing(QQ,'v');S=PolynomialRing(V,'t');t=S.gen();v=V.gen()
    q=S([V(row) for row in d['q_coefficients_t_then_v']]);scale=QQ(d['branch_squareclass_scale'])
    assert scale and q.degree()==6 and all(c.denominator()==1 for p in q.list() for c in p.list())
    assert gcd([ZZ(c) for p in q.list() for c in p.list()])==1
    A=[S(R(row)) for row in net['A']];B=[S(R(row)) for row in net['B']]
    f0,f1,f2=[B[i]+(u0+v*t)*A[i] for i in range(3)]
    h=R(d['trace_denominator_root']);assert h*h==C[0].denominator() and h(0)
    assert f2[0] and f2[0].degree()==0
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2;aa=-E.c4()/48
    # A direct rational-function version of the quartic identity, independent
    # of the constructor's denominator-cancelled numerator formula.
    T=PolynomialRing(V.fraction_field(),'t');F=T.fraction_field()
    def ev(f):return F(T(R(f.numerator())))/T(R(f.denominator()))
    slope=-F(T(f1))/T(f2)+ev(E.a1())/2
    H=slope**4-6*ev(cx)*slope**2-8*ev(cy)*slope-3*ev(cx)**2-4*ev(aa)
    assert T(f2)**4*H==T(h)**6*scale*T(q)
    assert max(p.degree() for p in q.list())==4
    assert all(q[i][j]==0 for i in range(7) for j in range(i+1,5))
    qa=R([q[i][4] for i in range(4,7)]);lim=R([q[j][j] for j in range(5)])
    assert qa==R(d['infinity_degeneration']['quadratic']) and lim==R(d['infinity_degeneration']['quartic'])
    assert qa.degree()==2 and qa(0) and qa.gcd(qa.derivative()).degree()==0
    assert lim.degree()==4 and lim.gcd(lim.derivative()).degree()==0
    # These conditions prove exact discriminant degree40-12=28 over Q.
    # Four roots have pairwise difference valuation1 in epsilon; the other
    # two tend to distinct nonzero roots. The sextic leading coefficient is
    # an epsilon-unit, so the product-of-root-differences formula applies.
    chosen=next(r for r in d['modular_checks'] if r['discriminant_degree']==28 and not r['roots'])
    p=ZZ(chosen['prime']);assert p==191 and p.is_prime()
    field=GF(p);vp=PolynomialRing(field,'v');xp=PolynomialRing(field,'t')
    coeff=[vp(row.list()) for row in q.list()];samples=[]
    for a in field:
        values=[c(a) for c in coeff]
        if not values[6]:continue
        f=xp(values);g=f.derivative();fc=list(reversed(f.list()));gc=list(reversed(g.list()))
        rows=[]
        for j in range(5):rows.append([field(0)]*j+fc+[field(0)]*(4-j))
        for j in range(6):rows.append([field(0)]*j+gc+[field(0)]*(5-j))
        syl=matrix(field,rows);assert syl.dimensions()==(11,11)
        value=-syl.det()/values[6];samples.append((a,value))
        if len(samples)==29:break
    assert len(samples)==29
    disc=vp.lagrange_polynomial(samples);assert disc.degree()==28
    assert disc==vp(chosen['coefficients'])
    assert disc.gcd(vp.gen()**p-vp.gen()).degree()==0
    return {'classification':'retrospective verified application and new deduction','status':'PASS_INDEPENDENT_ALL_RATIONAL_AFFINE_MEMBERS_HAVE_GENUS2',
        'exact_branch_identity':True,'global_discriminant_degree':28,
        'degree_proof':'Four-root epsilon cluster contributes12 to degree40 homogeneous discriminant scaling.',
        'prime':int(p),'independent_sylvester_determinants':29,'modular_discriminant':list(map(int,disc.list())),
        'modular_root_gcd_degree':0,'no_rational_discriminant_root':True,
        'affine_members_geometric_genus':2,'pencil_infinity':'Old bisection plus302 fibre; reducible and not a point construction.',
        'claim_boundary':'This diagnoses only the fixed next RR system through the known first chart witness. It is not a discovery rule or an exclusion of other centres or higher-degree systems.',
        'point_searches':0,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,np,bp]},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    d=verify();text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(d['status'],flush=True)

#!/usr/bin/env sage-python
"""Independent RR-net coefficient proof and retrospective frozen-chart bridge.

No point search. Generic construction is complete before this checker reads
the immutable search chart. Its returned coordinate is used only as a replay
witness, never to choose a net member or tune a parameter.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
CERT=ART/'det1092_first_centre_rr_net_v1.json';OUT=ART/'det1092_first_centre_rr_net_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for path,h in d['inputs'].items():assert sha(ROOT/path)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';parent=json.loads(pp.read_text())
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(v,ring=R):return ring.fraction_field()(ring(v['numerator']))/ring(v['denominator'])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);base=[E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    w=vector(ZZ,d['trace_word']);T=sum((n*P for n,P in zip(w,base)),E(0));G=matrix(QQ,parent['generic_height_gram']);assert w*G*w==10
    A=list(map(R,d['A']));B=list(map(R,d['B']))
    h=T[0].denominator().sqrt();terms=[R(h**3),R(T[0]*h**3),R(T[1]*h**3)]
    columns=[term*t**j for term,bound in zip(terms,[10,6,4]) for j in range(bound+1)]
    M=matrix(QQ,20,23,lambda i,j:columns[j][i]);assert M.rank()==20
    vectors=[]
    for row in [A,[t*f for f in A],B]:
        assert all(f.degree()<=bound for f,bound in zip(row,[10,6,4]))
        assert row[0]+row[1]*T[0]+row[2]*T[1]==0
        v=vector(QQ,[f[j] for f,bound in zip(row,[10,6,4]) for j in range(bound+1)])
        assert M*v==0;vectors.append(v)
    assert matrix(QQ,vectors).rank()==3
    evals=matrix(QQ,[[f(0) for f in row] for row in [A,[t*f for f in A],B]])
    assert evals.rank()==2 and evals[1]==0
    # Reconstruct elimination directly from the parent, not a constructor helper.
    def eliminate(curve,line,trace):
        field=curve.base_ring();X=PolynomialRing(field,'x');x=X.gen();f0,f1,f2=map(field,line)
        a1,a2,a3,a4,a6=curve.a_invariants();u=f0+f1*x
        pol=u*u-a1*x*u*f2-a3*u*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
        res,rem=pol.quo_rem(x-field(trace[0]));assert not rem and res.degree()==2
        return res
    res=eliminate(E,B,T);fixed=d['fixed_member'];assert res.list()==[dec(v) for v in fixed['residual_coefficients']]
    q=R(fixed['q_coefficients']);hs=dec(fixed['discriminant_square_factor']);assert res.discriminant()==hs*hs*q
    assert q.degree()==6 and q.gcd(q.derivative()).degree()==0
    assert R(E.discriminant()).gcd(q).degree()==0
    # Bridge the entire frozen chart, not just a fitted exceptional coordinate.
    ap=ART/'det1092_initial_unlock_construction_audit_v1.json';audit=json.loads(ap.read_text())['historical']
    assert list(-w)==audit['centre_word_in_generic17']
    Z=PolynomialRing(QQ,'z');F=Z.fraction_field();z=Z.gen()
    a,b,c,e=map(QQ,audit['map']['raw_slope_matrix']);ell=(a*z+b)/(c*z+e)
    # At zero the saved short-coordinate transport is X=36x+15,
    # Y=108(2y+x+1), so short slope = 6*literal slope+3.
    slope=(ell-3)/6
    u=-(B[1](0)+slope*B[2](0))/(A[1](0)+slope*A[2](0))
    line=[F(B[i](0))+u*A[i](0) for i in range(3)]
    E0=EllipticCurve(F,[v(0) for v in E.a_invariants()]);T0=E0([T[0](0),T[1](0)])
    res0=eliminate(E0,line,T0);cf=Z(audit['quartic_coefficients_low_to_high'])
    ratio=res0.discriminant()/cf
    assert ratio and ratio.is_square();scale=ratio.sqrt()
    n,den,root=map(ZZ,audit['homogeneous_square_witness']);z0=QQ(n)/den
    assert cf(z0)==(QQ(root)/den**2)**2
    # Construct the two fibre images from the replayed square, then certify one
    # independent gain without reading saved exceptional point coordinates.
    aa=res0[2](z0);bb=res0[1](z0);sq=scale(z0)*QQ(root)/den**2
    lp=[v(z0) for v in line]
    Es=EllipticCurve(QQ,[0,0,0,-27*E.c4()(0),-54*E.c6()(0)])
    def short(P):return Es([36*P[0]+15,108*(2*P[1]+P[0]+1)])
    Ezero=EllipticCurve(QQ,[v(0) for v in E.a_invariants()])
    x0=(-bb+sq)/(2*aa);y0=-(lp[0]+lp[1]*x0)/lp[2];P=Ezero([x0,y0])
    generic=[short(Ezero([p[0](0),p[1](0)])) for p in base]
    proof=runpy.run_path(str(CAS/'audit_det1092_initial_unlock.sage'))['finite_certificate'](Es,[*generic,short(P)])
    assert proof['rank']==18
    def rec(v):return {'numerator':list(map(str,v.numerator().list())),'denominator':list(map(str,v.denominator().list()))}
    return {'classification':'verified application and new deduction','status':'PASS_RR_NET_FULL_CHART_BRIDGE_AND_FIRST_GAIN',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,ap,CAS/'audit_det1092_initial_unlock.sage']},'checker_sha256':sha(Path(__file__)),
        'RR_rank':20,'kernel_dimension':3,'fibre_restriction_rank':2,
        'fixed_member_squarefree_degree':6,'fixed_member_geometric_genus':2,'branch_coprime_parent_discriminant':True,
        'chart_to_net_u':rec(u),'residual_discriminant_to_chart_square_scale':rec(scale),
        'witness_coordinate':str(z0),'reconstructed_point_literal302':list(map(str,P[:2])),
        'independence':proof,'point_searches':0,'target_selected_net_members':0,
        'boundary':'The whole chart is represented by a generic-input RR net. The successful coordinate is used only for retrospective verification. This does not predict square values, classify all covers, or explain why this fibre has a rank jump.'}
if __name__=='__main__':
    data=verify();text=json.dumps(data,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(data['status'],flush=True)

#!/usr/bin/env sage-python
"""Retrospective genus-drop audit of the first-witness RR pencil.

The witness selects a diagnostic pencil, not a new discovery construction.
One exact sextic, at most18 modular discriminants, no rational point search.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,lcm,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_first_witness_pencil_genus_gate_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json';bp=ART/'det1092_first_centre_rr_net_replay_v1.json';loader=CAS/'load_curve302_recovered_parent.sage'
    net=json.loads(np.read_text());bridge=json.loads(bp.read_text())
    E,base,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    R=E.base_ring().ring();C=-sum((n*P for n,P in zip(net['trace_word'],base)),E(0))
    ur=bridge['chart_to_net_u'];z0=QQ(bridge['witness_coordinate']);u0=R(ur['numerator'])(z0)/R(ur['denominator'])(z0)
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2;aa=-E.c4()/48
    h=R(C[0].denominator().sqrt());nx=R(cx*h*h);ny=R(cy*h**3)
    V=PolynomialRing(QQ,'v');v=V.gen();S=PolynomialRing(V,'t');t=S.gen()
    A=[S(R(row)) for row in net['A']];B=[S(R(row)) for row in net['B']]
    f1=B[1]+(u0+v*t)*A[1];f2=B[2]+(u0+v*t)*A[2]
    ell=-f1+S(R(E.a1()))*f2/2
    L,rem=f2.quo_rem(S(h));assert not rem
    raw=ell**4-6*S(nx)*ell**2*L**2-8*S(ny)*ell*L**3-(3*S(nx)**2+4*S(R(aa))*S(h)**4)*L**4
    q,rem=raw.quo_rem(S(h)**6);assert not rem and q.degree()==6
    q=S(q)
    coeff=[QQ(c) for p in q.list() for c in p.list()];den=lcm(c.denominator() for c in coeff);content=gcd(ZZ(c*den) for c in coeff)
    q=S(q*QQ(den)/content);scale=QQ(content)/den
    assert all(c.denominator()==1 for p in q.list() for c in p.list())
    assert max(p.degree() for p in q.list())==4
    assert all(q[i][j]==0 for i in range(7) for j in range(i+1,5))
    qa=R([q[i][4] for i in range(4,7)]);limit=R([q[j][j] for j in range(5)])
    assert qa.degree()==2 and qa(0) and qa.gcd(qa.derivative()).degree()==0
    assert limit.degree()==4 and limit.gcd(limit.derivative()).degree()==0
    rows=[];found=False
    for p in [17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]:
        vp=PolynomialRing(GF(p),'v');tp=PolynomialRing(vp.fraction_field(),'t')
        qp=tp([vp(c.list()) for c in q.list()]);dp=qp.discriminant()
        assert dp.denominator().degree()==0
        dp=vp(dp);roots=[int(a) for a,e in dp.roots()] if dp else []
        row={'prime':p,'discriminant_degree':int(dp.degree()),'roots':roots,'coefficients':list(map(int,dp.list()))}
        rows.append(row)
        if dp.degree()==28 and not roots:found=True;break
    return {'classification':'retrospective verified application','status':'PASS_NO_RATIONAL_GENUS_DROP_IN_FIRST_WITNESS_PENCIL' if found else 'BOUNDED_MODULAR_GATE_INCONCLUSIVE',
        'input_boundary':'The saved successful chart coordinate fixes this diagnostic pencil. It is not an oracle-free selection or construction.',
        'witness_coordinate':str(z0),'u0':str(u0),'q_coefficients_t_then_v':[[str(c) for c in p.list()] for p in q.list()],
        'branch_squareclass_scale':str(scale),'branch_identity':'f2^4*H(t,-f1/f2+a1/2)=h^6*scale*q(t,v)',
        'trace_denominator_root':list(map(str,h.list())),
        'infinity_degeneration':{'quadratic':list(map(str,qa.list())),'quartic':list(map(str,limit.list())),
            'discriminant_degree':28,'proof':'epsilon^4*q(t,1/epsilon) has four simple roots of order epsilon near zero and two simple separated roots; its discriminant has order12. Homogeneous sextic discriminant degree10 gives degree40-12=28.'},
        'modular_checks':rows,'no_rational_genus_drop':found,
        'limits':{'diagnostic_pencils':1,'modular_discriminants_max':18,'prime_bound':197,'global_discriminant_factorizations':0,'point_searches':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,np,bp,loader,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],[(r['prime'],r['discriminant_degree'],len(r['roots'])) for r in d['modular_checks']],flush=True)

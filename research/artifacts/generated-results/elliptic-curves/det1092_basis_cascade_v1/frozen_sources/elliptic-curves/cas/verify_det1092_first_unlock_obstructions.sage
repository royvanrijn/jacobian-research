#!/usr/bin/env sage-python
"""Independent coefficient replay; no RR constructor or rational point search.

Verified application. Check both fixed first-centre multisections, including
uniqueness of their RR lines and exact nonsplitting at zero. Completed search
evidence is consulted only after construction, to bind the centre labels.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,power_mod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_first_unlock_obstructions_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    pp=ART/'curve302_recovered_mw17_parent_v1.json'
    ap=ART/'det1092_initial_unlock_construction_audit_v1.json'
    parent=json.loads(pp.read_text());audit=json.loads(ap.read_text())
    R=PolynomialRing(QQ,'t');K=R.fraction_field();X=PolynomialRing(K,'x');x=X.gen()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    def old(v):return K(R(v['numerator_coefficients_low_to_high']))/R(v['denominator_coefficients_low_to_high'])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']])
    basis=[E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    G=matrix(QQ,parent['generic_height_gram']);rows=[];paths=[pp,ap]
    pairs=[('historical','det1092_historical_unlock_obstruction_v1.json'),
           ('autonomous_V3','det1092_initial_unlock_bisection_obstruction_v1.json')]
    for label,name in pairs:
        cp=ART/name;paths.append(cp);d=json.loads(cp.read_text())
        for path,digest in d['inputs'].items():assert sha(ROOT/path)==digest
        w=-vector(ZZ,d['trace_word']);assert list(w)==audit[label]['centre_word_in_generic17']
        assert w*G*w==10
        T=-sum((n*P for n,P in zip(w,basis)),E(0))
        f0,f1,f2=[R(d['line_coefficients'][k]) for k in ['f0','f1','f2']]
        assert [f0.degree(),f1.degree(),f2.degree()]==[9,5,3]
        assert f0.gcd(f1).gcd(f2).degree()==0
        assert f0+f1*T[0]+f2*T[1]==0
        # A separate rational-function matrix, clearing one common denominator,
        # checks that this is the unique RR member in the stated degree bounds.
        functions=[K(R.gen()**j)*v for v,bound in [(K(1),9),(T[0],5),(T[1],3)] for j in range(bound+1)]
        den=R(1)
        for v in functions:den=den.lcm(v.denominator())
        polys=[R(v*den) for v in functions];height=max(v.degree() for v in polys)
        A=matrix(QQ,height+1,20,lambda i,j:polys[j][i]);assert A.rank()==19
        residual=X([old(v) for v in d['residual_coefficients']]);c,b,a=residual.list()
        a1,a2,a3,a4,a6=E.a_invariants();u=f0+f1*x
        eliminated=u*u-a1*x*u*f2-a3*u*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
        assert eliminated==(x-T[0])*residual
        assert a(0) and f2(0) and E.discriminant()(0)
        disc=b*b-4*a*c;zero=QQ(disc(0));assert zero and not zero.is_square()
        # The general lift is valid by elimination; the residual roots are the
        # two non-trace intersections. No exceptional coordinates are needed.
        row={'case':label,'RR_rank':19,'zero_discriminant_nonsquare':True,
             'equation_and_trace_removal_verified':True,'no_rational_fibre_point_at_zero':True}
        if label=='historical':
            q=R(d['q_coefficients']);h=R(d['h_coefficients'])
            assert disc==q*h*h and h(0) and q.degree()==2
            assert q.gcd(q.derivative()).degree()==0
            z=QQ(q(0));assert z==QQ(d['zero_q'])
            p=ZZ(d['nonsquare_witness']['prime']);assert p.is_prime() and z.denominator()%p
            r=z.numerator()%p*(z.denominator()%p).inverse_mod(p)%p
            assert r==d['nonsquare_witness']['residue'] and power_mod(r,(p-1)//2,p)==p-1
            row.update({'q0_mod_prime':int(r),'prime':int(p),'geometric_genus':0})
        rows.append(row)
    return {'classification':'verified application','status':'PASS_INDEPENDENT_BOTH_FIRST_CENTRE_OBSTRUCTIONS',
        'cases':rows,'point_searches_run':0,'exceptional_point_construction_inputs':False,
        'claim_boundary':'The two specified genus-zero RR bisections cannot supply rational302 points at t=0. The actual arithmetic source of the first exceptional point and other multisections remain unresolved.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    d=verify();text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(d['status'],flush=True)

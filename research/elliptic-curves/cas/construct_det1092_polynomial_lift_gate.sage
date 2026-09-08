#!/usr/bin/env sage-python
"""Generic-only obstruction for all polynomial linear-norm lifts, degree<=24.

Two existing reduction primes, two exact one-parameter square identities,
nine unchanged addresses. No coefficient/parameter sampling or point search.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_polynomial_lift_gate_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True,default=int)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def degrees(parts):
    result={0}
    for n in parts:result|={k+n for k in result}
    return sorted(result)
def construct():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',
           ART/'curve302_parent_geometric_picard19_v1.json',
           ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent,geometry,roster=map(read,paths)
    protocol={'classification':'generic-only exact construction-family obstruction',
        'rule':'For alpha=f_t(k(t))*(k(t)-theta), classify polynomial k over Q of degree<=24 with no ramification at finite good parameter fibres. Freeze only the existing primes149 and151 and the two discriminant-root remainder families.',
        'limits':{'wall_seconds':25,'polynomial_degree_bound':24,'reduction_primes':[149,151],
            'boundary_templates':2,'existing_addresses':9,'sampled_coefficients':0,
            'sampled_parameters':0,'new_addresses':0,'point_searches':0,
            'number_fields':0,'class_groups':0,'Selmer_dimension_runs':0,
            'exceptional_points':0,'V3_inputs':0,'pilot_changes':0},
        'geometry_dependency':'EC-CURVE302-PARENT-GEOMETRIC-PICARD19: displayed17 rational sections are the full torsion-free geometric MW basis. Not recomputed here.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'protocol.json',protocol)
    assert geometry['status']=='PASS' and geometry['full_geometric_basis_is_displayed_rational_basis']
    assert geometry['geometric_generic_MW_rank']==17
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def rat(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(rat,parent['a_invariants']));assert ai[:3]==[1,1,1]
    a,b,c=R(5),R(16*ai[3]+8),R(64*ai[4]+16)
    discriminant=a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    Delta=discriminant.monic()
    assert Delta.degree()==24 and Delta.gcd(Delta.derivative())==1
    denominator=a*a-3*b
    r=((9*c-a*b)*(2*denominator).inverse_mod(Delta))%Delta
    s=(-a-2*r)%Delta
    assert r.degree()==s.degree()==23
    assert (r**3+a*r*r+b*r+c)%Delta==0
    assert (3*r*r+2*a*r+b)%Delta==0
    assert (s**3+a*s*s+b*s+c)%Delta==0
    assert (s-r).gcd(Delta)==1
    templates=[]
    for label,z,order in [('double',r,2),('simple',s,3)]:
        A=3*z+a;B=3*z*z+2*a*z+b;C=(z**3+a*z*z+b*z+c)//Delta
        if label=='simple':
            assert (B-A*A/4)%Delta==0
            h2=(B-A*A/4)//(2*Delta)
            assert h2.degree()<=22
        templates.append({'label':label,'root_remainder':list(map(str,z.list())),
            'A':list(map(str,A.list())),'B':list(map(str,B.list())),
            'C':list(map(str,C.list())),'known_residual_w_order':order})
    prime_rows=[];allowed=set(range(25))
    for p in protocol['limits']['reduction_primes']:
        F=GF(p);T=PolynomialRing(F,'t');dp=T(Delta)
        assert dp.degree()==24 and dp.is_squarefree()
        factors=dp.factor();assert all(e==1 for h,e in factors)
        parts=[int(h.degree()) for h,e in factors]
        possible=degrees(parts);allowed&=set(possible)
        row={'prime':p,'factor_degrees':parts,'possible_rational_factor_degrees':possible,
             'factors':[list(map(int,h.list())) for h,e in factors],'square_gates':[]}
        W=PolynomialRing(F,'w');w=W.gen();V=PolynomialRing(W,'t');t=V.gen()
        def cv(poly):return V(T(poly).list())
        for template in templates:
            A,B,C=[R(list(map(QQ,template[key]))) for key in ['A','B','C']]
            Q=cv(Delta)**2+w*cv(A*Delta)+w*w*cv(B)+w**3*cv(C)
            h=t**24
            for i in range(23,-1,-1):h+=(Q-h*h)[24+i]/F(2)*t**i
            residual=Q-h*h;order=template['known_residual_w_order']
            assert residual.degree()==23
            assert all(v%w**order==0 for v in residual.list())
            normalized=[residual[i]//w**order for i in range(24)]
            assert normalized[23].degree()==25-order
            bezout=[W(0)]*24;g=W(0)
            for i in range(23,-1,-1):
                gg,u,v=g.xgcd(normalized[i]);bezout=[u*b for b in bezout]
                bezout[i]+=v;g=gg
                if g==1:break
            assert sum(b*q for b,q in zip(bezout,normalized))==1
            row['square_gates'].append({'label':template['label'],
                'known_residual_w_order':order,
                'monic_square_root_truncation':[list(map(int,h[i].list())) for i in range(25)],
                'normalized_residuals':[list(map(int,q.list())) for q in normalized],
                'bezout_multipliers':[list(map(int,q.list())) for q in bezout],
                'degree_preserving_residual_t_power':23,
                'attained_w_degree_bound':25-order,'normalized_residual_gcd':1})
        prime_rows.append(row)
        retain(DIR/f'prime-{p}.json',row)
    assert allowed=={0,24}
    cases=[]
    for i,entry in enumerate(roster['cases']):
        tau=QQ(entry['parameter']);assert Delta(tau)
        cases.append({'index':i,'label':entry['label'],'parameter':entry['parameter'],
            'smooth_fibre':True,
            'uniform_conclusion':'Any specialization of a generically finite-good-unramified polynomial linear-norm lift of degree<=24 is inherited. This is not an exclusion of other seeds or isolated soluble fibres of ramified families.'})
    result={'status':'PASS_CANDIDATE_POLYNOMIAL_LIFT_OBSTRUCTION',
        'classification':'new deduction; independent replay required',
        'a':list(map(str,a.list())),'b':list(map(str,b.list())),'c':list(map(str,c.list())),
        'monic_discriminant':list(map(str,Delta.list())),
        'discriminant_factor_degree_intersection':sorted(allowed),
        'templates':templates,'primes':prime_rows,'cases':cases,
        'conclusion':'Every polynomial k in Q[t] of degree<=24 whose f(k)*(k-theta) class is unramified at every finite good parameter fibre has f(k) a rational polynomial square and is the Kummer class of a generic MW17 section.',
        'limitations':'Fixed original parameter and cubic model; not an obstruction for arbitrary rational abscissas, higher degrees, genuinely quadratic cubic-algebra representatives, ramified class families, or nonconstant base changes. No fibre rank upper bound.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'protocol.json',Path(__file__)]}}
    retain(DIR/'construction.json',result)
    print(result['status'],[r['factor_degrees'] for r in prime_rows],flush=True)
if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);construct()

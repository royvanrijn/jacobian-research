#!/usr/bin/env sage-python
"""One equation-only conic specialization, then a certified congruence family.

Fixed u=0; all odd primes<=1009; no parameter/point search or exceptional
input. The modulus is a proof output, not a score or tuned point selector.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,lcm,gcd,prime_range,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_conic_seed_progression_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True,default=int)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def record(f):return {'numerator':list(map(str,f.numerator().list())),
                      'denominator':list(map(str,f.denominator().list()))}
def residue(a,p):
    a=QQ(a);assert a.denominator()%p
    return int(a.numerator()*a.denominator().inverse_mod(p)%p)
def primitive_polynomials(functions,R):
    denominator=lcm([f.denominator() for f in functions])
    polys=[R(f*denominator) for f in functions]
    coefficient_den=lcm([a.denominator() for f in polys for a in f])
    polys=[f*coefficient_den for f in polys]
    common=gcd(polys);polys=[R(f/common) for f in polys]
    coefficient_den=lcm([a.denominator() for f in polys for a in f])
    polys=[f*coefficient_den for f in polys]
    content=gcd([ZZ(a) for f in polys for a in f])
    return [R(f/content) for f in polys]
def construct():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',
        ART/'det1092_orbit8044_rank18_base_change_v2.json',
        ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent,cover,roster=map(read,paths)
    protocol={'classification':'equation-only single-seed construction and specialization proof',
        'rule':'Reuse exactly the already constructed rational conic8044, evaluate its fixed parametrization at u=0, and retain all finite proof exposures at odd primes<=1009. Do not change u or the bound if proof fails. If rank18 is certified, extract independent proof rows and a congruence modulus preserving those reductions.',
        'limits':{'wall_seconds':25,'conics':1,'new_specialization_addresses':1,
            'conic_parameter':0,'finite_prime_bound':1009,'old_control_addresses':9,
            'parameter_sweeps':0,'point_searches':0,'exceptional_points':0,
            'V3_inputs':0,'pilot_changes':0,'class_groups':0,'Selmer_dimensions':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'protocol.json',protocol)
    for path,digest in cover['inputs'].items():assert sha(ROOT/path)==digest
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    U=PolynomialRing(QQ,'u');F=U.fraction_field()
    def decode(d,R,K):return K(R(d['numerator']))/R(d['denominator'])
    T=decode(cover['parametrization']['t_of_u'],U,F)
    S=decode(cover['parametrization']['s_of_u'],U,F)
    q=R(cover['curve_over_Q']['q_coefficients'])
    assert S*S==q(T) and T.denominator()(0)
    def pull(d):
        f=decode(d,R,K)
        return F(f.numerator()(T))/f.denominator()(T)
    ai=[pull(d) for d in parent['a_invariants']]
    generic=[[pull(d) for d in P] for P in parent['basis_weierstrass_coordinates']]
    new=[decode(cover['lift'][k],U,F) for k in ['x_of_u','y_of_u']]
    allpoints=generic+[new]
    projective=[primitive_polynomials([*P,F(1)],U) for P in allpoints]
    coefficient_pairs=[primitive_polynomials([a,F(1)],U) for a in ai]
    assert all(pair[1](0) for pair in coefficient_pairs)
    E=EllipticCurve(QQ,[a(0) for a in ai])
    points=[E([h(0) for h in P]) for P in projective]
    assert E.discriminant()
    input_data={'classification':'frozen equations before finite rank tests',
        'parameter':str(T(0)),'conic_parameter':0,
        'base_map':record(T),'conic_ordinate':record(S),
        'weierstrass_functions':[record(a) for a in ai],
        'point_functions':[[record(a) for a in P] for P in allpoints],
        'projective_point_polynomials':[[list(map(str,a.list())) for a in P] for P in projective],
        'coefficient_polynomial_pairs':[[list(map(str,a.list())) for a in P] for P in coefficient_pairs],
        'equation_at_zero':list(map(str,E.a_invariants())),
        'points_at_zero':[list(map(str,P)) for P in points],
        'generic_count':17,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'protocol.json']}}
    retain(DIR/'input.json',input_data)
    primes=list(prime_range(3,1010));rows=[];owners=[];exposures=[];no2=None
    b2,b4,b6=E.b2(),E.b4(),E.b6()
    for p0 in primes:
        p=int(p0)
        if any(a.denominator()%p==0 for a in E.a_invariants()):
            exposures.append({'p':p,'status':'SKIP_COEFFICIENT_DENOMINATOR'});continue
        if E.discriminant().valuation(p)!=0:
            exposures.append({'p':p,'status':'SKIP_BAD_REDUCTION'});continue
        # Yc=2*y+a1*x+a3 gives Yc^2=4*x^3+b2*x^2+2*b4*x+b6.
        # Use X=4*x and Y=4*Yc, a monic cubic with coefficients b2,8*b4,16*b6.
        aa,bb,cc=residue(b2,p),residue(8*b4,p),residue(16*b6,p)
        roots=[r for r in range(p) if (r*r*r+aa*r*r+bb*r+cc)%p==0]
        assert len(roots) in [0,1,3]
        if not roots and no2 is None:no2=p
        bits=[]
        for P in points:
            if P.is_zero():bits.append([0]*len(roots));continue
            Xp=4*P[0]
            if Xp.valuation(p)<0:bits.append([0]*len(roots));continue
            xx=residue(Xp,p);column=[]
            for r in roots:
                value=(xx-r)%p
                if not value:value=(3*r*r+2*aa*r+bb)%p
                assert value
                column.append(int(pow(value,(p-1)//2,p)==p-1))
            bits.append(column)
        block=matrix(GF(2),18,len(roots),[b for col in bits for b in col]).transpose()
        rr=[list(map(int,r)) for r in block.rows()]
        rows.extend(rr);owners.extend([p]*len(rr))
        exposures.append({'p':p,'status':'PASS_CUBIC_KUMMER_REDUCTION',
            'cubic':[cc,bb,aa,1],'rational_roots':roots,'rows':rr})
        retain(DIR/f'prime-{p}.json',exposures[-1])
    M=matrix(GF(2),len(rows),18,[v for row in rows for v in row])
    rank=int(M.rank());generic_rank=int(M[:,:17].rank())
    summary={'classification':'exact finite proof proposal; independent replay required',
        'rank':rank,'generic_rank':generic_rank,'no_rational_2_torsion_prime':no2,
        'exposures':exposures,'matrix_rows':rows,
        'status':'CANDIDATE_RANK18' if rank==18 and no2 else 'UNRESOLVED_FIXED_SPECIALIZATION',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'input.json',DIR/'protocol.json',Path(__file__)]}}
    retain(DIR/'finite-code.json',summary)
    print('CHECKPOINT_FIXED_U_ZERO',rank,generic_rank,no2,flush=True)
    if rank!=18 or no2 is None:return
    pivot_rows=list(M.transpose().pivots())
    selected=sorted(set([owners[i] for i in pivot_rows]+[no2]))
    # One p-adic residue ball preserves all coefficients and projective points.
    stability=[]
    Tpair=primitive_polynomials([T,F(1)],U)
    for p in selected:
        coefficient_orders=[int(pair[1](0).valuation(p)) for pair in coefficient_pairs]
        assert all(e>=0 for e in coefficient_orders)
        assert all(pair[0](0).valuation(p)>=e for pair,e in zip(coefficient_pairs,coefficient_orders))
        point_orders=[min(int(a(0).valuation(p)) for a in P if a(0)) for P in projective]
        t_order=int(Tpair[1](0).valuation(p))
        exponent=1+max(coefficient_orders+point_orders+[t_order])
        stability.append({'p':p,'exponent':exponent,
            'coefficient_denominator_orders':coefficient_orders,
            'point_projective_common_orders':point_orders,'base_map_denominator_order':t_order})
    modulus=prod(ZZ(r['p'])**r['exponent'] for r in stability)
    controls=[]
    for i,case in enumerate(roster['cases']):
        value=q(QQ(case['parameter']))
        numerator,denominator=value.numerator(),value.denominator()
        split=bool(value>=0 and numerator.is_square() and denominator.is_square())
        row={'index':i,'label':case['label'],'parameter':case['parameter'],'q_value':str(value),'split':split}
        if value>=0:
            row['numerator_floor_sqrt']=str(numerator.isqrt())
            row['denominator_floor_sqrt']=str(denominator.isqrt())
        controls.append(row)
    result={'status':'CANDIDATE_ORACLE_FREE_RANK18_ARITHMETIC_PROGRESSION',
        'classification':'new constructive deduction; independent replay required',
        'base_parameter_at_u_zero':str(T(0)),
        'progression':'u=B*n, n any rational integer; t=t(u), Q=the fixed conic lift.',
        'B':str(modulus),'B_bits':int(modulus.nbits()),
        'selected_proof_primes':selected,'pivot_row_indices':pivot_rows,
        'stability':stability,'base_map_polynomial_pair':[list(map(str,a.list())) for a in Tpair],
        'controls':controls,
        'boundary':'This is one extra direction and an infinite rational rank>=18 specialization family from the already fixed conic. The conic misses302 and does not explain its first seed. No amplifier or point search is run.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [DIR/'input.json',DIR/'finite-code.json',DIR/'protocol.json',Path(__file__)]}}
    retain(DIR/'progression.json',result)
    print(result['status'],'B_bits',modulus.nbits(),'primes',selected,flush=True)
if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);construct()

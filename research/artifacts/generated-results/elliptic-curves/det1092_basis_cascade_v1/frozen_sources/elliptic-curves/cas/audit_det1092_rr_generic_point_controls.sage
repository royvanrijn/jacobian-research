#!/usr/bin/env sage-python
"""Outcome-blind specificity controls using only a known generic section.

The RR centre is the historically calibrated generic word, unchanged. The
member selector and marked point use basis section0 only, never an exceptional
point or successful chart coordinate. Freeze nine cases; execute one per call.
"""
import argparse,hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,prime_range,lcm,gcd,power_mod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_rr_generic_point_controls_v1'
LOCAL=ROOT/'artifacts/local/elliptic-curves/det1092-rr-generic-point-controls-v1'
PROTOCOL=DIR/'protocol.json'
ROSTER=ROOT/'artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap-v2/roster.json'
TRANSPORT=ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def source():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json'
    up=ART/'det1092_universal_rr_descent_preflight_v1.json'
    parent,net,uni=[json.loads(p.read_text()) for p in [pp,np,up]]
    uni=uni['universal_genus2']
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field()
    def rat(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[rat(row) for row in parent['a_invariants']])
    basis=[E([rat(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    w=-vector(QQ,net['trace_word']);G=matrix(QQ,parent['generic_height_gram'])
    C=sum((n*P for n,P in zip(w,basis)),E(0));h=R(uni['h'])
    assert h*h==C[0].denominator()
    functions=[]
    for P in basis:
        functions.append(-(B[0]+B[1]*P[0]+B[2]*P[1])/(A[0]+A[1]*P[0]+A[2]*P[1]))
    return R,E,basis,A,B,w,G,C,h,uni,functions,[pp,np,up]
def freeze():
    R,E,basis,A,B,w,G,C,h,uni,functions,paths=source()
    roster=json.loads(ROSTER.read_text())
    transport=json.loads(TRANSPORT.read_text());a,b,c,d=map(QQ,transport['parameter_matrix'])
    cases=[]
    for row in roster['cases']:
        s=QQ(row['parameter']);t=(a*s+b)/(c*s+d)
        cases.append({'label':row['id'],'parameter':str(t),'reduced_parameter':str(s),
                      'elliptic_model':row['model']})
    cases.append({'label':'302-generic-section-control','parameter':'0',
                  'reduced_parameter':None,'elliptic_model':list(map(str,[v(0) for v in E.a_invariants()]))})
    assert len(cases)==9
    selector=functions[0];assert max(selector.numerator().degree(),selector.denominator().degree())==7
    protocol={'classification':'frozen generic-input specificity-control experiment',
              'selector':{'generic_basis_index':0,'u_numerator':list(map(str,selector.numerator().list())),
                          'u_denominator':list(map(str,selector.denominator().list())),'v':'0',
                          'rule':'u(t)=-B(S0(t))/A(S0(t)), v=0; mark the known point S0(t), with curve coordinate T=t.'},
              'universal_curve':'s^2=c*q(T;u(t),0), using the retained generic universal sextic',
              'basepoint':'Intersection with O; no exceptional point is used to choose its sign.',
              'cases':cases,'primes':list(map(int,list(prime_range(17,500))[:64])),
              'comparison_boundary':'Every marked elliptic point is the known generic basis section0. A non-generic Jacobian class would be a specificity counterexample, not an elliptic rank gain.',
              'limits':{'wall_seconds_per_case':25,'cases':9,'primes_per_case':64,'generic_section_divisors':17,
                        'point_searches':0,'global_Selmer_runs':0,'class_group_runs':0,'pilot_changes':0},
              'selection_inputs':'Generic equation, generic section0, calibrated generic RR net, frozen V4 roster parameters. No exceptional points, winning coordinates, scores or null outcomes enter selection.',
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[ROSTER,TRANSPORT,Path(__file__)]}}
    retain(PROTOCOL,protocol)
    print('PASS_FROZEN_GENERIC_POINT_CONTROLS',len(cases),flush=True)
def run(index):
    protocol=json.loads(PROTOCOL.read_text())
    for path,digest in protocol['inputs'].items():assert sha(ROOT/path)==digest
    case=protocol['cases'][index];tau=QQ(case['parameter'])
    R,E,basis,A,B,w,G,C,h,uni,functions,paths=source();T=R.gen()
    u=functions[0](tau)
    q=R([sum(QQ(row['coefficient'])*u**row['u'] for row in uni['sparse_q']
             if row['T']==i and row['v']==0) for i in range(7)])
    denominator=lcm(a.denominator() for a in q);content=gcd(ZZ(a*denominator) for a in q)
    q=R(q*denominator/content);scale=QQ(uni['scale'])*content/denominator
    assert q.degree()==6 and q.gcd(q.derivative())==1
    alpha=(A[2]/h)[0];b=R(B[2]/h);assert b.degree()==1 and alpha
    x0=-(b[0]+alpha*u)/b[1]
    assert scale*q(x0) and (scale*q(x0)).is_square()
    f0,f1,f2=[B[j]+u*A[j] for j in range(3)]
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    mnum=-f1+E.a1()*f2/2
    E0=EllipticCurve(QQ,[a(tau) for a in E.a_invariants()])
    assert E0.discriminant() and E0.j_invariant()==EllipticCurve(QQ,case['elliptic_model']).j_invariant()
    marked=E0([basis[0][0](tau),basis[0][1](tau)])
    polynomials=[];degrees=[];mod_s_hashes=[]
    for j,P in enumerate(basis):
        g=(functions[j].numerator()-u*functions[j].denominator()).monic()
        degree=ZZ(1+G[j,j]-(w*G)[j])
        assert g.degree()==degree and g.gcd(q)==g.gcd(g.derivative())==1
        sy=((2*(P[0]+E.b2()/12)+cx)*f2**2-mnum**2)/h**3
        assert sy.denominator().gcd(g)==1
        s=(sy.numerator()*sy.denominator().inverse_mod(g))%g
        assert (s*s-scale*q)%g==0
        polynomials.append(g);degrees.append(int(degree))
        mod_s_hashes.append(hashlib.sha256(json.dumps(list(map(str,s.list()))).encode()).hexdigest())
        if j==0:ystar=sy(tau)
    assert polynomials[0](tau)==0 and ystar*ystar==scale*q(tau) and ystar
    m=mnum(tau)/f2(tau);omega=h(tau)**3*ystar/f2(tau)**2
    xx=(m*m-cx(tau)+omega)/2;yy=m*(xx-cx(tau))-cy(tau)
    literal_x=xx-E.b2()(tau)/12;literal_y=yy-(E.a1()(tau)*literal_x+E.a3()(tau))/2
    assert E0([literal_x,literal_y])==marked
    rows=[];trials=[];six=None;five=None
    checkpoint=LOCAL/case['label'];checkpoint.mkdir(parents=True,exist_ok=True)
    for p in protocol['primes']:
        trial={'p':p};coefficients=list(q)+[x0,tau]+[a for g in polynomials for a in g]
        if any(a.denominator()%p==0 for a in coefficients):trial['status']='SKIP_DENOMINATOR'
        else:
            F=GF(p);S=PolynomialRing(F,'X');X=S.gen();f=S(q.list());gs=[S(g.list()) for g in polynomials]
            if f.degree()!=6 or f.gcd(f.derivative())!=1:trial['status']='SKIP_BAD_SEXTIC'
            elif f(F(x0))==0 or f(F(tau))==0 or any(g.gcd(f)!=1 for g in gs):trial['status']='SKIP_NONUNIT_DIVISOR'
            else:
                assert scale.valuation(p)%2==0
                factors=sorted([g.monic() for g,e in f.factor()],key=lambda a:(a.degree(),list(map(int,a))))
                bits=[]
                for a in factors:
                    exponent=(ZZ(p)**a.degree()-1)//2
                    anchor=power_mod(S(F(x0))-X,exponent,a)
                    values=[((-1)**g.degree()*g)%a for g in gs]+[S(F(tau))-X]
                    row=[]
                    for value,degree in zip(values,degrees+[1]):
                        sign=(power_mod(value,exponent,a)*anchor**degree)%a
                        assert sign in [S(1),S(-1)]
                        row.append(int(sign==S(-1)))
                    bits.append(row)
                raw=matrix(GF(2),bits)
                assert all(sum(raw[:,j].list(),GF(2)(0))==0 for j in range(18))
                quotient=matrix(GF(2),1,len(factors),[a.degree()%2 for a in factors]).right_kernel().basis_matrix()
                block=quotient*raw;rows.extend([list(map(int,z)) for z in block.rows()])
                pattern=[int(a.degree()) for a in factors]
                if pattern==[6] and six is None:six=p
                if pattern==[1,5] and five is None:five=p
                trial.update({'status':'PASS_LOCAL_KUMMER_BLOCK','degrees':pattern,
                              'factors':[list(map(int,a.list())) for a in factors],
                              'raw_rows':[list(map(int,z)) for z in raw.rows()],
                              'quotient_rows':[list(map(int,z)) for z in quotient.rows()],
                              'block_rows':[list(map(int,z)) for z in block.rows()]})
        trials.append(trial);retain(checkpoint/f'prime-{p:04d}.json',trial)
    M=matrix(GF(2),rows);H=M[:,:-1]
    separator=next((z for z in H.left_kernel().basis() if (z*M[:,-1])[0]),None)
    complete=(H.rank()==16 and M.rank()==17 and six is not None and five is not None)
    result={'classification':'generic-input specificity-control calculation; independent replay required',
            'status':'PASS_CANDIDATE_GENERIC_ELLIPTIC_POINT_NON_GENERIC_JACOBIAN_CLASS' if complete else 'INCONCLUSIVE_FROZEN_GENERIC_POINT_CONTROL',
            'case':case,'u':str(u),'v':'0','q':list(map(str,q.list())),'scale':str(scale),
            'base_x':str(x0),'marked_curve_point':[str(tau),str(ystar)],
            'marked_elliptic_point':list(map(str,marked[:2])),'elliptic_basis_word':[1]+[0]*16,
            'divisor_degrees':degrees,'divisor_s_mod_g_sha256':mod_s_hashes,
            'trials':trials,'matrix_rows':rows,'inherited_character_rank':int(H.rank()),
            'with_marked_character_rank':int(M.rank()),
            'separator':None if separator is None else list(map(int,separator)),
            'theta_primes':{'six_cycle':six,'one_plus_five':five},
            'elliptic_quotient_class':'ZERO: exactly the first inherited generic section, not an extra elliptic direction',
            'rank_conclusion':'PENDING_INDEPENDENT_REPLAY',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[PROTOCOL,Path(__file__)]}}
    retain(DIR/f'case-{index:02d}.json',result)
    print(case['label'],result['status'],'ranks',H.rank(),M.rank(),'theta',six,five,flush=True)
if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);LOCAL.mkdir(parents=True,exist_ok=True)
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','case']);parser.add_argument('--index',type=int,default=0)
    args=parser.parse_args()
    if args.mode=='freeze':freeze()
    else:run(args.index)

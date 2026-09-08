#!/usr/bin/env sage-python
"""Exact marked-point Kummer transport and symmetric-pair obstruction.

Ten already certified points, one symbolic identity, bounded small-prime
cubic-field witnesses. No new point, Selmer, or class-group search.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_marked_kummer_transport_v1.json'
CONTROL=ART/'det1092_rr_generic_point_controls_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def symbolic():
    V=PolynomialRing(QQ,['a2','a4','cx','cy','m']);a2,a4,cx,cy,m=V.gens()
    S=PolynomialRing(V,'X');X=S.gen()
    a6=cy*cy-cx**3-a2*cx*cx-a4*cx
    cubic=X**3+a2*X**2+a4*X+a6;line=m*(X-cx)-cy
    pair,rem=(cubic-line**2).quo_rem(X-cx);pair=S(pair)
    assert not rem and pair.degree()==2
    assert ((cx-X)*pair-line**2)%cubic==0
    return {'classification':'new deduction verified symbolically',
            'cubic':'X^3+a2*X^2+a4*X+cy^2-cx^3-a2*cx^2-a4*cx',
            'line':'m*(X-cx)-cy',
            'pair_polynomial':'X^2+(a2-m^2+cx)*X+a4+cx^2+a2*cx+m^2*cx+2*m*cy',
            'identity':'(cx-theta)*pair(theta)=line(theta)^2',
            'consequence':'The unordered pair has Kummer class delta(C), hence zero modulo the inherited subgroup.'}
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json'
    parent,net=[json.loads(p.read_text()) for p in [pp,np]]
    R=PolynomialRing(QQ,'T');K=R.fraction_field()
    def rat(v):return K(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(K,[rat(v) for v in parent['a_invariants']])
    basis=[E([rat(v) for v in P]) for P in parent['basis_weierstrass_coordinates']]
    w=-vector(QQ,net['trace_word']);C=sum((n*P for n,P in zip(w,basis)),E(0))
    A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    frp=ART/'det1092_first_centre_rr_net_replay_v1.json';first=json.loads(frp.read_text())
    fjp=ART/'det1092_rr_full_inherited_jacobian_v1.json';fd=json.loads(fjp.read_text())
    fcp=ART/'det1092_rr_residual_jacobian_class_v1.json';fc=json.loads(fcp.read_text())
    cases=[{'label':'302-first-unlock','tau':'0','u':fc['curve']['u'],
            'point':first['reconstructed_point_literal302'],'q':fc['curve']['q'],
            'theta_primes':{'six_cycle':83,'one_plus_five':179},
            'elliptic_relative_class':'NONZERO','strict_mod_generic':'NO: certified first class is outside M17+K_strict'}]
    paths=[pp,np,frp,fjp,fcp]
    for i in range(9):
        cp=CONTROL/f'case-{i:02d}.json';rp=CONTROL/f'case-{i:02d}-replay.json'
        c=json.loads(cp.read_text());v=json.loads(rp.read_text())
        assert v['status']=='PASS_GENERIC_ELLIPTIC_POINT_GIVES_NON_GENERIC_RATIONAL_JACOBIAN_CLASS'
        cases.append({'label':c['case']['label'],'tau':c['case']['parameter'],'u':c['u'],
                      'point':c['marked_elliptic_point'],'q':c['q'],'theta_primes':c['theta_primes'],
                      'elliptic_relative_class':'ZERO','strict_mod_generic':'ZERO relative class, already in MW17'})
        paths.extend([cp,rp])
    protocol=json.loads((CONTROL/'protocol.json').read_text());primes=protocol['primes']
    S=PolynomialRing(QQ,'X');X=S.gen();rows=[]
    for case in cases:
        tau,u=QQ(case['tau']),QQ(case['u'])
        E0=EllipticCurve(QQ,[a(tau) for a in E.a_invariants()])
        centre=E0([C[0](tau),C[1](tau)]);P=E0(case['point']);Q=centre-P
        if case['elliptic_relative_class']=='ZERO':assert P==E0([basis[0][0](tau),basis[0][1](tau)])
        def transform(P):return 4*P[0],8*P[1]+4*P[0]+4
        cx,cy=transform(centre);px,py=transform(P);qx,qy=transform(Q)
        cubic=X**3+5*X**2+(16*E0.a4()+8)*X+64*E0.a6()+16
        assert cubic(cx)==cy*cy and cubic(px)==py*py and cubic(qx)==qy*qy
        assert cy and py and qy and cubic.gcd(cubic.derivative())==1
        f0,f1,f2=[B[j](tau)+u*A[j](tau) for j in range(3)]
        slope=1-2*f1/f2;line=slope*X+4-8*f0/f2
        assert line(cx)==-cy and line(px)==py and line(qx)==qy
        pair,rem=(cubic-line**2).quo_rem(X-cx);assert not rem and pair==(X-px)*(X-qx)
        beta,beta_other,beta_centre=S([px,-1]),S([qx,-1]),S([cx,-1])
        assert ((beta*beta_other)*beta_centre-line*line)%cubic==0
        assert cubic.resultant(beta)==py*py and cubic.resultant(beta_centre)==cy*cy
        tests=[];irreducible=None
        for p in primes:
            if any(a.denominator()%p==0 for a in cubic):
                tests.append({'p':p,'status':'SKIP_DENOMINATOR'});continue
            F=PolynomialRing(GF(p),'x');fp=F(cubic.list())
            if fp.gcd(fp.derivative())!=1:
                tests.append({'p':p,'status':'SKIP_BAD_CUBIC'});continue
            degrees=sorted(int(g.degree()) for g,e in fp.factor())
            tests.append({'p':p,'status':'TESTED','factor_degrees':degrees})
            if degrees==[3]:irreducible=p;break
        assert irreducible is not None
        sextic=R(case['q']);p6=case['theta_primes']['six_cycle'];p5=case['theta_primes']['one_plus_five']
        patterns=[]
        for p,expected in [(p6,[6]),(p5,[1,5])]:
            F=PolynomialRing(GF(p),'x');f=F(sextic.list())
            assert f.degree()==6 and f.gcd(f.derivative())==1
            pattern=sorted(int(g.degree()) for g,e in f.factor());assert pattern==expected
            patterns.append({'p':p,'factor_degrees':pattern})
        rows.append({**case,'literal_complement':list(map(str,Q[:2])),
                     'elliptic_cubic':list(map(str,cubic.list())),
                     'beta_branch':list(map(str,beta.list())),
                     'beta_complement':list(map(str,beta_other.list())),
                     'beta_centre':list(map(str,beta_centre.list())),
                     'pair_polynomial':list(map(str,pair.list())),
                     'square_line':list(map(str,line.list())),
                     'branch_norm_square_root':str(py),'centre_norm_square_root':str(cy),
                     'cubic_field_tests':tests,'cubic_irreducibility_prime':irreducible,
                     'sextic_field_patterns':patterns,
                     'sextic_has_no_proper_Q_subfields':True,
                     'direct_sextic_to_cubic_field_norm':'NOT_DEFINED: no cubic subfield',
                     'symmetric_pair_relative_class':'ZERO',
                     'branch_choice_relative_class_is_sign_independent':True,
                     'Sha_image_of_rational_branch':'ZERO'})
    codes=[row for record in first['independence']['signatures'] for row in record['rows']]
    M=matrix(GF(2),codes);H=M[:,:-1]
    separator=next(z for z in H.left_kernel().basis() if (z*M[:,-1])[0])
    assert H.rank()==17 and M.rank()==18
    return {'classification':'verified application and new deduction',
            'status':'PASS_MARKED_CUBIC_KUMMER_TRANSPORT_AND_NORM_BRIDGE_OBSTRUCTIONS',
            'symbolic':symbolic(),'cases':rows,
            'first_relative_nonzero_certificate':{'matrix_rows':codes,'separator':list(map(int,separator)),
                                                  'generic_rank':17,'with_point_rank':18},
            'primitive_field_argument':'A6-cycle proves transitivity. A5-cycle cannot preserve a nontrivial block system of2 or3 blocks of size3 or2, so the sextic field has no proper subfields.',
            'tensor_norm_boundary':'For beta in the sextic algebra, Norm_{A tensor K/K}(beta)=Norm_{A/Q}(beta). A norm-square Jacobian representative therefore gives a square scalar, not a new cubic Kummer class, by this tensor-norm route.',
            'scope':'Marked-point transport is conditional on an actual split RR fibre. It is not a homomorphism from the genus2 Jacobian or an unevaluated Selmer class to the fixed elliptic fibre.',
            'limits':{'wall_seconds':25,'marked_cases':10,'symbolic_identities':1,'cubic_prime_cap_per_case':64,
                      'point_searches':0,'global_Selmer_runs':0,'class_group_runs':0,'pilot_changes':0},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[CONTROL/'protocol.json',Path(__file__)]}}
if __name__=='__main__':
    signal.alarm(25);data=build();payload=json.dumps(data,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(data['status'],flush=True)

#!/usr/bin/env sage-python
"""Independent exact replay of marked transport and two norm obstructions.

Consumes the completed first-gain independence/strict certificates. No new
point search, class group, Selmer computation, or control selection.
"""
import hashlib,json,signal
from itertools import combinations
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,VectorSpace
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_marked_kummer_transport_v1.json'
OUT=ART/'det1092_marked_kummer_transport_replay_v1.json'
CONTROL=ART/'det1092_rr_generic_point_controls_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def partitions(xs,size):
    if not xs:
        yield frozenset();return
    a=min(xs)
    for rest in combinations(sorted(xs-{a}),size-1):
        block=frozenset((a,)+rest)
        for tail in partitions(xs-block,size):yield tail|{block}
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    dependencies=[]
    def read(name):
        p=ART/name;dependencies.append(p);return json.loads(p.read_text())
    parent=read('curve302_recovered_mw17_parent_v1.json')
    net=read('det1092_first_centre_rr_net_v1.json')
    first=read('det1092_first_centre_rr_net_replay_v1.json')
    firstcurve=read('det1092_rr_residual_jacobian_class_v1.json')
    uni=read('det1092_universal_rr_descent_preflight_v1.json')
    fil=read('curve302_recovered_quotient_local_filtration_v1.json')
    audit=read('det1092_initial_unlock_construction_audit_v1.json')['historical']
    protocol=read('det1092_rr_generic_point_controls_v2/protocol.json')
    assert first['status']=='PASS_RR_NET_FULL_CHART_BRIDGE_AND_FIRST_GAIN'
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def rat(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[rat(row) for row in parent['a_invariants']])
    assert (E.a1(),E.a2(),E.a3())==(1,1,1)
    base=[E([rat(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    centre=-sum((ZZ(n)*P for n,P in zip(net['trace_word'],base)),E(0))
    A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    # Direct expansion, independently of the constructor's polynomial division.
    V=PolynomialRing(QQ,['a2','a4','cx','cy','m','x'])
    a2,a4,cx,cy,m,x=V.gens()
    F=x**3+a2*x*x+a4*x+cy*cy-cx**3-a2*cx*cx-a4*cx
    L=m*(x-cx)-cy
    G=x*x+(a2-m*m+cx)*x+a4+cx*cx+a2*cx+m*m*cx+2*m*cy
    assert (cx-x)*G-L*L==-F
    # All nontrivial block systems for a transitive degree-six action.
    perm={0:1,1:2,2:3,3:4,4:0,5:5};blocks=[]
    for size,count in [(2,15),(3,10)]:
        systems=list(partitions(frozenset(range(6)),size));assert len(systems)==count
        invariant=[b for b in systems if frozenset(frozenset(perm[j] for j in s) for s in b)==b]
        assert not invariant
        blocks.append({'block_size':size,'systems_checked':count,'fixed_by_five_cycle':0})
    S=PolynomialRing(QQ,'X');X=S.gen();rows=[]
    assert len(d['cases'])==10 and len(protocol['primes'])==64
    for i,row in enumerate(d['cases']):
        if i==0:
            assert row['label']=='302-first-unlock' and QQ(row['tau'])==0
            assert row['u']==firstcurve['curve']['u'] and row['q']==firstcurve['curve']['q']
            assert row['point']==first['reconstructed_point_literal302']
            assert row['elliptic_relative_class']=='NONZERO'
        else:
            control=read('det1092_rr_generic_point_controls_v2/case-%02d.json'%(i-1))
            assert row['label']==control['case']['label'] and row['tau']==control['case']['parameter']
            assert row['u']==control['u'] and row['q']==control['q']
            assert row['point']==control['marked_elliptic_point'] and row['theta_primes']==control['theta_primes']
            assert row['elliptic_relative_class']=='ZERO'
        t,u=QQ(row['tau']),QQ(row['u'])
        Et=EllipticCurve(QQ,[a(t) for a in E.a_invariants()])
        C=Et([centre[0](t),centre[1](t)]);P=Et(row['point']);Q=Et(row['literal_complement'])
        assert P+Q==C and not any(Z.is_zero() for Z in [P,Q,C])
        if i:
            assert P==Et([base[0][0](t),base[0][1](t)])
        else:
            oldx=(QQ(audit['point_short_model'][0])-15)/36
            oldy=(QQ(audit['point_short_model'][1])/108-oldx-1)/2
            old=Et([oldx,oldy]);assert P in [old,C-old]
        points=[(4*Z[0],8*Z[1]+4*Z[0]+4) for Z in [P,Q,C]]
        F=S(row['elliptic_cubic'])
        assert F==X**3+5*X**2+(16*Et.a4()+8)*X+64*Et.a6()+16
        assert F.gcd(F.derivative())==1
        for xx,yy in points:assert yy and F(xx)==yy*yy
        # A multiplication matrix computes the norm independently of resultants.
        alg=S.quotient(F,'theta');theta=alg.gen();bb=[]
        for key,(xx,yy) in zip(['beta_branch','beta_complement','beta_centre'],points):
            b=S(row[key]);assert b==xx-X;bb.append(alg(b))
            multiplication=matrix(QQ,3,3,lambda a,j:(bb[-1]*theta**j).lift()[a])
            assert multiplication.det()==yy*yy
        assert QQ(row['branch_norm_square_root'])==points[0][1]
        assert QQ(row['centre_norm_square_root'])==points[2][1]
        f=[B[j](t)+u*A[j](t) for j in range(3)];assert f[2]
        line=S(row['square_line']);pair=S(row['pair_polynomial'])
        assert line==(1-2*f[1]/f[2])*X+4-8*f[0]/f[2]
        assert all(line(xx)==yy for xx,yy in points[:2]) and line(points[2][0])==-points[2][1]
        assert pair==(X-points[0][0])*(X-points[1][0])
        assert F-line**2==(X-points[2][0])*pair
        assert bb[0]*bb[1]*bb[2]==alg(line)**2
        # Retain every prime tried before the first irreducible cubic.
        trials=row['cubic_field_tests'];assert 1<=len(trials)<=64
        assert [z['p'] for z in trials]==protocol['primes'][:len(trials)]
        for j,trial in enumerate(trials):
            p=trial['p']
            if any(a.denominator()%p==0 for a in F):
                assert trial['status']=='SKIP_DENOMINATOR';continue
            Fp=PolynomialRing(GF(p),'z')(F.list())
            if Fp.gcd(Fp.derivative())!=1:
                assert trial['status']=='SKIP_BAD_CUBIC';continue
            ds=sorted(int(g.degree()) for g,e in Fp.factor())
            assert trial=={'p':p,'status':'TESTED','factor_degrees':ds}
            assert (ds==[3])==(j==len(trials)-1)
        assert trials[-1]['factor_degrees']==[3] and trials[-1]['p']==row['cubic_irreducibility_prime']
        q=S(row['q']);assert q.degree()==6
        assert len(row['sextic_field_patterns'])==2
        for gate,expected,key in zip(row['sextic_field_patterns'],[[6],[1,5]],['six_cycle','one_plus_five']):
            p=gate['p'];assert p==row['theta_primes'][key] and p in protocol['primes']
            qp=PolynomialRing(GF(p),'z')(q.list())
            assert qp.degree()==6 and qp.gcd(qp.derivative())==1
            ds=sorted(int(g.degree()) for g,e in qp.factor());assert ds==expected==gate['factor_degrees']
        assert row['sextic_has_no_proper_Q_subfields'] and row['symmetric_pair_relative_class']=='ZERO'
        assert row['branch_choice_relative_class_is_sign_independent'] and row['Sha_image_of_rational_branch']=='ZERO'
        rows.append({'label':row['label'],'elliptic_relative_class':row['elliptic_relative_class'],
                     'cubic_prime':row['cubic_irreducibility_prime'],'cubic_trials':len(trials),
                     'sextic_patterns':row['sextic_field_patterns'],'symmetric_relative_class':'ZERO'})
    # Consume the independently certified finite-group maps, check their linear
    # consequence here. This is not a new point-cloud independence calculation.
    codes=[r for z in first['independence']['signatures'] for r in z['rows']]
    proof=d['first_relative_nonzero_certificate'];assert codes==proof['matrix_rows']
    M=matrix(GF(2),codes);H=M[:,:-1];z=vector(GF(2),proof['separator'])
    assert H.rank()==proof['generic_rank']==17 and M.rank()==proof['with_point_rank']==18
    assert z*H==0 and (z*M[:,-1])[0]==1
    diagnostic=uni['first_witness_diagnostic'];ambient=VectorSpace(GF(2),31)
    generic=ambient.subspace(matrix(GF(2),parent['basis_embedding_in_public_D']).columns())
    strict=ambient.subspace(fil['local_filtration_mod_2']['strict_kernel_public_words'])
    phi=vector(GF(2),diagnostic['separating_character']);w=vector(GF(2),diagnostic['public_word'])
    assert (generic.dimension(),strict.dimension(),(generic+strict).dimension())==(17,10,27)
    assert (generic+strict).basis_matrix()*phi==0 and phi*w==1
    # Bind the public-word/point equality to its completed independent replay.
    strict_replay=read('det1092_universal_rr_descent_preflight_replay_v1.json')
    assert strict_replay['status']=='PASS_INDEPENDENT_UNIVERSAL_RR_DESCENT_PREFLIGHT'
    for p,h in strict_replay['inputs'].items():assert sha(ROOT/p)==h
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_MARKED_KUMMER_TRANSPORT_AND_NORM_OBSTRUCTIONS',
            'symbolic_pair_identity':True,'block_system_test':blocks,'cases':rows,
            'first_relative_elliptic_class':'NONZERO; finite character separates it from MW17',
            'first_strict_mod_generic':False,'nine_control_relative_elliptic_classes':'ZERO: exactly generic section0',
            'direct_norm_obstruction':'Every sextic field is primitive and every elliptic cubic is a field; no cubic subfield, hence no direct relative field norm.',
            'tensor_norm_obstruction':'Base extension of the multiplication matrix leaves its rational determinant unchanged. Thus Norm_(A tensor K)/K(beta)=Norm_A/Q(beta); norm-square inputs give trivial cubic squareclasses.',
            'scope':'No Jacobian-to-fixed-elliptic-fibre homomorphism or unevaluated Selmer transport constructed. Other auxiliary algebras or correspondences are not excluded.',
            'limits':d['limits'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT]+dependencies},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(result['status'],flush=True)

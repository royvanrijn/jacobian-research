#!/usr/bin/env sage-python
"""Independent exact-divisor, finite-field and theta-characteristic replay.

Uses extension-field norms instead of the constructor's polynomial powering.
No constructor is imported. All64 retained prime outcomes are checked.
"""
import hashlib,json,signal,itertools
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_rr_full_inherited_jacobian_v1.json'
OUT=ART/'det1092_rr_full_inherited_jacobian_replay_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for path,h in d['inputs'].items():assert sha(ROOT/path)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json'
    dp=ART/'det1092_rr_residual_jacobian_class_v1.json';rp=ART/'det1092_first_centre_rr_net_replay_v1.json'
    parent,net,diag,first=[json.loads(p.read_text()) for p in [pp,np,dp,rp]]
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def rat(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[rat(row) for row in parent['a_invariants']])
    basis=[E([rat(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    w=-vector(QQ,net['trace_word']);G=matrix(QQ,parent['generic_height_gram'])
    C=sum((n*P for n,P in zip(w,basis)),E(0));h=R(C[0].denominator().sqrt())
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    u=QQ(diag['curve']['u']);q=R(diag['curve']['q']);c=QQ(diag['curve']['scale'])
    assert q.degree()==6 and q.gcd(q.derivative())==1
    x0=QQ(diag['inherited_x_values'][diag['base_pair_index']])
    assert c*q(x0)!=0 and (c*q(x0)).is_square()
    f0,f1,f2=[B[i]+u*A[i] for i in range(3)]
    mn=-f1+E.a1()*f2/2
    # Full exact divisor and surface transport checks.
    polynomials=[]
    for j,(generic,saved) in enumerate(zip(d['generic_r_functions'],d['divisors'])):
        P=basis[j]
        av=A[0]+A[1]*P[0]+A[2]*P[1];bv=B[0]+B[1]*P[0]+B[2]*P[1]
        numerator,denominator=R(generic['numerator']),R(generic['denominator'])
        assert numerator*av+denominator*bv==0
        assert numerator.gcd(denominator)==1
        expected=ZZ(1+G[j,j]-(w*G)[j]);assert generic['degree']==saved['degree']==expected
        g=R(saved['g']);assert g.is_monic() and g.degree()==expected
        assert numerator-u*denominator==g*(numerator-u*denominator).leading_coefficient()
        assert g.gcd(q)==g.gcd(g.derivative())==1
        sy=R(saved['s_mod_g']);assert (sy*sy-c*q)%g==0
        residual=f0+f1*P[0]+f2*P[1]
        assert residual.denominator().gcd(g)==1 and residual.numerator()%g==0
        expected_y=((2*(P[0]+E.b2()/12)+cx)*f2**2-mn**2)/h**3
        assert expected_y.denominator().gcd(g)==1
        assert (expected_y.numerator()-sy*expected_y.denominator())%g==0
        polynomials.append(g)
    # The retained first curve really contains the independently certified EC witness.
    sstar=(c*q(0)).sqrt();assert sstar and sstar*sstar==c*q(0)
    E0=EllipticCurve(QQ,[a(0) for a in E.a_invariants()]);mapped=[]
    m=mn(0)/f2(0)
    for sign in [1,-1]:
        omega=h(0)**3*sign*sstar/f2(0)**2
        xs=(m*m-cx(0)+omega)/2;ys=m*(xs-cx(0))-cy(0)
        xx=xs-E.b2()(0)/12;yy=ys-(E.a1()(0)*xx+E.a3()(0))/2
        P=E0([xx,yy]);mapped.append(list(map(str,P[:2])))
    assert first['reconstructed_point_literal302'] in mapped and first['independence']['rank']==18
    # Independently realize each saved finite-field factor, and use its norm.
    rows=[];checked={}
    assert len(d['trials'])==64
    for trial in d['trials']:
        p=trial['p'];coefficients=list(q)+[x0]+[a for g in polynomials for a in g]
        bad_den=any(a.denominator()%p==0 for a in coefficients)
        if bad_den:
            assert trial['status']=='SKIP_DENOMINATOR';continue
        F=GF(p);S=PolynomialRing(F,'x');f=S(q.list());gs=[S(g.list()) for g in polynomials]
        if f.degree()!=6 or f.gcd(f.derivative())!=1:
            assert trial['status']=='SKIP_BAD_SEXTIC';continue
        if f(F(x0))==0 or f(0)==0 or any(g.gcd(f)!=1 for g in gs):
            assert trial['status']=='SKIP_NONUNIT_DIVISOR';continue
        assert trial['status']=='PASS_LOCAL_KUMMER_BLOCK' and c.valuation(p)%2==0
        factors=[S(row) for row in trial['factors']]
        assert prod(factors)==f.monic() and all(a.is_irreducible() for a in factors)
        raw=[]
        for a in factors:
            degree=a.degree()
            if degree==1:L=F;theta=-a[0]
            else:L=GF(ZZ(p)**degree,name='theta',modulus=a);theta=L.gen()
            anchor=L(x0)-theta
            values=[(-1)**g.degree()*g(theta)/anchor**g.degree() for g in gs]+[-theta/anchor]
            row=[]
            for value in values:
                assert value
                norm=prod(value**(ZZ(p)**j) for j in range(degree))
                leg=F(norm)**((p-1)//2);assert leg in [F(1),F(-1)]
                row.append(int(leg==F(-1)))
            raw.append(row)
        Braw=matrix(GF(2),raw);assert raw==trial['raw_rows']
        quotient=matrix(GF(2),trial['quotient_rows'])
        parity=vector(GF(2),[a.degree()%2 for a in factors])
        assert quotient*parity==0 and quotient.rank()==len(factors)-int(bool(parity))
        block=quotient*Braw;assert [list(map(int,z)) for z in block.rows()]==trial['block_rows']
        rows.extend([list(map(int,z)) for z in block.rows()]);checked[p]=[int(a.degree()) for a in factors]
    assert rows==d['matrix_rows']
    M=matrix(GF(2),rows);H=M[:,:-1];z=vector(GF(2),d['separator'])
    assert H.rank()==d['inherited_character_rank']==16
    assert M.rank()==d['total_character_rank']==17 and z*H==0 and (z*M[:,-1])[0]==1
    # Frobenius gates for rational2-torsion and rational theta characteristics.
    assert checked[83]==[6] and checked[179]==[1,5]
    even_classes={min(mask,mask^63) for mask in range(64) if mask.bit_count()%2==0}
    def permute(mask,permutation):
        return sum(((mask>>j)&1)<<permutation[j] for j in range(6))
    cycle6=[1,2,3,4,5,0]
    assert sum(min(permute(mask,cycle6),permute(mask,cycle6)^63)==mask for mask in even_classes)==1
    triples={min(sum(1<<j for j in s),63^sum(1<<j for j in s)) for s in itertools.combinations(range(6),3)}
    cycle5=[0,2,3,4,5,1]
    assert len(triples)==10
    assert not any(min(permute(mask,cycle5),permute(mask,cycle5)^63)==mask for mask in triples)
    # NS restriction upper bound: normalize degree using O|C=P0, then kill C_min.
    NS=matrix(QQ,19,19);NS[0,0]=-2;NS[0,1]=NS[1,0]=1
    NS[2:,2:]=-G
    D=vector(QQ,[2,5]+list(w));Cmin=vector(QQ,[2,4]+list(w));O=vector(QQ,[1]+[0]*18)
    assert D*NS*Cmin==0 and D*NS*O==1 and O*NS*Cmin==0
    assert matrix(QQ,[O,Cmin]).rank()==2 and NS.rank()==19
    generators=[]
    for j in range(17):
        v=vector(QQ,[1,G[j,j]/2]+[int(i==j) for i in range(17)])
        generators.append(v-(D*NS*v)*O)
    fibre=vector(QQ,[0,1]+[0]*17);generators.append(fibre-2*O)
    N=matrix(QQ,generators)
    assert N.rank()==18 and N*(NS*D)==0 and Cmin in N.row_space()
    beta=(-t*(R([x0,-1]).inverse_mod(q)))%q
    assert q.resultant(beta)/q.leading_coefficient()**beta.degree()==q(0)/q(x0)
    assert (q(0)/q(x0)).is_square()
    return {'classification':'retrospective verified application and new deduction',
            'status':'PASS_INDEPENDENT_INHERITED_RANK17_PLUS_FIRST_WITNESS_JACOBIAN_RANK18',
            'curve_q':list(map(str,q.list())),'curve_scale':str(c),'base_x':str(x0),
            'first_curve_point':["0",str(sstar)],'first_point_maps_to_certified_EC_gain':True,
            'relative_Kummer_beta':list(map(str,beta.list())),
            'relative_Kummer_norm':str(q(0)/q(x0)),
            'inherited_fake_character_rank':16,'with_first_fake_character_rank':17,
            'Jacobian_rational_2_torsion_zero':True,'no_rational_theta_characteristic':True,
            'theta_primes':[83,179],'canonical_minus_2base_not_in_2J':True,
            'inherited_NS_restriction_rank':17,'rank_of_subgroup_with_first_witness':18,
            'whole_Jacobian_rank':'AT_LEAST_18; no upper bound on the whole Jacobian',
            'scope':'The curve is retrospectively selected through the first witness. This is not a blind selector, full Selmer computation, control-panel result, or MW16 cubic ideal-class identification.',
            'limits':{'wall_seconds':25,'prime_trials':64,'passing_blocks':len(checked),
                      'point_searches':0,'class_group_runs':0,'global_Selmer_runs':0,'pilot_changes':0},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,np,dp,rp]},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);r=verify();payload=json.dumps(r,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==payload
    else:OUT.write_text(payload)
    print(r['status'],flush=True)

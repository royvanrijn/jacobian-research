#!/usr/bin/env sage-python
"""Independent one-case replay for the frozen generic-point specificity panel."""
import argparse,hashlib,json,signal,itertools
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_rr_generic_point_controls_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(index):
    cp=DIR/f'case-{index:02d}.json';pp=DIR/'protocol.json'
    result,protocol=[json.loads(p.read_text()) for p in [cp,pp]]
    for doc in [protocol,result]:
        for path,digest in doc['inputs'].items():assert sha(ROOT/path)==digest
    old=json.loads((ART/'det1092_rr_generic_point_controls_v1/protocol.json').read_text())
    for key in ['selector','universal_curve','basepoint','cases','primes','limits']:
        assert protocol[key]==old[key]
    case=protocol['cases'][index];assert result['case']==case and len(protocol['cases'])==9
    parent_path=ART/'curve302_recovered_mw17_parent_v1.json';net_path=ART/'det1092_first_centre_rr_net_v1.json'
    parent,net=[json.loads(p.read_text()) for p in [parent_path,net_path]]
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def rat(row):return K(R(row['numerator']))/R(row['denominator'])
    E=EllipticCurve(K,[rat(row) for row in parent['a_invariants']])
    basis=[E([rat(row) for row in P]) for P in parent['basis_weierstrass_coordinates']]
    A,B=[[R(row) for row in net[key]] for key in ['A','B']]
    w=-vector(QQ,net['trace_word']);G=matrix(QQ,parent['generic_height_gram'])
    C=sum((n*P for n,P in zip(w,basis)),E(0));h=R(C[0].denominator().sqrt())
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2;aa=-E.c4()/48
    selector=K(R(protocol['selector']['u_numerator']))/R(protocol['selector']['u_denominator'])
    S0=basis[0]
    assert selector==-(B[0]+B[1]*S0[0]+B[2]*S0[1])/(A[0]+A[1]*S0[0]+A[2]*S0[1])
    assert max(selector.numerator().degree(),selector.denominator().degree())==7
    tau=QQ(case['parameter']);u=selector(tau)
    assert str(u)==result['u'] and result['v']=='0'
    E0=EllipticCurve(QQ,[a(tau) for a in E.a_invariants()])
    assert E0.discriminant() and E0.j_invariant()==EllipticCurve(QQ,case['elliptic_model']).j_invariant()
    if case['reduced_parameter'] is not None:
        transport=json.loads((ART/'det1092_reduced_parameter_chart_v1/generic-proof.json').read_text())
        a,b,c,d=map(QQ,transport['parameter_matrix']);red=QQ(case['reduced_parameter'])
        assert tau==(a*red+b)/(c*red+d)
    else:assert tau==0
    f0,f1,f2=[B[j]+u*A[j] for j in range(3)]
    mnum=R(-f1+E.a1()*f2/2);L=R(f2/h)
    nx,ny=R(cx*h*h),R(cy*h**3)
    raw=mnum**4-6*nx*mnum*mnum*L*L-8*ny*mnum*L**3-(3*nx*nx+4*R(aa)*h**4)*L**4
    q=R(result['q']);scale=QQ(result['scale'])
    assert raw==h**6*scale*q and q.degree()==6 and q.gcd(q.derivative())==1
    alpha=R(A[2]/h)[0];b=R(B[2]/h);x0=-(b[0]+alpha*u)/b[1]
    assert str(x0)==result['base_x'] and scale*q(x0) and (scale*q(x0)).is_square()
    marked=E0([S0[0](tau),S0[1](tau)])
    assert list(map(str,marked[:2]))==result['marked_elliptic_point']
    assert result['elliptic_basis_word']==[1]+[0]*16
    gs=[]
    for j,P in enumerate(basis):
        r=-(B[0]+B[1]*P[0]+B[2]*P[1])/(A[0]+A[1]*P[0]+A[2]*P[1])
        g=(r.numerator()-u*r.denominator()).monic()
        degree=ZZ(1+G[j,j]-(w*G)[j]);assert degree==result['divisor_degrees'][j]==g.degree()
        assert g.gcd(q)==g.gcd(g.derivative())==1
        y=((2*(P[0]+E.b2()/12)+cx)*f2**2-mnum**2)/h**3
        assert y.denominator().gcd(g)==1
        reduced=(y.numerator()*y.denominator().inverse_mod(g))%g
        assert (reduced**2-scale*q)%g==0
        assert hashlib.sha256(json.dumps(list(map(str,reduced.list()))).encode()).hexdigest()==result['divisor_s_mod_g_sha256'][j]
        incidence=f0+f1*P[0]+f2*P[1]
        assert incidence.denominator().gcd(g)==1 and incidence.numerator()%g==0
        gs.append(g)
        if j==0:ystar=y(tau)
    assert gs[0](tau)==0 and [str(tau),str(ystar)]==result['marked_curve_point']
    assert ystar and ystar**2==scale*q(tau)
    m=mnum(tau)/f2(tau);W=h(tau)**3*ystar/f2(tau)**2
    xx=(m*m-cx(tau)+W)/2;yy=m*(xx-cx(tau))-cy(tau)
    xx-=E.b2()(tau)/12;yy-=(E.a1()(tau)*xx+E.a3()(tau))/2
    assert E0([xx,yy])==marked
    rows=[];checked={}
    assert [r['p'] for r in result['trials']]==protocol['primes'] and len(result['trials'])==64
    for trial in result['trials']:
        p=trial['p'];coefficients=list(q)+[x0,tau]+[a for g in gs for a in g]
        if any(a.denominator()%p==0 for a in coefficients):
            assert trial['status']=='SKIP_DENOMINATOR';continue
        F=GF(p);S=PolynomialRing(F,'x');f=S(q.list());polys=[S(g.list()) for g in gs]
        if f.degree()!=6 or f.gcd(f.derivative())!=1:
            assert trial['status']=='SKIP_BAD_SEXTIC';continue
        if f(F(x0))==0 or f(F(tau))==0 or any(g.gcd(f)!=1 for g in polys):
            assert trial['status']=='SKIP_NONUNIT_DIVISOR';continue
        assert trial['status']=='PASS_LOCAL_KUMMER_BLOCK' and scale.valuation(p)%2==0
        factors=[S(row) for row in trial['factors']]
        assert prod(factors)==f.monic() and all(a.is_irreducible() for a in factors)
        bits=[]
        for a in factors:
            degree=a.degree()
            if degree==1:L=F;theta=-a[0]
            else:L=GF(ZZ(p)**degree,name='zeta',modulus=a);theta=L.gen()
            anchor=L(x0)-theta
            values=[(-1)**g.degree()*g(theta)/anchor**g.degree() for g in polys]+[(L(tau)-theta)/anchor]
            row=[]
            for value in values:
                assert value
                norm=prod(value**(ZZ(p)**i) for i in range(degree))
                leg=F(norm)**((p-1)//2);assert leg in [F(1),F(-1)]
                row.append(int(leg==F(-1)))
            bits.append(row)
        assert bits==trial['raw_rows'];raw=matrix(GF(2),bits)
        Q=matrix(GF(2),trial['quotient_rows']);parity=vector(GF(2),[a.degree()%2 for a in factors])
        assert Q*parity==0 and Q.rank()==len(factors)-int(bool(parity))
        block=Q*raw;assert [list(map(int,z)) for z in block.rows()]==trial['block_rows']
        rows.extend([list(map(int,z)) for z in block.rows()]);checked[p]=[int(a.degree()) for a in factors]
    assert rows==result['matrix_rows'];M=matrix(GF(2),rows);H=M[:,:-1]
    assert H.rank()==result['inherited_character_rank'] and M.rank()==result['with_marked_character_rank']
    six,five=result['theta_primes']['six_cycle'],result['theta_primes']['one_plus_five']
    complete=(H.rank()==16 and M.rank()==17 and six is not None and five is not None)
    if complete:
        z=vector(GF(2),result['separator']);assert z*H==0 and (z*M[:,-1])[0]==1
        assert checked[six]==[6] and checked[five]==[1,5]
        def moved(mask,p):return sum(((mask>>j)&1)<<p[j] for j in range(6))
        even={min(mask,mask^63) for mask in range(64) if mask.bit_count()%2==0}
        assert sum(min(moved(mask,[1,2,3,4,5,0]),moved(mask,[1,2,3,4,5,0])^63)==mask for mask in even)==1
        triples={min(sum(1<<j for j in v),63^sum(1<<j for j in v)) for v in itertools.combinations(range(6),3)}
        assert not any(min(moved(mask,[0,2,3,4,5,1]),moved(mask,[0,2,3,4,5,1])^63)==mask for mask in triples)
    NS=matrix(QQ,19,19);NS[0,0]=-2;NS[0,1]=NS[1,0]=1;NS[2:,2:]=-G
    D=vector(QQ,[2,5]+list(w));Cmin=vector(QQ,[2,4]+list(w));O=vector(QQ,[1]+[0]*18)
    assert NS.rank()==19 and D*NS*Cmin==0 and D*NS*O==1 and matrix(QQ,[O,Cmin]).rank()==2
    # Record genuine Selmer information, not an uncomputed group dimension.
    return {'classification':'verified application and new deduction; generic-input specificity control',
            'status':'PASS_GENERIC_ELLIPTIC_POINT_GIVES_NON_GENERIC_RATIONAL_JACOBIAN_CLASS' if complete else 'VERIFIED_INCONCLUSIVE_FROZEN_LOCAL_PANEL',
            'label':case['label'],'tau':str(tau),'generic_selector_degree':7,
            'inherited_character_rank':int(H.rank()),'with_marked_character_rank':int(M.rank()),
            'inherited_NS_image_rank':17 if complete else 'UNRESOLVED_BY_THIS_PANEL',
            'subgroup_with_marked_Jacobian_class_rank':18 if complete else 'UNRESOLVED_BY_THIS_PANEL',
            'Jacobian_class_outside_generic_NS_rational_span':complete,
            'elliptic_point_is_exactly_generic_basis_0':True,'elliptic_class_mod_MW17':'ZERO',
            'Jacobian_class_is_rational':True,'true_Kummer_class_in_Selmer':True,'image_in_Sha_2':'ZERO',
            'full_Selmer_group':'NOT_COMPUTED','other_Sha_classes':'NOT_COMPUTED',
            'strict_relative_elliptic_class':'ZERO modulo MW17; not a new strict exceptional direction',
            'passing_prime_blocks':len(checked),'theta_primes':result['theta_primes'],
            'interpretation':'If separated, non-generic Jacobian class is insufficient to imply an extra elliptic direction: this exact marked point is generic.',
            'global_family_corollary':'A passing specialization certifies rank at least18 over Q(t) for Jac(s^2=c*q(T;r0(t),0)); the original elliptic family has exact generic rank17.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [cp,pp,parent_path,net_path]},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=0);args=p.parse_args()
    r=verify(args.index);out=DIR/f'case-{args.index:02d}-replay.json';text=json.dumps(r,indent=2,sort_keys=True)+'\n'
    if out.exists():assert out.read_text()==text
    else:out.write_text(text)
    print(r['label'],r['status'],flush=True)

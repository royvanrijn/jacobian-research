#!/usr/bin/env python3
"""Universal polynomial identities and independent trace-form verification."""
import argparse
from pathlib import Path
import retrospective as r
import two_cover_pencil_geometry as pencils

OUTPUT=r.OUT/'rank_jump_two_cover_pencil_verification_v1.json'


def universal():
    from sage.all import QQ,PolynomialRing,matrix,vector
    R=PolynomialRing(QQ,['a','b','c','v0','v1','v2','lam','mu'])
    a,b,c,v0,v1,v2,lam,mu=R.gens()
    # Multiplication by theta, directly from theta^3=-a*theta^2-b*theta-c.
    H=matrix(R,[[0,0,-c],[1,0,-b],[0,1,-a]]);I=matrix.identity(R,3)
    T=H+a*I;M=v0*I+v1*H+v2*H**2
    S=matrix(R,[[0,0,1],[0,1,-a],[1,-a,a*a-b]])
    Q2=S*M;Q1=Q2*T
    assert Q2==Q2.transpose() and Q1==Q1.transpose() and S.det()==-1
    assert Q2.det()==-M.det()
    norm_poly=lam**3+2*a*lam**2*mu+(a*a+b)*lam*mu**2+(a*b-c)*mu**3
    assert (lam*I+mu*T).det()==norm_poly
    assert (lam*Q2+mu*Q1).det()==-M.det()*norm_poly
    e=vector(R,[1,0,0]);cyclic=matrix(R,[e,T*e,T*T*e]).transpose()
    assert cyclic.det()==1
    # Multiplication by an arbitrary element gives the labelled congruence.
    U=I+2*H+3*H**2
    assert U.transpose()*Q2*U==S*M*U**2
    assert U.transpose()*Q1*U==S*M*U**2*T
    return {'status':'PASS','coefficient_ring':str(R),
        'norm_pencil_polynomial':str(norm_poly),'cyclic_basis_determinant':'1',
        'first_conic_base_Gram':[[str(x) for x in row] for row in S.rows()],
        'identities':['Q2(beta)=Q2(1)*M_beta','Q1(beta)=Q2(beta)*(H+a*I)',
            'det(Q2(beta))=-Norm(beta)','det(lam*Q2+mu*Q1)=-Norm(beta)*Norm(lam+mu*(theta+a))',
            'M_u^T*Qj(beta)*M_u=Qj(beta*u^2), j=1,2 (explicit nonconstant u regression)']}


def compute():
    from sage.all import QQ,PolynomialRing,matrix,vector,GF
    inp=r.read(pencils.INPUT);data=r.read(pencils.OUTPUT)
    assert data['bindings']==pencils.bindings()
    for name,digest in inp['source_hashes'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    symbolic=universal();rows=[]
    for source,result in zip(inp['rows'],data['rows'],strict=True):
        assert result['status']=='PASS' and source['case_index']==result['case_index']
        c,b,a,_=map(QQ,source['cubic_ascending'])
        H=matrix(QQ,[[0,0,-c],[1,0,-b],[0,1,-a]]);I=matrix.identity(QQ,3);T=H+a*I
        Fprime=3*H**2+2*a*H+b*I;inv=Fprime.inverse()
        assert [list(map(str,x)) for x in T.rows()]==result['common_operator']
        R=PolynomialRing(QQ,'z');z=R.gen();f=z**3+a*z*z+b*z+c
        assert str(f.discriminant())==result['cubic_discriminant']
        irred=next((p for p in r.primes(199) if f.change_ring(GF(p)).is_irreducible()),None)
        assert irred is not None
        signatures=[]
        for item,rec in zip(source['classes'],result['classes'],strict=True):
            beta=list(map(QQ,item['beta']));M=beta[0]*I+beta[1]*H+beta[2]*H**2
            norm=M.det();assert norm==QQ(item['norm_root'])**2==QQ(rec['norm'])
            expected2=matrix(QQ,3,3,lambda i,j:(M*H**(i+j)*inv).trace())
            expected1=matrix(QQ,3,3,lambda i,j:(M*H**(i+j)*T*inv).trace())
            Q2=matrix(QQ,rec['Q2_Gram']);Q1=matrix(QQ,rec['Q1_Gram'])
            assert Q2==expected2 and Q1==expected1
            assert Q2.det()==-norm and Q2.inverse()*Q1==T
            # Five finite samples prove the claimed degree <=4 determinant identity;
            # the leading coefficient is also checked by the rank-three first quadric.
            for signed in rec['signs']:
                sigma=signed['sigma'];coeff=list(map(QQ,signed['pencil_determinant_coefficients']))
                assert coeff[0]==0 and coeff==[-sigma*norm*x for x in map(QQ,result['normalized_pencil_coefficients'])]
                for lam in [-2,-1,0,1,2]:
                    direct=sigma*(lam*Q2+Q1).det()
                    assert direct==sum(coeff[i]*QQ(lam)**(4-i) for i in range(5))
            signatures.append(item['kind'])
        rows.append({'case_index':source['case_index'],'id':source['id'],'status':'PASS',
            'class_count':len(signatures),'generic_class_count':signatures.count('generic'),
            'relative_class_count':signatures.count('relative'),'cubic_irreducibility_prime':irred,
            'trace_forms_verified':True,'both_signs_smooth':True,
            'normalized_pencil_coefficients':result['normalized_pencil_coefficients']})
        print('verified trace forms and pencil',source['case_index'],len(signatures),flush=True)
    assert sum(x['class_count'] for x in rows)==50
    # Existing exact rational/Sha small covers have the identical first conic.
    small=r.read(pencils.small.OUTPUT);small_result=data['rows'][6]
    for src,rec in zip(small['records'],small_result['classes'],strict=True):
        z=vector(QQ,src['twist_rational_point'][:3]);h=QQ(src['twist_rational_point'][3])
        Q1=matrix(QQ,rec['Q1_Gram']);Q2=matrix(QQ,rec['Q2_Gram'])
        assert (z*Q2*z)==0 and (z*Q1*z)==h*h and h
    return {'schema':'rank-jump.two-cover-pencil-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),pencils.INPUT,pencils.OUTPUT,pencils.PROTOCOL,pencils.small.OUTPUT)},
        'universal':symbolic,'rows':rows,'small_control_common_conics_split':True,
        'boundary':'Source rank/CT solubility certificates are retained inputs, not recomputed. This independently verifies the common trace geometry and explicit conic points.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS universal and independent trace geometry')
    else:r.write_new(OUTPUT,data)

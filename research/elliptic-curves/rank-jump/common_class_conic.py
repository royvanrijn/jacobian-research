#!/usr/bin/env python3
"""Universal common conic and twist double covers for one cubic Kummer class."""
import argparse
from pathlib import Path
import retrospective as r

OUTPUT=r.OUT/'rank_jump_common_class_conic_v1.json'


def compute():
    from sage.all import QQ,PolynomialRing,matrix
    R=PolynomialRing(QQ,names=('A','B','h0','h1','h2','z0','z1','z2','l','m','d','w'))
    A,B,h0,h1,h2,z0,z1,z2,l,m,d,w=R.gens();h=[h0,h1,h2];z=[z0,z1,z2]
    def mul(a,b):
        out=[R(0)]*5
        for i in range(3):
            for j in range(3):out[i+j]+=a[i]*b[j]
        for k in (4,3):out[k-2]-=A*out[k];out[k-3]-=B*out[k]
        return out[:3]
    def norm(a):
        return matrix(R,[mul(a,[R(int(i==j)) for i in range(3)]) for j in range(3)]).det()
    q=mul(h,mul(z,z));N=norm(h);Nz=norm(z)
    assert norm(q)==N*Nz*Nz
    H=[matrix(R,3,3,lambda i,j:q[k].derivative(z[i]).derivative(z[j])/2) for k in range(3)]
    assert all(sum(H[k][i,j]*z[i]*z[j] for i in range(3) for j in range(3))==q[k] for k in range(3))
    assert H[2].det()==-N
    pencil=(l*H[2]+m*H[1]).det()
    assert pencil==-N*(l**3+A*l*m*m-B*m**3)
    assert norm([l,-d*w*w,R(0)])==l**3+A*d*d*w**4*l+B*d**3*w**6
    encode=lambda M:[[str(M[i,j]) for j in range(3)] for i in range(3)]
    return {'schema':'rank-jump.common-class-conic.v1','status':'PASS',
        'cubic':'theta^3 + A*theta + B','kummer_class':'h0+h1*theta+h2*theta^2',
        'norm_h':str(N),'norm_z':str(Nz),'quadratic_forms':list(map(str,q)),
        'quadratic_matrices':list(map(encode,H)),
        'common_conic':'Q2(z0,z1,z2)=0',
        'twist_cover':'Q2(z)=0, Q1(z)+d*w^2=0 in P^3',
        'point_map':'x=Q0(z)/w^2, y=n*Norm(z)/w^3, where Norm(h)=n^2',
        'target':'y^2=x^3+d^2*A*x+d^3*B',
        'pencil_determinant':str(m*d*pencil),
        'conditions':'Characteristic0; d*Norm(h)*(-4*A^3-27*B^2) nonzero; Norm(h)=n^2. The class is nonzero when using the rational-point equivalence without an infinity exception.',
        'conclusion':'Every scalar twist cover of the same Kummer class has the same conic quotient. If the conic has a rational point it can be parametrized, giving quadratic twists of a single binary quartic. Conic solubility alone is insufficient.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),Path(r.__file__))}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS universal norm, common-conic and determinant-pencil identities')

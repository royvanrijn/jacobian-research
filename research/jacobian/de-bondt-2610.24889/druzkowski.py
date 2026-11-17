"""Direct integer matrix and rational collision for the 90-variable lift.

The construction uses only a cube decomposition and Sylvester's determinant
identity. It does not need a square matrix completion over a Hermite ring.
"""
from __future__ import annotations
import json
import sympy as s
from verify import X, ROOT, cubic12, homogeneous18, evaluate


def cube(v):return v.applyfunc(lambda a:a**3)


def main():
    data=json.loads((ROOT/'polarization.json').read_text())
    L=s.Matrix(data['forms'])
    B=s.Matrix([[s.Rational(c) for c in row] for row in data['B']])
    for i in range(17):
        for j in range(B.cols): B[i,j]*=6
    assert all(c.q==1 for c in B)
    e=s.zeros(18,1);e[17]=1
    padding=next(j for j in range(L.rows) if L.row(j)==e.T)
    assert B.row(17)==s.eye(L.rows).row(padding)
    shear=B[:,padding]-e
    L0=L.copy();L0[padding,:]=s.zeros(1,18)
    A=L0*B
    assert all(c.q==1 for c in A)
    assert A.row(padding)==s.zeros(1,A.cols)
    _,pivots=B.rref()
    assert len(pivots)==18
    C=s.zeros(B.cols,18)
    inv=B[:,list(pivots)].inv()
    for i,j in enumerate(pivots):C[j,:]=inv[i,:]
    assert B*C==s.eye(18)
    _,k,points=cubic12(2)
    h,ctail=homogeneous18(k)
    # Rational collisions for X+lambda H at lambda=6: tau=1/lambda.
    lifted=[s.Matrix(p+[-6*a for a in evaluate(ctail,p)]+[s.Rational(1,6)]) for p in points]
    def base(z):return z+6*s.Matrix(evaluate(h,list(z)))
    assert base(lifted[0])==base(lifted[1])
    def sheared_base(z):return base(z)-shear*z[17]**3
    u=C*lifted[0]
    trial=C*lifted[1]
    # Correct along ker(B), where the lift acts as the identity.
    image=u+cube(A*u)
    v=image-cube(A*trial)
    assert B*v==lifted[1]
    assert u!=v and u+cube(A*u)==v+cube(A*v)
    assert B*image==sheared_base(lifted[0])
    assert B*(trial+cube(A*trial))==sheared_base(lifted[1])
    # B*F(Y)=S(G_6(BY)), S(z)=z-shear*z18^3. The last entry
    # of shear is zero, so S is elementary with determinant one.
    # Sylvester reduces the 90-variable determinant to this sheared base.
    # The exact cube identity has already been checked in polarize.py.
    r=A.rank()
    powers=[];current=A
    for n in range(1,5):
        powers.append({'power':n,'zero':current.is_zero_matrix,'rank':current.rank()})
        if current.is_zero_matrix:break
        current=current*A
    payload={'dimension':A.rows,'coefficient_ring':'Z','A_rank':r,
             'padding_index_1_based':padding+1,'A_powers':powers,
             'A':[[int(c) for c in row] for row in A.tolist()],
             'u':[str(c) for c in u],'v':[str(c) for c in v],
             'image':[str(c) for c in image],
             'target_shear_vector':[str(c) for c in shear],
             'projection_identity':'B F = S(G6(BY)); S(z)=z-shear*z18^3',
             'lifted_base_points':[[str(c) for c in z] for z in lifted],
             'projection_B':[[int(c) for c in row] for row in B.tolist()]}
    (ROOT/'druzkowski90.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k in ['dimension','coefficient_ring','A_rank','padding_index_1_based','A_powers']},indent=2))
    print('EXACT 90-VARIABLE INTEGER MATRIX AND RATIONAL COLLISION VERIFIED')


if __name__=='__main__':main()

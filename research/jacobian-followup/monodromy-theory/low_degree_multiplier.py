"""Exact low-degree divisor identities and hyperplane normal forms."""
import sympy as S
from two_image_projective import conjugate
def taylor_data(q,v,xs):
    v=S.Matrix(v);xx=S.Matrix(xs);W=v-xx;B=v.jacobian(xs)
    dq=S.Matrix([[S.diff(q,x) for x in xs]]);H=S.hessian(q,xs)
    Wq=(dq*W)[0]
    W2q=(S.Matrix([[S.diff(Wq,x) for x in xs]])*W)[0]
    qv=q.subs(dict(zip(xs,v)),simultaneous=True)
    return dict(v=v,W=W,B=B,dq=dq,H=H,Wq=Wq,W2q=W2q,qv=qv,
      chain_identity=W2q-(W.T*H*W)[0]-(dq*(B-S.eye(len(xs)))*W)[0],
      taylor_identity=qv-q-Wq-(W.T*H*W)[0]/2)
def hyperplane_forms():
    x,Y,Z,a,b,c,d,e,f,g,k1,k2,k3=S.symbols("x Y Z a b c d e f g k1 k2 k3")
    xs=[x,Y,Z];M=S.Matrix([[a,0,0],[b,c,d],[e,f,g]]);kk=S.Matrix([k1,k2,k3])
    L=S.Matrix([[1,0,0],[0,1,0],[0,-g/d,1]])
    H_A=S.Matrix([0,x*x,Y*Y]);H_B=S.Matrix([0,x*x,x*Y])
    oldA=M*(S.Matrix(xs)+H_A)+kk;oldB=M*(S.Matrix(xs)+H_B)+kk
    vA=conjugate(oldA,L,xs);vB=conjugate(oldB,L,xs)
    return dict(xs=xs,parameters=[a,b,c,d,e,f,g,k1,k2,k3],M=M,L=L,r=f-g*c/d,
         oldA=oldA,oldB=oldB,vA=vA,vB=vB,WA=vA-S.Matrix(xs),WB=vB-S.Matrix(xs))

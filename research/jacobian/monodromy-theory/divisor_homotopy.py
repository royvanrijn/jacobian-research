"""Formal line-homotopy identities and an affine rank4 survivor."""
import sympy as S
def truncate(e,t,N):
    return S.Add(*[co*t**mon[0] for mon,co in S.Poly(S.expand(e),t).terms() if mon[0]<N])
def formal_inverse(W,xs,t,N):
    W=S.Matrix(W);xx=S.Matrix(xs);G=xx.copy()
    for i in range(N):
        G=(xx-t*W.subs(dict(zip(xs,G)),simultaneous=True)).applyfunc(lambda e:truncate(e,t,N))
    return G
def residue_moments(q,v,xs,count):
    v=S.Matrix(v);W=v-S.Matrix(xs);B=v.jacobian(xs)
    dq=S.Matrix([[S.diff(q,x) for x in xs]])
    return [(dq*B**j*W)[0].expand() for j in range(count)]
def affine_survivor():
    u,x,y,z,l,T=S.symbols("u x y z l T");ps=[u,x,y,z]
    modulus=S.Poly(l**4-5*l**3+10*l**2-10*l+5,l)
    def red(e):
        num,den=S.fraction(S.cancel(e))
        return S.cancel(S.rem(S.Poly(num,l),modulus).as_expr()/S.rem(S.Poly(den,l),modulus).as_expr())
    c2=1-2*l;c1=1-2*l+2*l*l;c0=1-2*l+2*l*l-l**3
    C=S.Matrix([[0,0,-c0],[1,0,-c1],[0,1,-c2]])
    M=S.diag(S.Integer(1),C)
    F=(x*S.eye(3)+y*C+z*C*C).det(method="domain-ge")
    F=S.expand(red(F));G=u*F
    X=S.Matrix([-l*G,*list(G*(M-l*S.eye(4))*S.Matrix(ps))])
    return dict(ps=ps,l=l,T=T,modulus=modulus,red=red,C=C,M=M,F=F,G=G,X=X,V=u*M*S.Matrix(ps))

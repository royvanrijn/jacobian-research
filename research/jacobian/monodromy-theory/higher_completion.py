"""Finite weighted-jet rectification and all-rank tautological completion."""
import sympy as S

def truncate(expr,t,D):
    poly=S.Poly(S.expand(expr),t)
    return S.Add(*[c*t**m[0] for m,c in poly.terms() if m[0]<D])

def rectify(G,variables):
    """G(0)=0,DG(0)=I. Return an actual triangular polynomial automorphism."""
    t=variables[0];us=list(variables[1:]);r=len(us);D=2*r
    if r<1:raise ValueError('at least two outputs required')
    origin={v:0 for v in variables}
    if G.subs(origin)!=S.zeros(r+1,1) or G.jacobian(variables).subs(origin)!=S.eye(r+1):
        raise ValueError('identity linear jet required')
    h=S.zeros(r,1)
    for iteration in range(D):
        value=G[1:,:].subs(dict(zip(us,h)),simultaneous=True)
        h=(h-value).applyfunc(lambda e:truncate(e,t,D))
    M=G[1:,:].jacobian(us).subs(dict(zip(us,h)),simultaneous=True)
    M=M.applyfunc(lambda e:truncate(e,t,D));U=M.copy();A=S.eye(r)
    modulus=S.Poly(t**D,t)
    for i in range(r-1,0,-1):
        for j in range(i):
            inverse=S.invert(S.Poly(U[i,i],t),modulus).as_expr()
            coefficient=truncate(U[i,j]*inverse,t,D)
            for row in range(r):
                U[row,j]=truncate(U[row,j]-coefficient*U[row,i],t,D)
                A[row,j]=truncate(A[row,j]-coefficient*A[row,i],t,D)
    phi=S.Matrix([t,*list(h+A*S.Matrix(us))])
    inverse=S.Matrix([t,*list(A.inv()*(S.Matrix(us)-h))])
    return dict(phi=phi,inverse=inverse,h=h,A=A,M=M,U=U,D=D,
                weights=[1]+list(range(r+1,2*r+1)))

def weighted_order(expr,variables,weights):
    poly=S.Poly(S.expand(expr),*variables)
    if poly.is_zero:return None
    return min(sum(a*b for a,b in zip(m,weights)) for m,c in poly.terms())

def weighted_stage(X,parameters,g,us):
    r=len(parameters);aa=list(range(r,2*r));ee=[1]+list(range(r+1,2*r+1))
    substitutions={p:g**a*u for p,a,u in zip(parameters,aa,us)}
    swept=X+g*S.Matrix([1,*parameters])
    E=S.Matrix([S.cancel(e.subs(substitutions,simultaneous=True)/g**b)
                for e,b in zip(swept,ee)])
    return dict(E=E,source_weights=aa,output_weights=ee,S=swept)

def companion_surface(parameters):
    r=len(parameters);tt=S.Dummy('q')
    polynomial=sum(p*tt**(i+1) for i,p in enumerate(parameters))
    coefficient=S.expand(sum((-1)**(j+1)*polynomial**j/S.Integer(j)
                             for j in range(1,r+2))).coeff(tt,r+1)
    return S.Matrix([*parameters,-coefficient])

def suspended_cubic(parameters):
    w=parameters[0]
    return S.Matrix([2*w+w*w,w*w+S.Rational(2,3)*w**3,*parameters[2:],0]) if len(parameters)>1 else S.Matrix([2*w+w*w,w*w+S.Rational(2,3)*w**3])

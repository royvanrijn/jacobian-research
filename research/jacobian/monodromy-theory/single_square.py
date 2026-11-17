"""Single-square directions, transverse fixed-locus and base-line gates."""
import sympy as S
def direction(M,ps):
    p=S.Matrix(ps);u,x,y,z=ps;a=S.Matrix(M)[:,3]
    return u*S.Matrix(M)*p+a*x*x
def wedge(v,w):
    return S.Matrix([v[i]*w[j]-v[j]*w[i] for i in range(len(v)) for j in range(i+1,len(v))])
def triple(v,w,t):
    from itertools import combinations
    return S.Matrix([S.Matrix.hstack(S.Matrix(v),S.Matrix(w),S.Matrix(t)).extract(ij,[0,1,2]).det() for ij in combinations(range(len(v)),3)])
def fixed_quadric_form(ps,kappa,c,row):
    p=S.Matrix(ps);u,x,*_=ps;c=S.Matrix(c)
    D=x*x+u*(S.Matrix([list(row)])*p)[0]
    M=kappa*S.eye(len(ps))+c*S.Matrix([list(row)])
    return dict(D=D,M=M,V=kappa*u*p+c*D)
def rank3_control(ps,beta):
    u,x,y,z=ps;lam=1+beta
    M=S.diag(1,lam,beta,lam);V=direction(M,ps);F=y
    X=S.Matrix([-lam*u*F,*list(F*(V-lam*u*S.Matrix(ps)))])
    P=X[1:,:].jacobian(ps)-S.Matrix(ps)*X[:1,:].jacobian(ps)
    return dict(M=M,V=V,F=F,lam=lam,X=X,P=P)

"""Uniform polar companions and fail-closed backward section certificates."""
import sympy as s

def parameter_product(n,gamma,x):
    return s.Poly(s.prod(gamma+2+3*j-x for j in range(n)),x)

def companions(B,n,gamma,z,y):
    J=s.Poly(sum(s.binomial(n,k)*s.prod(gamma+2+3*j-k for j in range(n))*y**k for k in range(n+1)),y)
    K=s.Poly(sum(s.ff(n,j)*B.nth(j)*z**j for j in range(n+1)),z)
    return J,K

def backward_certificate(parts,N,u):
    q={N-r:s.Poly(s.cancel(p.as_expr()/p.nth(0)),u,field=True) for r,p in enumerate(parts)}
    weights={};formal=False
    for j in range(N,2,-1):
        diff=q[j]-q[j-1]
        pivot=s.factor(diff.nth(1))
        if pivot==0:return {'status':'inconclusive_zero_pivot','at':j,'weights':weights,'q':q}
        formal=formal or bool(pivot.free_symbols)
        weights[j]=pivot
        q[j-3]=s.Poly(s.cancel(diff.as_expr()/(pivot*u)),u,field=True)
        assert q[j-3].degree()<=(j-3)//3 and q[j-3].nth(0)==1
    return {'status':'formal_nonzero_conditions' if formal else 'certified','weights':weights,'q':q}

def replay_certificate(certificate,N,u):
    if certificate['status']!='certified':return False
    q=certificate['q'];weights=certificate['weights']
    if any(q[j].as_expr()!=1 for j in [0,1,2]):return False
    for j in range(3,N+1):
        if weights[j]==0:return False
        if s.cancel(q[j].as_expr()-q[j-1].as_expr()-weights[j]*u*q[j-3].as_expr())!=0:return False
        if q[j].degree()>j//3 or q[j].nth(0)!=1:return False
    return True

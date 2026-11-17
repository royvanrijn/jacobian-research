"""Exact formulas for the signed-alternating affine-four-space family."""
import sympy as S

def family(m):
    if not isinstance(m,int) or m<2:raise ValueError("m must be an integer >=2")
    a,g,u,v,w,z,q,A,B=S.symbols("a g u v w z q A B")
    am=S.Rational(m*m,2*(2*m+1)*(m+1)**2)
    bm=S.Rational(m,2*(m+1)*(m+2))
    minimal=am*a*a+bm*a+S.Rational(1,24)
    def red(expr):
        return S.Poly(S.expand(expr),a).rem(S.Poly(minimal,a)).as_expr().expand()
    def inv(expr):return S.invert(red(expr),minimal,a)
    h=a*q**m+q
    b=a*S.Rational(m,m+1)*q**(m+1)+q*q/2
    d=a*a*S.Rational(m,2*m+1)*q**(2*m+1)+a*S.Rational(m+1,m+2)*q**(m+2)+q**3/3
    sigma=h.subs(q,w*w)*w-z
    X=S.Matrix([sigma,b.subs(q,w*w),d.subs(q,w*w)-sigma*sigma/2])
    Delta=S.Matrix([1,w,z]);sweep=X+g*Delta
    E=[S.cancel(sweep[i].subs({w:g*u,z:g**3*v},simultaneous=True)/g**e)
       for i,e in enumerate([1,2,4])]
    g0=-b.subs(q,1);gi=inv(g0)
    theta0=S.Matrix([g0,gi,red((h.subs(q,1)-b.subs(q,1))*gi**3)])
    base=dict(zip([g,u,v],theta0))
    def at(expr):return red(expr.subs(base))
    J=S.Matrix([[at(S.diff(e,t)) for t in [g,u,v]] for e in E])
    detJ=red(J.det())
    L=J.adjugate().applyfunc(red).applyfunc(lambda e:red(e*inv(detJ)))
    def directional(e,c):return sum(S.diff(e,t)*co for t,co in zip([g,u,v],c))
    c1,c2,c3=[L[:,i] for i in range(3)]
    H11=at(directional(directional(E[2],c1),c1))
    H12=at(directional(directional(E[2],c1),c2))
    H13=at(directional(directional(E[2],c1),c3))
    T111=at(directional(directional(directional(E[2],c1),c1),c1))
    f20=red(-H11/2);f11=red(-H12);f30=red(-T111/6-f20*H13)
    Q=S.expand((A-b)**2+2*q*h*(A-b)+2*q*d-2*q*B)
    leading=-a*a*S.Rational(m*m,(2*m+1)*(m+1)**2)
    return dict(m=m,n=2*m+2,a=a,g=g,u=u,v=v,w=w,z=z,q=q,A=A,B=B,
        minimal=minimal,red=red,inv=inv,h=h,b=b,d=d,sigma=sigma,
        X=X,Delta=Delta,sweep=sweep,E=E,theta0=theta0,base=base,
        J=J,L=L,f20=f20,f11=f11,f30=f30,H11=H11,H12=H12,H13=H13,
        T111=T111,Q=Q,leading=leading,g0=g0)

"""Primitive cubic wave families, including a curved denominator."""
import sympy as S

def reconstruct(A,b,c,w,z):
    lam=A*b*b;mu=A*b*c;nu=A*c*c
    x1z=-S.diff(lam,w)-S.diff(mu,z)
    x1w=S.diff(mu,w)+S.diff(nu,z)
    X1=S.expand(S.integrate(x1z,z)+S.integrate(x1w.subs(z,0),w))
    X2=S.expand(S.integrate(w*S.diff(X1,z)+lam,z)+S.integrate((w*S.diff(X1,w)-mu).subs(z,0),w))
    X3=S.expand(S.integrate(z*S.diff(X1,z)+mu,z)+S.integrate((-nu).subs(z,0),w))
    return S.Matrix([X1,X2,X3])

def curved_cubic(w,z,D=0):
    b=w*w-S.Rational(4,3)*z
    c=w**3-w*z+D
    return dict(b=b,c=c,X=reconstruct(S.Integer(1),b,c,w,z))

def curved_stage(g,u,v,D=1):
    w,z=S.symbols("w z")
    X=curved_cubic(w,z,D)["X"]
    source={w:g*g*u,z:g*v}
    return S.Matrix([
      S.cancel((X[0]+g).subs(source,simultaneous=True)/g),
      S.cancel((X[2]+g*z).subs(source,simultaneous=True)/g**2),
      S.cancel((X[1]+g*w).subs(source,simultaneous=True)/g**3)])

def marked_polynomials(t):
    P=-S.Rational(4,15)+t-S.Rational(4,3)*t*t+S.Rational(16,27)*t**3
    Q=-S.Rational(1,7)+S.Rational(2,5)*t-S.Rational(1,3)*t*t
    d=3*P
    G=S.Rational(3,5)-S.Rational(5,3)*t+S.Rational(4,3)*t*t+S.Rational(4,3)*d*t
    H=89600*t**6-336000*t**5+554400*t**4-509040*t**3+274050*t*t-82215*t+10854
    return dict(P=P,Q=Q,d=S.expand(d),G=S.expand(G),H=H)

def marked_stage_data():
    """Exact coefficients in Q[t]/H; gamma0 is a cube root of q0."""
    t,q,u,v=S.symbols("t q u v")
    marked=marked_polynomials(t);H=S.Poly(marked["H"],t,domain=S.QQ)
    def red(expr):
        num,den=S.fraction(S.cancel(expr))
        N=S.Poly(num,t,domain=S.QQ).rem(H)
        D=S.Poly(den,t,domain=S.QQ).rem(H)
        return (N*S.invert(D,H)).rem(H).as_expr()
    d,G=marked["d"],marked["G"]
    q0=red(-G**3/d**5);u0=red(d**3/G**2);v0=red(-t*d/G)
    g=S.symbols("gamma")
    E=curved_stage(g,u,v)
    Eq=S.Matrix([S.Add(*[co*q**(a//3)*u**b*v**c
                   for (a,b,c),co in S.Poly(expr,g,u,v).terms()])
                  for expr in E])
    bases=[q0,u0,v0]
    powers=[[S.Integer(1)] for _ in bases]
    def evaluate(expr):
        poly=S.Poly(expr,q,u,v)
        total=S.Integer(0)
        for exps,co in poly.terms():
            term=co
            for axis,power in enumerate(exps):
                while len(powers[axis])<=power:
                    powers[axis].append(red(powers[axis][-1]*bases[axis]))
                term=red(term*powers[axis][power])
            total=red(total+term)
        return total
    B=S.Matrix.hstack(3*q*Eq.diff(q),Eq.diff(u),Eq.diff(v)).applyfunc(evaluate)
    determinant=red(B.det())
    R=B.adjugate().applyfunc(lambda expr:red(expr/determinant))
    a,bb,cc=R[:,0]
    first_q=red(3*q0*a);second_q=red(6*q0*a*a)
    E3=Eq[2]
    second=red(first_q**2*evaluate(S.diff(E3,q,2))
      +second_q*evaluate(S.diff(E3,q))
      +2*first_q*bb*evaluate(S.diff(E3,q,u))
      +2*first_q*cc*evaluate(S.diff(E3,q,v))
      +bb**2*evaluate(S.diff(E3,u,2))
      +2*bb*cc*evaluate(S.diff(E3,u,v))
      +cc**2*evaluate(S.diff(E3,v,2)))
    return dict(t=t,q=q,u=u,v=v,H=H,red=red,evaluate=evaluate,
      q0=q0,u0=u0,v0=v0,E_q=Eq,B=B,R=R,determinant=determinant,
      H11=second,f20=red(-second/2),d=d,G=G)

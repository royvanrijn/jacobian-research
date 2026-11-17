"""Exact formulas for the minimal rational-pole sweep and its affine lift."""
import sympy as S
def surface(R,w,z,k):
    I=S.integrate(R,w);J=S.integrate(w*R,w);K0=S.integrate(R*R,w)
    X=S.Matrix([-(k+2)*w*z+w*R+2*k*I,
        k*k*w*w*z+w*w*R+2*(k-1)*J,-w*z*z+w*R*z-K0])
    return X.applyfunc(S.cancel),I,J,K0
def general_surface(r,R,w,z,k):
    if not isinstance(r,int) or r<2:raise ValueError("r must be an integer >=2")
    I=S.integrate(w**(r-2)*R,w)
    J=S.integrate(w**(r-1)*R,w)
    K0=S.integrate(w**(r-2)*R*R,w)
    X=S.Matrix([-(r+k)*w**(r-1)*z+w**(r-1)*R+2*k*I,
      -(r+k-1)*w**r*z+w**r*R+2*(k-1)*J,
      -S.Rational(r,2)*w**(r-1)*z*z+w**(r-1)*R*z-K0])
    return X.applyfunc(S.cancel),I,J,K0
def family(D):
    if not isinstance(D,int) or D<3:raise ValueError("D must be an integer >=3")
    a,k,g,u,v,w,z,Y1,Y2,Y3=S.symbols("alpha kappa gamma u v w z Y1 Y2 Y3")
    cd=S.Rational(D*D*(D-1)**2,(2*D+1)*(D+1)**2*(D+2)**2)
    bd=S.Rational(D*(D-1),3*(D+1)*(D+2)*(D+3))
    minimal=cd*a*a+bd*a+S.Rational(1,180)
    def red(e):
        e=S.Poly(S.expand(e),k).rem(S.Poly(k*k+k+1,k)).as_expr()
        return S.Poly(e,a).rem(S.Poly(minimal,a)).as_expr().expand()
    def inv_a(e):return S.invert(e,minimal,a)
    R=a*w**D+w*w
    X,I,J0,K0=surface(R,w,z,k);Delta=S.Matrix([1,w,z]);sweep=X+g*Delta
    Z=Y2-w*Y1+2*k*w*I-2*(k-1)*J0
    Gamma=w*Y1+(k+2)*Z-w*w*R-2*k*w*I
    P=S.expand(w**3*(Y3+K0)-w*(Y1-2*k*I)*Z-(k+1)*Z*Z)
    g0=-S.Rational(D*(D-1),(D+1)*(D+2))*a-S.Rational(1,6)
    z0=2*(k+D+1)*a/((D+1)*(D+2))+(k+3)/6
    gi=inv_a(g0)
    E=[S.cancel(sweep[i].subs({w:g*u,z:g*g*v},simultaneous=True)/g**e) for i,e in enumerate([1,2,3])]
    theta0=S.Matrix([g0,gi,red(z0*gi*gi)]);base=dict(zip([g,u,v],theta0))
    def at(e):return red(e.subs(base))
    J=S.Matrix([[at(S.diff(e,t)) for t in [g,u,v]] for e in E])
    L=J.adjugate().applyfunc(lambda e:red(e*g0))
    c1=L[:,0]
    h=red(-sum(S.diff(E[2],t,s)*c1[i]*c1[j] for i,t in enumerate([g,u,v]) for j,s in enumerate([g,u,v])).subs(base)/2)
    return dict(D=D,N=2*D+4,a=a,k=k,g=g,u=u,v=v,w=w,z=z,Y1=Y1,Y2=Y2,Y3=Y3,
       minimal=minimal,cd=cd,bd=bd,red=red,inv_a=inv_a,R=R,X=X,I=I,J0=J0,K0=K0,
       Delta=Delta,sweep=sweep,Z=Z,Gamma=Gamma,P=P,g0=g0,z0=z0,E=E,theta0=theta0,J=J,L=L,h=h)

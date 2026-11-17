"""Homogeneous signed sweeps and their independently curve-induced pair blocks."""
import sympy as S
def homogeneous_signed(ell,w,z,g,q=None):
    if not isinstance(ell,int) or ell<2:raise ValueError("ell>=2 integer")
    q=S.Symbol("q") if q is None else q
    b0=S.Rational(ell-1,ell);d0=S.Rational(ell-1,2*ell-1)
    sig=w**(2*ell-1)-z
    X=S.Matrix([sig,b0*w**(2*ell),d0*w**(4*ell-2)-sig**2/2])
    tail=S.div((ell-1)*q**(2*ell)-2*(2*ell-1)*q**ell+2*ell**2*q-(ell-1)*(2*ell-1),(q-1)**3,q)
    if tail[1]!=0:raise ArithmeticError("exact triple ramification factor")
    return dict(ell=ell,n=2*ell,N=4*ell,b0=b0,d0=d0,X=X,
      sweep=X+g*S.Matrix([1,w,z]),tail=S.expand(tail[0]),q=q,
      geometric_group="SL2(F3) on8" if ell==2 else f"2^{2*ell-1}:A{2*ell} on{4*ell}")
def block_curve(ell,u):
    c=-S.Rational((ell-1)**2,2*(2*ell-1)*ell**2)
    F=S.Rational(1,2)+u/ell+c*u*u
    return dict(F=F,c=c,phi=F**ell/u,degree=2*ell,
      critical=S.Rational(ell,ell-1),fiber_gcd_lcm=ell)
def centered_signed(ell,g,u,v,w=None,z=None):
    w=S.Symbol("w") if w is None else w;z=S.Symbol("z") if z is None else z
    data=homogeneous_signed(ell,w,z,g)
    old=data["X"].subs({w:1+w,z:w+z},simultaneous=True)
    b0,d0=data["b0"],data["d0"]
    X=S.Matrix([old[0]-1,old[1]-old[0]+1-b0,
       old[2]-old[1]+old[0]-d0-S.Rational(1,2)+b0]).applyfunc(S.expand)
    sweep=X+g*S.Matrix([1,w,z])
    E=S.Matrix([S.cancel(sweep[i].subs({w:g*u,z:g*g*v},simultaneous=True)/g**(i+1)) for i in range(3)])
    return dict(X=X,E=E,original=data)
def marked_signed_data(ell):
    g,u,v,r=S.symbols("gamma u v r")
    data=centered_signed(ell,g,u,v)
    original=data["original"];q=original["q"]
    minimum=S.Poly(original["tail"].subs(q,r*r),r,domain=S.QQ)
    def red(expr):
        num,den=S.fraction(S.cancel(expr))
        N=S.Poly(num,r,domain=S.QQ).rem(minimum);D=S.Poly(den,r,domain=S.QQ).rem(minimum)
        return (N*S.invert(D,minimum)).rem(minimum).as_expr()
    g0=red(original["b0"]*(1-r**(2*ell))/r)
    theta0=S.Matrix([g0,red((r-1)/g0),red((r**(2*ell-1)-r+g0)/g0**2)])
    E=data["E"];at=lambda expr:red(expr.subs(dict(zip([g,u,v],theta0)),simultaneous=True))
    J=E.jacobian([g,u,v]).applyfunc(at)
    det=red(J.det());L=J.adjugate().applyfunc(lambda expr:red(expr/det))
    Hess=S.hessian(E[2],[g,u,v]).applyfunc(at)
    f20=red(-(L[:,0].T*Hess*L[:,0])[0]/2)
    return dict(r=r,minimum=minimum,red=red,E=E,J=J,L=L,theta0=theta0,
       f20=f20,determinant=det,at=at,original=original)

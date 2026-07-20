"""General 3D weighted-lift validation controls."""
import sympy as sp

x,y,z,w=sp.symbols("x y z w")

def construct(seed,c=sp.Integer(1),b=sp.Integer(1)):
    seed=sp.expand(seed);c=sp.sympify(c);b=sp.sympify(b)
    assert seed.subs(w,0)==0 and seed.subs(w,1)==-c and sp.integrate(seed,(w,0,1))==0
    q=sp.integrate(w*sp.diff(seed,w)/c,w)
    kappa=sp.diff(seed,w).subs(w,1)/c
    assert kappa!=-2
    a=-(1+kappa)/(2+kappa)
    vv=x*y;ss=x**2*z;uu=1+vv;gamma=1+a*vv+b*ss;W=uu*gamma
    alpha=sp.cancel(uu+q.subs(w,W)/gamma**2)
    beta=sp.cancel(c+seed.subs(w,W)/gamma)
    F=tuple(sp.cancel(f) for f in (alpha/x**2,beta/x,x*gamma))
    assert all(sp.denom(f)==1 for f in F)
    assert sp.factor(sp.Matrix(F).jacobian((x,y,z)).det())==b*c
    return F


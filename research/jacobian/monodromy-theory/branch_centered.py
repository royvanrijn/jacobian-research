"""Branch-centered weighted completion and an explicit binary degree8 map."""
import sympy as S

def binary_sweeps(w,z,g):
    sigma=w**3-z
    original=S.Matrix([sigma,w**4/2,w**6/3-sigma**2/2])
    old=original.subs({w:w+1,z:w+z},simultaneous=True)
    centered=S.Matrix([old[0]-1,old[1]-old[0]+S.Rational(1,2),
       old[2]-old[1]+old[0]-S.Rational(1,3)]).applyfunc(S.expand)
    return dict(original_X=original,centered_X=centered,
       original=original+g*S.Matrix([1,w,z]),
       centered=centered+g*S.Matrix([1,w,z]))

def binary_stage(g,u,v):
    return S.Matrix([
      1+2*u+g*(3*u*u-v)+g*g*u**3,
      u+v+g*u**3+g*g*u**4/2,
      v+2*u*v-S.Rational(4,3)*u**3+
        g*(-2*u**4+3*u*u*v-v*v/2)+
        g*g*(-u**5+u**3*v)-g**3*u**6/6])

def binary_constants(r):
    theta0=S.Matrix([4*r/3,(r+3)/4,r/2])
    L=S.Matrix([[-8,-32*r/3,S.Rational(56,3)],
       [S.Rational(3,2)-r,2*r+4,5*r/2-S.Rational(7,2)],
       [r/4+S.Rational(15,4),5*r,-10]])
    return dict(theta0=theta0,L=L,f20=10*r-18)

def expanded_binary_map():
    """Sparse exact polynomial coordinates over Q(sqrt(-3)); no rational opens."""
    from sympy.polys.rings import ring
    K=S.QQ.algebraic_field(S.sqrt(-3))
    R,x,y,z,t=ring("x,y,z,t",K)
    rr=K.from_sympy(S.sqrt(-3))
    rat=lambda a,b=1:K(S.QQ(a,b))
    L=[[-rat(8),-rat(32,3)*rr,rat(56,3)],
       [rat(3,2)-rr,2*rr+rat(4),rat(5,2)*rr-rat(7,2)],
       [rr/rat(4)+rat(15,4),5*rr,-rat(10)]]
    theta0=[rat(4,3)*rr,(rr+rat(3))/rat(4),rr/rat(2)]
    f=10*rr-rat(18)
    vs=[x*y,x*x*z,x**3*t+f*x*x*y*y]
    g,u,v=[theta0[row]+sum(L[row][col]*vs[col] for col in range(3)) for row in range(3)]
    E1=1+2*u+g*(3*u*u-v)+g*g*u**3
    E2=u+v+g*u**3+rat(1,2)*g*g*u**4
    E3=v+2*u*v-rat(4,3)*u**3+g*(-2*u**4+3*u*u*v-rat(1,2)*v*v)+g*g*(-u**5+u**3*v)-rat(1,6)*g**3*u**6
    F=[x*g,E1.exquo(x),E2.exquo(x*x),E3.exquo(x**3)]
    return dict(K=K,ring=R,variables=[x,y,z,t],r=rr,gamma0=theta0[0],
       gamma=g,F=F,numerators=[E1,E2,E3])

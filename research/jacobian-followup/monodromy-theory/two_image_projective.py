"""Two quadratic normal forms and exact coupled companion curl certificate."""
import sympy as S
from primitive_multiplier import feedback_form
from rank_one_projective import homogenize
def nilpotent_normalizers(x,y,z,a,b,c,d):
    K=c/(2*d);B=(b-c*c/(4*d))/a
    square=S.Matrix([[1,0,0],[K/a,1/a,0],[0,-B/(d*a*a),1/(d*a*a)]])
    mixed=S.Matrix([[1,0,0],[0,1/a,0],[0,-b/(a*a*c),1/(a*c)]])
    return square,mixed
def conjugate(H,Smat,xs):
    inv=Smat.inv();xx=S.Matrix(xs)
    return (Smat*H.subs(dict(zip(xs,inv*xx)),simultaneous=True)).applyfunc(S.cancel)
def companion_certificate():
    u,x,y,z,lam,beta,gamma,tau,sig,delta=S.symbols("u x y z lam beta gamma tau sig delta")
    ps=[u,x,y,z];xs=[x,y,z];v=S.Matrix([x*x+y,z-sig*x+tau*y+beta*x*y+gamma*y*y,delta*x])
    V=homogenize(v,u,xs,2);f=feedback_form(V,ps,lam*u)
    D=S.cancel(f["detC"].subs(u,1)/delta);eta=f["eta"].subs(u,1)
    Nx=S.cancel(D*eta[1]);Ny=S.cancel(D*eta[2])
    Dp=S.Poly(D,*xs);Nxp=S.Poly(Nx,*xs);Nyp=S.Poly(Ny,*xs)
    # Work with polynomial numerators; no generic rational derivative simplification.
    R=Dp*(Nyp.diff(x)-Nxp.diff(y))-Nyp*Dp.diff(x)+Nxp*Dp.diff(y)
    return dict(ps=ps,parameters=[lam,beta,gamma,tau,sig,delta],v=v,V=V,feedback=f,
                D=D,eta=eta,Nx=Nx,Ny=Ny,curl_numerator=R)

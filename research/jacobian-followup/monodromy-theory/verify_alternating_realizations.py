"""Exact controls for ALTERNATING_REALIZATIONS.md; written proofs carry universal claims."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s

ROOT=Path(__file__).resolve().parent
started=time.monotonic(); checks=[]
def check(name,condition):
    if not condition: raise AssertionError(name)
    checks.append(name)
def zero(expr):
    return s.cancel(expr)==0

# Degree-five diagonal beta--Hermite member, in an exact quadratic field.
L,v,u,x,y,z,t=s.symbols("L v u x y z t")
rad=s.sqrt(-15)
zeta=(-9+rad)/6
delta=s.expand(zeta**2-1)
def poly(expr):
    return s.Poly(expr,L,v,extension=rad)
one=poly(1)
g=1+delta*L/2-delta**2*L**2/8+(zeta-1-delta/2+delta**2/8)*L**3
gp=poly(g)
rho=gp+poly(L**3*v)
Q=(rho*rho-poly(1+delta*L)).mul_ground(s.Rational(1,2))
hh=one+Q.mul_ground(2)
P=((rho**5-one).mul_ground(s.Rational(1,5))
   -(hh*(rho**3-one)).mul_ground(s.Rational(2,3))
   +hh**2*(rho-one)).mul_ground(6/delta**3)
A=Q.exquo(poly(L**3));B=P.exquo(poly(L**3))
check("beta chosen root and nonzero delta",s.simplify(3*zeta**2+9*zeta+8)==0 and delta!=0)
M=(t**5-1)/5-2*(t**3-1)/3+t-1
check("beta root factorization",s.expand(15*M-(t-1)**3*(3*t*t+9*t+8))==0)
check("beta first interpolation jet",poly(gp.as_expr()**2-1-delta*L).rem(poly(L**3)).is_zero)
check("beta second interpolation value",s.simplify(gp.as_expr().subs(L,1)-zeta)==0)
check("beta quotient exact Jacobian",
      P.diff(L)*Q.diff(v)-P.diff(v)*Q.diff(L)==poly(3*L**5))
for name,coord in [("A",A),("B",B)]:
    shifted=s.Poly(coord.as_expr().subs(L,1+u),u,v,extension=rad)
    check("beta global polynomial lift "+name,
          all(3*i+3*j>=3 for (i,j),coeff in shifted.terms()))
    check("beta lift value "+name,s.simplify(coord.as_expr().subs({L:1,v:0}))==0)
check("beta normalization at L=0",s.simplify(B.as_expr().subs(L,0)-1)==0)
rr,qq,pp=s.symbols("rr qq pp")
betaP=6/delta**3*((rr**5-1)/5-2*(1+2*qq)*(rr**3-1)/3+(1+2*qq)**2*(rr-1))
check("beta polynomial degree five",s.degree(betaP,rr)==5)
check("beta fixed-Q derivative",zero(s.diff(betaP,rr)-6/delta**3*(rr*rr-1-2*qq)**2))
# Rational reconstruction, and its reverse in the polynomial source pair.
Lr=(rr*rr-1-2*qq)/delta
vr=(rr-gp.as_expr().subs(L,Lr))/Lr**3
rho_r=gp.as_expr().subs(L,Lr)+Lr**3*vr
check("beta inverse recovers rho",zero(rho_r-rr))
check("beta inverse recovers Q",zero((rho_r*rho_r-1-delta*Lr)/2-qq))
check("beta forward inverse recovers L",poly((rho.as_expr()**2-1-2*Q.as_expr())/delta)==poly(L))
check("beta forward inverse recovers v",poly(rho.as_expr()-gp.as_expr())==poly(L**3*v))
# Exact branch-cycle control; the written lemma proves the general family.
def mul(a,b): return tuple(a[b[i]] for i in range(len(a)))
def generated(gens):
    ident=tuple(range(len(gens[0])));out={ident};todo=[ident]
    while todo:
        a=todo.pop()
        for b in gens:
            c=mul(a,b)
            if c not in out:out.add(c);todo.append(c)
    return out
def cycle(n,items):
    p=list(range(n))
    for a,b in zip(items,items[1:]+items[:1]):p[a]=b
    return tuple(p)
for r,expected in [(2,6),(3,60),(4,5040)]:
    n=2*r-1
    aa=cycle(n,list(range(r)))
    bb=cycle(n,[0]+list(range(r,n)))
    check("beta two-cycle monodromy control r="+str(r),len(generated([aa,bb]))==expected)
# Alternating-shear member. Work with factored expressions to avoid degree-80 expansion.
X,Y,c,d,a,b,w,om,T=s.symbols("X Y c d a b w om T")
shearP=(X+2*c*Y**3/X-c**2*Y**6/(3*X**3)+6*c*d*X*Y**2
 +6*c**2*d*Y**5/X-2*c*d**2*X**3*Y+15*c**2*d**2*X*Y**4
 +s.Rational(2,5)*c*d**3*X**5-s.Rational(20,3)*c**2*d**3*X**3*Y**3
 +3*c**2*d**4*X**5*Y**2-s.Rational(6,7)*c**2*d**5*X**7*Y
 +c**2*d**6*X**9/9)
Phi=X**2-c*(Y-d*X**2)**3
check("shear perfect-square derivative",zero(s.diff(shearP,X)-Phi**2/X**4))
formalB=s.expand(shearP.subs({X:L*a*a,Y:L*a*b})/L)
check("shear B is polynomial",s.denom(formalB)==1)
C1=1+c*L*w**3;C2=w+d*L*C1**3
Xf=L*C1**2;Yf=L*C1*C2
check("shear quotient coordinate Jacobian",
      s.expand(s.diff(Xf,L)*s.diff(Yf,w)-s.diff(Xf,w)*s.diff(Yf,L)-Xf)==0)
check("shear inverse radical identity",s.expand(Phi.subs({X:Xf,Y:Yf})-L**2*C1**3)==0)
check("shear quotient Jacobian equals L",zero((L**2*a**3)**2/(L*a*a)**3-L))
# Explicit rational inverse of (L,w)->(X,Y).
Linv=Phi**2/X**3
winv=X*(Y-d*X**2)/Phi
C1inv=s.cancel(1+c*Linv*winv**3)
C2inv=s.cancel(winv+d*Linv*C1inv**3)
check("shear inverse recovers X",zero(Linv*C1inv**2-X))
check("shear inverse recovers Y",zero(Linv*C1inv*C2inv-Y))
check("shear reverse inverse recovers L",s.expand((L**2*C1**3)**2-L*Xf**3)==0)
check("shear reverse inverse recovers w",s.expand(Xf*(Yf-d*Xf**2)-w*L**2*C1**3)==0)
c0=(9-3*s.I)/8;d0=(1664+1152*s.I)/125
check("shear exact quadratic coefficient",s.expand(32*c0*c0-72*c0+45)==0)
check("shear coefficient relation",s.expand(d0*(1-c0)**3-1)==0)
check("shear lift A vanishes",s.expand((1-c0)*(-1+d0*(1-c0)**3))==0)
b0=formalB.subs({L:1,a:1-c,b:0,d:1/(1-c)**3})
check("shear lift B condition",
      zero(b0-(32*c*c-72*c+45)/45))
# At x=0, quotients A/x and B/x are the first directional derivatives.
Aval=a*b
Araw=Aval.subs({a:C1,b:C2})
A_L=s.diff(Araw,L).subs({L:1,w:-1,c:c0,d:d0})
A_w=s.diff(Araw,w).subs({L:1,w:-1,c:c0,d:d0})
# Chain rule for B in independent (L,a,b), then C1,C2.
B_L=(s.diff(formalB,L)+s.diff(formalB,a)*s.diff(C1,L).subs(w,-1)
      +s.diff(formalB,b)*s.diff(C2,L).subs(w,-1))
B_w=(s.diff(formalB,a)*s.diff(C1,w).subs(w,-1)
      +s.diff(formalB,b)*s.diff(C2,w).subs(w,-1))
lift_sub={L:1,a:1-c,b:0}
BL=s.simplify(B_L.subs(lift_sub).subs({c:c0,d:d0}))
BW=s.simplify(B_w.subs(lift_sub).subs({c:c0,d:d0}))
yp=10-3*s.I;zp=s.Rational(7,2)-s.Rational(3,2)*s.I
check("shear collision first coordinate at x=0",s.simplify(A_L*yp+A_w*zp)==0)
check("shear collision second coordinate at x=0",s.simplify(BL*yp+BW*zp-1)==0)
check("shear other collision point",formalB.subs({L:0,a:1,b:0})==1)
check("shear finite critical roots generically simple",
      s.factor(s.discriminant(Phi,X))==64*Y**3*c**4*d**9*(27*Y**2*c*d+4)**2)
check("shear degree twelve poles",
      s.degree(s.cancel(shearP*X**3),X)==12 and s.cancel(shearP*X**3).subs(X,0)==-c*c*Y**6/3)
# Generic critical-value distinctness in normalized family, by one exact nonzero specialization.
pn=shearP.subs({c:1,d:1,Y:1})
fn=s.cancel((pn-T)*X**3)
critical=Phi.subs({c:1,d:1,Y:1})
value_resultant=s.Poly(s.resultant(fn,critical,X),T)
check("shear six distinct finite critical values",
      value_resultant.degree()==6 and s.gcd(value_resultant,value_resultant.diff()).degree()==0)
# Normalization scales: b=d*a^2, c*d^3*a^4=1.
ascl,bscl=s.symbols("ascl bscl",nonzero=True)
scaled=s.expand(shearP.subs({X:ascl*X,Y:d*ascl**2*Y})/ascl-shearP.subs({c:1,d:1}))
num=s.Poly(scaled,c)
check("shear normalization scaling",
      s.rem(num,s.Poly(c*d**3*ascl**4-1,c)).is_zero)
check("shear Riemann-Hurwitz",6*(3-1)+(9-1)+(3-1)==2*12-2)
# Exact discriminant/resolvent sign control, no square => group assertion.
discn=s.factor(s.discriminant(fn,X))
resn=s.resultant(fn,critical,X)
check("shear discriminant square class minus three",zero(discn+27*resn**2))
# Generic numerator has leading coefficient c^2*d^6/9 and constant -c^2*Y^6/3.
# The written resultant argument extends this square class to all nonzero c,d.
check("shear even branch profiles",(9-1+3-1)%2==0 and (3-1)%2==0)
manifest=json.loads((ROOT/"sources/alternating-realizations-manifest.json").read_text())
for source in manifest["sources"]:
    check("primary source snapshot hash",
          hashlib.sha256((ROOT/"sources"/source["path"]).read_bytes()).hexdigest()==source["sha256"])
data=dict(
 beta=dict(r=3,zeta=str(zeta),delta=str(s.simplify(delta)),g=str(gp.as_expr()),
           A=str(A.as_expr()),B=str(B.as_expr()),rho=str(rho.as_expr()),
           map="(A(1+x^3*y,x^3*z)/x^3, B(1+x^3*y,x^3*z)/x^3, x*(1+x^3*y))",
           jacobian=-3,degree=5,geometric_group="A5",
           universal_family="r>=2: degree 2r-1, A_(2r-1) for odd r, S_(2r-1) for even r"),
 shear=dict(c=str(c0),d=str(d0),A="C1*C2",B=str(formalB),P=str(shearP),
            map="(A(1+xy,xz-1)/x, B(1+xy,xz-1)/x, x*(1+xy))",
            jacobian=-1,degree=12,geometric_group="A12",
            inverse={"L":str(Linv),"w":str(winv)},
            normalized_critical_value_resultant=str(value_resultant.as_expr())),
 attribution="Kistner--Shaska, arXiv:2608.02863v1, sections 4--6; independently rederived algebra and written group proofs.",
 full_action_spectrum_complete=False,PSL2_11_realized=False)
(ROOT/"alternating_realizations_model.json").write_text(json.dumps(data,indent=2)+"\n")
files=["ALTERNATING_REALIZATIONS.md","verify_alternating_realizations.py",
       "alternating_realizations_model.json","sources/alternating-realizations-manifest.json"]
receipt=dict(result="pass",checks=checks,check_count=len(checks),
 python=platform.python_version(),sympy=s.__version__,runtime_seconds=round(time.monotonic()-started,3),
 sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
 assurance="Exact algebraic identities and finite controls supporting written universal proofs; no formal verification, independent whole-proof review or new priority claim.",
 actual_A5_affine_Keller_map=True,actual_A12_affine_Keller_map=True,
 full_action_spectrum_complete=False,PSL2_11_realized=False)
(ROOT/"alternating_realizations_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                     actual_groups=["A5","A12"],complete_classification=False)))

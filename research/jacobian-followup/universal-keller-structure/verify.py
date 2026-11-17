"""Small exact controls for written uniform theorems; not formal verification."""
from pathlib import Path
import hashlib
import json
import math
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
started=time.monotonic()
checks=[]
fibre_cases=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)

S,X,Y,r,B,C=s.symbols("S X Y r B C")
# Proof ingredients, including p|N and the small prime p=3.
for p,n in [(0,4),(3,3),(3,6),(3,9),(5,5),(5,10),(7,9)]:
    domain=s.QQ if not p else s.GF(p)
    h=S+S**3+S**n
    hp=s.diff(h,S)
    a=s.cancel((h.subs(S,X)-h.subs(S,Y))/(X-Y))
    rem=s.Poly(s.expand(a.subs(Y,-X)),X,domain=domain)
    check(f"char{p} degree{n} pair-incidence primitive coefficient",rem.nth(0)!=0)
    bfun=-hp.subs(S,r)/(2*r)
    cfun=h.subs(S,r)-r*hp.subs(S,r)/2
    identity=s.cancel(s.diff(cfun,r)-r**2*s.diff(bfun,r))
    check(f"char{p} degree{n} critical envelope identity",identity==0)
    numerator=s.Poly(s.expand(hp.subs(S,r)-r*s.diff(hp.subs(S,r),r)),r,domain=domain)
    check(f"char{p} degree{n} critical parameter separable",numerator.nth(0)!=0)
    check(f"char{p} degree{n} odd simple pole at critical origin",
          s.Poly(hp.subs(S,r),r,domain=domain).nth(0)!=0)
check("even-only countercontrol factors through S squared",
      s.expand((X**2)**2+B*X**2-C-(X**4+B*X**2-C))==0)
linearized=lambda a:a**8+B*a**2+a
check("characteristic-two linearized countercontrol",
      s.Poly(s.expand(linearized(X+Y)-linearized(X)-linearized(Y)),X,Y,B,modulus=2).is_zero)
check("degree-eight affine group is smaller than symmetric group",
      8*7*6*4==1344 and 1344<math.factorial(8))

class Ring:
    def __init__(self,P):
        self.P=P.monic();self.p=P.domain.characteristic();self.domain=P.domain
        self.one=self.poly(1)
    def poly(self,v):return s.Poly(v,S,domain=self.domain)
    def mod(self,v):
        if not isinstance(v,s.Poly):v=self.poly(v)
        return v.rem(self.P)
    def mul(self,a,b):return self.mod(a*b)
    def power(self,a,n):return self.mod(a**n)
    def inv(self,a):return s.invert(a,self.P)

def source_in_ring(R,xx,yy,zz,profiles):
    # Direct polynomial landing formula for the complete class in (7).
    P=R.poly
    m=R.mul
    pw=R.power
    tt=R.mod(R.one+m(xx,yy))
    qq=R.mod(m(pw(tt,2),zz)+m(pw(yy,2),R.one+P(3)*tt))
    pp=m(tt,qq)
    BB=R.mod(yy+P(3)*m(xx,qq))
    CC=R.mod(m(xx,P(5)-P(3)*tt)-m(pw(xx,3),zz))
    for j,(g0,g1,u,v) in profiles.items():
        regular=R.mod(P(g0)+P(g1)*pp)
        positive=m(regular,m(pw(tt,2),m(pw(xx,j-2),pw(qq,j))))
        simple=m(P(u),m(tt,m(pw(xx,j-2),pw(qq,j-1))))
        double=m(P(v),m(pw(xx,j-2),pw(qq,j-2)))
        BB=R.mod(BB+P(j)*(positive+simple+double))
        CC=R.mod(CC-P(j-2)*m(regular,pw(m(xx,qq),j)))
        if u or v:assert (j-2)%R.p==0
    return [pp,R.mod(P(-pow(2,-1,R.p))*BB),CC],tt,qq

def boundary_case(p,n,profiles,b,c,label):
    P=lambda v:s.Poly(v,S,modulus=p)
    phi=P(S)
    for j,(_,_,_,v) in profiles.items():
        phi-=P(v)*P((-S)**(j-2))
    phi=P(phi.as_expr())
    check(label+" boundary derivative one",phi.diff()==P(1))
    d=phi.degree()
    part=phi-P(b)
    if b%p==0:part=part.exquo(P(S))
    ranks=[]
    if part.degree()>0:
        R=Ring(part)
        yy=R.mod(R.poly(S))
        xx=-R.inv(yy)
        extra=R.poly(0)
        for j,(g0,_,_,_) in profiles.items():
            extra=R.mod(extra+R.poly((j-2)*g0)*R.power(-yy,j))
        zz=R.mul(R.power(yy,3),R.mod(R.poly(c)+R.poly(5)*R.inv(yy)+extra))
        out,tt,qq=source_in_ring(R,xx,yy,zz,profiles)
        check(label+" t-component whole algebra reconstructs",
              all(v==R.poly(a) for v,a in zip(out,[0,b%p,c%p]))
              and tt.is_zero and R.mul(qq,R.inv(qq))==R.one)
        check(label+" t-component recovers its generator",yy==R.mod(R.poly(S)))
        ranks.append(part.degree())
    Q=P(S+b*S**2-c*pow(2,-1,p))
    qrank=0
    if s.gcd(Q,Q.diff()).degree()==0:
        R=Ring(Q)
        D=R.mod(Q.diff())
        xx=R.mul(R.poly(S),R.inv(D))
        yy=R.poly(-2*b)
        zz=-R.mul(R.power(yy,2),R.mod(R.power(D,2)+R.poly(3)*D))
        out,tt,qq=source_in_ring(R,xx,yy,zz,profiles)
        check(label+" q-component whole algebra reconstructs",
              all(v==R.poly(a) for v,a in zip(out,[0,b%p,c%p]))
              and qq.is_zero and R.mul(tt,R.inv(tt))==R.one)
        check(label+" q-component recovers its generator",
              R.mul(xx,R.inv(tt))==R.mod(R.poly(S)))
        qrank=Q.degree()
    elif Q.degree()==2:
        check(label+" repeated quadratic component entirely deleted",
              s.gcd(Q,Q.diff()).degree()==1)
    expected=d if b%p==0 or (1+2*b*c)%p==0 else d+2
    check(label+" exact total boundary rank",sum(ranks)+qrank==expected)
    fibre_cases.append({"label":label,"characteristic":p,"generic_degree":n,
                        "boundary_degree":expected,"t_rank":sum(ranks),"q_rank":qrank})

for p in [3,5,7]:
    n=p+2
    profiles={n:(0,0,0,-1)}
    boundary_case(p,n,profiles,1,0,f"F{p} full-rank boundary")
    boundary_case(p,n,profiles,1,(-pow(2,-1,p))%p,f"F{p} field boundary")
    boundary_case(p,n,profiles,0,2,f"F{p} zero-b split boundary")
    P=s.Poly(S-S**p-1,S,modulus=p)
    check(f"F{p} field boundary Artin-Schreier irreducible",P.is_irreducible)
    check(f"F{p} b-zero boundary splits into p rational points",
          all(s.Poly(S-S**p,S,modulus=p).eval(a)==0 for a in range(p)))
boundary_case(3,5,{5:(0,0,1,0)},1,0,"F3 simple pole keeps old boundary")
boundary_case(3,8,{4:(1,1,0,0),5:(0,0,1,1),8:(0,0,2,2)},
              1,0,"F3 mixed double poles and regular profile")

# Complete pi-nonzero algebras, not just rational source points.
for p,n,profiles in [
    (3,5,{5:(0,0,0,-1)}),
    (5,7,{4:(1,1,0,0),7:(0,0,2,-1)}),
    (3,6,{5:(0,0,0,1),6:(1,0,0,0)})]:
    pi=1;b=1
    P=lambda v:s.Poly(v,S,modulus=p)
    coeff={j:(g0+g1*pi+u*pow(pi,-1,p)+v*pow(pi,-2,p))%p
           for j,(g0,g1,u,v) in profiles.items()}
    assert coeff[n]!=0, "chosen finite control loses its top degree"
    h=P(S+b*S**2+pi*S**3)
    for j,a in coeff.items():h+=P(a*pi**j*S**j)
    E=None
    for c in range(p):
        candidate=h-P(c*pow(2,-1,p))
        if s.gcd(candidate,candidate.diff()).degree()==0:
            E=candidate;break
    assert E is not None
    R=Ring(E)
    D=R.mod(E.diff())
    yy=R.poly(-2*b)-R.poly(3*pi)*R.poly(S)
    for j,a in coeff.items():yy-=R.poly(j*a*pi**j*S**(j-2))
    yy=R.mod(yy)
    xx=R.mul(R.poly(S),R.inv(D))
    zz=R.mod(R.poly(pi)*R.power(D,3)
             -R.mul(R.power(yy,2),R.mod(R.power(D,2)+R.poly(3)*D)))
    out,tt,qq=source_in_ring(R,xx,yy,zz,profiles)
    label=f"F{p} rank{n} pi-nonzero"
    check(label+" complete algebra reconstructs all targets",
          E.degree()==n and all(v==R.poly(a) for v,a in zip(out,[pi,b,c])))
    check(label+" inverse recovers root",R.mul(xx,R.inv(tt))==R.mod(R.poly(S)))
    fibre_cases.append({"label":label,"characteristic":p,"generic_degree":n,
                        "pi":pi,"b":b,"c":c,"rank":E.degree()})

# Normal-form controls: a nonconstant unimodular normal block, including
# a nonlinear invariant component which cannot create moving fibre points.
u1,u2,tau,target=s.symbols("u1 u2 tau target")
A=s.Matrix([[1,tau],[0,1]])
ff=s.Matrix([u1+tau*u2,u2,tau+u2**2])
check("equivariant affine-moving control is Keller",ff.jacobian([u1,u2,tau]).det()==1)
check("nonconstant normal matrix has polynomial inverse",
      A.det()==1 and A.inv()==s.Matrix([[1,-tau],[0,1]]))
check("normal-form matrix control in tame characteristic three",
      s.Poly(A.det()-1,tau,modulus=3).is_zero)
gauge={u1:u1-tau*u2}
normalized=ff.subs(gauge,simultaneous=True)
check("source gauge normalizes moving first jet",
      normalized==s.Matrix([u1,u2,tau+u2**2]))
check("invariant nonlinear component does not change fixed-target fibre",
      normalized.subs({u1:0,u2:0})==s.Matrix([0,0,tau]))
# The nonlinear cubic compiler permits a nontrivial C2 action on its fibre:
# this is a countercontrol against overextending the affine-moving exclusion.
x,y,z=s.symbols("x y z")
tt=1+x*y;qq=tt**2*z+y**2*(1+3*tt)
base=s.Matrix([tt*qq,-(y+3*x*qq)/2,x*(5-3*tt)-x**3*z])
flipped=base.subs({x:-x,y:-y},simultaneous=True)
check("cubic compiler C2 covariance",
      all(s.expand(v)==0 for v in flipped-s.Matrix([base[0],-base[1],-base[2]])))
points=[(0,0,-1),(-s.Rational(1,2),3,26),(s.Rational(1,2),-3,26)]
check("nonlinear moving cubic has fixed and exchanged fibre points",
      all(base.subs(dict(zip([x,y,z],pt)))==s.Matrix([-1,0,0]) for pt in points)
      and points[1][0]==-points[2][0] and points[1][1]==-points[2][1])
normal=s.Matrix([base[1],base[2]]).jacobian([x,y]).subs({x:0,y:0})
check("nonlinear cubic normal block remains globally unimodular",
      normal==s.Matrix([[-3*z/2,-s.Rational(1,2)],[2,0]]) and normal.det()==1)

manifest=json.loads((HERE/"sources.json").read_text())
for pin in manifest["local_sources"]:
    check("source pin "+pin["path"],
          hashlib.sha256((RESEARCH/pin["path"]).read_bytes()).hexdigest()==pin["sha256"])
receipt={
    "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
    "complete_fibre_controls":fibre_cases,
    "assurance":"Written uniform proofs with bounded exact controls, not formal verification or independent whole-proof replay. No finite sample is used to infer a generic Galois group or an all-degree statement.",
    "python":platform.python_version(),"sympy":s.__version__,
    "runtime_seconds":round(time.monotonic()-started,3),
    "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [HERE/"PROOF.md",HERE/"LITERATURE_AND_GAPS.md",Path(__file__),HERE/"sources.json"]}
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"fibre_controls":len(fibre_cases),
                  "seconds":receipt["runtime_seconds"]}))

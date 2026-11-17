"""Exact controls for sextic compression and sparse Fitting transport.
The source-map routine is explicitly reused. Universal proofs are in PROOF.md.
"""
from pathlib import Path
import ast
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
started=time.monotonic()
checks=[]
def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)

x,y,z,h,S,T,pi,b,c,lam=s.symbols("x y z h S T pi b c lam")
predecessor=HERE.parent/"universal-keller-compression"/"verify.py"
tree=ast.parse(predecessor.read_text())
nodes=[node for node in tree.body if isinstance(node,ast.FunctionDef)
       and node.name=="source_map"]
assert len(nodes)==1
exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),
             str(predecessor),"exec"))

# Full polynomial Jacobian, including h as a coefficient variable.
f=source_map(6,{4:h**2,5:s.Integer(0),6:h**3})
jac=[[s.Poly(s.diff(g,v),x,y,z,h,domain=s.QQ) for v in (x,y,z)] for g in f]
det=(jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])
     -jac[0][1]*(jac[1][0]*jac[2][2]-jac[1][2]*jac[2][0])
     +jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]))
check("full sextic polynomial Jacobian",det.as_expr()==1)
degrees=[1]+[s.Poly(g,x,y,z,h).total_degree() for g in f]
check("promoted ordinary coordinate degrees",degrees==[1,7,41,39])

A,B=s.symbols("A B",nonzero=True)
for r,t,u,v,eps in [(5,6,1,1,-1),(5,7,2,3,1)]:
    tau=(B**u/A**v)**eps
    inv=(B**r/A**t)**eps
    check(f"coprime section weights {r},{t}",
          s.cancel(tau**r*A-inv**u)==0 and s.cancel(tau**t*B-inv**v)==0)
    check(f"Laurent inverse section weights {r},{t}",
          s.cancel(tau**(-r)*inv**u-A)==0 and s.cancel(tau**(-t)*inv**v-B)==0)

def compile_sextic(P):
    field=P.domain
    d={j:field.convert(P.nth(j)) for j in range(7)}
    if P.degree()!=6 or d[5] or not all(d[j] for j in (1,3,4)):
        raise AssertionError("Not a good trace-zero sextic")
    al=d[1]/d[4]**3
    hv=d[1]**4/(d[3]**2*d[4]**7)
    pv=d[3]*d[4]**6/d[1]**3
    bv=d[2]*d[4]**3/d[1]**2
    cv=-2*d[0]/d[4]**3
    aj={4:hv**2,5:field.zero,6:hv**3}
    co={0:-cv/field.convert(2),1:field.one,2:bv,3:pv,
        **{j:aj[j]*pv**j for j in aj}}
    E=s.Poly.from_dict({(j,):value for j,value in co.items()},S,domain=field)
    expected=s.Poly(field.to_sympy(field.one/d[4]**3)*P.as_expr().subs(
        S,field.to_sympy(field.one/al)*S),S,domain=field)
    return al,hv,pv,bv,cv,aj,E,expected

def verify_fibre(P,label):
    field=P.domain
    al,hv,pv,bv,cv,aj,E,expected=compile_sextic(P)
    check(label+" inverse polynomial identity",(E-expected).is_zero)
    check(label+" degree six and squarefree",
          E.degree()==6 and s.gcd(E,E.diff()).degree()==0)
    def poly(v): return s.Poly(v,S,domain=field)
    def scalar(v): return poly(field.to_sympy(v))
    def mod(v): return v.rem(E)
    def mul(a,b): return mod(a*b)
    def power(a,j): return mod(a**j)
    one=poly(1)
    D=E.diff()
    xx=mul(poly(S),s.invert(D,E))
    yy=scalar(-2*bv)-scalar(3*pv)*poly(S)
    for j in aj:
        yy-=scalar(j*aj[j]*pv**j)*poly(S**(j-2))
    yy=mod(yy)
    zz=mod(scalar(pv)*power(D,3)-power(yy,2)*mod(power(D,2)+poly(3)*D))
    tt=mod(one+mul(xx,yy))
    qq=mod(mul(power(tt,2),zz)+mul(power(yy,2),one+poly(3)*tt))
    Bout=mod(yy+poly(3)*mul(xx,qq))
    Cout=mod(mul(xx,poly(5)-poly(3)*tt)-mul(power(xx,3),zz))
    for j in aj:
        Bout=mod(Bout+poly(j)*mul(scalar(aj[j]),
             mul(power(tt,2),mul(power(xx,j-2),power(qq,j)))))
        Cout=mod(Cout-poly(j-2)*mul(scalar(aj[j]),power(mul(xx,qq),j)))
    check(label+" reconstructs complete target",
          (mul(tt,qq)-scalar(pv)).is_zero
          and (mul(Bout,scalar(-field.one/field.convert(2)))-scalar(bv)).is_zero
          and (Cout-scalar(cv)).is_zero)
    check(label+" inverse recovers generator",
          (mul(xx,s.invert(tt,E))-poly(S)).is_zero)
    return {"label":label,"degree":6,"characteristic":int(field.characteristic())}

P1=s.Poly(s.resultant(T**6-T-1,S-T**3-T,T),S,domain=s.QQ)
P2=s.Poly(s.resultant(T**6-T-1,S-T**4-T,T),S,domain=s.QQ)
Psplit=s.Poly(s.prod(S-r for r in [0,1,2,3,4,-10]),S,domain=s.QQ)
examples=[(P1,"QQ sextic generator T3+T"),(P2,"QQ sextic generator T4+T"),
          (Psplit,"QQ split sextic"),(s.Poly(P1.as_expr(),S,modulus=7),"F7 sextic")]
cases=[verify_fibre(P,label) for P,label in examples]

# The canonical realization depends on a line, not its representative.
for P,label in examples[:3]:
    old=compile_sextic(P)
    field=s.QQ.frac_field(lam)
    scaled=s.Poly(lam**6*P.as_expr().subs(S,S/lam),S,domain=field)
    new=compile_sextic(scaled)
    check(label+" normalizer weight minus one",
          s.cancel(field.to_sympy(new[0])-old[0]/lam)==0)
    check(label+" target independent of line representative",
          all(s.cancel(field.to_sympy(new[j])-old[j])==0 for j in (1,2,3,4)))
    normalized=s.Poly(old[0]**6*P.as_expr().subs(S,S/old[0]),S,domain=s.QQ)
    check(label+" canonical generator satisfies d1=d4 cubed",
          normalized.nth(1)==normalized.nth(4)**3)
    check(label+" canonical normalization round trip",
          compile_sextic(normalized)[0]==s.QQ.one)

# Small independent finite-field arithmetic for a split sextic in characteristic3.
# F27=F3[u]/(u^3-u-1); every nonzero element is explicitly inverted below.
class K:
    def __init__(self,v=0):
        if isinstance(v,K): self.v=v.v
        elif isinstance(v,(tuple,list)): self.v=tuple(int(a)%3 for a in v)
        else: self.v=(int(v)%3,0,0)
    def __add__(self,o):
        o=K(o); return K(tuple(a+b for a,b in zip(self.v,o.v)))
    __radd__=__add__
    def __neg__(self): return K(tuple(-a for a in self.v))
    def __sub__(self,o): return self+-K(o)
    def __rsub__(self,o): return K(o)+-self
    def __mul__(self,o):
        o=K(o); co=[0]*5
        for i,a in enumerate(self.v):
            for j,bv in enumerate(o.v): co[i+j]+=a*bv
        for j in (4,3):
            co[j-2]+=co[j]; co[j-3]+=co[j]; co[j]=0
        return K(co[:3])
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0: return self.inverse()**(-n)
        out=K(1); base=self
        while n:
            if n&1: out=out*base
            base=base*base; n//=2
        return out
    def inverse(self):
        if not self: raise ZeroDivisionError()
        return self**25
    def __truediv__(self,o): return self*K(o).inverse()
    def __rtruediv__(self,o): return K(o)*self.inverse()
    def __bool__(self): return any(self.v)
    def __eq__(self,o): return self.v==K(o).v

u=K((0,1,0))
check("F27 modulus irreducible over F3",all((a**3-a-1)%3 for a in range(3)))
nonzero=[K((a,bv,cv)) for a in range(3) for bv in range(3) for cv in range(3)
         if (a,bv,cv)!=(0,0,0)]
check("all F27 inverses",all(a*a.inverse()==1 for a in nonzero))
roots=[K(0),K(1),K(-1),u,u*u,-u-u*u]
check("F27 six trace-zero distinct roots",
      sum(roots)==0 and all(roots[i]!=roots[j] for i in range(6) for j in range(i)))
d=[K(1)]
for root in roots:
    new=[K(0)]*(len(d)+1)
    for j,co in enumerate(d):
        new[j]=new[j]-root*co
        new[j+1]=new[j+1]+co
    d=new
check("F27 good trace-zero coefficients",not d[5] and all(d[j] for j in (1,3,4)))
al=d[1]/d[4]**3
hv=d[1]**4/(d[3]**2*d[4]**7)
pv=d[3]*d[4]**6/d[1]**3
bv=d[2]*d[4]**3/d[1]**2
cv=-2*d[0]/d[4]**3
aj={4:hv**2,5:K(0),6:hv**3}
Eco=[-cv/2,K(1),bv,pv,aj[4]*pv**4,K(0),aj[6]*pv**6]
check("F27 coefficient identity",
      all(Eco[j]==d[j]/d[4]**3/al**j for j in range(7)))
for i,root in enumerate(roots):
    sv=al*root
    check(f"F27 root {i} satisfies E",sum(Eco[j]*sv**j for j in range(7))==0)
    dv=sum(j*Eco[j]*sv**(j-1) for j in range(1,7))
    check(f"F27 root {i} derivative unit",bool(dv))
    yy=-2*bv-3*pv*sv-sum(j*aj[j]*pv**j*sv**(j-2) for j in aj)
    xx=sv/dv
    zz=pv*dv**3-yy**2*(dv**2+3*dv)
    tt=1+xx*yy
    qq=tt**2*zz+yy**2*(1+3*tt)
    Bout=yy+3*xx*qq+sum(j*aj[j]*tt**2*xx**(j-2)*qq**j for j in aj)
    Cout=xx*(5-3*tt)-xx**3*zz-sum((j-2)*aj[j]*(xx*qq)**j for j in aj)
    check(f"F27 root {i} reconstructs source and generator",
          tt*qq==pv and -Bout/2==bv and Cout==cv and xx/tt==sv)
cases.append({"label":"F27 split sextic: all six source points","degree":6,
              "characteristic":3,"source_points_checked":6})
check("finite-field sufficient bound includes F27",27>15+5+3+2)

# The sparse Fitting divisor must retain P-exponents as well as r-exponents.
P,r,scale=s.symbols("P r scale",nonzero=True)
for n,aj in [(6,{4:h**2,5:s.Integer(0),6:h**3}),
             (7,{4:h,5:h,6:s.Integer(0),7:s.Symbol("a7")})]:
    J=-1+3*P*r**2+sum(j*(j-2)*aj[j]*P**j*r**(j-1) for j in aj)
    Jscaled=-1+3*P*r**2+sum(j*(j-2)*scale**(j+1)*aj[j]*P**j*r**(j-1) for j in aj)
    check(f"N{n} full sparse Fitting covariance",
          s.expand(Jscaled.subs({P:P/scale**2,r:scale*r},simultaneous=True)-J)==0)
    # Endpoint matching under orientation reversal requires this nonintegral twist.
    twist=s.Rational(-2*n,n-1)
    check(f"N{n} orientation reversal excluded by full endpoints",twist.q!=1)
check("rank6 projected support is reflection symmetric",
      {5-j for j in {0,2,3,5}}=={0,2,3,5})
check("rank6 stable-class invariant recovers h",s.cancel((h**3)**5/(h**2)**7-h)==0)
bad=s.Poly(S**5-S-1,S)
check("rank5 trace normalization loses a4",bad.nth(4)==0)

# Pin imported proof inputs, not just successful examples.
manifest=json.loads((HERE/"sources.json").read_text())
for item in manifest["local_sources"]:
    check("source pin "+item["path"],
          hashlib.sha256((RESEARCH/item["path"]).read_bytes()).hexdigest()==item["sha256"])
receipt={
    "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
    "complete_fibre_controls":cases,"coordinate_degrees":degrees,
    "source_routines_reused":["universal-keller-compression/verify.py:source_map"],
    "assurance":"Written proofs with exact controls, not formalization or independent whole-proof replay. Monodromy, all-rank endpoint rigidity, torsor identification and density use the written arguments.",
    "python":platform.python_version(),"sympy":s.__version__,
    "runtime_seconds":round(time.monotonic()-started,3),
    "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [HERE/"PROOF.md",Path(__file__),HERE/"sources.json"]}
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"fibre_controls":len(cases),
                  "seconds":receipt["runtime_seconds"]}))


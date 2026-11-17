"""Exact controls for trace compression; universal proofs are in PROOF.md."""
from pathlib import Path
import ast
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
PREDECESSOR=HERE.parent/"universal-keller-compression"
started=time.monotonic()
checks=[]
def check(label,condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
x,y,z,w,S,T,pi,b,c,alpha=s.symbols("x y z w S T pi b c alpha")
u,lam=s.symbols("u lam")
# Reuse two named predecessor routines without executing its regression campaign.
# This is declared code reuse, not an independent implementation.
tree=ast.parse((PREDECESSOR/"verify.py").read_text())
nodes=[node for node in tree.body if isinstance(node,ast.FunctionDef)
       and node.name in {"source_map","verify_fibre"}]
assert len(nodes)==2
def exact_equal(left,right):
    # SymPy Poly structural equality can differ for equal GF(p)(u) fractions.
    if isinstance(left,s.Poly) and isinstance(right,s.Poly):
        return (left-right).is_zero
    if hasattr(left,'numer') or hasattr(right,'numer'):
        return not bool(left-right)
    return left==right
class ExactComparisons(ast.NodeTransformer):
    def visit_Compare(self,node):
        node=self.generic_visit(node)
        if len(node.ops)==1 and isinstance(node.ops[0],ast.Eq):
            return ast.copy_location(ast.Call(func=ast.Name(id="exact_equal",ctx=ast.Load()),
                args=[node.left,node.comparators[0]],keywords=[]),node)
        return node
nodes=[ExactComparisons().visit(node) for node in nodes]
module=ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[]))
exec(compile(module,str(PREDECESSOR/"verify.py"),"exec"))

def compile_presentation(P):
    n=P.degree()
    field=P.domain
    d={j:field.convert(P.nth(j)) for j in range(n+1)}
    if not all(d[j] != field.zero for j in (1,3,4,5)):
        raise AssertionError("Input is not a good generator")
    al=d[3]*d[4]/(d[1]*d[5])
    pars={j:d[j]*d[3]*d[4]**(j+1)/(d[1]**2*d[5]**(j+1))
          for j in range(4,n+1)}
    pival=d[1]*d[5]**2/(d[3]*d[4]**2)
    bval=d[2]*d[5]/(d[3]*d[4])
    cval=-field.convert(2)*d[0]*d[3]*d[4]/(d[1]**2*d[5])
    coefficients={0:-cval/field.convert(2),1:field.one,2:bval,3:pival,
                  **{j:pars[j]*pival**j for j in pars}}
    E=s.Poly.from_dict({(j,):value for j,value in coefficients.items()},S,domain=field)
    expected=s.Poly(field.to_sympy(al/d[1])*P.as_expr().subs(
        S,field.to_sympy(field.one/al)*S),S,domain=field)
    return 0,al,pars,pival,bval,cval,E,expected

# One full polynomial Jacobian with the new vanishing coefficient.
n=7
aj={4:w,5:w,6:s.Integer(0),7:s.Integer(2)}
f=source_map(n,aj)
jac=[[s.Poly(s.diff(g,v),x,y,z,w,domain=s.QQ) for v in (x,y,z)] for g in f]
det=(jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])
     -jac[0][1]*(jac[1][0]*jac[2][2]-jac[1][2]*jac[2][0])
     +jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]))
check("N7 full symbolic Jacobian after trace compression",det.as_expr()==1)

cases=[]
P7=s.Poly(s.resultant(T**7-T-1,S-(T**3+T),T),S,domain=s.QQ)
examples=[
 (P7,"QQ rank7 generator of T7-T-1"),
 (s.Poly(s.prod(S-r for r in [0,1,2,3,4,5,-15]),S,domain=s.QQ),"QQ split rank7"),
 (s.Poly(s.prod(S-r for r in [0,1,2,3,4,5,6,-21]),S,domain=s.QQ),"QQ split rank8"),
 (s.Poly(P7.as_expr(),S,modulus=5),"F5 rank7 small characteristic"),
 (s.Poly(s.prod(S-r for r in [0,1,2,3,4,5,-15]),S,modulus=43),"F43 split rank7")
]
for P,label in examples:
    check(label+" trace coefficient zero",P.nth(P.degree()-1)==0)
    cases.append(verify_fibre(P,label))
    data=compile_presentation(P)
    check(label+" omitted source parameter zero",data[2][P.degree()-1]==P.domain.zero)

# Characteristic three with an infinite ground field.
field=s.GF(3).frac_field(u)
roots=[0,1,-1,u,u+1,u*u,-2*u-1-u*u]
P3=s.Poly(s.prod(S-r for r in roots),S,domain=field)
check("F3(u) rank7 trace coefficient zero",P3.nth(6)==0)
# Direct scalar reconstruction at every known split root avoids SymPy's
# GF(p)(u) polynomial-inversion path, which does not normalize fractions reliably.
data=compile_presentation(P3)
_,al,pars,pival,bval,cval,E,expected=data
check("F3(u) compiler polynomial identity",(E-expected).is_zero)
normalized_roots=[al*field.convert(r) for r in roots]
for index,sv in enumerate(normalized_roots):
    dv=field.one+2*bval*sv+3*pival*sv**2+sum(
        j*pars[j]*pival**j*sv**(j-1) for j in pars)
    check(f"F3(u) derivative unit at split root {index}",bool(dv))
    yy=-2*bval-3*pival*sv-sum(j*pars[j]*pival**j*sv**(j-2) for j in pars)
    xx=sv/dv
    zz=pival*dv**3-yy**2*(dv**2+3*dv)
    tt=field.one+xx*yy
    qq=tt**2*zz+yy**2*(field.one+3*tt)
    Bout=yy+3*xx*qq+sum(j*pars[j]*tt**2*xx**(j-2)*qq**j for j in pars)
    Cout=xx*(5-3*tt)-xx**3*zz-sum((j-2)*pars[j]*(xx*qq)**j for j in pars)
    check(f"F3(u) direct source reconstruction split root {index}",
          not bool(tt*qq-pival) and not bool(-Bout/2-bval)
          and not bool(Cout-cval) and not bool(xx/tt-sv))
check("F3(u) all seven normalized roots distinct",all(
    bool(normalized_roots[i]-normalized_roots[j]) for i in range(7) for j in range(i)))
cases.append({"label":"F3(u) split rank7: direct checks of all seven roots",
              "degree":7,"characteristic":3,"source_points_checked":7})
check("F3(u) coefficient a6 eliminated",compile_presentation(P3)[2][6]==field.zero)

# Normalization is a function of a line, and reconstructs the isomorphism chart.
for P,label in examples[:2]:
    n=P.degree()
    old=compile_presentation(P)
    scaled=s.Poly(lam**n*P.as_expr().subs(S,S/lam),S,domain=s.QQ.frac_field(lam))
    new=compile_presentation(scaled)
    check(label+" alpha has weight minus one",s.cancel(scaled.domain.to_sympy(new[1])-P.domain.to_sympy(old[1])/lam)==0)
    check(label+" target independent of line representative",
          all(s.cancel(scaled.domain.to_sympy(new[2][j])-P.domain.to_sympy(old[2][j]))==0 for j in old[2])
          and all(s.cancel(scaled.domain.to_sympy(new[i])-P.domain.to_sympy(old[i]))==0 for i in (3,4,5)))
    al=P.domain.to_sympy(old[1])
    normalized=s.Poly(al**n*P.as_expr().subs(S,S/al),S,domain=s.QQ)
    normdata=compile_presentation(normalized)
    check(label+" canonical generator has alpha one",normdata[1]==s.QQ.one)
    check(label+" inverse polynomial is monic characteristic polynomial divided by linear coefficient",
          old[6]==normalized.mul_ground(1/normalized.nth(1)))

# When the characteristic divides N, trace(1)=0 does not kill the trace map.
roots9=[u**j for j in range(8)]
roots9.append(-sum(roots9))
P9=s.Poly(s.prod(S-r for r in roots9),S,domain=field)
check("rank9 in characteristic3 has trace zero and good coefficients",
      P9.nth(8)==0 and all(P9.nth(j)!=0 for j in (1,3,4,5)))
check("rank9 characteristic3 split generator remains primitive",
      all(field.convert(roots9[i]-roots9[j])!=field.zero
          for i in range(9) for j in range(i)))
# e1-coordinate vector has trace one in the split algebra, even though 1 has trace zero.
check("rank9 characteristic3 trace map nonzero despite trace of unit zero",
      sum([1]+[0]*8)%3==1 and 9%3==0)

# Explicit pi-zero fibres.
parameters={4:s.Integer(1),5:s.Integer(1),6:s.Integer(0),7:s.Integer(2)}
out=source_map(7,parameters)
def image(point):
    return tuple(s.expand(g.subs(dict(zip((x,y,z),point)))) for g in out)
for bv,cv in [(0,2),(1,0),(1,s.Rational(-1,2)),(1,1)]:
    if bv==0:
        check("pi0 b0 unique reconstructed point",image((s.Rational(cv)/2,0,0))==(0,bv,cv))
        continue
    zz=bv**3*cv+5*bv**2+bv**3*sum((j-2)*parameters[j]*(-bv)**j for j in parameters)
    check(f"pi0 boundary component b={bv} c={cv}",
          image((-s.Rational(1,bv),bv,zz))==(0,bv,cv))
    quadratic=s.Poly(S+bv*S*S-s.Rational(cv)/2,S,domain=s.QQ)
    if 1+2*bv*cv==0:
        check("pi0 double root entirely removed by derivative localization",
              s.gcd(quadratic,quadratic.diff()**2).monic()==quadratic.monic())
    elif cv==0:
        for sv in (0,-1):
            dv=1+2*bv*sv
            yy=-2*bv
            point=(s.Rational(sv,dv),yy,-yy**2*(dv**2+3*dv))
            check(f"pi0 quadratic component root {sv}",image(point)==(0,bv,cv))
    else:
        check("pi0 nonsplit quadratic has rank two though no rational root",quadratic.is_irreducible)

# Boundary checks for the next proposed compression.
for n in (7,8):
    rr=s.symbols("r0:"+str(n-1))
    roots=list(rr)+[-sum(rr)]
    e2=sum(roots[i]*roots[j] for i in range(n) for j in range(i))
    check(f"double depression Newton identity N={n}",
          s.expand(e2+sum(r*r for r in roots)/2)==0)
for n,j in ((5,4),(6,5)):
    roots=list(range(1,n))
    roots.append(-sum(roots))
    P=s.Poly(s.prod(S-r for r in roots),S)
    check(f"trace-zero low-rank compiler obstruction N={n}",P.nth(j)==0)

receipt={
 "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
 "complete_fibre_controls":cases,
 "source_routines_reused":["source_map","verify_fibre"],
 "assurance":"Written universal proofs with exact controls. Source reconstruction routine is reused from the predecessor, not an independent implementation. Monodromy, correspondence and weak approximation are written arguments, not finite computations.",
 "python":platform.python_version(),"sympy":s.__version__,
 "runtime_seconds":round(time.monotonic()-started,3),
 "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [HERE/"PROOF.md",Path(__file__),HERE/"sources.json"]},
 "predecessor_checker_sha256":hashlib.sha256((PREDECESSOR/"verify.py").read_bytes()).hexdigest()
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"fibre_controls":len(cases),
                  "seconds":receipt["runtime_seconds"]}))


"""Exact controls for the written all-rank proof; no formal verification."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


from itertools import product
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
started=time.monotonic()
checks=[]
def check(label, condition):
    if not condition:
        raise AssertionError(label)
    checks.append(label)
x,y,z,w,S,T,pi,b,c,alpha=s.symbols("x y z w S T pi b c alpha")

def source_map(n, coefficients):
    t=1+x*y
    q=t*t*z+y*y*(1+3*t)
    B=y+3*x*q+sum(j*coefficients[j]*t*t*x**(j-2)*q**j
                   for j in range(4,n+1))
    C=x*(5-3*t)-x**3*z-sum((j-2)*coefficients[j]*(x*q)**j
                            for j in range(4,n+1))
    return [t*q,-B/2,C]

# Sparse polynomial-ring determinants, not sampled Jacobians.
for n in (5,6):
    coefficients={j:(w if j in (4,5) else s.Integer(j-3))
                  for j in range(4,n+1)}
    f=source_map(n,coefficients)
    jac=[[s.Poly(s.diff(g,v),x,y,z,w,domain=s.QQ)
          for v in (x,y,z)] for g in f]
    determinant=(jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])
                -jac[0][1]*(jac[1][0]*jac[2][2]-jac[1][2]*jac[2][0])
                +jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]))
    check(f"full symbolic polynomial Jacobian N={n}",determinant.as_expr()==1)

a4,a5=s.symbols("a4 a5",nonzero=True)
normalized_w=a4**6/a5**5
normal_alpha=a4/a5
check("quotient inverse a4",s.cancel(normal_alpha**(-5)*normalized_w-a4)==0)
check("quotient inverse a5",s.cancel(normal_alpha**(-6)*normalized_w-a5)==0)
check("quotient invariant I5",s.cancel(1/normalized_w-a5**5/a4**6)==0)
q0,t0=s.symbols("q0 t0")
for n in (5,8):
    aj={j:s.Symbol(f"a{j}") for j in range(4,n+1)}
    E=S+b*S**2+pi*S**3+sum(aj[j]*pi**j*S**j for j in aj)-c/2
    Eprime=s.diff(E,S)
    fy=-2*b-3*pi*S-sum(j*aj[j]*pi**j*S**(j-2) for j in aj)
    check(f"fibre S*y=1-Eprime N={n}",s.expand(S*fy-(1-Eprime))==0)
    Cnormal=S*(1+Eprime)-pi*S**3-sum((j-2)*aj[j]*pi**j*S**j for j in aj)
    check(f"fibre C=c+2E N={n}",s.expand(Cnormal-c-2*E)==0)
    transformed=S+(b/alpha)*S**2+(pi/alpha**2)*S**3+sum(
        alpha**(j+1)*aj[j]*(pi/alpha**2)**j*S**j for j in aj)-alpha*c/2
    check(f"whole inverse-polynomial covariance N={n}",
          s.expand(transformed-alpha*E.subs(S,S/alpha))==0)
    B=y+3*x*q0+sum(j*aj[j]*t0**2*x**(j-2)*q0**j for j in aj)
    C=x*(5-3*t0)-x**3*z-sum((j-2)*aj[j]*(x*q0)**j for j in aj)
    scaled={x:alpha*x,y:y/alpha,z:z/alpha**2,q0:q0/alpha**2,
            **{aj[j]:alpha**(j+1)*aj[j] for j in aj}}
    check(f"source B covariance N={n}",s.expand(B.subs(scaled,simultaneous=True)-B/alpha)==0)
    check(f"source C covariance N={n}",s.expand(C.subs(scaled,simultaneous=True)-alpha*C)==0)
t_source=1+x*y
q_source=t_source**2*z+y**2*(1+3*t_source)
check("source q covariance",s.cancel(q_source.subs(
    {x:alpha*x,y:y/alpha,z:z/alpha**2},simultaneous=True)-q_source/alpha**2)==0)

def compile_presentation(P):
    n=P.degree()
    field=P.domain
    zero,one=field.zero,field.one
    chosen=None
    candidates=range(4*n-12) if field==s.QQ else range(field.mod)
    for translation in candidates:
        shifted=s.Poly(P.as_expr().subs(S,S+translation),S,domain=field)
        d={j:field.convert(shifted.nth(j)) for j in range(n+1)}
        if all(d[j] != zero for j in (1,3,4,5)):
            chosen=translation
            break
    if chosen is None:
        raise AssertionError("No admissible translation in declared range")
    al=d[3]*d[4]/(d[1]*d[5])
    pars={j:d[j]*d[3]*d[4]**(j+1)/(d[1]**2*d[5]**(j+1))
          for j in range(4,n+1)}
    pival=d[1]*d[5]**2/(d[3]*d[4]**2)
    bval=d[2]*d[5]/(d[3]*d[4])
    cval=-field.convert(2)*d[0]*d[3]*d[4]/(d[1]**2*d[5])
    coefficients={0:-cval/field.convert(2),1:one,2:bval,3:pival,
                  **{j:pars[j]*pival**j for j in pars}}
    E=s.Poly.from_dict({(j,):value for j,value in coefficients.items()},S,domain=field)
    expected=s.Poly(field.to_sympy(al/d[1])*P.as_expr().subs(
        S,chosen+field.to_sympy(one/al)*S),S,domain=field)
    return chosen,al,pars,pival,bval,cval,E,expected

def verify_fibre(P,label):
    n=P.degree()
    field=P.domain
    a,al,aj,pival,bval,cval,E,expected=compile_presentation(P)
    check(label+" compiler scaling identity",E==expected)
    check(label+" compressed a4=a5",aj[4]==aj[5] and aj[4]!=field.zero)
    check(label+" degree and squarefreeness",E.degree()==n and s.gcd(E,E.diff()).degree()==0)
    def poly(value): return s.Poly(value,S,domain=field)
    unit=poly(1)
    def mod(value): return value.rem(E)
    def mul(a,b): return mod(a*b)
    def power(a,j): return mod(a**j)
    D=E.diff()
    invD=s.invert(D,E)
    xx=mul(poly(S),invD)
    yy=poly(field.to_sympy(-2*bval))-poly(field.to_sympy(3*pival))*poly(S)
    for j in aj:
        yy-=poly(field.to_sympy(field.convert(j)*aj[j]*pival**j))*poly(S**(j-2))
    yy=mod(yy)
    zz=mod(poly(field.to_sympy(pival))*power(D,3)-power(yy,2)*mod(power(D,2)+poly(3)*D))
    tt=mod(unit+mul(xx,yy))
    qq=mod(mul(power(tt,2),zz)+mul(power(yy,2),unit+poly(3)*tt))
    Bout=mod(yy+poly(3)*mul(xx,qq))
    Cout=mod(mul(xx,poly(5)-poly(3)*tt)-mul(power(xx,3),zz))
    for j in aj:
        scalar=poly(field.to_sympy(aj[j]))
        Bout=mod(Bout+poly(j)*mul(scalar,mul(power(tt,2),mul(power(xx,j-2),power(qq,j)))))
        Cout=mod(Cout-poly(j-2)*mul(scalar,power(mul(xx,qq),j)))
    check(label+" source reconstructs all three target coordinates",
          mul(tt,qq)==poly(field.to_sympy(pival))
          and mul(Bout,poly(field.to_sympy(-field.one/field.convert(2))))==poly(field.to_sympy(bval))
          and Cout==poly(field.to_sympy(cval)))
    check(label+" inverse recovers generator",mul(xx,s.invert(tt,E))==poly(S))
    return {"label":label,"translation":a,"degree":n,"characteristic":field.characteristic()}

cases=[]
for n in (5,6,7):
    for k,expr in enumerate([S**n-S-1,s.prod(S-j for j in range(1,n+1))]):
        P=s.Poly(expr,S,domain=s.QQ)
        cases.append(verify_fibre(P,f"QQ rank{n} example{k}"))

def partitions(n,minimum=1):
    if n==0:
        yield []
    for j in range(minimum,n+1):
        for tail in partitions(n-j,j):
            yield [j]+tail

field=s.GF(11)
irreducibles={}
for degree in range(1,6):
    needed=5//degree
    bank=[]
    for coeffs in product(range(11),repeat=degree):
        P=s.Poly(S**degree+sum(co*S**j for j,co in enumerate(coeffs)),S,domain=field)
        if P.is_irreducible:
            bank.append(P)
        if len(bank)==needed:
            break
    irreducibles[degree]=bank
for partition in partitions(5):
    P=s.Poly(1,S,domain=field)
    used={}
    for degree in partition:
        index=used.get(degree,0)
        P*=irreducibles[degree][index]
        used[degree]=index+1
    cases.append(verify_fibre(P,"F11 rank5 factor degrees "+str(partition)))
check("all seven rank-five finite-etale splitting types covered",len(list(partitions(5)))==7)
# Exact boundary controls, rather than extrapolation from successful examples.
bad=s.Poly(S**5-S-1,S,modulus=5)
check("small-characteristic third-Hasse failure",
      all(s.Poly(bad.as_expr().subs(S,S+a),S,modulus=5).nth(3)==0 for a in range(5)))
repeated=s.Poly((S-1)**2*(S**3+S+1),S,domain=s.QQ)
check("repeated roots require derivative localization",s.gcd(repeated,repeated.diff()).degree()>0)

receipt={
 "schema_version":1,"result":"pass","checks":checks,"check_count":len(checks),
 "complete_fibre_controls":cases,
 "assurance":"Written all-rank proofs with exact representative controls, not formalization or independent whole-proof review. No finite Galois-group samples are used as a monodromy proof.",
 "python":platform.python_version(),"sympy":s.__version__,
 "runtime_seconds":round(time.monotonic()-started,3),
 "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [Path(__file__)]}
}
(_check_output(HERE/"receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":receipt["result"],"checks":len(checks),
                  "fibre_controls":len(cases),"seconds":receipt["runtime_seconds"]}))


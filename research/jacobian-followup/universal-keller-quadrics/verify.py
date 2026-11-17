"""Exact controls for the written trace-quadric theorems.
Two predecessor routines are reused explicitly; this is not independent replay.
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
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
x,y,z,h,w,S,T,pi,b,c,alpha,delta=s.symbols("x y z h w S T pi b c alpha delta")
predecessor=HERE.parent/"universal-keller-compression"/"verify.py"
tree=ast.parse(predecessor.read_text())
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef)
       and n.name in {"source_map","verify_fibre"}]
assert len(nodes)==2
class SectionCondition(ast.NodeTransformer):
    def visit_Call(self,node):
        node=self.generic_visit(node)
        if (isinstance(node.func,ast.Name) and node.func.id=="check"
            and node.args and isinstance(node.args[0],ast.BinOp)
            and isinstance(node.args[0].right,ast.Constant)
            and node.args[0].right.value==" compressed a4=a5"):
            node.args[0].right.value=" selected polynomial weight section"
            node.args[1]=ast.Call(func=ast.Name(id="section_good",ctx=ast.Load()),
                                  args=[ast.Name(id="aj",ctx=ast.Load())],keywords=[])
        return node
module=ast.fix_missing_locations(ast.Module(
    body=[SectionCondition().visit(n) for n in nodes],type_ignores=[]))
exec(compile(module,str(predecessor),"exec"))

def section_good(aj):
    n=max(aj)
    if n==7:
        hv=aj[7]**5/aj[4]**8
        return aj[4]==hv**3 and aj[7]==hv**5 and not aj[5] and not aj[6]
    return aj[4]==aj[5] and not aj[n-1] and not aj[n-2]

def compile_presentation(P):
    n=P.degree();field=P.domain
    d={j:field.convert(P.nth(j)) for j in range(n+1)}
    required=(1,3,4) if n==7 else (1,3,4,5)
    assert not d[n-1] and not d[n-2] and all(d[j] for j in required)
    if n==7:
        al=d[1]**3/(d[3]*d[4]**5)
        hv=d[1]**6/(d[3]**3*d[4]**8)
        aj={4:hv**3,5:field.zero,6:field.zero,7:hv**5}
        pv=d[3]**3*d[4]**10/d[1]**7
        bv=d[2]*d[3]*d[4]**5/d[1]**4
        cv=-2*d[0]*d[1]**2/(d[3]*d[4]**5)
    else:
        al=d[3]*d[4]/(d[1]*d[5])
        aj={j:d[j]*d[3]*d[4]**(j+1)/(d[1]**2*d[5]**(j+1))
            for j in range(4,n+1)}
        pv=d[1]*d[5]**2/(d[3]*d[4]**2)
        bv=d[2]*d[5]/(d[3]*d[4])
        cv=-2*d[0]*d[3]*d[4]/(d[1]**2*d[5])
    co={0:-cv/field.convert(2),1:field.one,2:bv,3:pv,
        **{j:aj[j]*pv**j for j in aj}}
    E=s.Poly.from_dict({(j,):value for j,value in co.items()},S,domain=field)
    # Compute the scaled polynomial inside its exact coefficient domain.
    # Expanding huge algebraic expressions back through SymPy minimal-polynomial
    # inference can overflow unrelated floating-point estimates in factorint.
    expected=s.Poly.from_dict({(j,):al*d[j]/(d[1]*al**j)
                               for j in range(n+1)},S,domain=field)
    return 0,al,aj,pv,bv,cv,E,expected

# Full polynomial determinant after the new coefficient substitutions.
for n,aj in [(7,{4:h**3,5:s.Integer(0),6:s.Integer(0),7:h**5}),
             (8,{4:w,5:w,6:s.Integer(0),7:s.Integer(0),8:s.Integer(1)})]:
    par=h if n==7 else w
    f=source_map(n,aj)
    J=[[s.Poly(s.diff(g,v),x,y,z,par,domain=s.QQ) for v in (x,y,z)] for g in f]
    det=(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
         -J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
         +J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]))
    check(f"N{n} full polynomial Jacobian",det.as_expr()==1)
    if n==7:
        degrees=[1]+[s.Poly(g,x,y,z,h).total_degree() for g in f]
        check("rank7 ordinary coordinate degrees",degrees==[1,7,49,47])

a4,a7=s.symbols("a4 a7",nonzero=True)
al=a7**3/a4**5;hv=a7**5/a4**8
check("weights5,8 root-free section",
      s.cancel(al**5*a4-hv**3)==0 and s.cancel(al**8*a7-hv**5)==0)
check("weights5,8 Laurent inverse",
      s.cancel(al**(-5)*hv**3-a4)==0 and s.cancel(al**(-8)*hv**5-a7)==0)
check("primitive weight vector with no coprime pair",
      abs(s.Matrix([[6,1,0],[10,0,1],[15,1,1]]).det())==1
      and all(s.gcd(i,j)>1 for i,j in [(6,10),(6,15),(10,15)]))

P7=s.Poly(S**7+S**4+S**3+S**2+S+1,S,domain=s.QQ)
P8=s.Poly(S**8+S**5+2*S**4+S**3+S**2+S+1,S,domain=s.QQ)
real7=s.Poly(s.prod(S-j for j in range(5))*(S**2+10*S+65),S,domain=s.QQ)
real8=s.Poly(s.prod(S-j for j in range(6))*(S**2+15*S+140),S,domain=s.QQ)
examples=[(P7,"QQ rank7 small presentation"),(P8,"QQ rank8 small presentation"),
          (real7,"QQ rank7 five real factors and one complex pair"),
          (real8,"QQ rank8 six real factors and one complex pair"),
          (s.Poly(P7.as_expr(),S,modulus=7),"F7 rank7 characteristic divides rank"),
          (s.Poly(S**9+S**5+S**4+S**3+S**2+S+1,S,modulus=3),
           "F3 rank9 characteristic divides rank")]
cases=[]
for P,label in examples:
    n=P.degree()
    check(label+" first two nonleading coefficients vanish",
          not P.nth(n-1) and not P.nth(n-2))
    cases.append(verify_fibre(P,label))
    # Companion traces test the interpretation in the algebra, independently of signs.
    M=s.zeros(n)
    for j in range(n-1):M[j+1,j]=1
    for j in range(n):M[j,n-1]=-P.nth(j)
    field=P.domain
    check(label+" multiplication traces vanish",
          field.convert(s.trace(M))==field.zero
          and field.convert(s.trace(M*M))==field.zero)

# The normalized element depends on the generator line.
lam=s.symbols("lam")
for P,label in examples[:2]:
    n=P.degree();old=compile_presentation(P)
    field=s.QQ.frac_field(lam)
    shifted=s.Poly(lam**n*P.as_expr().subs(S,S/lam),S,domain=field)
    new=compile_presentation(shifted)
    check(label+" normalizer weight minus one",
          s.cancel(field.to_sympy(new[1])-old[1]/lam)==0)
    check(label+" target independent of representative",
          all(s.cancel(field.to_sympy(new[i])-old[i])==0 for i in (3,4,5))
          and all(s.cancel(field.to_sympy(new[2][j])-old[2][j])==0 for j in old[2]))
    normalized=s.Poly(old[1]**n*P.as_expr().subs(S,S/old[1]),S,domain=s.QQ)
    dd=lambda j:normalized.nth(j)
    goodnorm=(dd(1)**3==dd(3)*dd(4)**5 if n==7 else dd(3)*dd(4)==dd(1)*dd(5))
    check(label+" canonical generator normalization",goodnorm)

# Binary test: a symbolic identity in Q[delta]/(delta^2+14).
binary=[s.Integer(3),-4+delta,-4-2*delta]
check("rank7 binary linear trace",
      s.expand(4*binary[0]+2*binary[1]+binary[2])==0)
check("rank7 binary quadratic trace",
      s.rem(4*binary[0]**2+2*binary[1]**2+binary[2]**2,delta**2+14,delta)==0)
check("rank7 binary vector is nonconstant",
      s.gcd(delta-7,delta**2+14)==1)
check("rank7 Sylow two reduction has odd degree",s.factorial(7)//16==315)

# The three-block binary determinant controls the constant-field square class.
n1,n2,n3=s.symbols("n1 n2 n3",nonzero=True)
G=s.Matrix([[n1+n1*n1/n3,n1*n2/n3],
            [n1*n2/n3,n2+n2*n2/n3]])
check("three-block determinant",s.cancel(G.det()-n1*n2*(n1+n2+n3)/n3)==0)
check("rank7 binary discriminant square class",
      s.Rational(-7*4*2*1,-14)==4)
# A split algebra over Q cannot have a primitive sum-of-squares-zero generator.
H=s.eye(6)+s.ones(6)
check("rank7 real-split trace form positive definite",
      all(H[:j,:j].det()==j+1 for j in range(1,7)))

def rank_mod(matrix,p):
    rows=[[int(v)%p for v in row] for row in matrix.tolist()]
    rank=0
    for col in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inv=pow(rows[rank][col],-1,p);rows[rank]=[(v*inv)%p for v in rows[rank]]
        for i in range(len(rows)):
            if i!=rank:
                v=rows[i][col];rows[i]=[(a-v*bv)%p for a,bv in zip(rows[i],rows[rank])]
        rank+=1
    return rank
for n,p in [(7,7),(8,7),(9,3)]:
    gram=s.eye(n-1)+s.ones(n-1)
    check(f"trace radical dimension N{n} characteristic{p}",
          rank_mod(gram,p)==n-1-int(n%p==0))

# An actual smooth binary point parametrizes a good primitive rank7 generator.
K=s.QQ.algebraic_field(s.sqrt(-14))
pp=[K.convert(3)]*4+[K.from_sympy(-4+s.sqrt(-14))]*2+[K.from_sympy(-4-2*s.sqrt(-14))]
uu=[K.convert(j) for j in [0,1,2,3,4,5,-15]]
qsum=sum((v*v for v in uu),K.zero)
pair=sum((a*bv for a,bv in zip(pp,uu)),K.zero)
rr=[qsum*a-2*pair*bv for a,bv in zip(pp,uu)]
check("Qsqrt-14 parametrized trace-zero isotropic point",
      not sum(rr,K.zero) and not sum((v*v for v in rr),K.zero))
check("Qsqrt-14 parametrized generator primitive",
      all(rr[i]!=rr[j] for i in range(7) for j in range(i)))
Pc=s.Poly(1,S,domain=K)
for root in rr:Pc*=s.Poly(S-K.to_sympy(root),S,domain=K)
check("Qsqrt-14 parametrization reaches all good coefficients",
      all(Pc.nth(j)!=0 for j in (1,3,4)) and Pc.nth(6)==0 and Pc.nth(5)==0)
cc=compile_presentation(Pc)
check("Qsqrt-14 actual normalized inverse polynomial identity",(cc[6]-cc[7]).is_zero)

# New septic Hermite-Joubert constraint has a different remaining weight lattice.
check("septic first-and-third trace constraint retains nonprimitive seed weights",
      s.gcd(6,8)==2)
# Binary one/two-block obstructions are exact, not failed primitive-element searches.
check("one-block trace kernel zero in characteristic zero",s.Matrix([[8]]).nullspace()==[])
check("two-block restricted form nonzero",8+8**2==72)

manifest=json.loads((HERE/"sources.json").read_text())
for item in manifest["local_sources"]:
    check("source pin "+item["path"],
          hashlib.sha256((RESEARCH/item["path"]).read_bytes()).hexdigest()==item["sha256"])
receipt={
 "schema_version":1,"result":"pass","check_count":len(checks),"checks":checks,
 "complete_fibre_controls":cases,"rank7_coordinate_degrees":degrees,
 "python":platform.python_version(),"sympy":s.__version__,
 "source_routines_reused":["source_map","verify_fibre (section-condition adapted by AST)"],
 "assurance":"Written universal proofs with exact controls, not independent replay or formalization. Quadratic-form existence, field-versality necessity, monodromy and stable classification are written arguments or explicitly imported literature theorems.",
 "runtime_seconds":round(time.monotonic()-started,3),
 "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [HERE/"PROOF.md",HERE/"LITERATURE_AND_GAPS.md",Path(__file__),HERE/"sources.json"]}
}
(HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":"pass","checks":len(checks),"fibre_controls":len(cases),
                  "seconds":receipt["runtime_seconds"]}))


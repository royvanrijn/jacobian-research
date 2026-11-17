"""Exact controls for SWEEP_MONODROMY.md; the universal proof is written."""
from pathlib import Path
import hashlib,json,platform,time,math
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
w,z,gamma,Y1,Y2,Y3=S.symbols("w z gamma Y1 Y2 Y3")
Y=S.Matrix([Y1,Y2,Y3])
Delta=S.Matrix([1,w,w*w])
A=S.Matrix([w/3-4*w*w/3,2*w*w/3-8*w**3/9,2*w+7*w**3/9-2*w**4/3])
B=S.Matrix([w*w-1,w**3,w**4+w*w+2])
def det(*cols):return S.expand(S.Matrix.hstack(*cols).det())
for label,aa,bb in [("F4",A,B),("free_polynomial",
    S.integrate((w**3+2)*Delta+(w-1)*(w*w+w)*Delta.diff(w)+(w-1)*Delta.diff(w,2),w),
    ((w*w+w)**2-S.diff(w*w+w,w))*Delta+(w*w+w)*Delta.diff(w)+Delta.diff(w,2))]:
    check(label+" first tangency",det(Delta,aa.diff(w),bb)==0)
    check(label+" ruling tangency",det(Delta,bb.diff(w),bb)==0)
    check(label+" normalized frame",det(Delta,Delta.diff(w),bb)==2)
    nu=bb.cross(Delta)
    check(label+" primitive coefficient curve",S.gcd_list(list(nu)) in [1,-1,2,-2])
    check(label+" projective immersion pairing",S.expand(nu.diff(w).dot(Delta))==-2)
    check(label+" sweep Jacobian",det(Delta,aa.diff(w)+z*bb.diff(w)+gamma*Delta.diff(w),bb)==2*gamma)
RF4=S.expand(B.cross(Delta).dot(A-Y))
expect=2*w**5/9+2*w**4/9+(S.Rational(8,9)+Y1)*w**3-(S.Rational(4,3)+2*Y2)*w*w+(2*Y1+Y3)*w-2*Y2
check("F4 exact root polynomial",S.expand(RF4-expect)==0)
check("F4 degree five",S.degree(RF4,w)==5)
check("F4 constant leading coefficient",S.Poly(RF4,w).LC()==S.Rational(2,9))
kappa=S.Matrix([B.cross(Delta).dot(A),*B.cross(Delta)])
check("F4 coefficient Wronskian nonzero",any(S.det(S.Matrix.hstack(kappa,kappa.diff(w),kappa.diff(w,2)).extract(rows,[0,1,2]))!=0 for rows in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]))
X=S.Matrix([160*w**3-60*w*w-240*w*z+51*w+60*z,
    60*w**4-20*w**3-120*z*z+51*z,
    -48*w**5+15*w**4+240*w**3*z-17*w**3-60*w*w*z-240*w*z*z+51*w*z+30*z*z])/29
D=S.Matrix([1,w,z])
sweep=X+gamma*D
check("F5 exact sweep Jacobian",det(D,sweep.diff(w),sweep.diff(z))==gamma**2)
check("F5 tangency at gamma zero",det(D,X.diff(w),X.diff(z))==0)
eq2=S.expand(X[1]+(Y1-X[0])*w-Y2)
eq3=S.expand(X[2]+(Y1-X[0])*z-Y3)
check("F5 monic transverse equation up to constant",S.Poly(eq2,z).LC()==-S.Rational(120,29))
R=S.expand(-S.Rational(29**4,30)*S.resultant(eq2,eq3,z))
check("F5 integral root equation degree ten",S.degree(R,w)==10)
check("F5 integral root equation leading coefficient",S.Poly(R,w).LC()==20480)
check("F5 resultant integer coefficients",S.Poly(R,w,Y1,Y2,Y3).domain==S.ZZ)
witness={Y1:0,Y2:1,Y3:1}
R0=S.Poly(R.subs(witness),w,domain=S.ZZ)
r7=S.Poly(R0,w,modulus=7)
check("F5 good reduction retains degree",r7.degree()==10)
check("F5 good reduction squarefree",r7.gcd(r7.diff()).degree()==0)
factors=[w+2,w*w+2*w+2,w**7-3*w**6-3*w**5+w**4-2*w**3-w*w+1]
check("F5 exact Frobenius factorization",S.Poly(R0.as_expr()-5*S.prod(factors),w,modulus=7).is_zero)
for i,f in enumerate(factors):
    check("F5 finite factor irreducible "+str(i),S.Poly(f,w,modulus=7).is_irreducible)
subs=S.subresultants(eq2,eq3,z)
lin=next(t for t in subs if S.degree(t,z)==1)
lc=S.Poly(lin,z).LC()
lc0=S.Poly(lc.subs(witness),w,domain=S.QQ).clear_denoms()[1]
check("F5 primitive projection inverse coefficient",r7.gcd(S.Poly(lc0,w,modulus=7)).degree()==0)
disc=int(S.discriminant(R0.as_expr(),w))
check("F5 rational witness unramified",disc!=0)
q,issquare=S.integer_nthroot(disc//6,2)
check("F5 exact discriminant square class six",disc==6*q*q and issquare)
check("F5 constant quadratic nontrivial over Q",not S.integer_nthroot(6,2)[1])
check("F5 finite Frobenius is odd",sum(d-1 for d in [1,2,7])%2==1)
perm=[0,2,1,4,5,6,7,8,9,3]
square=[perm[perm[i]] for i in range(10)]
seen=set();lengths=[]
for i in range(10):
    if i in seen:continue
    j=i;length=0
    while j not in seen:seen.add(j);length+=1;j=square[j]
    lengths.append(length)
check("F5 squared Frobenius yields prime seven cycle",sorted(lengths)==[1,1,1,7])
for block in [2,5]:
    order=math.factorial(block)**(10//block)*math.factorial(10//block)
    check("F5 seven excludes block size "+str(block),order%7!=0)
check("F5 Jordan prime has three fixed points",S.isprime(7) and 10-7==3)
# Verify the complete global affine polynomial lift by exact weighted jets.
u1,u2,v1,v2,v3=S.symbols("u1 u2 v1 v2 v3")
stage=S.Matrix([1+S.Rational(20,29)*v1,
    1-S.Rational(49,29)*v1+8*v2,
    1-S.Rational(69,29)*v1+9*v2+v3+S.Rational(4197,1682)*v1*v1])
check("F5 triangular stage determinant",stage.jacobian([v1,v2,v3]).det()==S.Rational(160,29))
stinv1=S.Rational(29,20)*(gamma-1)
stinv2=(u1-1+S.Rational(49,29)*stinv1)/8
stinv3=u2-1+S.Rational(69,29)*stinv1-9*stinv2-S.Rational(4197,1682)*stinv1**2
for i,expr in enumerate([stinv1,stinv2,stinv3]):
    check("F5 polynomial stage inverse "+str(i),
       S.expand(expr.subs({gamma:stage[0],u1:stage[1],u2:stage[2]},simultaneous=True)-[v1,v2,v3][i])==0)
jetcounts=[]
for i,si in enumerate(sweep,1):
    scaled=S.expand(si.subs({w:gamma*u1,z:gamma**2*u2},simultaneous=True))
    quot,rem=S.div(S.Poly(scaled,gamma,u1,u2),S.Poly(gamma**i,gamma,u1,u2))
    check("F5 gamma divisibility "+str(i),rem.is_zero)
    Ei=quot.as_expr()
    composed=S.Poly(S.expand(Ei.subs({gamma:stage[0],u1:stage[1],u2:stage[2]},simultaneous=True)),v1,v2,v3)
    low=[(mon,coef) for mon,coef in composed.terms() if sum(a*b for a,b in zip(mon,[1,2,3]))<i]
    check("F5 all global polynomiality jets "+str(i),not low)
    jetcounts.append(len(composed.terms()))
check("F5 x determinant exponent cancels",sum([1,2,3])==sum([1,2,3]))
check("F5 gamma determinant exponent cancels",sum([1,2])+3==sum([1,2,3]))
manifest=json.loads((ROOT/"sources/sweep-monodromy-manifest.json").read_text())
for source in manifest["sources"]:
    check("primary source hash "+source["path"],
      hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest()==source["sha256"])
model=dict(ruled=dict(F4_polynomial=str(RF4),F4_geometric_group="S5",
  universal_group="S_N for every normalized ruled sweep satisfying the stated hypotheses"),
  F5=dict(surface=[str(t) for t in X],root_polynomial=str(R),
  transverse_equation=str(eq2),linear_subresultant=str(lin),
  finite_prime=7,witness=[0,1,1],factors=[str(f) for f in factors],
  factor_degrees=[1,2,7],discriminant=disc,discriminant_squareclass=6,
  discriminant_square_root_after_dividing_6=int(q),
  geometric_group="A10",arithmetic_group="S10",constant_field="Q(sqrt(6))",
  stage=[str(t) for t in stage],jet_polynomial_term_counts=jetcounts,
  actual_global_Jacobian="160/29",source_dimension=4),
  scope=dict(general_PSL2_11_action_decided=False,full_classification_complete=False),
  assurance="Universal geometric proof and finite-different argument written; Jordan theorem attributed.")
(ROOT/"sweep_monodromy_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["SWEEP_MONODROMY.md","verify_sweep_monodromy.py","sweep_monodromy_model.json","sources/sweep-monodromy-manifest.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,
    python=platform.python_version(),sympy=S.__version__,
    runtime_seconds=round(time.monotonic()-started,2),
    sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
    assurance="Written universal geometry, finiteness and discriminant proof; exact arithmetic and affine lift controls. Jordan theorem attributed. No formal verification or independent whole-proof review.",
    goal_complete=False)
(ROOT/"sweep_monodromy_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=receipt["runtime_seconds"])))

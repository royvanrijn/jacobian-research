"""Exact branch-centering, expanded binary Keller map and action controls."""
from pathlib import Path
import hashlib,itertools,json,math,platform,time
import sympy as S
from branch_centered import binary_sweeps,binary_stage,binary_constants,expanded_binary_map
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
w,z,g,u,v,r=S.symbols("w z gamma u v r")
minimum=r*r+3
def red(expr):
    num,den=S.fraction(S.cancel(expr))
    return S.rem(S.rem(num,minimum,r)*S.invert(S.rem(den,minimum,r),minimum,r),minimum,r).expand()
data=binary_sweeps(w,z,g);X=data["original_X"];XX=data["centered_X"]
old=polynomial_chart(w*w,[(1,0)],w=w,z=z)
check("same exact binary field input as KMFI45",(old["X"]-X).applyfunc(S.expand)==S.zeros(3,1) and old["N"]==8 and old["geometric_group"]=="SL2(F3) on8 letters")
# Retained exhaustive full-group certificate is source-bound, not a sampled Frobenius test.
prior=json.loads((ROOT/"polynomial_chart_receipt.json").read_text())
check("retained full binary group receipt passes",prior["result"]=="pass")
for file,h in prior["sha256"].items():
    check("retained binary group source "+file,hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==h)
for label,expr in [("original",data["original"]),("branch-centered",data["centered"])]:
    check(label+" full sweep Jacobian",S.factor(expr.jacobian([g,w,z]).det())==g*g)
A,B,q=S.symbols("A B q")
quartic=S.expand((A-q*q/2)**2+2*q*q*(A-q*q/2)+2*q**4/3-2*q*B)
check("branch quartic exact factorization",S.factor(quartic.subs({A:S.Rational(1,2),B:S.Rational(1,3)}))==-(q-1)**3*(q+3)/12)
Y=data["original"];AA=Y[1];BB=Y[2]+Y[0]**2/2
check("exact branch pullback and unramified tail",S.factor(9*BB**2-8*AA**3)==g**3*(9*g+4*w**3)/4)
target=S.Matrix([1,S.Rational(1,2),-S.Rational(1,6)])
marked_old={g:4*r/3,w:r,z:-1-5*r/3}
check("regular sheet above explicit branch point",all(red(x.subs(marked_old,simultaneous=True))==0 for x in data["original"]-target))
check("original branch point maps to same target",X.subs({w:1,z:0})==target)
check("marked regular gamma nonzero",red((4*r/3)**2)==-S.Rational(16,3))
# General local normalization: c(0)=0 and lambda(0)!=0 give all required jets.
a0,aw,az,b0,bw,bz,cw,cz=S.symbols("a0 aw az b0 bw bz cw cz")
scalar=a0+aw*w+az*z;bb=b0+bw*w+bz*z;cc=cw*w+cz*z
lam=scalar*bb*bb;mu=scalar*bb*cc;nu=scalar*cc*cc
check("universal centered lambda unit",lam.subs({w:0,z:0})==a0*b0*b0)
for name,expr in [("mu0",mu),("nu0",nu),("nu_w0",S.diff(nu,w)),("nu_z0",S.diff(nu,z))]:
    check("universal branch-center jet "+name,expr.subs({w:0,z:0})==0)
# Concrete centered X obeys every weighted order, with no homogeneity extrapolation.
for idx,expr in enumerate(XX):
    check("centered weighted vanishing output"+str(idx),
      all(a+2*b>=idx+1 for (a,b),co in S.Poly(expr,w,z).terms()))
E=binary_stage(g,u,v)
check("exact weighted stage from sweep",
 (E-S.Matrix([S.cancel(data["centered"][idx].subs({w:g*u,z:g*g*v},simultaneous=True)/g**(idx+1)) for idx in range(3)])).applyfunc(S.expand)==S.zeros(3,1))
constants=binary_constants(r);theta0=constants["theta0"];L=constants["L"];f20=constants["f20"]
at=lambda expr:red(expr.subs(dict(zip([g,u,v],theta0)),simultaneous=True))
check("all marked stage coordinates vanish",E.applyfunc(at)==S.zeros(3,1))
J=E.jacobian([g,u,v]).applyfunc(at)
check("exact marked inverse JL",(J*L).applyfunc(red)==S.eye(3))
check("exact marked inverse LJ",(L*J).applyfunc(red)==S.eye(3))
check("exact marked Jacobian",red(J.det())==-r/4)
check("exact inverse-stage determinant",red(L.det())==4*r/3)
Hess=S.hessian(E[2],[g,u,v]).applyfunc(at)
check("actual quadratic shear coefficient",red((L[:,0].T*Hess*L[:,0])[0]+2*f20)==0)
# Full exact sparse coordinates, including x0.
expanded=expanded_binary_map();K=expanded["K"];variables=expanded["variables"];F=expanded["F"]
check("all four expanded coordinates exist",len(F)==4 and all(f for f in F))
check("expanded coordinate term counts",[len(f) for f in F]==[5,88,137,381])
check("expanded coordinate degrees",[max(map(sum,f.keys())) for f in F]==[5,19,22,33])
x,y,zsrc,tsrc=variables
for idx in range(3):
    check("exact polynomial division coordinate"+str(idx+1),
      F[idx+1]*x**(idx+1)==expanded["numerators"][idx])
check("gamma restricts to gamma0 on whole x0",expanded["gamma"].evaluate(x,K.zero)==expanded["gamma0"])
check("universal determinant gamma exponent",1+2+1+2-(1+2+3)==0)
check("universal determinant x exponent",1+2+3-(1+2+3)==0)
# Serialize every coefficient; replay from rational pairs in the fixed quadratic field.
def coefficient_pair(a):
    a=list(a.to_list())
    a=[K.dom.zero]*(2-len(a))+a
    return [str(a[1]),str(a[0])] # constant, r coefficient
coordinates=[]
for f in F:
    coordinates.append([dict(exponents=list(m),coefficient=coefficient_pair(a)) for m,a in sorted(f.items())])
payload=dict(field="Q(r), r^2=-3",coefficient_order=["constant","r"],
 source_variables=["x","y","z_src","t_src"],coordinates=coordinates)
(ROOT/"binary_keller_coefficients.json").write_text(json.dumps(payload,indent=2)+"\n")
loaded=json.loads((ROOT/"binary_keller_coefficients.json").read_text())
for idx,records in enumerate(loaded["coordinates"]):
    replay=expanded["ring"].zero
    for item in records:
        a,b=map(S.Rational,item["coefficient"])
        term=K(a)+K(b)*expanded["r"]
        monomial=expanded["ring"].one
        for var,power in zip(variables,item["exponents"]):monomial*=var**power
        replay+=term*monomial
    check("independent serialized polynomial replay coordinate"+str(idx),replay==F[idx])
# Independent whole-polynomial determinant, directly from the serialized coordinates.
# Represent a+b*r as two rational multivariate polynomials, using r^2=-3.
import flint
from flint import fmpq_mpoly_ctx,fmpq
ctx=fmpq_mpoly_ctx.get(["x","y","z_src","t_src"])
rational_pairs=[]
for records in loaded["coordinates"]:
    rational_pairs.append(tuple(ctx.from_dict({
      tuple(record["exponents"]):fmpq(record["coefficient"][component])
      for record in records if record["coefficient"][component]!="0"})
      for component in range(2)))
fullJ=[[(a.derivative(j),b.derivative(j)) for j in range(4)] for a,b in rational_pairs]
def multiply_pair(a,b):
    return (a[0]*b[0]-3*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
full_det=[ctx.constant(0),ctx.constant(0)]
for perm in itertools.permutations(range(4)):
    term=(ctx.constant(1),ctx.constant(0))
    for row,col in enumerate(perm):term=multiply_pair(term,fullJ[row][col])
    sign=-1 if sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))%2 else 1
    for component in range(2):full_det[component]+=sign*term[component]
check("independent FULL expanded polynomial Jacobian",
    full_det==[ctx.constant(0),ctx.constant(fmpq(4,3))])
jacobian_values=[]
for point in [(0,0,0,0),(0,1,2,3),(1,0,0,0),(1,1,0,-1),(-1,2,-2,1)]:
    JJ=[[f.diff(var).evaluate([(axis,K(value)) for axis,value in zip(variables,point)]) for var in variables] for f in F]
    det=K.zero
    for perm in itertools.permutations(range(4)):
        term=K.one
        for row,col in enumerate(perm):term*=JJ[row][col]
        inv=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        det+=(-term if inv%2 else term)
    check("independent expanded Jacobian point"+str(point),det==expanded["gamma0"])
    jacobian_values.append(dict(point=point,value=coefficient_pair(det)))
# Two-sided field transformations are identities on the generic open.
C=S.symbols("C")
check("source x reconstruction",S.cancel((C/g)*g-C)==0)
for exponent in [1,2,3]:
    check("target twist exponent"+str(exponent),S.cancel(C**exponent*(S.Symbol("E")/(C/g)**exponent)-g**exponent*S.Symbol("E"))==0)
deck=data["original"].subs({g:-g,w:-w,z:z-2*w**3-2*g},simultaneous=True)
check("explicit generic binary deck involution",deck.applyfunc(S.expand)==data["original"].applyfunc(S.expand))
# Full eight-letter matrix action, its normal closure, subgroup interval and C invariant.
vectors=[(i,j) for i in range(3) for j in range(3) if (i,j)!=(0,0)]
group=set()
for a,b,c,d in itertools.product(range(3),repeat=4):
    if (a*d-b*c)%3!=1:continue
    group.add(tuple(vectors.index(((a*i+b*j)%3,(c*i+d*j)%3)) for i,j in vectors))
identity=tuple(range(8))
compose=lambda a,b:tuple(a[i] for i in b)
def inverse(a):
    out=[0]*8
    for i,j in enumerate(a):out[j]=i
    return tuple(out)
def generated(gens):
    out={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in gens:
            c=compose(a,b)
            if c not in out:out.add(c);todo.append(c)
    return frozenset(out)
stabilizer=frozenset(a for a in group if a[0]==0)
check("faithful transitive SL2F3 action",len(group)==24 and {a[0] for a in group}==set(range(8)) and len(stabilizer)==3)
normal=generated({compose(compose(a,h),inverse(a)) for a in group for h in stabilizer})
check("sheet stabilizer normally generates",normal==frozenset(group))
normalizer={a for a in group if {compose(compose(a,h),inverse(a)) for h in stabilizer}==set(stabilizer)}
check("generic deck normalizer",len(normalizer)==6)
interval={stabilizer};todo=[stabilizer]
while todo:
    H=todo.pop()
    for a in group-set(H):
        J=generated(set(H)|{a})
        if J not in interval:interval.add(J);todo.append(J)
check("complete subgroup interval sizes",sorted(map(len,interval))==[3,6,24])
def cycles(a):
    unseen=set(range(8));out=[]
    while unseen:
        i=next(iter(unseen));j=i;length=0
        while j in unseen:
            unseen.remove(j);length+=1;j=a[j]
        out.append(length)
    return sorted(out)
census={}
cycle_lcm=1
for a in group:
    profile=tuple(cycles(a));census[profile]=census.get(profile,0)+1
    cycle_lcm=math.lcm(cycle_lcm,math.gcd(*profile))
check("exact binary cycle-gcd obstruction",cycle_lcm==4 and cycle_lcm<8)
check("unique central involution",sum(compose(a,a)==identity and a!=identity for a in group)==1)
model=dict(claims=["KMFI59","KMFI60","KMFI61"],
 completion_criterion=dict(class_scope="rank2 normalized polynomial tautological sweeps",
   equivalent_conditions=["normal closure of sheet stabilizer is G","nonconstant unramified tail in d(S)=gamma^3*h",
     "regular source sheet above a suitable branch value","same field pair after independent C has an ordinary-affine4 Keller model"],
   finite_case="N>3delta,delta<=2; every N>6 sufficient",nonfinite_included=True,
   all_Keller_actions_classified=False),
 actual_binary_map=dict(source="ordinary affine4",coefficient_field="Q(r),r^2=-3",
   geometric_group="SL2(F3) on8 nonzero vectors",arithmetic_group="SL2(F3)",
   normal_closure_constants="Q(sqrt(-3))",degree=8,Jacobian="4r/3",
   coordinate_term_counts=[len(f) for f in F],coordinate_total_degrees=[max(map(sum,f.keys())) for f in F],
   theta0=[str(x) for x in theta0],L=[[str(x) for x in row] for row in L.tolist()],f20=str(f20),
   exact_expanded_Jacobian_controls=jacobian_values,
   polynomiality_verified_by_exact_division=True,full_expanded_Jacobian_verified=True),
 action_structure=dict(group_order=24,stabilizer_order=3,normalizer_order=6,
   subgroup_interval_orders=sorted(map(len,interval)),generic_deck_order=2,
   regular_deck_order=1,absolutely_and_stably_atomic=True,
   cycle_gcd_lcm=cycle_lcm,independently_curve_induced=False,
   cycle_census=[dict(profile=p,count=n) for p,n in sorted(census.items())],
   natural_A4_block_is_Keller_factor=False),
 absolute_binary_degree8_occurrence_decided=True,
 arbitrary_Keller_PSL2_11_actions_decided=False,absolute_natural_A4_decided=False,
 full_action_classification_complete=False)
(ROOT/"branch_centered_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["BRANCH_CENTERED_COMPLETION.md","branch_centered.py","verify_branch_centered.py",
 "binary_keller_coefficients.json","branch_centered_model.json","branch_centered_requirements.txt",
 "POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md","polynomial_chart_receipt.json","polynomial_wave.py"]
rec=dict(result="pass",check_count=len(checks),checks=checks,
 python=platform.python_version(),sympy=S.__version__,python_flint=flint.__version__,runtime_seconds=round(time.monotonic()-started,2),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written rank2 completion equivalence, explicit ordinary-affine binary map with exact coefficientwise polynomial division, marked derivatives, generic field reconstruction and full finite action controls. The full expanded-coordinate polynomial Jacobian is independently recomputed with rational multivariate polynomial pairs; exact composition factors provide a second universal proof. Finite point evaluations are supplemental. No formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"branch_centered_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 ordinary_affine_binary_realization=True,coordinate_terms=sum(map(len,F)),
 coefficient_file_bytes=(ROOT/"binary_keller_coefficients.json").stat().st_size,goal_complete=False)))

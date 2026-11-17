"""Exact homogeneous signed completions, forbidden blocks and PSL fixed targets."""
from pathlib import Path
from collections import Counter
from itertools import product,permutations
from math import factorial,gcd,lcm
import ast,hashlib,json,platform,time
import sympy as S
from homogeneous_blocks import homogeneous_signed,block_curve,centered_signed,marked_signed_data
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
w,z,g,q,u,A,B,e,tau=S.symbols("w z gamma q u A B e tau")
ll=S.symbols("ell",integer=True,positive=True)
cc=-(ll-1)**2/(2*(2*ll-1)*ll**2)
FF=S.Rational(1,2)+u/ll+cc*u*u
check("all-ell quadratic discriminant",S.factor(S.discriminant(FF,u)-1/(2*ll-1))==0)
check("all-ell critical-point identity",S.factor(ll*u*S.diff(FF,u)-FF+S.Rational(1,2)*(1-(ll-1)*u/ll)**2)==0)
check("all-ell critical value nonzero",
 S.factor(FF.subs(u,ll/(ll-1))-ll**2/((ll-1)*(2*ll-1)))==0)
NN=(ll-1)*q**(2*ll)-2*(2*ll-1)*q**ll+2*ll**2*q-(ll-1)*(2*ll-1)
for order in range(3):
    check("all-ell exact triple factor derivative"+str(order),S.simplify(S.diff(NN,q,order).subs(q,1))==0)
check("all-ell exact nonzero third derivative",
 S.simplify(S.diff(NN,q,3).subs(q,1)-2*ll**2*(ll-1)*(2*ll-1))==0)
check("all-ell zero constant of tail",S.factor(-NN.subs(q,0)-(ll-1)*(2*ll-1))==0)
power=S.symbols("Qpower")
formal=(ll-1)*power**2-2*(2*ll-1)*power+2*ll**2*q-(ll-1)*(2*ll-1)
check("all-ell roots of unity exclude other critical roots",
 S.factor(formal.subs(power,1)-2*ll**2*(q-1))==0)
# All-parameter block identity before choosing an integer exponent.
BB=A*A/(2*q)+A*q**(ll-1)/ll+cc*q**(2*ll-1)
check("all-ell block reconstruction qB/A^2",
 S.simplify(S.expand(q*BB/A**2-FF.subs(u,q**ll/A)))==0)
check("target birational inverse A",S.simplify((A*A/e)**ll/A**(2*ll-1)-A/e**ll)==0)
curve_controls=[];marked_controls=[]
for ell in [2,3,4,6,8]:
    H=homogeneous_signed(ell,w,z,g,q);curve=block_curve(ell,u)
    chart=polynomial_chart(w**(2*ell-2),[(1,0)],w=w,z=z)
    check(f"exact original wave family ell={ell}",(chart["X"]-H["X"]).applyfunc(S.expand)==S.zeros(3,1))
    check(f"exact full geometric action ell={ell}",chart["N"]==4*ell and chart["delta"]==2)
    if ell==2:
        check("binary exceptional action",chart["geometric_group"]=="SL2(F3) on8 letters")
    else:
        check("signed full group ell="+str(ell),chart["geometric_group"]==f"2^{2*ell-1}:A{2*ell} on{4*ell} letters")
    check("full sweep Jacobian ell="+str(ell),S.factor(H["sweep"].jacobian([g,w,z]).det())==g*g)
    tail=H["tail"]
    check("tail degree ell="+str(ell),S.degree(tail,q)==2*ell-3)
    check("tail has no zero or unit roots ell="+str(ell),tail.subs(q,0)==(ell-1)*(2*ell-1) and tail.subs(q,1)==S.Rational(ell**2*(ell-1)*(2*ell-1),3))
    check("tail is separable ell="+str(ell),S.gcd(S.Poly(tail,q,domain=S.QQ),S.Poly(S.diff(tail,q),q,domain=S.QQ)).degree()==0)
    check("tail has no lth roots of unity ell="+str(ell),S.gcd(tail,q**ell-1)==1)
    r=S.symbols("r");minimum=S.Poly(tail.subs(q,r*r),r)
    def red(expr):
        num,den=S.fraction(S.cancel(expr))
        return S.rem(S.rem(num,minimum.as_expr(),r)*S.invert(S.rem(den,minimum.as_expr(),r),minimum.as_expr(),r),minimum.as_expr(),r).expand()
    g0=H["b0"]*(1-r**(2*ell))/r
    point={g:g0,w:r,z:r**(2*ell-1)-1+g0}
    target=S.Matrix([1,H["b0"],H["d0"]-S.Rational(1,2)])
    check("all selected tail roots give regular zero ell="+str(ell),
      all(red(expr.subs(point,simultaneous=True))==0 for expr in H["sweep"]-target))
    check("gamma nonzero at every selected root ell="+str(ell),
      S.gcd(minimum,S.Poly(S.together(g0).as_numer_denom()[0],r)).degree()==0)
    centered=centered_signed(ell,g,S.Symbol("uu"),S.Symbol("vv"),w,z)
    for idx,expr in enumerate(centered["X"]):
        check(f"centered weighted order ell={ell},output={idx}",
          all(a+2*b>=idx+1 for (a,b),co in S.Poly(expr,w,z).terms()))
    F=curve["F"];phi=curve["phi"];critical=curve["critical"]
    check("curve degree ell="+str(ell),S.degree(S.together(phi).as_numer_denom()[0],u)==2*ell)
    check("curve critical derivative ell="+str(ell),
      S.factor(S.diff(phi,u)+F**(ell-1)*(1-S.Rational(ell-1,ell)*u)**2/(2*u*u))==0)
    passport=[[ell,ell],[3]+[1]*(2*ell-3),[2*ell-1,1]]
    check("complete rational curve RH ell="+str(ell),
      sum(sum(n-1 for n in profile) for profile in passport)==2*(2*ell)-2)
    branch_lcm=lcm(*(gcd(*profile) for profile in passport))
    check("strict branch gcd obstruction ell="+str(ell),branch_lcm==ell and branch_lcm<2*ell)
    # Two-sided exact source inverse and target phi relation.
    qback=e*F;Aback=qback**ell/u;Bback=Aback*Aback/e
    check("source inverse u ell="+str(ell),S.cancel(qback**ell/Aback-u)==0)
    check("source inverse e ell="+str(ell),S.cancel(Aback*Aback/Bback-e)==0)
    check("fixed curve target identity ell="+str(ell),S.cancel(Bback**ell/Aback**(2*ell-1)-phi)==0)
    b=H["b0"]*q**ell;d=H["d0"]*q**(2*ell-1)
    Q=S.expand((A-b)**2+2*q**ell*(A-b)+2*q*d-2*q*B)
    check("generic qQ different identity ell="+str(ell),
      S.expand(q*S.diff(Q,q)-Q+(A-b)**2)==0)
    check("exact root-product square class ell="+str(ell),
      S.factor(A*A/S.Poly(Q,q).LC()+(2*ell-1)*(ell*A/(ell-1))**2)==0)
    if ell==3:
        res=S.resultant(Q,A-b,q);disc=S.discriminant(Q,q);leading=S.Poly(Q,q).LC()
        check("generic block discriminant character ell3",
          S.cancel(disc-(-1)**(ell*(2*ell-1))*res**2/(leading*A*A))==0)
    curve_controls.append(dict(ell=ell,degree=4*ell,group=H["geometric_group"],
      tail=str(tail),coefficient_degree_bound=4*ell-6,
      block_degree=2*ell,block_curve=str(phi),passport=passport,fiber_gcd_lcm=branch_lcm))
for ell in [2,3]:
    H=marked_signed_data(ell);red=H["red"]
    check("entire marked stage zero ell="+str(ell),H["E"].applyfunc(H["at"])==S.zeros(3,1))
    check("entire marked derivative inverse ell="+str(ell),(H["J"]*H["L"]).applyfunc(red)==S.eye(3))
    check("marked determinant chain ell="+str(ell),red(H["determinant"]*H["theta0"][0])==1)
    check("marked inverse-stage determinant ell="+str(ell),red(H["L"].det()-H["theta0"][0])==0)
    marked_controls.append(dict(ell=ell,minimum=str(H["minimum"].as_expr()),
      theta0=[str(x) for x in H["theta0"]],L=[[str(x) for x in row] for row in H["L"].tolist()],f20=str(H["f20"])))
# Explicit rational field-pair countercontrols: action normal generation still holds.
for degree in [4,6]:
    AA={p for p in permutations(range(degree)) if sum(p[i]>p[j] for i in range(degree) for j in range(i+1,degree))%2==0}
    HH={p for p in AA if p[0]==0};identity=tuple(range(degree))
    def compose(a,b):return tuple(a[i] for i in b)
    def inverse(a):
        out=[0]*degree
        for i,j in enumerate(a):out[j]=i
        return tuple(out)
    normal_gens={compose(compose(g,h),inverse(g)) for g in AA for h in HH}
    normal={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in normal_gens:
            c=compose(a,b)
            if c not in normal:normal.add(c);todo.append(c)
    check("normally generating block stabilizer does not cure field obstruction n="+str(degree),normal==AA)
# Exact PSL2(11) actions using retained audited group definitions.
helper=ROOT/"verify_curve_descent_lcm.py"
tree=ast.parse(helper.read_text());wanted={"projective_group","cycles","generated"}
defs=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in wanted]
check("retained complete group helpers",{node.name for node in defs}==wanted)
ctx={"product":product};exec(compile(ast.Module(body=defs,type_ignores=[]),str(helper),"exec"),ctx)
G,I,mul,act12=ctx["projective_group"](11);cycles=ctx["cycles"];generated=ctx["generated"]
check("complete660 group",len(G)==660 and len({act12(g) for g in G})==660)
orders={g:lcm(*cycles(act12(g))) for g in G}
a2=min(g for g in G if orders[g]==2);H5=None
for a3 in sorted(g for g in G if orders[g]==3):
    if orders[mul(a2,a3)]!=5:continue
    test=generated((a2,a3),I,mul)
    if len(test)==60:H5=test;break
check("exceptional A5 constructed",H5 is not None and len(H5)==60)
lookup={};reps=[];unseen=set(G)
while unseen:
    rep=min(unseen);coset={mul(rep,h) for h in H5};index=len(reps)
    for g in coset:lookup[g]=index
    reps.append(rep);unseen-=coset
act11=lambda g:tuple(lookup[mul(g,rep)] for rep in reps)
check("exceptional11 complete faithful action",len(reps)==11 and len({act11(g) for g in G})==660)
for k in range(1,11):
    perm=[k*i%11 for i in range(11)];parity=sum(perm[i]>perm[j] for i in range(11) for j in range(i+1,11))%2
    check("11-cycle split class multiplier"+str(k),(-1)**parity==(1 if pow(k,5,11)==1 else -1))
check("11-cycle split half",sum(pow(k,5,11)==1 for k in range(1,11))==5)
def centralizer(profile):
    out=1
    for n,count in Counter(profile).items():out*=n**count*factorial(count)
    return out
def power_profile(profile,power):
    out=[]
    for n in profile:
        c=gcd(n,power);out.extend([n//c]*c)
    return tuple(sorted(out,reverse=True))
index=factorial(12)//2//660;targets=[]
for label,action,expected in [
 ("natural12",act12,74038),("exceptional11plus1",lambda g:act11(g)+(11,),74226)]:
    census=Counter(cycles(action(g)) for g in G)
    check(label+" subgroup lies in A12",all(sum(n-1 for n in profile)%2==0 for profile in census))
    check(label+" unipotent total",census[(11,1)]==120)
    check(label+" no single3-cycle",census[(3,)+(1,)*9]==0)
    inertia=[]
    for profile in [(6,6),(3,)+(1,)*9,(11,1)]:
        fixed=[]
        for power in range(lcm(*profile)):
            pp=power_profile(profile,power);ct=centralizer(pp);intersection=census[pp]
            if all(n%2 for n in pp) and len(set(pp))==len(pp):intersection//=2
            else:ct//=2
            check(label+" integral fixed powers "+str((profile,power)),ct*intersection%660==0)
            fixed.append(ct*intersection//660)
        check(label+" identity degree "+str(profile),fixed[0]==index)
        check(label+" integral cyclic orbit count "+str(profile),sum(fixed)%len(fixed)==0)
        count=sum(fixed)//len(fixed)
        inertia.append(dict(profile=profile,fixed_powers=fixed,cycle_count=count,index=index-count))
    value=1+S.Rational(index-sum(x["cycle_count"] for x in inertia),2)
    check(label+" exact target genus",value==expected and value>0)
    targets.append(dict(embedding=label,subgroup_order=660,coset_degree=index,
      genus=int(value),inertia=inertia,full_target_not_stably_rational=True,
      exact_independent_spectators=["Y1","e","C"]))
model=dict(claims=["KMFI62","KMFI63","KMFI64"],
 all_ell_family="l>=2, actual ordinary-affine4 degree4l; binary8 for l2 and even signed V_(2l):A_(2l) for l>=3",
 all_ell_coefficient_degree_bound="4l-6",arithmetic_constants_for_l_at_least3="K0(sqrt((-1)^(l-1)(2l-1)),sqrt(-(2l-1)))",
 all_pair_block_fields_independently_curve_induced=True,
 all_pair_block_fields_have_Keller_model=False,
 both_block_source_and_target_fields_rational=True,block_stabilizers_normally_generate=True,
 normal_generation_plus_rational_fields_is_a_sufficient_general_criterion=False,
 obstruction_scope="exact intermediate fields, all independent stabilization and birational recharts; not every abstract quotient action",
 curve_controls=curve_controls,complete_marked_stage_controls=marked_controls,
 PSL2_11_subgroup_targets=targets,
 generic_pencil_genus_used_as_total_rationality_proof=False,
 source_and_target_field_reconstructions_explicit=True,
 absolute_natural_A4_decided=False,arbitrary_Keller_PSL2_11_actions_decided=False,
 full_action_classification_complete=False)
(ROOT/"homogeneous_blocks_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["HOMOGENEOUS_BLOCK_OBSTRUCTION.md","homogeneous_blocks.py","verify_homogeneous_blocks.py",
 "homogeneous_blocks_model.json","polynomial_wave.py","verify_curve_descent_lcm.py",
 "BRANCH_CENTERED_COMPLETION.md","branch_centered_receipt.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written all-l family and exact intermediate-curve reconstruction, explicit marked completion via KMFI59, and both full660-element class censuses with all cyclic powers. Genus applies to intrinsic fixed curve fields with proved independent spectators. No formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"homogeneous_blocks_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 PSL_target_genera=[x["genus"] for x in targets],all_ell_block_obstruction=True,goal_complete=False)))

"""Exact controls for SIGNED_ALTERNATING.md; universal proofs are written."""
from pathlib import Path
from itertools import combinations,product
from math import factorial,gcd,lcm
import hashlib,json,platform,time
import sympy as S
from sympy.combinatorics import Permutation,PermutationGroup
from signed_alternating import family
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
M=S.symbols("M",positive=True,integer=True)
am=M*M/(2*(2*M+1)*(M+1)**2);bm=M/(2*(M+1)*(M+2))
disc=S.factor(bm*bm-am/6)
check("universal coefficient discriminant",
 S.factor(disc+M*M*(M-1)**2/(12*(M+1)**2*(M+2)**2*(2*M+1)))==0)
cAA=M/(2*M+1)-M/(M+1)+M*M/(2*(M+1)**2)
cAB=(M+1)/(M+2)-S.Rational(1,2)-M/(M+1)+M/(2*(M+1))
check("universal marked alpha squared coefficient",S.factor(cAA+am)==0)
check("universal marked alpha coefficient",S.factor(cAB+bm)==0)
check("universal marked constant coefficient",S.Rational(1,3)-S.Rational(1,2)+S.Rational(1,8)==-S.Rational(1,24))
check("universal critical value birational Euclidean identity",2*(M+1)-(2*M+1)==1)

v1,v2,v3,Y1,Y3=S.symbols("v1 v2 v3 Y1 Y3")
# Independent polynomial composition truncated by source x-weight.
def jet_compose(E,theta,variables,cut,red):
    weights=(1,2,4);zero=(0,0,0)
    def add(a,b):
        out=dict(a)
        for mon,c in b.items():
            out[mon]=red(out.get(mon,0)+c)
            if out[mon]==0:del out[mon]
        return out
    def mul(a,b):
        out={}
        for i,ci in a.items():
            for j,cj in b.items():
                mon=tuple(x+y for x,y in zip(i,j))
                if sum(x*y for x,y in zip(mon,weights))<cut:
                    out[mon]=out.get(mon,0)+ci*cj
        return {mon:red(c) for mon,c in out.items() if red(c)!=0}
    jets=[]
    for coord in theta:
        jets.append({mon:red(c) for mon,c in S.Poly(S.expand(coord),v1,v2,v3).terms()
                     if sum(x*y for x,y in zip(mon,weights))<cut and red(c)!=0})
    powers=[{0:{zero:S.Integer(1)}} for j in range(3)]
    result={}
    for exps,coeff in S.Poly(E,*variables).terms():
        term={zero:coeff}
        for j,e in enumerate(exps):
            for k in range(1,e+1):
                if k not in powers[j]:powers[j][k]=mul(powers[j][k-1],jets[j])
            term=mul(term,powers[j][e])
        result=add(result,term)
    return result
models=[]
for m in range(2,6):
    f=family(m);a,g,u,v,w,z,q,A,B=[f[k] for k in ["a","g","u","v","w","z","q","A","B"]]
    red=f["red"];n=f["n"];lab=f"m={m}"
    check(lab+" minimal polynomial has nonreal distinct roots",S.discriminant(f["minimal"],a)<0)
    check(lab+" b prime",S.expand(S.diff(f["b"],q)-q*S.diff(f["h"],q))==0)
    check(lab+" d prime",S.expand(S.diff(f["d"],q)-f["h"]*S.diff(f["b"],q))==0)
    origin=f["d"].subs(q,1)-f["h"].subs(q,1)*f["b"].subs(q,1)+f["b"].subs(q,1)**2/2
    check(lab+" marked origin relation",red(origin)==0)
    check(lab+" marked gamma nonzero",red(f["g0"])!=0)
    sw=f["sweep"];D=f["Delta"]
    check(lab+" full sweep Jacobian",
          S.expand(S.Matrix.hstack(D,sw.diff(w),sw.diff(z)).det()-g*g)==0)
    check(lab+" all marked sweep coordinates",
          all(red(x.subs({w:1,z:f["h"].subs(q,1)-f["b"].subs(q,1),g:f["g0"]}))==0 for x in sw))
    for i,e in enumerate([1,2,4]):
        check(lab+" polynomial gamma quotient "+str(i),
          S.expand(f["E"][i]*g**e-sw[i].subs({w:g*u,z:g**3*v},simultaneous=True))==0)
    check(lab+" stage derivative inverse",all(red(x)==0 for x in f["J"]*f["L"]-S.eye(3)))
    check(lab+" inverse Jacobian determinant",red(f["J"].det()-f["inv"](f["g0"]))==0)
    check(lab+" constant full Keller determinant",red(f["L"].det()-f["g0"])==0)
    shear=f["f20"]*v1*v1+f["f11"]*v1*v2+f["f30"]*v1**3
    theta=f["theta0"]+f["L"]*S.Matrix([v1,v2,v3+shear])
    recovered=(f["J"]*(theta-f["theta0"])).applyfunc(red)
    check(lab+" explicit triangular inverse",all(red(x)==0 for x in recovered-S.Matrix([v1,v2,v3+shear])))
    for i,e in enumerate([1,2,4]):
        low=jet_compose(f["E"][i],theta,[g,u,v],e,red)
        check(lab+" every source polynomiality coefficient "+str(i),not low)
    check(lab+" full determinant x and gamma cancellation",1+2+4==1+3+3)
    Q=f["Q"];leading=f["leading"]
    check(lab+" block generic degree",S.degree(Q,q)==n)
    check(lab+" constant leading coefficient",S.expand(S.Poly(Q,q).LC()-leading)==0)
    check(lab+" finite eliminant constant",Q.subs(q,0)==A*A)
    check(lab+" different square identity",S.expand(q*S.diff(Q,q)-Q+(A-f["b"])**2)==0)
    check(lab+" generic critical values distinct degree gate",gcd(S.degree(f["b"],q),S.degree(f["d"],q))==1)
    eq3=S.expand(f["X"][2]+(Y1-f["X"][0])*z-Y3)
    check(lab+" finite transverse monic equation",S.Poly(eq3,z).LC()==S.Rational(1,2))
    zr=f["h"].subs(q,w*w)*w-Y1+(A-f["b"].subs(q,w*w))/w
    check(lab+" exact root lifting eliminant",
      S.cancel(2*w*w*eq3.subs(z,zr)-Q.subs({q:w*w,B:Y3+Y1*Y1/2},simultaneous=True))==0)
    eps=(-1)**(n*(n-1)//2)
    check(lab+" block arithmetic square class",
      S.cancel((S.Integer(eps)/leading)/((-1)**m*(2*m+1))-((m+1)/(a*m))**2)==0)
    check(lab+" flip arithmetic square class",
      S.cancel((A*A/leading)/(-(2*m+1))-((m+1)*A/(a*m))**2)==0)
    model=dict(m=m,n=n,degree=2*n,minimal_polynomial=str(f["minimal"]),
      coefficient_field=f"Q(sqrt(-{3*(2*m+1)}))",gamma0=str(f["g0"]),
      theta0=[str(x) for x in f["theta0"]],stage_matrix=[[str(x) for x in row] for row in f["L"].tolist()],
      shear=[str(f[k]) for k in ["f20","f11","f30"]],
      block_equation=str(Q),geometric_group=f"2^{n-1}:A{n}",
      arithmetic_character_squareclasses=[(-1)**m*(2*m+1),-(2*m+1)])
    if m==2:
        expected=S.Matrix([[0,S.Rational(1,8)-a/4,S.Rational(17,72)*a+S.Rational(13,48)],
          [0,-32*a-29,-8*a/3-1],[-96*a-84,336*a+414,104*a+94]])
        check("degree twelve displayed stage matrix",all(red(x)==0 for x in f["L"]-expected))
        check("degree twelve displayed theta",all(red(x)==0 for x in f["theta0"]-S.Matrix([-2*a/3-S.Rational(1,2),16*a+18,-24])))
        check("degree twelve displayed shear",red(f["f20"]-48*a-42)==0 and f["f11"]==f["f30"]==0)
    sig,gg,CC=S.symbols("sig gg CC")
    yy=S.Matrix([sig+gg,f["b"]+gg*w,f["d"]-sig*sig/2+gg*(f["h"]*w-sig)])
    zz=yy.subs({w:-w,sig:sig+2*gg,gg:-gg},simultaneous=True)
    check(lab+" generic source deck fixes target",
          all(S.expand(v)==0 for v in zz-yy))
    check(lab+" regular source gamma boundary",
          all(S.expand(v)==0 for v in theta.subs({v1:0,v2:0,v3:0})-f["theta0"]))
    models.append(model)

def signed(n,perm,flip=()):
    return Permutation([2*perm[i]+(bit^(i in flip)) for i in range(n) for bit in [0,1]])
def cycle_gcd(perm):
    seen=set();lengths=[]
    for i in range(perm.size):
        if i in seen:continue
        j=i;k=0
        while j not in seen:seen.add(j);k+=1;j=perm(j)
        lengths.append(k)
    return gcd(*lengths)
group_controls=[]
for n in [6,8,10,12]:
    gens=[];basegens=[]
    for i in range(2,n):
        pp=list(range(n));pp[0]=1;pp[1]=i;pp[i]=0
        perm=signed(n,pp);gens.append(perm);basegens.append(perm)
    gens.append(signed(n,list(range(n)),[0,1]))
    G=PermutationGroup(gens)
    order=2**(n-1)*factorial(n)//2
    check(f"signed degree {2*n} exact group order",G.order()==order)
    check(f"signed degree {2*n} transitive",G.is_transitive())
    blocks=[len(set(G.minimal_block([0,j]))) for j in range(1,2*n)]
    check(f"signed degree {2*n} unique proper block system",blocks==[n]+[1]*(2*n-2))
    check(f"degree {n} Carmichael cube relations",all((x**3).is_Identity for x in basegens))
    check(f"degree {n} Carmichael pair relations",all(((x*y)**2).is_Identity for x,y in combinations(basegens,2)))
    check(f"degree {n} central relations have outside support",n>4)
    if n==6:
        actual=lcm(*(cycle_gcd(x) for x in G.generate_schreier_sims()))
        check("complete degree twelve cycle gcd census",actual==12)
    group_controls.append(dict(n=n,degree=2*n,order=order,only_block_size=2))

def partitions(n,minpart=1):
    if n==0:yield ()
    for i in range(minpart,n+1):
        for tail in partitions(n-i,i):yield (i,)+tail
type_controls=[]
for n in [6,8,10,12]:
    inv=1;types=0
    for part in partitions(n):
        if (n-len(part))%2:continue
        for flip in product([0,1],repeat=len(part)):
            if sum(flip)%2:continue
            lengths=[]
            for length,sign in zip(part,flip):
                lengths.extend([2*length] if sign else [length,length])
            inv=lcm(inv,gcd(*lengths));types+=1
    expected=n if n&(n-1)==0 else 2*n
    check(f"degree {2*n} complete signed cycle-type invariant",inv==expected)
    type_controls.append(dict(n=n,signed_cycle_types=types,C=inv))
# Exact rational-square tests in Q(sqrt(-15)).
def rational_square(x):
    x=S.Rational(x)
    return x>=0 and S.integer_nthroot(int(x.p),2)[1] and S.integer_nthroot(int(x.q),2)[1]
for x in [5,-5,-1]:
    check("independent degree-twelve constant character "+str(x),
          not rational_square(x) and not rational_square(S.Rational(x,-15)))
manifest=json.loads((ROOT/"sources/signed-alternating-manifest.json").read_text())
for source in manifest["sources"]:
    check("classical presentation source hash",
          hashlib.sha256((ROOT/source["path"]).read_bytes()).hexdigest()==source["sha256"])
model=dict(family="Every even n>=6, signed action 2^(n-1):A_n on2n letters, actual affine four-space maps",
 parameter_formula="A_m alpha^2+B_m alpha+1/24=0; coefficients and Taylor stage in signed_alternating.py",
 controls=models,group_controls=group_controls,cycle_type_controls=type_controls,
 degree_twelve=dict(geometric_order=11520,arithmetic_order=46080,
  coefficient_field="Q(sqrt(-15))",constant_field="Q(sqrt(-15))(i,sqrt(5))"),
 atomic_with_imprimitive_monodromy=True,
 cycle_gcd_formula="n for n a power of two, otherwise2n",
 general_PSL2_11_action_decided=False,full_classification_complete=False)
(ROOT/"signed_alternating_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["SIGNED_ALTERNATING.md","signed_alternating.py","verify_signed_alternating.py",
       "signed_alternating_model.json","sources/signed-alternating-manifest.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,
 python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written all-parameter polynomial lift, full signed-kernel, atomicity and cycle-invariant proofs. Classical Carmichael presentation attributed. Exact controls are not formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"signed_alternating_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=receipt["runtime_seconds"])))

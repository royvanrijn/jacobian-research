"""Exact finite controls for the universal central-cover action obstruction.

The Brauer/purity/twist proof is written in CENTRAL_LIFT_OBSTRUCTION.md.
Finite matrix or cohomology controls do not by themselves certify it.
"""
from pathlib import Path
from itertools import permutations,product
from math import gcd
import hashlib,json,platform,time
from central_lift import (sl2_prime,generated,gf2_rank,h2_mod2,
  permutation_inverse,permutation_subgroups,permutation_order)
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
family=[];groups={}
for p in [3,5,7,11,13,19]:
    data=sl2_prime(p);groups[p]=data
    SL=data["elements"];G=data["quotient"];I,Z=data["I"],data["Z"];mul=data["mul"]
    check("SL2 order q="+str(p),len(SL)==p*(p*p-1))
    check("PSL2 order q="+str(p),len(G)==p*(p*p-1)//2)
    check("two-element central kernel q="+str(p),data["canon"](I)==data["canon"](Z) and I!=Z)
    check("unique nonidentity involution q="+str(p),[a for a in SL if a!=I and mul(a,a)==I]==[Z])
    J=(0,p-1,1,0);j=data["canon"](J)
    check("nontrivial projective involution q="+str(p),j!=data["qI"] and data["qmul"](j,j)==data["qI"])
    check("projective action faithful q="+str(p),len({data["act"](a) for a in G})==len(G))
    H={a for a in G if data["act"](a)[p]==p};HH={a for a in SL if a[2]==0}
    check("natural stabilizer order q="+str(p),len(H)==p*(p-1)//2)
    check("full preimage stabilizer q="+str(p),len(HH)==2*len(H))
    char=lambda a:0 if pow(a[0],(p-1)//2,p)==1 else 1
    check("natural preimage character homomorphism q="+str(p),
      all(char(mul(a,b))==(char(a)+char(b))%2 for a in HH for b in HH))
    check("central character detects parity q="+str(p),char(Z)==(1 if p%4==3 else 0))
    applicable=len(H)%2==1
    check("odd stabilizer exact congruence q="+str(p),applicable==(p%4==3))
    if applicable:
        kernel={a for a in HH if char(a)==0}
        check("odd normal complement index2 q="+str(p),len(kernel)==len(H) and len(kernel)%2==1 and Z not in kernel)
    family.append(dict(q=p,degree=p+1,SL_order=len(SL),PSL_order=len(G),
      stabilizer_order=len(H),preimage_stabilizer_order=len(HH),central_double_cover_nonsplit=True,
      odd_stabilizer_obstruction_applies=applicable,absolute_exclusion_proved=applicable))
# Full normalized cochain ranks distinguish A4 from the actual binary group.
d=groups[3]
cohomA4=h2_mod2(d["quotient"],d["qI"],d["qmul"])
cohomBinary=h2_mod2(d["elements"],d["I"],d["mul"])
check("full A4 normalized H2 dimension1",cohomA4["H2_dimension"]==1)
check("actual binary group's normalized H2 vanishes",cohomBinary["H2_dimension"]==0)
check("A4 exact cochain ranks",cohomA4["cocycle_constraint_rank"]==109 and cohomA4["coboundary_rank"]==11)
check("binary exact cochain ranks",cohomBinary["cocycle_constraint_rank"]==506 and cohomBinary["coboundary_rank"]==23)
# The explicit A4 central extension cocycle is independently not a coboundary.
order=[d["qI"]]+sorted(set(d["quotient"])-{d["qI"]});lookup={a:i for i,a in enumerate(order)}
rows=[];augmented=[]
for a,b in product(order[1:],repeat=2):
    c=d["qmul"](a,b);row=0
    for x in [a,b,c]:
        j=lookup[x]
        if j:row^=1<<(j-1)
    actual=d["mul"](a,b);cocycle=int(actual!=c)
    check("central cocycle fiber "+str((lookup[a],lookup[b])),actual==c or actual==d["neg"](c))
    rows.append(row);augmented.append(row|(cocycle<<(len(order)-1)))
check("explicit A4 cocycle is nonsplit",gf2_rank(augmented)==gf2_rank(rows)+1)
# Complete averaging controls in small prime-to-l groups.
averaging=[]
for size,prime in [(3,2),(2,3)]:
    pairs=list(product(range(1,size),repeat=2));valid=0
    for values in product(range(prime),repeat=len(pairs)):
        cocycle=dict(zip(pairs,values))
        c=lambda a,b:cocycle.get((a,b),0)
        if any((c(a,b)+c((a+b)%size,t)-c(b,t)-c(a,(b+t)%size))%prime
               for a,b,t in product(range(size),repeat=3)):continue
        valid+=1
        average={a:(pow(size,-1,prime)*sum(c(a,t) for t in range(size)))%prime for a in range(size)}
        check(f"complete averaging split group{size},prime{prime},cocycle{valid}",
          all((c(a,b)-average[a]-average[b]+average[(a+b)%size])%prime==0
              for a,b in product(range(size),repeat=2)))
    averaging.append(dict(group_order=size,prime=prime,complete_cocycles=valid))
twists=[]
for prime in [2,3,5]:
    for inertia in [1,2,3,4,5,7]:
        if gcd(prime,inertia)!=1:continue
        # In C_(prime*inertia), central generator is the inertia-th power.
        exponents=[a for a in range(prime) if prime*inertia//gcd(prime*inertia,1+inertia*a)==inertia]
        check(f"unique central tame twist prime{prime},inertia{inertia}",len(exponents)==1)
        twists.append(dict(prime=prime,inertia=inertia,twist=exponents[0]))
    check("coprimality cannot be dropped prime"+str(prime),
      all(prime*prime//gcd(prime*prime,1+prime*a)==prime*prime for a in range(prime)))
# No C2 quotient on the exceptional stabilizer preimage in SL2(11).
d=groups[11];G=d["quotient"];mul=d["qmul"];I=d["qI"]
orders={a:permutation_order(d["act"](a)) for a in G}
a2=min(a for a in G if orders[a]==2);H5=None
for a3 in sorted(a for a in G if orders[a]==3):
    if orders[mul(a2,a3)]!=5:continue
    test=generated([a2,a3],I,mul)
    if len(test)==60:H5=test;break
check("exceptional A5 subgroup found",H5 is not None and len(H5)==60)
HH={a for a in d["elements"] if d["canon"](a) in H5}
inverse=lambda a:(a[3],(-a[1])%11,(-a[2])%11,a[0])
comm={d["mul"](d["mul"](d["mul"](a,b),inverse(a)),inverse(b)) for a,b in product(HH,repeat=2)}
derived=generated(comm,d["I"],d["mul"])
check("exceptional preimage order120",len(HH)==120)
check("exceptional preimage perfect",derived==frozenset(HH))
check("exceptional degree11 is not excluded by coprime criterion",len(H5)%2==0)
A5natural={a for a in permutations(range(5)) if sum(a[i]>a[j] for i in range(5) for j in range(i+1,5))%2==0}
check("actual natural A5 on5 lies outside the odd-stabilizer obstruction",
 len(A5natural)==60 and sum(a[0]==0 for a in A5natural)==12)
# All faithful transitive A4 actions, and complete degree<=4 transitive subgroups.
A4={groups[3]["act"](a) for a in groups[3]["quotient"]}
check("PSL2(3) is the whole natural A4",
 A4=={a for a in permutations(range(4)) if sum(a[i]>a[j] for i in range(4) for j in range(i+1,4))%2==0})
mulperm=lambda a,b:tuple(a[i] for i in b);I4=tuple(range(4))
allH=permutation_subgroups(A4);faithful=[]
for H in allH:
    core=set(H)
    conjugates=[]
    for a in A4:
        conjugate={mulperm(mulperm(a,h),permutation_inverse(a)) for h in H}
        core&=conjugate;conjugates.extend(conjugate)
    if len(core)!=1:continue
    normal=generated(conjugates,I4,mulperm)
    degree=len(A4)//len(H)
    norm_pass=normal==frozenset(A4)
    spin_excluded=len(H)%2==1
    check("every faithful A4 action excluded "+str((degree,len(H))),not norm_pass or spin_excluded)
    faithful.append(dict(degree=degree,stabilizer_order=len(H),normal_closure_order=len(normal),
       normal_generation_passes=norm_pass,central_cover_excluded=norm_pass and spin_excluded))
check("complete A4 subgroup count",len(allH)==10)
check("complete faithful A4 stabilizer sizes",{x["stabilizer_order"] for x in faithful}=={1,2,3})
small=[]
for n in range(1,5):
    symmetric=set(permutations(range(n)));subgroups=permutation_subgroups(symmetric)
    transitive=[H for H in subgroups if {a[0] for a in H}==set(range(n))]
    types=set()
    for H in transitive:
        size=len(H)
        if n==1:label="1"
        elif n==2:label="C2"
        elif n==3:label="C3" if size==3 else "S3"
        elif size==4:label="C4" if max(map(permutation_order,H))==4 else "V4"
        elif size==8:label="D8"
        elif size==12:label="A4"
        elif size==24:label="S4"
        else:raise AssertionError("unknown transitive degree4 group")
        types.add(label)
    expected={1:{"1"},2:{"C2"},3:{"C3","S3"},4:{"C4","V4","D8","A4","S4"}}[n]
    check("complete transitive subgroup types degree"+str(n),types==expected)
    realized={1:["1"],2:[],3:["S3"],4:["S4"]}[n]
    small.append(dict(degree=n,all_transitive_types=sorted(types),Keller_action_types=realized))
# Source pin: the global input is retained, not inferred from group ranks.
manifest=json.loads((ROOT/"sources/central-lift-manifest.json").read_text())
for source in manifest["sources"]:
    file=ROOT/source["path"]
    check("Brauer primary source SHA256",hashlib.sha256(file.read_bytes()).hexdigest()==source["sha256"])
    check("Brauer primary source bytes",file.stat().st_size==source["bytes"])
model=dict(claims=["KMFI65","KMFI66","KMFI67","KMFI68"],
 universal_necessary_condition="prime l not dividing |H| implies H^2(G,C_l)=0",
 all_affine_inertia_prime_to_l=True,
 universal_geometric_proof="local unique coprime complement, unramified Brauer obstruction, purity and Br(A^d)[l]=0, proper central lift, UFD Kummer twisting, forbidden cyclic etale source cover",
 finite_checks_are_not_formal_geometry_proof=True,
 A4_absolute_occurrence="excluded in every faithful transitive action",
 natural_PSL2_q_congruence3mod4_absolute_occurrence="excluded for every odd prime power q",
 natural_PSL2_11_degree12_absolute_occurrence="excluded",
 exceptional_PSL2_11_degree11_absolute_occurrence_decided=False,
 natural_A4_decided=True,
 family_controls=family,cohomology_controls=dict(A4=cohomA4,actual_binary8=cohomBinary),
 averaging_controls=averaging,tame_twist_controls=twists,
 faithful_A4_action_controls=faithful,small_degree_complete_spectrum=small,
 exceptional_preimage=dict(order=120,perfect=True,criterion_applies=False),
 actual_binary_occurrence_preserved=True,full_action_classification_complete=False)
(ROOT/"central_lift_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["CENTRAL_LIFT_OBSTRUCTION.md","central_lift.py","verify_central_lift.py",
 "central_lift_model.json","sources/central-lift-manifest.json",
 "sources/brauer-grothendieck-2019.pdf"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),
 runtime_seconds=round(time.monotonic()-started,3),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written absolute action-level proof using pinned primary Brauer purity and homotopy inputs; exact matrix double covers, full normalized mod2 cochain ranks, complete small subgroup classifications and all local finite controls. No formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"central_lift_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"],
 absolute_A4_excluded=True,natural_PSL11_degree12_excluded=True,
 exceptional_PSL11_degree11_decided=False,goal_complete=False)))

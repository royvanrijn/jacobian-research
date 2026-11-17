"""Exact PSL2(11) subgroup-resolvent genus controls for SUBGROUP_REDUCTION.md."""
from pathlib import Path
from collections import Counter
from itertools import product
from math import factorial,gcd,lcm
import ast,hashlib,json,platform,time

ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)

# Load the retained pure group helper definitions without executing its checker.
helper=ROOT/"verify_curve_descent_lcm.py"
tree=ast.parse(helper.read_text())
names={"projective_group","cycles","generated"}
defs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
check("retained helper definitions found",{n.name for n in defs}==names)
ctx={"product":product}
exec(compile(ast.Module(body=defs,type_ignores=[]),str(helper),"exec"),ctx)
G,I,mul,act12=ctx["projective_group"](11)
cycles=ctx["cycles"];generated=ctx["generated"]
check("complete PSL2(11) matrix group",len(G)==660)
orders={g:lcm(*cycles(act12(g))) for g in G}
a2=min(g for g in G if orders[g]==2)
H=None
for a3 in sorted(g for g in G if orders[g]==3):
    if orders[mul(a2,a3)]!=5:continue
    candidate=generated((a2,a3),I,mul)
    if len(candidate)==60:
        H=candidate;break
check("exceptional A5 point stabilizer constructed",H is not None and len(H)==60)
lookup={};reps=[];unseen=set(G)
while unseen:
    r=min(unseen);coset={mul(r,h) for h in H}
    for g in coset:lookup[g]=len(reps)
    reps.append(r);unseen-=coset
act11=lambda g:tuple(lookup[mul(g,r)] for r in reps)
check("exceptional action has degree eleven",len(reps)==11)
B={g for g in G if act12(g)[11]==11}
check("natural degree-twelve stabilizer order",len(B)==55)
census11=Counter(cycles(act11(g)) for g in G)
census12=Counter(cycles(act12(g)) for g in G)
censusH=Counter(cycles(act11(g)) for g in H)
censusB=Counter(cycles(act12(g)) for g in B)
check("exceptional order-five class intersection",census11[(5,5,1)]==264)
check("exceptional order-eleven class total",census11[(11,)]==120)
check("exceptional involution profile",census11[(2,2,2,2,1,1,1)]==55)
check("exceptional order-three profile",census11[(3,3,3,1,1)]==110)
check("exceptional subgroup contains no single three-cycle",census11[(3,)+(1,)*8]==0)
check("A5 order-five class intersection",censusH[(5,5,1)]==24)
check("B55 unipotent class total",censusB[(11,1)]==10)
check("natural action order-three profile",census12[(3,3,3,3)]==110)
check("natural subgroup contains no single three-cycle",census12[(3,)+(1,)*9]==0)
# A cyclic order-eleven subgroup contains five powers in each split alternating class.
# Verify the sign of a multiplier on the eleven cyclic labels equals the Legendre symbol.
for multiplier in range(1,11):
    perm=[multiplier*i%11 for i in range(11)]
    parity=sum(perm[i]>perm[j] for i in range(11) for j in range(i+1,11))%2
    legendre=pow(multiplier,5,11)
    check("eleven-cycle split-class multiplier "+str(multiplier),
          (-1 if parity else 1)==(1 if legendre==1 else -1))
check("every unipotent subgroup splits its ten generators equally",
      sum(pow(i,5,11)==1 for i in range(1,11))==5)

def centralizer_sn(profile):
    out=1
    for length,count in Counter(profile).items():out*=length**count*factorial(count)
    return out
def split_in_an(profile):
    return all(v%2 for v in profile) and len(set(profile))==len(profile)
def power_profile(profile,power):
    parts=[]
    for length in profile:
        count=gcd(length,power)
        parts.extend([length//count]*count)
    return tuple(sorted(parts,reverse=True))
def coset_fixed(profile,ambient,subgroup_order,census):
    intersections=census[profile]
    centralizer=centralizer_sn(profile)
    if ambient=="A":
        if split_in_an(profile):
            # Only unipotent split classes have nonzero intersections in our cases.
            if intersections:
                check("supported split-class profile "+str(profile),profile in [(11,),(11,1)])
                check("split class total even "+str(profile),intersections%2==0)
                intersections//=2
        else:centralizer//=2
    numerator=centralizer*intersections
    check("coset fixed count integral "+str((ambient,profile,subgroup_order)),numerator%subgroup_order==0)
    return numerator//subgroup_order
def inertia_data(profile,ambient,n,subgroup_order,census):
    group_order=factorial(n)//(2 if ambient=="A" else 1)
    degree=group_order//subgroup_order
    order=lcm(*profile)
    fixed=[coset_fixed(power_profile(profile,k),ambient,subgroup_order,census) for k in range(order)]
    check("identity fixes whole coset space "+str((n,profile,subgroup_order)),fixed[0]==degree)
    check("Burnside cyclic orbit count integral "+str((n,profile,subgroup_order)),sum(fixed)%order==0)
    orbits=sum(fixed)//order
    # Recover full orbit-length profile independently of the Burnside average.
    exact_points={};orbit_counts={}
    for divisor in (d for d in range(1,order+1) if order%d==0):
        points=coset_fixed(power_profile(profile,divisor),ambient,subgroup_order,census)
        points-=sum(v for d,v in exact_points.items() if divisor%d==0)
        check("nonnegative integral inertia cycles "+str((n,profile,subgroup_order,divisor)),
              points>=0 and points%divisor==0)
        exact_points[divisor]=points
        if points:orbit_counts[divisor]=points//divisor
    check("inertia profile accounts for all sheets "+str((n,profile,subgroup_order)),
          sum(k*v for k,v in orbit_counts.items())==degree and sum(orbit_counts.values())==orbits)
    return dict(original_profile=profile,order=order,fixed_powers=fixed,
                coset_cycles=orbit_counts,index=degree-orbits)
def quotient_case(name,ambient,n,profiles,expected_target,expected_source,scope):
    censusJ,censusP=(census11,censusH) if n==11 else (census12,censusB)
    point_order=60 if n==11 else 55
    records={}
    for label,order,census,expected in [
       ("target",660,censusJ,expected_target),("source",point_order,censusP,expected_source)]:
        degree=factorial(n)//(2 if ambient=="A" else 1)//order
        rows=[inertia_data(profile,ambient,n,order,census) for profile in profiles]
        two_g_minus_two=-2*degree+sum(r["index"] for r in rows)
        check(name+" "+label+" Riemann-Hurwitz integral",two_g_minus_two%2==0)
        genus=(two_g_minus_two+2)//2
        check(name+" "+label+" exact genus",genus==expected and genus>0)
        records[label]=dict(subgroup_order=order,coset_degree=degree,
                            genus=genus,inertia=rows)
    relative_ram=(2*expected_source-2)-n*(2*expected_target-2)
    check(name+" relative ramification nonnegative",relative_ram>=0 and relative_ram%2==0)
    return dict(name=name,ambient=ambient+str(n),PSL_group_order=660,
                reduced_point_degree=n,scope=scope,relative_ramification_degree=relative_ram,**records)

cases=[]
cases.append(quotient_case(
 "A11_higher_power_radical","A",11,[(5,5,1),(3,)+(1,)*8,(11,)],
 5678,62493,"Fixed curve field; both target and source are not stably rational."))
cases.append(quotient_case(
 "S11_quadratic_beta","S",11,[(2,)*5+(1,),(6,)+(1,)*5,(11,)],
 7331,80641,"Fixed curve field after explicit two-parameter birational reconstruction; no affine target."))
cases.append(quotient_case(
 "S11_degree_ten_radical","S",11,[(10,1),(2,)+(1,)*9,(11,)],
 9343,102813,"Fixed curve field; no affine target."))
cases.append(quotient_case(
 "A12_level_two_cone","A",12,[(3,)+(1,)*9]*6+[(11,1)],
 527826,6333906,"Generic-pencil genera only; absolute rationality of the total surface is not decided."))
cases.append(quotient_case(
 "A12_alternating_shear","A",12,[(3,)+(1,)*9]*6+[(9,3)],
 524161,6289921,"Generic-pencil genera only; relative reduced twelve-sheet curve is unramified; total surface rationality open."))
check("A11 reduced cover ramification control",cases[0]["relative_ramification_degree"]==90)
check("A12 shear relative curve is etale",cases[-1]["relative_ramification_degree"]==0)
data=dict(cases=cases,
 subgroup_censuses={k:{str(p):v for p,v in c.items()} for k,c in
                    [("PSL11_exceptional",census11),("A5_exceptional_point",censusH),
                     ("PSL11_natural",census12),("B55_natural_point",censusB)]},
 fixed_curve_subgroup_reduction_to_affine_space="excluded for all three degree-eleven inputs",
 generic_pencil_genus_implies_total_nonrationality=False,
 PSL2_11_affine_Keller_map_constructed=False,full_classification_complete=False)
(ROOT/"subgroup_reduction_model.json").write_text(json.dumps(data,indent=2)+"\n")
files=["SUBGROUP_REDUCTION.md","verify_subgroup_reduction.py","subgroup_reduction_model.json",
       "verify_curve_descent_lcm.py"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),
 runtime_seconds=round(time.monotonic()-started,3),
 sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
 assurance="Complete 660-element actions, exact class intersections and all cyclic fixed powers; written genus and stable-rationality proof. No formal verification or independent whole-proof review.",
 fixed_curve_subgroup_reduction_excluded=True,general_PSL2_11_action_excluded=False,
 A12_total_surface_rationality_decided=False,full_classification_complete=False)
(ROOT/"subgroup_reduction_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),
 target_genera={c["name"]:c["target"]["genus"] for c in cases},
 runtime_seconds=receipt["runtime_seconds"],full_classification_complete=False)))

"""Exact finite controls for ORBIT_FORMS.md; the module proof is written there."""
from pathlib import Path
from itertools import combinations,product
from math import comb,factorial
from collections import Counter
import ast,json,hashlib,platform,re,time
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
# All m in the written module theorem satisfy every KMFI19 collision bound.
modules=[]
for family,start,canonical,extra in [("cone",2,18,0),("shear",3,12,8)]:
    for m in range(start,7):
        cross=m*(12-m);projective=cross-canonical-extra
        check(f"{family} m={m} projective degree",projective>=0)
        local=[("triple",3,1,2,0),("five-cusp",5,1,6,0)]
        if family=="cone":local.append(("origin",4,1,4,0))
        else:local.extend([("origin",4,3,4,8),("lost",3,2,0,4),
                          ("infinity-seven",7,4,30,3*projective)])
        bounds=[]
        for label,size,weight,volume,extra_order in local:
            cuts=[k*(size-k) for k in range(max(0,size-(12-m)),min(size,m)+1)]
            bound=max(cuts)*weight
            check(f"{family} m={m} {label} full cut bound",volume+extra_order>=bound)
            bounds.append(dict(place=label,required=bound,available=volume+extra_order))
        modules.append(dict(family=family,m=m,dimension=comb(12,m)-comb(12,m-1),
                            projective_factor=projective,local_bounds=bounds))
check("cone guaranteed form-space dimension",sum(x["dimension"] for x in modules if x["family"]=="cone")==912)
check("shear guaranteed form-space dimension",sum(x["dimension"] for x in modules if x["family"]=="shear")==858)
# Lagrange incidence relations at a fully split exact root fiber.
prime=991;omega=118
roots=[pow(omega,i,prime) for i in range(11)]+[0]
for m in range(2,7):
    sums={B:0 for B in combinations(range(12),m-1)}
    for A in combinations(range(12),m):
        den=1
        for i in A:
            for j in set(range(12))-set(A):den=den*(roots[i]-roots[j])%prime
        value=pow(den,-1,prime)
        for k in range(m):sums[A[:k]+A[k+1:]]=(sums[A[:k]+A[k+1:]]+value)%prime
    check(f"all Lagrange incidence relations m={m}",all(v==0 for v in sums.values()))

def cycles_perm(text):
    p=list(range(12))
    for body in re.findall(r"\(([\d,]+)\)",text):
        c=[int(v)-1 for v in body.split(",")]
        for a,b in zip(c,c[1:]+c[:1]):p[a]=b
    return tuple(p)
def read_gens(line):
    text=re.search(r",\[(\(.*)\]\],?$",line).group(1)
    depth=0;start=0;out=[]
    for i,c in enumerate(text):
        if c=="(":depth+=1
        if c==")":depth-=1
        if c=="," and depth==0:
            out.append(cycles_perm(text[start:i]));start=i+1
    out.append(cycles_perm(text[start:]))
    return out
def orbit_counts(gens):
    counts=[]
    for m in range(1,7):
        unseen=set(combinations(range(12),m));count=0
        while unseen:
            root=min(unseen);orbit={root};stack=[root]
            while stack:
                c=stack.pop()
                for g in gens:
                    d=tuple(sorted(g[i] for i in c))
                    if d not in orbit:orbit.add(d);stack.append(d)
            unseen-=orbit;count+=1
        counts.append(count)
    return counts
catalogue=(ROOT/"sources/primgrp-degree12.g").read_text()
# Exact retained classification is six entries, two explicit and four constructor records.
check("complete degree-twelve primitive excerpt",
      sum(line.startswith("[") for line in catalogue.splitlines())==6
      and all(token in catalogue for token in ['"M(11)"','"M(12)"','["PSL",2,11]','["PGL",2,11]','["Alt"]','["Sym"]']))
rows=[]
for line in catalogue.splitlines():
    if not (line.startswith("[1,") or line.startswith("[2,")):continue
    name=re.findall(r'"M\(\d+\)"',line)[0].strip('"')
    counts=orbit_counts(read_gens(line))
    expected=[1,1,1,2,2,3] if name=="M(11)" else [1,1,1,1,1,2]
    check(name+" full subset orbits",counts==expected)
    rows.append(dict(group=name,subset_orbits=counts,cone_pg_lower_bound=counts[5]-counts[0],
                     shear_pg_lower_bound=counts[5]-counts[1]))
# Retained PSL dictionary and independent orbit BFS using all exact group elements.
helper=ROOT/"verify_curve_descent_lcm.py";tree=ast.parse(helper.read_text());ctx={"product":product}
defs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ["projective_group","cycles","generated"]]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(helper),"exec"),ctx)
G,I,mul,act=ctx["projective_group"](11)
# Use all 660 permutations: tiny subset space, no large alternating enumeration.
counts=orbit_counts([act(g) for g in G])
check("PSL2(11) orbit census",counts==[1,1,1,2,2,6])
rows.append(dict(group="PSL2(11)",subset_orbits=counts,cone_pg_lower_bound=5,shear_pg_lower_bound=5))

# Imprimitive full block stabilizers. Their invariants also belong to any even subgroup.
wreaths=[]
for b in (2,3,4,6):
    nblocks=12//b
    a=list(range(12));a[0],a[1]=1,0
    rotate=list(range(12))
    for i in range(b):rotate[i]=(i+1)%b
    swap=list(range(12))
    for i in range(b):swap[i],swap[b+i]=b+i,i
    blockrotate=tuple(((i//b+1)%nblocks)*b+i%b for i in range(12))
    counts=orbit_counts([tuple(a),tuple(rotate),tuple(swap),blockrotate])
    check("block stabilizer has two pair orbits b="+str(b),counts[1]==2)
    check("block stabilizer positive high module b="+str(b),counts[5]>counts[1])
    wreaths.append(dict(block_size=b,subset_orbits=counts,form_lower_bound=counts[5]-counts[1]))

# Stabilizers for fields of unordered subsets and ordered root configurations.
configuration=[]
for k in range(1,7):
    r=lambda m:min(k,m)+1
    cone=r(6)-r(1);shear=r(6)-r(2)
    check("unordered subset source bounds k="+str(k),cone==k-1 and shear==max(0,k-2))
    configuration.append(dict(action="unordered subsets",k=k,degree=comb(12,k),
                              cone_pg_lower_bound=cone,shear_pg_lower_bound=shear))
    cone=2**k-(k+1);shear=2**k-1-k-comb(k,2)
    check("ordered configuration source bounds k="+str(k),cone>=0 and shear>=0)
    configuration.append(dict(action="ordered distinct roots",k=k,degree=factorial(12)//factorial(12-k),
                              cone_pg_lower_bound=cone,shear_pg_lower_bound=shear))
manifest=json.loads((ROOT/"sources/orbit-forms-manifest.json").read_text())
for source in manifest["sources"]:
    check("primary source hash "+source["path"],
          hashlib.sha256((ROOT/"sources"/source["path"]).read_bytes()).hexdigest()==source["sha256"])
raw=(ROOT/"sources/primgrp-gps1.g").read_bytes()
check("official primitive catalogue blob SHA",
      hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()=="78fe85c09f180028f3752535d044d65132e145a8")
check("catalogue completeness retained", "all primitive" in (ROOT/"sources/primgrp-README.md").read_text())
model=dict(modules=modules,subgroup_orbit_controls=rows,imprimitive_block_controls=wreaths,
 configuration_actions=configuration,
 all_proper_transitive_A12_subgroup_targets_of_both_audited_families_not_stably_rational=True,
 classification_inputs="Beaumont--Peterson Theorem 10 and official complete degree-twelve PrimGrp classification.",
 general_abstract_group_actions_excluded=False,full_classification_complete=False)
(ROOT/"orbit_forms_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["ORBIT_FORMS.md","verify_orbit_forms.py","orbit_forms_model.json",
       "SURFACE_FORMS.md","verify_curve_descent_lcm.py","sources/orbit-forms-manifest.json"]
receipt=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),
 runtime_seconds=round(time.monotonic()-started,3),
 sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
 assurance="Written regularity, incidence-module irreducibility and quotient proofs. Exact orbit controls; established external finite-group classification is attributed, not independently reproved. No formal verification or independent whole-proof review.",
 full_abstract_action_classification_complete=False)
(ROOT/"orbit_forms_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),runtime_seconds=receipt["runtime_seconds"],
                     subgroup_pg_bounds=rows,complete_classification=False)))

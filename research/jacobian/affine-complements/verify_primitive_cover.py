"""Finite controls for PRIMITIVE_QUARTIC_FIELD.md; purity is imported.

The normalization/DVR/purity proof is written, not formally certified.
All subgroups below are literal subgroups of S4, not conjugacy classes.
"""
from pathlib import Path
from itertools import permutations,combinations_with_replacement
from collections import Counter
from time import monotonic
import hashlib,json,platform

HERE=Path(__file__).resolve().parent
start=monotonic()
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)

# In a proper quartic tower the inner degree is two. These are every
# complete finite-DVR decomposition profile of that degree.
ef=[(e,f) for e in range(1,3) for f in range(1,3) if e*f<=2]
profiles=[]
for length in range(1,3):
    profiles.extend(p for p in combinations_with_replacement(ef,length)
                    if sum(e*f for e,f in p)==2)
check("complete quadratic DVR profile list",set(profiles)=={
 ((1,1),(1,1)),((1,2),),((2,1),)})
for p in profiles:
    composite=[(2*e,f) for e,f in p]
    check("no quadratic inner profile gives sole ramified e2 f1 prime "+str(p),
          composite!=[(2,1)])

points=range(4)
perms=tuple(permutations(points))
one=tuple(points)
def mul(p,q):return tuple(p[q[i]] for i in points)
def close(gens):
    gens=tuple(gens)
    found={one};todo=[one]
    while todo:
        p=todo.pop()
        for q in gens:
            product=mul(p,q)
            if product not in found:
                found.add(product);todo.append(product)
    return frozenset(found)
subgroups={frozenset({one})};todo=list(subgroups)
while todo:
    h=todo.pop()
    for g in perms:
        if g in h:continue
        k=close(tuple(h)+(g,))
        if k not in subgroups:
            subgroups.add(k);todo.append(k)
check("all literal S4 subgroups enumerated",len(subgroups)==30)
check("subgroup order census",dict(Counter(map(len,subgroups)))==
      {1:1,2:9,3:4,4:7,6:4,8:3,12:1,24:1})
pair_partitions={frozenset((frozenset((0,j)),frozenset(set(points)-{0,j})))
                 for j in (1,2,3)}
check("all possible nontrivial transitive block systems",len(pair_partitions)==3)
def transitive(h):return {g[0] for g in h}==set(points)
def invariant(h,partition):
    return all(frozenset(frozenset(g[i] for i in block) for block in partition)==partition
               for g in h)
def primitive(h):
    return transitive(h) and not any(invariant(h,p) for p in pair_partitions)
def is_transposition(g):return sum(g[i]!=i for i in points)==2
transitives=[h for h in subgroups if transitive(h)]
primitives=[h for h in transitives if primitive(h)]
with_transposition=[h for h in primitives if any(is_transposition(g) for g in h)]
check("nine literal transitive subgroups",len(transitives)==9)
check("two primitive groups: A4 and S4",sorted(map(len,primitives))==[12,24])
check("primitive plus transposition forces S4",len(with_transposition)==1 and
      len(with_transposition[0])==24)
check("transposition is six distinct conjugate edges",
      sum(is_transposition(g) for g in perms)==6)

def even(g):
    return sum(g[i]>g[j] for i in points for j in points if i<j)%2==0
alternating=frozenset(g for g in perms if even(g))
stabilizer=frozenset(g for g in perms if g[0]==0)
check("sign kernel has index two",len(alternating)==12)
check("quartic stabilizer has order six",len(stabilizer)==6)
check("quartic stabilizer not contained in sign kernel",not stabilizer<=alternating)
check("no subgroup strictly between quartic stabilizer and S4",
      not any(stabilizer<h and len(h)<24 for h in subgroups))
check("quadratic field occurs in closure but not in quartic field",
      alternating in subgroups and not stabilizer<=alternating)
manifest_path=HERE/"sources/primitive-cover-manifest.json"
manifest=json.loads(manifest_path.read_text())
for source in manifest["sources"]:
    data=(manifest_path.parent/source["file"]).read_bytes()
    check("purity source bytes "+source["file"],len(data)==source["bytes"])
    check("purity source sha256 "+source["file"],hashlib.sha256(data).hexdigest()==source["sha256"])
check("purity remains an external input",manifest["purity_theorem_is_external_input"] is True and
      manifest["whole_purity_proof_independently_reproved"] is False)
report=dict(as_of="2026-11-14",outcome="PASS",check_count=len(checks),checks=checks,
            quadratic_inner_dvr_profiles=profiles,
            s4_subgroups=dict(total=len(subgroups),transitive=len(transitives),
                             primitive_orders=sorted(map(len,primitives)),
                             primitive_with_transposition_orders=[len(h) for h in with_transposition]),
            written_result="Both actual quartic normalization fields have no proper intermediate field and have S4 Galois closure",
            assurance=dict(normalization_and_all_degree_dvr_argument="written proof",
                           purity="external Stacks 0BMB input, pinned source",
                           sign_field_in_galois_closure_excluded=False,
                           actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
                           independent_full_reproof=False,formal_verification=False,external_review=False),
            source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in [Path(__file__),HERE/"PRIMITIVE_QUARTIC_FIELD.md",manifest_path]},
            python=platform.python_version(),elapsed_seconds=round(monotonic()-start,6))
(HERE/"primitive_cover_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
print(f"PASS {len(checks)} finite controls; quartic field primitive, S4; purity imported")

#!/usr/bin/env python3
"""Literal ternary affine quotient certificates, with no cover-group input.

Standard library only. Every matrix and kernel word is checked directly.
Published bridge=Wirtinger and retained named-census identity remain imports.
"""
from collections import defaultdict
from copy import deepcopy
from itertools import product
from pathlib import Path
import hashlib,json,platform

HERE=Path(__file__).resolve().parent
INPUT=HERE/"ternary_affine_quotient_inputs.json"
NAMES=["13a_373","13a_375","13a_391","13a_392","13a_3264","13a_3266","13a_3275","13a_3276",
       "13n_1005","13n_1006","13n_1007","13n_1011","13n_1012","13n_1013","13n_1043","13n_1044",
       "13n_1045","13n_1046","13n_1047","13n_1048","13n_3842","13n_3844","13n_3855","13n_3856"]
ID=tuple(range(7))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rank(rows,n):
    a=[[v%3 for v in row] for row in rows];r=0
    for j in range(n):
        i=next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None:continue
        a[r],a[i]=a[i],a[r];c=pow(a[r][j],-1,3);a[r]=[c*v%3 for v in a[r]]
        for i in range(r+1,len(a)):
            if a[i][j]:
                c=a[i][j];a[i]=[(v-c*w)%3 for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r
def compose(p,q):return tuple(p[q[i]] for i in range(7))
def act(p,v):
    ans=[0]*7
    for i,a in enumerate(v):ans[p[i]]=a
    return ans
def affine_mul(a,b):
    p,v=a;q,w=b
    return compose(p,q),[(x+y)%3 for x,y in zip(v,act(p,w))]
def evaluate(word,gens):
    ans=(ID,[0]*7)
    for j in word:
        assert 0<=j<4
        ans=affine_mul(ans,gens[j])
    return ans
def homogeneous(p,v):
    return [[(int(p[j]==i)-int(p[6]==i))%3 for j in range(6)]+[v[i]] for i in range(6)]+[[0]*6+[1]]
def multiply(a,b):return [[sum(x*y for x,y in zip(row,col))%3 for col in zip(*b)] for row in a]
def difference(a):return [[(x-int(i==j))%3 for j,x in enumerate(row)] for i,row in enumerate(a)]
def matvec(rows,v):return sum(((a&v).bit_count()%2)<<i for i,a in enumerate(rows))
def gl3():
    ans={tuple(matvec(rows,v)-1 for v in range(1,8))
         for rows in product(range(1,8),repeat=3)
         if sorted(matvec(rows,v) for v in range(1,8))==list(range(1,8))}
    assert len(ans)==168
    return ans
def projection_group(gens):
    seen={ID};queue=[ID]
    for p in queue:
        for q in gens:
            r=compose(q,p)
            if r not in seen:seen.add(r);queue.append(r)
    return seen
def diagram(pd):
    ends=defaultdict(list);adj=defaultdict(set)
    for i,row in enumerate(pd):
        assert len(row)==4 and len(set(row))==4
        for h,e in enumerate(row):ends[e].append(4*i+h)
        adj[row[1]].add(row[3]);adj[row[3]].add(row[1])
    assert len(ends)==2*len(pd) and all(len(a)==2 for a in ends.values())
    twin={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
    seen=set();v=0
    for _ in range(2*len(pd)):
        assert v not in seen;seen.add(v)
        v=twin[4*(v//4)+(v+2)%4]
    assert v==0 and len({v//4 for v in seen})==len(pd)
    unseen=set(twin);faces=0
    while unseen:
        start=v=min(unseen)
        while True:
            unseen.remove(v);w=twin[v];v=4*(w//4)+(w+1)%4
            if v==start:break
        faces+=1
    assert faces==len(pd)+2
    unseen=set(ends);groups=[]
    while unseen:
        start=min(unseen);group={start};queue=[start]
        for v in queue:
            for w in adj[v]-group:group.add(w);queue.append(w)
        unseen-=group;groups.append(sorted(group))
    groups.sort();where={e:i for i,g in enumerate(groups) for e in g}
    assert len(groups)==len(pd)
    crossings=[(where[b],where[a],where[c]) for a,b,c,d in pd]
    return crossings,twin
def named_binding(case):
    _,a=diagram(case["pd"]);_,b=diagram(case["census_pd"])
    w=case["census_map_witness"];f=w["dart_images"]
    orientation,swap=w["orientation"],w["over_under_swap"]
    assert orientation in [-1,1] and swap in [0,1]
    assert sorted(f)==list(range(4*len(case["pd"])))
    for x,y in enumerate(f):
        assert f[a[x]]==b[y]
        assert ((x%2)^(y%2))==swap
        assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+orientation)%4
def check_case(case,G):
    named_binding(case);crossings,_=diagram(case["pd"])
    seeds=case["seed_arcs"];assert len(seeds)==len(set(seeds))==4
    known=set(seeds)
    for o,a,b in case["seed_trace"]:
        assert o in known and a in known and b not in known
        assert (o,a,b) in crossings or (o,b,a) in crossings
        known.add(b)
    assert known==set(range(13))
    perms=[tuple(matvec(rows,v)-1 for v in range(1,8)) for rows in case["meridian_matrices_F2"]]
    assert [list(p) for p in perms]==case["permutations"]
    assert all(p in G and compose(p,p)==ID and p!=ID for p in perms)
    translations=case["arc_translations_F3"]
    assert len(translations)==13
    for p,v in zip(perms,translations):
        assert len(v)==7 and all(x in [0,1,2] for x in v)
        assert sum(v)%3==0 and all((a+b)%3==0 for a,b in zip(v,act(p,v)))
    matrices=[homogeneous(p,v) for p,v in zip(perms,translations)]
    assert matrices==case["homogeneous_W_matrices_F3"],"literal homogeneous matrix changed"
    identity=[[int(i==j) for j in range(7)] for i in range(7)]
    assert all(multiply(m,m)==identity for m in matrices)
    assert all(rank(difference(m),7)==2 for m in matrices)
    for o,a,b in crossings:
        assert multiply(multiply(matrices[o],matrices[a]),matrices[o])==matrices[b],"Wirtinger matrix relation fails"
    assert rank([row for a in seeds for row in difference(matrices[a])],7)==7,"nonzero common fixed vector"
    seedgens=[(perms[a],translations[a]) for a in seeds]
    assert projection_group([p for p,v in seedgens])==G
    vectors=case["kernel_translation_basis_F3"];words=case["kernel_basis_words_on_seed_indices"]
    assert len(vectors)==len(words)==6
    assert all(sum(v)%3==0 for v in vectors)
    for v,w in zip(vectors,words):assert evaluate(w,seedgens)==(ID,v),"kernel word fails"
    assert rank(vectors,7)==6
    return {"name":case["name"],"image_order":122472,"kernel_translation_rank":6,
            "meridian_displacement_rank":2,"common_fixed_dimension":0,
            "meridional_lower":4,"bridge_upper":4,"mr_equals_bridge":4}
def full_group_control(case,ordinary_words):
    seeds=case["seed_arcs"];gens=[(tuple(case["permutations"][a]),case["arc_translations_F3"][a]) for a in seeds]
    proj=list(sorted(projection_group([p for p,v in gens])));where={p:i for i,p in enumerate(proj)}
    left=[[where[compose(p,q)] for q in proj] for p,v in gens]
    vectors=[[code//3**i%3 for i in range(6)] for code in range(729)]
    vectors=[v+[-sum(v)%3] for v in vectors]
    def encode(v):return sum(v[i]*3**i for i in range(6))
    maps=[[encode([(a+b)%3 for a,b in zip(v,act(p,w))]) for w in vectors] for p,v in gens]
    start=where[ID]*729;seen={start};queue=[start]
    for x in queue:
        g,c=divmod(x,729)
        for l,m in zip(left,maps):
            y=l[g]*729+m[c]
            if y not in seen:seen.add(y);queue.append(y)
    assert len(seen)==122472
    assert ordinary_words==[[0,2],[1,3]]
    ordinary=[evaluate(w,gens) for w in ordinary_words]
    left2=[[where[compose(p,q)] for q in proj] for p,v in ordinary]
    maps2=[[encode([(a+b)%3 for a,b in zip(v,act(p,w))]) for w in vectors] for p,v in ordinary]
    seen2={start};queue2=[start]
    for x in queue2:
        g,c=divmod(x,729)
        for l,m in zip(left2,maps2):
            y=l[g]*729+m[c]
            if y not in seen2:seen2.add(y);queue2.append(y)
    assert len(seen2)==122472
    assert compose(gens[0][0],gens[2][0])!=compose(gens[2][0],gens[0][0])
    first=(gens[0][0],tuple(gens[0][1]));orbit={first};queue=[first]
    for p,v in queue:
        for gen in gens:
            q,w=affine_mul(affine_mul(gen,(p,list(v))),gen);x=(q,tuple(w))
            if x not in orbit:orbit.add(x);queue.append(x)
    assert len(orbit)==189
    assert all((tuple(p),tuple(v)) in orbit for p,v in zip(case["permutations"],case["arc_translations_F3"]))
    return {"full_affine_closure_order":122472,"meridian_conjugacy_class_size":189,
            "ordinary_two_word_closure_order":122472,"ordinary_generating_rank":2}
def rejected(call):
    try:call()
    except (AssertionError,IndexError,KeyError):return True
    raise AssertionError("corrupted input accepted")
def main():
    data=json.loads(INPUT.read_text())
    assert [c["name"] for c in data["cases"]]==NAMES
    for path,pin in data["source_sha256"].items():assert digest(HERE/path)==pin
    old=json.loads((HERE/"seven_sheet_cover_inputs.json").read_text());old={c["name"]:c for c in old["cases"]}
    census=json.loads((HERE/"seven_sheet_census_verification.json").read_text())
    assert census["status"]=="PASS" and census["external_census_reread_count"]==29
    for case in data["cases"]:
        for k in ["pd","census_pd","census_map_witness","seed_arcs","seed_trace"]:assert case[k]==old[case["name"]][k]
    G=gl3();involutions={p for p in G if p!=ID and compose(p,p)==ID}
    assert len(involutions)==21
    mu=next(iter(involutions))
    inverse={p:tuple(p.index(i) for i in range(7)) for p in G}
    assert {compose(compose(p,mu),inverse[p]) for p in G}==involutions
    results=[check_case(c,G) for c in data["cases"]]
    finite=full_group_control(data["cases"][0],data["ordinary_generator_words_on_first_case_seeds"])
    bad=deepcopy(data["cases"][0]);bad["homogeneous_W_matrices_F3"][0][0][0]^=1
    bad2=deepcopy(data["cases"][0]);bad2["arc_translations_F3"]=[[0]*7 for _ in range(13)]
    bad2["homogeneous_W_matrices_F3"]=[homogeneous(tuple(p),v) for p,v in zip(bad2["permutations"],bad2["arc_translations_F3"])]
    bad3=deepcopy(data["cases"][0]);bad3["kernel_basis_words_on_seed_indices"][0]=[]
    bad4=deepcopy(data["cases"][0]);bad4["census_map_witness"]["dart_images"][0]=bad4["census_map_witness"]["dart_images"][1]
    controls={"changed_matrix":rejected(lambda:check_case(bad,G)),
              "principal_zero_translation":rejected(lambda:check_case(bad2,G)),
              "changed_kernel_word":rejected(lambda:check_case(bad3,G)),
              "changed_named_census_map":rejected(lambda:check_case(bad4,G))}
    files=[INPUT,Path(__file__)]+[HERE/p for p in data["source_sha256"]]
    receipt={"status":"PASS","python":platform.python_version(),"scope":data["scope"],
             "results":results,"group_control":finite,"projected_involution_class_size":21,
             "class_generating_rank":4,"negative_controls":controls,
             "input_sha256":{p.relative_to(HERE).as_posix():digest(p) for p in files},
             "assurance":"Direct exact finite matrix, word, group and diagram replay. No branched-cover group or rank bound is used. Named-census read is bound to its separate retained receipt; the published geometric bridge-Wirtinger theorem is imported. Internal replay only, without formal verification or external specialist review."}
    (HERE/"ternary_affine_quotient_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","full_Q_epimorphisms":len(results),"group_control":finite,
                      "class_generating_rank":4,"mr_equals_bridge":4,"negative_controls":controls},indent=2))
if __name__=="__main__":main()

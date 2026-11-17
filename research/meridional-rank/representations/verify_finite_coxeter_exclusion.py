
#!/usr/bin/env python3
"""Exact finite-Coxeter exclusion on 24 retained knot diagrams.

All H4 reflection colorings are enumerated, allowing repeated seeds.
Coordinates use Z[phi], phi^2=phi+1. No historical generating-set lists
or floating point group computations are imported.
"""
import json, itertools, hashlib
from pathlib import Path
from verify_seven_sheet_covers import diagram, census_bijection
D = Path(__file__).resolve().parent
zero=(0,0)
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def mul(a,b):return(a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def dot(a,b):
 out=zero
 for v,w in zip(a,b):out=add(out,mul(v,w))
 return out
def canonical(v):return min(v,tuple(neg(x) for x in v))
roots=set()
for i in range(4):
 for sg in [-1,1]:
  v=[zero]*4;v[i]=(2*sg,0);roots.add(tuple(v))
for signs in itertools.product([-1,1],repeat=4):roots.add(tuple((s,0) for s in signs))
for perm in itertools.permutations(range(4)):
 if sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))%2:continue
 for signs in itertools.product([-1,1],repeat=3):
  a=[zero,(signs[0],0),(0,signs[1]),(-signs[2],signs[2])]
  roots.add(tuple(a[perm[i]] for i in range(4)))
assert len(roots)==120 and all(dot(r,r)==(4,0) for r in roots)
C=sorted({canonical(r) for r in roots});assert len(C)==60
ix={r:i for i,r in enumerate(C)}
Q=[]
for u in C:
 row=[]
 for v in C:
  d=dot(u,v);assert d[0]%2==d[1]%2==0
  f=(d[0]//2,d[1]//2)
  w=tuple(add(a,neg(mul(f,b))) for a,b in zip(v,u))
  row.append(ix[canonical(w)])
 Q.append(row)
assert all(Q[i][Q[i][j]]==j for i in range(60) for j in range(60))
seen={0};todo=[0]
for x in todo:
 for row in Q:
  y=row[x]
  if y not in seen:seen.add(y);todo.append(y)
assert len(seen)==60
# Exact span determinant over Z[phi], only needed for accepted colorings.
perms4=list(itertools.permutations(range(4)))
sgns=[(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)) for p in perms4]
def rank4(seeds):
 out=zero
 for perm,sgn in zip(perms4,sgns):
  value=(sgn,0)
  for i,j in enumerate(perm):value=mul(value,C[seeds[i]][j])
  out=add(out,value)
 return out!=zero

def reflect(u,v):
    d=dot(u,v)
    assert d[0]%2==d[1]%2==0
    factor=(d[0]//2,d[1]//2)
    return tuple(add(a,neg(mul(factor,b))) for a,b in zip(v,u))


def simple_control():
    # Twice unit normals, with Cartan labels 5-3-3.
    simple=[((-2,0),zero,zero,zero),
            ((0,1),(-1,0),(-1,1),zero),
            (zero,(0,1),(1,-1),(-1,0)),
            (zero,zero,zero,(2,0))]
    assert all(v in roots for v in simple)
    assert rank4(tuple(ix[canonical(v)] for v in simple))
    target={(0,1):(0,-2),(1,2):(-2,0),(2,3):(-2,0)}
    allroots=sorted(roots); idx={v:i for i,v in enumerate(allroots)}
    perms=[tuple(idx[reflect(v,w)] for w in allroots) for v in simple]
    identity=tuple(range(120))
    for i,j in itertools.combinations(range(4),2):
        assert dot(simple[i],simple[j])==target.get((i,j),zero)
        product=tuple(perms[i][perms[j][k]] for k in range(120))
        value=identity;order=0
        while True:
            order+=1;value=tuple(product[k] for k in value)
            if value==identity:break
            assert order<=5
        assert order==({(0,1):5,(1,2):3,(2,3):3}.get((i,j),2))
    # These four reflections alone conjugate through every root direction.
    seen={0};queue=[0]
    for v in queue:
        for u in simple:
            w=Q[ix[canonical(u)]][v]
            if w not in seen:seen.add(w);queue.append(w)
    assert len(seen)==60
    # A known sharp pretzel cycle quotient P(5,3,3,2) accepts these roots.
    return {"simple_twice_roots":simple,"reflection_orders":[5,3,3,2],
            "root_directions":60,"signed_roots":120,"rank4_positive_control":True}


EXPECTED_H4 = {
"13a_373":121,"13a_375":1,"13a_391":1,"13a_392":121,
"13a_3264":121,"13a_3266":1,"13a_3275":1,"13a_3276":121,
"13n_1005":121,"13n_1006":121,"13n_1007":121,"13n_1011":1,
"13n_1012":1,"13n_1013":25,"13n_1043":1,"13n_1044":1,
"13n_1045":25,"13n_1046":121,"13n_1047":121,"13n_1048":121,
"13n_3842":121,"13n_3844":1,"13n_3855":1,"13n_3856":121}

EXPECTED_H4.update({'13a_304': 121, '13a_362': 121, '13a_623': 241, '13a_629': 241, '13a_3232': 1, '13a_3259': 1, '13a_3374': 121, '13a_3375': 121, '13n_834': 25, '13n_838': 1, '13n_839': 1, '13n_991': 25, '13n_992': 1, '13n_993': 1, '13n_1367': 121, '13n_1370': 121, '13n_1371': 121, '13n_1374': 121, '13n_1375': 121, '13n_1376': 121, '13n_3814': 1, '13n_3840': 1, '13n_3911': 121, '13n_3914': 121})

def verify_plan(case):
    cross,_=diagram(case["pd"])
    census_bijection(case["pd"],case["census_pd"],case["census_map_witness"])
    seen=set(case["seed_arcs"])
    assert len(seen)==len(case["seed_arcs"])==4
    for o,a,b in case["seed_trace"]:
        assert o in seen and a in seen and b not in seen
        assert (o,a,b) in cross or (o,b,a) in cross
        seen.add(b)
    assert seen==set(range(len(case["pd"])))
    return cross


def search(case,conjugations,full_rank=None):
    cross=verify_plan(case);count=len(conjugations)
    accepted=sharp=0
    for a,b,c in itertools.product(range(count),repeat=3):
        labels=[None]*len(case["pd"])
        for i,j in zip(case["seed_arcs"],(0,a,b,c)):labels[i]=j
        for o,x,y in case["seed_trace"]:
            labels[y]=conjugations[labels[o]][labels[x]]
        if not all(conjugations[labels[o]][labels[x]]==labels[y] for o,x,y in cross):
            continue
        accepted+=1
        if full_rank is not None:
            sharp+=full_rank((0,a,b,c))
    return {"normalized_trials":count**3,"accepted":accepted,"rank4_seed_colorings":sharp}


def symmetric_table(degree):
    pairs=list(itertools.combinations(range(degree),2))
    table=[]
    for p,q in pairs:
        def f(v):return q if v==p else p if v==q else v
        table.append([pairs.index(tuple(sorted((f(a),f(b))))) for a,b in pairs])
    return table


def main():
    import copy
    path=D/"seven_sheet_cover_inputs.json";inputs=json.loads(path.read_text())
    cases=[c for c in inputs["cases"] if c["actions"][0]["expected_H1_F3_dimension"]]
    cases += json.loads((D/"a6_affine_inputs.json").read_text())["cases"]
    assert [c["name"] for c in cases]==list(EXPECTED_H4)
    controls=simple_control()
    results=[]
    sym={r:symmetric_table(r) for r in (3,4,5)}
    for case in cases:
        row={"name":case["name"],"symmetric":{}}
        for degree,table in sym.items():
            result=search(case,table)
            assert result["accepted"]==1
            row["symmetric"][str(degree)]=result
        h4=search(case,Q,rank4)
        assert h4["accepted"]==EXPECTED_H4[case["name"]]
        assert h4["rank4_seed_colorings"]==0
        row["H4"]=h4
        results.append(row)
    bad=copy.deepcopy(cases[0])
    bad["seed_trace"][0][2]=bad["seed_arcs"][0]
    try:verify_plan(bad)
    except AssertionError:pass
    else:raise AssertionError("bad propagation plan accepted")
    # Repeating a simple root destroys the full-rank positive control.
    assert not rank4((0,0,0,0))
    files=[Path(__file__),path,D/"verify_seven_sheet_covers.py",D/"FINITE_COXETER_EXCLUSION.tex",D/"a6_affine_inputs.json"]
    receipt={"status":"PASS","scope":"All 48 retained affine sharp-four knots: no rank-four finite Coxeter reflection quotient. Complete first-seed-normalized transposition and H4 root-coloring searches; classification/projection reductions are written dependencies. Infinite Coxeter quotients are separately excluded on sixteen alternating cases; other cases remain outside this finite search.",
      "results":results,"controls":controls,"corruptions_rejected":["invalid seed propagation","dependent roots fail full-rank control"],
      "input_sha256":{p.relative_to(D).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
      "assurance":"Exact internal computation, without formal verification or external human review."}
    (D/"finite_coxeter_exclusion_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","knots":48,"H4_trials_per_knot":60**3,
                     "sharp_finite_Coxeter_hits":0,"symmetric_degrees":[3,4,5]},indent=2))


if __name__=="__main__":main()

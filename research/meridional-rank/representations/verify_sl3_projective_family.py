#!/usr/bin/env python3
"""Exact50-case projective-cover equality and known-family audit.

Requires pinned SnapPy/Spherogram/Regina and named census versions.
Replays full directed13-generator13-sheet chains with every meridian
cycle fill, plus actual double-cover decompositions into Seifert pieces.
"""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
from copy import deepcopy
from importlib.metadata import version
import json,hashlib,random,warnings
warnings.filterwarnings("ignore",message="Plink failed to import tkinter.*")
import snappy,spherogram,regina,snappy_15_knots
from spherogram import Link

HERE=Path(__file__).resolve().parent
INPUT=HERE/"sl3_projective_family_inputs.json"
EXPECTED=["13a_3167","13a_3168","13a_3169","13a_3170","13a_3188","13a_3189","13a_3190","13a_3191",
"13a_3293","13a_3294","13a_3295","13a_3298","13a_3347","13a_3348","13a_3351","13a_3402","13a_3403","13a_3411",
"13n_326","13n_329","13n_330","13n_341","13n_344","13n_345","13n_348","13n_351","13n_352","13n_355","13n_357","13n_358",
"13n_510","13n_513","13n_514","13n_518","13n_520","13n_521","13n_524","13n_525","13n_526","13n_529","13n_530","13n_531",
"13n_3739","13n_3742","13n_3745","13n_3748","13n_3770","13n_3773","13n_3776","13n_3779"]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(pd):
    endpoints=defaultdict(list);adj=defaultdict(set)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):endpoints[e].append(4*c+h)
        a,b=row[1],row[3];adj[a].add(b);adj[b].add(a)
    assert len(endpoints)==2*len(pd) and all(len(v)==2 for v in endpoints.values())
    twin={a:b for pair in endpoints.values() for a,b in [pair,pair[::-1]]}
    seen=set();v=0
    for _ in range(2*len(pd)):
        assert v not in seen;seen.add(v);v=twin[4*(v//4)+(v+2)%4]
    assert v==0 and len({v//4 for v in seen})==len(pd)
    unseen=set(twin);faces=0
    while unseen:
        start=v=min(unseen)
        while True:
            unseen.remove(v);x=twin[v];v=4*(x//4)+(x+1)%4
            if v==start:break
        faces+=1
    assert faces==len(pd)+2
    unseen=set(endpoints);components=[]
    while unseen:
        start=min(unseen);comp={start};queue=[start]
        for e in queue:
            for f in adj[e]-comp:comp.add(f);queue.append(f)
        unseen-=comp;components.append(sorted(comp))
    components.sort();where={e:i for i,c in enumerate(components) for e in c}
    cr=[(where[b],where[a],where[c]) for a,b,c,d in pd]
    assert len(components)==len(pd)
    return cr,twin
def full_chain(perms,cr):
    n=len(perms[0]);m=len(perms)
    edges=[(g,v) for g in range(m) for v in range(n)]
    parent=list(range(n))
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    tree=set()
    # Ascending full-generator tree differs from the compressed seed chain.
    for j,(g,a) in enumerate(edges):
        b=perms[g][a]
        if root(a)!=root(b):parent[root(a)]=root(b);tree.add(j)
    assert len(tree)==n-1
    chords=[j for j in range(len(edges)) if j not in tree]
    rows=[]
    for o,a,b in cr:
        for start in range(n):
            v=start;row=[0]*len(edges)
            for g in [o,a,o,b]:row[g*n+v]+=1;v=perms[g][v]
            assert v==start
            rows.append([row[j] for j in chords])
    cycles=[]
    for g,p in enumerate(perms):
        unused=set(range(n));cy=[]
        while unused:
            start=min(unused);v=start;row=[0]*len(edges);cycle=[]
            while True:
                unused.remove(v);cycle.append(v);row[g*n+v]+=1;v=p[v]
                if v==start:break
            rows.append([row[j] for j in chords]);cy.append(cycle)
        cycles.append(cy)
    assert all(sorted(map(len,c))==[1,1,1,1,1,2,2,2,2] for c in cycles)
    return rows,len(chords),tree,cycles
def rank_mod3(rows,columns):
    a=[[v%3 for v in row] for row in rows];r=0
    for j in range(columns):
        i=next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None:continue
        a[r],a[i]=a[i],a[r]
        factor=1 if a[r][j]==1 else 2;a[r]=[factor*v%3 for v in a[r]]
        for i in range(r+1,len(a)):
            if a[i][j]:
                factor=a[i][j];a[i]=[(v-factor*w)%3 for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r
def composition(p,q):return tuple(p[q[i]] for i in range(len(p)))
def group_order(perms):
    known={tuple(range(len(perms[0])))};queue=list(known)
    for p in queue:
        for q in perms:
            v=composition(q,p)
            if v not in known:known.add(v);queue.append(v)
    return len(known)
def map_witness(source,target):
    _,a=parse(source);_,b=parse(target);n=len(source)
    def rotate(x,e):return 4*(x//4)+(x+e)%4
    for orientation in [1,-1]:
        for first in range(4*n):
            mapping={0:first};queue=[0];used={first};valid=True
            for x in queue:
                y=mapping[x]
                for u,v in [(a[x],b[y]),(rotate(x,1),rotate(y,orientation))]:
                    if u in mapping:
                        if mapping[u]!=v:valid=False;break
                    else:
                        if v in used:valid=False;break
                        mapping[u]=v;used.add(v);queue.append(u)
                if not valid:break
            if valid and len(mapping)==4*n:
                f=[mapping[i] for i in range(4*n)];swap=f[0]%2
                if all((i%2)^(v%2)==swap for i,v in enumerate(f)):
                    assert all(f[a[x]]==b[f[x]] for x in a)
                    return {"dart_images":f,"orientation":orientation,"over_under_swap":swap}
    raise AssertionError("No exact crossing-preserving census map")

def mul3(a,b):
    return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3))%3 for i in range(3) for j in range(3))
def det3(a):
    return (a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]))%3
def projective_action(a,points):
    lookup={tuple(p):i for i,p in enumerate(points)};result=[]
    for p in points:
        q=[sum(a[3*i+j]*p[j] for j in range(3))%3 for i in range(3)]
        scale=pow(next(x for x in q if x),-1,3)
        result.append(lookup[tuple(scale*x%3 for x in q)])
    assert sorted(result)==list(range(13));return tuple(result)
def representation(row,points):
    perms=[tuple(p) for p in row["projective_permutations"]]
    assert len(perms)==len(row["pd"])==13
    identity=tuple(range(13));eye=(1,0,0,0,1,0,0,0,1)
    for a,p in zip(row["involution_matrices_F3"],perms):
        assert len(a)==9 and all(x in [0,1,2] for x in a)
        assert det3(a)==1 and mul3(a,a)==eye and tuple(a)!=eye
        assert projective_action(a,points)==p and composition(p,p)==identity
        unused=set(range(13));cycles=[]
        while unused:
            v=start=min(unused);length=0
            while True:
                unused.remove(v);length+=1;v=p[v]
                if v==start:break
            cycles.append(length)
        assert sorted(cycles)==[1,1,1,1,1,2,2,2,2]
    cr,_=parse(row["pd"])
    for o,a,b in cr:assert composition(composition(perms[o],perms[a]),perms[o])==perms[b]
    return perms,cr
def four_seed_upper(row,cr):
    known=set(row["seed_arcs"]);assert len(known)==4
    used=set()
    for o,a,b in row["seed_trace"]:
        assert o in known and a in known and b not in known
        assert (o,a,b) in cr or (o,b,a) in cr
        key=(o,)+tuple(sorted([a,b]));assert key not in used
        used.add(key);known.add(b)
    assert known==set(range(13)) and len(used)==9
def named_binding(row):
    L=Link("K"+row["name"].replace("_",""))
    pd=[list(r) for r in L.PD_code()]
    assert pd==row["census_pd"] and json.loads(json.dumps(L.DT_code()))==row["census_DT"]
    assert map_witness(row["pd"],pd)==row["census_map_witness"]
    assert row["census_map_witness"]["orientation"]==1 and row["census_map_witness"]["over_under_swap"]==0
def close_cover(row):
    t=row["double_cover"];seed=t["random_seed"]
    random.seed(seed);snappy.set_rand_seed(seed);regina.RandomEngine.reseedWithDefault()
    E=Link(row["pd"]).exterior()
    assert E.num_cusps()==1 and E.triangulation_isosig()==t["exterior_isosig_with_peripheral"]
    E.dehn_fill((2,0));covers=E.covers(2);assert len(covers)==1
    C=covers[0];assert C.cover_info()["degree"]==2 and C.cover_info()["type"]=="cyclic"
    assert C.num_cusps()==1 and C.cusp_info("filling")==[(1.0,0.0)]
    T=regina.Triangulation3(C.filled_triangulation()._to_string());T.simplify()
    assert T.isoSig()==t["closed_cover_isosig"] and T.size()==t["closed_tetrahedra"]
    assert T.isValid() and T.isClosed() and T.isConnected() and T.isOrientable()
    assert str(T.homology())==t["closed_H1"] and T.isIrreducible()
    return T
def full_vector(surface,T):
    return tuple(int(surface.triangles(j,k).longValue()) for j in range(T.size()) for k in range(4))+tuple(int(surface.quads(j,k).longValue()) for j in range(T.size()) for k in range(3))
def quad_vector(surface,T):
    return [int(surface.quads(j,k).longValue()) for j in range(T.size()) for k in range(3)]
def recognised_piece(U):
    assert U.isValid() and U.isConnected() and U.isOrientable()
    assert all(b.eulerChar()==0 for b in U.boundaryComponents())
    standard=regina.StandardTriangulation.recognise(U);assert standard is not None
    manifold=standard.manifold();assert isinstance(manifold,regina.SFSpace)
    return {"finite_isosig":U.isoSig(),"tetrahedra":U.size(),"boundary_components":U.countBoundaryComponents(),
            "recognition":{"triangulation_type":type(standard).__name__,"standard_name":standard.name(),
                           "manifold_type":type(manifold).__name__,"manifold":manifold.name()}}
def graph_cover(row):
    t=row["double_cover"];T=close_cover(row)
    all_surfaces=regina.NormalSurfaces(T,regina.NormalCoords.Quad,regina.NormalList.Vertex)
    surfaces=[all_surfaces[i] for i in t["disjoint_normal_torus_indices"]]
    assert len(surfaces) in [3,4]
    for surface,quad in zip(surfaces,t["disjoint_normal_quad_vectors"]):
        assert quad_vector(surface,T)==quad
        assert surface.isConnected() and surface.isOrientable() and surface.isTwoSided() and surface.eulerChar()==0
        assert surface.isIncompressible()
    assert len(t["disjoint_normal_quad_vectors"])==len(surfaces)
    assert all(a.disjoint(b) for a,b in combinations(surfaces,2))
    union=surfaces[0]
    for surface in surfaces[1:]:union=union+surface
    assert sorted(full_vector(s,T) for s in union.components())==sorted(full_vector(s,T) for s in surfaces)
    regina.RandomEngine.reseedWithDefault();cut=union.cutAlong();cut.simplify()
    pieces=[]
    for U in cut.triangulateComponents():U.simplify();pieces.append(recognised_piece(U))
    assert len(pieces)==len(surfaces)+1
    assert sorted(pieces,key=lambda p:p["finite_isosig"])==sorted(t["pieces"],key=lambda p:p["finite_isosig"])
    return {"closed_cover_isosig":T.isoSig(),"irreducible":True,"disjoint_incompressible_tori":len(surfaces),
            "Seifert_pieces":len(pieces),"prime_graph_manifold_double_cover":True,
            "within_BJW18_graph_cover_family":True}
def rejected(call):
    try:call()
    except (AssertionError,IndexError,KeyError,ValueError,RuntimeError,StopIteration):return True
    raise AssertionError("corrupted assertion accepted")
def main():
    data=json.loads(INPUT.read_text());assert [r["name"] for r in data["cases"]]==EXPECTED
    assert data["field"]==3 and data["degree"]==13 and data["projected_group_order"]==5616
    assert data["expected_H1_F3_dimension"]==2 and data["meridian_cycle_type"]==[1,1,1,1,1,2,2,2,2]
    versions={p:version(p) for p in ["snappy","spherogram","snappy_15_knots","snappy_manifolds"]};versions["regina"]=regina.versionString()
    assert versions=={"snappy":"3.3.2","spherogram":"2.4.1","snappy_15_knots":"1.2.1","snappy_manifolds":"1.4","regina":"7.4"}
    for name,pin in data["source_sha256"].items():assert sha(HERE/name)==pin
    provenance=data["census_provenance"];package=Path(snappy_15_knots.__file__).parent
    assert provenance["version"]==versions["snappy_15_knots"] and provenance["spherogram_version"]==versions["spherogram"]
    for pin in provenance["source_files"]:
        raw=(package/pin["path"]).read_bytes();assert len(raw)==pin["bytes"] and hashlib.sha256(raw).hexdigest()==pin["sha256"]
    panel={r["name"]:r["pd"] for r in json.loads((HERE/"joint_a5_panel.json").read_text())["cases"]}
    rows=[]
    for row in data["cases"]:
        assert row["pd"]==panel[row["name"]];named_binding(row)
        perms,cr=representation(row,data["projective_points"]);four_seed_upper(row,cr)
        assert group_order([perms[a] for a in row["seed_arcs"]])==5616
        matrix,ncols,tree,cycles=full_chain(perms,cr);rank=rank_mod3(matrix,ncols)
        assert ncols==157 and len(matrix)==286 and rank==155
        h=ncols-rank;assert h==2
        lower=(12+h+3)//4;assert lower==4
        topology=graph_cover(row)
        rows.append({"name":row["name"],"projected_group_order":5616,"degree":13,"meridian_delta":4,
                     "directed_generator_edges":169,"tree_edges":sorted(tree),"columns":157,
                     "lifted_crossing_rows":169,"meridian_cycle_fill_rows":117,"fixed_sheet_fill_rows":65,
                     "rank_F3":155,"H1_F3_dimension":2,"meridional_lower":4,"four_seed_bridge_upper":4,
                     "meridional_rank":4,"bridge_number":4,**topology})
        print("PASS",row["name"],"H1F3=2 mr=b=4 GRAPH",topology["Seifert_pieces"],flush=True)
    first=data["cases"][0]
    badname=deepcopy(first);badname["name"]="13a_373"
    badperm=deepcopy(first);badperm["projective_permutations"][0][0],badperm["projective_permutations"][0][1]=badperm["projective_permutations"][0][1],badperm["projective_permutations"][0][0]
    T=close_cover(first);surfaces=regina.NormalSurfaces(T,regina.NormalCoords.Quad,regina.NormalList.Vertex)
    chosen=surfaces[first["double_cover"]["disjoint_normal_torus_indices"][0]]
    badquad=list(first["double_cover"]["disjoint_normal_quad_vectors"][0]);badquad[0]+=1
    def same_quad():assert quad_vector(chosen,T)==badquad
    controls={"wrong_named_knot":rejected(lambda:named_binding(badname)),
              "changed_projective_permutation":rejected(lambda:representation(badperm,data["projective_points"])),
              "changed_normal_quad_vector":rejected(same_quad)}
    files=[INPUT,Path(__file__),HERE/"SL3_PROJECTIVE_FAMILY_AUDIT.tex"]+[HERE/p for p in data["source_sha256"]]
    receipt={"status":"PASS","scope":data["scope"],"versions":versions,"cases":rows,"negative_controls":controls,
             "input_sha256":{p.relative_to(HERE).as_posix():sha(p) for p in files},
             "assurance":"Exact internal named census, full cover-chain and actual graph-cover topology replay. Trusted SnapPy/Regina construction/recognition algorithms and written general-cover/BKVV/BJW imports. No new outside-family cases, worldwide priority, formal verification or external specialist review asserted."}
    (HERE/"sl3_projective_family_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS50 named equalities and50 existing BJW-family memberships",flush=True)
if __name__=="__main__":main()

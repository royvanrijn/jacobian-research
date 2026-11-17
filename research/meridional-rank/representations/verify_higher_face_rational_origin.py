#!/usr/bin/env python3
"""Exact finite checker for the bridge-eight rational origin of13a373.

Standard library only. BBK Theorem1 and the local band-to-twist
interpretation remain written mathematical imports. No MRC
counterexample or rank-preserving replacement is claimed.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import platform

HERE=Path(__file__).resolve().parent
INPUT=HERE/"higher_face_rational_origin_inputs.json"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def corners(pd):
    ends=defaultdict(list)
    for i,row in enumerate(pd):
        assert len(row)==4
        for h,e in enumerate(row):ends[e].append((i,h))
    assert all(len(v)==2 for v in ends.values())
    twin={a:b for pair in ends.values() for a,b in [pair,pair[::-1]]}
    unseen=set(twin);faces=[]
    while unseen:
        start=min(unseen);cycle=[];a=start
        while a not in cycle:
            assert a in unseen
            cycle.append(a);unseen.remove(a)
            i,h=twin[a];a=(i,(h-1)%4)
        assert a==start
        faces.append(cycle)
    assert len(faces)==len(pd)+2, "rotation system is not spherical"
    where={a:i for i,f in enumerate(faces) for a in f}
    adj={i:set() for i in range(len(faces))}
    for i in range(len(pd)):
        for h in range(4):
            a,b=where[i,h],where[i,(h+1)%4]
            adj[a].add(b);adj[b].add(a)
    colour={where[0,3]:0};queue=list(colour)
    for a in queue:
        for b in adj[a]:
            if b not in colour:colour[b]=1-colour[a];queue.append(b)
            assert colour[b]!=colour[a]
    assert len(colour)==len(faces), "disconnected projection"
    vertices=[i for i in adj if colour[i]==0]
    parity=[colour[where[i,0]] for i in range(len(pd))]
    # A G-coloured corner has index0 for parity0, index1 for parity1.
    edges=[(where[i,p],where[i,(p+2)%4]) for i,p in enumerate(parity)]
    return vertices,edges,parity

def strand_components(pd,pairings=None):
    parent={x:x for row in pd for x in row}
    def find(a):
        while parent[a]!=a:a=parent[a]
        return a
    def union(a,b):parent[find(a)]=find(b)
    for i,row in enumerate(pd):
        pairs=((0,2),(1,3)) if pairings is None or i not in pairings else pairings[i]
        for a,b in pairs:union(row[a],row[b])
    return len({find(a) for a in parent})

def expand_twist(row,n,p,nextlabel):
    assert n>=2, "base twist weight is below2"
    def recur(row,n,nextlabel):
        if n==1:return [row],nextlabel
        a,b,c,d=row;x,y=nextlabel,nextlabel+1
        if p==0:first,second=[a,x,y,d],[x,b,c,y]
        else:first,second=[a,b,x,y],[y,x,c,d]
        out,nextlabel=recur(first,n-1,nextlabel+2)
        return out+[second],nextlabel
    return recur(row,n,nextlabel)

def check(data):
    pd=data["target_pd"]
    assert len(pd)==13
    assert pd[3:5]==[[13,4,14,5],[5,14,6,15]]
    shadow=pd[:3]+[[13,4,6,15]]+pd[5:]
    assert shadow==data["box_shadow_pd"]
    vertices,edges,parity=corners(shadow)
    assert len(vertices)==8 and len(edges)==12
    assert all(a!=b for a,b in edges)
    assert len({tuple(sorted(e)) for e in edges})==12
    assert parity==data["G_corner_parity"]
    weights=data["absolute_twist_weights"]
    assert weights==[3,2,3,3,3,3,3,3,3,3,3,3]
    base=[];groups=[];nxt=27
    for row,n,p in zip(shadow,weights,parity):
        rows,nxt=expand_twist(row,n,p,nxt)
        groups.append(list(range(len(base),len(base)+len(rows))))
        base.extend(rows)
    assert base==data["base_pd"], "literal base diagram differs"
    assert len(base)==35 and strand_components(base)==1
    vv,ee,_=corners(base)
    assert len(vv)==8 and len(ee)==35
    collapsed=[]
    for group,n in zip(groups,weights):
        pairs={tuple(sorted(ee[i])) for i in group}
        assert len(group)==n and len(pairs)==1, "twist is not a parallel G bundle"
        a,b=next(iter(pairs));assert a!=b
        collapsed.append((a,b))
    assert len(set(collapsed))==12
    assert set(v for edge in collapsed for v in edge)==set(vv)
    ends=defaultdict(list)
    for row in base:
        for h,e in enumerate(row):ends[e].append(h)
    assert all(sorted(h%2 for h in pair)==[0,1] for pair in ends.values()), "base is not alternating"
    def graph_connected_without(removed):
        adjacent={v:set() for v in vv}
        for i,(a,b) in enumerate(ee):
            if i!=removed:adjacent[a].add(b);adjacent[b].add(a)
        reached={vv[0]};queue=list(reached)
        for a in queue:
            for b in adjacent[a]-reached:reached.add(b);queue.append(b)
        return reached==set(vv)
    assert all(a!=b for a,b in ee)
    assert all(graph_connected_without(i) for i in range(len(ee))), "G crossing graph has a bridge"
    # Reinsert exactly the original rational crossing boxes. The
    # exceptional box is its displayed two-crossing twist; all others
    # are single crossings, hence all are rational.
    recovered=shadow[:3]+pd[3:5]+shadow[4:]
    assert recovered==data["recovered_pd"]==pd, "target recovery differs"
    assert strand_components(pd)==1
    return {"base_crossings":35,"base_components":1,
            "simple_dual_vertices":8,"simple_dual_edges":12,
            "base_bridge_number_by_BBK_Theorem1":8,
            "all_base_twist_weights_at_least2":True,
            "literal_target_recovery":True,
            "base_twist_crossing_groups":groups,
            "parallel_G_twist_bundles":True,"alternating":True,
            "full_G_loopless":True,"full_G_bridge_free":True,
            "collapsed_base_dual_edges":[list(e) for e in collapsed]}

def rejected(call):
    try:call()
    except AssertionError:return True
    raise AssertionError("invalid input accepted")

def main():
    data=json.loads(INPUT.read_text())
    for path,pin in data["source_sha256"].items():assert digest(HERE/path)==pin
    retained=json.loads((HERE/"seed_cover_topology_inputs.json").read_text())
    target=next(c for c in retained["cases"] if c["name"]=="13a_373")
    assert target["pd"]==data["target_pd"]
    old=json.loads((HERE/"seed_cover_topology_verification.json").read_text())
    assert old["status"]=="PASS" and any(c["name"]=="13a_373" for c in old["cases"])
    result=check(data)
    bad=deepcopy(data);bad["base_pd"][0][0]+=1000
    bad2=deepcopy(data);bad2["recovered_pd"][0][0]+=1000
    negative={"changed_base_PD":rejected(lambda:check(bad)),
              "changed_recovered_PD":rejected(lambda:check(bad2)),
              "unit_base_twist":rejected(lambda:expand_twist([1,2,3,4],1,0,5))}
    paths=[INPUT,Path(__file__)]+[HERE/p for p in data["source_sha256"]]
    receipt={"status":"PASS","python":platform.python_version(),"scope":data["scope"],
             "finite_example":result,"negative_controls":negative,
             "input_sha256":{str(p.relative_to(HERE)):digest(p) for p in paths},
             "assurance":"Internal exact diagram, planar rotation, strand connectivity, parallel dual twist bundles, alternating/reduced base and literal rational replacement checks. BBK twisted-base interpretation and bridge theorem are written imports. The general universality proposition is a written local-gluing and circuit-merging argument, not inferred from this finite check. No formal verification or external specialist review."}
    (HERE/"higher_face_rational_origin_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS",**{k:v for k,v in result.items() if not isinstance(v,list)},"negative_controls":negative},indent=2))

if __name__=="__main__":main()

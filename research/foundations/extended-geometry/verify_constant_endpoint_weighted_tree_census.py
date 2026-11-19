#!/usr/bin/env python3
"""Independent small census via labelled Prufer trees, rotations and edge weights."""
if not __debug__:
    raise RuntimeError("Run without -O: exact assertion controls are required")
from itertools import product,permutations
from math import comb
import argparse,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def tree(code,n):
    degrees=[1]*n
    for v in code:degrees[v]+=1
    edges=[]
    for v in code:
        leaf=next(j for j in range(n)if degrees[j]==1)
        edges.append((leaf,v));degrees[leaf]-=1;degrees[v]-=1
    ends=[j for j in range(n)if degrees[j]==1];edges.append(tuple(ends))
    return edges

def weightings(d,a):
    for v in product(range(1,a+1),repeat=d):
        if sum(v)<=a:yield v

def census(a,r):
    n=r+1;classes=set();labelled=0
    for code in product(range(n),repeat=max(0,n-2)):
        edges=tree(code,n)
        adjacency=[[]for _ in range(n)]
        for eid,(u,v)in enumerate(edges):adjacency[u].append((v,eid));adjacency[v].append((u,eid))
        if len(adjacency[0])!=1 or adjacency[0][0][0]!=1:continue
        parent=[None]*n;parent[0]=-1;depth=[0]*n;order=[0]
        for v in order:
            for u,_ in adjacency[v]:
                if u!=parent[v]:parent[u]=v;depth[u]=depth[v]+1;order.append(u)
        blacks=[v for v in range(1,n)if depth[v]%2==0]
        if any(len(adjacency[v])>a for v in blacks):continue
        labelled+=1
        parentedge=adjacency[0][0][1]
        choices=[list(weightings(len(adjacency[v]),a))for v in blacks]
        rotations=[list(permutations([u for u,_ in adjacency[v]if u!=parent[v]]))for v in range(n)]
        for assignments in product(*choices):
            weights=[None]*r;weights[parentedge]=1
            for v,vals in zip(blacks,assignments):
                for(_,eid),w in zip(adjacency[v],vals):weights[eid]=w
            assert all(w is not None for w in weights)
            edgeid={frozenset((u,v)):i for i,(u,v)in enumerate(edges)}
            for children in product(*rotations):
                def key(v):return tuple((weights[edgeid[frozenset((v,u))]],key(u))for u in children[v])
                classes.add(key(0))
    expected,remainder=divmod(comb(a*r,r-1),r)
    assert remainder==0
    assert len(classes)==expected
    return {'a':a,'r':r,'labelled_topological_trees':labelled,'marked_oriented_weighted_classes':len(classes),'expected':expected}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output','--write-receipt',dest='output',type=Path,default=HERE/'constant_endpoint_weighted_tree_census_verification.json')
    args=parser.parse_args();rows=[]
    for a in (1,2,3,4):
        for r in range(1,6):rows.append(census(a,r))
    out={'passed':True,'method':'Independent small census of ALL labelled trees by Prufer words, marked leaf0/neighbor1, every cyclic rotation system and positive edge assignment satisfying black weighted-degree cap; canonical rooted ordered encodings quotient labels. Finite controls only.','results':rows,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'proof_sha256':hashlib.sha256((HERE/'CONSTANT_ENDPOINT_DAVENPORT_ZANNIER_CLASSIFICATION.md').read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))

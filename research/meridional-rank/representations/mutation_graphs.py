"""Exact small plane-graph and Whitney-switch certificate production."""

from collections import defaultdict
from itertools import combinations,permutations,product
import time


def components(n,edges,removed=()):
    remaining=set(range(n))-set(removed);result=[]
    while remaining:
        reached={min(remaining)};queue=list(reached)
        for v in queue:
            for a,b in edges:
                w=b if a==v else a if b==v else None
                if w is not None and w in remaining and w not in reached:reached.add(w);queue.append(w)
        remaining-=reached;result.append(sorted(reached))
    return result


def checkerboards(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    assert all(len(v)==2 for v in ends.values());pair={a:b for v in ends.values() for a,b in [v,v[::-1]]}
    current=(0,0);walk=[]
    for _ in range(2*len(pd)):
        assert current not in walk;walk.append(current);c,h=current;current=pair[c,(h+2)%4]
    assert current==(0,0)
    assert all(walk[i][1]%2!=walk[(i+1)%len(walk)][1]%2 for i in range(len(walk)))
    where={};faces=[]
    for start in sorted(pair):
        if start in where:continue
        current=start;face=[]
        while current not in where:
            where[current]=len(faces);face.append(current);c,h=pair[current];current=c,(h+1)%4
        assert current==start;faces.append(face)
    assert len(faces)==len(pd)+2
    colours={0:0};queue=[0]
    for f in queue:
        for dart in faces[f]:
            other=where[pair[dart]]
            if other not in colours:colours[other]=1-colours[f];queue.append(other)
            assert colours[other]!=colours[f]
    graphs=[]
    for colour in [0,1]:
        vertices=sorted(f for f in colours if colours[f]==colour);index={f:i for i,f in enumerate(vertices)}
        edges=[];parities=[]
        for c in range(len(pd)):
            h=0 if colours[where[c,0]]==colour else 1;parities.append(h)
            edges.append([index[where[c,h]],index[where[c,(h+2)%4]]])
        assert len(set(parities))==1
        n=len(vertices);assert all(a!=b for a,b in edges) and len(components(n,edges))==1
        assert all(len(components(n,edges[:i]+edges[i+1:]))==1 for i in range(len(edges)))
        graphs.append({'vertices':n,'edges':edges,'checkerboard_parity':parities[0],
            'face_indices':vertices,'two_vertex_connected':all(len(components(n,edges,[v]))==1 for v in range(n))})
    selected=min(range(2),key=lambda i:graphs[i]['vertices'])
    assert graphs[0]['vertices']!=graphs[1]['vertices']
    return {'faces':faces,'face_colours':[colours[i] for i in range(len(faces))],
        'graphs':graphs,'selected_colour':selected,'selected':graphs[selected]}


def canonical(n,edges):
    degrees=[sum(v==a or v==b for a,b in edges) for v in range(n)]
    groups=[[v for v in range(n) if degrees[v]==d] for d in sorted(set(degrees))]
    best=None;best_mapping=None
    for choices in product(*(permutations(g) for g in groups)):
        order=[v for choice in choices for v in choice];mapping={v:i for i,v in enumerate(order)}
        key=tuple(sorted(tuple(sorted((mapping[a],mapping[b]))) for a,b in edges))
        if best is None or key<best:best,best_mapping=key,[mapping[i] for i in range(n)]
    return (n,best),best_mapping


def final_maps(n,edges,target):
    key,a=canonical(n,edges);other,b=canonical(target['vertices'],target['edges']);assert key==other
    inv={v:i for i,v in enumerate(b)};vertices=[inv[a[i]] for i in range(n)]
    available=defaultdict(list)
    for i,(u,v) in enumerate(target['edges']):available[tuple(sorted((u,v)))].append(i)
    edge_map=[]
    for u,v in edges:edge_map.append(available[tuple(sorted((vertices[u],vertices[v])))].pop(0))
    return vertices,edge_map


def switches(n,edges):
    for u,v in combinations(range(n),2):
        pieces=components(n,edges,[u,v])
        if len(pieces)<2:continue
        for piece in pieces:
            side=[i for i,(a,b) in enumerate(edges) if a in piece or b in piece]
            if not side:continue
            left={x for i in side for x in edges[i]};right={x for i,e in enumerate(edges) if i not in side for x in e}
            if left&right!={u,v}:continue
            swap=lambda a:v if a==u else u if a==v else a
            new=[[swap(a),swap(b)] if i in side else [a,b] for i,(a,b) in enumerate(edges)]
            yield new,{'boundary':[u,v],'edge_side':side}


def search(source,target,limits):
    started=time.monotonic();n=source['vertices']
    if n!=target['vertices']:return {'status':'NO_CERTIFICATE','reason':'different graph ranks'}
    target_key,_=canonical(n,target['edges']);key,_=canonical(n,source['edges'])
    states=[source['edges']];parents=[None];moves=[None];seen={key:0}
    current=0;found=None
    while current<len(states):
        if time.monotonic()-started>limits['seconds_per_search']:break
        edges=states[current];key,_=canonical(n,edges)
        if key==target_key:found=current;break
        for new,move in switches(n,edges):
            key,_=canonical(n,new)
            if key in seen:continue
            if len(states)>=limits['maximum_states_per_search']:break
            seen[key]=len(states);states.append(new);parents.append(current);moves.append(move)
        if len(states)>=limits['maximum_states_per_search']:break
        current+=1
    if found is None:return {'status':'NO_CERTIFICATE','states':len(states),'expanded':current,'seconds':round(time.monotonic()-started,6)}
    trace=[];index=found
    while parents[index] is not None:trace.append(moves[index]);index=parents[index]
    trace.reverse();vertices,edges=final_maps(n,states[found],target)
    return {'status':'CERTIFIED_GRAPH_2_ISOMORPHISM','switches':trace,'final_vertex_map':vertices,
        'edge_map':edges,'states':len(states),'expanded':current,'seconds':round(time.monotonic()-started,6)}

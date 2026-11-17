#!/usr/bin/env python3
"""Independent corner/forest/binary-cycle replay; no producer graph imports."""
from collections import defaultdict
from copy import deepcopy
from itertools import combinations
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256((HERE/path).read_bytes()).hexdigest()


def regions_and_graph(pd, claimed):
    n=len(pd); ends=defaultdict(list)
    for c,row in enumerate(pd):
        assert len(row)==4
        for h,label in enumerate(row):ends[label].append((c,h))
    assert len(ends)==2*n and all(len(e)==2 for e in ends.values())
    mate={a:b for e in ends.values() for a,b in (e,e[::-1])}
    # Follow the knot, recording incoming crossing branches.
    walk=[]; current=(0,0)
    for _ in range(2*n):
        assert current not in walk; walk.append(current)
        c,h=current; current=mate[c,(h+2)%4]
    assert current==walk[0]
    assert all(walk[i][1]%2!=walk[(i+1)%(2*n)][1]%2 for i in range(2*n))
    assert all(sum(c==x for x,h in walk)==2 for c in range(n))
    # A corner h is between halfedges h and h+1. Glue corners across edges.
    parent=list(range(4*n))
    def root(a):
        while a!=parent[a]:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    def join(a,b):parent[root(a)]=root(b)
    for (c,h),(d,k) in ends.values():
        join(4*c+(h-1)%4,4*d+k)
        join(4*c+h,4*d+(k-1)%4)
    groups=defaultdict(list)
    for a in range(4*n):groups[root(a)].append(a)
    faces=sorted(groups.values()); assert len(faces)==n+2
    where={corner:i for i,face in enumerate(faces) for corner in face}
    adjacency={i:set() for i in range(len(faces))}
    for c in range(n):
        for h in range(4):
            a,b=where[4*c+h],where[4*c+(h+1)%4]
            assert a!=b;adjacency[a].add(b);adjacency[b].add(a)
    colour={0:0}; queue=[0]
    for a in queue:
        for b in adjacency[a]:
            if b not in colour:colour[b]=1-colour[a];queue.append(b)
            assert colour[b]!=colour[a]
    assert len(colour)==len(faces)
    counts=[sum(v==k for v in colour.values()) for k in (0,1)]
    assert counts[0]!=counts[1]
    chosen=counts.index(min(counts)); vertices=[v for v in range(len(faces)) if colour[v]==chosen]
    index={v:i for i,v in enumerate(vertices)}; edges=[]; parities=[]
    for c in range(n):
        h=0 if colour[where[4*c]]==chosen else 1
        parities.append(1-h)
        edges.append([index[where[4*c+h]],index[where[4*c+(h+2)%4]]])
    assert len(set(parities))==1
    assert claimed['vertices']==len(vertices) and claimed['checkerboard_parity']==parities[0]
    # Align vertices by their complete crossing-edge stars, not canonicalization.
    def stars(es,k):return [frozenset(i for i,e in enumerate(es) if v in e) for v in range(k)]
    actual=stars(edges,len(vertices)); expected=stars(claimed['edges'],len(vertices))
    assert len(set(actual))==len(actual) and set(actual)==set(expected)
    mapping=[expected.index(s) for s in actual]
    assert [sorted(mapping[v] for v in e) for e in edges]==[sorted(e) for e in claimed['edges']]
    trees=spanning_trees(claimed)
    assert trees and all(a!=b for a,b in edges)
    # Every edge occurs in a tree, and every edge is omitted by a tree.
    assert all(any(i in t for t in trees) and any(i not in t for t in trees) for i in range(n))
    return {'vertices':len(vertices),'edges':n,'regions':len(faces),'spanning_trees':len(trees),
        'checkerboard_parity':parities[0],'connected_reduced_alternating_knot':True}


def spanning_trees(graph):
    n=graph['vertices'];edges=graph['edges'];trees=[]
    for selected in combinations(range(len(edges)),n-1):
        parts=[{v} for v in range(n)]
        for i in selected:
            a,b=edges[i]; A=next(s for s in parts if a in s);B=next(s for s in parts if b in s)
            if A is B:break
            parts.remove(A);parts.remove(B);parts.append(A|B)
        else:
            assert len(parts)==1;trees.append(selected)
    return trees


def cycle_code(edges):
    parity=[0]*(1<<len(edges)); result=set()
    for mask in range(len(parity)):
        if mask:
            low=mask&-mask; a,b=edges[low.bit_length()-1]
            parity[mask]=parity[mask^low]^(1<<a)^(1<<b)
        if parity[mask]==0:result.add(mask)
    return result


def replay(source,target,certificate):
    assert certificate['status']=='CERTIFIED_GRAPH_2_ISOMORPHISM'
    n=source['vertices']; assert n==target['vertices']
    edges=deepcopy(source['edges']);m=len(edges);assert m==len(target['edges'])
    original=cycle_code(edges)
    assert len(original)==1<<(m-n+1)
    for move in certificate['switches']:
        u,v=move['boundary'];side=move['edge_side']
        assert 0<=u<n and 0<=v<n and u!=v and side and len(side)<m
        assert len(side)==len(set(side)) and all(0<=i<m for i in side)
        left={a for i in side for a in edges[i]}
        right={a for i,e in enumerate(edges) if i not in side for a in e}
        assert left&right=={u,v}
        exchange={u:v,v:u}
        edges=[[exchange.get(a,a),exchange.get(b,b)] if i in side else [a,b] for i,(a,b) in enumerate(edges)]
        assert cycle_code(edges)==original
    vm=certificate['final_vertex_map'];em=certificate['edge_map']
    assert sorted(vm)==list(range(n)) and sorted(em)==list(range(m))
    for i,(a,b) in enumerate(edges):assert sorted((vm[a],vm[b]))==sorted(target['edges'][em[i]])
    transported={sum(1<<em[i] for i in range(m) if mask>>i&1) for mask in original}
    assert transported==cycle_code(target['edges'])
    return {'switches':len(certificate['switches']),'cycle_vectors':len(original),
        'edge_subsets_checked':1<<m,'edge_map':em}


def main():
    result=json.loads((HERE/'mutation_cover_results.json').read_text())
    plan=json.loads((HERE/'mutation_cover_inputs.json').read_text())
    panel={c['name']:c for c in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    assert result['status']=='COMPLETE_CERTIFICATES'
    for path,digest in result['input_sha256'].items():assert sha(path)==digest,path
    selection={'391':['13a_129','13a_136','13a_194','13a_197'],
        '401':['13a_135','13a_137','13a_196','13a_198']}
    assert plan['selection']==selection
    assert set(result['diagrams'])=={name for names in selection.values() for name in names}
    diagrams={}
    for block,names in selection.items():
        for name in names:
            data=result['diagrams'][name];assert data['pd']==panel[name]['pd']
            diagrams[name]=regions_and_graph(data['pd'],data['selected'])
            assert diagrams[name]['spanning_trees']==int(block)
    expected=[(block,names[0],name) for block,names in selection.items() for name in names[1:]]
    assert [(r['block'],r['source'],r['target']) for r in result['pairs']]==expected
    pairs=[]
    for row in result['pairs']:
        a=result['diagrams'][row['source']]['selected'];b=result['diagrams'][row['target']]['selected']
        assert a['checkerboard_parity']==b['checkerboard_parity']
        assert row['checkerboard_parities']==[a['checkerboard_parity'],b['checkerboard_parity']]
        assert row['requires_mirror_for_same_convention'] is False
        pairs.append({'source':row['source'],'target':row['target'],**replay(a,b,row['certificate'])})
    assert len(result['controls'])==2
    controls=[]
    for control in result['controls']:
        if 'diagrams' in control:
            for d in control['diagrams']:regions_and_graph(d['pd'],d['selected'])
        controls.append({'name':control['name'],**replay(control['source'],control['target'],control['certificate'])})
    assert controls[0]['switches']>0
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    first=result['pairs'][0]; a=result['diagrams'][first['source']]['selected'];b=result['diagrams'][first['target']]['selected']
    bad=deepcopy(first['certificate']);bad['switches'][0]['boundary']=[0,1]
    reject('wrong actual switch boundary',lambda:replay(a,b,bad))
    bad=deepcopy(first['certificate']);bad['edge_map'][0],bad['edge_map'][5]=bad['edge_map'][5],bad['edge_map'][0]
    reject('bijective wrong crossing-edge transport',lambda:replay(a,b,bad))
    bad=deepcopy(first['certificate']);bad['final_vertex_map'][0],bad['final_vertex_map'][1]=bad['final_vertex_map'][1],bad['final_vertex_map'][0]
    reject('wrong final vertex transport',lambda:replay(a,b,bad))
    badpd=deepcopy(result['diagrams'][first['source']]['pd']);badpd[0]=badpd[0][1:]+badpd[0][:1]
    reject('single crossing change breaks alternating convention',lambda:regions_and_graph(badpd,a))
    files=['verify_mutation_covers.py', 'mutation_cover_inputs.json', 'mutation_cover_results.json', 'parabolic_panel.json', 'MUTATION_COMMON_DOUBLE_COVERS.tex']
    receipt={'status':'PASS independent mutation/common-cover inputs','diagrams':diagrams,'pairs':pairs,'controls':controls,
        'corruptions_rejected':rejected,'blocks':selection,'actual_rank3_proved':False,
        'assurance':'Independent corner unions, tree enumeration and all8192 edge-subsets per13-edge graph. No producer graph/search imports. Actual mutation/common-cover theorem additionally uses read Greene Cor4.5 and Morton Thm3, explicitly imported in the written proof; no bridge or meridional-rank mutation invariance assumed.',
        'input_sha256':{path:sha(path) for path in files}}
    (HERE/'mutation_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'diagrams':len(diagrams),'pairs':len(pairs),'controls':len(controls),'corruptions':len(rejected)}))


if __name__=='__main__':main()

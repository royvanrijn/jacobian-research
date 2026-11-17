#!/usr/bin/env python3
"""Independent input replay; the uniform rank/additivity theorem is written."""
from collections import defaultdict
from copy import deepcopy
from itertools import permutations
import hashlib
import json
from pathlib import Path

from verify_integer_lifts import replay

HERE=Path(__file__).resolve().parent


def traced_arcs(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    assert all(len(v)==2 for v in ends.values());visited=set();groups=[]
    for c in range(len(pd)):
        for h in [0,2]:
            if (c,h) in visited:continue
            visited.add((c,h));current=(c,h);edges=[]
            for _ in range(2*len(pd)+1):
                e=pd[current[0]][current[1]];edges.append(e);pair=ends[e]
                d,k=pair[1] if pair[0]==current else pair[0]
                if k%2==0:visited.add((d,k));break
                current=(d,(k+2)%4)
            else:raise AssertionError('arc did not terminate')
            groups.append(sorted(edges))
    assert len(visited)==2*len(pd) and len(groups)==len(pd)
    assert sorted(e for g in groups for e in g)==sorted(ends)
    return sorted(groups)


def knot_and_sphere(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    pair={a:b for v in ends.values() for a,b in [v,v[::-1]]}
    start=(0,0);current=start;seen=set()
    for _ in range(2*len(pd)):
        assert current not in seen;seen.add(current);c,h=current;current=pair[c,(h+2)%4]
    assert current==start
    unseen=set(pair);faces=0
    while unseen:
        start=next(iter(unseen));current=start
        while current in unseen:
            unseen.remove(current);c,h=pair[current];current=(c,(h+1)%4)
        assert current==start;faces+=1
    assert faces==len(pd)+2


def source_factor(name,panel):
    if name=='trefoil':return [[1,4,2,5],[3,6,4,1],[5,2,6,3]],1
    if name=='unknot':return [[1,2,2,1]],0
    return panel[name],2


def check_case(row,expected,panel):
    assert row['factors']==expected['factors'] and row['palettes']==expected['palettes']
    palettes=row['palettes'];R=max(s for p in palettes for s in p)+1
    assert set(range(R))==set(s for p in palettes for s in p)
    assert all(p[:2]==[0,1] and len(set(p))==len(p) for p in palettes)
    pd=row['pd'];knot_and_sphere(pd);labels={int(e):tuple(p) for e,p in row['edge_transpositions'].items()}
    assert all(0<=a<b<R for a,b in labels.values())
    # Undo cuts in each literal source factor; all changed ports belong to edge1.
    offset=0;first=0;marked={};ranks=[]
    for i,name in enumerate(row['factors']):
        source,rank=source_factor(name,panel);ranks.append(rank);palette=palettes[i]
        relabel=row['source_sheet_relabels'][i];assert sorted(relabel)==list(range(len(palette)))
        for local,crossing in enumerate(source):
            for h,e in enumerate(crossing):
                actual=pd[first+local][h]
                if e==1:marked[first+local,h]=i
                else:assert actual==offset+e
                assert set(labels[actual])<=set(palette)
        offset+=2*len(source);first+=len(source)
    assert first==len(pd) and len(marked)==2*len(row['factors'])
    assert all(labels[pd[c][h]]==(0,1) for c,h in marked)
    # Reconstruct lifted gluing graph from sheet equivalence classes, with no rank formula.
    vertices=[];which={}
    for i,p in enumerate(palettes):
        partition=[set(p)]+[{s} for s in range(R) if s not in p]
        for component in partition:
            index=len(vertices);vertices.append([i,'main'] if component==set(p) else [i,next(iter(component))])
            for s in component:which[i,s]=index
    lifted=[]
    for i in range(len(palettes)-1):
        for sphere in [[0,1]]+[[s] for s in range(2,R)]:lifted.append([which[i,sphere[0]],which[i+1,sphere[0]],sphere])
    assert row['gluing_graph']=={'degree':R,'vertices':vertices,'edges':lifted,'cycle_rank':expected['expected_free']}
    # Each main-component H1 is finite in these controls; graph cycles give all free H1.
    parent=list(range(len(vertices)));cycles=0
    def find(v):
        while parent[v]!=v:v=parent[v]
        return v
    for a,b,_ in lifted:
        a,b=find(a),find(b)
        if a==b:cycles+=1
        else:parent[a]=b
    assert len({find(v) for v in range(len(vertices))})==1 and cycles==expected['expected_free']
    groups=traced_arcs(pd);cover=row['cover'];assert groups==cover['arc_edges']
    where={e:i for i,g in enumerate(groups) for e in g};pairs=[labels[g[0]] for g in groups]
    assert all(all(labels[e]==pairs[i] for e in g) for i,g in enumerate(groups))
    assert [list(p) for p in pairs]==cover['arc_transpositions']
    # Check the supplied gauge tree and reconstruct inverse lifted paths pointwise.
    tree=cover['tree_arc_indices'];forest=list(range(R))
    def root(v):
        while forest[v]!=v:v=forest[v]
        return v
    for g in tree:
        a,b=pairs[g];assert root(a)!=root(b);forest[root(a)]=root(b)
    assert len(tree)==R-1
    retained=[i for i in range(len(groups)) if i not in tree];assert retained==cover['retained_arc_indices']
    index={g:i for i,g in enumerate(retained)};matrix=[];words=[]
    for a,b,c,d in pd:
        assert where[b]==where[d];o,u,v=where[b],where[a],where[c]
        # Direct conjugation of transposition endpoints checks the monodromy relation.
        x,y=pairs[o]
        moved=lambda s:y if s==x else x if s==y else s
        assert tuple(sorted(moved(s) for s in pairs[u]))==pairs[v]
        for start in range(R):
            sheet=start;word=[];row_vector=[0]*len(retained)
            for g in [v,o,u,o]:
                p,q=pairs[g]
                if sheet not in [p,q]:continue
                sign=1 if sheet==p else -1
                if g in index:word.append(sign*(index[g]+1));row_vector[index[g]]+=sign
                sheet=q if sheet==p else p
            assert sheet==start
            matrix.append([-x for x in row_vector]);words.append([-g for g in word[::-1]])
    assert matrix==cover['homology_matrix'] and words==cover['lifted_relators']
    diagonal=replay(matrix,cover['smith'])
    free=len(retained)-len(diagonal);torsion=[d for d in diagonal if d>1]
    assert free==cover['H1_free_rank']==expected['expected_free']
    assert torsion==cover['H1_torsion']==expected['expected_torsion']
    rank=sum(ranks)+cycles;lower=R-1+rank
    assert rank==row['rank_from_uniform_topology_and_Grushko']==expected['expected_cover_rank']
    assert lower==row['meridional_lower']==expected['expected_lower']
    return {'id':row['id'],'degree':R,'crossings':len(pd),'sphere_graph_cycles':cycles,
            'H1_free_rank':free,'H1_torsion':torsion,'conditional_cover_rank':rank,'meridional_lower':lower}


def check_torus(row):
    def permutation(a):
        return [sum((sum(a[i][j]*((v>>j)&1) for j in range(3))%2)<<i for i in range(3)) for v in range(8)]
    x,y,t,z=[permutation(row[k]) for k in ['X3','Y3','meridian3','Z3']]
    assert all(sorted(p)==list(range(8)) for p in [x,y,t,z])
    def at(word,start,images):
        for g in word[::-1]:start=images[abs(g)-1][start] if g>0 else images[-g-1].index(start)
        return start
    assert all(at([1]*4,i,[x])==at([1]*7,i,[y])==at([1]*2,i,[t])==i for i in range(8))
    assert t==[at([-1,2,2],i,[x,y]) for i in range(8)]
    assert z==[at([1,2,1],i,[t,y]) for i in range(8)]
    assert any(y[z[i]]!=z[y[i]] for i in range(8))
    relator=row['orbifold_relator'];relations=[]
    for coset in [0,1]:
        state=coset;word=[]
        for g in relator:
            if abs(g)==2:state^=1
            else:word.append((state+1)*(1 if g>0 else -1))
        assert state==coset;relations.append(word)
    assert relator==[1,1,2]*4+[-1]*7 and relations==row['double_cover_relators']
    assert all(at(w,i,[y,z])==i for w in relations for i in range(8))
    matrix=[[sum((1 if g>0 else -1) for g in w if abs(g)==j) for j in [1,2]] for w in relations]
    assert matrix==row['exponent_matrix']==[[-3,4],[4,-3]]
    assert abs(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0])==7
    counts={}
    for r in [3,4,5]:counts[r]=sum(all(at([1]*7,i,[p])==i for i in range(r)) for p in permutations(range(r)))
    assert counts=={3:1,4:1,5:1}
    assert row['double_cover_rank']==2 and row['maximum_simple_cover_rank_bound']==3 and row['known_meridional_rank_and_bridge']==4
    return {'double_cover_rank':2,'maximum_simple_cover_rank_bound':3,'meridional_rank':4,'seventh_power_identity_counts':counts,'nonabelian_matrix_witness_checked':True}


def main():
    data=json.loads((HERE/'cover_sum_results.json').read_text());plan=json.loads((HERE/'cover_sum_inputs.json').read_text())
    for p,h in data['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r['pd'] for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    assert len(data['controls'])==len(plan['controls'])==8
    results=[check_case(row,expected,panel) for row,expected in zip(data['controls'],plan['controls'])]
    boundary=check_torus(data['T47_boundary'])
    # Mathematically relevant failures: incompatible meridian palette, lost graph cycle, false Smith step, false torus image.
    bad=[]
    r=deepcopy(data['controls'][1]);e=deepcopy(plan['controls'][1]);r['palettes'][1]=[0,2,3];e['palettes']=deepcopy(r['palettes']);bad.append((r,e))
    r=deepcopy(data['controls'][1]);e=deepcopy(plan['controls'][1]);r['gluing_graph']['cycle_rank']=0;e['expected_free']=0;bad.append((r,e))
    r=deepcopy(data['controls'][0]);step=next(s for s in r['cover']['smith']['steps'] if s[1]=='add');step[4]+=1;bad.append((r,plan['controls'][0]))
    for row,expected in bad:
        try:check_case(row,expected,panel)
        except AssertionError:pass
        else:raise AssertionError('invalid sum input accepted')
    r=deepcopy(data['T47_boundary']);r['Z3']=[[int(i==j) for j in range(3)] for i in range(3)]
    try:check_torus(r)
    except AssertionError:pass
    else:raise AssertionError('invalid torus witness accepted')
    paths=['verify_cover_sums.py', 'cover_sum_inputs.json', 'cover_sum_results.json', 'joint_a5_panel.json', 'verify_integer_lifts.py', 'integer_lifts.py', 'SIMPLE_COVER_RANK_ADDITIVITY.tex']
    report={'status':'PASS independent simple-cover sum inputs','controls':results,'control_count':8,'T47_boundary':boundary,
            'rejected_corruptions':4,'scope':'Independent finite diagram/monodromy/lifted-word/Smith and torus-boundary inputs. The all-knots additivity theorem requires the written normal-generation, lifted-sphere topology and imported Grushko theorem; finite controls do not prove those uniform statements. General meridional additivity, worldwide novelty and full MRC remain open.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'cover_sum_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'controls':results,'T47_boundary':boundary,'rejected_corruptions':4},indent=2))


if __name__=='__main__':main()

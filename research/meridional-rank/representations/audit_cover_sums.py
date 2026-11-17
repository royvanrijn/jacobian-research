#!/usr/bin/env python3
"""Exact diagram controls for the written simple-cover connected-sum theorem."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import time

from integer_lifts import smith_certificate
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent
ROOTS=[(0,1),(0,2),(1,2)]


def arc_data(pd):
    parent={e:e for c in pd for e in c}
    def root(e):
        while parent[e]!=e:e=parent[e]
        return e
    for _,a,_,b in pd:parent[root(a)]=root(b)
    groups=defaultdict(list)
    for e in sorted(parent):groups[root(e)].append(e)
    groups=sorted(groups.values());where={e:i for i,g in enumerate(groups) for e in g}
    return groups,[(where[b],where[a],where[c]) for a,b,c,d in pd]


def source(name,panel,old):
    if name=='trefoil':
        pd=[[1,4,2,5],[3,6,4,1],[5,2,6,3]];return pd,{i:(0,1) for i in range(1,7)},1,[0,1]
    if name=='unknot':return [[1,2,2,1]],{1:(0,1),2:(0,1)},0,[0,1]
    pd=panel[name]['pd'];cr=oriented_diagram(pd);labels=old[name]['arc_labels'];edges={}
    for (a,b,c,d),(o,i,j,_) in zip(pd,cr):
        for e,g in [(a,i),(b,o),(c,j),(d,o)]:
            pair=ROOTS[labels[g]]
            if e in edges:assert edges[e]==pair
            edges[e]=pair
    u,v=edges[1];other=next(s for s in range(3) if s not in [u,v])
    relabel={u:0,v:1,other:2}
    edges={e:tuple(sorted(relabel[s] for s in p)) for e,p in edges.items()}
    assert edges[1]==(0,1)
    return pd,edges,2,[relabel[i] for i in range(3)]


def join(pd,labels,other,other_labels):
    offset=max(labels);shifted=[[e+offset for e in c] for c in other]
    shifted_labels={e+offset:p for e,p in other_labels.items()};a,b=1,offset+1
    assert labels[a]==shifted_labels[b]==(0,1)
    left=[c[:] for c in pd];right=[c[:] for c in shifted]
    l=[(i,h) for i,c in enumerate(left) for h,e in enumerate(c) if e==a]
    r=[(i,h) for i,c in enumerate(right) for h,e in enumerate(c) if e==b]
    assert len(l)==len(r)==2
    i,h=l[1];left[i][h]=b;i,h=r[1];right[i][h]=a
    return left+right,{**labels,**shifted_labels}


def graph_of_pieces(palettes):
    R=max(max(p) for p in palettes)+1;vertices=[]
    for i,p in enumerate(palettes):
        vertices.append([i,'main'])
        vertices.extend([i,s] for s in range(R) if s not in p)
    lookup={tuple(v):i for i,v in enumerate(vertices)}
    def vertex(i,s):return lookup[(i,'main' if s in palettes[i] else s)]
    edges=[]
    for i in range(len(palettes)-1):
        edges.append([vertex(i,0),vertex(i+1,0),[0,1]])
        edges.extend([vertex(i,s),vertex(i+1,s),[s]] for s in range(2,R))
    seen={0};queue=[0]
    for v in queue:
        for a,b,_ in edges:
            w=b if a==v else a if b==v else None
            if w is not None and w not in seen:seen.add(w);queue.append(w)
    assert len(seen)==len(vertices)
    beta=len(edges)-len(vertices)+1
    assert beta==sum(len(p)-2 for p in palettes)-(R-2)
    return {'degree':R,'vertices':vertices,'edges':edges,'cycle_rank':beta}


def folded_cover(pd,edge_labels,R):
    groups,crossings=arc_data(pd);pairs=[]
    for g in groups:
        p=edge_labels[g[0]];assert all(edge_labels[e]==p for e in g);pairs.append(p)
    parent=list(range(R));tree=[]
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    for i,(a,b) in enumerate(pairs):
        if root(a)!=root(b):parent[root(a)]=root(b);tree.append(i)
    assert len(tree)==R-1
    retained=[i for i in range(len(groups)) if i not in tree];index={g:i for i,g in enumerate(retained)}
    words=[];matrix=[]
    for o,a,b in crossings:
        for start in range(R):
            v=start;word=[];row=[0]*len(retained)
            for g in [o,a,o,b]:
                x,y=pairs[g]
                if v==x or v==y:
                    sign=1 if v==x else -1
                    if g in index:word.append(sign*(index[g]+1));row[index[g]]+=sign
                    v=y if v==x else x
            assert v==start;words.append(word);matrix.append(row)
    smith=smith_certificate(matrix,len(retained))
    return {'arc_edges':groups,'arc_transpositions':pairs,'tree_arc_indices':tree,'retained_arc_indices':retained,
            'lifted_relators':words,'homology_matrix':matrix,'smith':smith,
            'H1_free_rank':len(retained)-len(smith['diagonal']),'H1_torsion':[d for d in smith['diagonal'] if d>1]}


def matmul(a,b):return [[sum(x*y for x,y in zip(r,c))%2 for c in zip(*b)] for r in a]
def matpow(a,n):
    out=[[int(i==j) for j in range(len(a))] for i in range(len(a))]
    for _ in range(n):out=matmul(out,a)
    return out


def torus_boundary():
    previous=json.loads((HERE/'gl3f2_verification.json').read_text())['T47_control']
    x=[r[:3] for r in previous['homogeneous_X'][:3]];y=[r[:3] for r in previous['homogeneous_Y'][:3]]
    t=matmul(matpow(x,3),matpow(y,2));z=matmul(matmul(t,y),t);identity=matpow(x,0)
    assert matpow(x,4)==matpow(y,7)==matpow(t,2)==identity
    assert matmul(y,z)!=matmul(z,y)
    # t is a sheet swap; y is even. Rewrite the one orbifold relation from both cosets.
    orbifold=[1,1,2]*4+[-1]*7;relations=[]
    for start in [0,1]:
        state=start;word=[]
        for g in orbifold:
            if abs(g)==2:state=1-state
            else:word.append((1 if state==0 else 2)*(1 if g>0 else -1))
        assert state==start;relations.append(word)
    matrices={1:y,2:z};inverses={1:matpow(y,6),2:matpow(z,6)}
    for w in relations:
        value=identity
        for g in w:value=matmul(value,matrices[g] if g>0 else inverses[-g])
        assert value==identity
    exponents=[[sum((1 if g>0 else -1) for g in w if abs(g)==j) for j in [1,2]] for w in relations]
    assert exponents==[[-3,4],[4,-3]]
    return {'torus_group_relator':[1]*4+[-2]*7,'meridian_word':[-1,2,2],
            'orbifold_generators':['y','t'],'orbifold_relator':orbifold,'double_cover_generators':['y','t*y*t'],
            'double_cover_relators':relations,'exponent_matrix':exponents,'X3':x,'Y3':y,'meridian3':t,'Z3':z,
            'noncommuting_pair':['Y3','Z3'],'double_cover_rank':2,'maximum_simple_cover_rank_bound':3,
            'known_meridional_rank_and_bridge':4,
            'higher_degree_exclusion':'Any onto S_r simple coloring has r<=5 from upper4. Central x^4=y^7 maps to1 for r>=3, hence y^7=1; no S_r with r<=5 has an element of order7, so y=1 and the image would be cyclic.'}


def main():
    started=time.monotonic();plan=json.loads((HERE/'cover_sum_plan.json').read_text())
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    old={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    rows=[]
    for case in plan['controls']:
        pd=[];labels={};ranks=[];transports=[]
        for name,palette in zip(case['factors'],case['palettes']):
            base,colours,rank,relabel=source(name,panel,old);ranks.append(rank);transports.append(relabel)
            assert len(palette)==len(relabel)
            assert palette[:2]==[0,1] and len(set(palette))==len(palette)
            embedded={e:tuple(sorted(palette[s] for s in p)) for e,p in colours.items()}
            if not pd:pd=[c[:] for c in base];labels=embedded
            else:pd,labels=join(pd,labels,base,embedded)
        topology=graph_of_pieces(case['palettes']);R=topology['degree'];cover=folded_cover(pd,labels,R)
        assert cover['H1_free_rank']==case['expected_free'] and cover['H1_torsion']==case['expected_torsion']
        rank=sum(ranks)+topology['cycle_rank'];lower=R-1+rank
        assert rank==case['expected_cover_rank'] and lower==case['expected_lower']
        rows.append({'id':case['id'],'factors':case['factors'],'palettes':case['palettes'],'source_sheet_relabels':transports,'pd':pd,
                     'edge_transpositions':labels,'gluing_graph':topology,'cover':cover,
                     'rank_from_uniform_topology_and_Grushko':rank,'meridional_lower':lower,
                     'finite_check_boundary':'Homology and finite gluing graph checked; actual group-rank formula requires the written topology and classical Grushko theorem.'})
        assert time.monotonic()-started<plan['limits']['seconds_total']
    paths=['cover_sum_plan.json','audit_cover_sums.py','integer_lifts.py','parabolic_lifts.py','joint_a5_panel.json',
           'three_sheet_nonabelian_results.json','three_sheet_splice_verification.json','periodic_twist_family_verification.json','gl3f2_verification.json']
    result={'status':'CONSTRUCTED connected-sum cover inputs','controls':rows,'T47_boundary':torus_boundary(),
            'scope':'Finite input controls for a separate uniform all-colorings additivity proof; no full MRC or general meridional-additivity conclusion.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'cover_sum_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'controls':[{k:r[k] for k in ['id','rank_from_uniform_topology_and_Grushko','meridional_lower']}|{'degree':r['gluing_graph']['degree'],'H1_free':r['cover']['H1_free_rank'],'H1_torsion':r['cover']['H1_torsion']} for r in rows],'T47_C':3},indent=2))


if __name__=='__main__':main()

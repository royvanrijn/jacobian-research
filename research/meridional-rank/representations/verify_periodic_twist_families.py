#!/usr/bin/env python3
"""Independent finite inputs for the written uniform periodic-twist theorem."""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
import hashlib
import json
from pathlib import Path
import re

from parabolic_lifts import oriented_diagram  # old certificate indices only

HERE=Path(__file__).resolve().parent


def word_at(word,point,images):
    for g in word[::-1]:
        p=images[abs(g)-1]
        point=p[point] if g>0 else p.index(point)
    return point


def components(pd,removed=None):
    graph={e:set() for row in pd for e in row}
    for c,(_,b,_,d) in enumerate(pd):
        if c!=removed:graph[b].add(d);graph[d].add(b)
    unseen=set(graph);groups=[]
    while unseen:
        queue=[min(unseen)];found=set(queue)
        for e in queue:
            for f in graph[e]-found:found.add(f);queue.append(f)
        unseen-=found;groups.append(sorted(found))
    return sorted(groups)


def trace_check(pd,row):
    c=row['selected_crossing_zero_based'];axis=row['axis'];assert (c,axis)==(0,0)
    groups=components(pd,c);assert groups==row['outside_arc_edges']
    where={e:i for i,g in enumerate(groups) for e in g}
    boundary=[where[e] for e in pd[c]];assert boundary==row['boundary_arc_indices']
    known=set(row['four_outside_seeds']);assert len(known)==4 and known<=set(range(len(groups)))
    tangle=0
    for step in row['uniform_upper_trace']:
        if step[0]=='integer_tangle':
            _,ports,added=step;assert ports in [[0,1],[2,3]] and tangle==0
            assert all(boundary[h] in known for h in ports)
            assert added==sorted(set(boundary)-known);known.update(boundary);tangle+=1
        else:
            kind,j,o,a,b=step;assert kind=='crossing' and j!=c
            x,y,z,w=pd[j];assert o==where[y]==where[w]
            assert {a,b}=={where[x],where[z]} and o in known and a in known and b not in known
            known.add(b)
    assert known==set(range(len(groups))) and tangle==1
    return groups


def corner_graph(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    assert all(len(v)==2 for v in ends.values())
    parent={(c,h):(c,h) for c in range(len(pd)) for h in range(4)}
    def find(a):
        while a!=parent[a]:a=parent[a]
        return a
    def union(a,b):parent[find(a)]=find(b)
    for (c,h),(d,k) in ends.values():
        union((c,(h-1)%4),(d,k));union((c,h),(d,(k-1)%4))
    faces={find(a) for a in parent};assert len(faces)==len(pd)+2
    adj={f:set() for f in faces}
    for c in range(len(pd)):
        for h in range(4):
            a,b=find((c,h)),find((c,(h+1)%4));adj[a].add(b);adj[b].add(a)
    colour={find((0,3)):0};queue=list(colour)
    for a in queue:
        for b in adj[a]:
            if b not in colour:colour[b]=1-colour[a];queue.append(b)
            assert colour[b]!=colour[a]
    vertices=sorted(a for a in faces if colour[a]==0);index={a:i for i,a in enumerate(vertices)};edges=[]
    for c in range(len(pd)):
        p=0 if colour[find((c,0))]==0 else 1
        edges.append([index[find((c,p))],index[find((c,(p+2)%4))]])
    return list(range(len(vertices))),edges


def incidence(vertices,edges):
    return sorted(tuple(i for i,(a,b) in enumerate(edges) for v in [a,b] if v==u) for u in vertices)


def tree_count(vertices,edges):
    count=0
    for choice in combinations(edges,len(vertices)-1):
        parent={u:u for u in vertices}
        def find(u):
            while u!=parent[u]:u=parent[u]
            return u
        for a,b in choice:
            a,b=find(a),find(b)
            if a==b:break
            parent[a]=b
        else:count+=1
    return count


def rational_det(matrix):
    a=[list(map(Fraction,row)) for row in matrix];value=Fraction(1)
    for j in range(len(a)):
        k=next((k for k in range(j,len(a)) if a[k][j]),None)
        if k is None:return 0
        if k!=j:a[k],a[j]=a[j],a[k];value=-value
        pivot=a[j][j];value*=pivot
        for i in range(j+1,len(a)):
            factor=a[i][j]/pivot
            for k in range(j,len(a)):a[i][k]-=factor*a[j][k]
    assert value.denominator==1
    return abs(value.numerator)


def fox_determinant(pd):
    groups=components(pd);assert len(groups)==len(pd);where={e:i for i,g in enumerate(groups) for e in g}
    matrix=[]
    for a,b,c,d in pd:
        assert where[b]==where[d];v=[0]*len(groups);v[where[b]]+=2;v[where[a]]-=1;v[where[c]]-=1;matrix.append(v)
    return rational_det([r[:-1] for r in matrix[:-1]])


def alternating_knot(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    start=(0,0);current=start;seen=[]
    for _ in range(2*len(pd)):
        assert current not in seen;seen.append(current);c,h=current
        outgoing=(c,(h+2)%4);pair=ends[pd[c][outgoing[1]]];current=pair[1] if pair[0]==outgoing else pair[0]
    assert current==start and len(seen)==2*len(pd)
    assert all(seen[i][1]%2!=seen[(i+1)%len(seen)][1]%2 for i in range(len(seen)))


def check_family(row,pd):
    images=row['knot_meridian_images'];assert len(images)==13
    identity=list(range(75));block_graph=[set() for _ in range(3)]
    for p in images:
        assert sorted(p)==identity and all(p[p[i]]==i for i in identity)
        blocks=[{p[25*s+i]//25 for i in range(25)} for s in range(3)]
        assert all(len(b)==1 for b in blocks);sigma=[next(iter(b)) for b in blocks]
        assert sorted(sigma)==[0,1,2] and sum(sigma[s]!=s for s in range(3))==2
        assert all(p[i]==i for i in identity if sigma[i//25]==i//25)
        for s in range(3):block_graph[s].add(sigma[s])
    reached={0};queue=[0]
    for s in queue:
        for t in block_graph[s]-reached:reached.add(t);queue.append(t)
    assert reached=={0,1,2}
    crossings=oriented_diagram(pd);edge_images={}
    for (a,b,c,d),(o,x,y,_) in zip(pd,crossings):
        for e,g in [(a,x),(b,o),(c,y),(d,o)]:
            if e in edge_images:assert edge_images[e]==images[g]
            edge_images[e]=images[g]
        assert all(word_at([o+1,x+1,o+1,y+1],i,images)==i for i in identity)
    # Every listed meridian occurs outside the selected crossing.
    c=row['selected_crossing_zero_based'];outer={e for j,r in enumerate(pd) if j!=c for e in r}
    assert all(p in [edge_images[e] for e in outer] for p in images)
    witness=row['stabilizer_noncommutator'];wa,wb=witness['stabilizer_words']
    assert all(word_at(w,i,images)<25 for w in [wa,wb] for i in range(25))
    comm=wa+wb+[-g for g in wa[::-1]]+[-g for g in wb[::-1]]
    assert witness['commutator_word']==comm
    assert word_at(comm,witness['moved_point'],images)==witness['image']!=witness['moved_point']
    groups=trace_check(pd,row);axis=row['axis'];assert row['source_crossing_ports']==pd[c]
    x,y=[edge_images[pd[c][h]] for h in [axis,axis+1]];m=row['product_order']
    assert m>0 and all(word_at([1,2]*m,i,[x,y])==i for i in identity)
    assert m==3 and any(word_at([1,2],i,[x,y])!=i for i in identity)
    vertices,edges=corner_graph(pd)
    assert incidence(vertices,edges)==incidence(row['checkerboard_vertices'],row['checkerboard_edges'])
    u,v=edges[c];others=edges[:c]+edges[c+1:]
    deleted=tree_count(vertices,others)
    contracted=tree_count([w for w in vertices if w!=v],[[u if a==v else a,u if b==v else b] for a,b in others])
    assert deleted==row['deleted_edge_tree_count'] and contracted==row['contracted_edge_tree_count']
    assert row['determinant_constant']==deleted+contracted==fox_determinant(pd)
    assert row['determinant_slope_in_twist_length']==contracted>0
    assert row['determinant_slope_in_N']==2*m*contracted
    assert [sample['N'] for sample in row['samples']]==[0,1,2]
    samples=[]
    for sample in row['samples']:
        N=sample['N'];n=1+2*m*N;new=sample['pd'];assert sample['twist_length']==n
        assert new[:c]==pd[:c] and new[c+n:]==pd[c+1:] and len(new)==12+n==sample['crossings']
        p=list(pd[c]);inside=dict(edge_images);left,right=x[:],y[:]
        ls=[p[0]]+[26+i for i in range(1,n)]+[p[3]]
        rs=[p[1]]+[26+n-1+i for i in range(1,n)]+[p[2]]
        assert new[c:c+n]==[[ls[i],rs[i],rs[i+1],ls[i+1]] for i in range(n)]
        for j in range(n):
            a,b,z,d=new[c+j]
            assert inside.get(a,left)==left and inside.get(b,right)==right
            inside[a]=left;inside[b]=right;new_right=[right[left[right[i]]] for i in identity]
            assert inside.get(z,new_right)==new_right and inside.get(d,right)==right
            inside[z]=new_right;inside[d]=right;left,right=right,new_right
        assert left==edge_images[p[3]] and right==edge_images[p[2]]
        assert new[c][0:2]==p[0:2] and new[c+n-1][2:4]==p[2:4]
        for a,b,z,d in new:
            assert inside[b]==inside[d]
            assert all(inside[z][i]==inside[b][inside[a][inside[b][i]]] for i in identity)
        alternating_knot(new)
        vs,es=corner_graph(new);expected=edges[:c]+[edges[c]]*n+edges[c+1:]
        assert incidence(vs,es)==incidence(vertices,expected)  # parallel-edge expansion
        det=fox_determinant(new);assert det==sample['determinant']==deleted+n*contracted
        arcs=components(new);where={e:i for i,g in enumerate(arcs) for e in g}
        known={where[groups[s][0]] for s in row['four_outside_seeds']};changed=True
        while changed:
            changed=False
            for a,b,z,d in new:
                o,i,j=where[b],where[a],where[z]
                if o in known and bool({i,j}&known):
                    added={i,j}-known;known.update(added);changed|=bool(added)
        assert len(known)==len(arcs)
        samples.append({'N':N,'crossings':len(new),'determinant':det,'all_meridian_relators_and_four_seeds_pass':True})
    return {'name':row['name'],'period_added_crossings':2*m,'meridional_lower':4,'bridge_upper':4,
            'determinant_constant':deleted+contracted,'determinant_slope_in_N':2*m*contracted,
            'stabilizer_moved_point':witness['moved_point'],'stabilizer_image':witness['image'],'samples':samples}


def main():
    data=json.loads((HERE/'periodic_twist_family_results.json').read_text())
    for p,h in data['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r['pd'] for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    names=['13a_3145','13a_3147','13a_3148','13a_3149'];assert [r['name'] for r in data['families']]==names
    proof=(HERE/'PERIODIC_TWIST_FAMILIES.tex').read_text()
    literal=proof.split('All four PDs have the common seven-tuple prefix',1)[1]
    prefix_text,suffix_text=literal.split('\\begin{center}',1)
    tuples=lambda text:[list(map(int,w.split(','))) for w in re.findall(r'\((\d+,\d+,\d+,\d+)\)',text)]
    prefix=tuples(prefix_text);suffix=tuples(suffix_text.split('\\bottomrule',1)[0])
    assert len(prefix)==7 and len(suffix)==24
    assert all(prefix+suffix[i::4]==panel[name] for i,name in enumerate(names))
    results=[check_family(r,panel[r['name']]) for r in data['families']]
    original=data['families'][0];bad=[]
    r=deepcopy(original);r['knot_meridian_images'][0][0],r['knot_meridian_images'][0][1]=r['knot_meridian_images'][0][1],r['knot_meridian_images'][0][0];bad.append(r)
    r=deepcopy(original);r['product_order']=1;bad.append(r)
    r=deepcopy(original);r['uniform_upper_trace'][0][1]=[0,3];bad.append(r)
    r=deepcopy(original);r['determinant_slope_in_N']+=1;bad.append(r)
    for r in bad:
        try:check_family(r,panel[r['name']])
        except AssertionError:pass
        else:raise AssertionError('corrupted family certificate accepted')
    paths=['periodic_twist_family_plan.json','audit_periodic_twist_families.py','periodic_twist_family_results.json',
           'verify_periodic_twist_families.py','PERIODIC_TWIST_FAMILIES.tex','joint_a5_panel.json','finite_cover_quotient_verification.json','parabolic_lifts.py']
    report={'status':'PASS independent periodic twist family inputs','families':results,'family_count':4,
            'explicit_expanded_diagrams':12,'knot_action_degree':75,'blocks':3,'block_size':25,'rejected_corruptions':len(bad),
            'uniform_proof_boundary':'Finite checks certify specific input permutations, diagrams, uniform outside traces and determinant coefficients. All-N conclusions require the written spanning-tree gauge, Hurwitz-period and integer-tangle induction proofs, not interpolation from the12 samples. Named seed provenance retained; family mutual disjointness and worldwide novelty unclaimed. Unformalized and externally unreviewed.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'periodic_twist_family_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'families':results,'rejected_corruptions':len(bad)},indent=2))


if __name__=='__main__':main()

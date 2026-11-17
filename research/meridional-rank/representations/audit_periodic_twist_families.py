#!/usr/bin/env python3
"""Construct finite knot actions and uniform four-seed twist-family templates."""
from collections import defaultdict, deque
from itertools import combinations
from math import lcm
import hashlib
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent
ROOTS=[(0,1),(0,2),(1,2)]


def compose(a,b):return [a[b[i]] for i in range(len(a))]
def inverse(a):return [a.index(i) for i in range(len(a))]
def word_value(word,images):
    out=list(range(len(images[0])))
    for g in word:out=compose(out,images[abs(g)-1] if g>0 else inverse(images[-g-1]))
    return out
def order(p):
    unseen=set(range(len(p)));answer=1
    while unseen:
        start=next(iter(unseen));i=start;length=0
        while i in unseen:unseen.remove(i);length+=1;i=p[i]
        assert i==start;answer=lcm(answer,length)
    return answer


def induced(old,full,crossings):
    identity=list(range(25));voltages={int(g):p for g,p in full['original_cover_generator_images'].items()}
    for sign in [1,-1]:
        images=[]
        for arc,label in enumerate(old['arc_labels']):
            u,v=ROOTS[label];p=voltages.get(arc+1,identity)
            if sign<0:p=inverse(p)
            inv=inverse(p);image=[]
            for sheet in range(3):
                for point in range(25):
                    image.append(25*v+p[point] if sheet==u else 25*u+inv[point] if sheet==v else 25*sheet+point)
            images.append(image)
        if all(word_value([o+1,a+1,o+1,b+1],images)==list(range(75)) for o,a,b,_ in crossings):break
    else:raise AssertionError('neither voltage convention satisfies the source relators')
    paths={0:[]};queue=[0]
    for u in queue:
        for g in old['tree_generators']:
            a,b=ROOTS[old['arc_labels'][g-1]]
            if u in [a,b]:
                v=b if u==a else a
                if v not in paths:paths[v]=[g]+paths[u];queue.append(v)
    assert set(paths)=={0,1,2}
    words=[]
    for g in old['remaining_generators']:
        u,v=ROOTS[old['arc_labels'][g-1]]
        words.append(paths[v][::-1]+[g]+paths[u])
    comm=words[0]+words[1]+[-g for g in words[0][::-1]]+[-g for g in words[1][::-1]]
    values=[word_value(w,images) for w in words]
    assert all(all(v[i]<25 for i in range(25)) for v in values)
    c=word_value(comm,images);point=next(i for i in range(25) if c[i]!=i)
    return images,{'voltage_sign':sign,'stabilizer_words':words,'commutator_word':comm,'moved_point':point,'image':c[point]}


def outside_arcs(pd,removed):
    edges=sorted({e for c in pd for e in c});parent={e:e for e in edges}
    def root(e):
        while e!=parent[e]:e=parent[e]
        return e
    for i,(_,b,_,d) in enumerate(pd):
        if i!=removed:parent[root(b)]=root(d)
    groups=defaultdict(list)
    for e in edges:groups[root(e)].append(e)
    groups=sorted(groups.values());where={e:i for i,g in enumerate(groups) for e in g}
    crossings=[(i,where[b],where[a],where[c]) for i,(a,b,c,d) in enumerate(pd) if i!=removed]
    return groups,crossings,[where[e] for e in pd[removed]]


def propagate(count,crossings,boundary,axis,seeds):
    known=set(seeds);trace=[];used=False;changed=True
    pairs=[(axis,(axis+1)%4),((axis+2)%4,(axis+3)%4)]
    while changed:
        changed=False
        for c,o,a,b in crossings:
            if o in known:
                if a in known and b not in known:known.add(b);trace.append(['crossing',c,o,a,b]);changed=True
                if b in known and a not in known:known.add(a);trace.append(['crossing',c,o,b,a]);changed=True
        if not used:
            pair=next((p for p in pairs if all(boundary[i] in known for i in p)),None)
            if pair is not None:
                added=sorted(set(boundary)-known);known.update(boundary);used=True;changed=True
                trace.append(['integer_tangle',list(pair),added])
    return trace if len(known)==count and used else None


def checkerboard(pd):
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    assert all(len(v)==2 for v in ends.values());pair={a:b for v in ends.values() for a,b in [v,v[::-1]]}
    face={};faces=[]
    for d in pair:
        if d in face:continue
        current=d;path=[]
        while current not in face:
            face[current]=len(faces);path.append(current);c,h=pair[current];current=(c,(h+1)%4)
        assert current==d;faces.append(path)
    assert len(faces)==len(pd)+2
    colours={0:0};queue=[0]
    for f in queue:
        for d in faces[f]:
            other=face[pair[d]]
            if other not in colours:colours[other]=1-colours[f];queue.append(other)
            assert colours[other]!=colours[f]
    vertices=sorted(f for f in colours if colours[f]==0);edges=[]
    for c in range(len(pd)):
        parity=0 if colours[face[c,0]]==0 else 1
        edges.append([face[c,parity],face[c,(parity+2)%4]])
    return vertices,edges


def determinant(matrix):
    a=[r[:] for r in matrix];n=len(a);previous=1;sign=1
    if not n:return 1
    for k in range(n-1):
        if not a[k][k]:
            row=next((i for i in range(k+1,n) if a[i][k]),None)
            if row is None:return 0
            a[k],a[row]=a[row],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*pivot-a[i][k]*a[k][j];assert value%previous==0;a[i][j]=value//previous
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def trees(vertices,edges):
    index={v:i for i,v in enumerate(vertices)};a=[[0]*len(vertices) for v in vertices]
    for u,v in edges:
        if u==v:continue
        i,j=index[u],index[v];a[i][i]+=1;a[j][j]+=1;a[i][j]-=1;a[j][i]-=1
    return determinant([r[:-1] for r in a[:-1]])


def twist(pd,crossing,axis,length):
    assert length>=1 and length%2==1;ports=pd[crossing];next_edge=max(e for row in pd for e in row)+1
    left=[ports[axis]]+list(range(next_edge,next_edge+length-1))+[ports[(axis+3)%4]]
    next_edge+=length-1
    right=[ports[(axis+1)%4]]+list(range(next_edge,next_edge+length-1))+[ports[(axis+2)%4]]
    pieces=[]
    for j in range(length):
        pieces.append([left[j],right[j],right[j+1],left[j+1]] if axis==0 else [right[j],right[j+1],left[j+1],left[j]])
    return pd[:crossing]+pieces+pd[crossing+1:]


def main():
    plan=json.loads((HERE/'periodic_twist_family_plan.json').read_text())
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    old={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    full={r['name']:r for r in json.loads((HERE/'finite_cover_full_results.json').read_text())['cases']}
    results=[]
    for name in plan['targets']:
        pd=panel[name]['pd'];cr=oriented_diagram(pd);images,witness=induced(old[name],full[name],cr)
        arc_of={e:o for row,(o,a,b,_) in zip(pd,cr) for e in [row[1],row[3]]}
        assert len(arc_of)==26
        candidates=[]
        for c in range(13):
            groups,outside,boundary=outside_arcs(pd,c)
            for axis in [0,1]:
                for seeds in combinations(range(len(groups)),4):
                    trace=propagate(len(groups),outside,boundary,axis,seeds)
                    if trace is not None:
                        pair=[images[arc_of[pd[c][h]]] for h in [axis,axis+1]]
                        m=order(compose(*pair));candidates.append((m,c,axis,list(seeds),trace,groups,outside,boundary));break
        assert candidates
        m,c,axis,seeds,trace,groups,outside,boundary=min(candidates)
        vertices,edges=checkerboard(pd);total=trees(vertices,edges);u,v=edges[c]
        deletion=trees(vertices,edges[:c]+edges[c+1:])
        contraction=trees([w for w in vertices if w!=v],[[u if a==v else a,u if b==v else b] for i,(a,b) in enumerate(edges) if i!=c])
        assert total==deletion+contraction and deletion>0 and contraction>0
        sample=[]
        for N in range(3):
            length=1+2*m*N;new_pd=twist(pd,c,axis,length);vs,es=checkerboard(new_pd);sample.append({'N':N,'twist_length':length,'crossings':len(new_pd),'determinant':trees(vs,es),'pd':new_pd})
        slope=(sample[1]['determinant']-total)//(2*m)
        assert slope in [deletion,contraction] and all(r['determinant']==total+2*m*slope*r['N'] for r in sample)
        results.append({'name':name,'knot_meridian_images':images,'stabilizer_noncommutator':witness,
                        'selected_crossing_zero_based':c,'source_crossing_ports':pd[c],'axis':axis,'product_order':m,
                        'twist_length_formula':f'1+{2*m}*N','four_outside_seeds':seeds,'outside_arc_edges':groups,
                        'outside_crossings':outside,'boundary_arc_indices':boundary,'uniform_upper_trace':trace,
                        'successful_local_templates':len(candidates),'checkerboard_vertices':vertices,'checkerboard_edges':edges,
                        'deleted_edge_tree_count':deletion,'contracted_edge_tree_count':contraction,
                        'determinant_constant':total,'determinant_slope_in_twist_length':slope,
                        'determinant_slope_in_N':2*m*slope,'samples':sample})
    paths=['periodic_twist_family_plan.json','audit_periodic_twist_families.py','joint_a5_panel.json','three_sheet_nonabelian_results.json','finite_cover_full_results.json','finite_cover_quotient_verification.json','parabolic_lifts.py']
    report={'status':'CONSTRUCTED periodic sharp-four family inputs','families':results,
            'scope':'Finite input construction; uniform proofs and independent replay required before a theorem claim. Distinctness uses the determinant, not approximate recognition. No worldwide novelty or family exclusion inferred.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'periodic_twist_family_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'families':[{k:r[k] for k in ['name','selected_crossing_zero_based','axis','product_order','successful_local_templates','determinant_constant','determinant_slope_in_N']} for r in results]},indent=2))


if __name__=='__main__':main()

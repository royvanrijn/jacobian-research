#!/usr/bin/env python3
"""Independent pointwise and literal-diagram replay of the ternary control.

No producer matrix arithmetic, root equations, group enumeration or diagram
joins are imported. The old parser supplies only legacy arc names, checked
against independent geometric paths and orientation.
"""

from collections import defaultdict
from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram
from verify_double_transpositions import geometric_crossings

HERE=Path(__file__).resolve().parent
POINTS=list(product(range(3),repeat=3))
POINT_INDEX={v:i for i,v in enumerate(POINTS)}


def act(matrix,point):
    n=len(point)
    return tuple(sum(matrix[n*i+j]*point[j] for j in range(n))%3 for i in range(n))


def affine(matrix,translation,point):
    return tuple((x+y)%3 for x,y in zip(act(matrix,point),translation))


def mat_rank(rows,n):
    a=[[v%3 for v in row] for row in rows];pivot=0
    for c in range(n):
        k=next((i for i in range(pivot,len(a)) if a[i][c]),None)
        if k is None:continue
        a[pivot],a[k]=a[k],a[pivot]
        scale=1 if a[pivot][c]==1 else 2
        a[pivot]=[v*scale%3 for v in a[pivot]]
        for i in range(len(a)):
            if i!=pivot:
                value=a[i][c];a[i]=[(x-value*y)%3 for x,y in zip(a[i],a[pivot])]
        pivot+=1
        if pivot==len(a):break
    return pivot


def inverse_matrix(a):
    n=3 if len(a)==9 else 4
    # Recover inverse columns by exhaustive vectors, not cofactor/Gaussian code.
    points=list(product(range(3),repeat=n));preimages={act(a,v):v for v in points}
    assert len(preimages)==3**n
    columns=[preimages[tuple(int(i==j) for i in range(n))] for j in range(n)]
    return tuple(columns[j][i] for i in range(n) for j in range(n))


def signed_crossings(pd):
    triples=geometric_crossings(pd)
    ends=defaultdict(list)
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends[e].append((c,h))
    incoming=set();current=(0,0)
    for _ in range(2*len(pd)):
        assert current not in incoming;incoming.add(current)
        c,h=current;outgoing=(c,(h+2)%4);pair=ends[pd[outgoing[0]][outgoing[1]]]
        current=pair[1] if pair[0]==outgoing else pair[0]
    assert current==(0,0) and all((c,0) in incoming for c in range(len(pd)))
    return [(*triple,int((c,1) in incoming)) for c,triple in enumerate(triples)]


def word_for_crossing(row):
    o,a,b,positive=row
    return [-(o+1),a+1,o+1,-(b+1)] if positive else [o+1,a+1,-(o+1),-(b+1)]


def evaluate_word(matrices,word,point):
    inverses={abs(i)-1:inverse_matrix(matrices[abs(i)-1]) for i in word if i<0}
    for letter in reversed(word):
        index=abs(letter)-1
        point=act(matrices[index] if letter>0 else inverses[index],point)
    return point


def root_from_points(a):
    differences={tuple((x-y)%3 for x,y in zip(act(a,v),v)) for v in POINTS}
    assert len(differences)==3
    return min(v for v in differences if any(v))


def pointwise_root_matrix(crossings,matrices):
    n=len(matrices);roots=[root_from_points(a) for a in matrices]
    inverses=[inverse_matrix(a) for a in matrices];rows=[]
    for crossing in crossings:
        word=word_for_crossing(crossing);block=[[],[],[]]
        for variable in range(n):
            point=(0,0,0)
            for letter in reversed(word):
                i=abs(letter)-1
                if letter>0:
                    point=affine(matrices[i],roots[i] if i==variable else (0,0,0),point)
                else:
                    # Inverse translation is computed from the actual inverse.
                    translation=act(inverses[i],roots[i]) if i==variable else (0,0,0)
                    point=affine(inverses[i],tuple(-v%3 for v in translation),point)
            for j in range(3):block[j].append(point[j])
        rows.extend(block)
    return rows,roots


def check_control(report):
    matrices=report['chain_matrices'];points4=list(product(range(3),repeat=4))
    assert len(matrices)==4
    for a in matrices:
        assert len({tuple((x-y)%3 for x,y in zip(act(a,v),v)) for v in points4})==3
        assert all(act(a,act(a,act(a,v)))==v for v in points4)
    assert [v for v in points4 if all(act(a,v)==v for a in matrices)]==[(0,0,0,0)]
    for i in range(3):
        word=[i+1,i+2,i+1,-(i+2),-(i+1),-(i+2)]
        assert all(evaluate_word(matrices,word,v)==v for v in points4)
    # Restore all three literal source trefoils; only their outgoing edge1
    # ports have been cycled by the two connected sums.
    trefoil=[[1,4,2,5],[3,6,4,1],[5,2,6,3]]
    unnumbered=[[e+6*i for e in row] for i in range(3) for row in trefoil]
    for i in range(3):unnumbered[3*i+1][3]=6*((i+1)%3)+1
    names={int(k):v for k,v in report['edge_renumbering'].items()}
    assert set(names)==set(names.values())==set(range(1,19))
    pd=report['pd'];assert [[names[e] for e in row] for row in unnumbered]==pd
    crossings=signed_crossings(pd)
    assert crossings==[tuple(row) for row in oriented_diagram(pd)]
    full=report['arc_matrices'];assert len(full)==9
    for crossing in crossings:
        assert all(evaluate_word(full,word_for_crossing(crossing),v)==v for v in points4)
    seeds=report['seeds'];seen=set(seeds);assert len(seen)==4
    valid={(o,a,b) for o,a,b,_ in crossings}|{(o,b,a) for o,a,b,_ in crossings}
    for o,a,b in report['seed_trace']:
        assert (o,a,b) in valid and o in seen and a in seen and b not in seen;seen.add(b)
    assert seen==set(range(9))
    linear=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in full]
    rows,roots=pointwise_root_matrix(crossings,linear)
    assert rows==report['root_matrix'] and 9-mat_rank(rows,9)==report['root_dimension']==4
    basis=report['root_basis'];assert len(basis)==mat_rank(basis,9)==4
    assert all(sum(x*y for x,y in zip(row,v))%3==0 for row in rows for v in basis)
    coefficients=report['arc_root_coefficients']
    assert all(tuple(m[4*i+3] for i in range(3))==tuple(c*x%3 for x in r) for m,r,c in zip(full,roots,coefficients))
    # Independent permutation group model on all27 vectors, multiplying on
    # the opposite side from the producer's dense-matrix closure.
    first_three=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in matrices[:3]]
    generators=[tuple(POINT_INDEX[act(a,v)] for v in POINTS) for a in first_three]
    group={tuple(range(27))};queue=list(group)
    for p in queue:
        for g in generators:
            q=tuple(g[p[i]] for i in range(27))
            if q not in group:assert len(group)<5616;group.add(q);queue.append(q)
    assert len(group)==report['linear_first_three_image_order']==5616
    return {'homogeneous_points_checked':81,'chain_braid_relators':3,'literal_PD_crossings':9,
        'root_dimension':4,'linear_image_order':5616,'meridional_rank':4,'bridge_number':4}


def main():
    report=json.loads((HERE/'ternary_transvection_control.json').read_text())
    for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    result=check_control(report);rejects=[]
    for kind in ['false_functional_rank','false_crossing_image','lost_kernel_vector']:
        bad=deepcopy(report)
        if kind=='false_functional_rank':bad['chain_matrices'][3]=bad['chain_matrices'][2]
        elif kind=='false_crossing_image':bad['arc_matrices'][0][0]=(bad['arc_matrices'][0][0]+1)%3
        else:bad['root_basis']=bad['root_basis'][:-1]
        try:check_control(bad)
        except (AssertionError,KeyError,StopIteration):rejects.append(kind)
        else:raise AssertionError('accepted corrupted control')
    files=['ternary_transvection_control.json','verify_ternary_transvection_control.py','verify_double_transpositions.py','parabolic_lifts.py']
    receipt={'status':'PASS independent ternary triple-trefoil control','checks':result,'rejected_corruptions':rejects,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Known composite knot; pointwise matrices, direct affine word coefficients, literal trefoil-source restoration and geometric orientation. Source parser only supplies checked legacy arc names. Full affine order uses the separate linear-plus-translation proof, not a direct151632-element enumeration.'}
    (HERE/'ternary_transvection_control_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','checks','rejected_corruptions']},indent=2))


if __name__=='__main__':main()

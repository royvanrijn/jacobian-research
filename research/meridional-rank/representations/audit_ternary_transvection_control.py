#!/usr/bin/env python3
"""Construct the frozen known triple-trefoil order3 affine control."""

import hashlib
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram
from search_lift_panel import seed_certificate
from ternary_transvections import (identity,image,linear_group,lift_matrix,
    model,multiply,power,rank)

HERE=Path(__file__).resolve().parent
TREFOIL=[[1,4,2,5],[3,6,4,1],[5,2,6,3]]


def invert4(a):
    # Gaussian inverse, used only for the explicit rebase words.
    rows=[list(a[4*i:4*i+4])+[int(i==j) for j in range(4)] for i in range(4)]
    for c in range(4):
        k=next(i for i in range(c,4) if rows[i][c]);rows[c],rows[k]=rows[k],rows[c]
        inv=pow(rows[c][c],-1,3);rows[c]=[v*inv%3 for v in rows[c]]
        for i in range(4):
            if i!=c:
                f=rows[i][c];rows[i]=[(x-f*y)%3 for x,y in zip(rows[i],rows[c])]
    return tuple(v for row in rows for v in row[4:])


def source(a,b):
    third=multiply(multiply(power(a,2),b),a)
    arc_values=[a,b,third];edges={}
    for row,(o,u,v,_) in zip(TREFOIL,oriented_diagram(TREFOIL)):
        for e,g in [(row[0],u),(row[1],o),(row[2],v),(row[3],o)]:
            if e in edges:assert edges[e]==arc_values[g]
            edges[e]=arc_values[g]
    return [r[:] for r in TREFOIL],edges


def join(pd,edges,other,other_edges):
    offset=max(edges);right=[[e+offset for e in row] for row in other];left=[r[:] for r in pd]
    assert edges[1]==other_edges[1]
    l=[(i,h) for i,row in enumerate(left) for h,e in enumerate(row) if e==1]
    r=[(i,h) for i,row in enumerate(right) for h,e in enumerate(row) if e==offset+1]
    i,h=l[1];left[i][h]=offset+1;i,h=r[1];right[i][h]=1
    return left+right,{**edges,**{e+offset:p for e,p in other_edges.items()}}


def renumber(pd,edges):
    ends={}
    for c,row in enumerate(pd):
        for h,e in enumerate(row):ends.setdefault(e,[]).append((c,h))
    current=(0,0);start=current;names={}
    while True:
        c,h=current;e=pd[c][h];assert e not in names;names[e]=len(names)+1
        outgoing=(c,(h+2)%4);other_edge=pd[outgoing[0]][outgoing[1]]
        pair=ends[other_edge];current=pair[1] if pair[0]==outgoing else pair[0]
        if current==start:break
    assert len(names)==2*len(pd)
    return [[names[e] for e in row] for row in pd],{names[e]:p for e,p in edges.items()},names


def main():
    started=time.monotonic();plan=json.loads((HERE/'ternary_transvection_control_plan.json').read_text())
    matrices=[]
    for v,f in zip(plan['directions'],plan['homogeneous_functionals']):
        u=v+[0];assert sum(a*b for a,b in zip(u,f))%3==0
        a=tuple((int(i==j)+u[i]*f[j])%3 for i in range(4) for j in range(4))
        assert power(a,3)==identity(4)
        assert rank([[(a[4*i+j]-int(i==j))%3 for j in range(4)] for i in range(4)],4)==1
        matrices.append(a)
    assert rank(plan['homogeneous_functionals'],4)==4
    for a,b in zip(matrices,matrices[1:]):assert multiply(multiply(a,b),a)==multiply(multiply(b,a),b)
    a,b,c,d=matrices
    g=multiply(a,b);h=invert4(g)
    k=invert4(multiply(multiply(b,c),g))
    partners=[b,multiply(multiply(h,c),g),multiply(multiply(k,d),invert4(k))]
    assert multiply(multiply(h,b),g)==a
    assert multiply(multiply(k,c),invert4(k))==a
    pd,edges=source(a,partners[0])
    for partner in partners[1:]:
        other,other_edges=source(a,partner);pd,edges=join(pd,edges,other,other_edges)
    pd,edges,renumbering=renumber(pd,edges);crossings=oriented_diagram(pd)
    full=[]
    for row,(o,u,v,_) in zip(pd,crossings):
        while len(full)<=max(o,u,v):full.append(None)
        for edge,index in [(row[0],u),(row[1],o),(row[2],v),(row[3],o)]:
            if full[index] is not None:assert full[index]==edges[edge]
            full[index]=edges[edge]
    for o,u,v,sign in crossings:
        left=multiply(multiply(power(full[o],2),full[u]),full[o]) if sign else multiply(multiply(full[o],full[u]),power(full[o],2))
        assert left==full[v]
    elements,_,_,centralizer,_=model();indices={m:i for i,m in enumerate(elements)}
    labels=[indices[tuple(m[4*i+j] for i in range(3) for j in range(3))] for m in full]
    matrix,basis,roots=lift_matrix(crossings,labels,elements);assert len(basis)==4
    coefficients=[]
    for m,r in zip(full,roots):
        position=next(i for i,v in enumerate(r) if v);coefficient=m[4*position+3]
        assert all(m[4*i+3]==coefficient*r[i]%3 for i in range(3));coefficients.append(coefficient)
    assert all(sum(x*y for x,y in zip(row,coefficients))%3==0 for row in matrix)
    seeds,trace,trials=seed_certificate(9,[r[:3] for r in crossings],4,126,time.monotonic()+20)
    assert seeds is not None
    first_three_linear=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in matrices[:3]]
    group=linear_group(first_three_linear)
    assert len(group)==5616
    files=['ternary_transvection_control_plan.json','audit_ternary_transvection_control.py','ternary_transvections.py','parabolic_lifts.py','search_lift_panel.py','verify_integer_lifts.py']
    report={'status':'CONSTRUCTED known triple-trefoil order3 control','directions':plan['directions'],
        'functionals':plan['homogeneous_functionals'],'chain_matrices':matrices,
        'chain_braid_relators_checked':3,'functional_rank':4,'linear_first_three_image_order':len(group),
        'linear_group':'SL3(F3)','affine_group_order_from_linear_and_translation_proof':27*len(group),
        'rebased_partners':partners,'rebase_conjugators':[identity(4),h,k],
        'pd':pd,'edge_renumbering':renumbering,'arc_matrices':full,'arc_linear_labels':labels,'arc_root_coefficients':coefficients,
        'root_matrix':matrix,'root_basis':basis,'root_dimension':len(basis),
        'seeds':seeds,'seed_trace':trace,'seed_trials':trials,'meridional_rank':4,'bridge_number':4,
        'transvection_class_size':len(elements),'centralizer_size':len(centralizer),
        'seconds':round(time.monotonic()-started,6),
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Known connected sum of three trefoils, no new knot closure. Full rank-one common-fixed-space and root-nullity witnesses outside order2 bases. Target search not executed.'}
    (HERE/'ternary_transvection_control.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','root_dimension','functional_rank','linear_first_three_image_order','transvection_class_size','centralizer_size','seconds']},indent=2))


if __name__=='__main__':main()

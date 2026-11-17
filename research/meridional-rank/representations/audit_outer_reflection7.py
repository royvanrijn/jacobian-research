#!/usr/bin/env python3
"""Use the outer deck action to construct ternary reflection bases over F7."""
import hashlib
from itertools import product
import json
from pathlib import Path
from gl3f2 import group_model
from parabolic_lifts import oriented_diagram
from verify_integer_lifts import nullspace

HERE=Path(__file__).resolve().parent
P=7


def mm(a,b):return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%7 for i in range(2) for j in range(2))
def det(a):return (a[0]*a[3]-a[1]*a[2])%7
def norm(a):
    s=pow(next(x for x in a if x),-1,7)
    return tuple(x*s%7 for x in a)
def inverse(a):return norm((a[3],-a[1]%7,-a[2]%7,a[0]))


def projective_model():
    elements=sorted({norm(a) for a in product(range(7),repeat=4) if det(a)})
    identity=(1,0,0,1);elements.remove(identity);elements.insert(0,identity)
    assert len(elements)==336
    index={a:i for i,a in enumerate(elements)}
    table=[[index[norm(mm(a,b))] for b in elements] for a in elements]
    inverses=[index[inverse(a)] for a in elements]
    even=[i for i,a in enumerate(elements) if pow(det(a),3,7)==1];assert len(even)==168
    orders=[]
    for i in range(336):
        x=i;n=1
        while x:x=table[x][i];n+=1
        orders.append(n)
    return elements,table,inverses,orders,even


def identify(oldtable,stdpair,table,orders,even):
    a,b=stdpair;mapping=None
    for u in even:
        if orders[u]!=2:continue
        for v in even:
            if orders[v]!=3 or orders[table[u][v]]!=7:continue
            trial={0:0};queue=[0];failed=False
            for x in queue:
                for g,h in [(a,u),(b,v)]:
                    y,z=oldtable[x][g],table[trial[x]][h]
                    if y in trial:
                        if trial[y]!=z:failed=True;break
                    else:trial[y]=z;queue.append(y)
                if failed:break
            if not failed and len(trial)==len(set(trial.values()))==168:
                mapping=[trial[i] for i in range(168)];break
        if mapping is not None:break
    assert mapping is not None
    assert all(mapping[oldtable[i][j]]==table[mapping[i]][mapping[j]] for i in range(168) for j in range(168))
    return mapping


def adjoint_twist(a):
    x,y,z,w=a;scale=pow(det(a),3,7)*pow(det(a),-1,7)%7
    return tuple(c*scale%7 for c in [x*w+y*z,-x*z,y*w,-2*x*y,x*x,-y*y,2*z*w,-z*z,w*w])
def mult3(a,b):return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3))%7 for i in range(3) for j in range(3))
def vec3(a,v):return tuple(sum(a[3*i+j]*v[j] for j in range(3))%7 for i in range(3))
I=(1,0,0,0,1,0,0,0,1)


def roots(images,words):
    rs=[]
    for a in images:
        d=[(x-I[i])%7 for i,x in enumerate(a)]
        assert len(nullspace([d[3*i:3*i+3] for i in range(3)],3,7))==2
        r=next(tuple(d[3*i+j] for i in range(3)) for j in range(3) if any(d[3*i+j] for i in range(3)))
        rs.append(r)
    rows=[]
    for word in words:
        prefix=I;block=[[0]*13 for _ in range(3)]
        for l in word:
            j=abs(l)-1;v=vec3(prefix,rs[j])
            for i in range(3):block[i][j]=(block[i][j]+v[i])%7
            prefix=mult3(prefix,images[j])
        assert prefix==I;rows.extend(block)
    basis=nullspace(rows,13,7)
    return {'roots':rs,'matrix':rows,'kernel_basis':basis,'dimension':len(basis),'sharp_four':len(basis)>=4}


def main():
    plan=json.loads((HERE/'outer_reflection7_plan.json').read_text());assert plan['frozen_before_execution']
    source=json.loads((HERE/'gl3f2_results.json').read_text());assert source['status']=='COMPLETE'
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases'][:16]}
    _,oldtable,_,_,_,_,_,_,std=group_model()
    elements,table,inverses,orders,even=projective_model()
    identification=identify(oldtable,std,table,orders,even)
    reps=[adjoint_twist(a) for a in elements]
    assert len(set(reps))==336 and reps[0]==I
    assert all(mult3(reps[i],reps[j])==reps[table[i][j]] for i in range(336) for j in range(336))
    results=[]
    for row in source['results']:
        cr=oriented_diagram(panel[row['name']]['pd'])
        words=[[o+1,a+1,o+1,b+1] for o,a,b,_ in cr]
        for mapping in row['maps']:
            arcs=[identification[i] for i in mapping['arc_images']]
            deck=[t for t in range(336) if orders[t]==2 and all(table[table[t][g]][inverses[t]]==inverses[g] for g in arcs)]
            assert len(deck)==1 and deck[0] not in even
            meridians=[table[g][deck[0]] for g in arcs]
            assert all(orders[m]==2 and m not in even for m in meridians)
            images=[reps[m] for m in meridians]
            root=roots(images,words)
            result={'name':row['name'],'orbit_index':mapping['orbit_index'],'deck_projective_index':deck[0],
                    'projective_meridians':meridians,'meridian_matrices':images,**root}
            results.append(result)
            print(json.dumps({'name':row['name'],'orbit_index':mapping['orbit_index'],'dimension':root['dimension'],'sharp_four':root['sharp_four']}),flush=True)
    assert len(results)==10
    paths=['audit_outer_reflection7.py','outer_reflection7_plan.json','gl3f2.py','gl3f2_results.json','joint_a5_panel.json','parabolic_lifts.py','verify_integer_lifts.py']
    report={'status':'COMPLETE','projective_elements':elements,'GL_to_PSL_identification':identification,
            'results':results,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'outer_reflection7_results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

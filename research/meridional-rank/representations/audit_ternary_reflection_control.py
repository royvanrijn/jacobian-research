#!/usr/bin/env python3
"""Known order2 ternary affine control on the retained triple-trefoil PD."""

import hashlib
import json
from pathlib import Path

from audit_ternary_transvection_control import TREFOIL,invert4,join,renumber
from parabolic_lifts import oriented_diagram
from ternary_reflections import lift_matrix,model
from ternary_transvections import identity,linear_group,multiply,power,rank

HERE=Path(__file__).resolve().parent


def source(a,b):
    third=multiply(multiply(a,b),a);values=[a,b,third];edges={}
    for row,(o,u,v,_) in zip(TREFOIL,oriented_diagram(TREFOIL)):
        for edge,i in [(row[0],u),(row[1],o),(row[2],v),(row[3],o)]:
            if edge in edges:assert edges[edge]==values[i]
            edges[edge]=values[i]
    return [r[:] for r in TREFOIL],edges


def main():
    plan=json.loads((HERE/'ternary_reflection_plan.json').read_text());control=plan['control']
    matrices=[]
    for v,f in zip(control['directions'],control['functionals']):
        u=v+[0];assert sum(a*b for a,b in zip(u,f))%3==1
        m=tuple((int(i==j)+u[i]*f[j])%3 for i in range(4) for j in range(4))
        assert power(m,2)==identity(4);matrices.append(m)
    assert rank(control['functionals'],4)==4
    for a,b in zip(matrices,matrices[1:]):assert multiply(multiply(a,b),a)==multiply(multiply(b,a),b)
    a,b,c,d=matrices;g=multiply(a,b);h=invert4(g);k=invert4(multiply(multiply(b,c),g))
    partners=[b,multiply(multiply(h,c),g),multiply(multiply(k,d),invert4(k))]
    assert multiply(multiply(h,b),g)==a and multiply(multiply(k,c),invert4(k))==a
    pd,edges=source(a,partners[0])
    for p in partners[1:]:other,e=source(a,p);pd,edges=join(pd,edges,other,e)
    pd,edges,names=renumber(pd,edges);cr=oriented_diagram(pd)
    old=json.loads((HERE/'ternary_transvection_control.json').read_text())
    assert pd==old['pd'] and {int(k):v for k,v in old['edge_renumbering'].items()}==names
    full=[None]*9
    for row,(o,u,v,_) in zip(pd,cr):
        for edge,i in [(row[0],u),(row[1],o),(row[2],v),(row[3],o)]:
            if full[i] is not None:assert full[i]==edges[edge]
            full[i]=edges[edge]
    for o,u,v,_ in cr:assert multiply(multiply(full[o],full[u]),full[o])==full[v]
    elements,_,centralizer,_=model();indices={m:i for i,m in enumerate(elements)}
    linear=[tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in full]
    labels=[indices[m] for m in linear];matrix,basis,roots=lift_matrix(cr,labels,elements);assert len(basis)==4
    coefficients=[]
    for m,r in zip(full,roots):
        i=next(j for j,v in enumerate(r) if v);c=m[4*i+3]
        assert all(m[4*j+3]==c*r[j]%3 for j in range(3));coefficients.append(c)
    group=linear_group([tuple(m[4*i+j] for i in range(3) for j in range(3)) for m in matrices[:3]])
    assert len(group)==24
    linear_fourth=tuple(matrices[3][4*i+j] for i in range(3) for j in range(3));assert linear_fourth in group
    files=['ternary_reflection_plan.json','audit_ternary_reflection_control.py','ternary_reflections.py',
        'ternary_transvections.py','audit_ternary_transvection_control.py','ternary_transvection_control.json','parabolic_lifts.py']
    report={'status':'CONSTRUCTED known triple-trefoil ternary reflection control','directions':control['directions'],
        'functionals':control['functionals'],'chain_matrices':matrices,'functional_rank':4,
        'linear_first_three_image_order':24,'affine_group_order_from_linear_and_translation_proof':648,
        'rebased_partners':partners,'rebase_conjugators':[identity(4),h,k],
        'pd':pd,'edge_renumbering':names,'arc_matrices':full,'arc_linear_labels':labels,
        'arc_root_coefficients':coefficients,'root_matrix':matrix,'root_basis':basis,'root_dimension':4,
        'seeds':old['seeds'],'seed_trace':old['seed_trace'],'meridional_rank':4,'bridge_number':4,
        'reflection_class_size':117,'centralizer_size':len(centralizer),
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Known composite control; independent replay required. Affine order uses the separately written24-element linear-copy and translation-orbit span proof.'}
    (HERE/'ternary_reflection_control.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','functional_rank','root_dimension','linear_first_three_image_order','reflection_class_size','centralizer_size']},indent=2))


if __name__=='__main__':main()

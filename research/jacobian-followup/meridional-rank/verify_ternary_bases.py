#!/usr/bin/env python3
"""Independent full ternary3 rank-one base and root-space verification."""

from collections import Counter
from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path
import time

from verify_ternary_transvection_control import (POINTS,POINT_INDEX,act,mat_rank,
    pointwise_root_matrix,root_from_points,signed_crossings)

HERE=Path(__file__).resolve().parent


def permutation(a):return tuple(POINT_INDEX[act(a,v)] for v in POINTS)


def inverse_permutation(p):
    q=[0]*len(p)
    for i,j in enumerate(p):q[j]=i
    return tuple(q)


def catalogue_data(catalogue,order):
    elements=[tuple(m) for m in catalogue['elements']];identity=tuple(range(27))
    # Filter all3^9 matrices by independent row reduction and rank-one defect,
    # rather than producing canonical direction-functional pairs.
    gl=[];expected=set()
    for a in product(range(3),repeat=9):
        rows=[a[3*i:3*i+3] for i in range(3)]
        if mat_rank(rows,3)!=3:continue
        gl.append(a)
        difference=[[(a[3*i+j]-int(i==j))%3 for j in range(3)] for i in range(3)]
        if mat_rank(difference,3)!=1:continue
        p=permutation(a);power=identity
        for _ in range(order):power=tuple(p[power[i]] for i in range(27))
        if power==identity:expected.add(a)
    assert len(gl)==11232 and len(elements)==len(set(elements))==len(expected)
    assert set(elements)==expected
    perms=[permutation(a) for a in elements];indices={p:i for i,p in enumerate(perms)}
    centralizer=[];actions=[];tau=perms[0]
    for g in gl:
        p=permutation(g)
        if any(p[tau[i]]!=tau[p[i]] for i in range(27)):continue
        centralizer.append(list(g));pinv=inverse_permutation(p)
        actions.append([indices[tuple(p[q[pinv[i]]] for i in range(27))] for q in perms])
    assert centralizer==catalogue['centralizer']
    size=len(elements);marked=bytearray(size**3);weights=0
    for row in catalogue['representatives']:
        a,b,c,d=row['tuple'];assert a==0 and all(0<=j<size for j in [b,c,d])
        orbit={(g[b]*size+g[c])*size+g[d] for g in actions}
        assert len(orbit)==row['raw_weight'] and len(orbit)*row['stabilizer']==len(actions)
        for encoded in orbit:assert not marked[encoded];marked[encoded]=1
        weights+=len(orbit)
    assert all(marked) and weights==size**3
    assert weights==catalogue['counts']['raw_normalized_tuples']
    assert len(catalogue['representatives'])==catalogue['counts']['four_tuple_orbits']
    assert {int(k):v for k,v in catalogue['counts']['stabilizer_histogram'].items()}==dict(Counter(r['stabilizer'] for r in catalogue['representatives']))
    forward=[];backward=[]
    for p in perms:
        inv=inverse_permutation(p)
        forward.append([indices[tuple(inv[q[p[i]]] for i in range(27))] for q in perms])
        backward.append([indices[tuple(p[q[inv[i]]] for i in range(27))] for q in perms])
    return elements,perms,actions,forward,backward


def trace_steps(crossings,seeds,trace):
    seen=set(seeds);assert len(seen)==4;steps=[]
    for o,source,target in trace:
        match=next((r for r in crossings if r[0]==o and {r[1],r[2]}=={source,target}),None)
        assert match is not None and o in seen and source in seen and target not in seen
        _,a,b,sign=match
        use_forward=bool(sign) if source==a else not bool(sign)
        steps.append((o,source,target,use_forward));seen.add(target)
    assert seen==set(range(len(crossings)))
    return steps


def labels_from_trace(assignment,seeds,steps,n,forward,backward):
    labels=[-1]*n
    for seed,label in zip(seeds,assignment):labels[seed]=label
    for o,a,b,direct in steps:labels[b]=(forward if direct else backward)[labels[o]][labels[a]]
    assert -1 not in labels
    return labels


def group_order(generators,cap=11232):
    reached={tuple(range(27))};queue=list(reached)
    for p in queue:
        for q in generators:
            value=tuple(q[p[i]] for i in range(27))
            if value not in reached:assert len(reached)<cap;reached.add(value);queue.append(value)
    return len(reached)


def check_kernel(base,crossings,elements):
    matrices=[elements[i] for i in base['arc_labels']]
    rows,roots=pointwise_root_matrix(crossings,matrices);n=len(matrices)
    assert rows==base['root_matrix'] and [list(r) for r in roots]==base['roots']
    dimension=n-mat_rank(rows,n);basis=base['root_basis']
    assert dimension==base['root_dimension']==len(basis)==mat_rank(basis,n)
    assert all(sum(a*b for a,b in zip(row,v))%3==0 for row in rows for v in basis)
    return dimension


def main():
    started=time.monotonic();panel=json.loads((HERE/'parabolic_panel.json').read_text())['cases']
    receipts=[];rejects=[]
    for prefix,order,expected_orbits in [('ternary_transvection',3,21632),('ternary_reflection',2,34084)]:
        catalogue=json.loads((HERE/(prefix+'_catalogue.json')).read_text())
        report=json.loads((HERE/(prefix+'_results.json')).read_text())
        assert report['status']=='COMPLETE' and report['sharp_names']==[]
        for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
        assert hashlib.sha256((HERE/(prefix+'_catalogue.json')).read_bytes()).hexdigest()==report['catalogue_sha256']
        elements,perms,actions,forward,backward=catalogue_data(catalogue,order)
        assert len(catalogue['representatives'])==expected_orbits
        control=json.loads((HERE/(prefix+'_control.json')).read_text())
        end=report['end_to_end_control'];crossings=signed_crossings(control['pd'])
        steps=trace_steps(crossings,control['seeds'],control['seed_trace'])
        labels=labels_from_trace(end['canonical_seed_tuple'],control['seeds'],steps,9,forward,backward)
        assert labels==end['arc_labels']
        assert all((forward if s else backward)[labels[o]][labels[a]]==labels[b] for o,a,b,s in crossings)
        rows,_=pointwise_root_matrix(crossings,[elements[i] for i in labels])
        assert 9-mat_rank(rows,9)==end['root_dimension']==4
        assert [r['name'] for r in report['results']]==[c['name'] for c in panel]
        results=[]
        for result,case in zip(report['results'],panel):
            crossings=signed_crossings(case['pd']);n=len(crossings)
            steps=trace_steps(crossings,result['seeds'],result['seed_trace'])
            actual={}
            for row in catalogue['representatives']:
                assignment=tuple(row['tuple'])
                labels=labels_from_trace(assignment,result['seeds'],steps,n,forward,backward)
                if all((forward if s else backward)[labels[o]][labels[a]]==labels[b] for o,a,b,s in crossings):actual[assignment]=labels
            stored={tuple(b['seed_assignment']):b for b in result['bases']}
            assert set(stored)==set(actual) and result['seed_orbits_tested']==expected_orbits
            verified=[]
            for assignment,labels in actual.items():
                base=stored[assignment];assert labels==base['arc_labels']
                dimension=check_kernel(base,crossings,elements)
                image_order=group_order([perms[labels[s]] for s in result['seeds']])
                principal=[]
                for label in labels:
                    a=elements[label];r=root_from_points(a);i=next(j for j,v in enumerate(r) if v)
                    principal.append([(int(i==j)-a[3*i+j])%3 for j in range(3)])
                principal_dimension=mat_rank(principal,3)
                assert principal_dimension==dimension
                assert all(sum(row[i]*principal[i][j] for i in range(n))%3==0 for row in base['root_matrix'] for j in range(3))
                verified.append({'seed_assignment':assignment,'root_dimension':dimension,'knot_image_order':image_order,
                    'constant_base':len(set(labels))==1,'principal_dimension':principal_dimension,
                    'all_admissible_affine_lifts_principal':True})
            results.append({'name':case['name'],'base_orbits':verified,'root_dimension_histogram':dict(Counter(b['root_dimension'] for b in verified))})
        # Alter an actual kernel and an actual canonical orbit, without relying
        # on input-hash rejection to detect either mathematical corruption.
        base=deepcopy(report['results'][0]['bases'][0]);base['root_basis'][0]=[0]*13
        try:check_kernel(base,signed_crossings(panel[0]['pd']),elements)
        except AssertionError:rejects.append(prefix+':false_kernel_basis')
        else:raise AssertionError('accepted false kernel basis')
        bad=deepcopy(catalogue);bad['representatives'][1]=deepcopy(bad['representatives'][0])
        try:catalogue_data(bad,order)
        except AssertionError:rejects.append(prefix+':overlapping_tuple_orbits')
        else:raise AssertionError('accepted overlapping tuple orbits')
        receipts.append({'meridian_order':order,'class_size':len(elements),'tuple_orbits':expected_orbits,
            'raw_tuples':len(elements)**3,'control_root_dimension':4,'targets':results})
        print(json.dumps({'order':order,'status':'PASS full ternary class/orbit/geometry/pointwise-root replay',
            'image_order_histogram':dict(Counter(b['knot_image_order'] for r in results for b in r['base_orbits']))}),flush=True)
    files=['verify_ternary_bases.py','verify_ternary_transvection_control.py','verify_double_transpositions.py',
        'ternary_transvection_results.json','ternary_transvection_catalogue.json',
        'ternary_reflection_results.json','ternary_reflection_catalogue.json','parabolic_panel.json']
    receipt={'status':'PASS independent complete ternary3 rank-one base gate','classes':receipts,
        'rejected_corruptions':rejects,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Exactly16 diagrams and all rank-one GL3(F3) meridian bases, both possible nonidentity eigenvalues. Full centralizer raw orbit marking; separate geometric signed trace propagation, affine point kernel and permutation image enumeration. No all-field or higher-dimensional exclusion; no target closure. Representation orbits are not asserted to be distinct quotient kernels.'}
    (HERE/'ternary_base_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'rejected_corruptions':rejects,'seconds':round(time.monotonic()-started,6)},indent=2))


if __name__=='__main__':main()

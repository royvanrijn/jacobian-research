#!/usr/bin/env python3
"""Unrestricted ordinary cover cocycles through the new ternary quotients."""

import hashlib
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram
from ternary_transvections import I3,identity,inverse,multiply,rank
from verify_integer_lifts import nullspace

HERE=Path(__file__).resolve().parent


def core_words(pd):
    return [[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in oriented_diagram(pd)]+[[1]]


def ordinary_matrix(matrices,words):
    n=len(matrices);rows=[];inverses=[inverse(a) for a in matrices]
    for word in words:
        block=[[0]*(3*n) for _ in range(3)];prefix=I3
        for letter in word:
            i=abs(letter)-1
            coefficient=prefix if letter>0 else tuple(-v%3 for v in multiply(prefix,inverses[i]))
            for j in range(3):
                for k in range(3):block[j][3*i+k]=(block[j][3*i+k]+coefficient[3*j+k])%3
            prefix=multiply(prefix,matrices[i] if letter>0 else inverses[i])
        assert prefix==I3;rows.extend(block)
    basis=nullspace(rows,3*n,3)
    fixed_rank=rank([[(a[3*i+j]-int(i==j))%3 for j in range(3)] for a in matrices for i in range(3)],3)
    return {'coefficient_matrix':rows,'kernel_basis':basis,'Z1_dimension':len(basis),
        'fixed_space_dimension':3-fixed_rank,'H1_dimension':len(basis)-fixed_rank,
        'rank_lower_bound':(len(basis)+2)//3}


def main():
    started=time.monotonic();plan=json.loads((HERE/'ternary_cover_module_plan.json').read_text())
    checked=json.loads((HERE/'ternary_base_verification.json').read_text())
    assert checked['status']=='PASS independent complete ternary3 rank-one base gate'
    for p,digest in checked['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    source=json.loads((HERE/'ternary_reflection_results.json').read_text())
    catalogue=json.loads((HERE/'ternary_reflection_catalogue.json').read_text())
    elements=[tuple(m) for m in catalogue['elements']]
    panel={c['name']:c for c in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    control=json.loads((HERE/'ternary_reflection_control.json').read_text())
    positive=ordinary_matrix([I3]*9,core_words(control['pd']));assert positive['Z1_dimension']==9 and positive['rank_lower_bound']==3
    files=['ternary_cover_module_plan.json','audit_ternary_cover_modules.py','ternary_transvections.py',
        'verify_integer_lifts.py','parabolic_lifts.py','parabolic_panel.json','ternary_reflection_results.json',
        'ternary_reflection_catalogue.json','ternary_base_verification.json','ternary_reflection_control.json']
    report={'status':'RUNNING','positive_control':positive,'results':[],
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    path=HERE/'ternary_cover_module_results.json'
    for row in source['results']:
        for index,base in enumerate(row['bases']):
            if len(set(base['arc_labels']))==1:continue
            assert time.monotonic()-started<plan['limits']['seconds_total']
            meridians=[elements[i] for i in base['arc_labels']]
            matrices=[multiply(a,meridians[0]) for a in meridians];assert matrices[0]==I3
            dual=[tuple(inverse(a)[3*j+i] for i in range(3) for j in range(3)) for a in matrices]
            words=core_words(panel[row['name']]['pd'])
            modules=[{'module':'natural3','matrices':matrices,**ordinary_matrix(matrices,words)},
                     {'module':'dual3','matrices':dual,**ordinary_matrix(dual,words)}]
            assert all(m['fixed_space_dimension']==0 for m in modules)
            result={'name':row['name'],'representation_index':index,'seed_assignment':base['seed_assignment'],
                'core_words':words,'modules':modules,'sharp_four':any(m['Z1_dimension']>=7 for m in modules)}
            report['results'].append(result);path.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'name':row['name'],'representation':index,
                'Z1':[m['Z1_dimension'] for m in modules],'H1':[m['H1_dimension'] for m in modules],
                'sharp_four':result['sharp_four']}),flush=True)
    assert len(report['results'])==plan['limits']['representation_orbits']==16
    report['status']='COMPLETE';report['sharp_names']=sorted({r['name'] for r in report['results'] if r['sharp_four']})
    report['seconds']=round(time.monotonic()-started,6);path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'sharp_names':report['sharp_names'],'seconds':report['seconds']},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Check both restricted Ext spaces for the frozen stabilized H3 bases."""
import hashlib
import json
from pathlib import Path
import time
from audit_meridian_tangents import exact_certificate,ma,mn
from parabolic_lifts import oriented_diagram
from quadratic_integer_lifts import IDENTITY,ZERO,root_matrix

HERE=Path(__file__).resolve().parent


def transpose(a):
    return tuple(a[3*j+i] for i in range(3) for j in range(3))


def eigenroot(a):
    diff=ma(IDENTITY,mn(a))
    return next(tuple(diff[3*i+j] for i in range(3)) for j in range(3)
                if any(diff[3*i+j]!=ZERO for i in range(3)))


def main():
    plan=json.loads((HERE/'stabilization_plan.json').read_text())
    tangent=json.loads((HERE/'meridian_tangent_results.json').read_text())
    assert plan['frozen_before_execution'] and tangent['status']=='COMPLETE'
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    results=[];start=time.monotonic()
    for row in tangent['results']:
        assert time.monotonic()-start<plan['limits']['seconds_total']
        matrices=[tuple(map(tuple,m)) for m in row['matrices']]
        words=[[o+1,a+1,o+1,b+1] for o,a,b,_ in oriented_diagram(panel[row['name']]['pd'])]
        modules=[]
        for label,images in [('W',matrices),('Wdual',[transpose(m) for m in matrices])]:
            roots=[eigenroot(m) for m in images]
            matrix=root_matrix(images,roots,words)
            cert=exact_certificate(matrix)
            assert cert['nullity']==3
            modules.append({'module':label,'roots':roots,'matrix':matrix,'certificate':cert})
        results.append({'name':row['name'],'orbit_index':row['orbit_index'],'modules':modules})
    paths=['audit_stabilization.py','stabilization_plan.json','audit_meridian_tangents.py',
           'meridian_tangent_results.json','quadratic_integer_lifts.py','parabolic_panel.json','parabolic_lifts.py']
    result={'status':'COMPLETE','results':results,
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'stabilization_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'COMPLETE','bases':len(results),'restricted_modules':sum(len(r['modules']) for r in results),'all_complex_nullities':3},indent=2))


if __name__=='__main__':
    main()

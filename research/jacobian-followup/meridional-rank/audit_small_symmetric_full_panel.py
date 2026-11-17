#!/usr/bin/env python3
"""Complete small symmetric scalar lifts on the frozen293-name panel."""
from collections import Counter
from itertools import combinations, permutations, product
import hashlib
import json
from pathlib import Path
import time

from integer_lifts import relator_matrix, smith_certificate
from joint_a5 import group_model
from parabolic_lifts import oriented_diagram
from search_lift_panel import conjugation_table, extend_colouring
from signed_lifts import sign_potentials

HERE=Path(__file__).resolve().parent


def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))


def seed_orbits(r):
    edges,table=conjugation_table(r);index={e:i for i,e in enumerate(edges)}
    actions=[tuple(index[tuple(sorted((s[a],s[b])))] for a,b in edges)
             for s in permutations(range(r))]
    stabilizer=[a for a in actions if a[0]==0];representatives=[]
    for tail in product(range(len(edges)),repeat=3):
        values=(0,)+tail
        if sign_potentials(r,[edges[i] for i in values],0)[0]!=1:continue
        if values==min(tuple(a[i] for i in values) for a in stabilizer):representatives.append(values)
    return edges,table,representatives


def main():
    plan=json.loads((HERE/'small_symmetric_full_panel_plan.json').read_text())
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases']
    covers=json.loads((HERE/'joint_a5_results.json').read_text())
    elements,mul,inv,orders,autos,_=group_model();sheets=list(permutations(range(5)))
    identities=tuple(range(5));models={r:seed_orbits(r) for r in [2,3,4,5]}
    inputs=['small_symmetric_full_panel_plan.json','audit_small_symmetric_full_panel.py','joint_a5_panel.json',
            'joint_a5_results.json','joint_a5.py','parabolic_lifts.py','search_lift_panel.py','signed_lifts.py','integer_lifts.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    path=HERE/'small_symmetric_full_panel_results.json'
    report=json.loads(path.read_text()) if path.exists() else {'status':'IN_PROGRESS','input_sha256':hashes,'results':[]}
    assert report['input_sha256']==hashes
    done={r['name'] for r in report['results']};started=time.monotonic()
    for case,old in zip(panel,covers['results']):
        if case['name'] in done:continue
        if time.monotonic()-started>plan['limits']['seconds_total']:break
        tick=time.monotonic();assert case['name']==old['name'] and old['status']=='COMPLETE'
        crossings=[c[:3] for c in oriented_diagram(case['pd'])];n=len(crossings)
        words=[[o+1,a+1,o+1,b+1] for o,a,b in crossings]
        deck=[]
        for m in old['epimorphism_representatives']:
            labels=m['arc_images'];implementers=[j for j,a in enumerate(autos) if all(a[i]==inv[i] for i in labels)]
            assert len(implementers)<=1
            item={'orbit_index':m['orbit_index'],'deck_type':'absent','implementer':None}
            if implementers:
                d=sheets[implementers[0]];parity=sum(d[i]>d[j] for i in range(5) for j in range(i+1,5))%2
                assert compose(d,d)==identities
                meridians=[compose(elements[i],d) for i in labels]
                assert all(compose(a,a)==identities for a in meridians)
                item.update({'deck_type':'outer' if parity else 'inner','implementer':d,'meridians':meridians})
            deck.append(item)
        result={'name':case['name'],'status':'COMPLETE','seeds':old['seeds'],'seed_trace':old['seed_trace'],
                'bridge_upper_bound':4,'A5_deck_extensions':deck,'degrees':[]}
        for r,(edges,table,representatives) in models.items():
            certificates=[]
            for assignment in representatives:
                assert time.monotonic()-tick<plan['limits']['seconds_per_target']
                labels=extend_colouring(n,crossings,old['seeds'],assignment,table)
                if labels is None:continue
                matrix=relator_matrix(r,[edges[i] for i in labels],words)
                smith=smith_certificate(matrix,n)
                assert smith['maximum_dimension_all_characteristics']<=4
                certificates.append({'seed_assignment':assignment,'arc_labels':labels,'matrix':matrix,'smith':smith})
            result['degrees'].append({'degree':r,'seed_orbits_tested':len(representatives),'certificates':certificates})
        # S5 conjugacy classes are exactly the outer-invariant A5 kernel classes.
        assert len(result['degrees'][-1]['certificates'])==sum(d['deck_type']=='outer' for d in deck)
        result['maximum_bound_all_degrees_fields']=max((c['smith']['maximum_dimension_all_characteristics']
            for d in result['degrees'] for c in d['certificates']),default=0)
        result['sharp_four']=result['maximum_bound_all_degrees_fields']==4
        result['seconds']=round(time.monotonic()-tick,6);report['results'].append(result)
        path.write_text(json.dumps(report,indent=2)+'\n')
        if result['sharp_four'] or len(report['results'])%50==0:
            print(json.dumps({'completed':len(report['results']),'name':case['name'],'bases':[len(d['certificates']) for d in result['degrees']], 'sharp':result['sharp_four']}),flush=True)
    report['status']='COMPLETE' if len(report['results'])==293 else 'CAPPED_INCOMPLETE'
    report['seed_orbit_counts']={r:len(v[2]) for r,v in models.items()}
    report['deck_histogram']=dict(Counter(d['deck_type'] for row in report['results'] for d in row['A5_deck_extensions']))
    report['maximum_bound_histogram']=dict(Counter(row['maximum_bound_all_degrees_fields'] for row in report['results']))
    report['sharp_names']=[r['name'] for r in report['results'] if r['sharp_four']]
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','seed_orbit_counts','deck_histogram','maximum_bound_histogram','sharp_names']},indent=2))


if __name__=='__main__':main()

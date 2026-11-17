#!/usr/bin/env python3
"""Frozen original16 order3 transvection base and affine-root search."""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import time

from double_transpositions import tuple_representatives
from parabolic_lifts import extend,oriented_diagram
from ternary_transvections import determinant,inverse,lift_matrix,model,multiply

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic();plan=json.loads((HERE/'ternary_transvection_plan.json').read_text())
    positive=json.loads((HERE/'ternary_transvection_control_verification.json').read_text())
    assert positive['status']=='PASS independent ternary triple-trefoil control'
    for p,digest in positive['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    inputs=['ternary_transvection_plan.json','audit_ternary_transvections.py','ternary_transvections.py',
        'double_transpositions.py','parabolic_lifts.py','verify_integer_lifts.py','parabolic_panel.json',
        'parabolic_panel_results.json','ternary_transvection_control.json','ternary_transvection_control_verification.json']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    path=HERE/'ternary_transvection_results.json'
    report={'status':'RUNNING','input_sha256':hashes,'results':[]}
    if path.exists():
        report=json.loads(path.read_text());assert report['input_sha256']==hashes
    elements,forward,backward,centralizer,actions=model()
    representatives,levels=tuple_representatives(actions);index={m:i for i,m in enumerate(elements)}
    assert len(representatives)<=plan['limits']['tuple_orbit_cap']
    assert sum(w for _,_,w in representatives)==plan['limits']['raw_normalized_tuples']
    report['tuple_catalogue']={'class_size':104,'centralizer_size':108,'prefix_counts':levels,
        'four_tuple_orbits':len(representatives),'raw_normalized_tuples':sum(w for _,_,w in representatives),
        'stabilizer_histogram':dict(sorted(Counter(s for _,s,_ in representatives).items()))}
    catalogue={'elements':elements,'centralizer':centralizer,
        'representatives':[{'tuple':t,'stabilizer':s,'raw_weight':w} for t,s,w in representatives],
        'counts':report['tuple_catalogue']}
    catalogue_path=HERE/'ternary_transvection_catalogue.json'
    catalogue_path.write_text(json.dumps(catalogue,separators=(',',':'))+'\n')
    report['catalogue_sha256']=hashlib.sha256(catalogue_path.read_bytes()).hexdigest()
    # End-to-end positive on the exact independently checked control PD.
    control=json.loads((HERE/'ternary_transvection_control.json').read_text())
    seed_matrices=[elements[control['arc_linear_labels'][s]] for s in control['seeds']]
    g=next(g for g in product(range(3),repeat=9) if determinant(g) and
        multiply(multiply(g,seed_matrices[0]),inverse(g))==elements[0])
    assignment=tuple(index[multiply(multiply(g,a),inverse(g))] for a in seed_matrices)
    canonical=min(tuple(action[i] for i in assignment) for action in actions)
    assert canonical in {t for t,_,_ in representatives}
    cr=oriented_diagram(control['pd'])
    labels=extend(9,cr,control['seeds'],canonical,forward,backward);assert labels is not None
    _,basis,_=lift_matrix(cr,labels,elements);assert len(basis)==4
    report['end_to_end_control']={'normalizing_matrix':g,'normalized_seed_tuple':assignment,
        'canonical_seed_tuple':canonical,'arc_labels':labels,'root_dimension':len(basis)}
    panel=json.loads((HERE/'parabolic_panel.json').read_text())['cases']
    old=json.loads((HERE/'parabolic_panel_results.json').read_text())['results']
    assert [r['name'] for r in panel]==plan['targets']
    done={r['name'] for r in report['results']}
    for case,previous in zip(panel,old):
        if case['name'] in done:continue
        if time.monotonic()-started>plan['limits']['seconds_total']:break
        assert case['name']==previous['name'] and previous['status']=='COMPLETE'
        tick=time.monotonic();crossings=oriented_diagram(case['pd']);valid=[];tested=0
        for assignment,stabilizer,weight in representatives:
            if time.monotonic()-tick>plan['limits']['seconds_per_target']:break
            tested+=1
            labels=extend(len(crossings),crossings,previous['seeds'],assignment,forward,backward)
            if labels is None:continue
            matrix,basis,roots=lift_matrix(crossings,labels,elements)
            assert len(basis)<=4,'root bound contradicts certified upper4'
            valid.append({'seed_assignment':assignment,'arc_labels':labels,'stabilizer':stabilizer,
                'raw_weight':weight,'root_dimension':len(basis),'root_basis':basis,'root_matrix':matrix,'roots':roots})
        row={'name':case['name'],'status':'COMPLETE' if tested==len(representatives) else 'CAPPED_INCOMPLETE',
            'seeds':previous['seeds'],'seed_trace':previous['seed_trace'],'bridge_upper_bound':4,
            'seed_orbits_tested':tested,'bases':valid,'root_dimension_histogram':dict(Counter(r['root_dimension'] for r in valid)),
            'maximum_root_dimension':max((r['root_dimension'] for r in valid),default=0),
            'sharp_four':any(r['root_dimension']==4 for r in valid),'seconds':round(time.monotonic()-tick,6)}
        report['results'].append(row);path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:row[k] for k in ['name','status','seed_orbits_tested','root_dimension_histogram','maximum_root_dimension','seconds']}),flush=True)
    report['status']='COMPLETE' if len(report['results'])==16 and all(r['status']=='COMPLETE' for r in report['results']) else 'CAPPED_INCOMPLETE'
    report['sharp_names']=[r['name'] for r in report['results'] if r['sharp_four']]
    report['maximum_dimension_histogram']=dict(Counter(r['maximum_root_dimension'] for r in report['results']))
    report['seconds_total']=round(time.monotonic()-started,6)
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','tuple_catalogue','sharp_names','maximum_dimension_histogram','seconds_total']},indent=2))


if __name__=='__main__':main()

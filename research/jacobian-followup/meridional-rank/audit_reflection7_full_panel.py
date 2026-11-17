#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
from audit_outer_reflection7 import projective_model,identify,adjoint_twist,roots,det
from gl3f2 import group_model
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'reflection7_full_panel_plan.json').read_text());assert plan['frozen_before_execution']
    source=json.loads((HERE/'gl3f2_full_panel_results.json').read_text());assert source['status']=='COMPLETE'
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    _,oldtable,_,_,_,_,_,_,std=group_model()
    elements,table,inverses,orders,even=projective_model();identification=identify(oldtable,std,table,orders,even)
    # This list gives -Ad(A), not a homomorphism on PGL. Meridian signs come
    # from the knot's independent parity character; all Wirtinger words are checked.
    reps=[tuple((-pow(det(a),3,7)*x)%7 for x in adjoint_twist(a)) for a in elements]
    results=[]
    for row in source['results']:
        words=[[o+1,a+1,o+1,b+1] for o,a,b,_ in oriented_diagram(panel[row['name']]['pd'])]
        for mapping in row['maps']:
            arcs=[identification[i] for i in mapping['arc_images']]
            deck=[t for t in range(336) if orders[t]==2 and all(table[table[t][g]][inverses[t]]==inverses[g] for g in arcs)]
            assert len(deck)==1
            assert (deck[0] in even)==(mapping['deck_type']=='inner')
            meridians=[table[g][deck[0]] for g in arcs];assert all(orders[m]==2 for m in meridians)
            images=[reps[m] for m in meridians]
            result={'name':row['name'],'orbit_index':mapping['orbit_index'],'deck_type':mapping['deck_type'],
                    'deck_projective_index':deck[0],'projective_meridians':meridians,'meridian_matrices':images,
                    **roots(images,words)}
            results.append(result)
    assert len(results)==191
    paths=['reflection7_full_panel_plan.json','audit_reflection7_full_panel.py','audit_outer_reflection7.py',
           'gl3f2.py','gl3f2_full_panel_results.json','joint_a5_panel.json','parabolic_lifts.py','verify_integer_lifts.py']
    report={'status':'COMPLETE','projective_elements':elements,'GL_to_PSL_identification':identification,
            'results':results,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'reflection7_full_panel_results.json').write_text(json.dumps(report,indent=2)+'\n')
    from collections import Counter
    print(json.dumps({'status':'COMPLETE','maps':len(results),'root_dimension_histogram':dict(Counter(r['dimension'] for r in results)),
                      'sharp':[{'name':r['name'],'orbit_index':r['orbit_index'],'dimension':r['dimension']} for r in results if r['sharp_four']]},indent=2))


if __name__=='__main__':main()

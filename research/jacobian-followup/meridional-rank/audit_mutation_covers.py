#!/usr/bin/env python3
"""Frozen Tait-graph certificates for the two proposed common-cover blocks."""
import hashlib
import json
from pathlib import Path
import time

from mutation_graphs import checkerboards,search

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic();plan=json.loads((HERE/'mutation_cover_plan.json').read_text())
    panel={c['name']:c for c in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    files=['mutation_cover_plan.json','mutation_graphs.py','audit_mutation_covers.py','parabolic_panel.json',
        'projective13_source_audit.json','sources/knotinfo/2026.11.2_projective13_geometry_subset.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'diagrams':{},'controls':[],'pairs':[]}
    path=HERE/'mutation_cover_results.json'
    for names in plan['selection'].values():
        for name in names:
            pd=panel[name]['pd'];data=checkerboards(pd);assert len(pd)==13 and data['selected']['vertices']<=7
            assert data['selected']['two_vertex_connected']
            report['diagrams'][name]={'pd':pd,**data}
    # Exercise a genuine switch: its endpoint degree multiset changes.
    synthetic=[[0,2],[0,2],[2,3],[3,1],[0,3],[0,4],[4,5],[5,1],[5,1],[4,1],[0,5],[0,1],[0,1]]
    side=[0,1,2,3,4];swapped=[[1 if a==0 else 0 if a==1 else a,1 if b==0 else 0 if b==1 else b] if i in side else [a,b] for i,(a,b) in enumerate(synthetic)]
    control={'vertices':6,'edges':synthetic};target={'vertices':6,'edges':swapped}
    certificate=search(control,target,plan['limits']);assert certificate['status'].startswith('CERTIFIED') and certificate['switches']
    report['controls'].append({'name':'synthetic_one_switch','source':control,'target':target,'certificate':certificate})
    from spherogram import RationalTangle
    pretzels=[]
    for slopes in [[2,3,3,5],[2,3,5,3]]:
        tangle=RationalTangle(0)
        for q in slopes:tangle=tangle+RationalTangle(1,q)
        link=tangle.numerator_closure();pd=[list(c) for c in link.PD_code()]
        assert len(link.link_components)==1 and len(pd)==13
        pretzels.append({'slopes':slopes,'pd':pd,**checkerboards(pd)})
    certificate=search(pretzels[0]['selected'],pretzels[1]['selected'],plan['limits'])
    assert certificate['status'].startswith('CERTIFIED')
    report['controls'].append({'name':'positive_pretzel_pair','diagrams':pretzels,'source':pretzels[0]['selected'],'target':pretzels[1]['selected'],'certificate':certificate})
    for label,names in plan['selection'].items():
        source=report['diagrams'][names[0]]['selected']
        for name in names[1:]:
            assert time.monotonic()-started<plan['limits']['seconds_total']
            target=report['diagrams'][name]['selected'];certificate=search(source,target,plan['limits'])
            same=source['checkerboard_parity']==target['checkerboard_parity']
            row={'block':label,'source':names[0],'target':name,'checkerboard_parities':[source['checkerboard_parity'],target['checkerboard_parity']],
                'requires_mirror_for_same_convention':not same,'certificate':certificate}
            report['pairs'].append(row);path.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'source':names[0],'target':name,'status':certificate['status'],'switches':len(certificate.get('switches',[])),
                'states':certificate.get('states'),'requires_mirror':not same}),flush=True)
    report['status']='COMPLETE_CERTIFICATES' if all(r['certificate']['status'].startswith('CERTIFIED') for r in report['pairs']) else 'PARTIAL_CERTIFICATES'
    report['seconds']=round(time.monotonic()-started,6)
    report['assurance']='Finite graph certificates only; mutation/common-cover conclusion also uses read primary graph-to-mutation and torus-lift arguments. No actual group rank3 or target equality is asserted.'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'seconds':report['seconds'],'pairs':len(report['pairs'])},indent=2))


if __name__=='__main__':main()

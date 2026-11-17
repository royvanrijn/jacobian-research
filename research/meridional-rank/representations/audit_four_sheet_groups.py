#!/usr/bin/env python3
"""Exact perfect four-sheet cover group presentations and a bounded A5 gate."""
import hashlib
from itertools import combinations,product
import json
from pathlib import Path
import time
from audit_three_sheet_nonabelian import eliminate,finite_group
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'four_sheet_group_plan.json').read_text())
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    source=json.loads((HERE/'small_symmetric_full_panel_results.json').read_text());edges=list(combinations(range(4),2))
    elements,table,inv=finite_group('A5');rows=[];started=time.monotonic()
    for case in source['results']:
        for certificate in case['degrees'][2]['certificates']:
            assert time.monotonic()-started<plan['limits']['seconds_total'];tick=time.monotonic()
            name=case['name'];labels=certificate['arc_labels'];cr=[c[:3] for c in oriented_diagram(panel[name]['pd'])]
            parent=list(range(4));tree=[]
            def root(i):
                while parent[i]!=i:i=parent[i]
                return i
            for i,label in enumerate(labels,1):
                a,b=edges[label]
                if root(a)!=root(b):parent[root(a)]=root(b);tree.append(i)
            assert len(tree)==3
            relators=[]
            for o,a,b in cr:
                for start in range(4):
                    sheet=start;word=[]
                    for i in [o,a,o,b]:
                        p,q=edges[labels[i]]
                        if sheet in (p,q):
                            if i+1 not in tree:word.append(i+1 if sheet==p else -(i+1))
                            sheet=q if sheet==p else p
                    assert sheet==start;relators.append(word)
            generators=[i+1 for i in range(len(labels)) if i+1 not in tree]
            remaining,words,steps=eliminate(generators,relators)
            record={'name':name,'arc_labels':labels,'tree_generators':tree,'generators':generators,'relators':relators,
                    'eliminations':steps,'remaining_generators':remaining,'remaining_relators':words}
            if not remaining:
                assert not words;record['status']='PROVED_TRIVIAL_BY_ELIMINATION'
            elif len(remaining)<=plan['limits']['A5_max_generators']:
                witness=None;tested=0
                for values in product(range(60),repeat=len(remaining)):
                    tested+=1
                    if time.monotonic()-tick>plan['limits']['seconds_per_cover']:break
                    if not any(values):continue
                    assignment=dict(zip(remaining,values));valid=True
                    for word in words:
                        value=0
                        for a in word:value=table[value][assignment[a] if a>0 else inv[assignment[-a]]]
                        if value:valid=False;break
                    if not valid:continue
                    for step in reversed(steps):
                        value=0
                        for a in step['replacement']:value=table[value][assignment[a] if a>0 else inv[assignment[-a]]]
                        assignment[step['generator']]=value
                    witness={'A5_elements':elements,'original_generator_images':assignment};break
                record.update({'status':'NONTRIVIAL_WITNESS' if witness else 'A5_EXHAUSTED_SOURCE_UNRESOLVED' if tested==60**len(remaining) else 'CAPPED_INCOMPLETE',
                               'A5_assignments_tested':tested,'witness':witness})
            else:record['status']='GENERATOR_CAP_SOURCE_UNRESOLVED'
            rows.append(record)
            report={'status':'IN_PROGRESS','results':rows,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['four_sheet_group_plan.json','audit_four_sheet_groups.py','audit_three_sheet_nonabelian.py','parabolic_lifts.py','joint_a5_panel.json','small_symmetric_full_panel_results.json']}}
            (HERE/'four_sheet_group_results.json').write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'name':name,'remaining_generators':len(remaining),'status':record['status']}),flush=True)
    assert len(rows)==21
    report['status']='COMPLETE';(HERE/'four_sheet_group_results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

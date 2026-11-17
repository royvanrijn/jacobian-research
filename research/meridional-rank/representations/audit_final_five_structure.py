#!/usr/bin/env python3
"""Exact source audit of the final five original targets, with three controls."""
import hashlib
import json
from pathlib import Path
import time
from compact_cover_groups import free,simplify
from mutation_graphs import checkerboards
from parabolic_lifts import oriented_diagram
from probe_outer_tangles import diagram
from audit_remaining_cover_structure import candidate_cuts

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'final_five_structure_plan.json').read_text());started=time.monotonic()
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    files=['final_five_structure_plan.json','audit_final_five_structure.py','parabolic_panel.json',
        'compact_cover_groups.py','mutation_graphs.py','parabolic_lifts.py','probe_outer_tangles.py',
        'audit_remaining_cover_structure.py','three_piece_cover_verification.json','three_piece_source_audit.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'results':[]}
    path=HERE/'final_five_structure_results.json';limits={**plan['limits'],'seconds_per_target':plan['limits']['seconds_per_group']}
    for name in plan['controls']+plan['targets']:
        assert time.monotonic()-started<plan['limits']['seconds_total']
        source_started=time.monotonic();pd=panel[name]['pd'];triples=oriented_diagram(pd)
        rels=[free([x for x in (o+1,-a-1,o+1,-b-1) if abs(x)!=1]) for o,a,b,_ in triples]
        row={'name':name,'pd':pd,'checkerboards':checkerboards(pd),'full_raw_core_relators':rels,
            'compact_core':simplify(list(range(2,14)),rels,limits),'cuts':[]}
        cuts=candidate_cuts(pd);row['candidate_cuts']=cuts
        for cut in cuts[:plan['limits']['fragment_cuts_reduced_each']]:
            if time.monotonic()-source_started>plan['limits']['seconds_per_source'] or time.monotonic()-started>plan['limits']['seconds_total']:
                row['cut_status']='TIME_CAP';break
            parts=[]
            for field in ('inside_crossings','outside_crossings'):
                part=diagram([pd[i] for i in cut[field]]);assert part['boundary_labels']==cut['boundary']
                parts.append({**part,'formal_core_reduction':simplify(part['generators'],part['raw_relators'],limits)})
            row['cuts'].append({**cut,'fragments':parts})
        else:row['cut_status']='ALL_CANDIDATES' if len(cuts)<=plan['limits']['fragment_cuts_reduced_each'] else 'CANDIDATE_CAP'
        row['seconds']=round(time.monotonic()-source_started,6);report['results'].append(row)
        path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'name':name,'full':row['compact_core']['final'],
            'cuts':[{'inside':c['inside_crossings'],'boundary':c['boundary'],
                'fragments':[f['formal_core_reduction']['final'] for f in c['fragments']]} for c in row['cuts']],
            'cut_status':row['cut_status'],'seconds':row['seconds']}),flush=True)
    report['status']='COMPLETE_FROZEN_SELECTION';report['seconds']=round(time.monotonic()-started,6)
    report['assurance']='Exact finite source/presentation inputs only; independent replay and explicit peripheral/topology argument required. No new knot closure, rank bound or group transfer from this audit alone.'
    path.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

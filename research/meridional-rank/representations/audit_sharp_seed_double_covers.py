#!/usr/bin/env python3
"""Frozen double-cover source/fragment audit of the four sharp simple3cover seeds."""
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
    plan=json.loads((HERE/'sharp_seed_double_cover_plan.json').read_text());started=time.monotonic()
    targets={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    controls={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    provenance={r['name']:r for r in json.loads((HERE/'splice_knot_provenance.json').read_text())['results']}
    files=['sharp_seed_double_cover_plan.json','audit_sharp_seed_double_covers.py','joint_a5_panel.json',
        'parabolic_panel.json','splice_knot_provenance.json','compact_cover_groups.py','mutation_graphs.py',
        'parabolic_lifts.py','probe_outer_tangles.py','audit_remaining_cover_structure.py','graph_cover_family_scope_verification.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'results':[]}
    path=HERE/'sharp_seed_double_cover_results.json';limits={**plan['limits'],'seconds_per_target':plan['limits']['seconds_per_group']}
    for name in plan['controls']+plan['targets']:
        assert time.monotonic()-started<plan['limits']['seconds_total']
        source_started=time.monotonic();pd=(controls if name in plan['controls'] else targets)[name]['pd']
        if name in plan['targets']:assert pd==[[e+1 for e in row] for row in provenance[name]['census_pd']]
        triples=oriented_diagram(pd)
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
        print(json.dumps({'name':name,'full':row['compact_core']['final'],'cut_status':row['cut_status'],
            'cuts':[{'inside':c['inside_crossings'],'boundary':c['boundary'],
                'fragments':[f['formal_core_reduction']['final'] for f in c['fragments']]} for c in row['cuts']],
            'seconds':row['seconds']}),flush=True)
    report['status']='COMPLETE_FROZEN_SELECTION';report['seconds']=round(time.monotonic()-started,6)
    report['assurance']='Exact source/presentation inputs only. Independent replay and a written primitive peripheral/manifold argument required before actual-cover or known-family classification. Existing sharp knot equalities are not re-proved or expanded by this audit alone.'
    path.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Frozen small source/partition audit with exact presentation traces."""
from collections import Counter,defaultdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import time
from compact_cover_groups import free,simplify,graph_core
from mutation_graphs import checkerboards
from parabolic_lifts import oriented_diagram
from probe_outer_tangles import diagram

HERE=Path(__file__).resolve().parent


def candidate_cuts(pd):
    endpoints=defaultdict(list)
    for c,row in enumerate(pd):
        for edge in row:endpoints[edge].append(c)
    assert all(len(ends)==2 for ends in endpoints.values())
    def connected(vertices):
        reached={min(vertices)};queue=list(reached)
        for c in queue:
            for ends in endpoints.values():
                if c in ends:
                    for d in ends:
                        if d in vertices and d not in reached:reached.add(d);queue.append(d)
        return reached==vertices
    result=[]
    for size in (6,7):
        for tail in combinations(range(1,13),size-1):
            left={0,*tail};right=set(range(13))-left
            boundary=sorted(e for e,ends in endpoints.items() if (ends[0] in left)!=(ends[1] in left))
            if len(boundary)!=4 or not connected(left) or not connected(right):continue
            result.append({'inside_crossings':sorted(left),'outside_crossings':sorted(right),'boundary':boundary})
    return sorted(result,key=lambda r:(r['inside_crossings']!=list(range(7)),r['inside_crossings']))


def main():
    plan=json.loads((HERE/'remaining_cover_structure_plan.json').read_text());started=time.monotonic()
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    files=['remaining_cover_structure_plan.json','audit_remaining_cover_structure.py','parabolic_panel.json',
        'compact_cover_groups.py','mutation_graphs.py','parabolic_lifts.py','probe_outer_tangles.py',
        'four_piece_cover_verification.json','four_piece_source_audit.json']
    report={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},'results':[]}
    path=HERE/'remaining_cover_structure_results.json'
    limits={**plan['limits'],'seconds_per_target':plan['limits']['seconds_per_group']}
    for name in plan['controls']+plan['representatives']:
        assert time.monotonic()-started<plan['limits']['seconds_total']
        source_started=time.monotonic();pd=panel[name]['pd'];graphs=checkerboards(pd)
        triples=oriented_diagram(pd);rels=[free([x for x in (o+1,-a-1,o+1,-b-1) if abs(x)!=1]) for o,a,b,_ in triples]
        row={'name':name,'pd':pd,'checkerboards':graphs,'full_raw_core_relators':rels,
            'compact_core':simplify(list(range(2,14)),rels,limits),'cuts':[]}
        try:row['K4_core']=graph_core(graphs['selected'])
        except AssertionError:row['K4_core']=None
        cuts=candidate_cuts(pd);row['candidate_cuts']=cuts
        for cut in cuts[:plan['limits']['fragment_cuts_reduced_each']]:
            if time.monotonic()-source_started>plan['limits']['seconds_per_source'] or time.monotonic()-started>plan['limits']['seconds_total']:
                row['cut_status']='TIME_CAP';break
            parts=[]
            for field in ['inside_crossings','outside_crossings']:
                part=diagram([pd[i] for i in cut[field]])
                assert part['boundary_labels']==cut['boundary']
                parts.append({**part,'formal_core_reduction':simplify(part['generators'],part['raw_relators'],limits)})
            row['cuts'].append({**cut,'fragments':parts})
        else:row['cut_status']='ALL_CANDIDATES' if len(cuts)<=plan['limits']['fragment_cuts_reduced_each'] else 'CANDIDATE_CAP'
        row['seconds']=round(time.monotonic()-source_started,6);report['results'].append(row)
        path.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'name':name,'full_generators':len(row['compact_core']['final']['generators']),
            'full_relator_lengths':list(map(len,row['compact_core']['final']['relators'])),'K4':row['K4_core'] is not None,
            'cuts_total':len(cuts),'cuts_checked':len(row['cuts']),
            'fragment_generator_relator_counts':[[[len(p['formal_core_reduction']['final']['generators']),len(p['formal_core_reduction']['final']['relators'])] for p in c['fragments']] for c in row['cuts']],
            'seconds':row['seconds']}),flush=True)
    report['status']='COMPLETE_FROZEN_SELECTION';report['seconds']=round(time.monotonic()-started,6)
    report['assurance']='Exact finite source/presentation inputs; independent replay required. No new knot, rank lower bound, manifold model, common-cover or rationality claim.'
    path.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

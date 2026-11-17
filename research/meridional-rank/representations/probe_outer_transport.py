#!/usr/bin/env python3
"""Exact boundary-generation gate, with a fail-closed stalled result."""
import hashlib
import json
from pathlib import Path
from compact_cover_groups import free,inv,substitute
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'outer_tangle_transport_plan.json').read_text())
    source=json.loads((HERE/'mutation_cover_results.json').read_text())['diagrams']
    fragments=json.loads((HERE/'outer_tangle_results.json').read_text())['fragments'];results=[]
    for row in fragments:
        pd=source[row['name']]['pd'];crossings=oriented_diagram(pd);edge_arcs={}
        for (o,a,b,_),record in zip(crossings,pd):
            for edge,arc in [(record[0],a),(record[1],o),(record[2],b),(record[3],o)]:
                assert edge not in edge_arcs or edge_arcs[edge]==arc;edge_arcs[edge]=arc
        transported=[]
        for edges in row['arc_edge_classes']:
            names={edge_arcs[e] for e in edges};assert len(names)==1;transported.append(names.pop())
        reduction=row['formal_core_reduction'];mapping={int(g):w for g,w in reduction['original_generators_in_final'].items()};mapping[1]=[]
        known={};initial=[]
        for i,arc in enumerate(transported):
            initial.append({'outside_arc':i+1,'full_arc':arc,'word':mapping[i+1]})
            if arc not in known:known[arc]=mapping[i+1]
        trace=[]
        while True:
            changed=False
            for i,(o,a,b,_) in enumerate(crossings):
                if o not in known:continue
                if a in known and b not in known:
                    known[b]=free(known[o]+inv(known[a])+known[o]);trace.append({'crossing':i,'over':o,'under_known':a,'under_new':b});changed=True
                if b in known and a not in known:
                    known[a]=free(known[o]+inv(known[b])+known[o]);trace.append({'crossing':i,'over':o,'under_known':b,'under_new':a});changed=True
            if not changed:break
        missing=sorted(set(range(13))-set(known))
        results.append({'name':row['name'],'outside_to_full_arcs':transported,'outside_base_full_arc':transported[0],
            'initial':initial,'trace':trace,'known_words':known,'unresolved_full_arcs':missing,
            'status':'CANDIDATE_TWO_GENERATOR_UPPER' if not missing else 'BOUNDARY_GENERATION_STALLED',
            'assurance':'A stall concerns only this explicit propagation, not actual group rank or every boundary-generation method.'})
    files=['outer_tangle_transport_plan.json','probe_outer_transport.py','outer_tangle_results.json','mutation_cover_results.json','parabolic_lifts.py','compact_cover_groups.py']
    result={'status':'COMPLETE_GATE','results':results,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'outer_tangle_transport.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ['name','status','outside_base_full_arc','unresolved_full_arcs']} for r in results]))


if __name__=='__main__':main()

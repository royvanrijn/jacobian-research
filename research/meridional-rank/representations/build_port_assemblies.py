#!/usr/bin/env python3
"""Exact formal core assembly with four retained boundary identifications."""
import hashlib
import json
from pathlib import Path
from compact_cover_groups import free,inv,substitute,simplify

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'port_assembly_plan.json').read_text())
    inner=json.loads((HERE/'inner_tangle_results.json').read_text())['fragments']
    outer=json.loads((HERE/'outer_tangle_results.json').read_text())['fragments']
    limits=json.loads((HERE/'compact_cover_plan.json').read_text())['limits'];rows=[]
    for I,O in zip(inner,outer):
        assert I['name']==O['name'] and I['boundary_labels']==O['boundary_labels']
        pieces=[]
        for data,newgens in [(I,[1,2]),(O,[3,4])]:
            reduced=data['formal_core_reduction'];gens=reduced['final']['generators'];assert len(gens)==2
            rename={g:[h] for g,h in zip(gens,newgens)}
            original={int(g):substitute(w,rename) for g,w in reduced['original_generators_in_final'].items()};original[1]=[]
            pieces.append({'generators':newgens,'relators':[substitute(w,rename) for w in reduced['final']['relators']],
                'original_arc_words':original})
        ports=[]
        for edge in I['boundary_labels']:
            i=next(j+1 for j,s in enumerate(I['arc_edge_classes']) if edge in s)
            o=next(j+1 for j,s in enumerate(O['arc_edge_classes']) if edge in s)
            u=pieces[0]['original_arc_words'][i];v=pieces[1]['original_arc_words'][o]
            ports.append({'edge':edge,'inner_arc':i,'outer_arc':o,'inner_word':u,'outer_word':v,
                'relation':free(u+[5]+inv(v))})
        rels=pieces[0]['relators']+pieces[1]['relators']+[p['relation'] for p in ports]
        rows.append({'name':I['name'],'pieces':pieces,'ports':ports,
            'assembly_reduction':simplify([1,2,3,4,5],rels,limits)})
    files=['port_assembly_plan.json','build_port_assemblies.py','inner_tangle_results.json','outer_tangle_results.json','compact_cover_groups.py','compact_cover_plan.json']
    report={'status':'COMPLETE_ASSEMBLIES','targets':rows,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'assurance':'Core-presentation base change and port gluing only; no injective peripheral amalgam or topological piece classification.'}
    (HERE/'port_assembly_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{'name':r['name'],'ports':r['ports'],'final':r['assembly_reduction']['final']} for r in rows]))


if __name__=='__main__':main()

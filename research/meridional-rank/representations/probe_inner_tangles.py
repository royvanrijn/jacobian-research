#!/usr/bin/env python3
"""Fixed complementary7crossing formal presentations and port data."""
import hashlib
import json
from pathlib import Path
from compact_cover_groups import simplify
from probe_outer_tangles import diagram

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'inner_tangle_plan.json').read_text());sources=json.loads((HERE/'mutation_cover_results.json').read_text())['diagrams']
    limits=json.loads((HERE/'compact_cover_plan.json').read_text())['limits'];rows=[]
    for name in plan['targets']:
        data=diagram(sources[name]['pd'][:7]);reduced=simplify(data['generators'],data['raw_relators'],limits)
        rows.append({'name':name,**data,'formal_core_reduction':reduced})
    files=['inner_tangle_plan.json','probe_inner_tangles.py','probe_outer_tangles.py','compact_cover_groups.py','compact_cover_plan.json','mutation_cover_results.json']
    result={'status':'COMPLETE_FORMAL_FRAGMENTS','fragments':rows,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'assurance':'Formal complementary presentations and boundary labels only; no target group amalgam, JSJ or rank conclusion.'}
    (HERE/'inner_tangle_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{'name':r['name'],'boundary':r['boundary_labels'],'final':r['formal_core_reduction']['final']} for r in rows]))


if __name__=='__main__':main()

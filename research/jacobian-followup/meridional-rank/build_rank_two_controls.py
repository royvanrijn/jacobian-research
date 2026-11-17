#!/usr/bin/env python3
"""Freeze independent actual diagrams for known pretzel controls."""
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
from spherogram import RationalTangle
from audit_diagram_identity import pd_walk

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/'rank_two_equivalence_plan.json').read_text());assert version('spherogram')=='2.4.1'
    rows=[]
    for parameters in plan['finite_controls']['pretzel_parameters']:
        tangle=RationalTangle(0)
        for n in parameters:tangle=tangle+RationalTangle(1,n)
        link=tangle.numerator_closure();assert len(link.link_components)==1
        pd=[list(c) for c in link.PD_code()]
        code=[(c+1)*(1 if over else -1) for c,over,_ in pd_walk(pd)]
        rows.append({'parameters':parameters,'pd':pd,'gauss_code':code,'bridge_upper':len(parameters),
                     'construction':'Numerator closure of the sum of reciprocal rational tangles; known pretzel family'})
    report={'status':'CONSTRUCTED','spherogram_version':version('spherogram'),'controls':rows,
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['build_rank_two_controls.py','rank_two_equivalence_plan.json','audit_diagram_identity.py']}}
    (HERE/'rank_two_control_diagrams.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'CONSTRUCTED','crossing_counts':[len(r['pd']) for r in rows]}))


if __name__=='__main__':main()

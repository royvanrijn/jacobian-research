#!/usr/bin/env python3
"""Targeted exact conversion to two trefoil peripheral-amalgam models."""
import hashlib
import json
from pathlib import Path
from audit_three_sheet_structure import inv,subst,normalize,shorten,measure

HERE=Path(__file__).resolve().parent


def models(kind,e):
    a=(1,);b=(2,);ap=(1,)*(3*e) if e>0 else (-1,)*(-3*e)
    if kind=='A':return normalize([(2,)*3+(1,)*2+(2,)*3+(-1,),ap+b+ap+(-2,)*2])
    return normalize([(2,)*2+a+(2,)*2+(-1,)*2,(ap+b)*3+(-2,)*2])


def main():
    plan=json.loads((HERE/'three_sheet_splice_plan.json').read_text());source=json.loads((HERE/'three_sheet_structural_results.json').read_text());rows=[]
    for old in source['results']:
        if old['name'] not in plan['model_assignments']:continue
        kind,e=plan['model_assignments'][old['name']];words=normalize([tuple(w) for w in old['final_relators']]);steps=[]
        if old['name']=='13a_3148':
            forward={1:(1,),2:(2,-1,-1)};backward={1:(1,),2:(2,1,1)}
            steps.append({'operation':'nielsen','forward':forward,'inverse':backward,'before_measure':measure(words)})
            words=normalize([subst(w,forward) for w in words]);steps[-1]['after_measure']=measure(words)
        for _ in range(plan['limits']['steps_per_case']):
            candidate=shorten(words)
            if candidate is None:break
            new,step=candidate;step.update({'before_measure':measure(words),'after_measure':measure(new)});steps.append(step);words=new
        target=models(kind,e)
        rows.append({'name':old['name'],'model':kind,'epsilon':e,'initial_relators':old['final_relators'],
                     'steps':steps,'final_relators':words,'target_relators':target,'exact_match':words==target})
        print(json.dumps({'name':old['name'],'model':kind,'epsilon':e,'exact_match':words==target,'final_relators':words,'target':target,'steps':len(steps)}),flush=True)
    result={'status':'COMPLETE','results':rows,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['three_sheet_splice_plan.json','audit_three_sheet_splices.py','audit_three_sheet_structure.py','three_sheet_structural_results.json','three_sheet_structure_verification.json']}}
    (HERE/'three_sheet_splice_results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()

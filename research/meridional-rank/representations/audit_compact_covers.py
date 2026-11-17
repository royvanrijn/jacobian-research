#!/usr/bin/env python3
"""Exact compact core presentations for two certified common-cover blocks."""
import hashlib
import json
from pathlib import Path
import time
from compact_cover_groups import free,simplify,graph_core
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic();plan=json.loads((HERE/'compact_cover_plan.json').read_text())
    mutation=json.loads((HERE/'mutation_cover_results.json').read_text())
    files=['compact_cover_plan.json','compact_cover_groups.py','audit_compact_covers.py','parabolic_lifts.py',
        'mutation_cover_results.json','mutation_cover_verification.json','core_generator_results.json']
    result={'status':'RUNNING','input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'controls':[],'targets':[]}
    for label,gens,rels in [
        ('C3*C3*C3',[1,2,3],[[1]*3,[2]*3,[3]*3]),
        ('C3*C3',[1,2],[[1]*3,[2]*3]),
        ('rebased_C3*C3_extra_generator',[1,2,3],[[1,2]*3,[2]*3,[3,-2,-1]])]:
        row={'name':label,**simplify(gens,rels,plan['limits'])}
        assert len(row['final']['generators'])==(3 if label=='C3*C3*C3' else 2)
        result['controls'].append(row)
    path=HERE/'compact_cover_results.json'
    for name in plan['targets']:
        assert time.monotonic()-started<plan['limits']['seconds_total']
        data=mutation['diagrams'][name];pd=data['pd'];crossings=oriented_diagram(pd)
        rels=[free([x for x in [o+1,-(a+1),o+1,-(b+1)] if abs(x)!=1]) for o,a,b,_ in crossings]
        row={'name':name,'pd':pd,'base_arc':0,'crossings':[list(r) for r in crossings],
            'raw_relators':rels,'graph_core':graph_core(data['selected']),
            **simplify(list(range(2,len(pd)+1)),rels,plan['limits'])}
        result['targets'].append(row);path.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'name':name,'status':row['status'],'generators':row['final']['generators'],
            'relator_lengths':list(map(len,row['final']['relators'])),'steps':len(row['trace']),'seconds':row['seconds']}),flush=True)
    result['status']='COMPLETE';result['seconds']=round(time.monotonic()-started,6)
    result['assurance']='Every accepted presentation move is exact; independent replay required. A local minimum/cap proves no rank lower bound. Abstract graph K4 reduction does not imply a geometric splitting.'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'seconds':result['seconds']}))


if __name__=='__main__':main()

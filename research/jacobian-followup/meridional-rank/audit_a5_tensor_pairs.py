#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import time
from audit_a5_all_modules import small_cocycles
from quadratic_integer_lifts import residue_fields

HERE=Path(__file__).resolve().parent


def tensor(a,b,field):
    return [field.mul(a[2*i+j],b[2*k+l]) for i in range(2) for k in range(2) for j in range(2) for l in range(2)]


def main():
    plan=json.loads((HERE/'a5_tensor_pair_plan.json').read_text());assert plan['frozen_before_execution']
    source=json.loads((HERE/'a5_all_modules_results.json').read_text());assert source['status']=='COMPLETE'
    grouped={}
    for row in source['results']:grouped.setdefault(row['name'],[]).append(row)
    selected={name:rows for name,rows in grouped.items() if len(rows)==2};assert len(selected)==4
    field=residue_fields(2)[0];results=[];start=time.monotonic()
    for name,rows in selected.items():
        assert rows[0]['words']==rows[1]['words']
        for i in range(2):
            for j in range(2):
                assert time.monotonic()-start<30
                a,b=rows[0]['binary_simple_modules'][i],rows[1]['binary_simple_modules'][j]
                images=[tensor(x,y,field) for x,y in zip(a['images'],b['images'])]
                result=small_cocycles(images,rows[0]['words'],4,field)
                result.update({'name':name,'orbit_indices':[r['orbit_index'] for r in rows],
                               'frobenius_choices':[i,j]})
                results.append(result)
    paths=['audit_a5_tensor_pairs.py','a5_tensor_pair_plan.json','a5_all_modules_results.json','audit_a5_all_modules.py']
    report={'status':'COMPLETE','results':results,
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'a5_tensor_pair_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{'name':r['name'],'twists':r['frobenius_choices'],'Z1':r['Z1_dimension'],'H1':r['H1_dimension'],'sharp_four':r['sharp_four']} for r in results],indent=2))


if __name__=='__main__':main()

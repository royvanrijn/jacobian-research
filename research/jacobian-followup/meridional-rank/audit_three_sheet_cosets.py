#!/usr/bin/env python3
"""Bounded coset actions, accepted only as directly evaluated quotient witnesses."""
import hashlib
import json
from pathlib import Path
import signal
import sympy
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup,coset_enumeration_r

HERE=Path(__file__).resolve().parent


def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):return tuple(a.index(i) for i in range(len(a)))
def evaluate(word,assignment):
    identity=tuple(range(len(next(iter(assignment.values())))));value=identity
    for a in word:value=mul(value,assignment[a] if a>0 else inverse(assignment[-a]))
    return value
def timeout(*_):raise TimeoutError('frozen case wall-time cap')


def main():
    assert sympy.__version__=='1.14.0'
    plan=json.loads((HERE/'three_sheet_coset_plan.json').read_text())
    source=json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())
    F,x,y=free_group('x,y');basis=[x,y];rows=[];signal.signal(signal.SIGALRM,timeout)
    for row in source['results']:
        if len(row['remaining_generators'])!=2:continue
        names=row['remaining_generators'];dictionary=dict(zip(names,basis));relators=[]
        for word in row['remaining_relators']:
            value=F.identity
            for a in word:value*=dictionary[abs(a)]**(1 if a>0 else -1)
            relators.append(value)
        record={'name':row['name']};signal.alarm(plan['limits']['seconds_per_case'])
        try:
            group=FpGroup(F,relators)
            table=coset_enumeration_r(group,[],max_cosets=plan['limits']['maximum_cosets'])
            assert table.is_complete();table.compress();table.standardize()
            assert list(table.A)==[x,x**-1,y,y**-1]
            assignment={names[j]:tuple(r[2*j+1] for r in table.table) for j in range(2)}
            degree=len(table.table);identity=tuple(range(degree))
            assert all(evaluate(w,assignment)==identity for w in row['remaining_relators'])
            nonabelian=mul(assignment[names[0]],assignment[names[1]])!=mul(assignment[names[1]],assignment[names[0]])
            for step in reversed(row['eliminations']):assignment[step['generator']]=evaluate(step['replacement'],assignment)
            assert all(evaluate(w,assignment)==identity for w in row['relators'])
            record.update({'status':'NONABELIAN_WITNESS' if nonabelian else 'CYCLIC_IMAGE_ONLY','degree':degree,
                           'original_generator_permutations':assignment,'coset_table':table.table,
                           'scope':'An actual source-group cyclicity/order claim is not accepted from this action alone.'})
        except (ValueError,TimeoutError) as error:record.update({'status':'CAPPED_INCOMPLETE','reason':str(error)})
        finally:signal.alarm(0)
        rows.append(record)
        report={'status':'IN_PROGRESS','sympy_version':sympy.__version__,'results':rows,
                'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['three_sheet_coset_plan.json','audit_three_sheet_cosets.py','three_sheet_nonabelian_results.json']}}
        (HERE/'three_sheet_coset_results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:record[k] for k in ['name','status','degree'] if k in record}),flush=True)
    report['status']='COMPLETE';(HERE/'three_sheet_coset_results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

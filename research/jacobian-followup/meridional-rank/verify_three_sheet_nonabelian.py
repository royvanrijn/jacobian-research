#!/usr/bin/env python3
"""Replay folded branched-cover words, exact Tietze moves and finite tests."""
from collections import Counter
from functools import reduce
from itertools import combinations,product
from math import gcd
import hashlib
import json
from pathlib import Path
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def inverse(w):return [-x for x in reversed(w)]
def free(w):
    out=[]
    for x in w:
        if out and out[-1]==-x:out.pop()
        else:out.append(x)
    return out
def cyclic(w):
    out=free(w)
    while len(out)>1 and out[0]==-out[-1]:out=out[1:-1]
    return out
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))


def group(name):
    if name=='S3':gens=[(1,0,2),(1,2,0)]
    elif name=='A4':gens=[(1,2,0,3),(1,0,3,2)]
    elif name=='A5':gens=[(1,2,0,3,4),(1,2,3,4,0)]
    else:
        n=int(name[1:])//2;gens=[tuple((i+1)%n for i in range(n)),tuple(-i%n for i in range(n))]
    identity=tuple(range(len(gens[0])));elements=[identity];seen={identity}
    for a in elements:
        for b in gens:
            c=mul(a,b)
            if c not in seen:seen.add(c);elements.append(c)
    expected={'S3':6,'A4':12,'D10':10,'D14':14,'A5':60}[name];assert len(elements)==expected
    table=[[elements.index(mul(a,b)) for b in elements] for a in elements]
    inv=[next(j for j in range(expected) if table[i][j]==0) for i in range(expected)]
    assert any(table[a][b]!=table[b][a] for a in range(expected) for b in range(expected))
    return table,inv


def main():
    report=json.loads((HERE/'three_sheet_nonabelian_results.json').read_text());assert report['status']=='COMPLETE'
    for p,h in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    models={name:group(name) for name in ['S3','A4','D10','D14','A5']};summaries=[];cyclic_orders=Counter()
    for row in report['results']:
        oriented=oriented_diagram(panel[row['name']]['pd']);labels=row['arc_labels'];n=len(labels)
        edges=list(combinations(range(3),2));tree=row['tree_generators'];seen={0}
        assert len(tree)==2 and len(set(tree))==2
        for _ in range(2):
            for i in tree:
                a,b=edges[labels[i-1]]
                if a in seen or b in seen:seen.update([a,b])
        assert seen=={0,1,2}
        transports=[]
        for i,label in enumerate(labels,1):
            a,b=edges[label];perm=list(range(3));perm[a],perm[b]=b,a
            words=[[] for _ in range(3)]
            if i not in tree:words[a]=[i];words[b]=[-i]
            transports.append((perm,words))
        words=[]
        for o,a,b,sign in oriented:
            relator=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            for start in range(3):
                sheet=start;w=[]
                for letter in relator:
                    permutation,steps=transports[abs(letter)-1];target=permutation[sheet]
                    w.extend(steps[sheet] if letter>0 else inverse(steps[target]));sheet=target
                assert sheet==start;words.append(w)
        assert words==row['relators']
        generators=set(row['generators']);assert generators==set(range(1,n+1))-set(tree)
        words=[cyclic(w) for w in words if w]
        for step in row['eliminations']:
            j=step['relator_index'];w=words[j];g=step['generator'];replacement=step['replacement']
            assert w==step['relator'] and sum(abs(a)==g for a in w)==1 and g in generators
            substituted=[]
            for a in w:substituted.extend(replacement if a==g else inverse(replacement) if a==-g else [a])
            assert free(substituted)==[] and all(abs(a) in generators-{g} for a in replacement)
            words.pop(j);updated=[]
            for w in words:
                value=[]
                for a in w:value.extend(replacement if a==g else inverse(replacement) if a==-g else [a])
                value=cyclic(value)
                if value:updated.append(value)
            words=updated;generators.remove(g)
        assert words==row['remaining_relators'] and generators==set(row['remaining_generators'])
        if len(generators)<=1:
            order=1 if not generators else reduce(gcd,(abs(sum(1 if a>0 else -1 for a in w)) for w in words),0)
            assert row['cyclic_group_order']==order;cyclic_orders[order]+=1
            summaries.append({'name':row['name'],'cyclic_order':order})
        else:
            generators=sorted(generators);assert len(generators)==2
            for test in row['finite_tests']:
                assert test['status']=='EXHAUSTED';table,inv=models[test['group']];order=len(table);nonabelian=0
                for values in product(range(order),repeat=2):
                    if table[values[0]][values[1]]==table[values[1]][values[0]]:continue
                    assignment=dict(zip(generators,values));valid=True
                    for w in words:
                        value=0
                        for a in w:value=table[value][assignment[a] if a>0 else inv[assignment[-a]]]
                        if value:valid=False;break
                    nonabelian+=valid
                assert nonabelian==0 and test['tested']==order**2
            summaries.append({'name':row['name'],'remaining_generators':2,'nonabelianity':'UNRESOLVED','tested_nonabelian_images':0})
    assert len(summaries)==26 and sum(cyclic_orders.values())==19
    cosets=json.loads((HERE/'three_sheet_coset_results.json').read_text());assert cosets['status']=='COMPLETE'
    for p,h in cosets['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    originals={r['name']:r for r in report['results']};coset_status=Counter()
    for row in cosets['results']:
        coset_status[row['status']]+=1
        if row['status']=='CAPPED_INCOMPLETE':continue
        assert row['status']=='CYCLIC_IMAGE_ONLY';degree=row['degree']
        images={int(k):tuple(v) for k,v in row['original_generator_permutations'].items()}
        assert all(sorted(p)==list(range(degree)) for p in images.values())
        for word in originals[row['name']]['relators']:
            result=tuple(range(degree))
            for letter in word:
                p=images[abs(letter)]
                if letter<0:p=tuple(p.index(i) for i in range(degree))
                result=mul(result,p)
            assert result==tuple(range(degree))
        assert all(mul(a,b)==mul(b,a) for a in images.values() for b in images.values())
    assert coset_status=={'CAPPED_INCOMPLETE':4,'CYCLIC_IMAGE_ONLY':3}
    paths=['verify_three_sheet_nonabelian.py','three_sheet_nonabelian_plan.json','three_sheet_nonabelian_results.json',
           'rank_two_equivalence_verification.json','joint_a5_panel.json','parabolic_lifts.py','three_sheet_coset_results.json']
    receipt={'status':'PASS','covers':26,'cyclic_orders':dict(cyclic_orders),'remaining_two_generator_groups':7,
             'finite_nonabelian_targets':['S3','A4','D10','D14','A5'],'new_sharp_cases':0,'results':summaries,
             'coset_followup':dict(coset_status),'coset_boundary':'Three returned cyclic actions are replayed as quotients only; source cyclicity is not inferred. Four enumerations hit the frozen cap; all seven source ranks remain unresolved.',
             'scope':'Independent signed edge-transport reconstruction and exact Tietze substitution replay; finite targets independently generated from permutations and all noncommuting pairs evaluated. PD parser and retained colouring inputs shared. Seven actual group ranks/cyclicity unresolved; no inferred three-bridge presentation.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'three_sheet_nonabelian_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','covers','cyclic_orders','remaining_two_generator_groups','new_sharp_cases']},indent=2))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact presentations and a frozen nonabelianity gate for all26 S3 covers."""
from itertools import permutations,product
from math import gcd
import hashlib
import json
from pathlib import Path
import time
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def reduce_word(word):
    stack=[]
    for a in word:
        if stack and stack[-1]==-a:stack.pop()
        else:stack.append(a)
    while len(stack)>1 and stack[0]==-stack[-1]:stack=stack[1:-1]
    return stack


def inverse(word):return [-a for a in reversed(word)]


def eliminate(generators,relators):
    current=generators[:];words=[reduce_word(w) for w in relators if w];steps=[]
    while True:
        choice=next(((j,abs(a),a) for j,w in enumerate(words) for a in w if sum(abs(b)==abs(a) for b in w)==1),None)
        if choice is None:break
        j,g,a=choice;word=words[j];i=word.index(a);rest=word[i+1:]+word[:i]
        replacement=inverse(rest) if a>0 else rest
        steps.append({'generator':g,'relator_index':j,'relator':word,'replacement':replacement})
        words.pop(j);new=[]
        for word in words:
            expanded=[]
            for letter in word:expanded.extend(replacement if letter==g else inverse(replacement) if letter==-g else [letter])
            value=reduce_word(expanded)
            if value:new.append(value)
        words=new;current.remove(g)
    return current,words,steps


def presentation(crossings,labels):
    edges=[(0,1),(0,2),(1,2)];parent=list(range(3));tree=[]
    def root(x):
        while parent[x]!=x:x=parent[x]
        return x
    for i,label in enumerate(labels):
        a,b=edges[label]
        if root(a)!=root(b):parent[root(a)]=root(b);tree.append(i+1)
    assert len(tree)==2
    words=[]
    for o,a,b in crossings:
        for start in range(3):
            sheet=start;word=[]
            for i in [o,a,o,b]:
                p,q=edges[labels[i]]
                if sheet in (p,q):
                    if i+1 not in tree:word.append(i+1 if sheet==p else -(i+1))
                    sheet=q if sheet==p else p
            assert sheet==start
            words.append(word)
    generators=[i+1 for i in range(len(labels)) if i+1 not in tree]
    return generators,words,tree


def finite_group(name):
    if name in ['S3','A4','A5']:
        n=3 if name=='S3' else int(name[1:]);elements=list(permutations(range(n)))
        if name!='S3':elements=[p for p in elements if sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2==0]
    else:
        n=int(name[1:])//2;elements=[tuple((sign*i+a)%n for i in range(n)) for sign in [1,-1] for a in range(n)]
    identity=tuple(range(len(elements[0])));elements.remove(identity);elements.insert(0,identity)
    index={p:i for i,p in enumerate(elements)}
    table=[[index[tuple(a[b[i]] for i in range(len(a))) ] for b in elements] for a in elements]
    inv=[next(j for j in range(len(elements)) if table[i][j]==0) for i in range(len(elements))]
    return elements,table,inv


def main():
    plan=json.loads((HERE/'three_sheet_nonabelian_plan.json').read_text());panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases']
    colours=json.loads((HERE/'small_symmetric_full_panel_results.json').read_text())['results']
    inputs=['three_sheet_nonabelian_plan.json','audit_three_sheet_nonabelian.py','joint_a5_panel.json','small_symmetric_full_panel_results.json','parabolic_lifts.py']
    report={'status':'IN_PROGRESS','results':[],'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    path=HERE/'three_sheet_nonabelian_results.json';started=time.monotonic();models={name:finite_group(name) for name in plan['targets']}
    for case,row in zip(panel,colours):
        assert case['name']==row['name'];cr=[c[:3] for c in oriented_diagram(case['pd'])]
        for cert in row['degrees'][1]['certificates']:
            assert time.monotonic()-started<plan['limits']['seconds_total'];tick=time.monotonic()
            gens,rels,tree=presentation(cr,cert['arc_labels']);remaining,words,steps=eliminate(gens,rels)
            record={'name':case['name'],'arc_labels':cert['arc_labels'],'tree_generators':tree,'generators':gens,'relators':rels,
                    'eliminations':steps,'remaining_generators':remaining,'remaining_relators':words,'finite_tests':[]}
            if len(remaining)<=1:
                order=0
                for w in words:order=gcd(order,abs(sum(1 if a>0 else -1 for a in w)))
                record['cyclic_group_order']=1 if not remaining else order
            else:
                for name,(elements,table,inv) in models.items():
                    cap=6 if name=='S3' else 3 if name=='A5' else 4
                    if len(remaining)>cap:
                        record['finite_tests'].append({'group':name,'status':'SKIPPED_GENERATOR_CAP'});continue
                    tested=0;witness=None
                    for values in product(range(len(elements)),repeat=len(remaining)):
                        tested+=1
                        if time.monotonic()-tick>plan['limits']['seconds_per_target']:break
                        if all(table[a][b]==table[b][a] for a in values for b in values):continue
                        assignment=dict(zip(remaining,values));ok=True
                        for word in words:
                            x=0
                            for a in word:x=table[x][assignment[a] if a>0 else inv[assignment[-a]]]
                            if x:ok=False;break
                        if not ok:continue
                        for step in reversed(steps):
                            x=0
                            for a in step['replacement']:x=table[x][assignment[a] if a>0 else inv[assignment[-a]]]
                            assignment[step['generator']]=x
                        witness={'original_generator_images':assignment,'permutations':elements,'remaining_images':values};break
                    record['finite_tests'].append({'group':name,'status':'NONABELIAN_WITNESS' if witness else 'EXHAUSTED' if tested==len(elements)**len(remaining) else 'CAPPED_INCOMPLETE','tested':tested,'witness':witness})
                    if witness:record['sharp_four']=True;break
            report['results'].append(record);path.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({'name':case['name'],'generators':len(remaining),'cyclic_order':record.get('cyclic_group_order'),'sharp':record.get('sharp_four',False)}),flush=True)
    report['status']='COMPLETE';report['new_sharp_names']=[r['name'] for r in report['results'] if r.get('sharp_four')]
    path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':'COMPLETE','covers':len(report['results']),'sharp':report['new_sharp_names']}))


if __name__=='__main__':main()

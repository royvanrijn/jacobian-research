#!/usr/bin/env python3
"""Exact length-decreasing presentation rewrites, with full replay trace."""
import hashlib
import json
from math import gcd
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def inv(w):return tuple(-a for a in reversed(w))
def free(w):
    out=[]
    for a in w:
        if out and out[-1]==-a:out.pop()
        else:out.append(a)
    return tuple(out)
def cyclic(w):
    w=free(w)
    while len(w)>1 and w[0]==-w[-1]:w=w[1:-1]
    return w
def rotations(w):return [w[i:]+w[:i] for i in range(len(w))]
def canonical(w):
    w=cyclic(w)
    return min(rotations(w)+rotations(inv(w))) if w else ()
def normalize(words):return sorted({canonical(w) for w in words}-{()})
def subst(w,mapping):return free(a for x in w for a in (mapping[x] if x>0 else inv(mapping[-x])))
def measure(words):return sum(map(len,words))


def shorten(words):
    best=None
    for i,r in enumerate(words):
        replacements={}
        for signed in [r,inv(r)]:
            for relation in rotations(signed):
                for k in range(len(r)//2+1,len(r)+1):
                    pattern,complement=relation[:k],relation[k:]
                    replacements.setdefault(pattern,(inv(complement),relation))
        for j,w in enumerate(words):
            if i==j:continue
            for rotated in rotations(w):
                for k in range(1,len(rotated)+1):
                    match=replacements.get(rotated[:k])
                    if match is None:continue
                    replacement,relation=match
                    new=normalize(words[:j]+words[j+1:]+[replacement+rotated[k:]])
                    key=(measure(new),new)
                    if best is None or key<best[0]:
                        best=(key,new,{'operation':'relator_shortening','used_relator':list(r),'used_signed_rotation':list(relation),
                                      'target_relator':list(w),'target_rotation':list(rotated),'pattern_length':k})
    return None if best is None else best[1:]


def nielsen(words,active,exponents):
    if len(active)!=2:return None
    a,b=active;best=None
    for x,y in [(a,b),(b,a)]:
        for n in exponents:
            power=(y,)*n if n>0 else (-y,)*(-n)
            reverse=inv(power)
            for side in ['left','right']:
                mapping={a:(a,),b:(b,)};inverse={a:(a,),b:(b,)}
                mapping[x]=power+(x,) if side=='left' else (x,)+power
                inverse[x]=reverse+(x,) if side=='left' else (x,)+reverse
                new=normalize([subst(w,mapping) for w in words])
                if measure(new)>=measure(words):continue
                key=(measure(new),new)
                if best is None or key<best[0]:
                    best=(key,new,{'operation':'nielsen','forward':mapping,'inverse':inverse})
    return None if best is None else best[1:]


def main():
    plan=json.loads((HERE/'three_sheet_structural_plan.json').read_text())
    source=json.loads((HERE/'three_sheet_nonabelian_results.json').read_text());rows=[];started=time.monotonic()
    for old in source['results']:
        if len(old['remaining_generators'])!=2:continue
        names=old['remaining_generators'];translate={names[0]:1,names[1]:2}
        words=normalize([[translate[abs(a)]*(1 if a>0 else -1) for a in w] for w in old['remaining_relators']])
        initial=[list(w) for w in words];active=[1,2];images={1:(1,),2:(2,)};steps=[]
        for _ in range(plan['limits']['steps_per_case']):
            assert time.monotonic()-started<plan['limits']['seconds_total']
            eliminate=next(((w,abs(a),a) for w in words for a in w if sum(abs(b)==abs(a) for b in w)==1),None)
            if eliminate is not None:
                w,g,a=eliminate;i=w.index(a);rest=w[i+1:]+w[:i];replacement=inv(rest) if a>0 else rest
                mapping={x:(x,) for x in active};mapping[g]=replacement
                new=normalize([subst(v,mapping) for v in words if v!=w]);step={'operation':'eliminate','relator':list(w),'generator':g,'replacement':list(replacement)}
                images={x:subst(v,mapping) for x,v in images.items()};active.remove(g)
            else:
                candidate=shorten(words)
                if candidate is None:candidate=nielsen(words,active,plan['limits']['nielsen_exponents'])
                if candidate is None:break
                new,step=candidate
                if step['operation']=='nielsen':images={x:subst(v,step['forward']) for x,v in images.items()}
            assert measure(new)<measure(words) or len(active)<2
            assert all(len(w)<=plan['limits']['maximum_intermediate_word_length'] for w in new)
            step['before_measure']=measure(words);step['after_measure']=measure(new);steps.append(step);words=new
        record={'name':old['name'],'original_generator_names':names,'initial_canonical_relators':initial,'steps':steps,
                'active_generators':active,'final_relators':words,'original_generator_words':images,'status':'EXACT_REWRITE_PROPOSAL'}
        if len(active)<=1:
            order=1 if not active else 0
            if active:
                for w in words:order=gcd(order,abs(sum(1 if a>0 else -1 for a in w)))
            record['cyclic_order']=order
        rows.append(record)
        report={'status':'IN_PROGRESS','results':rows,'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['three_sheet_structural_plan.json','audit_three_sheet_structure.py','three_sheet_nonabelian_results.json']}}
        (HERE/'three_sheet_structural_results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({'name':old['name'],'active':active,'final_relators':words,'cyclic_order':record.get('cyclic_order'),'steps':len(steps)}),flush=True)
    report['status']='COMPLETE';(HERE/'three_sheet_structural_results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Replay only exact presentation isomorphisms, without producer heuristics."""
from functools import reduce
from math import gcd
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def invert(word):return [-a for a in word[::-1]]
def evaluate(word,substitution):
    value=[]
    for a in word:
        for b in substitution[a] if a>0 else invert(substitution[-a]):
            if value and value[-1]==-b:value.pop()
            else:value.append(b)
    return value
def rotations(word):return [word[i:]+word[:i] for i in range(len(word))]
def cyclic(word):
    value=evaluate(word,{abs(a):[abs(a)] for a in word})
    while len(value)>1 and value[0]==-value[-1]:value=value[1:-1]
    return value
def canonical(word):
    value=cyclic(word)
    return min(rotations(value)+rotations(invert(value))) if value else []
def normalized(words):return [list(w) for w in sorted({tuple(canonical(w)) for w in words}-{()})]


def main():
    report=json.loads((HERE/'three_sheet_structural_results.json').read_text());assert report['status']=='COMPLETE'
    for p,h in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    previous=json.loads((HERE/'three_sheet_nonabelian_verification.json').read_text());assert previous['status']=='PASS'
    for p,h in previous['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    originals={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    rows=[]
    for row in report['results']:
        old=originals[row['name']];a,b=old['remaining_generators'];translate={a:1,b:2}
        words=normalized([[translate[abs(x)]*(1 if x>0 else -1) for x in w] for w in old['remaining_relators']])
        assert words==row['initial_canonical_relators'];active={1,2};images={1:[1],2:[2]}
        for step in row['steps']:
            before=sum(map(len,words));op=step['operation']
            if op=='relator_shortening':
                r=step['used_relator'];target=step['target_relator'];signed=step['used_signed_rotation'];rotated=step['target_rotation'];k=step['pattern_length']
                assert r in words and target in words and r!=target
                assert signed in rotations(r)+rotations(invert(r)) and rotated in rotations(target)
                assert k>len(r)//2 and k<=len(r) and rotated[:k]==signed[:k]
                replacement=invert(signed[k:]);new=words[:];new.remove(target);new.append(replacement+rotated[k:]);words=normalized(new)
            elif op=='nielsen':
                forward={int(k):v for k,v in step['forward'].items()};backward={int(k):v for k,v in step['inverse'].items()}
                assert set(forward)==set(backward)==active
                for g in active:assert evaluate(forward[g],backward)==evaluate(backward[g],forward)==[g]
                words=normalized([evaluate(w,forward) for w in words]);images={g:evaluate(w,forward) for g,w in images.items()}
            elif op=='eliminate':
                r=step['relator'];g=step['generator'];replacement=step['replacement']
                assert r in words and g in active and sum(abs(a)==g for a in r)==1
                substitution={a:[a] for a in active};substitution[g]=replacement
                assert evaluate(r,substitution)==[] and all(abs(a) in active-{g} for a in replacement)
                new=words[:];new.remove(r);words=normalized([evaluate(w,substitution) for w in new])
                images={a:evaluate(w,substitution) for a,w in images.items()};active.remove(g)
            else:raise AssertionError(op)
            assert before==step['before_measure'] and sum(map(len,words))==step['after_measure']
        assert words==row['final_relators'] and active==set(row['active_generators'])
        assert images=={int(k):v for k,v in row['original_generator_words'].items()}
        result={'name':row['name'],'final_relators':words,'original_generator_words':images,'verified_steps':len(row['steps'])}
        if len(active)==1:
            order=reduce(gcd,(abs(sum(1 if a>0 else -1 for a in w)) for w in words),0)
            assert order==row['cyclic_order'];result['proved_cyclic_order']=order
        else:
            assert len(words)==len(active)==2
            matrix=[[sum((1 if x>0 else -1) for x in w if abs(x)==g) for g in [1,2]] for w in words]
            determinant=abs(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0])
            assert reduce(gcd,(abs(x) for r in matrix for x in r),0)==1
            assert determinant in [35,37];result['abelianization_order']=determinant;result['exponent_matrix']=matrix
        rows.append(result)
    assert len(rows)==7 and [r.get('proved_cyclic_order') for r in rows[-3:]]==[2,2,4]
    paths=['verify_three_sheet_structure.py','three_sheet_structural_plan.json','three_sheet_structural_results.json','three_sheet_nonabelian_results.json','three_sheet_nonabelian_verification.json']
    receipt={'status':'PASS','presentations':7,'new_proved_cyclic':3,'total_proved_cyclic_covers':22,'remaining_groups':4,
             'results':rows,'scope':'Exact normal-closure-preserving relator substitutions, invertible Nielsen maps and explicit Tietze eliminations independently replayed. Three actual source groups now proved cyclic, beyond the earlier quotient images. Four groups remain unresolved; no knot3bridge presentation follows.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'three_sheet_structure_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','new_cyclic':3,'remaining':[(r['name'],r['abelianization_order']) for r in rows if 'abelianization_order' in r]}))


if __name__=='__main__':main()

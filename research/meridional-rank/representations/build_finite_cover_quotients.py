#!/usr/bin/env python3
"""Explicit degree25 finite witnesses from regular peripheral translations."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def compose(a,b):return [a[b[i]] for i in range(len(a))]
def inverse(a):return [a.index(i) for i in range(len(a))]
def power(a,n):
    if n<0:return power(inverse(a),-n)
    result=list(range(len(a)))
    for _ in range(n):result=compose(result,a)
    return result
def evaluate(word,assignment):
    result=list(range(len(next(iter(assignment.values())))))
    for g in word:result=compose(result,assignment[g] if g>0 else inverse(assignment[-g]))
    return result


def main():
    plan=json.loads((HERE/'finite_cover_quotient_plan.json').read_text())
    originals={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    structural={r['name']:r for r in json.loads((HERE/'three_sheet_structural_results.json').read_text())['results']}
    splices=json.loads((HERE/'three_sheet_splice_verification.json').read_text());assert splices['status']=='PASS exact sharp-four evidence'
    points=[(x,y) for y in range(5) for x in range(5)]
    p=[3,0,2,1,4];q=[4,1,3,2,0]
    alpha=[points.index((p[x],(y+2)%5)) for x,y in points]
    delta=[points.index((q[x],(y+3)%5)) for x,y in points]
    identity=list(range(25));fibre=power(delta,2);meridian=compose(delta,inverse(alpha))
    assert power(alpha,3)==fibre==[points.index((x,(y+1)%5)) for x,y in points]
    assert meridian==[points.index(((x+1)%5,(y+1)%5)) for x,y in points]
    cases=[];full=[]
    for row in splices['results']:
        name=row['name'];kind=row['model'];epsilon=row['epsilon']
        linear=[points.index(((x-y)%5,((1+epsilon)*x-y)%5)) for x,y in points]
        if kind=='B':linear=[points.index((-points[i][0]%5,-points[i][1]%5)) for i in linear]
        model_pair={1:alpha,2:compose(compose(linear,alpha if kind=='A' else delta),inverse(linear))}
        model=next(r for r in json.loads((HERE/'three_sheet_splice_results.json').read_text())['results'] if r['name']==name)
        assert all(evaluate(w,model_pair)==identity for w in model['target_relators'])
        # The stored maps send each previous source generator into the new basis.
        intermediate={g:evaluate(w,model_pair) for g,w in ((int(k),v) for k,v in row['splice_basis_transport'].items())}
        source_pair={g:evaluate(w,intermediate) for g,w in ((int(k),v) for k,v in structural[name]['original_generator_words'].items())}
        old=originals[name];names=old['remaining_generators'];rename={names[0]:1,names[1]:2}
        words=[[rename[abs(a)]*(1 if a>0 else -1) for a in w] for w in old['remaining_relators']]
        assert all(evaluate(w,source_pair)==identity for w in words)
        commutator=evaluate([1,2,-1,-2],source_pair);moved=next(i for i,v in enumerate(commutator) if i!=v)
        cases.append({'name':name,'degree':25,'generator_names':names,'relators':words,
                      'A':source_pair[1],'B':source_pair[2],
                      'commutator_witness':{'point':moved,'image':commutator[moved]}})
        assignment={names[0]:source_pair[1],names[1]:source_pair[2]}
        for step in reversed(old['eliminations']):assignment[step['generator']]=evaluate(step['replacement'],assignment)
        assert all(evaluate(w,assignment)==identity for w in old['relators'])
        full.append({'name':name,'original_cover_generator_images':assignment,
                     'original_cover_relators':old['relators'],'noncommuting_original_generators':names,
                     'model':kind,'epsilon':epsilon,'model_A':model_pair[1],'model_B':model_pair[2],
                     'linear_gluing_permutation':linear})
    certificate={'schema':1,'permutation_convention':'zero-based image arrays; compose(P,Q)[i]=P[Q[i]]',
                 'commutator_word':[1,2,-1,-2],'cases':cases}
    (HERE/'finite_cover_certificates.json').write_text(json.dumps(certificate,indent=2)+'\n')
    inputs=['finite_cover_quotient_plan.json','build_finite_cover_quotients.py','three_sheet_nonabelian_results.json',
            'three_sheet_structural_results.json','three_sheet_splice_results.json','three_sheet_splice_verification.json','finite_cover_certificates.json']
    report={'status':'CONSTRUCTED','degree':25,'cases':full,'trefoil_vertex_permutations':{'alpha':alpha,'delta':delta,'meridian':meridian,'fibre':fibre},
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/'finite_cover_full_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'CONSTRUCTED','certificates':len(cases),'degree':25,
                      'commutator_moved_points':[(r['name'],r['commutator_witness']) for r in cases]},indent=2))


if __name__=='__main__':main()

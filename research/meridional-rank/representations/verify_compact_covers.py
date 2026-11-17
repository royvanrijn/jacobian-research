#!/usr/bin/env python3
"""Replay exact presentation changes without importing the simplifier."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations,product
import hashlib
import json
from math import gcd
from pathlib import Path
from verify_double_transpositions import geometric_crossings

HERE=Path(__file__).resolve().parent


def reduce_word(w):
    result=[]
    for letter in w:
        if result[-1:]==[-letter]:result.pop()
        else:result.append(letter)
    return result


def backwards(w):return list(map(lambda x:-x,reversed(w)))


def relator_key(w):
    word=reduce_word(w)
    while word and word[0]==-word[-1]:word=word[1:-1]
    if not word:return ()
    return min(tuple(orientation[k:]+orientation[:k]) for orientation in (word,backwards(word)) for k in range(len(word)))


def clean(words):return [list(w) for w in sorted(set(map(relator_key,words))-{()})]


def replace(w,mapping):
    result=[]
    for letter in w:result=reduce_word(result+(mapping[abs(letter)] if letter>0 else backwards(mapping[-letter])))
    return result


def keys(d):return {int(g):w for g,w in d.items()}


def replay(row):
    gens=row['initial']['generators'][:];rels=clean(row['initial']['relators'])
    assert row['initial']['relators']==rels and len(gens)==len(set(gens))
    forward={g:[g] for g in gens};back={g:[g] for g in gens}
    for step in row['trace']:
        kind=step['kind'];old_cost=sum(map(len,rels))
        if kind=='eliminate':
            g=step['generator'];r=rels[step['relator']];assert g in gens
            positions=[i for i,x in enumerate(r) if abs(x)==g];assert len(positions)==1
            k=positions[0];rest=r[k+1:]+r[:k]
            word=backwards(rest) if r[k]>0 else rest
            assert word==step['solved_word'] and all(abs(x) in set(gens)-{g} for x in word)
            mapping={h:[h] for h in gens};mapping[g]=word
            assert not replace(r,mapping)
            rels=clean([replace(w,mapping) for i,w in enumerate(rels) if i!=step['relator']])
            forward={h:replace(w,mapping) for h,w in forward.items()};gens.remove(g);back.pop(g)
        elif kind=='automorphism':
            mapping=keys(step['images']);inverse=keys(step['inverse_images'])
            assert set(mapping)==set(inverse)==set(gens)
            assert all(abs(x) in gens for w in list(mapping.values())+list(inverse.values()) for x in w)
            for g in gens:
                assert replace(mapping[g],inverse)==[g] and replace(inverse[g],mapping)==[g]
            rels=clean([replace(w,mapping) for w in rels]);assert sum(map(len,rels))<old_cost
            forward={h:replace(w,mapping) for h,w in forward.items()}
            back={h:replace(inverse[h],back) for h in gens}
        elif kind=='relator_dehn':
            i,j=step['target'],step['source'];assert i!=j
            S=rels[j][:] if step['source_sign']==1 else backwards(rels[j]);assert step['source_sign'] in (1,-1)
            sr,tr=step['source_rotation'],step['target_rotation'];assert 0<=sr<len(S) and 0<=tr<len(rels[i])
            S=S[sr:]+S[:sr];T=rels[i][tr:]+rels[i][:tr];k=step['matched']
            assert len(S)//2<k<=min(len(S),len(T)) and S[:k]==T[:k]
            shortened=backwards(S[k:])+T[k:]
            # The replacement is precisely S^-1*T after free cancellation.
            assert reduce_word(backwards(S)+T)==reduce_word(shortened)
            rels=clean([w if n!=i else shortened for n,w in enumerate(rels)])
            assert sum(map(len,rels))<old_cost
        else:raise AssertionError('unknown presentation operation')
    assert row['final']=={'generators':gens,'relators':rels}
    assert keys(row['original_generators_in_final'])==forward
    assert keys(row['final_generators_in_original'])==back
    return gens,rels


def determinant(m):
    if len(m)==1:return m[0][0]
    return sum((-1)**j*m[0][j]*determinant([r[:j]+r[j+1:] for r in m[1:]]) for j in range(len(m)))


def abelianization(gens,rels):
    matrix=[[sum((1 if x>0 else -1) for x in w if abs(x)==g) for g in gens] for w in rels]
    n=len(gens);assert len(rels)==n
    divisors=[]
    for k in range(1,n+1):
        value=0
        for rows in combinations(range(n),k):
            for cols in combinations(range(n),k):value=gcd(value,abs(determinant([[matrix[i][j] for j in cols] for i in rows])))
        divisors.append(value)
    diagonal=[divisors[0]]+[divisors[i]//divisors[i-1] for i in range(1,n)]
    assert all(v>0 for v in diagonal)
    return diagonal


def graph_replay(graph,certificate):
    edges={i:list(e) for i,e in enumerate(graph['edges'])};vertices=set(range(graph['vertices']))
    # T and F count trees and separating spanning forests of each2terminal network.
    networks={i:{'T':1,'F':1,'crossings':[i]} for i in edges}
    for step in certificate['trace']:
        ids=step['edges'];assert len(ids)==len(set(ids)) and all(i in edges for i in ids)
        children=[networks[i] for i in ids]
        if step['kind']=='parallel':
            assert len(ids)>=2 and all(sorted(edges[i])==sorted(step['ends']) for i in ids)
            T,F=children[0]['T'],children[0]['F']
            for child in children[1:]:T,F=T*child['F']+child['T']*F,F*child['F']
        else:
            assert step['kind']=='series' and len(ids)==2
            v=step['vertex'];assert v in vertices and {i for i,e in edges.items() if v in e}==set(ids)
            endpoints=[next(x for x in edges[i] if x!=v) for i in ids]
            assert len(set(endpoints))==2 and sorted(endpoints)==sorted(step['ends'])
            A,B=children;T=A['T']*B['T'];F=A['T']*B['F']+B['T']*A['F'];vertices.remove(v)
        new=step['new_edge'];assert new not in edges
        leaves=sorted(x for child in children for x in child['crossings']);assert len(leaves)==len(set(leaves))
        for i in ids:edges.pop(i);networks.pop(i)
        edges[new]=step['ends'];networks[new]={'T':T,'F':F,'crossings':leaves}
    assert sorted(vertices)==certificate['vertices'] and edges=={int(i):e for i,e in certificate['edges'].items()}
    assert len(vertices)==4 and {tuple(sorted(e)) for e in edges.values()}==set(combinations(sorted(vertices),2))
    assert sorted(x for network in networks.values() for x in network['crossings'])==list(range(13))
    variable=next(i for i,d in networks.items() if len(d['crossings'])==7)
    double=next(i for i,d in networks.items() if len(d['crossings'])==2)
    assert len(set(edges[variable])&set(edges[double]))==1
    assert all(d['T']==d['F']==1 for i,d in networks.items() if i not in (variable,double))
    assert (networks[double]['T'],networks[double]['F'])==(2,1)
    coefficients=[0,0]
    for chosen in combinations(edges,3):
        parts=[{v} for v in vertices]
        for i in chosen:
            a,b=edges[i];A=next(s for s in parts if a in s);B=next(s for s in parts if b in s)
            if A is B:break
            parts.remove(A);parts.remove(B);parts.append(A|B)
        else:
            assert len(parts)==1
            value=1
            for i,d in networks.items():
                if i!=variable:value*=d['T'] if i in chosen else d['F']
            coefficients[0 if variable in chosen else 1]+=value
    assert coefficients==[11,13]
    data=networks[variable]
    return {'tetrahedral_vertices':sorted(vertices),'variable_branch_T_F':[data['T'],data['F']],
        'variable_branch_crossings':data['crossings'],'tree_count_coefficients':coefficients,
        'tree_count':11*data['T']+13*data['F'],
        'scope':'Exact abstract two-terminal graph network; T/F is conductance, not a certified rational-tangle slope.'}


def quotient(row,source):
    # Independently act on projective lines; no producer matrix arithmetic.
    vectors=list(product(range(3),repeat=3))
    canon=lambda v:min(v,tuple(-x%3 for x in v))
    lines=sorted({canon(v) for v in vectors if any(v)})
    matrices=next(m['matrices'] for m in source['modules'] if m['module']=='natural3')
    permutations={i+1:tuple(lines.index(canon(tuple(sum(M[3*a+b]*v[b] for b in range(3))%3 for a in range(3)))) for v in lines) for i,M in enumerate(matrices)}
    identity=tuple(range(13));assert permutations[1]==identity
    def evaluate(w,images):
        result=identity
        for x in reversed(w):
            p=images[abs(x)]
            if x<0:p=tuple(p.index(i) for i in range(13))
            result=tuple(p[result[i]] for i in range(13))
        return result
    assert all(evaluate(w,permutations)==identity for w in row['raw_relators'])
    images={int(g):evaluate(w,permutations) for g,w in row['final_generators_in_original'].items()}
    assert all(evaluate(w,images)==identity for w in row['final']['relators'])
    assert all(evaluate(w,images)==permutations[int(g)] for g,w in row['original_generators_in_final'].items())
    seen={identity};queue=[identity]
    for p in queue:
        for q in images.values():
            r=tuple(q[p[i]] for i in range(13))
            if r not in seen:seen.add(r);queue.append(r)
        assert len(seen)<=5616
    assert len(seen)==5616
    return {'order':len(seen),'generator_permutations':images,'all_original_and_final_relators_pass':True}


def main():
    result=json.loads((HERE/'compact_cover_results.json').read_text())
    mutation=json.loads((HERE/'mutation_cover_results.json').read_text())
    source=json.loads((HERE/'ternary_cover_module_results.json').read_text())['results']
    assert result['status']=='COMPLETE' and [r['name'] for r in result['targets']]==['13a_129','13a_135']
    for path,digest in result['input_sha256'].items():assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest,path
    controls=[]
    for row in result['controls']:
        gens,rels=replay(row);diagonal=abelianization(gens,rels)
        assert diagonal==[3]*len(gens)
        controls.append({'name':row['name'],'generators':len(gens),'F3_betti':len(diagonal),'steps':len(row['trace'])})
    targets=[]
    for row,D in zip(result['targets'],[391,401]):
        original=mutation['diagrams'][row['name']];assert row['pd']==original['pd'] and row['base_arc']==0
        triples=geometric_crossings(row['pd']);rels=[]
        assert [list(t) for t in triples]==[r[:3] for r in row['crossings']]
        for o,a,b in triples:rels.append(reduce_word([x for x in (o+1,-a-1,o+1,-b-1) if abs(x)!=1]))
        assert rels==row['raw_relators'] and clean(rels)==row['initial']['relators']
        assert row['initial']['generators']==list(range(2,14))
        gens,rels=replay(row);assert len(gens)==3 and len(rels)==3
        diagonal=abelianization(gens,rels);assert diagonal==[1,1,D]
        graph=graph_replay(original['selected'],row['graph_core']);assert graph['tree_count']==D
        quotient_data=quotient(row,next(s for s in source if s['name']==row['name']))
        targets.append({'name':row['name'],'presentation':row['final'],'trace_steps':dict(Counter(s['kind'] for s in row['trace'])),
            'abelianization_smith':diagonal,'nonabelian_quotient':quotient_data,'graph_network':graph,'actual_rank_interval':[2,3]})
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    bad=deepcopy(result['targets'][0]);step=next(s for s in bad['trace'] if s['kind']=='eliminate');step['solved_word']+=[step['generator']]
    reject('wrong solved Tietze word',lambda:replay(bad))
    bad=deepcopy(result['targets'][0]);step=next(s for s in bad['trace'] if s['kind']=='automorphism');key=next(iter(step['inverse_images']));step['inverse_images'][key]=[]
    reject('noninverse alleged free automorphism',lambda:replay(bad))
    bad=deepcopy(result['targets'][0]);step=next(s for s in bad['trace'] if s['kind']=='relator_dehn');step['matched']=0
    reject('invalid relator subword replacement',lambda:replay(bad))
    bad=deepcopy(result['targets'][0]);bad['final']['relators'][0][0]*=-1
    reject('altered compact final relator',lambda:replay(bad))
    bad=deepcopy(result['targets'][0]['graph_core']);bad['trace'][0]['ends']=[0,1]
    reject('wrong parallel network boundary',lambda:graph_replay(mutation['diagrams']['13a_129']['selected'],bad))
    files=['verify_compact_covers.py','compact_cover_results.json','compact_cover_plan.json','mutation_cover_results.json',
        'mutation_cover_verification.json','ternary_cover_module_results.json','ternary_cover_module_verification.json',
        'verify_double_transpositions.py','parabolic_lifts.py']
    receipt={'status':'PASS independent compact common-cover presentations','targets':targets,'controls':controls,
        'corruptions_rejected':rejected,'rank3_proved':False,'new_target_closures':0,
        'assurance':'Every accepted Tietze, inverse free automorphism and relator replacement is replayed without simplifier imports. Geometric arcs checked via legacy parser only for label transport, disclosed. Core-group identification is the retained classical STW input; mutation transfers only ordinary cover rank. Graph conductances do not certify rational-tangle slopes or JSJ data.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'compact_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'targets':len(targets),'controls':len(controls),'corruptions':len(rejected),
        'relator_lengths':[[len(w) for w in r['presentation']['relators']] for r in targets]}))


if __name__=='__main__':main()

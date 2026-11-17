#!/usr/bin/env python3
"""Independent pointwise replay of the finite cover-group certificates."""
from collections import Counter
from copy import deepcopy
from itertools import permutations
import hashlib
import json
from pathlib import Path

from check_finite_cover_certificates import check
from verify_three_sheet_splices import planar_arcs
from verify_three_sheet_structure import canonical
from parabolic_lifts import oriented_diagram  # certificate index transport only

HERE=Path(__file__).resolve().parent


def inverses(images):
    out={}
    for g,p in images.items():
        q=[None]*len(p)
        for i,j in enumerate(p):q[j]=i
        out[g]=q
    return out


def point_word(word,start,images,inv):
    point=start
    for g in word[::-1]:point=(images[g] if g>0 else inv[-g])[point]
    return point


def perm_word(word,images):
    inv=inverses(images)
    return [point_word(word,i,images,inv) for i in range(len(next(iter(images.values()))))]


def generated(gens):
    identity=tuple(range(len(gens[0])));seen={identity};queue=[identity]
    for p in queue:
        for q in gens:
            r=tuple(q[p[i]] for i in range(len(p)))
            if r not in seen:seen.add(r);queue.append(r)
    return seen


def main():
    core=json.loads((HERE/'finite_cover_certificates.json').read_text());assert core['schema']==1
    full=json.loads((HERE/'finite_cover_full_results.json').read_text());assert full['status']=='CONSTRUCTED'
    for p,h in full['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    source={r['name']:r for r in json.loads((HERE/'three_sheet_nonabelian_results.json').read_text())['results']}
    panel={r['name']:r for r in json.loads((HERE/'joint_a5_panel.json').read_text())['cases']}
    names=['13a_3145','13a_3147','13a_3148','13a_3149'];assert [r['name'] for r in core['cases']]==names
    assert [r['name'] for r in full['cases']]==names
    identity=list(range(25));results=[];edges=[(0,1),(0,2),(1,2)]
    p=[3,0,2,1,4];q=[4,1,3,2,0]
    alpha=[p[i%5]+5*((i//5+2)%5) for i in range(25)]
    delta=[q[i%5]+5*((i//5+3)%5) for i in range(25)]
    for case,row in zip(core['cases'],full['cases']):
        assert case['name']==row['name'];old=source[case['name']]
        assert case['generator_names']==old['remaining_generators'] and case['degree']==25
        rename=dict(zip(old['remaining_generators'],[1,2]))
        expected=[[rename[abs(a)]*(1 if a>0 else -1) for a in w] for w in old['remaining_relators']]
        assert case['relators']==expected
        # Reconstruct the formulas printed in the proof, including the old basis.
        epsilon=row['epsilon'];kind=row['model']
        assert (kind,epsilon)==dict(zip(names,[('A',1),('A',-1),('B',1),('B',-1)]))[case['name']]
        linear=[((i%5-i//5)%5)+5*(((1+epsilon)*(i%5)-i//5)%5) for i in range(25)]
        if kind=='B':linear=[(-j%5)+5*(-(j//5)%5) for j in linear]
        model_b=perm_word([1,2,-1],{1:linear,2:alpha if kind=='A' else delta})
        assert row['linear_gluing_permutation']==linear and row['model_A']==alpha and row['model_B']==model_b
        source_b=model_b
        if case['name']=='13a_3148':source_b=perm_word([-1,-1,2,-1,-1],{1:alpha,2:model_b})
        if case['name']=='13a_3149':source_b=perm_word([1,1,2],{1:alpha,2:model_b})
        assert case['A']==alpha and case['B']==source_b
        check(case)
        pair={1:case['A'],2:case['B']};pair_inv=inverses(pair)
        assert all(point_word(w,i,pair,pair_inv)==i for w in expected for i in range(25))
        commutator=perm_word([1,2,-1,-2],pair);witness=case['commutator_witness']
        assert commutator[witness['point']]==witness['image']!=witness['point']
        images={int(g):p for g,p in row['original_cover_generator_images'].items()}
        assert set(images)==set(old['generators']) and all(sorted(p)==identity for p in images.values())
        assert [images[g] for g in old['remaining_generators']]==[case['A'],case['B']]
        assert row['original_cover_relators']==old['relators']
        inv=inverses(images)
        assert all(point_word(w,i,images,inv)==i for w in old['relators'] for i in range(25))
        # Rebuild the actual branched words from independent geometric arc traces.
        pd=panel[case['name']]['pd'];cr,_=planar_arcs(pd);transport=[c[2] for c in oriented_diagram(pd)]
        labels=[old['arc_labels'][j] for j in transport];assert len(set(labels))>=2
        actual=[]
        for o,a,b in cr:
            for start in range(3):
                sheet=start;word=[]
                for arc in [o,a,o,b]:
                    p,q=edges[labels[arc]];g=transport[arc]+1
                    if sheet in [p,q]:
                        if g not in old['tree_generators']:word.append(g if sheet==p else -g)
                        sheet=q if sheet==p else p
                assert sheet==start;actual.append(word)
        assert Counter(tuple(canonical(w)) for w in actual)==Counter(tuple(canonical(w)) for w in old['relators'])
        assert all(point_word(w,i,images,inv)==i for w in actual for i in range(25))
        results.append({'name':case['name'],'degree':25,'original_two_generator_relators_checked':len(expected),
                        'original_full_cover_relators_checked':len(actual),'full_cover_generators':len(images),
                        'commutator_witness':witness,'all_relators_fix_every_point':True})
    # Verify the structured vertex action without trusting the producer's composition.
    vertex=full['trefoil_vertex_permutations'];p=[3,0,2,1,4];q=[4,1,3,2,0]
    assert vertex['alpha']==[p[i%5]+5*((i//5+2)%5) for i in range(25)]
    assert vertex['delta']==[q[i%5]+5*((i//5+3)%5) for i in range(25)]
    assert perm_word([1,1],{1:vertex['delta']})==perm_word([1,1,1],{1:vertex['alpha']})==vertex['fibre']
    assert perm_word([1,-2],{1:vertex['delta'],2:vertex['alpha']})==vertex['meridian']
    assert len(generated([tuple(vertex['alpha']),tuple(vertex['delta'])]))==300
    projection=generated([tuple(p),tuple(q)])
    even={g for g in permutations(range(5)) if sum(g[i]>g[j] for i in range(5) for j in range(i+1,5))%2==0}
    assert projection==even
    translations={point_word([1]*a+[2]*b,0,{1:vertex['meridian'],2:vertex['fibre']},inverses({1:vertex['meridian'],2:vertex['fibre']})) for a in range(5) for b in range(5)}
    assert translations==set(range(25))
    # Meaningful rejection controls: malformed, false relation, commuting and wrong witness.
    bad=[]
    c=deepcopy(core['cases'][0]);c['A'][0]=c['A'][1];bad.append(c)
    c=deepcopy(core['cases'][0]);c['relators'].append([1]);bad.append(c)
    c=deepcopy(core['cases'][0]);c['A']=identity[:];c['B']=identity[:];bad.append(c)
    c=deepcopy(core['cases'][0]);c['commutator_witness']['image']=c['commutator_witness']['point'];bad.append(c)
    for c in bad:
        try:check(c)
        except ValueError:pass
        else:raise AssertionError('invalid certificate accepted')
    cyclic=json.loads((HERE/'three_sheet_structure_verification.json').read_text());assert cyclic['status']=='PASS'
    controls=[(r['name'],r['proved_cyclic_order']) for r in cyclic['results'] if 'proved_cyclic_order' in r]
    assert controls==[('13a_3291',2),('13a_3346',2),('13n_4637',4)]
    provenance=json.loads((HERE/'splice_knot_provenance.json').read_text());upper=json.loads((HERE/'three_sheet_splice_verification.json').read_text())
    assert provenance['status']=='PASS exact named PD provenance' and upper['status']=='PASS exact sharp-four evidence'
    assert [r['name'] for r in upper['results']]==names and all(r['meridional_lower']==r['bridge_upper']==4 for r in upper['results'])
    for receipt in [cyclic,provenance,upper]:
        for path,h in receipt['input_sha256'].items():assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==h
    paths=['finite_cover_quotient_plan.json','build_finite_cover_quotients.py','finite_cover_certificates.json',
           'check_finite_cover_certificates.py','finite_cover_full_results.json','verify_finite_cover_quotients.py',
           'three_sheet_nonabelian_results.json','three_sheet_structure_verification.json','three_sheet_splice_verification.json',
           'splice_knot_provenance.json','joint_a5_panel.json','verify_three_sheet_splices.py','verify_three_sheet_structure.py','parabolic_lifts.py']
    report={'status':'PASS independent finite commutator certificates','cases':results,'finite_certificate_count':4,
            'degree':25,'vertex_image_order':300,'peripheral_translation_group_order':25,'rejected_corruptions':len(bad),
            'exact_cyclic_controls':controls,'new_named_cases_beyond_retained_four':0,
            'assurance':'All original two-generator and full branched-cover relators checked pointwise, independently of producer composition. Finite positive witnesses need no inference from the amalgam normal-form theorem or residual finiteness. Named PD/upper and cover-rank criteria retain separate checked dependencies; not a formal theorem-prover verification or global priority audit.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'finite_cover_quotient_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'degree':25,'cases':results,'cyclic_controls':controls},indent=2))


if __name__=='__main__':main()

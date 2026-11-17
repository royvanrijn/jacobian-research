#!/usr/bin/env python3
"""Independent unrestricted affine-point cover-cocycle reconstruction."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path

from verify_ternary_transvection_control import (POINTS,act,inverse_matrix,mat_rank,
    signed_crossings)

HERE=Path(__file__).resolve().parent


def check_module(module,words):
    matrices=module['matrices'];n=len(matrices);inverses=[inverse_matrix(a) for a in matrices];rows=[]
    for word in words:
        block=[[],[],[]]
        # Give each generator each of the three independent translation
        # coordinates, then evaluate the affine word rightmost first at zero.
        for variable in range(3*n):
            generator,coordinate=divmod(variable,3);point=(0,0,0)
            unit=tuple(int(i==coordinate) for i in range(3))
            for letter in reversed(word):
                i=abs(letter)-1;translation=unit if i==generator else (0,0,0)
                if letter>0:point=tuple((a+b)%3 for a,b in zip(act(matrices[i],point),translation))
                else:point=act(inverses[i],tuple((a-b)%3 for a,b in zip(point,translation)))
            for j in range(3):block[j].append(point[j])
        rows.extend(block)
        for start in POINTS:
            point=start
            for letter in reversed(word):point=act(matrices[abs(letter)-1] if letter>0 else inverses[abs(letter)-1],point)
            assert point==start
    assert rows==module['coefficient_matrix']
    dimension=3*n-mat_rank(rows,3*n);basis=module['kernel_basis']
    assert dimension==module['Z1_dimension']==len(basis)==mat_rank(basis,3*n)
    assert all(sum(x*y for x,y in zip(row,v))%3==0 for row in rows for v in basis)
    fixed=[v for v in POINTS if all(act(a,v)==v for a in matrices)]
    assert len(fixed)==3**module['fixed_space_dimension']
    assert module['H1_dimension']==dimension-3+module['fixed_space_dimension']
    assert module['rank_lower_bound']==(dimension+2)//3
    return {'Z1_dimension':dimension,'H1_dimension':module['H1_dimension'],'fixed_vectors':len(fixed),
        'cover_rank_lower_bound':module['rank_lower_bound']}


def main():
    report=json.loads((HERE/'ternary_cover_module_results.json').read_text())
    assert report['status']=='COMPLETE' and report['sharp_names']==[]
    for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    control=json.loads((HERE/'ternary_reflection_control.json').read_text())
    positive={**report['positive_control'],'matrices':[[1,0,0,0,1,0,0,0,1]]*9}
    def words(pd):return [[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in signed_crossings(pd)]+[[1]]
    checked_control=check_module(positive,words(control['pd']));assert checked_control['Z1_dimension']==9
    bases=json.loads((HERE/'ternary_reflection_results.json').read_text())['results']
    catalogue=json.loads((HERE/'ternary_reflection_catalogue.json').read_text())['elements']
    lookup={(r['name'],i):b for r in bases for i,b in enumerate(r['bases']) if len(set(b['arc_labels']))>1}
    panel={c['name']:c for c in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    result=[]
    for row in report['results']:
        base=lookup[row['name'],row['representation_index']];meridians=[catalogue[i] for i in base['arc_labels']]
        natural=row['modules'][0]['matrices'];dual=row['modules'][1]['matrices']
        assert row['core_words']==words(panel[row['name']]['pd'])
        for a,m in zip(natural,meridians):
            assert all(act(a,v)==act(m,act(meridians[0],v)) for v in POINTS)
        units=[tuple(int(i==j) for i in range(3)) for j in range(3)]
        for a,b in zip(natural,dual):
            for u in units:
                for v in units:assert sum(x*y for x,y in zip(act(a,u),act(b,v)))%3==sum(x*y for x,y in zip(u,v))%3
        tests=[check_module(m,row['core_words']) for m in row['modules']]
        assert all(m['Z1_dimension']==4 and m['H1_dimension']==1 and m['fixed_vectors']==1 for m in tests)
        result.append({'name':row['name'],'representation_index':row['representation_index'],'modules':tests})
    assert len(result)==16
    bad=deepcopy(report['results'][0]['modules'][0]);bad['coefficient_matrix'][0][0]=(bad['coefficient_matrix'][0][0]+1)%3
    try:check_module(bad,report['results'][0]['core_words'])
    except AssertionError:bad_coefficient=True
    else:raise AssertionError('accepted false ordinary coefficient')
    bad=deepcopy(report['results'][0]['modules'][0]);bad['kernel_basis']=bad['kernel_basis'][:-1]
    try:check_module(bad,report['results'][0]['core_words'])
    except AssertionError:bad_basis=True
    else:raise AssertionError('accepted incomplete ordinary kernel')
    # Finite coset inputs for the next13-point homology gate, on the same
    # known knot with cover group C3*C3*C3. Kurosh topology is written separately.
    source_control=json.loads((HERE/'ternary_transvection_control.json').read_text())
    def line(v):
        first=next(x for x in v if x)
        return tuple(x*(1 if first==1 else 2)%3 for x in v)
    lines=sorted({line(v) for v in POINTS if any(v)});assert len(lines)==13
    projective=[];cycles=[]
    for m in source_control['chain_matrices'][:3]:
        a=tuple(m[4*i+j] for i in range(3) for j in range(3))
        p=[lines.index(line(act(a,v))) for v in lines];projective.append(p)
        unseen=set(range(13));lengths=[]
        while unseen:
            start=min(unseen);point=start;length=0
            while point in unseen:unseen.remove(point);length+=1;point=p[point]
            assert point==start;lengths.append(length)
        assert sorted(lengths)==[1,1,1,1,3,3,3];cycles.append(sorted(lengths))
    reached={0};queue=[0]
    for point in queue:
        for p in projective:
            if p[point] not in reached:reached.add(p[point]);queue.append(p[point])
    assert len(reached)==13
    vertices=13+sum(len(c) for c in cycles);edges=3*13;free=edges-vertices+1
    factors=sum(c.count(1) for c in cycles);assert (vertices,edges,free,factors)==(34,39,6,12)
    projective_control={'points':13,'generator_cycle_lengths':cycles,'graph_vertices':vertices,
        'graph_edges':edges,'free_rank_from_graph':free,'C3_vertex_factors':factors,
        'H1_F3_dimension_from_written_Kurosh_argument':free+factors,'base_H1_F3_dimension':3,
        'homology_gain':free+factors-3,'sharp_gain_threshold':13,
        'assurance':'Finite permutations/cycles/connectivity checked; subgroup structure F6*(*12 C3) uses the separately written graph-of-groups/Kurosh argument. Target13point covers have not been tested.'}
    files=['verify_ternary_cover_modules.py','verify_ternary_transvection_control.py',
        'ternary_cover_module_results.json','ternary_base_verification.json','ternary_reflection_results.json',
        'ternary_reflection_catalogue.json','ternary_transvection_control.json','parabolic_panel.json']
    receipt={'status':'PASS independent natural-dual ternary cover modules','control':checked_control,
        'representation_orbits':16,'module_tests':32,'results':result,
        'false_coefficient_rejected':bad_coefficient,'incomplete_kernel_rejected':bad_basis,
        'projective_index13_control':projective_control,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Only natural3 and dual3 overF3 through the new16 representation orbits. Full ordinary cover cocycles computed without extension/root conditions. Z1=4<7 on all32, so no sharp4; other SL3F3 modules/fields and actual cover ranks remain open.'}
    (HERE/'ternary_cover_module_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','control','representation_orbits','module_tests','false_coefficient_rejected','incomplete_kernel_rejected']},indent=2))


if __name__=='__main__':main()

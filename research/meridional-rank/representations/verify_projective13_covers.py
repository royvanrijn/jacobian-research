#!/usr/bin/env python3
"""Independent thirteen-sheet Fox/Smith and augmentation-coordinate replay."""

from collections import Counter
from copy import deepcopy
import hashlib
from itertools import product
import json
from pathlib import Path

from verify_core_subgroup_homology import fox_cover_matrix
from verify_integer_lifts import replay
from verify_ternary_transvection_control import act,signed_crossings

HERE=Path(__file__).resolve().parent


def line(v):return min(tuple(v),tuple(-x%3 for x in v))
LINES=sorted({line(v) for v in product(range(3),repeat=3) if any(v)})
PLANES=[{j for j,v in enumerate(LINES) if sum(x*y for x,y in zip(f,v))%3==0} for f in LINES]


def point_action(matrices):return [tuple(LINES.index(line(act(a,v))) for v in LINES) for a in matrices]


def plane_action(points):
    return [tuple(PLANES.index({p[i] for i in plane}) for plane in PLANES) for p in points]


def rank_mod(matrix,columns,prime):
    a=[[v%prime for v in row] for row in matrix];pivot=0
    for j in range(columns):
        k=next((i for i in range(pivot,len(a)) if a[i][j]),None)
        if k is None:continue
        a[pivot],a[k]=a[k],a[pivot];scale=pow(a[pivot][j],-1,prime)
        a[pivot]=[v*scale%prime for v in a[pivot]]
        for i in range(pivot+1,len(a)):
            value=a[i][j]
            if value:a[i]=[(x-value*y)%prime for x,y in zip(a[i],a[pivot])]
        pivot+=1
        if pivot==len(a):break
    return pivot


def factors(value):
    result=[];d=2
    while d*d<=value:
        if value%d==0:
            result.append(d)
            while value%d==0:value//=d
        d+=1
    if value>1:result.append(value)
    return result


def check_homology(row):
    images=[tuple(p) for p in row['generator_images']]
    matrix,columns=fox_cover_matrix(images,row['presentation_words'],row['tree_edges'])
    diagonal=replay(matrix,row['smith'])
    free=columns-len(diagonal);torsion=[d for d in diagonal if d>1]
    assert free==row['H1_free_rank'] and torsion==row['H1_torsion_factors']
    assert len(matrix)==row['rows'] and columns==row['columns']
    return free,torsion


def augmentation_word_matrix(images,words):
    n=len(images);inverse=[tuple(p.index(i) for i in range(13)) for p in images];rows=[]
    for word in words:
        block=[[] for _ in range(12)]
        for variable in range(12*n):
            generator,coordinate=divmod(variable,12);value=[0]*13
            for letter in reversed(word):
                i=abs(letter)-1
                unit=[int(j==coordinate)-int(j==12) for j in range(13)] if i==generator else [0]*13
                if letter<0:value=[x-y for x,y in zip(value,unit)]
                permutation=images[i] if letter>0 else inverse[i];new=[0]*13
                for j in range(13):new[permutation[j]]=value[j]
                value=new
                if letter>0:value=[x+y for x,y in zip(value,unit)]
            assert sum(value)==0
            for j in range(12):block[j].append(value[j])
        rows.extend(block)
    return rows


def check_augmentation(data,images,words):
    matrix=augmentation_word_matrix(images,words);n=12*len(images);p=data['prime']
    assert hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest()==data['matrix_sha256']
    dimension=n-rank_mod(matrix,n,p);basis=data['kernel_basis']
    assert dimension==data['Z1_dimension']==len(basis)==rank_mod(basis,n,p)
    assert all(sum(x*y for x,y in zip(row,v))%p==0 for row in matrix for v in basis)
    # Transitive permutation modules have only constant invariants; constants
    # lie in the sum-zero module exactly when p divides13.
    fixed=int(13%p==0);assert data['fixed_dimension']==fixed
    assert data['H1_dimension']==dimension-12+fixed
    assert data['rank3_gate_passed']==(dimension>=25)
    return matrix,dimension


def paired_group(points,planes):
    identity=tuple(range(13));seen={identity:identity};queue=[identity]
    generators=list(dict.fromkeys(zip(map(tuple,points),map(tuple,planes))))
    for a in queue:
        b=seen[a]
        for p,q in generators:
            x=tuple(p[a[i]] for i in range(13));y=tuple(q[b[i]] for i in range(13))
            if x in seen:assert seen[x]==y
            else:assert len(seen)<5616;seen[x]=y;queue.append(x)
    assert len(seen)==5616
    stabilizer=[b for a,b in seen.items() if a[0]==0];assert len(stabilizer)==432
    # An equivariant map from points to planes would send point0 to a plane
    # fixed by its full point stabilizer. There is no such plane.
    assert not any(all(p[j]==j for p in stabilizer) for j in range(13))
    return len(seen)


def check_alias(row,previous):
    assert row['name']==previous['name'] and row['source_group']==previous['source_group']
    assert row['presentation_words']==previous['presentation_words']
    h=row['sheet_intertwiner'];assert sorted(h)==list(range(13))
    for a,b in zip(row['generator_images'],previous['generator_images']):
        assert all(h[a[i]]==b[h[i]] for i in range(13))
    assert row['maximum_gain_all_primes_except13']==previous['maximum_gain_all_primes_except13']
    assert row['char13_Z1_dimension']==previous['augmentation_F13']['Z1_dimension']


def main():
    report=json.loads((HERE/'projective13_cover_results.json').read_text())
    assert report['status']=='COMPLETE' and len(report['results'])==32 and report['sharp_names']==[]
    for p,digest in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    source=json.loads((HERE/'ternary_cover_module_results.json').read_text())['results']
    lookup={(r['name'],r['representation_index']):r for r in source}
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    base={};verified=[];group_checks={}
    for name,row in report['base_homology'].items():
        words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in signed_crossings(panel[name]['pd'])]+[[1]]
        matrix=[]
        for word in words:
            v=[0]*13
            for letter in word:v[abs(letter)-1]+=1 if letter>0 else -1
            matrix.append(v)
        diagonal=replay(matrix,row['smith']);assert len(diagonal)==13
        torsion=[d for d in diagonal if d>1];assert torsion==row['H1_torsion_factors'] and len(torsion)==1
        base[name]=torsion[0]
    for i,row in enumerate(report['results']):
        source_row=lookup[row['name'],row['representation_index']]
        points=point_action(source_row['modules'][0]['matrices']);planes=plane_action(points)
        expected=points if row['action']=='point' else planes
        assert [list(p) for p in expected]==row['generator_images']
        expected_words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in signed_crossings(panel[row['name']]['pd'])]+[[1]]
        assert row['presentation_words']==expected_words
        if 'equivalent_to' in row:
            assert row['equivalent_to']<i
            check_alias(row,report['results'][row['equivalent_to']]);continue
        if row['name'] not in group_checks:group_checks[row['name']]=paired_group(points,planes)
        free,torsion=check_homology(row);assert free==2 and len(torsion)==2 and torsion[0]==2
        D=base[row['name']];primes=sorted(set(factors(D)+[p for d in torsion for p in factors(d)]))
        checks=[]
        for p in primes+[3]:
            if p==13:continue
            bg=int(D%p==0);bh=free+sum(d%p==0 for d in torsion)
            checks.append({'prime':p,'b1_G':bg,'b1_H':bh,'gain':bh-bg,'augmentation_Z1':12+bh-bg})
        maximum=max(c['gain'] for c in checks)
        assert maximum==row['maximum_gain_all_primes_except13']==4
        assert 2 in primes and D%2 and all(d%2==0 for d in torsion)
        module,at13=check_augmentation(row['augmentation_F13'],row['generator_images'],row['presentation_words']);assert at13==14
        assert 12*13-rank_mod(module,12*13,2)==16
        verified.append({'name':row['name'],'representation_index':row['representation_index'],'action':row['action'],
            'H1_free_rank':free,'H1_torsion_factors':torsion,'base_D':D,'prime_checks':checks,
            'maximum_Betti_gain_except13':maximum,'augmentation_Z1_F13':at13,'maximum_augmentation_Z1_all_fields':16})
    assert len(verified)==report['distinct_source_actions']==16 and len(group_checks)==8
    for name in base:
        rows=[r for r in report['results'] if r['name']==name]
        assert sorted({r['representation_index'] for r in rows})==[1,2]
        aliases=[r for r in rows if 'equivalent_to' in r];assert len(aliases)==2
        for alias in aliases:
            earlier=report['results'][alias['equivalent_to']]
            assert alias['representation_index']==2 and earlier['representation_index']==1
            assert alias['action']!=earlier['action']
    controls=[]
    for row in report['controls']:
        free,torsion=check_homology(row);assert free==6 and torsion==[3]*12
        _,at3=check_augmentation(row['augmentation_F3'],row['generator_images'],row['presentation_words']);assert at3==27
        _,at13=check_augmentation(row['augmentation_F13'],row['generator_images'],row['presentation_words']);assert at13==18
        assert free+len(torsion)-3==row['gain_F3']==15
        controls.append({'action':row['action'],'H1_free_rank':free,'H1_torsion_factors':torsion,'gain_F3':15,'Z1_F3':at3,'Z1_F13':at13})
    rejects=[];first=next(r for r in report['results'] if 'equivalent_to' not in r)
    for kind in ['altered_Smith_step','false_free_rank']:
        bad=deepcopy(first)
        if kind=='altered_Smith_step':next(s for s in bad['smith']['steps'] if s[1]=='add')[4]+=1
        else:bad['H1_free_rank']=1
        try:check_homology(bad)
        except AssertionError:rejects.append(kind)
        else:raise AssertionError('accepted altered integer homology')
    bad=deepcopy(first['augmentation_F13']);bad['kernel_basis'][0]=[0]*156
    try:check_augmentation(bad,first['generator_images'],first['presentation_words'])
    except AssertionError:rejects.append('false_characteristic13_kernel')
    else:raise AssertionError('accepted false characteristic13 kernel')
    bad=deepcopy(next(r for r in report['results'] if 'equivalent_to' in r));bad['sheet_intertwiner'][0],bad['sheet_intertwiner'][1]=bad['sheet_intertwiner'][1],bad['sheet_intertwiner'][0]
    try:check_alias(bad,report['results'][bad['equivalent_to']])
    except AssertionError:rejects.append('false_alias_intertwiner')
    else:raise AssertionError('accepted false action alias')
    degree_control=json.loads((HERE/'characteristic_degree_control_verification.json').read_text())
    assert degree_control['status']=='PASS independent positive degree-characteristic control'
    for p,digest in degree_control['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==digest
    assert degree_control['gain']==23 and degree_control['augmentation_Z1']==36 and degree_control['split_prediction']==35
    files=['projective13_cover_plan.json','projective13_cover_results.json','verify_projective13_covers.py',
        'verify_core_subgroup_homology.py','verify_integer_lifts.py','verify_ternary_transvection_control.py',
        'ternary_cover_module_results.json','parabolic_panel.json','characteristic_degree_control_verification.json']
    receipt={'status':'PASS independent projective13 all-field cover test','source_actions':32,'distinct_source_covers':16,
        'knots':8,'group_checks':group_checks,'results':verified,'controls':controls,'rejected_corruptions':rejects,
        'maximum_Betti_gain_except13':4,'maximum_augmentation_Z1_all_fields':16,'required_augmentation_Z1':25,
        'retained_quotient_kernels_per_source_among_these_maps':1,'new_target_closures':0,
        'positive_degree_characteristic_control':{'degree':13,'characteristic':13,'gain':23,'augmentation_Z1':36,'split_prediction':35,'rejected_corruptions':2},
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'scope':'Only exact32 point/dual actions from16 retained representation orbits on8 original knots;16 inequivalent labelled Gsets, one faithful SL3 quotient kernel per knot among these maps. Smith/UCT and split augmentation give all fields except13; direct module kernel covers13. Does not exhaust all index13 subgroups, SL3 maps/modules, fields for other modules or actual cover ranks. Uniform arguments remain written, unformalized and externally unreviewed.'}
    (HERE/'projective13_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','source_actions','distinct_source_covers','maximum_Betti_gain_except13','maximum_augmentation_Z1_all_fields','rejected_corruptions','new_target_closures']},indent=2))


if __name__=='__main__':main()

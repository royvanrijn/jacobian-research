#!/usr/bin/env python3
"""Independent integral/affine replay and simple-module coverage audit.

Uses the old independent integer H3 model and Smith-step replay. The
producer's quadratic-field elimination and cocycle builder are not used.
Simple-module completeness is explicitly imported from ATLAS.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
from verify_icosahedral_lifts import build_group,expand,expand_rectangular,matrix_product,compose,inverse
from verify_integer_lifts import nullspace,replay
from verify_six_point_covers import action_matrix

HERE=Path(__file__).resolve().parent


def rank(a,n,p):return n-len(nullspace(a,n,p))
def mm(a,b,p=None):
    out=matrix_product(a,b)
    return [[x%p for x in r] for r in out] if p else out
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def invmod(a,p):
    n=len(a);rows=[list(r)+e for r,e in zip(a,eye(n))]
    for j in range(n):
        k=next(k for k in range(j,n) if rows[k][j]%p)
        rows[j],rows[k]=rows[k],rows[j];c=pow(rows[j][j],-1,p)
        rows[j]=[x*c%p for x in rows[j]]
        for k in range(n):
            if k!=j:
                c=rows[k][j];rows[k]=[(x-c*y)%p for x,y in zip(rows[k],rows[j])]
    return [r[n:] for r in rows]


def affine_matrix(images,inverses,words,p=None):
    # Evaluate affine words from right to left on zero, one translation at a time.
    n=len(images[0]);out=[]
    for word in words:
        block=[[0]*(n*len(images)) for _ in range(n)]
        for variable in range(n*len(images)):
            point=[0]*n
            for letter in reversed(word):
                i=abs(letter)-1;a=images[i] if letter>0 else inverses[i]
                point=[sum(r[k]*point[k] for k in range(n)) for r in a]
                if i==variable//n:
                    point=[v+(int(j==variable%n) if letter>0 else -inverses[i][j][variable%n]) for j,v in enumerate(point)]
                if p:point=[x%p for x in point]
            for j,x in enumerate(point):block[j][variable]=x
        out.extend(block)
    return out


def principal(images):
    n=len(images[0]);return [[a[i][j]-int(i==j) for j in range(n)] for a in images for i in range(n)]


def algebra_dimension(images,p,scalars=None):
    n=len(images[0]);basis={};queue=[eye(n)]
    def insert(a):
        v=[x%p for r in a for x in r]
        for j,row in sorted(basis.items()):
            c=v[j];v=[(x-c*y)%p for x,y in zip(v,row)]
        j=next((j for j,x in enumerate(v) if x),None)
        if j is None:return False
        c=pow(v[j],-1,p);basis[j]=[x*c%p for x in v];return True
    insert(queue[0]);generators=images+([scalars] if scalars is not None else [])
    for a in queue:
        for g in generators:
            b=mm(a,g,p)
            if insert(b):queue.append(b)
    return len(basis)


def tau_scalar(n):
    return [[int(i//2==j//2)*[[0,1],[1,1]][i%2][j%2] for j in range(2*n)] for i in range(2*n)]


def quotient_check(images,item):
    line=item['invariant_line'];pivot=item['pivot'];free=[i for i in range(3) if i!=pivot]
    q=[]
    for i in free:
        q.append([((int(i==j)-(line[i][0] if j==pivot else 0))%2,
                   (-line[i][1] if j==pivot else 0)%2) for j in range(3)])
    projection=expand_rectangular(q)
    assert rank(projection,6,2)==4
    small=[expand(a,2) for a in item['images']]
    for old,new in zip(images,small):assert mm(projection,old,2)==mm(new,projection,2)
    return small


def finite_check(item):
    n=item['dimension'];images=[expand(a,n) for a in item['images']]
    inverses=[invmod(a,2) for a in images]
    return images,inverses


def field_dimensions(row,expanded,invexpanded,expected):
    p,t=expected['prime'],expected['tau_root']
    if t is None:
        z=78-rank(expanded,78,p);v=6-rank(invexpanded,6,p)
        assert z%2==v%2==0;z//=2;v//=2
    else:
        a=[[(x+t*y)%p for x,y in r] for r in row['cocycle_matrix']]
        b=[[(x+t*y)%p for x,y in r] for r in row['principal_matrix']]
        z=39-rank(a,39,p);v=3-rank(b,3,p)
    assert (z,v,z-3+v)==(expected['Z1_dimension'],expected['invariants'],expected['H1_dimension'])
    assert z<=4
    return {'prime':p,'tau_root':t,'Z1_dimension':z,'invariants':v}


def joint_order(a,b):
    unit=tuple(range(5));seen={(unit,unit)};queue=[(unit,unit)]
    for x,y in queue:
        for u,v in zip(a,b):
            z=(compose(x,u),compose(y,v))
            if z not in seen:seen.add(z);queue.append(z)
    return len(seen)


def main():
    report=json.loads((HERE/'a5_all_modules_results.json').read_text())
    tensor=json.loads((HERE/'a5_tensor_pair_results.json').read_text())
    assert report['status']==tensor['status']=='COMPLETE'
    for receipt in [report,tensor]:
        for p,h in receipt['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    group,_,_,unit=build_group(report['H3_group_certificate'])
    catalogue=json.loads((HERE/'joint_a5_catalogue.json').read_text())
    original=json.loads((HERE/'joint_a5_results.json').read_text())['results'][:16]
    maps={(r['name'],m['orbit_index']):m for r in original for m in r['epimorphism_representatives']}
    six=json.loads((HERE/'six_point_cover_verification.json').read_text())
    old={(r['name'],r['orbit_index']):r for r in six['results']}
    aug=json.loads((HERE/'augmentation_cover_verification.json').read_text())
    assert aug['all_index_five_subgroups_covered'] and len(maps)==16
    results=[];by_key={};controls=[]
    for row in report['results']:
        key=row['name'],row['orbit_index'];mapping=maps[key]
        elements=[tuple(catalogue['elements'][i]) for i in mapping['arc_images']]
        images=[group[(p,0)] for p in elements];inverses=[group[(inverse(p),0)] for p in elements]
        assert images==[expand(m,3) for m in row['images']]
        expanded=affine_matrix(images,inverses,row['words'])
        assert expanded==expand_rectangular(row['cocycle_matrix'])
        invexpanded=principal(images);assert invexpanded==expand_rectangular(row['principal_matrix'])
        diagonal=replay(expanded,row['integer_smith']);invdiag=replay(invexpanded,row['invariant_smith'])
        assert diagonal in [[1]*71+[5],[1]*70+[7,35]] and invdiag==[1]*4+[2,2]
        fields=[field_dimensions(row,expanded,invexpanded,r) for r in row['field_dimensions']]
        assert row['generic_Z1']==3 and row['generic_invariants']==0
        assert row['candidate_primes']==([2,3,5,7] if 7 in diagonal else [2,3,5])
        assert old[key]['H1_H_torsion_factors'][0] % 7==0 if 7 in diagonal else old[key]['H1_H_torsion_factors'][0] % 7!=0
        binary=[]
        first=quotient_check(images,row['binary_simple_modules'][0])
        # An explicit nonsplit module is included in the composition-series scope.
        transpose_images=[list(map(list,zip(*a))) for a in images]
        assert rank(principal(transpose_images),6,2)==6
        for j,item in enumerate(row['binary_simple_modules']):
            small,invs=finite_check(item)
            if j:
                assert item['images']==[[[(a+b)%2,b%2] for a,b in matrix] for matrix in row['binary_simple_modules'][0]['images']]
            coc=affine_matrix(small,invs,row['words'],2)
            assert coc==[[x%2 for x in r] for r in expand_rectangular(item['cocycle_matrix'])]
            z=4*13-rank(coc,4*13,2);v=4-rank(principal(small),4,2)
            assert z==6 and v==0 and item['Z1_dimension']==3 and item['H1_dimension']==1
            assert algebra_dimension(small,2,tau_scalar(2))==8
            binary.append({'simple_dimension':2,'Z1_dimension':z//2,'H1_dimension':1})
            if not controls:
                for free_rank in [2,3]:
                    free_z=(4*free_rank-rank(affine_matrix(small[:free_rank],invs[:free_rank],[],2),4*free_rank,2))//2
                    assert free_z==2*free_rank
                    controls.append({'group':'F'+str(free_rank),'module_dimension':2,'Z1_dimension':free_z,'rank_lower_bound':(free_z+1)//2})
        assert any((a[0][1]+a[3][1])%2 for a in row['binary_simple_modules'][0]['images'])
        # Absolute simplicity of the remaining modular models needed for coverage.
        assert algebra_dimension([[[x%3 for x in r] for r in a] for a in images[1:]],3,tau_scalar(3))==18
        assert any(sum(a[2*i][2*i+1] for i in range(3))%3 for a in images)
        p5=[[[((row['images'][g][3*i+j][0]+3*row['images'][g][3*i+j][1])%5) for j in range(3)] for i in range(3)] for g in range(13)]
        assert algebra_dimension(p5,5)==9
        # Reuse all-field V4/V5 bounds, rather than re-running their cover searches.
        assert old[key]['V4_characteristic5']['Z1_dimension']==4
        assert [m['Z1_dimension'] for m in old[key]['V5_nonsplit_characteristics']]==[6,5]
        assert max(r['gain'] for r in old[key]['field_checks'])<=1
        assert any(r['prime']==5 and r['gain']==0 for r in old[key]['field_checks'])
        assert any(r.get('generic') and r['gain']==0 for r in old[key]['field_checks'])
        results.append({'name':key[0],'orbit_index':key[1],
                        'W3_integer_nonunit_factors':[d for d in diagonal if d!=1],
                        'W3_fields':fields,'binary_simples':binary,
                        'all_A5_module_ratio_bound':'3/2',
                        'ordinary_rank_bound_at_most':2})
        by_key[key]=row
    # Universal modular model identification for V4 at2/3 and V5 at5.
    a,b=(1,0,3,2,4),(2,1,4,3,0)
    four=[action_matrix(g) for g in [a,b]]
    assert algebra_dimension(four,2)==algebra_dimension(four,3)==16
    # Obtain a coherent six-point action from the earlier retained coefficient model.
    oldsix=json.loads((HERE/'six_point_cover_results.json').read_text())['A5_six_point_action']
    perms=[tuple(p) for p in catalogue['elements']]
    five=[action_matrix(oldsix[perms.index(g)]) for g in [a,b]]
    assert algebra_dimension(five,5)==25
    pairs=[];orders={}
    for row in tensor['results']:
        key0=(row['name'],row['orbit_indices'][0]);key1=(row['name'],row['orbit_indices'][1])
        if row['name'] not in orders:
            x,y=[maps[k]['arc_images'] for k in [key0,key1]]
            orders[row['name']]=joint_order([tuple(catalogue['elements'][i]) for i in x],[tuple(catalogue['elements'][i]) for i in y])
            assert orders[row['name']]==3600
        i,j=row['frobenius_choices'];aa=by_key[key0]['binary_simple_modules'][i]['images'];bb=by_key[key1]['binary_simple_modules'][j]['images']
        def fmul(x,y):return [(x[0]*y[0]+x[1]*y[1])%2,(x[0]*y[1]+x[1]*y[0]+x[1]*y[1])%2]
        expected=[[fmul(a[2*i+j],b[2*k+l]) for i in range(2) for k in range(2) for j in range(2) for l in range(2)] for a,b in zip(aa,bb)]
        assert row['images']==expected
        images,inverses=finite_check(row)
        coc=affine_matrix(images,inverses,by_key[key0]['words'],2)
        assert coc==[[x%2 for x in r] for r in expand_rectangular(row['cocycle_matrix'])]
        assert (8*13-rank(coc,8*13,2))//2==row['Z1_dimension']==6
        assert rank(principal(images),8,2)==8 and row['H1_dimension']==2
        assert algebra_dimension(images,2,tau_scalar(4))==32
        pairs.append({'name':row['name'],'frobenius_choices':row['frobenius_choices'],'Z1_dimension':6,'H1_dimension':2,'threshold_Z1':9})
    assert len(pairs)==16
    source=(HERE/'sources/atlas-a5.html').read_text()
    assert 'All faithful irreducibles in characteristic 2' in source and 'Both faithful irreducibles in characteristic 5' in source
    paths=['verify_a5_all_modules.py','SIMPLE_MODULE_BOUNDS_AND_A5_EXHAUSTION.tex','a5_all_modules_results.json','a5_tensor_pair_results.json',
           'a5_all_modules_plan.json','a5_tensor_pair_plan.json','joint_a5_results.json','joint_a5_catalogue.json',
           'augmentation_cover_verification.json','six_point_cover_verification.json','six_point_cover_results.json',
           'verify_icosahedral_lifts.py','verify_integer_lifts.py','sources/atlas-a5.html']
    receipt={'status':'PASS','kernel_classes':16,'new_W3_matrices':16,'new_binary_simple_actions':32,
             'results':results,'tensor_pairs':pairs,'joint_image_orders':orders,
             'simple_module_completeness':'Imported from primary ATLAS A5 page; modular model identification checked by full matrix-algebra spans. Scalar-extension and composition-series proofs written separately.',
             'all_A5_module_ratio_bound':'3/2','all_A5_ordinary_cocycle_rank_bounds_at_most':2,
             'new_sharp_cases':0,'rational_A5_induced_cover_H1':'zero, including the60-sheet kernel covers',
             'controls':controls,'nonsplit_W3_characteristic2_extension_checked':True,
             'scope':'All finite-dimensional coefficient modules through one A5 quotient at a time, all fields; natural2x2 product tensors only in characteristic2. Does not establish actual group rank2 or exclude other product modules/quotients/restricted meridional tests.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'a5_all_modules_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','kernel_classes','new_W3_matrices','new_binary_simple_actions','all_A5_module_ratio_bound','new_sharp_cases']},indent=2))


if __name__=='__main__':main()

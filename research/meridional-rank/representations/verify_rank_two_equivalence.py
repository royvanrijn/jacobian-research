#!/usr/bin/env python3
"""Independent symmetric-orbit, affine-point and rank-two control replay."""
from collections import Counter
from itertools import combinations,permutations,product
from math import factorial
import hashlib
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram
from search_lift_panel import seed_certificate
from verify_integer_lifts import replay,nullspace,crossings_from_walk
from verify_joint_a5 import matrix_model

HERE=Path(__file__).resolve().parent


def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):return tuple(a.index(i) for i in range(len(a)))


def closure(generators):
    n=len(generators[0]);identity=tuple(range(n));seen={identity};todo=[identity]
    for a in todo:
        for b in generators:
            c=mul(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return seen


def orbit_catalogue(degree,arity=4):
    group=list(permutations(range(degree)));edges=list(combinations(range(degree),2));trans=[]
    for i,j in edges:
        p=list(range(degree));p[i],p[j]=p[j],p[i];trans.append(tuple(p))
    actions=[[trans.index(mul(mul(g,t),inverse(g))) for t in trans] for g in group]
    seen=set();out=[]
    for values in product(range(len(trans)),repeat=arity):
        if values in seen:continue
        orbit={tuple(a[i] for i in values) for a in actions};seen.update(orbit)
        if len(closure([trans[i] for i in values]))==factorial(degree):out.append(min(orbit))
    assert len(seen)==len(trans)**arity
    return edges,trans,sorted(out)


def transport_a5(elements):
    source=[[elements.index(mul(a,b)) for b in elements] for a in elements]
    psl,orders,autos,_=matrix_model()
    source_orders=[]
    for i in range(60):
        x=i;n=1
        while x:x=source[x][i];n+=1
        source_orders.append(n)
    a,b=next((a,b) for a in range(60) for b in range(60) if source_orders[a]==2 and source_orders[b]==3 and source_orders[source[a][b]]==5)
    for c,d in ((c,d) for c in range(60) for d in range(60) if orders[c]==2 and orders[d]==3 and orders[psl[c][d]]==5):
        image={0:0};todo=[0];ok=True
        for g in todo:
            for x,y in [(a,c),(b,d)]:
                h,k=source[g][x],psl[image[g]][y]
                if h in image:
                    if image[h]!=k:ok=False;break
                else:image[h]=k;todo.append(h)
            if not ok:break
        if ok and len(image)==len(set(image.values()))==60:break
    else:raise AssertionError('no A5/PSL transport')
    identification=[image[i] for i in range(60)]
    assert all(identification[source[a][b]]==psl[identification[a]][identification[b]] for a in range(60) for b in range(60))
    inv=[next(j for j in range(60) if psl[i][j]==0) for i in range(60)]
    inner={tuple(psl[psl[g][i]][inv[g]] for i in range(60)):g for g in range(60)}
    return identification,psl,inv,autos,inner


def matrix_from_points(crossings,degree,labels):
    edges=list(combinations(range(degree),2));n=len(labels);rows=[]
    for o,a,b in crossings:
        word=[o,a,o,b];block=[[0]*n for _ in range(degree)]
        for variable in range(n):
            point=[0]*degree
            for i in reversed(word):
                u,v=edges[labels[i]];point[u],point[v]=point[v],point[u]
                if i==variable:point[u]+=1;point[v]-=1
            for j,x in enumerate(point):block[j][variable]=x
        rows.extend(r for r in block if any(r))
    return rows


def labels_from_trace(crossings,seeds,trace,assignment,trans):
    index={p:i for i,p in enumerate(trans)};labels={s:trans[i] for s,i in zip(seeds,assignment)}
    for o,a,b in trace:
        assert o in labels and a in labels and b not in labels
        labels[b]=mul(mul(labels[o],labels[a]),labels[o])
    assert len(labels)==len(crossings)
    if any(mul(mul(labels[o],labels[a]),labels[o])!=labels[b] for o,a,b in crossings):return None
    return [index[labels[i]] for i in range(len(crossings))]


def B(x,y):return ((x&1)*((y>>1)&1)^((x>>1)&1)*(y&1))&1
def Q(x):return (x&1)^((x>>1)&1)^((x&1)*((x>>1)&1))
def rank(rows,columns):return columns-len(nullspace(rows,columns,2))


def affine_generators(roots,z_arrays,epsilon):
    q=len(z_arrays);size=1<<(q+2);out=[]
    for i,v in enumerate(roots):
        u=v|sum(a[i]<<(j+2) for j,a in enumerate(z_arrays))
        out.append(tuple(x^u if B(x,v)^epsilon[i] else x for x in range(size)))
    return out


def check_relators(generators,crossings):
    identity=tuple(range(len(generators[0])))
    assert all(mul(g,g)==identity for g in generators)
    assert all(mul(mul(mul(generators[o],generators[a]),generators[o]),generators[b])==identity for o,a,b in crossings)


def positive_controls():
    diagrams=json.loads((HERE/'rank_two_control_diagrams.json').read_text());results=[]
    for source in diagrams['controls']:
        crossings=crossings_from_walk(source['gauss_code']);n=len(crossings);upper=source['bridge_upper']
        seeds,trace,_=seed_certificate(n,crossings,upper,10000,time.monotonic()+15);assert seeds
        edges,trans,reps=orbit_catalogue(3,upper);best=None
        for assignment in reps:
            labels=labels_from_trace(crossings,seeds,trace,assignment,trans)
            if labels is None:continue
            rows=matrix_from_points(crossings,3,labels);basis=nullspace(rows,n,2)
            if best is None or len(basis)>len(best[2]):best=(labels,rows,basis,assignment)
        assert best is not None
        labels,rows,basis,assignment=best;roots=[sum(1<<j for j in edges[i] if j<2) for i in labels]
        principal=[[B(s,v) for v in roots] for s in [1,2]];assert rank(principal,n)==2
        span=principal[:];classes=[]
        for vector in basis:
            if rank(span+[vector],n)>len(span):classes.append(vector);span.append(vector)
        assert len(basis)==upper and len(classes)==upper-2
        # The uniform radical-lifting lemma predicts all basis pairs work,
        # not merely a specially selected maximal independent set.
        for z in basis:
            for e in basis:check_relators(affine_generators(roots,[z],e),crossings)
        epsilon=classes[0];z_arrays=classes[1:]
        generators=affine_generators(roots,z_arrays,epsilon);check_relators(generators,crossings)
        group=closure(generators);size=1<<(upper-1)
        assert len(group)==24*8**(upper-3) and {g[0] for g in group}==set(range(size))
        # A dependent radical/translation pair has the explicit invariant.
        dependent=affine_generators(roots,[epsilon],epsilon)
        invariant=[((x>>2)&1)^Q(x&3) for x in range(8)]
        assert len(set(invariant))==2 and all(invariant[g[x]]==invariant[x] for g in dependent for x in range(8))
        # Independent four-point action identifies an actual S4 lift.
        plane=affine_generators(roots,[],epsilon);plane_group=closure(plane)
        assert len(plane_group)==24 and all(sum(g[i]!=i for i in range(4))==2 for g in plane)
        # The deleted permutation module's root space agrees with S3's.
        e4=list(combinations(range(4),2));labels4=[e4.index(tuple(i for i in range(4) if g[i]!=i)) for g in plane]
        rows4=matrix_from_points(crossings,4,labels4)
        assert len(nullspace(rows4,n,2))==len(basis)
        results.append({'parameters':source['parameters'],'seeds':seeds,'seed_trace':trace,'S3_assignment':assignment,
                        'S3_arc_labels':labels,'binary_lift_dimension':len(basis),'admissible_classes':classes,
                        'principal_basis':principal,'all_basis_pair_relators_checked':len(basis)**2,
                        'directions':[v|sum(a[i]<<(j+2) for j,a in enumerate(z_arrays)) for i,v in enumerate(roots)],
                        'offsets':epsilon,'affine_permutations':generators,'dimension':upper-1,'orbit_size':size,
                        'image_order':len(group),'dependent_quadratic_invariant':invariant,
                        'S4_lift_dimension':len(basis),'known_family_only':True})
    return results


def binary3_images():
    catalogue=json.loads((HERE/'gl3f2_catalogue.json').read_text())
    matrices=catalogue['matrices'];perms=[]
    def act(a,v):return sum(((row&v).bit_count()%2)<<i for i,row in enumerate(a))
    for a in matrices:perms.append(tuple(act(a,v)-1 for v in range(1,8)))
    assert len(set(perms))==168 and all(sorted(g)==list(range(7)) for g in perms)
    assert {tuple(a) for a in matrices}=={a for a in product(range(8),repeat=3) if rank([[(r>>j)&1 for j in range(3)] for r in a],3)==3}
    index={p:i for i,p in enumerate(perms)};table=[[index[mul(a,b)] for b in perms] for a in perms]
    inv=[index[inverse(a)] for a in perms]
    def generated(gens):
        seen={0};todo=[0]
        for x in todo:
            for y in gens:
                z=table[x][y]
                if z not in seen:seen.add(z);todo.append(z)
        return frozenset(seen)
    subgroups={frozenset(s) for s in catalogue['subgroups']};assert len(subgroups)==179
    assert frozenset([0]) in subgroups and frozenset(range(168)) in subgroups
    classified=[]
    for h in subgroups:
        gens=[];growing=frozenset([0])
        for x in sorted(h):
            if x not in growing:gens.append(x);growing=generated(gens)
        assert growing==h
        for x in range(168):
            if x not in h:assert generated(gens+[x]) in subgroups
        normal_generators=[]
        for t in h-{0}:
            if table[t][t]!=0:continue
            conjugates={table[table[g][t]][inv[g]] for g in h}
            if generated(conjugates)==h:normal_generators.append(t)
        if not normal_generators:continue
        assert len(h) in [2,6,24,168]
        fixed=[v for v in range(1,8) if all(act(matrices[g],v)==v for g in gens)]
        covectors=[f for f in range(1,8) if all(all(((f&act(matrices[g],v)).bit_count()-(f&v).bit_count())%2==0 for v in range(8)) for g in gens)]
        directions=[]
        for t in normal_generators:
            differences={act(matrices[t],v)^v for v in range(8)}-{0}
            assert len(differences)==1;directions.append(next(iter(differences)))
        vector_rows=lambda vs:[[(v>>j)&1 for j in range(3)] for v in vs]
        if len(h)==6:
            assert len(fixed)==1 and rank(vector_rows(directions),3)==2 and rank(vector_rows(directions+fixed),3)==3
            assert len(set(directions))==3
            assert len({tuple(sorted(set(directions)).index(act(matrices[g],v)) for v in sorted(set(directions))) for g in h})==6
        if len(h)==24:
            sylow3={frozenset([0,t,table[t][t]]) for t in h-{0} if table[table[t][t]][t]==0}
            assert len(sylow3)==4;sylow3=sorted(sylow3,key=lambda a:tuple(sorted(a)))
            actions={g:tuple(sylow3.index(frozenset(table[table[g][t]][inv[g]] for t in a)) for a in sylow3) for g in h}
            assert len(set(actions.values()))==24
            assert all(sum(actions[t][i]!=i for i in range(4))==2 for t in normal_generators)
            assert (len(fixed),len(covectors)) in [(1,0),(0,1)]
            if fixed:
                assert fixed[0] not in directions
                cosets=sorted({min(v,v^fixed[0]) for v in range(1,8)}-{0})
                assert len(cosets)==3
                assert len({tuple(cosets.index(min(act(matrices[g],v),act(matrices[g],v)^fixed[0])) for v in cosets) for g in h})==6
            else:
                assert all((covectors[0]&v).bit_count()%2==0 for v in directions)
                plane=[v for v in range(1,8) if (covectors[0]&v).bit_count()%2==0]
                assert len({tuple(plane.index(act(matrices[g],v)) for v in plane) for g in h})==6
        classified.append({'subgroup_order':len(h),'fixed_dimension':(len(fixed)+1).bit_length()-1,
                           'fixed_covector_dimension':(len(covectors)+1).bit_length()-1,
                           'normal_generating_involutions':len(normal_generators)})
    old=json.loads((HERE/'gl3f2_verification.json').read_text());assert old['status']=='PASS' and old['binary_root_checks']==58
    for p,h in old['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    return {'subgroups_exhausted':179,'eligible_image_subgroups_by_order':dict(Counter(r['subgroup_order'] for r in classified)),
            'module_fixed_space_types':dict(Counter(f"{r['subgroup_order']}:{r['fixed_dimension']}:{r['fixed_covector_dimension']}" for r in classified)),
            'full_GL3_binary_root_actions_retained':58,'all_binary3_root_bounds_at_most':3}


def main():
    source=json.loads((HERE/'small_symmetric_full_panel_results.json').read_text());assert source['status']=='COMPLETE'
    for receipt in [source,json.loads((HERE/'rank_two_control_diagrams.json').read_text())]:
        for p,h in receipt['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    cases=json.loads((HERE/'joint_a5_panel.json').read_text())['cases']
    cover=json.loads((HERE/'joint_a5_results.json').read_text())['results']
    elements=[tuple(p) for p in json.loads((HERE/'joint_a5_catalogue.json').read_text())['elements']]
    identification,table,inverses,automorphisms,inner=transport_a5(elements)
    models={r:orbit_catalogue(r) for r in [2,3,4,5]};rows=[];decks=Counter();base_counts=Counter();binary3=Counter()
    for case,old,result in zip(cases,cover,source['results']):
        assert case['name']==old['name']==result['name'];cr=[c[:3] for c in oriented_diagram(case['pd'])];n=len(cr)
        seeds=result['seeds'];trace=result['seed_trace'];assert seeds==old['seeds'] and len(set(seeds))==4
        seen=set(seeds)
        for o,a,b in trace:assert o in seen and a in seen and b not in seen;seen.add(b)
        assert len(seen)==n
        for before,extension in zip(old['epimorphism_representatives'],result['A5_deck_extensions']):
            images=[identification[i] for i in before['arc_images']]
            implementers=[tuple(a) for a in automorphisms if all(a[i]==inverses[i] for i in images)]
            assert len(implementers)<=1
            dtype='absent' if not implementers else 'inner' if implementers[0] in inner else 'outer'
            assert extension['deck_type']==dtype;decks[dtype]+=1
            if dtype=='inner':assert identification[elements.index(tuple(extension['implementer']))]==inner[implementers[0]]
        maximum=0;s3dimension=None
        for d in result['degrees']:
            r=d['degree'];edges,trans,reps=models[r];assert len(reps)==d['seed_orbits_tested']
            found={}
            for assignment in reps:
                labels=labels_from_trace(cr,seeds,trace,assignment,trans)
                if labels is not None:found[assignment]=labels
            assert found=={tuple(c['seed_assignment']):c['arc_labels'] for c in d['certificates']}
            base_counts[r]+=len(found)
            for certificate in d['certificates']:
                matrix=matrix_from_points(cr,r,certificate['arc_labels']);assert matrix==certificate['matrix']
                diagonal=replay(matrix,certificate['smith']);bound=n-diagonal.count(1);maximum=max(maximum,bound)
                if r==3:
                    ell=n-sum(x%2!=0 for x in diagonal);binary3[ell]+=1;s3dimension=ell
        assert maximum==result['maximum_bound_all_degrees_fields']<4
        rows.append({'name':case['name'],'maximum_symmetric_bound_all_fields':maximum,'S3_binary_lift_dimension':s3dimension})
        if len(rows)%75==0:print('PASS independent symmetric panel',len(rows),flush=True)
    assert decks=={'inner':170} and base_counts=={2:293,3:26,4:21,5:0}
    assert binary3=={2:5,3:21}
    controls=positive_controls();binary3_classification=binary3_images()
    # Explicit known T45 outer-deck control: x^4=y^5=1, m=x*y^-1.
    x=(1,2,3,0,4);y=(1,2,3,4,0);m=mul(x,inverse(y))
    assert len(closure([x,y]))==120 and sum(m[i]!=i for i in range(5))==2
    assert mul(m,m)==tuple(range(5)) and m not in elements
    assert len(closure([mul(x,x),y]))==60
    paths=['verify_rank_two_equivalence.py','small_symmetric_full_panel_plan.json','small_symmetric_full_panel_results.json',
           'rank_two_equivalence_plan.json','rank_two_control_diagrams.json','joint_a5_results.json','joint_a5_panel.json',
           'joint_a5_catalogue.json','parabolic_lifts.py','verify_joint_a5.py','verify_integer_lifts.py','search_lift_panel.py',
           'gl3f2_catalogue.json','gl3f2_verification.json']
    report={'status':'PASS','targets':len(rows),'seed_orbits':{str(r):len(m[2]) for r,m in models.items()},
            'deck_histogram':dict(decks),'symmetric_base_counts':dict(base_counts),'S3_binary_lift_dimensions':dict(binary3),
            'results':rows,'positive_known_pretzel_controls':controls,'T45_outer_control':{'x':x,'y':y,'meridian':m,'image_order':120},
            'binary3_linear_image_classification':binary3_classification,
            'new_open_cases':0,'original_common_form_no_representation_targets':272,'original_common_form_maximum_bound':3,
            'scope':'Finite evidence for separately written uniform radical-lifting/rank2 equivalence proofs, not an all-knot proof. Independent full tuple-orbit marking and fixed-trace permutation closure, reverse affine-point integer matrices and Smith replay; source PD parser and existing Smith replay shared. Known-family controls only.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'rank_two_equivalence_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','targets':293,'base_counts':dict(base_counts),'binary3':dict(binary3),
                      'positive_control_orders':[r['image_order'] for r in controls],'new_open_cases':0},indent=2))


if __name__=='__main__':main()

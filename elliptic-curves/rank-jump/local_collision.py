#!/usr/bin/env python3
"""Retrospective root-collision incidence and common CT-annihilator audit."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import retrospective as r
import blocks

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'LOCAL_COLLISION_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_local_collision_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_local_collision_v1.json'
BASE='661246f'


def canonical(vs):
    p=r.basis(vs)
    for k in sorted(p):
        for j in p:
            if j!=k and p[j]>>k&1:p[j]^=p[k]
    return [p[k] for k in sorted(p)]


def orthogonal(vs,n):
    p={v.bit_length()-1:v for v in canonical(vs)}
    return canonical([(1<<j)^sum(((v>>j)&1)<<k for k,v in p.items()) for j in range(n) if j not in p])


def intersection(*spaces,n=20):
    return orthogonal([v for s in spaces for v in orthogonal(s,n)],n)


def lift(word,base):
    v=0
    for i,b in enumerate(base):
        if word>>i&1:v^=b
    return v


def coordinates(v,base):
    p={}
    for i,b in enumerate(base):
        word=1<<i
        while b:
            k=b.bit_length()-1
            if k not in p:p[k]=(b,word);break
            b^=p[k][0];word^=p[k][1]
        else:raise ValueError('dependent basis')
    word=0
    for k in sorted(p,reverse=True):
        if v>>k&1:v^=p[k][0];word^=p[k][1]
    if v:raise ValueError('outside basis span')
    return word


def capture():
    bindings={}
    def source(path):
        raw=subprocess.check_output(['git','show',BASE+':'+path],cwd=r.ROOT)
        bindings[path]=r.digest(raw);return json.loads(raw)
    local=source(r.read(PROTOCOL)['local_source'])
    old=source('artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json')
    ct=source('artifacts/generated-results/elliptic-curves/rank_jump_block_inputs_v1.json')['ct']
    old_by_u={int(x['parameter_u']):[b['mask'] for b in x['W_u_basis']] for x in old['runs']}
    for x in local['runs']:
        u=int(x['parameter_u']);base=[b['mask'] for b in x['W_u_basis']]
        if u in old_by_u:assert base==old_by_u[u]
        c=next(c for c in ct if c['u']==u)
        if 'admissible_anchor_masks' in c:assert base==c['admissible_anchor_masks']
        assert len(c['matrix'])==len(base)
    rows=[]
    for x in local['runs']:
        rows.append({k:x[k] for k in ('parameter_u','W_u_dimension','W_u_basis','all_local_kummer_images_complete','combined_condition_rank',
          'discriminant_multiplier','newly_bad_primes_relative_to_anchor','finite_local_conditions','real_local_condition')})
        rows[-1]['W_u_basis']=[b['mask'] for b in x['W_u_basis']]
        for p in rows[-1]['finite_local_conditions']:p.pop('basis_witnesses',None)
    r.write_new(INPUT,{'schema':'rank-jump.local-collision-inputs.v1','bindings':bindings,
       'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'anchor':local['anchor'],'rows':rows,
       'ct':[{'u':c['u'],'matrix':c['matrix']} for c in ct]})


def pairing(w,z,B):
    return sum(((r.pack(row)&z).bit_count()&1) for i,row in enumerate(B) if w>>i&1)&1


def compare_local(row,model,points):
    u=r.F(row['parameter_u']);A,B=map(r.F,model[3:]);D=1+A*u*u+B*u**3
    assert str(D)==row['discriminant_multiplier']
    old=[];new=[];checks=[];place_rows=[]
    assert row['all_local_kummer_images_complete']
    for pinfo in row['finite_local_conditions']:
        p=pinfo['prime'];assert pinfo['point_kummer_basis_complete']
        columns=r.transpose(pinfo['known_span_quotient_rows'])
        actual=canonical(columns)
        isnew=p in row['newly_bad_primes_relative_to_anchor']
        (new if isnew else old).extend(actual)
        place_rows.append({'prime':p,'new':isnew,'constraint_basis':actual,'codimension':len(actual)})
        if not isnew:continue
        n=D.numerator;exp=0
        while n%p==0:n//=p;exp+=1
        eligible=p>2 and D.denominator%p and u.numerator%p and u.denominator%p and exp==1
        if not eligible:
            checks.append({'prime':p,'eligible':False,'valuation_D':exp});continue
        assert (4*r.mod(A,p)**3+27*r.mod(B,p)**2)%p
        root=pow(r.mod(u,p),-1,p)
        assert (root**3+r.mod(A,p)*root+r.mod(B,p))%p==0
        predicted=r.pack([r.point_signature(model,P,[(p,[root])]) for P in points])
        checks.append({'prime':p,'eligible':True,'valuation_D':exp,'noncolliding_anchor_root_mod_p':root,
          'predicted_constraint':predicted,'actual_constraint_basis':actual,'matches':canonical([predicted])==actual})
    real=canonical(r.transpose(row['real_local_condition']['known_span_quotient_rows']))
    all_constraints=canonical(old+new+real)
    assert len(all_constraints)==row['combined_condition_rank']
    W=orthogonal(all_constraints,20)
    assert W==canonical(row['W_u_basis'])
    assert len(W)==row['W_u_dimension']
    return {'u':int(u),'dimension':len(W),'basis':W,'constraint_basis':all_constraints,
      'inherited_finite_constraint_rank':r.rank(old),'new_prime_constraint_rank':r.rank(new),'real_constraint_rank':r.rank(real),
      'inherited_new_constraint_intersection_dimension':r.rank(old)+r.rank(new)-r.rank(old+new),
      'local_places':place_rows,'root_character_checks':checks}


def build(check=False):
    inp=r.read(INPUT);assert inp['protocol_sha256']==r.digest(PROTOCOL.read_bytes())
    model=inp['anchor']['short_model_ainvariants'];points=inp['anchor']['known_points_on_short_model'];r.short(model,points)
    local=[compare_local(row,model,points) for row in inp['rows']]
    indexed={row['u']:row for row in local};cts={c['u']:c['matrix'] for c in inp['ct']}
    bases={int(row['parameter_u']):row['W_u_basis'] for row in inp['rows']}
    radicals={u:canonical([lift(w,bases[u]) for w in orthogonal(map(r.pack,B),len(B))]) for u,B in cts.items()}
    pairs=[]
    for u,v in combinations(sorted(indexed),2):
        C=intersection(indexed[u]['basis'],indexed[v]['basis']);rad_common=intersection(radicals[u],radicals[v])
        details=[]
        for t in (u,v):
            words=[coordinates(c,bases[t]) for c in C];B=cts[t]
            rectangular=[r.pack([pairing(w,1<<j,B) for j in range(len(B))]) for w in words]
            self_form=[r.pack([pairing(w,z,B) for z in words]) for w in words]
            cap=len(C)-r.rank(rectangular);weak=len(C)-r.rank(self_form)
            assert cap==len(intersection(C,radicals[t]))
            details.append({'u':t,'rectangular_obstruction_rank':r.rank(rectangular),'necessary_soluble_dimension_cap':cap,
              'common_space_self_pairing_rank':r.rank(self_form),'weaker_self_radical_dimension':weak})
        pairs.append({'u':u,'v':v,'common_local_dimension':len(C),'common_local_basis':C,
           'pairings_on_common_space':details,'common_necessary_soluble_dimension':len(rad_common),
           'common_necessary_soluble_basis':rad_common})
    deformations=[u for u in indexed if u]
    jointlocal=intersection(*(indexed[u]['basis'] for u in deformations))
    jointrad=intersection(*(radicals[u] for u in deformations))
    eligible=[c for row in local for c in row['root_character_checks'] if c['eligible']]
    # Enumerate only subsets of the six already computed spaces (at most 63).
    minimal_killers=[]
    for size in range(1,7):
        for us in combinations(sorted(deformations),size):
            if not intersection(*(radicals[u] for u in us)) and not any(set(x)<=set(us) for x in minimal_killers):minimal_killers.append(list(us))
    result={'schema':'rank-jump.local-collision.v1','input_sha256':r.digest(INPUT.read_bytes()),'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
      'script_sha256':r.digest(Path(__file__).read_bytes()),'local_rows':local,'embedded_restricted_radicals':{str(u):x for u,x in radicals.items()},
      'pair_comparisons':pairs,'joint_local_basis':jointlocal,'joint_local_dimension':len(jointlocal),
      'joint_necessary_soluble_basis':jointrad,'joint_necessary_soluble_dimension':len(jointrad),
      'minimal_tested_deformation_sets_with_zero_common_radical':minimal_killers,
      'summary':{'eligible_new_prime_events':len(eligible),'matched_events':sum(c['matches'] for c in eligible),
        'primary_endpoint':'PASS_ON_ALL_ELIGIBLE_RETAINED_EVENTS' if all(c['matches'] for c in eligible) else 'HYPOTHESIS_REFUTED',
        'pairs_with_strictly_weaker_common_self_test':sum(any(x['weaker_self_radical_dimension']>x['necessary_soluble_dimension_cap'] for x in a['pairings_on_common_space']) for a in pairs)},
      'boundary':'Local conditions and CT arithmetic inherited from pinned complete-subspace certificates; exact new linear algebra and independently evaluated root-character prediction. No full Selmer or rank bound and no point search.'}
    if check:
        if r.read(OUTPUT)!=result:raise ValueError('local-collision replay mismatch')
        print('PASS deterministic local-collision replay')
    else:r.write_new(OUTPUT,result)
    print(json.dumps(result['summary']));print('joint local',len(jointlocal),'joint radical',len(jointrad))
    for x in local:print('u',x['u'],'local dim',x['dimension'],'old/new cuts',x['inherited_finite_constraint_rank'],x['new_prime_constraint_rank'])
    print('minimal radical killers',minimal_killers)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('capture','build','check'));a=p.parse_args()
    capture() if a.mode=='capture' else build(a.mode=='check')

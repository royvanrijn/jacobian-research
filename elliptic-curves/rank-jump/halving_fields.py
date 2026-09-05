#!/usr/bin/env python3
"""Maximal joint E[2]-torsor extensions and explicit retrospective quartics."""
import argparse
from itertools import combinations,product
from pathlib import Path
import retrospective as r
import local_collision as lc
import ct_variation as cv
from cubic_bridge import Cubic

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'HALVING_FIELD_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_halving_fields_v1.json'


def compose(A,B):
    return tuple(cv.transpose([cv.apply(A,c) for c in cv.transpose(B,2)],2))


def subspaces(n):
    """All reduced binary subspaces, without enumerating generating lists."""
    for d in range(n+1):
        for pivots in combinations(range(n),d):
            slots=[(i,j) for i,k in enumerate(pivots) for j in range(k) if j not in pivots]
            for mask in range(1<<len(slots)):
                rows=[1<<p for p in pivots]
                for bit,(i,j) in enumerate(slots):
                    if mask>>bit&1:rows[i]|=1<<j
                yield rows


def module_certificate():
    I=(1,2);G=[I]+[A for A in product(range(4),repeat=2) if r.rank(A)==2 and A!=I]
    index={A:i for i,A in enumerate(G)};table=[[index[compose(A,B)] for B in G] for A in G]
    cocycles=[]
    for tail in product(range(4),repeat=5):
        c=(0,)+tail
        if all(c[table[i][j]]==c[i]^cv.apply(A,c[j]) for i,A in enumerate(G) for j in range(6)):
            cocycles.append(list(c))
    boundaries=[[cv.apply(A,v)^v for A in G] for v in range(4)]
    assert sorted(cocycles)==sorted(boundaries) and len(cocycles)==4
    centralizer=[A for A in product(range(4),repeat=2) if all(compose(A,B)==compose(B,A) for B in G)]
    assert centralizer==[(0,0),I]
    assert r.rank([A[0]|A[1]<<2 for A in G])==4
    checks=[]
    for copies in (1,2,3):
        n=2*copies;count=0;stable=[]
        for H in subspaces(n):
            count+=1;piv=r.basis(H)
            def act(A,v):return sum(cv.apply(A,(v>>(2*i))&3)<<(2*i) for i in range(copies))
            if not all(r.reduce(act(A,h),piv)==0 for A in G for h in H):continue
            relations=[]
            for lam in range(1<<copies):
                def evaluate(h):
                    z=0
                    for i in range(copies):
                        if lam>>i&1:z^=(h>>(2*i))&3
                    return z
                if all(evaluate(h)==0 for h in H):relations.append(lam)
            relbase=lc.canonical(relations)
            assert len(H)==2*(copies-len(relbase))
            if len(H)<n:assert len(relbase)>0
            stable.append({'basis':H,'coordinate_relation_basis':relbase})
        assert count==(5,67,2825)[copies-1]
        assert len(stable)==(2,5,16)[copies-1]
        affine_counts=[]
        for A in G:
            count=0
            for offsets in product(range(4),repeat=copies):
                if all(any(cv.apply(A,x)^v==x for x in range(4)) for v in offsets):count+=1
            rank=r.rank([a^(1<<i) for i,a in enumerate(A)])
            assert count==2**(rank*copies)
            affine_counts.append({'linear_part':list(A),'rank_of_g_minus_I':rank,
                                  'offset_tuples_fixing_a_point_in_every_factor':count})
        numerator=sum(x['offset_tuples_fixing_a_point_in_every_factor'] for x in affine_counts)
        density=r.F(numerator,6*4**copies)
        assert density==r.F(1,3)+r.F(1,2**(copies+1))+r.F(1,6*4**copies)
        checks.append({'copies':copies,'all_subspace_count':(5,67,2825)[copies-1],
          'stable_subspaces':stable,'affine_fixed_point_counts':affine_counts,
          'joint_fixed_point_fraction':str(density)})
    return {'group_matrices':[list(A) for A in G],'multiplication_table':table,
      'normalized_maps_tested':4**5,'one_cocycles':cocycles,'one_coboundaries':boundaries,
      'H1_dimension':0,'centralizer':[list(A) for A in centralizer],'action_algebra_dimension':4,
      'small_multiplicity_checks':checks}


def trim(a):
    while a and a[-1]==0:a.pop()
    return a


def remainder(a,b,p):
    a=trim([x%p for x in a]);b=trim([x%p for x in b])
    if not b:raise ZeroDivisionError
    inv=pow(b[-1],-1,p)
    while len(a)>=len(b):
        d=len(a)-len(b);c=a[-1]*inv%p
        for j,x in enumerate(b):a[d+j]=(a[d+j]-c*x)%p
        trim(a)
    return a


def gcd(a,b,p):
    while b:a,b=b,remainder(a,b,p)
    if not a:return []
    return [x*pow(a[-1],-1,p)%p for x in a]


def polynomial_power(a,e,h,p):
    def mul(a,b):
        out=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
        return remainder(out,h,p)
    out=[1]
    while e:
        if e&1:out=mul(out,a)
        a=mul(a,a);e//=2
    return out


def quartic_witnesses(h):
    found={}
    for p in r.primes(1999):
        try:hp=[r.mod(c,p) for c in h]
        except ValueError:continue
        if gcd(hp,[(i*hp[i])%p for i in range(1,5)],p)!=[1]:continue
        roots=[x for x in range(p) if sum(c*pow(x,i,p) for i,c in enumerate(hp))%p==0]
        if len(roots)==1 and '1,3' not in found:
            found['1,3']={'prime':p,'linear_roots':roots,'reduction':hp}
        if not roots and '4' not in found:
            xp2=polynomial_power([0,1],p*p,hp,p)+[0,0]
            xp2[1]=(xp2[1]-1)%p
            if gcd(hp,trim(xp2),p)==[1]:
                found['4']={'prime':p,'linear_roots':[],'reduction':hp,'gcd_with_x_p_squared_minus_x':[1]}
        if len(found)==2:return found
    raise ValueError('bounded Frobenius certificate incomplete; no group claim')


def quartic(row,index,role):
    model,points=r.short(row['model'],row['generic_points']+row['points'])
    A,B=map(r.F,model[3:]);p,q=map(r.F,points[index]);assert q
    h=[A*A-4*B*p,-4*A*p-8*B,-2*A,-4*p,r.F(1)]
    a,b,c,d=h[3],h[2],h[1],h[0]
    resolvent=[4*b*d-a*a*d-c*c,a*c-4*d,-b,r.F(1)]
    K=Cubic(A,B);rho=K.add(K.scalar(2*A),K.add(K.scale(K.theta,4*p),K.scale(K.square(K.theta),4)))
    result=K.scalar(0)
    for coeff in reversed(resolvent):result=K.add(K.mul(result,rho),K.scalar(coeff))
    assert result==K.scalar(0)
    return {'id':row['id'],'role':role,'independent_input_index':index,'short_model':model,
      'point':[str(p),str(q)],'halving_quartic_ascending':list(map(str,h)),
      'cubic_resolvent_ascending':list(map(str,resolvent)),
      'resolvent_generator_in_cubic':list(map(str,rho)),
      'quartic_discriminant_from_identity':str(4096*q**4*(-4*A**3-27*B*B)),
      'frobenius_witnesses':quartic_witnesses(h),'quartic_galois_group':'S4',
      'splitting_field_degree_over_two_division_field':4}


def build(check=False):
    inp=r.read(r.INPUT);prior=r.read(r.RESULT);local=r.read(lc.INPUT)
    assert prior['input_sha256']==r.digest(r.INPUT.read_bytes())
    indexed={x['id']:x for x in inp['rows']};characterized={};rows=[]
    for row in inp['rows']:
        out,_,_=r.characterize(row)
        assert out==next(x for x in prior['panel'] if x['id']==row['id'])
        assert out['galois']['galois_group']=='S3' and out['epsilon_M']==0
        characterized[row['id']]=out
        m=out['generic_kummer_dimension_exact'];n=out['certified_independent_subgroup_rank_exact'];q=n-m
        assert q==out['certified_independent_quotient_rank_exact']
        rows.append({'id':row['id'],'family':row['family'],'origin':row['origin'],
          'generic_mod2_dimension':m,'witness_mod2_dimension':n,'witness_quotient_mod2_dimension':q,
          'independent_input_indices':out['independent_input_indices'],
          'equation_galois_certificate':out['galois'],
          'rank_input_scope':'inherited exact global fingerprints' if 'supplied_global_kummer_fingerprints' in row else 'point membership and finite Kummer signatures recomputed',
          'generic_halving_field_degree_over_L':str(4**m),
          'witness_halving_field_degree_over_L':str(4**n),
          'relative_degree_over_generic_halving_field':str(4**q),
          'joint_galois_group':'(F2^2)^'+str(n)+' semidirect S3 (diagonal natural action)',
          'extension_scope':'Only the specified independent witness subgroup. Whole-curve halving field UNKNOWN.'})
    quartics=[]
    for pair in inp['pairs']:
        for side in ('high','low'):
            row=indexed[pair[side]];out=characterized[row['id']]
            quartics.append(quartic(row,out['independent_input_indices'][0],side+'_first_generic'))
            if side=='high':
                idx=out['independent_input_indices'][out['generic_kummer_dimension_exact']]
                quartics.append(quartic(row,idx,'high_first_quotient'))
    # Same class space and finite torsor fields, opposite solubility endpoints.
    W=next(x for x in local['rows'] if int(x['parameter_u'])==-1)
    B=next(x['matrix'] for x in local['ct'] if x['u']==-1)
    fixed={'parameters':[0,-1],'common_class_basis_in_anchor_coordinates':W['W_u_basis'],
      'dimension':len(B),'finite_torsor_field_degree_over_L_at_both_parameters':str(4**len(B)),
      'finite_torsor_field_is_the_same_under_labelled_E2_identification':True,
      'anchor_rationally_soluble_dimension':len(B),'deformed_restricted_CT_rank':r.rank(map(r.pack,B)),
      'deformed_necessary_soluble_dimension_cap':len(B)-r.rank(map(r.pack,B)),
      'boundary':'This is the finite E[2]-torsor field. The genus-one cover changes with the embedding E[2] into E_u.'}
    result={'schema':'rank-jump.halving-fields.v1','input_sha256':r.digest(r.INPUT.read_bytes()),
      'panel_report_sha256':r.digest(r.RESULT.read_bytes()),'local_input_sha256':r.digest(lc.INPUT.read_bytes()),
      'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'script_sha256':r.digest(Path(__file__).read_bytes()),
      'module_certificate':module_certificate(),'panel':rows,'paired_quartics':quartics,'fixed_class_counterexample':fixed,
      'summary':{'panel_rows':len(rows),'point_signature_replays':sum(x['rank_input_scope'].startswith('point') for x in rows),
        'inherited_signature_rows':sum(x['rank_input_scope'].startswith('inherited') for x in rows),
        'halving_quartics':len(quartics),'S4_quartics_with_bounded_Frobenius_witnesses':len(quartics),
        'joint_extension_collapse':'RULED_OUT_FOR_ALL_SPECIFIED_INDEPENDENT_SUBGROUPS',
        'solubility_prediction_from_finite_torsor_field_alone':'REFUTED_BY_FIXED_CLASS_COUNTEREXAMPLE'},
      'boundary':'No new curve, point, rank or prospective selector. Maximal degree is a consequence of known H1 independence, equally valid for insoluble classes. Arithmetic constructions of several rational points remain possible and unexplained.'}
    if check:
        if r.read(OUTPUT)!=result:raise ValueError('halving-field replay mismatch')
        print('PASS halving-field replay')
    else:r.write_new(OUTPUT,result)
    print(result['summary'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('build','check'));a=p.parse_args();build(a.mode=='check')

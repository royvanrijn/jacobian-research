#!/usr/bin/env python3
"""Supplemental retrospective local-support, CT block and height diagnostics."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import subprocess
import zipfile
import retrospective as r

INPUT=r.OUT/'rank_jump_block_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_block_report_v1.json'


def capture():
    bindings={}; heights=[];ct=[]
    archive=r.OUT/'rank25_followup_diagnostics_v2.zip'
    bindings[str(archive.relative_to(r.ROOT))]=r.digest(archive.read_bytes())
    with zipfile.ZipFile(archive) as z:
        for i in range(3):
            name=f'artifacts/local/elliptic-curves/specialized-rank25-parity-sample-v1/input-{i}.json'
            raw=z.read(name);d=json.loads(raw)
            heights.append({'id':d['family']+':'+d['parameter'],'gram':d['metric_gram'],
              'generic_rank':d['generic_dimension'],'model':d['curve'],
              'points':d['initial_state']['state']['reductions']['points'],
              'source_member':name,'source_sha256':r.digest(raw)})
    for u in [-3,-2,0,1,2,3]:
        tag='um'+str(-u) if u<0 else 'u'+str(u)
        name='artifacts/generated-results/elliptic-curves/fixed_field_comparison_'+tag+'_v1.json.gz'
        raw=subprocess.check_output(['git','show',r.BASE+':'+name],cwd=r.ROOT)
        bindings[name]=r.digest(raw);d=json.loads(gzip.decompress(raw))
        ct.append({'u':u,'matrix':d['ct']['matrix'],'admissible_anchor_masks':d['admissible_masks'],
          'summary':d['summary'],'source':name})
    name='artifacts/generated-results/elliptic-curves/fixed_cubic_u_minus1_cassels_tate_v1.json'
    raw=subprocess.check_output(['git','show',r.BASE+':'+name],cwd=r.ROOT)
    bindings[name]=r.digest(raw);d=json.loads(raw)
    ct.append({'u':-1,'matrix':d['arithmetic']['matrix'],'source':name})
    legacy=[]
    for filename in ('fermigier_rank_jump_fingerprints_v1.json','nagao_section7_rank_jump_fingerprint_v1.json'):
        name='artifacts/generated-results/elliptic-curves/'+filename
        raw=subprocess.check_output(['git','show',r.BASE+':'+name],cwd=r.ROOT)
        bindings[name]=r.digest(raw)
        for f in json.loads(raw)['fingerprints']:
            legacy.append({k:f[k] for k in ('label','certified_rank_lower_bound','generic_rank','global_or_canonical_model','exact_generic_embedding','quotient_structure','degree_visibility')})
            legacy[-1]['parameter']=f.get('canonical_parameter_u',f.get('constructor_parameter_T'))
            legacy[-1]['source']=name
    name='artifacts/generated-results/elliptic-curves/half_lattice_heldout_273_302_verification_v1.json'
    raw=subprocess.check_output(['git','show',r.BASE+':'+name],cwd=r.ROOT);bindings[name]=r.digest(raw)
    heldout=[{k:f[k] for k in ('curve','configuration','dimension','displayed_public_rank','exact_heldout_quotient_rank_recovered','heldout_dimension')} for f in json.loads(raw)['results']]
    r.write_new(INPUT,{'schema':'rank-jump.block-inputs.v1','bindings':bindings,'heights':heights,'ct':ct,'legacy':legacy,'heldout':heldout})


def ct_normal_form(matrix):
    n=len(matrix);assert all(len(row)==n for row in matrix)
    assert all(matrix[i][i]==0 and matrix[i][j]==matrix[j][i] for i in range(n) for j in range(n))
    rows=[r.pack(row) for row in matrix]
    def pair(x,y):
        return sum(((rows[i]&y).bit_count()&1) for i in range(n) if x>>i&1)&1
    remaining=[1<<i for i in range(n)];hyperbolic=[]
    while True:
        ij=next(((i,j) for i,x in enumerate(remaining) for j,y in enumerate(remaining[:i]) if pair(x,y)),None)
        if ij is None:break
        i,j=ij;a,b=remaining[i],remaining[j];hyperbolic.append([a,b])
        remaining=[v^(a if pair(v,b) else 0)^(b if pair(v,a) else 0) for k,v in enumerate(remaining) if k not in (i,j)]
    words=[v for h in hyperbolic for v in h]+remaining
    assert r.rank(words)==n
    for i,x in enumerate(words):
        for j,y in enumerate(words):assert pair(x,y)==int(i<len(hyperbolic)*2 and j==(i^1))
    return {'dimension':n,'pairing_rank':2*len(hyperbolic),'hyperbolic_pairs_in_input_coordinates':hyperbolic,
      'restricted_radical_basis_in_input_coordinates':remaining,'restricted_radical_dimension':len(remaining),
      'nonzero_classes_certifiably_obstructed':2**n-2**len(remaining),
      'soluble_dimension_upper_bound_in_this_subspace':len(remaining),
      'block_scope':'Exact symplectic normal form, not canonical arithmetic components or a full Selmer decomposition.'}


def exact_inverse(a):
    n=len(a);m=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k=next(k for k in range(j,n) if m[k][j]);m[j],m[k]=m[k],m[j]
        s=m[j][j];m[j]=[x/s for x in m[j]]
        for i in range(n):
            if i!=j:
                c=m[i][j];m[i]=[x-c*y for x,y in zip(m[i],m[j])]
    return [row[n:] for row in m]


def smith_diagonal(matrix):
    """Integer Smith reduction using only unimodular row/column operations."""
    a=[list(map(int,row)) for row in matrix];n=len(a);m=len(a[0])
    for k in range(min(n,m)):
        candidates=[(abs(a[i][j]),i,j) for i in range(k,n) for j in range(k,m) if a[i][j]]
        if not candidates:break
        _,i,j=min(candidates);a[k],a[i]=a[i],a[k]
        for row in a:row[k],row[j]=row[j],row[k]
        while True:
            for i in range(k+1,n):
                if a[i][k]:
                    q=a[i][k]//a[k][k];a[i]=[x-q*y for x,y in zip(a[i],a[k])]
                    if a[i][k]:a[i],a[k]=a[k],a[i]
                    break
            else:
                for j in range(k+1,m):
                    if a[k][j]:
                        q=a[k][j]//a[k][k]
                        for row in a:row[j]-=q*row[k]
                        if a[k][j]:
                            for row in a:row[j],row[k]=row[k],row[j]
                        break
                else:
                    bad=next(((i,j) for i in range(k+1,n) for j in range(k+1,m) if a[i][j]%a[k][k]),None)
                    if bad is None:break
                    a[k]=[x+y for x,y in zip(a[k],a[bad[0]])]
        if a[k][k]<0:a[k]=[-v for v in a[k]]
    return [a[i][i] for i in range(min(n,m))]


def height_profile(d):
    # Exact rational algebra on stored decimal approximations; no claim of rigorous heights.
    H=[[F(v) for v in row] for row in d['gram']];n=len(H);g=d['generic_rank'];q=n-g
    asym=max(abs(H[i][j]-H[j][i]) for i in range(n) for j in range(n))
    H=[[(H[i][j]+H[j][i])/2 for j in range(n)] for i in range(n)]
    inv=exact_inverse([row[:g] for row in H[:g]])
    coeff=[[sum(inv[i][k]*H[k][g+j] for k in range(g)) for j in range(q)] for i in range(g)]
    S=[[H[g+i][g+j]-sum(H[g+i][k]*coeff[k][j] for k in range(g)) for j in range(q)] for i in range(q)]
    assert all(S[i][i]>0 for i in range(q))
    # Squared correlations avoid platform-dependent sqrt rounding.
    corr=[[S[i][j]**2/(S[i][i]*S[j][j]) for j in range(q)] for i in range(q)]
    def components(threshold):
        sets=[{i} for i in range(q)]
        for i in range(q):
            for j in range(i):
                if corr[i][j]>threshold:
                    a=next(s for s in sets if i in s);b=next(s for s in sets if j in s)
                    if a is not b:a.update(b);sets.remove(b)
        return [sorted(s) for s in sets]
    return {'id':d['id'],'layer':'visibility','dimension':q,'height_source':'retained decimal canonical-height approximations; no interval bounds',
      'gram_asymmetry':str(asym),'schur_complement_decimal':[[format(float(x),'.12g') for x in row] for row in S],
      'quotient_to_raw_height_ratios':[format(float(S[i][i]/H[g+i][g+i]),'.12g') for i in range(q)],
      'exact_zero_cross_entries':sum(S[i][j]==0 for i in range(q) for j in range(i)),
      'maximum_squared_correlation':format(float(max(corr[i][j] for i in range(q) for j in range(i))),'.12g'),
      'components_abs_correlation_gt_half':components(F(1,4)),
      'components_abs_correlation_gt_three_quarters':components(F(9,16)),
      'scope':'Schur complement removes generic projections. Component thresholds are descriptive and basis-dependent, not arithmetic factors or rank predictors.'}


def local_support(row):
    if 'points' not in row:return None
    model=row['model'];g=row['generic_rank'];pts=row['generic_points']+row['points'];local=[]
    for p in row['prime_list']:
        roots=r.roots_at(model[3],model[4],p)
        if not roots:continue
        v=[r.point_signature(model,P,[(p,roots)]) for P in pts]
        a=r.rank(v[:g]);b=r.rank(v)
        local.append({'prime':p,'rational_cubic_root_count':len(roots),'generic_image_rank':a,'combined_image_rank':b,'new_local_image_dimension':b-a})
    return {'id':row['id'],'layer':'incidence','local_rows':local,
      'places_with_new_local_image':sum(x['new_local_image_dimension']>0 for x in local),
      'scope':'Selected good auxiliary primes only; these are residue Kummer images, not bad-place admissibility or cubic ideal support.'}


def build(check=False):
    d=r.read(INPUT);panel=r.read(r.INPUT);pr=r.read(r.RESULT)
    ct=[dict(u=a['u'],**ct_normal_form(a['matrix'])) for a in d['ct']]
    heights=[height_profile(h) for h in d['heights']]
    local=[x for row in panel['rows'] if (x:=local_support(row))]
    groups={}
    for row in panel['rows']:
        A,B=map(F,row['model'][3:]);j=6912*A**3/(4*A**3+27*B**2)
        groups.setdefault(str(j),[]).append(row['id'])
    # Equal j identifies geometric isomorphism. Certify Q-isomorphism by rational square scaling.
    byid={row['id']:row for row in panel['rows']};duplicates=[]
    for ids in groups.values():
        if len(ids)<2:continue
        for i,x in enumerate(ids):
            for y in ids[:i]:
                a,b=map(F,byid[x]['model'][3:]);c,e=map(F,byid[y]['model'][3:]);u2=b*c/(a*e)
                sq=u2>0 and r.isqrt(u2.numerator)**2==u2.numerator and r.isqrt(u2.denominator)**2==u2.denominator
                assert sq and a==u2*u2*c and b==u2**3*e
        duplicates.append(ids)
    # Hamming proximity vs chart subspace overlap is exploratory and presentation-dependent.
    chart_geometry=[]
    for cm in pr['chart_maps']:
        bins={}
        for x in cm['chart_pair_intersections']:
            b=bins.setdefault(x['parity_hamming_distance'],[0,0]);b[0]+=1;b[1]+=x['intersection_dimension']
        chart_geometry.append({'id':cm['id'],'hamming_bins':[{'distance':k,'pair_count':v[0],'intersection_dimension_sum':v[1]} for k,v in sorted(bins.items())],
          'scope':'Descriptive in the frozen marked basis; Hamming distance is not invariant under GL(r,Z). No causal/statistical test.'})
    legacy=[]
    for f in d['legacy']:
        mat=list(map(list,zip(*f['exact_generic_embedding']['columns'])))
        diagonal=smith_diagonal(mat)
        assert diagonal==f['exact_generic_embedding']['smith_invariant_factors']
        model,_=r.short(f['global_or_canonical_model'],[])
        legacy.append({'id':f['label'],'family':'Fermigier/Mestre MW12' if f['label'].startswith('fermigier') else 'Nagao section 7 MW12',
          'parameter':f['parameter'],'generic_rank':f['generic_rank'],'curve_rank_lower_bound':f['certified_rank_lower_bound'],
          'free_quotient_rank_inside_displayed_L':len(mat)-sum(v!=0 for v in diagonal),'smith_diagonal':diagonal,
          'two_torsion_dimension_of_L_over_M':sum(v%2==0 for v in diagonal),
          'tensor_mod2_dimension_of_L_over_M':len(mat)-r.rank([r.pack([int(v)%2 for v in col]) for col in f['exact_generic_embedding']['columns']]),
          'full_MW_saturation':'UNKNOWN; Smith form concerns displayed L only',
          'galois':r.galois(model),'source':f['source']})
    out={'schema':'rank-jump.block-report.v1','input_sha256':r.digest(INPUT.read_bytes()),
       'panel_input_sha256':r.digest(r.INPUT.read_bytes()),'panel_report_sha256':r.digest(r.RESULT.read_bytes()),
       'script_sha256':r.digest(Path(__file__).read_bytes()),'ct_blocks':ct,'height_profiles':heights,'local_profiles':local,
       'distinct_curves_in_main_panel':len(groups),'repeated_curve_rows':duplicates,'chart_geometry':chart_geometry,
       'legacy_comparison_panel':legacy,'heldout_controls_with_unknown_family':d['heldout'],
       'summary':{'profiles':len(local),'profiles_with_new_selected_good_place_image':sum(x['places_with_new_local_image']>0 for x in local)}}
    if check:
        if r.read(OUTPUT)!=out:raise ValueError('block report mismatch')
        print('PASS block replay')
    else:r.write_new(OUTPUT,out)
    print(json.dumps(out['summary']));print('distinct curves',len(groups))
    for x in ct:print('CT',x['u'],x['dimension'],x['pairing_rank'],x['restricted_radical_dimension'])
    for h in heights:print('HEIGHT',h['id'],h['maximum_squared_correlation'],h['components_abs_correlation_gt_half'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('capture','build','check'));a=p.parse_args()
    capture() if a.mode=='capture' else build(a.mode=='check')

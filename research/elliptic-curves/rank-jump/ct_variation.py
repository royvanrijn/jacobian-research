#!/usr/bin/env python3
"""Exact common-space CT pencils and rational normalization gauge audit."""
import argparse
import gzip
import json
from itertools import combinations
from pathlib import Path
import subprocess
import retrospective as r
import local_collision as lc
from reciprocity import hilbert

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'CT_VARIATION_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_ct_variation_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_ct_variation_v1.json'


def apply(A,v):
    return r.pack([(a&v).bit_count()%2 for a in A])


def transpose(A,n):
    return [r.pack([(a>>j)&1 for a in A]) for j in range(n)]


def solve(A,b,n):
    """Particular solution of rows A times x=b; reject inconsistent systems."""
    pivots={}
    for i,a in enumerate(A):
        z=(b>>i)&1
        while a:
            k=a.bit_length()-1
            if k not in pivots:pivots[k]=(a,z);break
            q,t=pivots[k];a^=q;z^=t
        else:
            if z:raise ValueError('inconsistent system')
    x=0
    for k,(a,z) in sorted(pivots.items()):
        if z^((a&x).bit_count()%2):x|=1<<k
    assert x<1<<n and apply(A,x)==b
    return x


def bilinear(A,v,w):
    return (v&apply(A,w)).bit_count()%2


def capture():
    protocol=r.read(PROTOCOL);entries=[];bindings={}
    for tag in protocol['gauge_sources']:
        path='artifacts/generated-results/elliptic-curves/fixed_field_comparison_'+tag+'_v1.json.gz'
        raw=subprocess.check_output(['git','show',protocol['gauge_sources_commit']+':'+path],cwd=r.ROOT)
        bindings[path]=r.digest(raw);d=json.loads(gzip.decompress(raw))
        covers={c['mask']:c for c in d['covers']}
        for i,p in enumerate(d['ct']['pairs'][:3]):
            entries.append({'u':int(d['parameter_u']),'pair_index':i,'masks':p['masks'],
              'first_quartic':covers[p['masks'][0]]['quartic'],
              'second_quartic_leading':covers[p['masks'][1]]['quartic'][4],
              **{k:p[k] for k in ('gamma','local_terms','value','support_factors')}})
    r.write_new(INPUT,{'schema':'rank-jump.ct-variation-inputs.v1','protocol_sha256':r.digest(PROTOCOL.read_bytes()),
       'bindings':bindings,'entries':entries})


def gauge(entry,scalars):
    a=r.F(entry['second_quartic_leading']);terms=entry['local_terms'];places=[x['place'] for x in terms]
    assert len(set(places))==len(places) and all(p in places for p in (2,3,5,'infinity'))
    # This support assertion needs no fresh factorization. At an omitted odd
    # place a and each selected c are units, so (a,c)_p=1.
    for n in (abs(a.numerator),a.denominator):
        for p in places:
            if p=='infinity':continue
            while n%p==0:n//=p
        assert n==1
    product=1
    for term in terms:
        x=r.F(term['x']);g=sum(r.F(c)*x**i for i,c in enumerate(entry['gamma']))
        q=sum(r.F(c)*x**i for i,c in enumerate(entry['first_quartic']))
        assert g==r.F(term['gamma_value']) and q==r.F(term['q_value'])
        h=hilbert(a,g,term['place']);assert h==term['hilbert_symbol'];product*=h
    assert int(product==-1)==entry['value']
    experiments=[]
    for c in scalars:
        flips=[];negative=[];ratio=1;changed_product=1
        for term in terms:
            p=term['place'];h=hilbert(a,c,p);ratio*=h
            changed=hilbert(a,c*r.F(term['gamma_value']),p)
            assert changed==h*term['hilbert_symbol'];changed_product*=changed
            if h==-1:flips.append(p)
            if changed==-1:negative.append(p)
        assert ratio==1 and changed_product==product
        experiments.append({'scalar':c,'flipped_places':flips,'negative_places_after':negative,
          'global_pairing_value':int(changed_product==-1)})
    return {'u':entry['u'],'pair_index':entry['pair_index'],'masks':entry['masks'],
       'global_pairing_value':entry['value'],'negative_places_before':[t['place'] for t in terms if t['hilbert_symbol']==-1],
       'rescalings':experiments}


def pencil(A,B,base):
    """Certificate for this retained nonsingular-difference pencil, not a classifier."""
    n=len(base);D=[a^b for a,b in zip(A,B)];assert r.rank(D)==n
    T=transpose([solve(D,col,n) for col in transpose(A,n)],n)
    assert all(apply(D,apply(T,1<<i))==apply(A,1<<i) for i in range(n))
    powers=[[1<<i for i in range(n)]]
    for _ in range(6):powers.append(transpose([apply(T,c) for c in transpose(powers[-1],n)],n))
    N=lc.orthogonal(powers[5],n);U=lc.orthogonal([t^(1<<i) for i,t in enumerate(T)],n)
    assert len(N)==10 and len(U)==2 and r.rank(N+U)==n
    assert all(not bilinear(H,v,w) for H in (A,B,D) for v in N for w in U)
    assert all(not apply(powers[5],v) for v in N)
    assert any(apply(powers[4],v) for v in N)
    # Find two length-five cyclic vectors, normalized into paired chains.
    e=next(v for v in N if apply(powers[4],v));chain_e=[apply(P,e) for P in powers[:5]]
    constraints=[r.pack([bilinear(D,v,w) for w in N]) for v in chain_e]
    f=lc.lift(solve(constraints,1<<4,len(N)),N);chain_f=[apply(P,f) for P in powers[:5]]
    assert r.rank(chain_e+chain_f)==10
    for i,v in enumerate(chain_e):
        for j,w in enumerate(chain_f):
            assert bilinear(D,v,w)==int(i+j==4)
            assert bilinear(A,v,w)==int(i+j==3)
    for chain in (chain_e,chain_f):
        assert all(not bilinear(H,v,w) for H in (A,B,D) for v in chain for w in chain)
    kernels=[len(N)-r.rank([apply(P,v) for v in N]) for P in powers[:6]]
    assert kernels==[0,2,4,6,8,10]
    assert r.rank([r.pack([bilinear(D,v,w) for w in N]) for v in N])==10
    assert r.rank([r.pack([bilinear(D,v,w) for w in U]) for v in U])==2
    return {'common_basis_in_anchor_coordinates':base,'difference_matrix_packed_rows':D,'operator_T_packed_rows':T,
       'nilpotent_primary_basis_in_anchor_coordinates':[lc.lift(v,base) for v in N],
       'eigenvalue_one_basis_in_anchor_coordinates':[lc.lift(v,base) for v in U],
       'nilpotent_chain_e_in_anchor_coordinates':[lc.lift(v,base) for v in chain_e],
       'nilpotent_chain_f_in_anchor_coordinates':[lc.lift(v,base) for v in chain_f],
       'nilpotent_kernel_dimensions':kernels,'primary_dimensions':[10,2],
       'minimal_polynomial':'t^5*(t+1)','characteristic_polynomial':'t^10*(t+1)^2',
       'paired_chain_formula':'D(e_i,f_j)=[i+j=4]; A(e_i,f_j)=[i+j=3]; B=A+D; both chains individually isotropic for all three forms.',
       'scope':'Simultaneous CT obstruction components; these are not certified soluble MW components.'}


def build(check=False):
    inp=r.read(lc.INPUT);g=r.read(INPUT);protocol=r.read(PROTOCOL)
    assert g['protocol_sha256']==r.digest(PROTOCOL.read_bytes())
    bases={int(x['parameter_u']):x['W_u_basis'] for x in inp['rows']};forms={x['u']:x['matrix'] for x in inp['ct']}
    rows=[];pencils=[]
    for u,v in combinations(sorted(bases),2):
        C=lc.intersection(bases[u],bases[v]);restricted=[]
        for t in (u,v):
            words=[lc.coordinates(c,bases[t]) for c in C]
            restricted.append([r.pack([lc.pairing(a,b,forms[t]) for b in words]) for a in words])
        A,B=restricted;D=[a^b for a,b in zip(A,B)];rank=r.rank(D)
        assert all(not (d>>i&1) for i,d in enumerate(D))
        assert D==transpose(D,len(C)) and rank%2==0
        rows.append({'u':u,'v':v,'common_dimension':len(C),'common_basis_in_anchor_coordinates':C,
           'first_form_packed_rows':A,'second_form_packed_rows':B,'first_rank':r.rank(A),'second_rank':r.rank(B),
           'difference_rank':rank,'difference_radical_in_anchor_coordinates':[lc.lift(w,C) for w in lc.orthogonal(D,len(C))],
           'minimum_class_feature_dimension_for_bilinear_factorization':rank,
           'minimum_rank_two_summands_for_difference':rank//2})
        if u and v and rank==len(C):pencils.append({'u':u,'v':v,**pencil(A,B,C)})
    assert [(x['u'],x['v']) for x in pencils]==[(-1,1)]
    gauges=[gauge(e,protocol['gauge_scalars']) for e in g['entries']]
    result={'schema':'rank-jump.ct-variation.v1','matrix_input_sha256':r.digest(lc.INPUT.read_bytes()),
      'gauge_input_sha256':r.digest(INPUT.read_bytes()),'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
      'script_sha256':r.digest(Path(__file__).read_bytes()),'comparisons':rows,'nondegenerate_difference_pencils':pencils,'gauge_audits':gauges,
      'summary':{'comparison_count':len(rows),'nonzero_deformation_difference_rank_range':[min(x['difference_rank'] for x in rows if x['u'] and x['v']),max(x['difference_rank'] for x in rows if x['u'] and x['v'])],
       'gauge_entries':len(gauges),'rescalings':sum(len(x['rescalings']) for x in gauges),
       'rescalings_changing_local_attribution':sum(bool(s['flipped_places']) for x in gauges for s in x['rescalings']),
       'global_pairing_changes_under_rescaling':0},
      'boundary':'Pinned restricted CT matrices, not full Selmer forms. Linear-algebra components need an additional arithmetic interpretation. No rational point, positive rank predictor, or parameter-selection rule is produced.'}
    if check:
        if result!=r.read(OUTPUT):raise ValueError('CT variation replay mismatch')
        print('PASS CT variation replay')
    else:r.write_new(OUTPUT,result)
    print(json.dumps(result['summary']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('capture','build','check'));a=p.parse_args()
    capture() if a.mode=='capture' else build(a.mode=='check')

#!/usr/bin/env python3
"""Independent geometric/valuation witnesses for the multi-fibre comparison."""
import argparse
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path
from sage.all import ZZ
import retrospective as r
import fibre_discrimination as census
from verify_collision_defect import gcd as modgcd, evaluate as ev

HERE=Path(__file__).resolve().parent
SOURCE=r.OUT/'rank_jump_fibre_discrimination_v1.json'
COLLISION=r.OUT/'rank_jump_fibre_collision_panel_v1.json'
OUTPUT=r.OUT/'rank_jump_fibre_discrimination_verification_v1.json'
PROTOCOL=HERE/'FIBRE_DISCRIMINATION_VERIFICATION_PROTOCOL.json'


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def power(a,n):
    v=[F(1)]
    for _ in range(n):v=mul(v,a)
    return v


def transport():
    path=r.OUT/'compact_six_r17_atlas_v1.json'
    fam=next(f for f in r.read(path)['families'] if f['family']=='11952')
    rawpath=r.ROOT/fam['source'];assert census.hash_file(rawpath)==fam['source_sha256']
    raw=r.read(rawpath)['weierstrass_model'];a,b,c,d=map(F,fam['base_matrix_a_b_c_d']);u=F(fam['total_scale_from_literal_source'])
    for key,degree in [('A',8),('B',12)]:
        ans=[F(0)]*(degree+1)
        for i,x in enumerate(map(F,raw[key+'_coefficients_low_to_high'])):
            for j,y in enumerate(mul(power([b,a],i),power([d,c],degree-i))):ans[j]+=x*y/u**(degree//2)
        assert ans==list(map(F,fam[key+'_coefficients_low_to_high']))
    return {'status':'PASS_ALL_COEFFICIENTS','matrix':list(map(str,[a,b,c,d])),
            'scaling_u':str(u),'raw_source_sha256':census.hash_file(rawpath)}


def roots(g,p):
    if len(g)==1:return []
    if len(g)==2:return [(-g[0]*pow(g[1],-1,p))%p]
    assert len(g)==3
    if p==2:return [i for i in range(2) if ev(g,i)%2==0]
    delta=(g[1]*g[1]-4*g[0]*g[2])%p
    if delta==0:return [(-g[1]*pow(2*g[2],-1,p))%p]
    if pow(delta,(p-1)//2,p)!=1:return []
    # Tonelli--Shanks. The returned roots are independently substituted.
    odd=p-1;s=0
    while odd%2==0:odd//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1:z+=1
    c=pow(z,odd,p);x=pow(delta,(odd+1)//2,p);t=pow(delta,odd,p);m=s
    while t!=1:
        i=0;v=t
        while v!=1:v=v*v%p;i+=1
        assert i<m
        b=pow(c,1<<(m-i-1),p);x=x*b%p;c=b*b%p;t=t*c%p;m=i
    assert x*x%p==delta
    return sorted({((-g[1]+e*x)*pow(2*g[2],-1,p))%p for e in [-1,1]})


def determinant(q,w):
    c,b,a=q;f,e,d=w
    M=[[a,b,c,0],[0,a,b,c],[d,e,f,0],[0,d,e,f]]
    answer=0
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        answer+=sign*prod(M[i][perm[i]] for i in range(4))
    return answer


def arithmetic(block,row):
    qs=[x['primitive_form'] for x in block['compatible_forms']];n=len(qs)
    assert all(x['content_is_square'] and ZZ(x['removed_content']).is_square() for x in block['compatible_forms'])
    assert row['full_support_status']=='COMPLETE'
    assert len(row['pair_factorizations'])==len(block['pair_resultants'])
    primes=set();incident={}
    for pair,fac in zip(block['pair_resultants'],row['pair_factorizations'],strict=True):
        assert pair['indices']==fac['indices'] and pair['resultant']==fac['resultant']
        i,j=pair['indices'];assert determinant(qs[i],qs[j])==int(pair['resultant'])!=0
        assert prod(int(p)**e for p,e in fac['factors'])==abs(int(pair['resultant']))
        for p,e in fac['factors']:
            p=int(p);assert ZZ(p).is_prime(proof=True);primes.add(p);incident.setdefault(p,[]).append((i,j))
    assert sorted(primes)==list(map(int,row['certified_support_primes'])) and len(primes)==row['collision_support_size']
    basis=[];witnesses=[]
    for p in sorted(primes):
        if p<5:continue
        rr=set()
        for i,j in incident[p]:
            rr.update(roots(modgcd(qs[i],qs[j],p),p))
        for chart,root in [('finite',x) for x in sorted(rr)]+[('infinity',0)]:
            forms=qs if chart=='finite' else [q[::-1] for q in qs]
            inds=[i for i,q in enumerate(forms) if ev(q,root)%p==0]
            if len(inds)!=2:continue
            mask=sum(1<<i for i in inds)
            if r.rank(basis+[mask])==len(basis):continue
            a=[ev(forms[i],root)//p%p for i in inds]
            b=[(forms[i][1]+2*forms[i][2]*root)%p for i in inds]
            if not all(b) or (a[0]*b[1]-a[1]*b[0])%p==0:continue
            for s in range(min(p,16)):
                xx=root+p*s;vals=[ev(q,xx) for q in forms]
                units=[(v//p if i in inds else v)%p for i,v in enumerate(vals)]
                unit=prod(units)%p
                if not unit or pow(unit,(p-1)//2,p)!=1:continue
                assert all(v%p==0 and v%(p*p)!=0 if i in inds else v%p!=0 for i,v in enumerate(vals))
                assert prod(vals)//(p*p)%p==unit
                basis.append(mask);witnesses.append({'prime':str(p),'chart':chart,'root':str(root),'s':s,'mask':mask,'product_unit':str(unit)})
                break
            if len(basis)==n-1:break
        if len(basis)==n-1:break
    return {'block_key':row['block_key'],'compatible_cover_count':n,'collision_support_size':len(primes),
            'collision_primes':list(map(str,sorted(primes))), 'pair_count':comb(n,2),
            'realizable_defect_span_interval':[len(basis),n-1],
            'realizable_defect_span_exact':len(basis) if len(basis)==n-1 else 'UNKNOWN',
            'independent_local_witnesses':witnesses,
            'simultaneous_carrier_degree':2**n,'simultaneous_carrier_genus':1+2**(n-1)*(n-2)}


def compute():
    data=r.read(SOURCE);collision=r.read(COLLISION);inputs=r.read(census.INPUT)
    for obj,key in [(data,'bindings'),(collision,'bindings'),(inputs,'sources')]:
        for path,sha in obj[key].items():assert census.hash_file(r.ROOT/path)==sha
    assert census.compute()==data
    model_transport=transport()
    rows=[arithmetic(data['blocks'][row['block_key']],row) for row in collision['rows']]
    assert len(rows)==sum(b['compatible_cover_count']>=2 for b in data['blocks'].values())
    # Recheck stored rank matrix algebra on newer source cohorts; point proofs remain inherited.
    matrix_rows=0
    for fn in [census.SOURCES[1],*census.SOURCES[3:]]:
        old=r.read(r.OUT/fn)
        for curve in old.get('rows',old.get('curves')):
            proof=curve['rank_certificate'];mats=[r.pack(v) for s in proof['signatures'] for v in s['rows']]
            assert r.rank(mats)==proof['rank_lower_bound']
            assert r.rank([v&((1<<17)-1) for v in mats])==17
            matrix_rows+=1
    byid={x['observation_id']:x for x in data['rows']};accepted=[];unknown=[]
    for pair in data['height_matched_pairs']:
        hi,lo=byid[pair['high']],byid[pair['low']]
        if hi['cohort'] is None or hi['boxes'] is None:
            unknown.append(pair);continue
        accepted.append(pair)
    scored=[x for x in data['rows'] if x['source']==census.SOURCES[1]]
    score_pairs=[]
    for hi in scored:
        if hi['retained_quotient_rank']<5:continue
        lows=[lo for lo in scored if lo['cohort']==hi['cohort'] and lo['boxes']==hi['boxes'] and lo['retained_quotient_rank']==0]
        if not lows:continue
        lo=min(lows,key=lambda x:(abs(x['score_units']-hi['score_units']),x['observation_id']))
        score_pairs.append({'high':hi['observation_id'],'low':lo['observation_id'],
            'gain':hi['retained_quotient_rank'],'cover_counts':[hi['compatible_cover_count'],lo['compatible_cover_count']],
            'absolute_score_gap':abs(hi['score_units']-lo['score_units']),
            'relative_score_gap':str(F(abs(hi['score_units']-lo['score_units']),abs(hi['score_units']))),
            'height_matched':False,'coverage_notes':[hi['coverage_note'],lo['coverage_note']]})
    catalogue_context=[]
    for fn in census.SOURCES[3:]:
        for old in r.read(r.OUT/fn)['curves']:
            if old['rank_lower_bound']>=26:catalogue_context.append({'observation_id':fn+':'+old['id'],'catalogue_matches':old['icarm_matches'],'previous_matches':old['previous_matches']})
    return {'schema':'rank-jump.fibre-discrimination-verification.v1','status':'PASS','rows':rows,
        'coefficient_transport':model_transport,'source_rank_matrices_replayed':matrix_rows,
        'accepted_height_pairs':accepted,'unknown_exposure_pairs_not_treated_as_matched':unknown,
        'score_matched_pairs':score_pairs,'newer_high_observation_catalogue_context':catalogue_context,
        'bindings':{str(p.relative_to(r.ROOT)):census.hash_file(p) for p in (SOURCE,COLLISION,PROTOCOL,Path(__file__),HERE/'verify_collision_defect.py',HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS',len(d['rows']),'arithmetic blocks',d['source_rank_matrices_replayed'],'rank matrices',flush=True)
    for row in d['rows']:print(row['block_key'],row['compatible_cover_count'],row['collision_support_size'],row['realizable_defect_span_interval'],flush=True)

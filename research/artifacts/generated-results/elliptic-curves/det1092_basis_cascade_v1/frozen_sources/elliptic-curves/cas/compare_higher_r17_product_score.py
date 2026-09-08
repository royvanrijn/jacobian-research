#!/usr/bin/env python3
"""Frozen paired score comparison on saved higher-height traces; no CAS trace calls."""
import argparse
from pathlib import Path
from math import log
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import benchmark_periodic_nagao_scanner as scan
import extend_higher_r17_stratified as parent
import compare_bounded_prime_selectors as formula
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint,digest
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=parent.parent.LOCAL/'higher-r17-product-score-v1';EXCLUDED=parent.parent.LOCAL/'higher24-r17-pari-v1/protocol.json'
PRIMES=[q for q in range(5,4094) if _is_prime(q)]
def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(formula.__file__).resolve(),Path(scan.__file__).resolve(),Path(parent.engine.__file__).resolve(),Path(parent.scoring.__file__).resolve(),spec.ATLAS)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve paired score protocol')
    data=cert.read(parent.D/'result.json');gate=cert.read(parent.D/'replay.supervisor.json');excluded=cert.read(EXCLUDED)
    if data['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or len(data['rows'])!=6144 or gate['outcome']!='completed' or gate['returncode']!=0:raise ArithmeticError('complete saved6144 trace replay required')
    addresses=[{'family':r['family'],'retained_index':r['retained_index'],'parameter':r['parameter']} for r in excluded['rows']]
    if len(addresses)!=24 or len({(r['family'],r['retained_index']) for r in addresses})!=24:raise ArithmeticError('fixed prior24 address exclusions required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-r17-product-score.v1','sources':sources(),'parent_protocol_sha256':cert.hashed(parent.D/'protocol.json'),'parent_result_sha256':cert.hashed(parent.D/'result.json'),'parent_replay_sha256':cert.hashed(parent.D/'replay.supervisor.json'),'excluded_protocol_sha256':cert.hashed(EXCLUDED),'excluded_addresses':addresses,'short_primes':PRIMES,'selection_maximum_prime':32749,'validation_interval':[32771,65521],'quantization':'Reuse compare_bounded_prime_selectors.contributions: round each binary-math.log contribution to nearest1e-12 then sum integers. This is a reproducible numerical selector, not an exact transcendental ordering or a rank certificate.','product_formula':'sum over nonsingular displayed reductions of log((p+1-a_p)/p); same omitted-prime policy as the original score.','selection':'Exclude all24 prior addresses uniformly, without reading their results. Per family choose top two remaining product scores and top two remaining original S1 scores. Tie order: good-prime count descending, denominator ascending, signed numerator ascending. Merge overlaps into one candidate with both arms, no refill. All primes above32749 are validation only, excluded even from ties.','gate':'The completed higher-height population retained6144 addresses and all traces. The original fixed24 point batch has weaker observed point recovery than its lower-height predecessor. Literature and the existing small ordinary-fibre diagnostic supply an independently defined product-score alternative; test it on this saved population. Current24 outcomes are known to the researcher but neither their ranks nor points enter this implementation. Untested candidates have no point outcomes in this comparison.','references':['https://arxiv.org/pdf/2003.00077 sections2,7','https://www.dpmms.cam.ac.uk/~taf1000/papers/rankcongr.pdf section5','elliptic-curves/cas/compare_bounded_prime_selectors.py'],'future_point_scope':'After complete score replay, a separate fixed protocol may expose at most24 distinct addresses, at most two per family per arm, to identical generic17-only43/49-chart point searches. No point search is launched or implied by this score script.','limits':{'wall_seconds':300,'rss_bytes':1610612736,'workers':1,'checkpoint_rows':128},'boundaries':'No new parameter, finite-field trace or point search. The entire pool was previously truncated by short-prime S1; this cannot recover discarded addresses. Reweighting, validation association or a later bounded point outcome is not a rank predictor theorem, superiority guarantee, rank upper bound or universal-novelty claim.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['parent_protocol_sha256']!=cert.hashed(parent.D/'protocol.json') or p['parent_result_sha256']!=cert.hashed(parent.D/'result.json') or p['excluded_protocol_sha256']!=cert.hashed(EXCLUDED):raise ArithmeticError('frozen comparison inputs changed')
    return p

def short_scores(rows):
    atlas={f['family']:f for f in cert.read(spec.ATLAS)['families']};values={};bindings={}
    for family in sorted(atlas):
        pool=[r for r in rows if r['family']==family];scores=[[0,0,0] for r in pool];f=atlas[family];model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
        if len(pool)!=1024:raise ArithmeticError('fixed1024 frame roster differs')
        for q in PRIMES:
            path=scan.old.PARENT/family/'trace-tables'/f'{q}.json';t=cert.read(path)
            if t['input']!={'family':family,'model_hash':digest(model),'prime':q}:raise ArithmeticError('short canonical model binding differs')
            bindings[str(path.relative_to(ROOT))]=cert.hashed(path)
            for row,v in zip(pool,scores):
                x=row['numerator']*pow(row['denominator'],-1,q)%q if row['denominator']%q else q
                if t['good'][x]:
                    a,b=formula.contributions(t['traces'][x],q);v[0]+=a;v[1]+=b;v[2]+=1
        for row,v in zip(pool,scores):
            if v[0]!=row['score_units'] or v[2]!=row['good_primes']:raise ArithmeticError('original562-prime scores differ')
            values[(family,row['retained_index'])]=v[1]
        print('PRODUCT SHORT FRAME',family,len(pool),flush=True)
    return values,bindings

def selection(rows,p):
    excluded={(r['family'],r['retained_index']) for r in p['excluded_addresses']};result={};roster={}
    for f in sorted({r['family'] for r in rows}):
        available=[r for r in rows if r['family']==f and (f,r['retained_index']) not in excluded];result[f]={}
        for arm,key in [('product','product_selection_units'),('original_s1','combined_selection_units')]:
            chosen=sorted(available,key=lambda r:(-r[key],-r['combined_good'],r['denominator'],r['numerator']))[:2];result[f][arm]=[r['retained_index'] for r in chosen]
            for r in chosen:
                k=(f,r['retained_index'])
                if k not in roster:roster[k]={key:r[key] for key in ('family','retained_index','parameter','numerator','denominator','product_selection_units','combined_selection_units','combined_good')};roster[k]['arms']=[]
                roster[k]['arms'].append(arm)
    candidates=[roster[k] for k in sorted(roster)]
    if not 12<=len(candidates)<=24:raise ArithmeticError('paired candidate count differs')
    return result,candidates

def compute(check=False):
    p=protocol();out=D/'result.json'
    if not check and out.exists():raise FileExistsError('preserve paired score results')
    source=cert.read(parent.D/'result.json');short,bindings=short_scores(source['rows']);result={'schema':'elliptic-curves.higher-r17-product-score-result.v1','status':'RUNNING','protocol_hash':digest(p),'short_table_hashes':bindings,'rows':[]}
    for row in source['rows']:
        path=parent.D/row['family']/f"candidate-{row['retained_index']:04}"/'raw.json';raw=cert.read(path)
        if cert.hashed(path)!=row['raw_sha256'] or raw['program']!=parent.engine.program(row['model']) or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('saved trace binding or status differs')
        traces,ms=parent.engine.parse(raw['stdout'],row['model']);s=parent.scoring.sums(traces)
        if any(row[k]!=v for k,v in s.items()) or row['combined_selection_units']!=row['score_units']+s['extension_selection_units'] or row['combined_good']!=row['good_primes']+s['extension_good']:raise ArithmeticError('original extension scores differ')
        low=high=0
        for q,t in traces:
            if t is None:continue
            v=round(log((q+1-t)/q)*10**12)
            if q<=32749:low+=v
            else:high+=v
        result['rows'].append({k:row[k] for k in ('family','retained_index','parameter','numerator','denominator','combined_selection_units','combined_good','validation_units','validation_good','raw_sha256')}|{'product_short_units':short[(row['family'],row['retained_index'])],'product_selection_units':short[(row['family'],row['retained_index'])]+low,'product_validation_units':high})
        if len(result['rows'])%128==0 and not check:checkpoint(out,result)
        if len(result['rows'])%512==0:print('PRODUCT EXTENDED',len(result['rows']),'of6144',flush=True)
    result['selection'],result['prospective_candidates']=selection(result['rows'],p);result['status']='COMPLETE_FIXED_PAIRED_SCORE_COMPARISON'
    if check:
        if cert.read(out)!=result:raise ArithmeticError('full6144 score comparison replay differs')
    else:checkpoint(out,result)
    print('PRODUCT COMPARISON PASS',len(result['prospective_candidates']),'distinct paired candidates',result['selection'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();prepare() if a.stage=='prepare' else compute(a.stage=='replay')

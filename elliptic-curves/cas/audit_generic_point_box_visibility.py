#!/usr/bin/env python3
"""Exact retrospective generic-point visibility in two completed24-curve cohorts."""
import argparse
from pathlib import Path
from collections import Counter
from fractions import Fraction
import certify_compact_r17_candidates as cert
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/generic-point-box-visibility-v1';OUT=ART/'generic_point_box_visibility_v1.json';INPUTS=[ART/'retention24_r17_results_v1.json',ART/'higher24_r17_results_v1.json']
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),ROOT/'elliptic-curves/cas/search_observability.py',*INPUTS]}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve generic visibility protocol')
    for p in INPUTS:
        if len(cert.read(p)['curves'])!=24:raise ArithmeticError('fixed24 cohort required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.generic-point-box-visibility.v1','sources':sources(),'cohorts':[str(p.relative_to(ROOT)) for p in INPUTS],'scope':'All17 original generic sections, both signs, on every completed initial chart of both fixed24 cohorts. Check exact chart identities and primitive coordinates, including tangent/endpoints and infinity. Compare actual recorded points with declared completed boxes. No new point, trace, centre or search.','endpoint':'Retrospective conditional visibility of known points already used in centre construction; this is not a masked control or calibration of exceptional-point visibility. Any within-completed-box omission remains an explicit audit discrepancy pending investigation.','bounds':{'maximum_curves':48,'maximum_charts':2160,'maximum_point_chart_observations':73440,'wall_seconds':300,'rss_bytes':1610612736,'workers':1}})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen generic visibility inputs differ')
    return p

def compute(check=False):
    p=protocol();out=D/'result.json'
    if not check and (out.exists() or OUT.exists()):raise FileExistsError('preserve exact visibility audit')
    d={'schema':'elliptic-curves.generic-point-box-visibility-result.v1','status':'RUNNING','protocol_hash':digest(p),'rows':[],'sources':sources()};total=0
    for input_path in INPUTS:
        for row in cert.read(input_path)['curves']:
            path=ROOT/row['discovery_witness']['path'];raw=cert.read(path)
            if cert.hashed(path)!=row['discovery_witness']['sha256'] or len(raw['generic_points'])!=17:raise ArithmeticError('fixed witness or generic roster differs')
            d['sources'][str(path.relative_to(ROOT))]=cert.hashed(path);best=[None]*17;seen=set();counts=Counter();discrepancies=[];finite_minima=[];observation_hashes=[]
            for i,c in enumerate(raw['charts']):
                r=c['search']
                if r['status']!='bounded_search_complete' or r['height_bound']!=125000 or r['timeout_seconds']!=10:raise ArithmeticError('identical completed point budgets required')
                record={**r,'completed_denominator':r['height_bound']}
                for k,P in enumerate(raw['generic_points']):
                    x,y=map(cert.F,P)
                    for sign in (-1,1):
                        v=point_visibility(record,(x,sign*y));entry={'section':k,'sign':sign,'chart':i,'visibility':v};observation_hashes.append(digest(entry));counts[v['status']]+=1;total+=1
                        if v['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(entry)
                        if v.get('in_completed_box') or v['status']=='KNOWN_POINTED_ENDPOINT':seen.add(k)
                        H=v.get('minimum_affine_height')
                        key=(0 if v.get('at_parameter_infinity') or v['status']=='KNOWN_POINTED_ENDPOINT' else H if H is not None else 10**1000,i,sign)
                        if best[k] is None or key<best[k][0]:best[k]=(key,entry)
            minima=[b[1] for b in best]
            if any(b is None for b in best):raise ArithmeticError('missing generic-section observation')
            result={'cohort':input_path.stem,'id':row['id'],'family':row['family'],'parameter':row['parameter'],'initial_rank_lower_bound':row['rank_lower_bound'],'charts':len(raw['charts']),'point_chart_observations':len(observation_hashes),'ordered_observation_digest':digest(observation_hashes),'status_counts':dict(counts),'generic_sections_with_a_visible_sign':len(seen),'section_minima':minima,'discrepancies':discrepancies};d['rows'].append(result)
            if not check:checkpoint(out,d)
            print('GENERIC VISIBILITY',input_path.stem,row['id'],len(seen),'of17',len(discrepancies),'discrepancies',flush=True)
    if len(d['rows'])!=48 or total!=73440:raise ArithmeticError('fixed48/73440 visibility product differs')
    d['observations_checked']=total;d['status']='COMPLETE_EXACT_VISIBILITY_AUDIT';d['discrepancy_count']=sum(len(r['discrepancies']) for r in d['rows']);d['cohort_summary']={name:{'visible_section_counts':{str(k):v for k,v in sorted(Counter(r['generic_sections_with_a_visible_sign'] for r in d['rows'] if r['cohort']==name).items())},'discrepancies':sum(len(r['discrepancies']) for r in d['rows'] if r['cohort']==name)} for name in sorted({r['cohort'] for r in d['rows']})};d['claim_boundary']=p['endpoint']
    if check:
        if cert.read(out)!=d or cert.read(OUT)!=d:raise ArithmeticError('exact generic visibility replay differs')
    else:checkpoint(out,d);checkpoint(OUT,d)
    print('GENERIC VISIBILITY ALL73440',d['cohort_summary'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();prepare() if a.stage=='prepare' else compute(a.stage=='replay')

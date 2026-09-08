#!/usr/bin/env python3
"""Freeze and replay matched score strata from existing corrected MW16 survivors."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
from fractions import Fraction
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'retained-score-stratification-v1'
SCAN=LOCAL/'corrected-mw16-higher-annuli-v1'
ORIGINAL=ART/'corrected_mw16_higher_selection_v1.json'
OUT=ART/'retained_mw16_score_strata_selection_v1.json'
SEED='user-retained-score-strata-family-height-v1'
ARMS=('top','moderate','lower')

def digest(text):return hashlib.sha256(text.encode()).hexdigest()
def sources():
    paths=[Path(__file__).resolve(),Path(cert.__file__).resolve(),spec.ATLAS,
           SCAN/'protocol.json',SCAN/'result.json',SCAN/'replay.json',ORIGINAL,
           D/'cancelled-queued-r17-sweep.json']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve stratification design')
    if cert.read(SCAN/'replay.json')['status']!='PASS':raise ArithmeticError('completed retained pool required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-score-strata-protocol.v1',
        'sources':sources(),'seed':SEED,'blocks':20,'triplets_per_block':1,'curves_per_arm':20,
        'strata':{'top':[1,128,128],'moderate':[1025,8192,1024],'lower':[16385,65536,2048]},
        'block_fields':['family','band','sign'],'block_retained_count':65536,
        'j_height_ratio_maximum':4,'parameter_height_ratio_maximum':2,
        'selection_wall_seconds':1800,'selection_rss_bytes':2147483648,
        'point_exposure':{'generic_rank':16,'charts_per_curve':43,'height':125000,
                          'seconds_per_chart':10,'seconds_per_curve':600,'workers':2,
                          'rank_stop':None,'maximum_boxes':2580},
        'prerequisites':['corrected60-mw16-pari-v1/post-batch/ledger.json',
                         'corrected60_mw16_point_portable_replay_v1.json'],
        'score':'Corrected combined score through32749 for every retained candidate; order by descending score, good count, denominator, signed numerator within each block. No65521 score or withheld validation enters strata.',
        'matching':'Fixed SHA256 samples within the declared rank windows; seeded top-anchor order, then closest exact rational j-height ratio with seeded tie break. Enforce pairwise j-height ratio<=4 and parameter-height ratio<=2. Exclude Q-isomorphs of the original60 prospective models and previously chosen curves. One complete triplet per family/band/sign; no relaxed caliper or enlarged pool on failure.',
        'exposure_order':'Hash-order the20 blocks, cycle through the six arm permutations across blocks, and dispatch each matched triplet in that order with two identical workers.',
        'scope':'After the original corrected experiment and its182 isolated point-proof stages finish unchanged, select60 matched candidates from existing retained scores. No new parameter enumeration, new scalar traces, point outcomes, catalogue or validation-prime input. This file selects and replays a roster; it launches no point search. Rank gains, actual computation, completion and censoring are measured by the subsequent fixed2580-box experiment.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['seed']!=SEED:raise ArithmeticError('frozen stratification inputs changed')
    return p

def completion_gate():
    a=cert.read(LOCAL/'corrected60-mw16-pari-v1/post-batch/ledger.json')
    b=cert.read(ART/'corrected60_mw16_point_portable_replay_v1.json')
    if a['status']!='PASS' or b['status']!='PASS' or b['logical_stages']!=182:
        raise ArithmeticError('original corrected experiment must finish unchanged first')

def model(f,row):
    n,d=row['numerator'],row['denominator']
    coefficients=[]
    for key,w in [('A_coefficients_low_to_high',8),('B_coefficients_low_to_high',12)]:
        value=sum(Fraction(c)*n**i*d**(w-i) for i,c in enumerate(f[key]))
        if value.denominator!=1:raise ArithmeticError('integral homogeneous model required')
        coefficients.append(value.numerator)
    a,b=coefficients;delta=-16*(4*a**3+27*b**2)
    if not delta:raise ArithmeticError('singular retained candidate')
    j=Fraction((-48*a)**3,delta)
    return [0,0,0,a,b],j,max(abs(j.numerator),j.denominator)

def compatible(rows,p):
    for key,bound in [('j_height',p['j_height_ratio_maximum']),('parameter_height',p['parameter_height_ratio_maximum'])]:
        values=[int(r[key]) for r in rows]
        if max(values)>bound*min(values):return False
    return True

def ratio(a,b):return Fraction(max(a,b),min(a,b))

def expected():
    p=protocol();completion_gate()
    population=cert.read(SCAN/'result.json')
    if population['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or len(population['shards'])!=320:
        raise ArithmeticError('all320 saved shards required')
    families={r['fibration_id']:r for r in cert.read(spec.ATLAS)['families']}
    excluded={}
    for row in cert.read(ORIGINAL)['selected']:
        m=list(map(Fraction,row['model']));v=cert.weierstrass_invariants(m);j=v['c4']**3/v['discriminant']
        excluded.setdefault(j,[]).append(m)
    blocks=list(itertools.product(sorted(families),(2,3),(-1,1)))
    blocks.sort(key=lambda b:digest(SEED+'|block|'+repr(b)))
    records=[];failures=[];execution=[];permutations=list(itertools.permutations(ARMS))
    for block_index,(family,band,sign) in enumerate(blocks):
        block=f'{family}-b{band}-s{sign}';pool=[];files=[]
        for manifest in population['shards']:
            if (manifest['family'],manifest['band'],manifest['sign'])!=(family,band,sign):continue
            path=ROOT/manifest['retained_path']
            if cert.hashed(path)!=manifest['retained_sha256']:raise ArithmeticError('retained shard changed')
            pool.extend({**r,'slice_id':manifest['id']} for r in cert.read(path)['rows']);files.append(manifest['retained_sha256'])
        if len(files)!=16 or len(pool)!=65536 or len({r['parameter'] for r in pool})!=65536:
            raise ArithmeticError('complete signed retained block required')
        pool.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        arms={};sampling={}
        for arm,(lo,hi,size) in p['strata'].items():
            chosen=[{**r,'retained_rank':i+1,'sampling_hash':digest(SEED+'|'+block+'|'+arm+'|'+r['parameter'])}
                    for i,r in enumerate(pool) if lo<=i+1<=hi]
            chosen.sort(key=lambda r:r['sampling_hash']);chosen=chosen[:size]
            sampling[arm]={'sample_size':len(chosen),'parameters_sha256':digest(json.dumps([(r['parameter'],r['retained_rank']) for r in chosen],separators=(',',':')))}
            arms[arm]=[]
            for r in chosen:
                m,j,h=model(families[family],r)
                if any(cert.isomorphic(m,old) for old in excluded.get(j,[])):continue
                arms[arm].append({**r,'arm':arm,'block':block,'family':family,'band':band,'sign':sign,
                    'model':list(map(str,m)),'j_invariant':str(j),'j_height':str(h),
                    'parameter_height':max(abs(r['numerator']),r['denominator'])})
        selected=None
        for top in arms['top']:
            close=lambda r:(ratio(int(top['j_height']),int(r['j_height'])),r['sampling_hash'])
            moderate=sorted((r for r in arms['moderate'] if compatible([top,r],p)),key=close)
            lower=sorted((r for r in arms['lower'] if compatible([top,r],p)),key=close)
            for mid in moderate:
                for low in lower:
                    triple=[top,mid,low]
                    if not compatible(triple,p):continue
                    if any(cert.isomorphic(a['model'],b['model']) for a,b in itertools.combinations(triple,2)):continue
                    selected=triple;break
                if selected:break
            if selected:break
        if selected is None:
            failures.append({'block':block,'reason':'NO_MATCH_IN_FIXED_SAMPLED_POOLS','sampling':sampling});continue
        for r in selected:
            r['id']='strata-'+block+'-'+r['arm']
            excluded.setdefault(Fraction(r['j_invariant']),[]).append(list(map(Fraction,r['model'])))
        by_arm={r['arm']:r for r in selected}
        order=permutations[block_index%6]
        execution.extend(by_arm[a] for a in order)
        records.append({'block':block,'sampling':sampling,'execution_arm_order':list(order),
                        'j_height_ratio':str(ratio(max(int(r['j_height']) for r in selected),min(int(r['j_height']) for r in selected))),
                        'parameter_height_ratio':str(Fraction(max(r['parameter_height'] for r in selected),min(r['parameter_height'] for r in selected)))})
        print('MATCHED SCORE TRIPLET',block,flush=True)
    return {'schema':'elliptic-curves.retained-score-strata-selection.v1',
        'status':'PASS_FROZEN60_MATCHED_SELECTION' if not failures and len(execution)==60 else 'MATCHING_INCOMPLETE_NO_POINT_SEARCH',
        'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),
        'selected':execution,'blocks':records,'matching_failures':failures,
        'claim_boundary':p['scope']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','select','check']);args=a.parse_args()
    if args.stage=='prepare':prepare()
    else:
        data=expected()
        if args.stage=='check':
            if cert.read(OUT)!=data:raise ArithmeticError('matched selection replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve matched roster')
            checkpoint(OUT,data)
        print('RETAINED SCORE STRATA',data['status'],flush=True)

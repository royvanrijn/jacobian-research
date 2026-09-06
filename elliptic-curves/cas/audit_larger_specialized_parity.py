#!/usr/bin/env python3
"""Fixed two-case geometry-only enlargement of the specialized parity sample."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/larger-specialized-parity-v1';OUT=ART/'larger_specialized_parity_v1.json'
CASES=[{'id':'native28','rank':28,'maps':'artifacts/local/elliptic-curves/native28-specialized-parity-control-v1/maps.json'},
       {'id':'new-20260906-186','rank':27,'maps':'artifacts/local/elliptic-curves/full11952-specialized-followup-v1/new-20260906-186/maps.json'}]
def sources():
    paths=[Path(__file__).resolve(),CAS/'prepare_larger_specialized_parity.sage',CAS/'prospective_half_lattice_v2.sage',*[ROOT/r['maps'] for r in CASES]]
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve two-case parity audit')
    for case in CASES:
        d=cert.read(ROOT/case['maps'])
        if d['status']!='COMPLETE_DECLARED_MAPS' or len(d['sample'])!=2048 or len(d['centres'])!=49 or len(d['rounded_gram'])!=case['rank']:raise ArithmeticError('completed original geometry required')
    checkpoint(D/'protocol.json',{'sources':sources(),'rows':CASES,'sample_size':65536,'original_sample_size':2048,
        'domain':'larger-specialized-parity-v1','generic_prefix':17,'selected_centres':49,
        'workers':2,'seconds_per_case':300,'replay_seconds':300,'rss_bytes':1610612736,
        'gate':'The known28 specialized control recovered29 with the old2048-mask policy, while the own27 ID186 stayed27 after the same policy and all49 million-height boxes. Sampling only2048 masks leaves most specialized parity classes untested. Measure a finite enlarged sample in the unchanged frozen metric before proposing any new point exposure.',
        'scope':'Reuse each original2048 sample exactly and append deterministic distinct masks to65536, with nonzero quotient above the17 generic bits. Keep the original rounded metric and unimodular reduction fixed. Compute numerical CVP representatives, exactly verify every parity and integer rounded norm, and compare the49 largest computed norms to the old49. No rational chart maps, point search, new parameter, covering/optimality theorem, rank claim or automatic further enlargement.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen larger parity inputs differ')
    return p
def masks(p,case):
    old=cert.read(ROOT/case['maps']);result=[r['parity'] for r in old['sample']];seen=set(result);i=0
    if len(seen)!=2048:raise ArithmeticError('original distinct masks required')
    while len(result)<p['sample_size']:
        m=int(digest([p['domain'],case['id'],i]),16)%(1<<case['rank']);i+=1
        if m>>p['generic_prefix'] and m not in seen:result.append(m);seen.add(m)
    return result
def expected():
    p=protocol();rows=[]
    for case in p['rows']:
        old=cert.read(ROOT/case['maps']);path=D/case['id']/'sample.json';d=cert.read(path)
        if d['status']!='COMPLETE_FIXED_SAMPLE' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or d['case']!=case or d['sample'][:2048]!=old['sample'] or [r['parity'] for r in d['sample']]!=masks(p,case):raise ArithmeticError('sample provenance differs')
        g=old['rounded_gram'];rank=case['rank']
        for r in d['sample']:
            word=r['representative']
            if len(word)!=rank or any(type(c)!=int or (c-(r['parity']>>j))%2 for j,c in enumerate(word)):raise ArithmeticError('exact parity differs')
            norm=sum(word[i]*g[i][j]*word[j] for i in range(rank) for j in range(rank))
            if norm!=r['metric_norm'] or norm<0:raise ArithmeticError('exact rounded norm differs')
        top=sorted(d['sample'],key=lambda r:(-r['metric_norm'],r['parity']))[:49]
        if d['centres']!=top:raise ArithmeticError('largest computed norm ordering differs')
        oldmasks={r['parity'] for r in old['centres']};newmasks={r['parity'] for r in top}
        rows.append({'id':case['id'],'dimension':rank,'sampled_masks':len(d['sample']),
            'old_maximum_computed_norm':old['centres'][0]['metric_norm'],'new_maximum_computed_norm':top[0]['metric_norm'],
            'old_49th_computed_norm':old['centres'][-1]['metric_norm'],'new_49th_computed_norm':top[-1]['metric_norm'],
            'retained_old_top49_masks':len(oldmasks&newmasks),'new_top49_masks':len(newmasks-oldmasks),
            'sample_sha256':cert.hashed(path)})
    return {'schema':'elliptic-curves.larger-specialized-parity.v1','status':'PASS_FIXED_GEOMETRY_AUDIT','sources':sources(),
        'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'claim_boundary':p['scope']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('larger parity summary differs')
        else:
            if OUT.exists():raise FileExistsError('preserve parity summary')
            checkpoint(OUT,d)
        print('LARGER PARITY SAMPLE',d['rows'],flush=True)

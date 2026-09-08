#!/usr/bin/env python3
"""Bind the three positive-ratio parent proofs and exact Q-distinctness."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compute():
    paths={Path(__file__).resolve()};stages=[]
    jobs=[('kihara-positive-parents-v1','backend_failure'),('kihara-positive-parents-v2','completed'),('kihara-positive-parents-standalone-v1','backend_failure'),('kihara-positive-parents-standalone-v2','completed')]
    for folder,outcome in jobs:
        d=LOCAL/folder;p=read(d/'protocol.json');s=read(d/'supervisor.json');assert s['outcome']==outcome
        log=Path(s['log']);assert sha(log)==s['log_sha256'];paths.update(d.glob('*.json'));paths.add(log)
        for name,h in p['sources'].items():f=ROOT/name;assert sha(f)==h;paths.add(f)
        stages.append({'stage':folder,'outcome':outcome,'wall_seconds':s['wall_seconds'],'supervisor':str((d/'supervisor.json').relative_to(ROOT))})
    for folder,name in [('kihara-positive-parents-standalone-v1','verify_kihara_positive_parents.sage'),('kihara-positive-parents-standalone-v2','verify_kihara_positive_parents_v2.sage')]:
        d=LOCAL/folder;canonical=ROOT/'elliptic-curves/cas'/name;assert (d/name).read_bytes()==canonical.read_bytes();paths.update([d/name,canonical])
        b=ART/'kihara_positive_parent_bundle_v1.json';assert (d/b.name).read_bytes()==b.read_bytes();paths.update([d/b.name,b])
    rows=[read(LOCAL/'kihara-positive-parents-standalone-v1/parent0.json')]+[read(LOCAL/'kihara-positive-parents-standalone-v2'/('parent'+str(i)+'.json')) for i in (1,2)]
    assert [r['ratio'] for r in rows]==['1','2','3'] and [r['generic_section_span_rank'] for r in rows]==[7,9,12]
    inputs=read(ART/'kihara_positive_parent_bundle_v1.json')['rows']
    for r,i in zip(rows,inputs):r['raw_A']=i['raw_A'];r['raw_B']=i['raw_B'];r['label']='kihara-ratio-'+r['ratio']
    old_k=ART/'kihara_five_parent_distinctness_v1.json';old_m=ART/'mestre_parent_portfolio_intake_v1.json';paths.update([old_k,old_m]);mk=read(old_m)
    old=[dict(r,label='kihara-path-'+r['path_parameter']) for r in read(old_k)['rows']]+[dict(r,label=r['id']) for r in mk['rows']+[mk['reference_parent']]]
    def fingerprint(row):return {s['prime']:s['surface_point_count'] for s in row['surface_counts'] if s.get('status','PASS')=='PASS'}
    pairs=[]
    for i,r in enumerate(rows):
        f=fingerprint(r)
        for other in old+rows[i+1:]:
            g=fingerprint(other);w=[p for p in (131,239,251) if p in f and p in g and f[p]!=g[p]];assert w,(r['label'],other['label'])
            pairs.append({'left':r['label'],'right':other['label'],'separating_good_primes':w,'status':'PROVED_Q_DISTINCT'})
    assert len(pairs)==39
    return {'schema':'kihara-positive-parent-expansion.v1','status':'PASS','rows':rows,'pairwise_separations':pairs,'new_Q_distinct_parents':3,'total_proved_Q_distinct_parents_at_least':14,
        'point_search_boxes':0,'inventory_additions':0,'stages':stages,'total_supervised_seconds':sum(s['wall_seconds'] for s in stages),'sources':{str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)},
        'scope':'Three positive-ratio Kihara parents outside the retained rank14 path are Q-nonisomorphic to every prior retained parent and to one another. New surface counts are independently computed at131,239,251; prior counts reuse their existing certificates. Direct generic section identities and exact height Grams prove spans7,9,12. These span ranks are not full generic-rank upper bounds. The configurations are8I1+6I2+I4,14I1+2I2+I6,20I1+I4. Full NS lattices, geometric parent inequivalence and literature novelty remain UNKNOWN. The first standalone run completed ratio1 before discovering I6 at ratio2; V2 checks only the remaining rows with the corrected component-killing multiplier. All failures are charged and all completed inputs preserved. No parameter ranking, point-search exposure or near-record curve was produced.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=compute();out=ART/'kihara_positive_parent_expansion_v1.json'
    if a.check:assert r==read(out)
    else:
        with out.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print('PASS three new Q-distinct parents; section spans7,9,12;',r['total_supervised_seconds'],'supervised seconds')

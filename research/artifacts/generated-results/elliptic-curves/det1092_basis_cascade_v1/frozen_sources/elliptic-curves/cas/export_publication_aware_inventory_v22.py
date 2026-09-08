#!/usr/bin/env python3
"""Preserve discovery IDs and make current catalogue/public-point provenance explicit."""
import argparse,csv,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import memory_rank_certificate as memory
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OLD=ART/'new_high_rank_curve_index_v21.json'
OLD_REPLAY=ART/'new_high_rank_curve_index_v21_memory_replay_v1.json'
PUBLIC=ART/'inventory188_public28_reproduction_v1.json'
PUBLIC_REPLAY=ART/'inventory188_public28_sage_replay_v1.json'
FRESH=ART/'inventory200_current_catalogue_comparison_v1.json'
FRESH_REPLAY=ART/'inventory200_current_catalogue_sage_replay_v1.json'
COHORT_REPLAY=ART/'nearcut60_current_catalogue_sage_v1.json'
WRAPPER=ART/'inventory188_public28_point_export_v1.json'
OUT=ART/'new_high_rank_curve_index_v22.json'
CHECK=ART/'new_high_rank_curve_index_v22_replay_v1.json'


def bound(path):
    d=cert.read(path)
    if d.get('status')!='PASS':raise ArithmeticError('completed input required: '+str(path))
    for n,h in d['sources'].items():
        if cert.hashed(ROOT/n)!=h:raise ArithmeticError('input changed: '+n)
    return d


def expected_wrapper():
    public=bound(PUBLIC);bound(PUBLIC_REPLAY)
    row={k:public[k] for k in ('family','parameter','curve','rank_certificate','rank_lower_bound')}
    row.update(id=public['local_id'],points=public['independent_points'],
               rank_provenance='PUBLIC_POINT_REPRODUCTION',catalogue_id=619)
    return {'schema':'elliptic-curves.public28-point-export.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),PUBLIC,PUBLIC_REPLAY)},
            'curves':[row],'claim_boundary':public['claim_boundary']}


def expected():
    old=cert.read(OLD);replay=bound(OLD_REPLAY);public=bound(PUBLIC)
    fresh=bound(FRESH);bound(FRESH_REPLAY);cohort=bound(COHORT_REPLAY)
    if replay['curves_checked']!=201 or not replay['csv_checked'] or len(old['curves'])!=201:
        raise ArithmeticError('complete201-curve prior replay required')
    if cert.read(WRAPPER)!=expected_wrapper():raise ArithmeticError('public point export differs')
    raw=ROOT/'artifacts/local/elliptic-curves/inventory200-current-catalogue-v1/database.json'
    if cert.hashed(raw)!=fresh['catalogue_intake']['sha256']:raise ArithmeticError('current catalogue changed')
    equations=[{k:r[k] for k in ('id','ainvs','rank_lower_bound','conductor')} for r in cert.read(raw)['curves']]
    byj={}
    for q in equations:
        v=cert.weierstrass_invariants(tuple(map(cert.F,q['ainvs'])))
        byj.setdefault(v['c4']**3/v['discriminant'],[]).append(q)
    known={r['id']:r['q_isomorphism_matches'] for r in fresh['inventory_comparisons']}
    cohort_known={r['id']:r['q_isomorphism_matches'] for r in cohort['comparisons']}
    rows=[]
    for original in old['curves']:
        r=dict(original);r['local_search_rank_lower_bound']=r['rank_lower_bound']
        r['rank_provenance']='LOCAL_POINT_SEARCH'
        if r['id']==public['local_id']:
            r['original_local_source']={k:r[k] for k in ('source_certificate','source_curve_index','rank_lower_bound')}
            r.update(points=public['independent_points'],rank_certificate=public['rank_certificate'],
                     rank_lower_bound=public['rank_lower_bound'],rank_provenance='PUBLIC_POINT_REPRODUCTION',
                     source_certificate=WRAPPER.name,source_curve_index=0)
        model=tuple(map(cert.F,r['curve']));v=cert.weierstrass_invariants(model);j=v['c4']**3/v['discriminant']
        if str(j)!=r['j_invariant']:raise ArithmeticError('inventory j differs')
        matches=[q['id'] for q in byj.get(j,[]) if cert.isomorphic(model,q['ainvs'])]
        if r['id'] in known:expected_matches=known[r['id']]
        else:
            incoming=cert.read(ART/original['source_certificate'])['curves'][original['source_curve_index']]
            expected_matches=cohort_known[incoming['id']]
        if matches!=expected_matches:raise ArithmeticError('independent current-catalogue comparison differs')
        r['current_catalogue_matches']=matches
        r['novelty_status']='CATALOGUE_MATCH' if matches else 'UNMATCHED_IN_PINNED620'
        rows.append(r)
    rows.sort(key=lambda r:(-r['rank_lower_bound'],r['family'],r['parameter']))
    paths=[Path(__file__).resolve(),Path(cert.__file__),Path(memory.__file__),OLD,OLD_REPLAY,PUBLIC,PUBLIC_REPLAY,FRESH,FRESH_REPLAY,COHORT_REPLAY,WRAPPER,raw]
    sources={str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}
    names=set(old['source_certificate_hashes'])|{WRAPPER.name}
    for n,h in old['source_certificate_hashes'].items():
        if cert.hashed(ART/n)!=h:raise ArithmeticError('historical point source changed')
    return {'schema':'elliptic-curves.new-high-rank-index.v22','sources':sources,
            'previous_index_sha256':cert.hashed(OLD),'curves':rows,
            'source_certificate_hashes':{n:cert.hashed(ART/n) for n in sorted(names)},
            'catalogue':{'url':'https://elliptic-rank.icarm.cloud/database.json','curve_count':620,
                         'raw_sha256':fresh['catalogue_intake']['sha256'],'equations':equations},
            'unmatched_curve_count':sum(not r['current_catalogue_matches'] for r in rows),
            'catalogue_matched_ids':[r['id'] for r in rows if r['current_catalogue_matches']],
            'claim_boundary':'Research inventory with stable discovery IDs and explicit current publication status. All201 previous point bases were replayed in V21; the only changed rank/basis is ID188, now28 by independently reproduced public points. Its original local search bound remains27. ID12 and ID188 are current catalogue matches;199 equations remain unmatched in the pinned620 snapshot. No new rank28 discovery, exact rank, first-discovery priority or universal novelty.'}


def csv_rows(d):
    keys=['id','rank_lower_bound','local_search_rank_lower_bound','rank_provenance','novelty_status','family','parameter']
    header=keys+['a1','a2','a3','a4','a6','current_catalogue_matches','source_certificate','source_curve_index']
    return [header]+[[str(v) for v in [*(r[k] for k in keys),*r['curve'],','.join(map(str,r['current_catalogue_matches'])),r['source_certificate'],r['source_curve_index']]] for r in d['curves']]


def build():
    if any(p.exists() for p in (WRAPPER,OUT,OUT.with_suffix('.csv'),CHECK)):raise FileExistsError('preserve publication-aware inventory')
    checkpoint(WRAPPER,expected_wrapper());d=expected();checkpoint(OUT,d)
    with OUT.with_suffix('.csv').open('x',newline='') as f:csv.writer(f).writerows(csv_rows(d))


def check():
    d=expected()
    if cert.read(OUT)!=d:raise ArithmeticError('publication-aware inventory extraction differs')
    if list(csv.reader(OUT.with_suffix('.csv').open(newline='')))!=csv_rows(d):raise ArithmeticError('publication-aware CSV differs')
    old={r['id']:r for r in cert.read(OLD)['curves']};changed=[]
    for r in d['curves']:
        before=old[r['id']]
        if any(r[k]!=before[k] for k in ('curve','points','rank_certificate','rank_lower_bound')):
            changed.append(r['id']);proof=r['rank_certificate'];model=tuple(map(cert.F,r['curve']));points=[tuple(map(cert.F,p)) for p in r['points']]
            actual=memory.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
            if json.dumps(actual,sort_keys=True)!=json.dumps(proof,sort_keys=True):raise ArithmeticError('changed28-point proof differs')
    if changed!=['new-20260906-188'] or d['unmatched_curve_count']!=199:raise ArithmeticError('unexpected changed proof or publication count')
    evidence={'schema':'elliptic-curves.inventory-v22-provenance-replay.v1','status':'PASS',
              'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),OUT,OUT.with_suffix('.csv'),OLD_REPLAY,PUBLIC_REPLAY)},
              'curves_checked':201,'unchanged_point_bases_bound_to_full_v21_replay':200,
              'changed_basis_independently_replayed_ids':changed,'csv_checked':True,
              'current_catalogue_count':620,'unmatched_curves':199,
              'claim_boundary':'Every unchanged basis is identical to its hash-bound fully replayed V21 source. The changed public28 basis is rechecked by exact finite quotients and its separate Sage group-enumeration certificate. Current catalogue comparisons and every CSV entry are checked. No new search or claim that a published curve is a new discovery.'}
    if CHECK.exists():
        if cert.read(CHECK)!=evidence:raise ArithmeticError('inventory replay evidence differs')
    else:checkpoint(CHECK,evidence)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
    check() if a.check else build()
    print('PUBLICATION-AWARE INVENTORY201;199 UNMATCHED; PUBLIC28 EXPLICIT')

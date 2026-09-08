#!/usr/bin/env python3
"""Fixed point exposure in49 previously frozen own27 charts of the public28 control."""
import argparse
import sys
from pathlib import Path
import full11952_specialized_followup as engine
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits

ROOT, CAS, ART = engine.ROOT, engine.CAS, engine.ART
OLD = ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1'
D = ROOT/'artifacts/local/elliptic-curves/inventory188-fixed49-point-control-v1'
OUT = ART/'inventory188_fixed49_point_control_v1.json'
MOD2 = ART/'inventory188_fixed49_point_control_mod2_v1.json'
MODL = ART/'inventory188_fixed49_point_control_modl_v1.json'


def configure():
    engine.D = D
    engine.SEED = D/'seed.json'
    engine.ROW = {'initial_rank': 27}
    engine.protocol = protocol


def source_hashes():
    paths = [Path(__file__).resolve(), OLD/'seed.json', OLD/'maps.json',
             CAS/'full11952_specialized_followup.py', CAS/'half_lattice_pointed_sieve.py',
             CAS/'pointed_quartic_search.py', CAS/'audit_recorded_point_mod2_rank_v3.py',
             CAS/'audit_retained_cloud_modl.py']
    return {**engine.sources(), **{str(p.relative_to(ROOT)): engine.cert.hashed(p) for p in paths}}


def protocol():
    p = engine.cert.read(D/'protocol.json')
    assert p['sources'] == source_hashes(), 'fixed control sources changed'
    assert engine.cert.hashed(D/'seed.json') == p['seed_sha256']
    original = engine.cert.read(OLD/'maps.json')
    original['protocol_hash'] = digest(p)
    assert engine.cert.read(D/'maps.json') == original, 'only metadata binding may differ'
    return p


def freeze():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve fixed point control')
    seed = engine.cert.read(OLD/'seed.json')
    assert len(seed['points']) == 27 and len(seed['generic_points']) == 17
    assert seed['parameter'] == '110314/102227'
    assert seed['points'][:17] == seed['generic_points']
    proof = seed['rank_certificate']
    assert digest(engine.checked_rank(tuple(map(engine.cert.F,seed['curve'])),
        [tuple(map(engine.cert.F,P)) for P in seed['points']],
        [r['prime'] for r in proof['signatures']], proof['no_rational_2_torsion_prime'])) == digest(proof)
    maps = engine.cert.read(OLD/'maps.json')
    assert maps['status'] == 'COMPLETE_DECLARED_MAPS' and len(maps['rows']) == 49
    checkpoint(D/'seed.json',seed)
    p = {'schema':'elliptic-curves.inventory188-fixed49-point-control.v1',
         'sources':source_hashes(), 'seed_sha256':engine.cert.hashed(D/'seed.json'),
         'sample_size':2048,'sample_domain':'full11952-specialized-followup-v1',
         'charts':49,'height':125000,'seconds_per_chart':10,'target_rank':28,
         'rank_stop':False,'gp_sha256':engine.cert.hashed(Path('/usr/bin/gp')),
         'worker_seconds':900,'replay_seconds':600,'rss_bytes':1610612736,'maximum_workers':1,
         'gate':'A retrospectively selected known28 curve has only27 old search points. Direct and fixed translated public representatives remain outside the saved own27 boxes, but that does not measure recovery of other representatives. Execute these49 already fixed maps exactly once using only the old27 seed, to measure actual quotient recovery.',
         'scope':'Known-control detector calibration, not a blind curve discovery or a new parameter population. All49 charts receive125000 height and ten seconds each, without rank stopping, new maps, adaptive centres or following search. Public points are absent from the point worker input.',
         'following_campaign':None}
    checkpoint(D/'protocol.json',p)
    maps['protocol_hash']=digest(p);checkpoint(D/'maps.json',maps)
    print('FROZEN49 existing maps / initial27 / no rank stop')


def worker():
    p=protocol();maps=engine.cert.read(D/'maps.json');out=D/'result.json'
    if out.exists():raise FileExistsError('preserve49 point attempt')
    cache=engine.ReductionCache(engine.MemoryFactStore());seed,state=engine.initial(cache)
    model=tuple(map(engine.cert.F,seed['curve']))
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')}
    data.update(protocol_hash=digest(p),maps_sha256=engine.cert.hashed(D/'maps.json'),
        initial_state=state.record(),initial_dimension=27,centres=maps['centres'],metric_gram=maps['metric_gram'],
        charts=[],status='RUNNING',rank_lower_bound=27,final_state=state.record(),arithmetic_facts=cache.store.snapshot())
    checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=engine.rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive)
        rep=m['centre']['representative']+[0]*(state.rank-27)
        search=engine.PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy'])
        r,points=engine.backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256'])
        compression=engine.compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record()
        data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),
            'archive_sha256':engine.cert.hashed(ap),'search':r,'admission_compression':compression,
            'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank})
        data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot())
        checkpoint(out,data);print('CONTROL49',i+1,r['status'],'rank',state.rank,flush=True)
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)


def jobs(p):
    return [('points',[sys.executable,str(Path(__file__).resolve()),'worker'],p['worker_seconds']),
        ('replay',[sys.executable,str(Path(__file__).resolve()),'replay'],p['replay_seconds']),
        ('mod2-build',[sys.executable,str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--input',str(D/'result.json'),
                      '--input-sha256',None,'--output',str(MOD2),'--prime-bound','997'],300),
        ('mod2-check',[sys.executable,str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(MOD2)],300),
        ('modl-build',[sys.executable,str(CAS/'audit_retained_cloud_modl.py'),'--input',str(MOD2),'--output',str(MODL)],300),
        ('modl-check',[sys.executable,str(CAS/'audit_retained_cloud_modl.py'),'--check',str(MODL)],300)]


def launch():
    p=protocol();ledger_path=D/'ledger.json'
    if ledger_path.exists():raise FileExistsError('preserve fixed control ledger')
    ledger={'status':'RUNNING','stages':[]};checkpoint(ledger_path,ledger)
    for name,argv,seconds in jobs(p):
        argv=[engine.cert.hashed(D/'result.json') if x is None else x for x in argv]
        s=run(argv,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),
              checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['stages'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
        checkpoint(ledger_path,ledger);print(name,s['outcome'],s['returncode'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(ledger_path,ledger);return
    ledger['status']='PASS';checkpoint(ledger_path,ledger)


def report(check=False):
    p=protocol();ledger=engine.cert.read(D/'ledger.json');data=engine.cert.read(D/'result.json')
    assert ledger['status']=='PASS' and [r['name'] for r in ledger['stages']]==[j[0] for j in jobs(p)]
    assert len(data['charts'])==49 and data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT'
    for r in ledger['stages']:
        s=r['supervision'];assert r['status']=='PASS' and s['outcome']=='completed' and s['returncode']==0
        assert engine.cert.hashed(D/(r['name']+'.log'))==s['log_sha256']
    mod2=engine.cert.read(MOD2);modl=engine.cert.read(MODL)
    paths=[Path(__file__).resolve(),D/'protocol.json',D/'ledger.json',D/'result.json',D/'maps.json',MOD2,MODL]
    r={'schema':'elliptic-curves.inventory188-fixed49-outcome.v1','status':'PASS',
       'sources':{str(q.relative_to(ROOT)):engine.cert.hashed(q) for q in paths},
       'attempted_charts':49,'completed_charts':sum(c['search']['status']=='bounded_search_complete' for c in data['charts']),
       'initial_rank_lower_bound':27,'rank_lower_bound':mod2['rank_lower_bound'],
       'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in modl['audits']},
       'retained_points':len(mod2['points']),'point_worker_seconds':ledger['stages'][0]['supervision']['wall_seconds'],
       'total_supervised_seconds':sum(x['supervision']['wall_seconds'] for x in ledger['stages']),
       'claim_boundary':p['scope']}
    assert r['rank_lower_bound']==data['rank_lower_bound']
    if check:assert r==engine.cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve control report')
        checkpoint(OUT,r)
    print(r)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['freeze','launch','worker','replay','report','check'])
    args=parser.parse_args();configure()
    if args.mode=='replay':engine.replay()
    elif args.mode=='check':report(True)
    else:globals()[args.mode]()

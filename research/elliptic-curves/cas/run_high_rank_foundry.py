#!/usr/bin/env python3
"""Prepare, launch, inspect, stop or resume a detached elliptic-curve foundry."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import fcntl
from fractions import Fraction as F
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import sqlite3
import subprocess
import sys
import time
import traceback
import zipfile

from v3_warm_support import read, sha, atomic as _atomic, same_process, process_info, require
from high_rank_foundry_policy import DEFAULTS, FAMILIES, apply_result, choose_lane, utility, next_bank, effective_batch_calls
from high_rank_foundry_intake import Intake, jkey, matches, supported

CAS=Path(__file__).resolve().parent
ROOT=CAS.parents[1]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
DEFAULT_FOLDER=LOCAL/'high-rank-foundry-v3'


def save(path,data,immutable=False):
    _atomic(path,json.loads(json.dumps(data)),immutable=immutable)


def exclusive(path):
    path.parent.mkdir(parents=True,exist_ok=True)
    fd=os.open(path,os.O_CREAT|os.O_RDWR,0o600)
    try: fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BaseException:
        os.close(fd);raise RuntimeError('another process owns '+str(path))
    return fd


def connect(folder):
    db=sqlite3.connect(folder/'ledger.sqlite',timeout=30)
    db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=FULL')
    db.executescript('CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);'
                     'CREATE TABLE IF NOT EXISTS curves (id TEXT PRIMARY KEY, j TEXT NOT NULL, record TEXT NOT NULL);'
                     'CREATE TABLE IF NOT EXISTS jobs (id INTEGER PRIMARY KEY, cid TEXT NOT NULL, state TEXT NOT NULL, record TEXT NOT NULL);'
                     'CREATE INDEX IF NOT EXISTS jobs_state ON jobs(state);')
    return db


def meta(db,key,default=None):
    row=db.execute('SELECT value FROM meta WHERE key=?',(key,)).fetchone()
    return json.loads(row[0]) if row else default


def putmeta(db,key,value):
    db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)',(key,json.dumps(value)))


def curves(db):
    return [json.loads(r[0]) for r in db.execute('SELECT record FROM curves ORDER BY id')]


def putcurve(db,c):
    db.execute('INSERT OR REPLACE INTO curves VALUES (?,?,?)',(c['id'],c['j_key'],json.dumps(c)))


def putjob(db,j):
    db.execute('INSERT OR REPLACE INTO jobs VALUES (?,?,?,?)',(j['id'],j['cid'],j['state'],json.dumps(j)))


def jobs(db,states=None):
    if states is None:
        query,args='SELECT record FROM jobs ORDER BY id',()
    else:
        query='SELECT record FROM jobs WHERE state IN ('+','.join('?' for _ in states)+') ORDER BY id'
        args=tuple(states)
    return [json.loads(r[0]) for r in db.execute(query,args)]


def guard(folder):
    manifest=read(folder/'manifest.json')
    frozen=Path(manifest['frozen_root'])
    require(sha(folder/'config.json')==manifest['config_sha256'],'foundry configuration changed')
    for name,digest in manifest['files'].items():
        require(sha(frozen/name)==digest,'frozen source/input changed: '+name)
    for name,digest in manifest['executables'].items():
        require(sha(Path(name))==digest,'arithmetic executable changed: '+name)
    return frozen


def fresh_record(row):
    return {**row,'rank':17,'state':'NEW','total_calls':0,'stale_calls':0,'batches':0,
            'gaining_batches':0,'last_batch_gain':0,'bank_index':0,'head':None,
            'packet':None,'packet_sha256':None,'history':[],'revivals':0}


def prepare(folder,args):
    require(not folder.exists(),'preserve existing foundry; use resume or a new folder')
    folder.mkdir(parents=True)
    frozen=folder/'runtime/research';frozen.mkdir(parents=True)
    # Copy source bytes, not symlinks into a moving checkout. No large raw searches.
    paths=[]
    for directory in (ROOT/'elliptic-curves/cas',ROOT/'elliptic-curves/ecsearch'):
        paths += [p for p in directory.rglob('*') if p.is_file() and p.suffix in ('.py','.sage','.gp','.c','.cpp','.h','.sh') and '__pycache__' not in p.parts]
    paths += list((ROOT/'elliptic-curves').glob('*.py'))
    for path in sorted(set(paths)):
        dest=frozen/path.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    inputs=(ART/'compact_six_r17_atlas_v1.json',ART/'r17_exact_maximum_parity_classes_v1.json')
    for path in inputs:
        dest=frozen/path.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    data=frozen/'foundry-inputs';data.mkdir()
    inherited_config={}
    if args.inherit_state:
        require(args.catalogue is None,'inherited campaigns retain their pinned catalogue')
        from high_rank_foundry_migration import freeze
        inherited_config=freeze(args.inherit_state,frozen)
        if args.unlimited_daily_budget:
            inherited_config.update(daily_point_calls=None,daily_worker_seconds=None,
                                    batch_growth_step=25,batch_growth_every=100,
                                    max_batch_calls=300)
    else:
        projection=read(LOCAL/'r17-60-panel-v1/selection-input.json')
        public_meta=read(ROOT/'elliptic-curves/data/icarm_current.json')
        catalogue_path=args.catalogue
        if catalogue_path:
            catalogue=read(catalogue_path)
        else:
            catalogue_path=ROOT/public_meta['snapshot']
            catalogue=json.loads(gzip.decompress(catalogue_path.read_bytes()))
        require(catalogue['count']==len(catalogue['curves']) and catalogue['count']>=630,'invalid pinned public catalogue')
        inventory_path=ROOT/'elliptic-curves/data/research_curves/database.json'
        inventory=read(inventory_path)
        equations=projection['excluded_models']+[r['ainvs'] for r in catalogue['curves']]+[r['ainvs'] for r in inventory['curves']]
        prior_foundry=[]
        for path in sorted(ART.glob('high-rank-foundry*/job-*.json')):
            r=read(path)
            if r.get('status')=='PASS_CERTIFIED_SEARCH':
                prior_foundry.append((path,r))
                equations.append(r['packet']['curve'])
        reserved=projection['reserved_addresses']+[[r['family'],str(F(r['parameter']))] for r in inventory['curves'] if r.get('family') in FAMILIES and r.get('parameter') is not None]
        snapshot={'pool':projection['pool'],'excluded_models':equations,'reserved_addresses':reserved,
                  'provenance':{'selection_projection_sha256':sha(LOCAL/'r17-60-panel-v1/selection-input.json'),
                                'inventory_sha256':sha(inventory_path),'catalogue_sha256':sha(catalogue_path),
                                'catalogue_file':str(catalogue_path),'frozen_at':time.time()},
                  'boundary':'Exact exclusions relative to pinned repository and public snapshots. No score is a mathematical exclusion.'}
        snapshot['prior_foundry_certificates']={str(p):sha(p) for p,r in prior_foundry}
        snapshot['reserved_addresses'] += [[r['family'],r['parameter']] for p,r in prior_foundry]
        save(data/'intake.json',snapshot,immutable=True);save(data/'catalogue.json',catalogue,immutable=True)
        panel_path=ART/'r17_60_panel_results_v1.json';panel=read(panel_path)
        warm=[r for r in panel['results'] if r['rank_lower_bound']>=25 or r['id'] in panel['late_gain_candidates']]
        for path,r in prior_foundry:
            if r['rank_lower_bound']<=17: continue
            req_path=path.parent/'config.json'
            oldroot=Path(read(req_path)['frozen_root'])
            req=read(oldroot/Path(r['packet_path']).parent/'request.json')
            previous=req['candidate']
            c=apply_result(previous,r)
            imported={'id':r['id'],'family':r['family'],'parameter':r['parameter'],'packet':r['packet'],
                      'rank_lower_bound':r['rank_lower_bound'],'total_calls':c['total_calls'],
                      'complement_calls_since_last_gain':c['stale_calls'],'seed_calls':0,'gain_timeline':c['history']}
            aliases=[w for w in warm if jkey(w['packet']['curve'])==jkey(r['packet']['curve']) and matches(r['packet']['curve'],[w['packet']['curve']])]
            if aliases and all(w['rank_lower_bound']>r['rank_lower_bound'] for w in aliases): continue
            warm=[w for w in warm if w not in aliases]+[imported]
        save(data/'warm.json',{'source':str(panel_path),'sha256':sha(panel_path),'rows':warm},immutable=True)
    config={**DEFAULTS,**inherited_config,'workers':args.workers,'sage':str(Path.home()/'.local/bin/sage'),
            'frozen_root':str(frozen),'export':str(ART/folder.name),
            'intake':'foundry-inputs/intake.json','catalogue':'foundry-inputs/catalogue.json',
            'created_at':time.time(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
            'campaign':'Continuous, resource-budgeted search in renewable daily epochs; no model calls or scheduled AI turns.'}
    require(1<=config['workers']<=8,'worker count outside declared maximum')
    save(folder/'config.json',config,immutable=True)
    probe=subprocess.check_output([config['sage'],'-python','-c',
        'import json,sys,sage.version,numpy; from sage.all import pari; print(json.dumps(dict(python=sys.executable,sage=sage.version.version,pari=str(pari.version()),numpy=numpy.__version__)))'],env=runtime_env(config))
    save(frozen/'foundry-inputs/software.json',json.loads(probe.decode().splitlines()[-1]),immutable=True)
    source_files={str(p.relative_to(frozen)):sha(p) for p in sorted(frozen.rglob('*')) if p.is_file()}
    executables={str(p):sha(p) for p in (Path(config['sage']),Path('/usr/bin/gp'),Path(read(frozen/'foundry-inputs/software.json')['python']))}
    manifest={'schema':'high-rank-foundry-frozen-runtime.v1','frozen_root':str(frozen),
              'files':source_files,'executables':executables,'config_sha256':sha(folder/'config.json')}
    save(folder/'manifest.json',manifest,immutable=True)
    out=Path(config['export']);out.mkdir(parents=True,exist_ok=False)
    for name in ('manifest.json','config.json'):
        shutil.copyfile(folder/name,out/name)
    with zipfile.ZipFile(out/'frozen-runtime.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for name in source_files: archive.write(frozen/name,'research/'+name)
    db=connect(folder)
    with db:
        putmeta(db,'status','PREPARED_PENDING_PREFLIGHT');putmeta(db,'fresh_completed',0)
    db.close()
    # Run preflight from the snapshot, so absent dependencies fail before launch.
    command=[sys.executable,str(frozen/'elliptic-curves/cas/run_high_rank_foundry.py'),'preflight','--folder',str(folder)]
    p=subprocess.run(command,cwd=frozen,env=runtime_env(config),check=False)
    require(p.returncode==0,'frozen runtime preflight failed; preserve this preparation')
    print(json.dumps({'status':'PREPARED','folder':str(folder),'export':str(out),'sources':len(source_files)}))


def runtime_env(config):
    env={k:v for k,v in os.environ.items() if not any(s in k.upper() for s in ('API_KEY','ACCESS_TOKEN','AUTH_TOKEN'))}
    env.update(PATH=str(Path(config['sage']).parent)+os.pathsep+env.get('PATH',''),
               PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
               PYTHONPATH=str(Path(config['frozen_root'])/'elliptic-curves/cas'))
    return env


def preflight(folder):
    frozen=guard(folder);config=read(folder/'config.json');db=connect(folder)
    if (folder/'preflight.json').exists():
        require(read(folder/'preflight.json')['manifest_sha256']==sha(folder/'manifest.json'),'preflight manifest changed')
        db.close();print('FOUNDRY_PREFLIGHT_ALREADY_PASSED');return
    require(meta(db,'status')=='PREPARED_PENDING_PREFLIGHT','preflight cannot reset an active ledger')
    import high_rank_foundry_arithmetic as arithmetic
    import pari_pointed_backend
    pari_pointed_backend.sources()  # All transitive chart bindings must exist.
    intake=read(frozen/config['intake'])
    require(len(intake['pool'])==6144,'retained population count differs')
    warm=read(frozen/'foundry-inputs/warm.json')['rows']
    for r in warm:
        p=r['packet'];arithmetic.verify_packet(p,{**r,'model':p['curve']})
        path=frozen/'foundry-warm'/f"{r['id']}.json";save(path,p,immutable=True)
        c=fresh_record({'id':'warm-'+r['id'],'family':r['family'],'parameter':r['parameter'],
                       'model':p['curve'],'j_key':jkey(p['curve']),'lane':'historical_watchlist',
                       'parameter_height':max(abs(F(r['parameter']).numerator),F(r['parameter']).denominator)})
        c.update(rank=r['rank_lower_bound'],total_calls=r['total_calls'],stale_calls=r['complement_calls_since_last_gain'],
                 state='READY',packet=str(path.relative_to(frozen)),packet_sha256=sha(path),
                 bank_index=0,last_batch_gain=1 if r['complement_calls_since_last_gain']<=35 else 0,
                 gaining_batches=1,history=[{**e,'call':e['call']+(r['seed_calls'] if e['phase'].startswith('complement') else 0)} for e in r['gain_timeline']],
                 historical_baseline=True,revivals=0)
        # These old packets contain exact points, but import does not count as discovery.
        with db: putcurve(db,c)
    inherited_packets=0
    if (frozen/'foundry-inputs/inherited-state.json').exists():
        from high_rank_foundry_migration import restore
        inherited_packets=restore(frozen,db)
    # Exercise a real fresh generic packet and full seed-geometry preparation,
    # without calling the point backend. This catches missing runtime/layout inputs.
    from high_rank_foundry_job import step
    row,_=Intake(intake).next({})
    probe=frozen/'foundry-preflight';probe.mkdir(exist_ok=True)
    save(probe/'request.json',{'candidate':row,'config':config,'bank_index':0},immutable=True)
    step(probe,'generic');step(probe,'seed-prepare')
    with db: putmeta(db,'status','READY')
    save(folder/'preflight.json',{'status':'PASS_FROZEN_NATIVE_INPUTS','warm_packets':len(warm),'inherited_packets':inherited_packets,'point_searches':0,
                                'manifest_sha256':sha(folder/'manifest.json')},immutable=True)
    db.close()
    print('FOUNDRY_PREFLIGHT_PASS',len(warm)+inherited_packets,'certificates; zero point searches',flush=True)


def export_result(folder,j,result):
    config=read(folder/'config.json');frozen=Path(config['frozen_root']);job=frozen/j['path']
    require(result['request_sha256']==sha(job/'request.json'),'job request binding changed')
    for name,digest in result['evidence'].items():
        require(sha(job/name)==digest,'job evidence changed: '+name)
    if result['status']=='PASS_CERTIFIED_SEARCH':
        packet=read(frozen/result['packet_path'])
        require(sha(frozen/result['packet_path'])==result['packet_sha256'],'result packet binding changed')
        receipt=read(job/'packet-verified.json')
        require(receipt['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS' and receipt['packet_sha256']==result['packet_sha256'],'independent packet replay missing')
        require(result['rank_lower_bound']==packet['rank_lower_bound']==receipt['rank_lower_bound']==len(packet['points']),
                'rank header differs from verified points')
        request=read(job/'request.json')['candidate']
        require(result['id']==j['cid'] and all(result[k]==request[k] for k in ('family','parameter')) and
                packet['curve']==request['model'],'certificate candidate binding differs')
        portable={**result,'packet':packet,'certificate_replay':receipt,
                  'runtime_manifest_sha256':sha(folder/'manifest.json')}
    elif result['status']=='UNRESOLVED_GENERIC_SEED':
        request=read(job/'request.json')
        gate=read(job/'generic-seed-gate.json')
        require(request['kind']=='fresh' and not request['candidate'].get('packet'),
                'only an uncertified fresh input can miss the seed gate')
        require(gate==result['seed_gate'] and sha(job/'generic-seed-gate.json')==result['seed_gate_sha256'],
                'generic seed gate binding changed')
        require(gate['status']=='UNRESOLVED_GENERIC_SEED' and 0<=gate['admitted_columns']<=17 and
                (gate['admitted_columns']<17 or gate['no_rational_2_torsion_prime'] is None) and
                gate['point_searches']==result['calls']==0 and result['rank_lower_bound'] is None,
                'invalid unresolved seed outcome')
        require(all(gate[k]==result[k]==request['candidate'][k] for k in ('id','family','parameter')) and
                gate['curve']==request['candidate']['model'],'unresolved seed candidate differs')
        require(not any(job.glob('point-invocations/*/attempt-*.json')),'seed miss after point invocation')
        portable={**result,'runtime_manifest_sha256':sha(folder/'manifest.json')}
    else: portable=result
    path=Path(config['export'])/f"job-{j['id']:07d}-{j['cid']}.json"
    save(path,portable,immutable=True)
    return str(path)


def summary(folder,db):
    config=read(folder/'config.json');cc=curves(db);jj=jobs(db)
    complete=[j for j in jj if j['state']=='DONE']
    stats=defaultdict(lambda:{'attempted':0,'certified':0,'seeded':0,'rank_sum':0,'calls':0,'unresolved':0})
    for c in cc:
        if c['state']=='NEW' or c.get('historical_baseline'): continue
        key=c['family']+'/'+c['lane']
        g=stats[key];g['attempted']+=1;g['unresolved']+=c['state'].startswith('QUARANTINED') or c['state']=='UNCERTIFIED_SEED'
        if c.get('packet'):
            g['certified']+=1;g['seeded']+=c['rank']>17;g['rank_sum']+=c['rank'];g['calls']+=c['total_calls']
    best=sorted((c for c in cc if c.get('packet')),key=lambda c:(-c['rank'],-utility(c),c['id']))[:20]
    result={'status':meta(db,'status'),'updated_at':time.time(),'guardian':meta(db,'guardian'),
            'inherited_from':meta(db,'inherited_from'),
            'controller':meta(db,'controller'),'workers':config['workers'],
            'jobs':dict(Counter(j['state'] for j in jj)),
            'active_jobs':[{k:j.get(k) for k in ('id','cid','kind','pid','token','path')} for j in jj if j['state']=='RUNNING'],
            'new_curves':sum(c.get('packet') is not None and not c.get('historical_baseline') for c in cc),
            'certified_rank_counts':dict(sorted(Counter(c['rank'] for c in cc if c.get('packet') and not c.get('historical_baseline')).items())),
            'point_calls':sum(j.get('logical_calls',j.get('charged_calls',0)) for j in complete),
            'charged_point_budget':sum(j.get('charged_calls',0) for j in jj if j['state'] in ('DONE','FAILED','SKIPPED')),
            'queued_curves':dict(Counter(c['state'] for c in cc)),
            'worker_seconds':sum(j.get('charged_seconds',0) for j in jj if j['state'] in ('DONE','FAILED','SKIPPED')),
            'groups':dict(stats),'best':[{k:c.get(k) for k in ('id','family','parameter','rank','state','total_calls','stale_calls','packet','certificate','historical_baseline')} for c in best],
            'claim_boundary':'Adaptively selected, variably censored lower bounds. Group summaries do not establish a causal law, exact ranks, world novelty, or conductor records.'}
    save(folder/'status.json',result)
    save(Path(config['export'])/'LIVE_STATUS.json',result)
    from high_rank_foundry_report import write_report
    write_report(Path(config['export'])/'REPORT.md',result,config)
    return result


def spawn_job(folder,db,j):
    config=read(folder/'config.json');frozen=Path(config['frozen_root']);job=frozen/j['path']
    job.mkdir(parents=True,exist_ok=True)
    if 'request' in j: save(job/'request.json',j['request'],immutable=True)
    lease=read(job/'lease.json') if (job/'lease.json').exists() else {}
    if same_process(lease.get('pid'),lease.get('token')):
        j.update(state='RUNNING',pid=lease['pid'],token=lease['token'],launched_at=lease['started_at'])
        with db: putjob(db,j)
        return
    with (job/'job.log').open('ab',buffering=0) as log:
        proc=subprocess.Popen([sys.executable,str(frozen/'elliptic-curves/cas/high_rank_foundry_job.py'),str(job)],
             cwd=frozen,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
             start_new_session=True,close_fds=True,env=runtime_env(config))
    token=process_info(proc.pid)
    j.update(state='RUNNING',pid=proc.pid,token=token['start_token'] if token else None,launched_at=time.time())
    with db: putjob(db,j)


def consume(folder,db,j):
    config=read(folder/'config.json');frozen=Path(config['frozen_root']);path=frozen/j['path']/'result.json'
    guard(folder)
    result=read(path);j['result_sha256']=sha(path)
    j['export']=export_result(folder,j,result)
    j['charged_seconds']=max(time.time()-j['created_at'],result['wall_seconds'])
    j['charged_calls']=result.get('point_invocation_intents',result['calls']) if result['calls'] is not None else j['reserved_calls']
    j['logical_calls']=result['calls']
    j['finished_at']=time.time()
    c=json.loads(db.execute('SELECT record FROM curves WHERE id=?',(j['cid'],)).fetchone()[0])
    failures=meta(db,'consecutive_failures',0)
    if result['status']=='PASS_CERTIFIED_SEARCH':
        c=apply_result(c,result)
        c['certificate']=j['export']
        if result['unresolved_cloud']: c['state']='QUARANTINED_CLOUD'
        c['novelty']=result['novelty'];c['conductor_gate']=result['conductor_gate'];j['state']='DONE'
        j['rank_before']=j['initial_rank'];j['rank_after']=c['rank'];failures=0
        if j['kind']=='fresh': putmeta(db,'fresh_completed',meta(db,'fresh_completed',0)+1)
        if c['rank']>=28 and c['rank']>j['initial_rank']:
            reported=c['novelty'].get('highest_reported_lower_bound')
            label='HIT' if reported is None or c['rank']>reported else 'REPRODUCTION'
            save(Path(config['export'])/f"{label}-rank-{c['rank']}-{j['id']:07d}.json",
                 {'id':c['id'],'rank_lower_bound':c['rank'],'certificate':j['export'],'novelty':c['novelty'],
                  'claim_boundary':'Certified lower bound. Public novelty is only relative to the pinned catalogue.'},immutable=True)
    elif result['status']=='UNRESOLVED_GENERIC_SEED':
        c['state']='UNCERTIFIED_SEED';c['rank']=None;c['seed_gate_export']=j['export'];j['state']='SKIPPED'
        # Neutral for the engineering failure streak: neither increments nor
        # resets it. Several malformed jobs cannot hide behind input misses.
    elif result['status']=='PASS_CONDUCTOR_REPLAY':
        c['conductor']=result['conductor'];c['conductor_comparison']=result['comparison'];j['state']='DONE';failures=0
        if result['comparison']=='BEATS_PINNED_REPORTED_THRESHOLD':
            save(Path(config['export'])/f"CONDUCTOR-HIT-{j['id']:07d}.json",{'certificate':j['export'],'comparison':result['comparison']},immutable=True)
    else:
        c['state']='QUARANTINED';c['failure_export']=j['export'];j['state']='FAILED';failures+=1
    with db:
        putcurve(db,c);putjob(db,j);putmeta(db,'consecutive_failures',failures)
        if failures>=config['max_consecutive_failures']: putmeta(db,'status','HALTED_REPEATED_FAILURE')
    print(json.dumps({'job':j['id'],'curve':c['id'],'state':j['state'],'rank':c['rank'],'calls':j['charged_calls']}),flush=True)


def dispatch(folder,db,intake):
    config=read(folder/'config.json');frozen=Path(config['frozen_root']);cc=curves(db);jj=jobs(db)
    busy={j['cid'] for j in jj if j['state'] in ('RUNNING','PENDING')}
    available=[c for c in cc if c['state']=='READY' and c['id'] not in busy]
    exposure={'fresh':0,'exploit':0,'conductor':0}
    for j in jj:
        exposure[j['kind']]+=j.get('charged_seconds',min(config['job_seconds'],max(1,time.time()-j['created_at'])))
    fresh_count=meta(db,'fresh_completed',0)
    revival=False;rejected=[]
    # Revival gets a bounded minority of exploit slots, never steals fresh share.
    cooled=[c for c in cc if c['state']=='COOLED' and c['id'] not in busy and c.get('revivals',0)<config['max_revivals']]
    revive_mark=meta(db,'last_revival',0)
    revival_due=bool(cooled) and fresh_count-revive_mark>=config['revival_every']
    kind=choose_lane(exposure['fresh'],exposure['exploit'],config,has_exploit=bool(available) or revival_due)
    if kind=='exploit' and revival_due:
        c=max(cooled,key=lambda c:(utility(c),c['id']));revival=True
        c['revivals']=c.get('revivals',0)+1;c['bank_index']=max(3,c['bank_index']+1)
        putmeta(db,'last_revival',fresh_count)
    elif kind=='fresh':
        seen=defaultdict(list)
        for c in cc: seen[c['j_key']].append(c['model'])
        row,rejected=intake.next(seen);c=fresh_record(row)
    else:
        c=max(available,key=lambda c:(utility(c),c['id']))
    # Exact factorization gets at most 2% of charged/reserved worker time.
    eligible=[c for c in cc if c.get('conductor_gate',{}).get('eligible') and c.get('packet') and
              c['id'] not in busy and c.get('conductor_attempted_rank',0)<c['rank']]
    total=sum(exposure.values())
    if kind=='exploit' and not revival and eligible and total>600 and exposure['conductor']+150<=config['conductor_share']*total:
        kind='conductor';c=max(eligible,key=lambda c:(c['rank'],c['id']));c['conductor_attempted_rank']=c['rank']
    bank=c['bank_index'] if kind=='fresh' or revival else next_bank(c)
    use_cached=bool(kind=='exploit' and not revival and c.get('head') and bank==c['bank_index'])
    allowance=0 if kind=='conductor' else effective_batch_calls(config,fresh_count)
    jid=(db.execute('SELECT max(id) FROM jobs').fetchone()[0] or 0)+1
    path=frozen/'foundry-jobs'/f'job-{jid:07d}'
    request={'schema':'foundry-job.v1','kind':kind,'candidate':c,'bank_index':bank,'use_cached':use_cached,
             'revival':revival,'allowance':allowance,'config':config,'catalogue':config['catalogue'],
             'conductor_gate':c.get('conductor_gate'),'decision':{'utility':utility(c),'worker_exposure_seconds':exposure,
             'fresh_share':config['fresh_share'],'rejected':rejected,'intake_cursor_after':intake.cursor}}
    j={'id':jid,'cid':c['id'],'state':'PENDING','kind':kind,'path':str(path.relative_to(frozen)),
       'created_at':time.time(),'day':int(time.time()//86400),'initial_rank':c['rank'],'retries':0,
       'reserved_calls':(allowance+98+config['max_crash_retries']) if kind=='fresh' or (revival and c['rank']==17) else 0 if kind=='conductor' else allowance+config['max_crash_retries'],
       'reserved_seconds':150 if kind=='conductor' else config['job_seconds'],'request':request}
    with db:
        putcurve(db,c);putjob(db,j);putmeta(db,'intake_cursor',intake.cursor)
    spawn_job(folder,db,j)


def controller(folder,job_limit=None):
    fd=exclusive(folder/'controller.lock');frozen=guard(folder);config=read(folder/'config.json');db=connect(folder)
    token=process_info(os.getpid())
    with db:
        putmeta(db,'controller',{'pid':os.getpid(),'token':token['start_token']});putmeta(db,'status','RUNNING')
    intake=Intake(read(frozen/config['intake']),meta(db,'intake_cursor'))
    children={};stop=False
    signal.signal(signal.SIGTERM,lambda *_:(folder/'STOP').touch())
    try:
        while True:
            for j in jobs(db,('RUNNING','PENDING')):
                job=frozen/j['path']
                lease=read(job/'lease.json') if (job/'lease.json').exists() else {}
                if same_process(lease.get('pid'),lease.get('token')) and j.get('pid')!=lease['pid']:
                    j.update(state='RUNNING',pid=lease['pid'],token=lease['token'],launched_at=lease['started_at'])
                    with db: putjob(db,j)
                if (job/'result.json').exists():
                    consume(folder,db,j)
                    try: os.waitpid(j.get('pid',-1),os.WNOHANG)
                    except (ChildProcessError,ProcessLookupError): pass
                elif same_process(j.get('pid'),j.get('token')):
                    if time.time()-j['launched_at']>config['job_seconds']:
                        os.kill(j['pid'],signal.SIGKILL)
                elif j['state']=='PENDING' or j['retries']<config['max_crash_retries']:
                    if j['state']!='PENDING': j['retries']+=1
                    guard(folder);spawn_job(folder,db,j)
                else:
                    save(job/'result.json',{'status':'UNRESOLVED_JOB_FAILURE','error':'Repeated process interruption; raw checkpoints retained.',
                        'calls':None,'wall_seconds':time.time()-j['created_at'],'request_sha256':sha(job/'request.json'),
                        'evidence':{}},immutable=True)
            active=jobs(db,('RUNNING','PENDING'))
            stop=(folder/'STOP').exists() or meta(db,'status')=='HALTED_REPEATED_FAILURE'
            count=len(jobs(db))
            if job_limit is not None and count>=job_limit: stop=True
            if stop and not active:
                with db:
                    if meta(db,'status')!='HALTED_REPEATED_FAILURE': putmeta(db,'status','STOPPED')
                summary(folder,db);return
            if not stop and len(active)<config['workers']:
                day=int(time.time()//86400);today=[j for j in jobs(db) if j['state'] in ('RUNNING','PENDING') or int(j.get('finished_at',j['created_at'])//86400)==day]
                calls=sum(j.get('charged_calls',j['reserved_calls']) for j in today)
                wall=sum(j.get('charged_seconds',j['reserved_seconds']) for j in today)
                disk=shutil.disk_usage(folder).free;inodes=os.statvfs(folder).f_favail
                wait_reason=None
                point_cap=config.get('daily_point_calls')
                worker_cap=config.get('daily_worker_seconds')
                if ((point_cap is not None and calls+198+config['max_crash_retries']>point_cap) or
                    (worker_cap is not None and wall+config['job_seconds']>worker_cap)):
                    wait_reason='WAIT_DAILY_BUDGET'
                elif disk<config['min_free_gib']*1024**3 or inodes<config['min_free_inodes']:
                    wait_reason='WAIT_DISK_RESERVE'
                with db: putmeta(db,'status',wait_reason or 'RUNNING')
                if wait_reason is None:
                    guard(folder);dispatch(folder,db,intake)
            summary(folder,db)
            time.sleep(3)
    except BaseException:
        with db: putmeta(db,'last_controller_error',traceback.format_exc())
        raise
    finally:
        db.close();os.close(fd)


def guardian(folder,job_limit=None):
    fd=exclusive(folder/'guardian.lock');config=read(folder/'config.json');frozen=guard(folder)
    db=connect(folder)
    with db: putmeta(db,'guardian',{'pid':os.getpid(),'token':process_info(os.getpid())['start_token']})
    db.close();crashes=0
    while True:
        guard(folder)
        with (folder/'controller.log').open('ab',buffering=0) as log:
            proc=subprocess.Popen([sys.executable,str(frozen/'elliptic-curves/cas/run_high_rank_foundry.py'),
                 'controller','--folder',str(folder)]+(['--job-limit',str(job_limit)] if job_limit is not None else []),cwd=frozen,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                 close_fds=True,env=runtime_env(config))
            started=time.monotonic();rc=proc.wait()
        if time.monotonic()-started>600: crashes=0
        if rc==0:
            db=connect(folder);state=meta(db,'status');db.close()
            if (folder/'STOP').exists() or state in ('STOPPED','HALTED_REPEATED_FAILURE','HALTED_CONTROLLER_FAILURE'):
                break
            # A bounded controller run ending cleanly is an opportunity to
            # continue the autonomous campaign, not a reason to leave it idle.
            time.sleep(1)
            continue
        crashes+=1
        save(folder/'guardian-error.json',{'controller_returncode':rc,'crashes':crashes,'at':time.time()})
        if crashes>config['max_crash_retries']:
            db=connect(folder)
            with db: putmeta(db,'status','HALTED_CONTROLLER_FAILURE')
            summary(folder,db);db.close();break
        time.sleep(min(30,2**crashes))
    os.close(fd)


def launch(folder,job_limit=None):
    config=read(folder/'config.json');frozen=guard(folder)
    require(read(folder/'preflight.json')['status']=='PASS_FROZEN_NATIVE_INPUTS','preflight required')
    fd=exclusive(folder/'launch.lock')
    try:
        db=connect(folder)
        for role in ('guardian','controller'):
            p=meta(db,role,{}) or {}
            require(not same_process(p.get('pid'),p.get('token')),'foundry already running')
        db.close();(folder/'STOP').unlink(missing_ok=True)
        with (folder/'guardian.log').open('ab',buffering=0) as log:
            proc=subprocess.Popen([sys.executable,str(frozen/'elliptic-curves/cas/run_high_rank_foundry.py'),
                 'guardian','--folder',str(folder)]+(['--job-limit',str(job_limit)] if job_limit is not None else []),cwd=frozen,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                 start_new_session=True,close_fds=True,env=runtime_env(config))
        info=process_info(proc.pid)
        require(info is not None,'guardian failed before process identity could be recorded')
        db=connect(folder)
        with db:putmeta(db,'guardian',{'pid':proc.pid,'token':info['start_token']})
        db.close()
        save(folder/'launch-receipt.json',{'status':'LAUNCHED','guardian_pid':proc.pid,'token':info['start_token'],
             'at':time.time(),'job_limit':job_limit,'manifest_sha256':sha(folder/'manifest.json')})
        print(json.dumps({'status':'LAUNCHED','guardian_pid':proc.pid,'folder':str(folder)}))
    finally: os.close(fd)


def status(folder):
    db=connect(folder);result=summary(folder,db);db.close()
    for role in ('guardian','controller'):
        p=result.get(role) or {};result[role+'_alive']=same_process(p.get('pid'),p.get('token'))
    if result['status'] in ('RUNNING','WAIT_DAILY_BUDGET','WAIT_DISK_RESERVE') and not result['controller_alive']:
        result['observed_status']='INTERRUPTED_OR_RESTARTING'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','preflight','launch','resume','guardian','controller','status','stop','verify'))
    p.add_argument('--folder',type=Path,default=DEFAULT_FOLDER);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--catalogue',type=Path);p.add_argument('--job-limit',type=int)
    p.add_argument('--unlimited-daily-budget',action='store_true')
    p.add_argument('--inherit-state',type=Path,help='Preserve a stopped campaign ledger and intake in a new source snapshot')
    a=p.parse_args();folder=a.folder.resolve()
    if a.mode=='prepare':prepare(folder,a)
    elif a.mode=='preflight':preflight(folder)
    elif a.mode in ('launch','resume'):launch(folder,a.job_limit)
    elif a.mode=='guardian':guardian(folder,a.job_limit)
    elif a.mode=='controller':controller(folder,a.job_limit)
    elif a.mode=='stop':(folder/'STOP').touch();print('Graceful stop requested; current jobs finish replay.')
    elif a.mode=='verify':guard(folder);print('FROZEN_MANIFEST_PASS')
    else:print(json.dumps(status(folder),indent=2))


if __name__=='__main__':main()

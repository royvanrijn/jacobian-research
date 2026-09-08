#!/usr/bin/env python3
"""Retrospective two-seed V3 amplifier experiment for curve 302.

The only oracle use is seed construction: recover the exact displayed directions
`recovered-strict-02` and `recovered-strict-03` from the immutable M31 diagnostic.
Each seed is independently certified with the generic M17 prefix to give rank 18.
After `seed-proof.json` is frozen, the search context contains only that rank-18
basis, the generic degree-two orbit table, and the unchanged V3 numerical policy.

Commands (standard Python): launch | resume | status | diagnose | stop
The detached worker is launched under Sage and runs both seeds sequentially.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from fractions import Fraction as F

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
D = LOCAL/'curve302-seeded-v3-amplifier-v1'
STATE = D/'state.json'
LOG = D/'worker.log'
SELF = Path(__file__).resolve()
V3 = LOCAL/'adaptive-visibility-cascade-v3'
ORBITS = ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
VISIBILITY = ART/'curve302_residual_visibility_geometry_v1.json'
M24 = ART/'curve302_recovered_followup_wave_03_mod2_v1.json'
SEEDS = ('recovered-strict-02','recovered-strict-03')
INITIAL_RANK = 18
TARGET_RANK = 31


def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name('.'+path.name+'.tmp-'+str(os.getpid()))
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)

def process_alive(pid):
    if not isinstance(pid,int) or pid<=0:return False
    try: os.kill(pid,0)
    except (ProcessLookupError,PermissionError): return False
    try:
        stat=Path(f'/proc/{pid}/stat').read_text().split()
        return len(stat)>2 and stat[2] != 'Z'
    except Exception: return True

def sage_launcher():
    for value in (os.environ.get('V3_SAGE'),shutil.which('sage'),str(Path.home()/'.local/bin/sage'),'/usr/bin/sage'):
        if value and Path(value).is_file() and os.access(value,os.X_OK): return str(Path(value).resolve())
    raise RuntimeError('Sage launcher not found; set V3_SAGE=/absolute/path/to/sage')

def tail(path,n=10000):
    path=Path(path)
    if not path.exists(): return ''
    data=path.read_bytes(); return data[-n:].decode(errors='replace')

def update(status,**extra):
    old=read(STATE) if STATE.exists() else {}
    old.update(extra,status=status,updated_unix=time.time())
    if status.startswith('RUNNING') or status in ('PREPARING','REPLAYING'):
        old['pid']=os.getpid()
    atomic(STATE,old)
    print('CURVE302_SEEDED_V3',status,json.dumps(extra,sort_keys=True),flush=True)


def construct_seed(seed_id):
    """Oracle stage only: reconstruct one displayed M31 direction, then seal M18."""
    from importlib.machinery import SourceFileLoader
    from v3_warm_support import atomic as immutable_atomic, point_tuple, curve_tuple, require
    from v3_warm_engine import certified_state
    from memory_rank_certificate import checked_rank

    folder=D/seed_id; folder.mkdir(parents=True,exist_ok=True)
    seed_path,proof_path,protocol_path=folder/'seed-input.json',folder/'seed-proof.json',folder/'protocol.json'
    if seed_path.exists() and proof_path.exists() and protocol_path.exists():
        return

    diag=SourceFileLoader('curve302_seed_diag',str(CAS/'audit_curve302_exceptional_subgroup_landscape.sage')).load_module()
    group=SourceFileLoader('curve302_seed_group',str(CAS/'half_lattice_pointed_sieve.py')).load_module()
    visibility=read(VISIBILITY); m24=read(M24)
    require(visibility['status']=='PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC','visibility oracle not passed')
    model=tuple(F(v) for v in m24['curve'])
    base17=tuple(tuple(F(v) for v in p) for p in m24['independent_points'][:17])
    entry=next((row for row in visibility['directions'] if row['id']==seed_id),None)
    require(entry is not None,'requested seed direction absent')
    target=diag.primary_target(entry,base17,model,group)
    points=(*base17,target)

    parent_proof=m24['rank_certificate']
    primes=[int(s['prime']) for s in parent_proof['signatures']]
    torsion=int(parent_proof['no_rational_2_torsion_prime'])
    proof=checked_rank(model,points,primes,torsion)
    require(proof['rank_lower_bound']==18,'seed did not certify rank 18')
    state=certified_state(model,points,proof)
    require(state.rank==18 and tuple(state.basis)==points,'certified M18 seed changed order')

    engine=SourceFileLoader('curve302_seed_policy_engine',str(CAS/'adaptive_visibility_cascade_v3.sage')).load_module()
    original=read(V3/'protocol.json')
    require(engine.sources()==original['sources'],'frozen V3 numerical sources changed')
    policy=dict(original)
    policy.update(schema='curve302-seeded-v3-amplifier.v1',initial_rank=18,target_rank=31,
                  oracle_seed_direction=seed_id,
                  scope=('Retrospective amplifier test. Known M31 data are used only to construct and certify '
                         'the displayed M18 seed. Search thereafter reuses the unchanged frozen V3 numerical policy.'))
    source_bindings={
        str(VISIBILITY.relative_to(ROOT)):sha(VISIBILITY),str(M24.relative_to(ROOT)):sha(M24),
        str(ORBITS.relative_to(ROOT)):sha(ORBITS),str((V3/'protocol.json').relative_to(ROOT)):sha(V3/'protocol.json'),
        str(SELF.relative_to(ROOT)):sha(SELF),
    }
    seed={'family':'curve302-det1092','seed_direction':seed_id,'curve':list(map(str,model)),
          'points':[list(map(str,p)) for p in points],'initial_rank':18,'generic_rank':17,
          'oracle_boundary':'The 18th point is recovered from immutable known-M31 diagnostic data; no later search step may read that diagnostic.'}
    immutable_atomic(seed_path,seed,immutable=True)
    immutable_atomic(proof_path,proof,immutable=True)
    policy['seed_inputs']={str(seed_path.relative_to(ROOT)):sha(seed_path),str(proof_path.relative_to(ROOT)):sha(proof_path)}
    policy['oracle_inputs']=source_bindings
    immutable_atomic(protocol_path,policy,immutable=True)
    print('SEEDED_PREPARED',seed_id,'rank18',flush=True)


def context(seed_id):
    from importlib.machinery import SourceFileLoader
    from v3_warm_engine import Context, certified_state
    from v3_warm_support import point_tuple, curve_tuple, require, sha as vsha
    folder=D/seed_id; policy=read(folder/'protocol.json')
    seed,proof=read(folder/'seed-input.json'),read(folder/'seed-proof.json')
    model=curve_tuple(seed['curve']); points=point_tuple(seed['points'])
    require(seed['seed_direction']==seed_id and len(points)==18,'wrong sealed seed')
    require(vsha(folder/'seed-input.json')==policy['seed_inputs'][str((folder/'seed-input.json').relative_to(ROOT))],'seed input changed')
    require(vsha(folder/'seed-proof.json')==policy['seed_inputs'][str((folder/'seed-proof.json').relative_to(ROOT))],'seed proof changed')
    state=certified_state(model,points,proof); require(state.rank==18,'M18 seed replay failed')
    engine=SourceFileLoader('curve302_seeded_numeric_'+seed_id.replace('-','_'),str(CAS/'adaptive_visibility_cascade_v3.sage')).load_module()
    require(engine.sources()==read(V3/'protocol.json')['sources'],'V3 numerical engine changed')
    engine.D=folder; engine.ORBITS=ORBITS; engine.v1.D=folder; engine.v1.ORBITS=ORBITS
    # Deliberately do not call v1.guard(): its D-specific generic-input contract belongs
    # to the original M17 calibration. All numerical source bytes are hash-checked above.
    sources={**engine.sources(),str(SELF.relative_to(ROOT)):sha(SELF)}
    return Context(seed_id,folder,ROOT,policy,engine,model,state,sources)


def verify_case(ctx):
    """Independent replay adapted from det1092_v3_replay, with initial rank 18."""
    from v3_warm_support import (atomic as immutable_atomic, bindings, check_chart, curve_tuple,
        indexed_paths, point_tuple, require, sha as vsha, terminal_structure, within)
    from v3_warm_engine import certified_state, restore_state
    from v3_warm_replay import verify_landscape
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl

    folder,engine,model,state,policy=ctx.folder,ctx.engine,ctx.model,ctx.state,ctx.policy
    terminal=terminal_structure(folder); require(terminal is not None,'no sealed terminal')
    require(terminal['initial_rank']==18,'terminal did not start from M18')
    require(terminal.get('execution_sources')==ctx.sources,'execution source binding differs')
    vd=folder/'seeded-verification'; vd.mkdir(exist_ok=True)
    mapper=engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000,silent=True)
    tested=set(); reports=[]
    for stage in terminal['stages']:
        wd=folder/f"replay-M17/epoch-{stage['epoch']:02d}"
        basis=point_tuple(state.basis); selection=read(wd/'selection.json')
        require(point_tuple(selection['basis'])==basis,'replay subgroup differs')
        metrics=verify_landscape(engine,model,basis,tested,wd,policy)
        charts=indexed_paths(wd,expected=stage['charts']); require(charts,'stage without charts')
        state=restore_state(read(wd/'epoch-seed-state.json'),model,basis)
        prefix=[]; expected=list(basis); seen={(x,abs(y)) for x,y in basis}
        for j,path in enumerate(charts):
            chart=read(path); check_chart(chart,selection,j); prefix.append(chart)
            require(mapper.mapping(model,basis,chart['centre'])==chart['mapping'],'quartic map differs')
            search=PointedQuarticSearch(state=state,centre={'coefficients':chart['centre']['representative']},
                                        coordinate_policy=chart['mapping']['coordinate_policy'])
            backend.replay(search,chart['mapping'],chart['search'])
            for raw in chart['search']['finite_curve_points']:
                p=F(raw['x']),F(raw['y']); key=p[0],abs(p[1])
                if key not in seen: seen.add(key); expected.append(p)
            cloud_path,audit_path=wd/f'cloud-{j:03d}.json',wd/f'mod2-{j:03d}.json'
            snap,audit=read(cloud_path),read(audit_path)
            require(snap['charts']==prefix and snap['final_state']==state.record(),'cumulative snapshot differs')
            require(audit['input_sha256']==vsha(cloud_path) and point_tuple(audit['points'])==tuple(expected),'mod2 audit cloud differs')
            mod2.check(audit_path); tested.add(tuple(chart['centre']['point']))
            if j < len(charts)-1: require(audit['rank_lower_bound']==len(basis),'stale chart after gain')
        last=read(wd/stage['audit']); enlarged=point_tuple(last['independent_points'])
        require(enlarged[:len(basis)]==basis and len(enlarged)==stage['after'],'stage rank/prefix differs')
        state=certified_state(model,enlarged,last['rank_certificate'])
        modl.check(wd/'modl.json'); odd=read(wd/'modl.json')
        ranks={str(a['modulus']):a['finite_column_rank'] for a in odd['audits']}
        require(ranks=={'3':state.rank,'5':state.rank},'odd-prime ranks differ')
        reports.append({'epoch':stage['epoch'],'before':len(basis),'after':state.rank,
                        'charts_replayed':len(charts),'landscape':metrics,'modl_ranks':ranks})
        print('SEEDED_REPLAY',ctx.case,'epoch',stage['epoch'],'rank',state.rank,flush=True)
    require(state.rank==terminal['final_rank_lower_bound'],'terminal/replay rank differs')
    result={'schema':'curve302-seeded-v3-amplifier-result.v1','status':'PASS_INDEPENDENT_SEEDED_V3_REPLAY',
            'seed_direction':ctx.case,'initial_rank':18,'rank_lower_bound':state.rank,'gain':state.rank-18,
            'charts':terminal['charts'],'stop_reason':terminal['stop_reason'],'stages':reports,
            'terminal_sha256':vsha(folder/'replay-M17/terminal.json'),'protocol_sha256':vsha(folder/'protocol.json'),
            'claim_boundary':'Retrospective known-seed amplifier experiment; not a prospective selector or new rank claim.'}
    path=folder/'seeded-verified.json'
    if path.exists(): require(read(path)==result,'existing seeded verification differs')
    else: immutable_atomic(path,result,immutable=True)
    print('SEEDED_VERIFIED',ctx.case,'18->'+str(state.rank),terminal['stop_reason'],flush=True)
    return result


def worker():
    from v3_warm_support import require
    import det1092_v3_worker as searcher
    update('PREPARING')
    results=[]
    for seed_id in SEEDS:
        update('PREPARING',case=seed_id); construct_seed(seed_id)
        ctx=context(seed_id)
        verified=ctx.folder/'seeded-verified.json'
        if verified.exists():
            result=read(verified); results.append(result); continue
        terminal=ctx.folder/'replay-M17/terminal.json'
        if not terminal.exists():
            update('RUNNING_SEARCH',case=seed_id,initial_rank=18)
            searcher.run_search(ctx)
        update('REPLAYING',case=seed_id)
        result=verify_case(ctx); results.append(result)
    summary={'status':'COMPLETE_TWO_SEED_AMPLIFIER','results':results,
             'best_rank_lower_bound':max(r['rank_lower_bound'] for r in results),
             'any_full_M31':any(r['rank_lower_bound']>=31 for r in results)}
    atomic(D/'summary.json',summary); update('COMPLETE_TWO_SEED_AMPLIFIER',results=results)


def launch(resume=False):
    D.mkdir(parents=True,exist_ok=True)
    old=read(STATE) if STATE.exists() else {}
    if process_alive(old.get('pid')): raise RuntimeError('seeded V3 worker is already alive')
    if old and not resume and old.get('status') not in (None,'NOT_LAUNCHED'):
        raise RuntimeError('existing experiment; use resume')
    sage=sage_launcher()
    with LOG.open('ab',buffering=0) as log:
        proc=subprocess.Popen([sage,'-python','-u',str(SELF),'worker'],cwd=ROOT,stdin=subprocess.DEVNULL,
                              stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                              env={**os.environ,'PYTHONUNBUFFERED':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    atomic(STATE,{'status':'LAUNCHED','pid':proc.pid,'sage':sage,'updated_unix':time.time()})
    print('Launched seeded curve302 V3 amplifier pid='+str(proc.pid))
    print('Status: python3 elliptic-curves/cas/run_curve302_seeded_v3_amplifier.py status')

def status():
    row=read(STATE) if STATE.exists() else {'status':'NOT_LAUNCHED'}
    row['process_alive']=process_alive(row.get('pid'))
    row['progress_tail']=tail(LOG,12000).splitlines()[-8:]
    cases={}
    for seed in SEEDS:
        folder=D/seed
        if (folder/'seeded-verified.json').exists(): cases[seed]=read(folder/'seeded-verified.json')
        elif (folder/'replay-M17/terminal.json').exists(): cases[seed]={'status':'SEARCH_SEALED_REPLAY_PENDING',**read(folder/'replay-M17/terminal.json')}
        elif folder.exists(): cases[seed]={'status':'PREPARED_OR_SEARCHING'}
        else: cases[seed]={'status':'PENDING'}
    row['cases']={k:{x:v.get(x) for x in ('status','initial_rank','rank_lower_bound','gain','charts','stop_reason')} for k,v in cases.items()}
    print(json.dumps(row,indent=2,sort_keys=True))

def diagnose(): status(); print('\n--- worker log tail ---\n'+tail(LOG,30000))
def stop():
    row=read(STATE) if STATE.exists() else {}
    pid=row.get('pid')
    if not process_alive(pid): print('No live worker'); return
    os.kill(pid,signal.SIGTERM); print('Stop requested; checkpoints retained')

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('action',choices=('launch','resume','status','diagnose','stop','worker')); a=p.parse_args()
    if a.action=='worker': worker()
    elif a.action=='launch': launch(False)
    elif a.action=='resume': launch(True)
    elif a.action=='status': status()
    elif a.action=='diagnose': diagnose()
    else: stop()

if __name__=='__main__': main()

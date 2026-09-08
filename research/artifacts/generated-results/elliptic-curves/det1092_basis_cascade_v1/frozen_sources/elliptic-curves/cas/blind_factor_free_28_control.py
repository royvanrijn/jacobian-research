#!/usr/bin/env python3
"""Full-roster factor-free recovery regression; only original27 seed data.

The curve is a known public28 control. No public point, old winning chart,
old map, or outcome certificate enters geometry or point execution.
"""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
D=ROOT/'artifacts/local/elliptic-curves/blind-factor-free-28-control-v1'
ORIGINAL=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1/seed.json'
SEED=D/'seed.json'
READS=set()

def guard(event,args):
    if event!='open' or not isinstance(args[0],(str,bytes)):return
    path=Path(args[0]).resolve()
    if not path.is_relative_to(ROOT/'artifacts'):return
    if not (path.is_relative_to(D) or path==ORIGINAL):
        raise PermissionError('control data isolation rejected '+str(path.relative_to(ROOT)))
    READS.add(str(path.relative_to(ROOT)))

sys.addaudithook(guard)
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROW={'initial_rank':27,'generic_dimension':17}
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in
            [Path(__file__),CAS/'prepare_blind_factor_free_28.sage']}}

def protocol():
    p=cert.read(D/'protocol.json')
    assert p['sources']==sources()
    assert cert.hashed(SEED)==p['seed_sha256'] and cert.hashed(ORIGINAL)==p['original_seed_sha256']
    return p

def masks(p):
    result=[];i=0
    while len(result)<2048:
        m=int(digest([p['sample_domain'],i]),16)%(1<<27);i+=1
        if m>>17 and m not in result:result.append(m)
    return result

def configure():
    engine.D,engine.SEED,engine.ROW=D,SEED,ROW
    engine.protocol,engine.masks=protocol,masks

def freeze():
    assert not (D/'protocol.json').exists()
    seed=cert.read(ORIGINAL)
    assert len(seed['points'])==27 and len(seed['generic_points'])==17
    assert seed['points'][:17]==seed['generic_points']
    proof=seed['rank_certificate']
    actual=checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],
                        [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    assert digest(actual)==digest(proof)
    checkpoint(SEED,seed)
    checkpoint(D/'protocol.json',dict(schema='blind-factor-free-28-control.v1',sources=sources(),
        seed_sha256=cert.hashed(SEED),original_seed_sha256=cert.hashed(ORIGINAL),
        sample_domain='full11952-specialized-followup-v1',sample_size=2048,charts=49,
        height=125000,seconds_per_chart=10,target_rank=28,rank_stop=False,
        gp_sha256=cert.hashed(Path('/usr/bin/gp')),rss_bytes=2147483648,
        geometry_wall_seconds=180,worker_wall_seconds=600,replay_wall_seconds=300,
        maximum_workers=1,maximum_point_boxes=49,
        selection='Regenerate all2048 masks, rounded384-bit metric, LLL/CVP and49 largest computed norms from original27 seed. No chart index or public exceptional point input.',
        success='All49 declared boxes complete; returned point cloud independently certifies rank at least28 from the original27 seed.',
        boundary='Known-curve regression with public-point data blinded. Prior knowledge of a one-chart success is disclosed; it does not enter this fixed full-roster selection. No new discovery claim.',
        isolation='Python audit hook rejects all repository artifact reads outside this run directory and the original27 seed, starting before research module imports.'))
    checkpoint(D/'freeze-data-access.json',sorted(READS));print('FROZEN49 control boxes; original27-only input')

def launch():
    p=protocol();out=D/'ledger.json';assert not out.exists()
    ledger={'status':'RUNNING','stages':[]};checkpoint(out,ledger)
    jobs=[('geometry',[SAGE,str(CAS/'prepare_blind_factor_free_28.sage')],180),
          ('worker',[sys.executable,str(Path(__file__)),'worker'],600),
          ('replay',[sys.executable,str(Path(__file__)),'replay'],300)]
    for name,cmd,limit in jobs:
        s=run(cmd,limits=Limits(limit,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
        print(name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
        if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed control stage')
    ledger['status']='PASS';checkpoint(out,ledger)
    result=cert.read(D/'result.json')
    print('CONTROL completed',len(result['charts']),'rank',result['rank_lower_bound'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','launch','worker','replay']);args=parser.parse_args()
    configure()
    if args.mode in ('worker','replay'):
        getattr(engine,args.mode)();checkpoint(D/(args.mode+'-data-access.json'),sorted(READS))
    else:globals()[args.mode]()

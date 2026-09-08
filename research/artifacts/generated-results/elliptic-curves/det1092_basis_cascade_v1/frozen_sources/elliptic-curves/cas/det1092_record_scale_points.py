#!/usr/bin/env python3
"""Prospective deployment of calibrated adaptive point exposure on48 fixed fibres."""
import argparse,sys,shutil
from pathlib import Path
import det1092_record_scale_selection as selection
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT,CAS,ART,SAGE=selection.ROOT,selection.CAS,selection.ART,selection.SAGE
cert=selection.cert;BATCH=ROOT/'artifacts/local/elliptic-curves/det1092-record-scale-points-v1'

def sources():
    names=['det1092_record_scale_points.py','prepare_det1092_record_scale_seed.sage',
        'prepare_det1092_record_scale_maps.sage','certify_factor_free_exposure_v3.py','verify_factor_free_rank.sage']
    return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def freeze():
    assert not (BATCH/'protocol.json').exists();q=selection.protocol()
    result=cert.read(selection.D/'selection-result.json');assert result['status']=='PASS'
    assert cert.read(selection.D/'ledger.json')['status']=='PASS'
    assert len(result['selected'])==48 and len({r['id'] for r in result['selected']})==48
    paths=[selection.D/'protocol.json',selection.D/'selection-result.json',selection.SOURCE,selection.GATE]
    checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.det1092-record-scale-points.v1',sources=sources(),
        inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},rows=result['selected'],maximum_curves=48,
        maximum_waves_per_curve=5,maximum_point_boxes=11760,charts=49,sample_size=2048,
        sample_domain='full11952-specialized-followup-v1',height=125000,seconds_per_chart=10,
        rank_stop=False,target_rank=32,maximum_workers=1,rss_bytes=2147483648,
        seed_seconds=180,geometry_wall_seconds=180,worker_wall_seconds=1200,replay_wall_seconds=1200,
        gp_sha256=cert.hashed(selection.scalar.GP),
        gate='The original17-only302 calibration plus guarded recovered-only follow-ups recovers24 and seven of the14 exceptional directions. Exact rational identities are checked only afterward. The million independent addresses and48 score-stratified record-scale choices are now frozen.',
        adaptive_policy=q['point_policy'],
        stop_policy='Every active wave finishes all49 boxes and independent cloud proofs. Continue only on certified rank gain, up to five total waves. A no-gain wave or a certified lower bound at least32 ends that curve. All48 fixed curves are attempted in selection order; no replacement or automatic larger population. A failed/censored stage stops this driver with preserved evidence.',
        prefix_policy='Every initial seed is the exact17 specialized generic sections. At wave1 admit any nonzero parity. Later require parity outside the subgroup used as input to the previous wave, with all older basis points retained as a prefix. The original generic17 is unchanged in metadata.',
        selection='Only the frozen48 equations from the640..703 and704..767-bit j-numerator strata. No public exceptional points,302 coordinate, known rank label, validation prime or catalogue input. Public comparisons occur only after discovery.',
        boundary='Prospective point and rank-lower-bound experiment. A low or unchanged finite bound is not an upper rank bound. No generic specialization injectivity is assumed: each seed must pass exact finite independence. Raw maps need no complete discriminant factorization.'))
    print('FROZEN48 record-scale fibres, at most11760 calibrated boxes',flush=True)

def campaign():
    p=cert.read(BATCH/'protocol.json');assert p['sources']==sources()
    assert all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items());return p

def configure(index,wave):
    global D,SEED,ROW
    candidate=campaign()['rows'][index];D=BATCH/candidate['id']/('wave-'+str(wave+1).zfill(2));SEED=D/'seed.json'
    ROW=cert.read(D/'protocol.json')['rows'][0]
    engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def protocol():
    p=cert.read(D/'protocol.json');top=campaign()
    assert p['campaign_sha256']==cert.hashed(BATCH/'protocol.json') and p['sources']==top['sources']
    assert p['seed_sha256']==cert.hashed(SEED);return p

def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
        if m>>ROW['mask_floor'] and m not in result:result.append(m)
    return result

def launch():
    p=campaign();out=BATCH/'ledger.json';assert not out.exists()
    ledger=dict(status='RUNNING',rows=[],completed_boxes=0);checkpoint(out,ledger)
    for index,candidate in enumerate(p['rows']):
        folder=BATCH/candidate['id'];folder.mkdir(exist_ok=False)
        entry=dict(id=candidate['id'],status='RUNNING',height_bin=candidate['height_bin'],stratum=candidate['stratum'],stages=[],waves=[])
        ledger['rows'].append(entry);checkpoint(out,ledger)
        def stage(name,cmd,seconds,dest=folder,stages=None,cwd=ROOT):
            stages=entry['stages'] if stages is None else stages
            s=run(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=dest/(name+'.log'),checkpoint_path=dest/(name+'.supervisor.json'),cwd=cwd)
            ok=s['outcome']=='completed' and s['returncode']==0
            stages.append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
            print(candidate['id'],dest.name,name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
            if not ok:
                ledger['status']=entry['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed/censored record-scale point stage')
        stage('initial-seed',[SAGE,str(CAS/'prepare_det1092_record_scale_seed.sage'),'--index',str(index)],p['seed_seconds'])
        original=cert.read(folder/'initial-seed.json');seed=original;floor=0
        for wave in range(p['maximum_waves_per_curve']):
            wdir=folder/('wave-'+str(wave+1).zfill(2));wdir.mkdir(exist_ok=False);rank=len(seed['points'])
            assert rank>floor and seed['points'][:17]==original['points']
            proof=seed['rank_certificate'];actual=checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],
                [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
            assert digest(actual)==digest(proof)
            checkpoint(wdir/'seed.json',seed)
            checkpoint(wdir/'seed-cloud.json',dict(curve=seed['curve'],points=seed['points'],signatures=proof['signatures'],rank_certificate=proof,rank_lower_bound=rank,independent_column_indices=list(range(rank))))
            checkpoint(wdir/'protocol.json',dict(p,campaign_sha256=cert.hashed(BATCH/'protocol.json'),seed_sha256=cert.hashed(wdir/'seed.json'),
                rows=[dict(id=wdir.name,initial_rank=rank,generic_dimension=17,mask_floor=floor)],candidate_id=candidate['id'],wave=wave+1))
            w=dict(id=wdir.name,status='RUNNING',initial_rank=rank,mask_floor=floor,stages=[]);entry['waves'].append(w);checkpoint(out,ledger)
            fresh=wdir/'seed-replay';fresh.mkdir(exist_ok=False)
            for src in [wdir/'seed-cloud.json',CAS/'verify_factor_free_rank.sage']:shutil.copy2(src,fresh/src.name)
            checkpoint(fresh/'protocol.json',dict(files={q.name:cert.hashed(q) for q in fresh.iterdir()},seconds=120,rss_bytes=p['rss_bytes']))
            stage('seed-replay',[SAGE,str(fresh/'verify_factor_free_rank.sage'),'--input',str(fresh/'seed-cloud.json')],120,wdir,w['stages'],fresh)
            stage('geometry',[SAGE,str(CAS/'prepare_det1092_record_scale_maps.sage'),'--index',str(index),'--wave',str(wave)],p['geometry_wall_seconds'],wdir,w['stages'])
            for name in ['worker','replay']:
                stage(name,[sys.executable,str(Path(__file__).resolve()),name,'--index',str(index),'--wave',str(wave)],p[name+'_wall_seconds'],wdir,w['stages'])
            prefix='det1092_record_scale_'+candidate['id'].replace('-','_')+'_'+wdir.name.replace('-','_')
            stage('certificates',[sys.executable,str(CAS/'certify_factor_free_exposure_v3.py'),'--run',str(wdir),'--protocol',str(wdir/'protocol.json'),'--prefix',prefix],1800,wdir,w['stages'])
            c=cert.read(wdir/'certification-ledger.json');assert c['status']=='PASS'
            assert all(v==c['rank_lower_bound'] for v in c['odd_modulus_ranks'].values())
            cloudpath=ART/(prefix+'_mod2_v1.json');cloud=cert.read(cloudpath)
            assert cloud['independent_points'][:rank]==seed['points']
            w.update(status='PASS',rank_lower_bound=c['rank_lower_bound'],completed_boxes=49,cloud_path=str(cloudpath.relative_to(ROOT)),cloud_sha256=cert.hashed(cloudpath))
            ledger['completed_boxes']+=49;entry['rank_lower_bound']=c['rank_lower_bound'];checkpoint(out,ledger)
            if c['rank_lower_bound']==rank:entry['stop_reason']='NO_CERTIFIED_GAIN';break
            if c['rank_lower_bound']>=32:entry['stop_reason']='CERTIFIED_AT_LEAST32';break
            floor=rank;seed=dict(original,points=cloud['independent_points'],rank_certificate=cloud['rank_certificate'])
        else:entry['stop_reason']='FIVE_WAVE_LIMIT'
        entry['status']='PASS';checkpoint(out,ledger)
        print('CURVE COMPLETE',candidate['id'],'rank >=',entry['rank_lower_bound'],entry['stop_reason'],flush=True)
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','launch','worker','replay']);ap.add_argument('--index',type=int);ap.add_argument('--wave',type=int);a=ap.parse_args()
    if a.stage in ('worker','replay'):configure(a.index,a.wave);getattr(engine,a.stage)()
    else:globals()[a.stage]()

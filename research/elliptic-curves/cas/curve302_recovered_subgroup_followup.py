#!/usr/bin/env python3
"""Bounded adaptive follow-up from the19 independently recovered302 directions."""
import argparse,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
BATCH=ROOT/'artifacts/local/elliptic-curves/curve302-recovered-subgroup-followup-v1'
ORIGINAL=ROOT/'artifacts/local/elliptic-curves/curve302-focused-point-exposure-v2/curve302-generic17/seed.json'
CLOUD=ART/'curve302_focused_curve302_generic17_mod2_v1.json'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
READS=set()

def install_guard():
    def guard(event,args):
        if event!='open' or not isinstance(args[0],(str,bytes)):return
        path=Path(args[0]).resolve()
        if not path.is_relative_to(ROOT/'artifacts'):return
        if not (path.is_relative_to(BATCH) or path in (ORIGINAL,CLOUD)):
            raise PermissionError('recovered-only control rejected artifact read '+str(path.relative_to(ROOT)))
        READS.add(str(path.relative_to(ROOT)))
    sys.addaudithook(guard)

# Geometry installs the guard before importing its arithmetic modules; command
# workers install it here, before importing any research implementation.
if __name__=='__main__' and len(sys.argv)>1 and sys.argv[1] in ('worker','replay'):
    install_guard()
import certify_compact_r17_candidates as cert
import retained_native19_trial_v3 as engine
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits

def sources():
    names=['curve302_recovered_subgroup_followup.py','prepare_curve302_recovered_followup.sage',
           'certify_factor_free_exposure_v3.py','verify_factor_free_rank.sage']
    return {**engine.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def freeze():
    assert not (BATCH/'protocol.json').exists()
    old=cert.read(ORIGINAL);cloud=cert.read(CLOUD)
    assert cloud['status']=='COMPLETE_DECLARED_FINITE_AUDIT' and cloud['rank_lower_bound']==19
    assert cloud['independent_points'][:17]==old['points'] and len(old['points'])==17
    checkpoint(BATCH/'protocol.json',dict(schema='elliptic-curves.curve302-recovered-subgroup-followup.v1',
        sources=sources(),inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in [ORIGINAL,CLOUD]},
        maximum_waves=4,charts=49,maximum_point_boxes=196,sample_size=2048,
        sample_domain='full11952-specialized-followup-v1',height=125000,seconds_per_chart=10,
        target_rank=32,rank_stop=False,geometry_wall_seconds=180,worker_wall_seconds=1200,
        replay_wall_seconds=1200,rss_bytes=2147483648,maximum_workers=1,
        gp_sha256=cert.hashed(Path('/usr/bin/gp')),
        gate='The generic-only302 control independently recovered19 from17. Exact finite independence implies centre parities outside the prior subgroup give distinct degree-two classes. Test newly available geometry before enlarging the parameter population.',
        adaptive_policy='Start from the complete19-point finite-certified basis of the generic-only control. Each wave freezes2048 sampled parities and49 largest computed norm centres before points. Require parity outside the subgroup used to generate the preceding wave (initially original17). Certify the complete cloud modulo2/3/5 with independent Sage rank and geometry. Continue only if its certified independent basis grows; preserve the old basis as a prefix. Stop after a no-gain wave, a certified bound at least31, or four waves. No partial-wave rank stop or timeout retry.',
        selection='Only points actually found from the generic17 input may enter later seeds. No public exceptional point, public31 seed, retrospective visibility coordinate, score, validation prime or new fibre.',
        isolation='Geometry and point/history workers reject repository artifact reads outside this campaign and the original generic17 seed/recovered19 cloud. The fixed source manifest hashes implementation files. Seed and final arithmetic proof processes are separately recorded.',
        boundary='Known-curve recovery control. No rational-span upper bound, exhaustive parity claim, guarantee of reaching31, or automatic sibling population. Public31 union checks, if performed, follow terminal point waves and never affect selection.'))
    print('FROZEN at most4 waves/196 boxes from recovered19',flush=True)

def campaign():
    p=cert.read(BATCH/'protocol.json')
    assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items())
    return p

def protocol():
    p=cert.read(D/'protocol.json');top=campaign()
    assert p['campaign_sha256']==cert.hashed(BATCH/'protocol.json')
    assert cert.hashed(SEED)==p['seed_sha256'] and p['sources']==top['sources']
    return p

def configure(index):
    global ROW,D,SEED
    D=BATCH/('wave-'+str(index+1).zfill(2));SEED=D/'seed.json'
    ROW=cert.read(D/'protocol.json')['rows'][0]
    engine.ROW,engine.D,engine.SEED=ROW,D,SEED;engine.protocol,engine.masks=protocol,masks

def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
        if m>>ROW['mask_floor'] and m not in result:result.append(m)
    return result

def launch():
    p=campaign();out=BATCH/'ledger.json';assert not out.exists()
    ledger=dict(status='RUNNING',waves=[]);checkpoint(out,ledger)
    cloud=cert.read(CLOUD);old=cert.read(ORIGINAL);floor=17
    for i in range(p['maximum_waves']):
        folder=BATCH/('wave-'+str(i+1).zfill(2));folder.mkdir(exist_ok=False)
        pts=cloud['independent_points'];rank=len(pts);proof=cloud['rank_certificate']
        actual=checked_rank(tuple(map(cert.F,cloud['curve'])),[tuple(map(cert.F,P)) for P in pts],
            [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        assert digest(actual)==digest(proof) and pts[:17]==old['points'] and rank>floor
        seed=dict(family='det1092',parameter='0',curve=cloud['curve'],points=pts,
                  generic_points=old['points'],rank_certificate=proof)
        checkpoint(folder/'seed.json',seed)
        checkpoint(folder/'seed-cloud.json',dict(curve=seed['curve'],points=pts,signatures=proof['signatures'],
            rank_certificate=proof,rank_lower_bound=rank,independent_column_indices=list(range(rank))))
        row=dict(id=folder.name,initial_rank=rank,generic_dimension=17,mask_floor=floor)
        checkpoint(folder/'protocol.json',dict(p,campaign_sha256=cert.hashed(BATCH/'protocol.json'),
            seed_sha256=cert.hashed(folder/'seed.json'),rows=[row],wave=i+1))
        configure(i);entry=dict(id=folder.name,status='RUNNING',initial_rank=rank,mask_floor=floor,stages=[])
        ledger['waves'].append(entry);checkpoint(out,ledger)
        def stage(name,cmd,seconds,cwd=ROOT):
            s=run(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=folder/(name+'.log'),
                checkpoint_path=folder/(name+'.supervisor.json'),cwd=cwd)
            ok=s['outcome']=='completed' and s['returncode']==0
            entry['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger)
            print(folder.name,name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
            if not ok:
                entry['status']=ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger)
                raise ArithmeticError('preserve failed adaptive stage')
        fresh=folder/'seed-replay';fresh.mkdir(exist_ok=False)
        for src in [folder/'seed-cloud.json',CAS/'verify_factor_free_rank.sage']:shutil.copy2(src,fresh/src.name)
        checkpoint(fresh/'protocol.json',dict(files={q.name:cert.hashed(q) for q in fresh.iterdir()},seconds=120,rss_bytes=p['rss_bytes']))
        stage('seed-replay',[SAGE,str(fresh/'verify_factor_free_rank.sage'),'--input',str(fresh/'seed-cloud.json')],120,fresh)
        stage('geometry',[SAGE,str(CAS/'prepare_curve302_recovered_followup.sage'),'--index',str(i)],p['geometry_wall_seconds'])
        for name in ['worker','replay']:
            stage(name,[sys.executable,str(Path(__file__).resolve()),name,'--index',str(i)],p[name+'_wall_seconds'])
        prefix='curve302_recovered_followup_'+folder.name.replace('-','_')
        stage('certificates',[sys.executable,str(CAS/'certify_factor_free_exposure_v3.py'),'--run',str(folder),
            '--protocol',str(folder/'protocol.json'),'--prefix',prefix],1800)
        c=cert.read(folder/'certification-ledger.json');assert c['status']=='PASS'
        assert all(v==c['rank_lower_bound'] for v in c['odd_modulus_ranks'].values())
        cloud=cert.read(ART/(prefix+'_mod2_v1.json'))
        assert cloud['independent_points'][:rank]==pts
        entry.update(status='PASS',rank_lower_bound=c['rank_lower_bound'],completed_boxes=49,
            cloud_path=str((ART/(prefix+'_mod2_v1.json')).relative_to(ROOT)),cloud_sha256=cert.hashed(ART/(prefix+'_mod2_v1.json')))
        checkpoint(out,ledger)
        if c['rank_lower_bound']==rank:ledger['stop_reason']='NO_CERTIFIED_GAIN';break
        if c['rank_lower_bound']>=31:ledger['stop_reason']='RECOVERED_AT_LEAST31';break
        floor=rank
    else:ledger['stop_reason']='FOUR_WAVE_LIMIT'
    ledger.update(status='PASS',rank_lower_bound=cloud['rank_lower_bound'],completed_boxes=49*len(ledger['waves']))
    checkpoint(out,ledger);print('TERMINAL',ledger['stop_reason'],ledger['rank_lower_bound'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','launch','worker','replay']);ap.add_argument('--index',type=int);a=ap.parse_args()
    if a.stage in ('worker','replay'):
        configure(a.index);getattr(engine,a.stage)();checkpoint(D/(a.stage+'-data-access.json'),sorted(READS))
    else:globals()[a.stage]()

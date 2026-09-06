#!/usr/bin/env python3
"""Run the frozen retained-pool comparison only after the completed comparison and retained selection."""
import argparse
import sys
import time
from pathlib import Path
import nearcut60v2_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits
ROOT=batch.ROOT; CAS=batch.CAS; D=batch.extension.D/'point-controller-v2'


def sources():
    names=['finish_nearcut60v2_mw16_experiment.py','nearcut60v2_mw16_pari_batch.py',
           'prepare_nearcut60v2_mw16_pari_batch.sage','verify_nearcut60v2_mw16_points.py',
           'report_nearcut60v2_mw16_experiment.py','replay_corrected60_mw16_geometry.py',
           'audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py']
    paths=[*(CAS/n for n in names),batch.extension.D/'controller/protocol.json',batch.extension.D/'protocol.json']
    return {**batch.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve comparison controller')
    batch.extension.protocol()
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.nearcut60v2-mw16-controller.v1',
        'sources':sources(),'wait_seconds':100000,'rss_bytes':4294967296,
        'stages':[
            ['freeze','nearcut60v2_mw16_pari_batch.py',['freeze'],120],
            ['maps','nearcut60v2_mw16_pari_batch.py',['maps'],4000],
            ['baselines','nearcut60v2_mw16_pari_batch.py',['baselines'],4000],
            ['points','nearcut60v2_mw16_pari_batch.py',['batch'],20000],
            ['verify','verify_nearcut60v2_mw16_points.py',['all'],43000],
            ['report','report_nearcut60v2_mw16_experiment.py',[],180],
            ['report-check','report_nearcut60v2_mw16_experiment.py',['--check'],180]],
        'scope':'The preceding corrected and matched trials are complete. Require frozen near-finalist selection and replay before this separate same-size retained-pool trial. All60 maps and independent16 baseline certificates precede all point attempts. Exactly43 generic charts per curve with identical125000 height,10sec/chart,600sec/curve,two workers. No rank stop/retry/refill. Exact histories, exposure geometry, full retained clouds mod2/3/5, provenance and censor-aware measured-yield reporting. Per-curve proof failures remain unresolved and do not discard allocated rows. No validation-driven adaptation or new parameter scan.'})


def launch():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('comparison sources changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve controller attempt')
    ledger={'status':'WAITING_FOR_ORIGINAL_PROOFS_AND_MATCHING','rows':[]};checkpoint(D/'ledger.json',ledger)
    deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            state=cert.read(batch.extension.D/'controller/ledger.json')['status']
            if state=='PASS_FROZEN60_RETAINED_SELECTION':break
            if state=='SELECTION_INCOMPLETE_NO_POINT_SEARCH':
                ledger['status']=state;checkpoint(D/'ledger.json',ledger);return
            if state not in ('WAITING_FOR_UNCHANGED_CORRECTED_TRIAL','RUNNING') or time.monotonic()>deadline:
                raise ArithmeticError('original/matching prerequisite failed or finite wait elapsed')
            time.sleep(5)
        if p['sources']!=sources():raise ArithmeticError('comparison sources changed during wait')
        batch.extension.completion_gate()
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for label,script,args,seconds in p['stages']:
            s=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),cwd=ROOT,
                  log_path=D/(label+'.log'),checkpoint_path=D/(label+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':label,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            print('NEARCUT60V2 CONTROLLER',label,ledger['rows'][-1]['status'],flush=True)
            if not ok:
                if label in ('maps','baselines'):
                    # These gates have accounted for all60 allocations but prohibit points.
                    for reporting_args in ([],['--check']):
                        suffix='gate-report-check' if reporting_args else 'gate-report'
                        r=run([sys.executable,str(CAS/'report_nearcut60v2_mw16_experiment.py'),*reporting_args],limits=Limits(180,p['rss_bytes']),cwd=ROOT,log_path=D/(suffix+'.log'),checkpoint_path=D/(suffix+'.supervisor.json'))
                        ledger['rows'].append({'name':suffix,'status':'PASS' if r['outcome']=='completed' and r['returncode']==0 else 'FAILED_OR_CENSORED','supervision':r});checkpoint(D/'ledger.json',ledger)
                raise ArithmeticError('fixed comparison stage failed/censored; retained without retry')
        ledger['status']='COMPLETE_FIXED_RETAINED_TRIAL_AND_ACCOUNTING';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

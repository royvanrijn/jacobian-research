#!/usr/bin/env python3
"""Bounded independent accounting audit after the frozen comparison terminates."""
import argparse,sys,time
from pathlib import Path
import audit_strata60_mw16_accounting as audit
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=audit.ROOT;D=audit.LOCAL/'retained-score-stratification-v1/accounting-controller'
UPSTREAM=audit.LOCAL/'retained-score-stratification-v1/point-controller'


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve independent audit protocol')
    paths=[Path(__file__).resolve(),Path(audit.__file__).resolve(),UPSTREAM/'protocol.json',
           ROOT/'elliptic-curves/cas/report_strata60_mw16_experiment.py']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':172800,'seconds_per_stage':300,'rss_bytes':2147483648,
        'scope':'After the frozen comparison terminates and its report replay passes, independently reconstruct every allocated row, measured cost, completed exposure, certified-gain binding and policy criterion twice from retained files. No point search, score calculation, cohort change, retry or follow-on sweep.'})


def launch():
    p=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve audit launch')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('audit source changed')
    ledger={'status':'WAITING_FOR_TERMINAL_COMPARISON','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            up=cert.read(UPSTREAM/'ledger.json');state=up['status']
            if state=='COMPLETE_FIXED_COMPARISON_AND_ACCOUNTING':break
            if state=='MATCHING_INCOMPLETE_NO_POINT_SEARCH':
                ledger['status']='NOT_RUN_MATCHING_INCOMPLETE';checkpoint(D/'ledger.json',ledger);return
            if state=='FAILED_OR_CENSORED':
                if any(r['name']=='gate-report-check' and r['status']=='PASS' for r in up['rows']):break
                raise ArithmeticError('upstream terminated without replayed accounting')
            if state not in ('WAITING_FOR_ORIGINAL_PROOFS_AND_MATCHING','RUNNING') or time.monotonic()>deadline:
                raise ArithmeticError('upstream failed or fixed wait expired')
            time.sleep(5)
        if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('audit source changed during wait')
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,args in [('build',[]),('check',['--check'])]:
            s=run([sys.executable,str(Path(audit.__file__).resolve()),*args],limits=Limits(p['seconds_per_stage'],p['rss_bytes']),cwd=ROOT,
                log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            if not ok:raise ArithmeticError('independent audit failed; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

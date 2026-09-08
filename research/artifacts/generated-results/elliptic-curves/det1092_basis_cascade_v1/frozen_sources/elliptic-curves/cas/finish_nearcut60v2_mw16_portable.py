#!/usr/bin/env python3
"""Package and replay terminal comparison point proofs after independent accounting."""
import argparse,sys,time
from pathlib import Path
import nearcut60v2_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.extension.D/'portable-controller-v2'
UPSTREAM=batch.extension.D/'accounting-controller-v2'


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve portable comparison protocol')
    preflight=batch.ART/'nearcut60v2_mw16_portability_preflight_v1.json'
    if cert.read(preflight)['status']!='PASS':raise ArithmeticError('isolated prospective input preflight required')
    paths=[Path(__file__).resolve(),CAS/'package_nearcut60v2_mw16_points.py',CAS/'verify_nearcut60v2_mw16_points_portable.py',
           CAS/'package_recorded_mod2_audit.py',UPSTREAM/'protocol.json',preflight]
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':172800,'rss_bytes':2147483648,
        'stages':[['package','package_nearcut60v2_mw16_points.py',300],['isolated-replay','verify_nearcut60v2_mw16_points_portable.py',12000]],
        'scope':'After the fixed comparison and independent accounting terminate successfully, preserve all60 allocated outcomes in a standalone archive and replay exact histories, geometry and mod2/3/5 point proofs for every certified row, at most240 isolated stages. Failed/unresolved rows stay explicitly unresolved. No new point search, budget change, retry, score, catalogue input or follow-on sweep. Accounting evidence is embedded, not recomputed by isolated point replay.'})


def launch():
    p=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve portable comparison launch')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('portable sources changed')
    ledger={'status':'WAITING_FOR_COMPARISON_AND_ACCOUNTING','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            state=cert.read(UPSTREAM/'ledger.json')['status']
            if state=='PASS':break
            if state=='NOT_RUN_MATCHING_INCOMPLETE':
                ledger['status']='NOT_RUN_MATCHING_INCOMPLETE';checkpoint(D/'ledger.json',ledger);return
            if state not in ('WAITING_FOR_TERMINAL_COMPARISON','RUNNING') or time.monotonic()>deadline:raise ArithmeticError('accounting failed or fixed wait elapsed')
            time.sleep(5)
        if cert.read(batch.extension.D/'point-controller-v2/ledger.json')['status']!='COMPLETE_FIXED_RETAINED_TRIAL_AND_ACCOUNTING':
            ledger['status']='NOT_RUN_PRESEARCH_GATE_FAILURE';checkpoint(D/'ledger.json',ledger);return
        if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('portable sources changed during wait')
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,script,seconds in p['stages']:
            s=run([sys.executable,str(CAS/script)],limits=Limits(seconds,p['rss_bytes']),cwd=ROOT,
                log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            if not ok:raise ArithmeticError('portable comparison failed; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

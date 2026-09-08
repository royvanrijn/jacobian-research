#!/usr/bin/env python3
"""Wait for the unchanged corrected trial, then freeze/replay retained-score matches."""
import argparse,sys,time
from pathlib import Path
import select_retained_mw16_score_strata as selection
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=selection.ROOT;D=selection.D/'controller'
UPSTREAM=selection.LOCAL/'corrected60-mw16-point-portable-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve post-corrected comparison controller')
    selection.protocol()
    paths=[Path(__file__).resolve(),Path(selection.__file__).resolve(),selection.D/'protocol.json',UPSTREAM/'protocol.json']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-score-strata-controller.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':90000,'stages':['select','check'],'seconds_per_stage':1800,
        'rss_bytes':2147483648,'scope':'Let the original corrected MW16 experiment and its182 isolated point-only proof stages finish unchanged. Then select and replay exactly twenty matched score-stratum triplets from existing retained candidates. No new parameter scan, scalar trace, validation input or point search. A matching failure is recorded without expanded pools or relaxed calipers.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('comparison controller changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve comparison controller attempt')
    ledger={'status':'WAITING_FOR_UNCHANGED_CORRECTED_TRIAL','rows':[]};checkpoint(D/'ledger.json',ledger)
    deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            state=cert.read(UPSTREAM/'ledger.json')['status']
            if state=='PASS':break
            if state not in ['WAITING_FOR_POINT_FINALIZATION','RUNNING'] or time.monotonic()>deadline:
                raise ArithmeticError('original trial failed/censored or fixed wait elapsed')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for stage in p['stages']:
            s=run([sys.executable,str(Path(selection.__file__).resolve()),stage],
                  limits=Limits(p['seconds_per_stage'],p['rss_bytes']),cwd=ROOT,
                  log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':stage,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
            checkpoint(D/'ledger.json',ledger)
            if not ok:raise ArithmeticError('fixed matching stage failed; no retry')
        ledger['status']=cert.read(selection.OUT)['status'];checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);args=p.parse_args();globals()[args.stage]()

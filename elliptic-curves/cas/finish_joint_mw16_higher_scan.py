#!/usr/bin/env python3
"""Wait for exact five-family cache gates, then scan the frozen higher slices."""
import argparse,sys,time
from pathlib import Path
import scan_joint_mw16_higher_annuli as scan
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=scan.ROOT;CAS=scan.CAS;D=scan.D/'controller'
UPSTREAM=ROOT/'artifacts/local/elliptic-curves/mw16-joint-cache-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher-slice controller')
    paths=[Path(__file__).resolve(),Path(scan.__file__),Path(scan.old.__file__),
           scan.old.D/'protocol.json',Path(scan.binary.__file__),UPSTREAM/'protocol.json',
           CAS/'verify_periodic_nagao_scanner.py',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    rows=scan.roster()
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'rows':rows,'primitive_population':sum(r['primitive_population'] for r in rows),
        'wait_seconds':18000,'rss_bytes':2147483648,
        'stages':[('prepare',180),('benchmark',300),('run',1200),('replay',600)],
        'scope':'Freeze twenty untouched higher MW16 slices before the complete cache and all40960 saved-score gates finish. After all gates pass, score every address through32749 before retaining4096 per slice. Require twenty complete signed-frame regressions, first-slice-per-band runtime gates and all81920 independent retained-score replays. Exactly191215782 primitive addresses in16384<H<=65536 and65536<H<=262144; no compact-region rescan or prior outer slice reuse. No public-record input, adaptive refill or point search.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()) or p['rows']!=scan.roster():raise ArithmeticError('frozen higher-slice scope changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve higher-slice launch')
    ledger={'status':'WAITING_FOR_FIVE_CACHE_PROOFS','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            status=cert.read(UPSTREAM/'ledger.json')['status']
            if status=='PASS':break
            if status not in ('RUNNING','WAITING_FOR_CACHE_BUILD') or time.monotonic()>deadline:raise ArithmeticError('cache gates failed/censored or fixed wait expired')
            time.sleep(5)
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,seconds in p['stages']:
            s=run([sys.executable,str(Path(scan.__file__)),name],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            print('HIGHER MW16 SCAN',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('higher-slice stage failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','launch']);a=parser.parse_args();globals()[a.stage]()

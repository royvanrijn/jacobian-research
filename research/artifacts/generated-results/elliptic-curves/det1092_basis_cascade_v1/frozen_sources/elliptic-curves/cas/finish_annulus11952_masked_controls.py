#!/usr/bin/env python3
"""Four fixed new-fibre masked diagnostics, with oracles opened after all searches."""
import argparse,sys
from pathlib import Path
import annulus11952_masked_controls as control
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=control.ROOT;CAS=control.CAS;D=control.D/'controller'
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve four-control controller')
    paths=[Path(__file__).resolve(),Path(control.__file__),CAS/'prepare_annulus11952_masked_controls.sage',CAS/'audit_annulus11952_masked_relations.sage',control.D/'protocol.json',CAS/'research_runtime/supervisor.py',CAS/'research_runtime/store.py']
    jobs=[('maps',control.SAGE,'prepare_annulus11952_masked_controls.sage',[],600),('masked-points',sys.executable,'annulus11952_masked_controls.py',['launch'],1800),('relations',control.SAGE,'audit_annulus11952_masked_relations.sage',[],600),('relations-check',control.SAGE,'audit_annulus11952_masked_relations.sage',['--check'],300)]
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'jobs':jobs,'rss_bytes':2147483648,'scope':'Exactly four frozen new11952 cohort indices0,16,32,48. Independently certify17 generic points, withhold section0, prepare twelve sampled deep-centre maps using only the sixteen-point metric block, execute48 fixed125000/10second boxes and replay their geometry. Open all four recorded oracles only after all masked attempts/replays complete. Bounded finite-quotient relation proposals count recovery only after exact rational equality with a nonzero withheld coefficient. This diagnostic adds no new parameter or rank result and cannot alter the running original64-fibre batch.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('fixed masked-controller source changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve masked-controller attempt')
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    for name,exe,script,args,seconds in p['jobs']:
        s=run([exe,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('NEW11952 MASKED CONTROL',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('fixed diagnostic failed/censored; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

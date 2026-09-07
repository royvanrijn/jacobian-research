#!/usr/bin/env python3
"""Supervised post-search certification, without further point exposure."""
import sys
from pathlib import Path
import mestre_parent_calibration as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

def main():
    d=batch.BATCH;out=d/'post-verification-ledger.json'
    if out.exists():raise FileExistsError('preserve independent certification ledger')
    if cert.read(d/'ledger.json')['status']!='PASS':raise ArithmeticError('terminal frozen search required')
    stages=[('cloud-audits',sys.executable,'audit_mestre_parent_calibration.py',[],900),
      ('exact-geometry',batch.SAGE,'replay_mestre_parent_calibration_geometry.sage',[],300),
      ('independent-rank',batch.SAGE,'verify_mestre_parent_calibration_rank.sage',[],120),
      ('input-span',batch.SAGE,'certify_mestre_input_span.sage',[],120),
      ('input-span-check',batch.SAGE,'certify_mestre_input_span.sage',['--check'],120),
      ('input-odd-ranks',sys.executable,'audit_mestre_seed_clouds.py',[],180)]
    ledger={'status':'RUNNING','stages':[]};checkpoint(out,ledger)
    for name,python,source,args,seconds in stages:
        s=run([python,str(batch.CAS/source),*args],limits=Limits(seconds,1610612736),
          log_path=d/(name+'.log'),checkpoint_path=d/(name+'.supervisor.json'),cwd=batch.ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['stages'].append({'name':name,'source':str((batch.CAS/source).relative_to(batch.ROOT)),
          'source_sha256':cert.hashed(batch.CAS/source),'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
        checkpoint(out,ledger);print(name,s['outcome'],s['returncode'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger)
            raise ArithmeticError('independent certification failed or censored')
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':main()

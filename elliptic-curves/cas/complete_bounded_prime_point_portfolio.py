#!/usr/bin/env python3
"""Certify every completed row of the fixed prospective point portfolio."""
import sys
from pathlib import Path
import bounded_prime_point_portfolio as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

def main():
 p=batch.protocol();d=batch.BATCH;out=d/'post-verification-ledger.json';assert cert.read(d/'ledger.json')['status']=='PASS' and not out.exists()
 checkpoint(d/'post-verification-protocol.json',dict(source_sha256=cert.hashed(Path(__file__)),verifier_sha256=cert.hashed(batch.CAS/'certify_factor_free_exposure_v3.py'),protocol_sha256=cert.hashed(d/'protocol.json'),maximum_workers=1,seconds_per_curve=1100,rss_bytes=2147483648))
 ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
 for row in p['rows']:
  if next(r for r in cert.read(d/'ledger.json')['rows'] if r['id']==row['id'])['status']=='SKIPPED_NO_PROVED_NEW_HISTORY_BOX':continue
  folder=d/row['id'];prefix='bounded_prime_point_'+row['id'].replace('-','_')
  cmd=[sys.executable,str(batch.CAS/'certify_factor_free_exposure_v3.py'),'--run',str(folder),'--protocol',str(d/'protocol.json'),'--prefix',prefix]
  s=run(cmd,limits=Limits(1100,p['rss_bytes']),log_path=folder/'certification-driver.log',checkpoint_path=folder/'certification-driver.supervisor.json',cwd=batch.ROOT)
  ok=s['outcome']=='completed' and s['returncode']==0
  ledger['rows'].append(dict(id=row['id'],status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger);print(row['id'],'all proofs',s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed or censored certification')
 ledger['status']='PASS';checkpoint(out,ledger);print('PASS complete prospective portfolio certification',flush=True)
if __name__=='__main__':main()

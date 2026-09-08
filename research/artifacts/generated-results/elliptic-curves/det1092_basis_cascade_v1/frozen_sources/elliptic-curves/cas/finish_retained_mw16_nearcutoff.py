#!/usr/bin/env python3
"""Bind and supervise read-only replay of the already frozen nearcutoff roster."""
import argparse,sys
from pathlib import Path
import select_retained_mw16_nearcutoff as selection
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=selection.ROOT;D=selection.D/'controller'


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve retained selection check protocol')
    selection.protocol();selection.completion_gate()
    paths=[Path(__file__).resolve(),Path(selection.__file__),selection.D/'protocol.json',selection.OUT]
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wall_seconds':120,'rss_bytes':1073741824,
        'scope':'Exact read-only replay of the frozen60 retained near-finalist selection. The first build and read-only replay already ran before this binding stage. No new candidate generation, score, point search or cohort modification.'})


def launch():
    p=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve bound selection replay')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen selection check source changed')
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    s=run([sys.executable,str(Path(selection.__file__).resolve()),'check'],limits=Limits(p['wall_seconds'],p['rss_bytes']),cwd=ROOT,log_path=D/'check.log',checkpoint_path=D/'check.supervisor.json')
    ok=s['outcome']=='completed' and s['returncode']==0
    ledger['rows']=[{'name':'check','status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s}]
    ledger['status']=cert.read(selection.OUT)['status'] if ok else 'FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger)
    if not ok:raise ArithmeticError('selection check failed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

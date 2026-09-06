#!/usr/bin/env python3
"""Relocate frozen prospective inputs and imports without reading point outcomes."""
import sys,shutil
from pathlib import Path
import strata60_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.extension.D/'portability-preflight'
OUT=batch.ART/'strata60_mw16_portability_preflight_v1.json'


def main():
    if (D/'protocol.json').exists() or OUT.exists():raise FileExistsError('preserve relocation preflight')
    p=batch.protocol();workspace=D/'workspace'
    paths={ROOT/n for n in p['sources']}
    paths.update([batch.D/'protocol.json',batch.extension.OUT,batch.extension.D/'protocol.json',
        batch.extension.D/'controller/ledger.json',CAS/'verify_strata60_mw16_points.py',
        CAS/'package_strata60_mw16_points.py',CAS/'verify_strata60_mw16_points_portable.py'])
    paths=dependencies(paths)
    hashes={str(q.relative_to(ROOT)):cert.hashed(q) for q in sorted(paths)}
    checkpoint(D/'protocol.json',{'sources':hashes,'script_sha256':cert.hashed(Path(__file__).resolve()),
        'wall_seconds':120,'rss_bytes':1073741824,
        'scope':'Read only frozen prospective point inputs and source dependency closure. Copy them to a separate root and validate relocated imports and all60 protocol bindings under isolated Python (-I). No point result, rank outcome, search, allocation change or proof-execution claim.'})
    for q in sorted(paths):
        target=workspace/q.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(q,target)
    code='import sys;sys.path.insert(0,'+repr(str(workspace/'elliptic-curves/cas'))+');import strata60_mw16_pari_batch as b;import verify_strata60_mw16_points; p=b.protocol();assert len(p["rows"])==60;assert p["stop_rank"] is None;print("RELOCATED FROZEN60 INPUTS AND PROOF IMPORTS PASS")'
    s=run([sys.executable,'-I','-c',code],limits=Limits(120,1073741824),cwd=workspace,
        log_path=D/'preflight.log',checkpoint_path=D/'preflight.supervisor.json')
    if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('relocation preflight failed')
    for name,h in hashes.items():
        if cert.hashed(workspace/name)!=h or cert.hashed(ROOT/name)!=h:raise ArithmeticError('immutable preflight member changed')
    checkpoint(OUT,{'schema':'elliptic-curves.strata60-portability-preflight.v1','status':'PASS',
        'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),D/'protocol.json',D/'preflight.supervisor.json')},
        'files_relocated':len(paths),'allocated_curves':60,'supervision':s,
        'claim_boundary':'Isolated Python imports and frozen prospective60 protocol bindings succeed in a separate workspace with no point outcomes copied. No point proof, matching rerun, score or search was executed.'})
    print('STRATA60 PROSPECTIVE PORTABILITY PREFLIGHT PASS',len(paths),'FILES',flush=True)


if __name__=='__main__':main()

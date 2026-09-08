#!/usr/bin/env python3
"""Bounded certificate-only union/visibility supplement after the focused run."""
import sys,shutil
from pathlib import Path
import curve302_focused_point_exposure_v2 as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

ROOT,CAS,ART=batch.ROOT,batch.CAS,batch.ART
D=batch.BATCH/'union-replay';cert=batch.cert

def main():
    assert cert.read(batch.BATCH/'run-ledger.json')['status']=='PASS'
    assert not (D/'protocol.json').exists()
    names=['audit_curve302_focused_union.py','certify_curve302_focused_union.py',
        'audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py',
        'verify_factor_free_rank.sage','search_observability.py']
    checkpoint(D/'protocol.json',dict(sources={str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names},
        input_run_ledger_sha256=cert.hashed(batch.BATCH/'run-ledger.json'),
        maximum_workers=1,rss_bytes=2147483648,seconds_per_stage=180,point_searches=0))
    ledger=dict(status='RUNNING',stages=[]);checkpoint(D/'ledger.json',ledger)
    m2=ART/'curve302_focused_union_mod2_v1.json';ml=ART/'curve302_focused_union_modl_v1.json'
    jobs=[('build',[sys.executable,str(CAS/names[0])]),
        ('mod2-check',[sys.executable,str(CAS/names[2]),'--check',str(m2)]),
        ('modl-build',[sys.executable,str(CAS/names[3]),'--input',str(m2),'--output',str(ml)]),
        ('modl-check',[sys.executable,str(CAS/names[3]),'--check',str(ml)])]
    def stage(name,cmd,cwd):
        s=run(cmd,limits=Limits(180,2147483648),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=cwd)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(D/'ledger.json',ledger)
        print('UNION',name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('preserve failed union audit')
    for name,cmd in jobs:stage(name,cmd,ROOT)
    fresh=D/'standalone';fresh.mkdir(exist_ok=False)
    for src in [m2,CAS/names[4]]:shutil.copy2(src,fresh/src.name)
    stage('standalone',[batch.SAGE,str(fresh/names[4]),'--input',str(fresh/m2.name)],fresh)
    c=cert.read(m2);odd=cert.read(ml)
    assert all(a['finite_column_rank']==c['rank_lower_bound'] for a in odd['audits'])
    ledger.update(status='PASS',rank_lower_bound=c['rank_lower_bound'],point_count=len(c['points']),
        files={str(p.relative_to(ROOT)):cert.hashed(p) for p in [m2,ml,ART/'curve302_focused_visibility_v1.json',fresh/m2.name,fresh/names[4]]})
    checkpoint(D/'ledger.json',ledger)

if __name__=='__main__':main()

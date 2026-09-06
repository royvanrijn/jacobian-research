#!/usr/bin/env python3
"""Seal the fixed endpoint trial only after all original completion/proof jobs pass."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BASE=ROOT/'artifacts/local/elliptic-curves/endpoint-specialized-parity-v1';D=BASE/'finalization-v1'

def prepare():
    out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve endpoint finalization protocol')
    names=['finalize_endpoint_point_evidence.py','replay_endpoint_specialized_geometry.sage','report_endpoint_specialized_trial.py','package_endpoint_point_trial.py','verify_endpoint_point_trial_bundle.py','report_endpoint_point_trial_portable.py','certify_new27_small_prime_saturation.py']
    paths=[*(CAS/n for n in names),BASE/'protocol.json',BASE/'completion-v1/protocol.json',BASE/'cloud-verification-v2/protocol.json',ART/'endpoint_section_spans_evidence_v1.json']
    checkpoint(out,{'schema':'elliptic-curves.endpoint-point-finalization.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_deadline_seconds':19000,'geometry_seconds':300,'aggregate_seconds':180,'package_seconds':1200,'portable_seconds':30000,'rss_bytes':3221225472,'scope':'Wait for both fixed21-curve completion and cloud-proof ledgers. Proceed only if both are PASS with all21 rows. Run exact geometry and aggregate build/replay, render the canonical finite outcome note, package immutable evidence, then perform all66 isolated checks and bind their result. No new point search or retry is launched, and no mathematical status or inventory is promoted automatically.'})

def write_note(d):
    out=ROOT/'elliptic-curves/notes/ENDPOINT_POINT_TRIAL_2026-09-06.md'
    high=sum(r['rank_lower_bound']>=22 for r in d['rows']);best=max(r['rank_lower_bound'] for r in d['rows'])
    summary=('No endpoint reaches the inventory threshold of 22.' if not high else f'{high} endpoints reach the inventory threshold of 22 and require separate promotion.')
    lines=['# Completed bounded point trial on the 21 omitted endpoints','',f"All 21 endpoint attempts and their exact history and point checks are complete. Of {d['declared_point_boxes']} declared boxes, {d['attempted_point_boxes']} were attempted and {d['completed_point_boxes']} completed their search bounds. The strongest certified lower bound is {best}. {summary}",'','The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/endpoint_specialized_trial_v1.json) binds the frozen roster, effective replay supervision, raw point records and full retained-cloud certificates modulo 2, 3 and 5.','', '| Endpoint | Initial section-span rank | Certified point bound | Completed / attempted boxes |','| --- | ---: | ---: | ---: |']
    for r in d['rows']:lines.append(f"| {r['id']} | {r['initial_rank']} | {r['rank_lower_bound']} | {r['completed_boxes']} / {r['attempted_charts']} |")
    lines += ['', '## Fixed geometry and exposure','', 'All 21 nonsingular endpoints were selected in the exact endpoint-audit order. Each used its separately certified independent subset of 11–17 section values. The first 256 distinct nonzero SHA256-derived parities supplied a fixed numerical sample; 12 representatives were chosen by computed specialized norms, with exact parity and rounded-norm transports. All 21 map files preceded every point search.','', 'The point budget was height 125,000 and ten seconds per PARI chart, with one point worker and a 300-second worker cap per curve. A curve could stop at 22 pending independent replay. No extra mask, parameter, chart, refill or larger point box was introduced.','', '## Preserved replay timeout','', 'The second endpoint returned 18,966 raw point records. Its original exact history replay exceeded 180 seconds. That timeout and the failed original launcher remain retained. A separately frozen completion protocol permitted one 600-second history replay for each remaining or previously timed-out completed worker. The second history then passed in 205.643 seconds. All point-search limits stayed fixed; no point worker was repeated.','', 'The full endpoint protocols and raw records remain under `artifacts/local/elliptic-curves/endpoint-specialized-parity-v1/`. Its `completion-v1/` and `cloud-verification-v2/` directories preserve the effective histories and point proofs.','', '## Claim boundary and replay','', 'The [exact section-span proof](ENDPOINT_SECTION_SPANS_2026-09-06.md) concerns only the transported generic sections. Its five proper known-section indices did not change this already frozen trial retrospectively. Trial point clouds contain the chosen independent seed plus every returned finite point; bounded outcomes do not prove whole-curve ranks, saturation or point absence.','', 'The same evidence supplement includes the [six own27-point subgroup saturation proofs at 2, 3 and 5](NEW27_SMALL_PRIME_SATURATION_2026-09-06.md). Those specific prime-saturation statements neither prove full saturation nor rule out independent directions beyond 27.','', 'The [evidence manifest](../../artifacts/generated-results/elliptic-curves/endpoint_point_trial_evidence_v1.json) and `verify_endpoint_point_trial_bundle.py` support 66 isolated checks. The [replay report](../../artifacts/generated-results/elliptic-curves/endpoint_point_trial_portable_replay_v1.json), when present, records whether all stages passed. No new point searches run during isolated replay.','']
    out.write_text('\n'.join(lines))

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen endpoint finalization sources differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve finalization execution')
    ledger={'status':'WAITING_FOR_FIXED_JOBS','rows':[]};checkpoint(out,ledger);deadline=time.monotonic()+p['wait_deadline_seconds']
    while True:
        states=[cert.read(BASE/f/'ledger.json') for f in ('completion-v1','cloud-verification-v2')]
        if any(d['status'] not in ('RUNNING','PASS') for d in states):raise ArithmeticError('upstream failure/censor retained; finalization stopped')
        if all(d['status']=='PASS' and len(d['rows'])==21 for d in states):break
        if time.monotonic()>deadline:raise TimeoutError('declared upstream wait deadline')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(out,ledger)
    jobs=[('geometry',['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'replay_endpoint_specialized_geometry.sage')],p['geometry_seconds'],BASE/'geometry.supervisor.json'),('aggregate-build',[sys.executable,str(CAS/'report_endpoint_specialized_trial.py')],p['aggregate_seconds'],D/'aggregate-build.supervisor.json'),('aggregate-check',[sys.executable,str(CAS/'report_endpoint_specialized_trial.py'),'--check'],p['aggregate_seconds'],D/'aggregate-check.supervisor.json'),('package',[sys.executable,str(CAS/'package_endpoint_point_trial.py')],p['package_seconds'],D/'package.supervisor.json'),('isolated-replay',[sys.executable,str(CAS/'verify_endpoint_point_trial_bundle.py')],p['portable_seconds'],D/'isolated-replay.supervisor.json'),('portable-report',[sys.executable,str(CAS/'report_endpoint_point_trial_portable.py')],p['aggregate_seconds'],D/'portable-report.supervisor.json')]
    for name,command,seconds,path in jobs:
        if name=='package':write_note(cert.read(ART/'endpoint_specialized_trial_v1.json'))
        if path.exists():raise FileExistsError('preserve single finalization stage')
        s=run(command,limits=Limits(seconds,p['rss_bytes']),log_path=path.with_suffix('.log'),checkpoint_path=path,cwd=ROOT);ok=s['outcome']=='completed' and s['returncode']==0
        ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger);print('ENDPOINT FINALIZATION',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:raise ArithmeticError('finalization stage failed or censored; no retry')
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()

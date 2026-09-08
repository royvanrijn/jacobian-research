#!/usr/bin/env python3
"""One bounded cached-score benchmark after the full11952 cache replay closes."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
import build_extended_projective_trace_cache_11952 as parent
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=parent.ROOT;CAS=parent.CAS;D=ROOT/'artifacts/local/elliptic-curves/retained-extended-cache-benchmark-controller-v1'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve one cached-score controller')
    p=parent.protocol();paths=[Path(__file__).resolve(),CAS/'benchmark_retained_extended_cache.py',CAS/'newfamily/score_retained_projective_cache.cpp',ROOT/'elliptic-curves/tests/test_retained_projective_cache_scorer.py',parent.D/'protocol.json']
    if len(p['primes'])!=2948 or p['projective_rows']!=52989620:raise ArithmeticError('fixed full11952 source required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-extended-cache-controller.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_deadline_seconds':6000,'rss_bytes':4294967296,'jobs':[{'name':n,'seconds':s} for n,s in [('prepare',60),('encode',600),('encoding-check',600),('compile',60),('score',60),('score-check',60)]],'scope':'Wait only for the already frozen2948-prime11952 cache and its exact replay. If it passes, freeze a binary-cache benchmark on exactly967 previously scored outer candidates, encode and byte-replay the cache, compile the scorer, and run/check one exact score comparison. Any failure or censoring stops without retries. No full parameter scan, new finalist, point search or automatic next campaign.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen score-controller source changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve cached-score controller')
    ledger={'status':'WAITING_FOR_CACHE_REPLAY','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_deadline_seconds']
    while True:
        for label in ['run','check']:
            path=parent.D/(label+'.supervisor.json')
            if path.exists():
                s=cert.read(path)
                if s['outcome']!='running' and (s['outcome']!='completed' or s['returncode']!=0):raise ArithmeticError('upstream cache failure retained; score benchmark stopped')
        path=parent.D/'check.supervisor.json'
        if path.exists() and cert.read(path)['outcome']=='completed':break
        if time.monotonic()>deadline:raise TimeoutError('fixed cache wait deadline')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
    for job in p['jobs']:
        name=job['name'];s=run([sys.executable,str(CAS/'benchmark_retained_extended_cache.py'),name],limits=Limits(job['seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('CACHED SCORE BENCHMARK',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('cached score benchmark stopped; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','launch']);v=a.parse_args();globals()[v.stage]()

#!/usr/bin/env sage-python
"""Checkpointed top-64 full-prime continuation; at most four bounded workers.

Reuses completed measurements by exact parameter and retains all run outcomes.
Public equations only exclude already-catalogued isomorphism classes after the
fullscore order is frozen; public ranks/points do not select fresh candidates.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys,time

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
full=SourceFileLoader('batch_fullscore',str(CAS/'compact_r17_fullscore.sage')).load_module()
m=full.compact
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
from certify_compact_r17_candidates import isomorphic, DATABASE
from elliptic_candidate_record import weierstrass_invariants

DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-top64-v1'


def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    old=m.read(full.DIRECTORY/'population.json')
    protocol={**m.read(full.DIRECTORY/'protocol.json'),
        'schema':'elliptic-curves.compact-r17-top64-protocol.v1','sources':m.provenance(),
        'selection_source':{str(p.relative_to(ROOT)):m.hashed(p) for p in
            (Path(__file__).resolve(),CAS/'compact_r17_fullscore.sage',CAS/'compact_r17_wide.sage',full.wide.SCANNER)},
        'finalists':64,'parent_population_sha256':m.hashed(full.DIRECTORY/'population.json'),
        'admission_prime_bound':251,'max_workers':4,'worker_wall_limit_seconds':300,
        'worker_rss_limit_bytes':1610612736,'global_submission_wall_limit_seconds':3600,
        'novelty_filter':'After freezing the top-64 order, skip exact Q-isomorphism matches in the pinned 584-curve equation snapshot. Never read public points or ranks for selection.',
        'stop_after_certified_rank_lower_bound':28,
        'gate':'Full-prime selector and fixed deep-chart detector independently recover known ranks 25 and 26; fresh compact specializations have certified ranks 22,23,24. Short-prefix hard selection is removed.'}
    pp=directory/'protocol.json'
    if pp.exists() and m.read(pp)!=protocol:raise ArithmeticError('batch protocol changed')
    checkpoint(pp,protocol)
    pop={'protocol_hash':m.identity(protocol),'candidate_count':old['candidate_count'],
        'retained_candidates':old['retained_candidates'],'finalists':old['retained_candidates'][:64]}
    checkpoint(directory/'population.json',pop)
    equations=[{'id':r['id'],'ainvs':r['ainvs']} for r in m.read(DATABASE)['curves']]
    byj={}
    for r in equations:
        inv=weierstrass_invariants(r['ainvs']);j=inv['c4']**3/inv['discriminant'];byj.setdefault(j,[]).append(r)
    family=m.load_q12o5867_data(m.MODEL,m.SECTIONS);prior={}
    for d in (m.DIRECTORY,full.wide.DIRECTORY,full.DIRECTORY):
        for path in d.glob('*/result.json'):
            r=m.read(path)
            if r['status']=='COMPLETE':prior[r['parameter']]=path
    roster=[]
    for i,r in enumerate(pop['finalists']):
        t=m.F(r['parameter']);s=m.evaluate_projective_specialization(family,t.numerator,t.denominator)
        inv=weierstrass_invariants(s.model);j=inv['c4']**3/inv['discriminant']
        matches=[x['id'] for x in byj.get(j,[]) if isomorphic(s.model,x['ainvs'])]
        path=prior.get(r['parameter'])
        roster.append({'index':i,'parameter':r['parameter'],'snapshot_matches':matches,
            'prior_complete_measurement':str(path.relative_to(ROOT)) if path else None,
            'prior_measurement_sha256':m.hashed(path) if path else None})
    checkpoint(directory/'roster.json',{'protocol_hash':m.identity(protocol),'database_sha256':m.hashed(DATABASE),'rows':roster})
    print('ROSTER',len(roster),'known',sum(bool(r['snapshot_matches']) for r in roster),'prior',sum(bool(r['prior_complete_measurement']) for r in roster),flush=True)


def run(directory):
    protocol=m.read(directory/'protocol.json');roster=m.read(directory/'roster.json')
    if roster['protocol_hash']!=m.identity(protocol):raise ArithmeticError('roster binding changed')
    start=time.monotonic();ledger=[];stop=directory/'TARGET28.json'
    def job(row):
        out={'index':row['index'],'parameter':row['parameter']}
        if row['snapshot_matches']:return {**out,'status':'ALREADY_IN_PINNED_DATABASE','matches':row['snapshot_matches']}
        if row['prior_complete_measurement']:
            path=ROOT/row['prior_complete_measurement']
            if m.hashed(path)!=row['prior_measurement_sha256']:raise ArithmeticError('prior completed measurement changed')
            r=m.read(path)
            return {**out,'status':'REUSED_COMPLETE_MEASUREMENT','rank_lower_bound':r['rank_lower_bound'],'path':str(path.relative_to(ROOT))}
        if stop.exists() or time.monotonic()-start>=3600:return {**out,'status':'NOT_DISPATCHED_BUDGET_OR_TARGET_STOP'}
        path=directory/f'candidate-{row["index"]:02}'/'result.json'
        if path.exists() and m.read(path)['status']=='COMPLETE':
            r=m.read(path);return {**out,'status':'COMPLETE','rank_lower_bound':r['rank_lower_bound'],'path':str(path.relative_to(ROOT))}
        try:
            result=capture(['sage','-python',str(CAS/'compact_r17_fullscore.sage'),'run','--directory',str(directory),'--index',str(row['index'])],
                limits=Limits(300,1610612736),log_path=directory/f'worker-{row["index"]:02}.log')
            out['supervision']=result.supervision;out['status']='COMPLETE'
        except Exception as ex:out.update(status='WORKER_STOPPED',error=type(ex).__name__+': '+str(ex))
        if path.exists():
            r=m.read(path);out.update(rank_lower_bound=r['rank_lower_bound'],completed_chart_records=len(r['charts']),path=str(path.relative_to(ROOT)))
            if r['rank_lower_bound']>=28:checkpoint(stop,out)
        print('BATCH',out['index'],out['parameter'],out.get('rank_lower_bound'),out['status'],flush=True)
        return out
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending=[pool.submit(job,row) for row in roster['rows']]
        for future in as_completed(pending):
            ledger.append(future.result());checkpoint(directory/'ledger.json',{'protocol_hash':m.identity(protocol),
                'status':'RUNNING','rows':sorted(ledger,key=lambda r:r['index'])})
    checkpoint(directory/'ledger.json',{'protocol_hash':m.identity(protocol),'status':'COMPLETE_DECLARED_BATCH',
        'rows':sorted(ledger,key=lambda r:r['index']),'claim_boundary':'Worker stops, skipped known curves, and unrun rows are not rank exclusions.'})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);a=p.parse_args()
    (prepare if a.stage=='prepare' else run)(a.directory)


if __name__=='__main__':main()

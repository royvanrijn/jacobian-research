#!/usr/bin/env python3
"""Bind the completed fresh outer selection and fixed60 point outcomes."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import fresh60_mw16_pari_batch as batch
from research_runtime.store import checkpoint
ROOT=batch.ROOT;ART=batch.ART;D=batch.D;CONTROL=ROOT/'artifacts/local/elliptic-curves/fresh60-mw16-controller-v1'
INPUT=ART/'fresh60_mw16_results_v1.json';OUT=ART/'fresh60_mw16_experiment_v1.json'
def expected():
    p=batch.protocol();ledger=cert.read(CONTROL/'ledger.json');source=cert.read(INPUT);odd=cert.read(D/'odd-cloud-audit/ledger.json')
    if ledger['status']!='PASS' or odd['status']!='PASS' or len(source['curves'])!=60 or len(odd['rows'])!=60:raise ArithmeticError('all selection, maps, points, histories and clouds required')
    required=['freeze','maps','points','verify','odd','geometry']
    if [r['name'] for r in ledger['rows']]!=required:raise ArithmeticError('complete controller stages differ')
    paths={Path(__file__).resolve(),CONTROL/'protocol.json',CONTROL/'ledger.json',INPUT,D/'protocol.json',D/'ledger.json',D/'verification-ledger.json',D/'odd-cloud-audit/ledger.json',batch.extension.OUT,batch.extension.D/'fresh-validation.json'}
    for r in ledger['rows']:
        s=r['supervision'];path=Path(s['log']);paths.add(path)
        if r['status']!='PASS' or s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(path)!=s['log_sha256']:raise ArithmeticError('successful stage transcript differs')
    rows=[];stronger=[];validation={r['id']:r for r in cert.read(batch.extension.D/'fresh-validation.json')['rows']}
    for r,o in zip(source['curves'],odd['rows']):
        if r['id']!=o['id']:raise ArithmeticError('cloud roster differs')
        cloudpath=ROOT/o['input'];oddpath=ROOT/o['output'];cloud=cert.read(cloudpath);paths.update([cloudpath,oddpath])
        if cert.hashed(cloudpath)!=o['input_sha256'] or cert.hashed(oddpath)!=o['output_sha256'] or cloud['rank_lower_bound']!=r['rank_lower_bound']:raise ArithmeticError('bound source differs')
        if max(o['odd_lower_bounds'].values())>r['rank_lower_bound']:stronger.append(r['id'])
        rows.append({'id':r['id'],'family':r['family'],'band':r['band'],'parameter':r['parameter'],'rank_lower_bound':r['rank_lower_bound'],
            'odd_modulus_lower_bounds':o['odd_lower_bounds'],'retained_points':len(cloud['points']),
            'completed_boxes':r['completed_boxes'],'attempted_charts':r['attempted_charts'],
            'search_status':r['search_status'],'icarm_matches':r['icarm_matches'],'previous_matches':r['previous_matches'],
            'fresh_validation_units':validation[r['id']]['validation_units']})
    return {'schema':'elliptic-curves.fresh60-mw16-experiment.v1','status':'PASS_STRONGER_ODD_BOUNDS_REQUIRE_REVIEW' if stronger else 'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'rows':rows,
        'stronger_odd_prime_bounds':stronger,'completed_point_boxes':sum(r['completed_boxes'] for r in rows),
        'attempted_charts':sum(r['attempted_charts'] for r in rows),'maximum_declared_boxes':2580,
        'retained_point_witnesses':sum(r['retained_points'] for r in rows),
        'claim_boundary':'Ten fresh signed denominator slices in4096<H<=16384, disjoint from the earlier outer cohort, five MW16 families with equal retention and point budgets. All40960 retained extended scores precede60 distinct finalists; wholly disjoint validation never changes selection. Up to2580 generic16-only point charts with all maps before points. Catalogue and977 prior-equation comparison follows terminal exact point proofs. Matching finite bounds do not prove exact rank, saturation or absence. No universal novelty, rank density, calibrated sensitivity or isolated portability claim.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('fresh60 aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve fresh60 aggregate')
        checkpoint(OUT,d)
    print('FRESH60 OUTCOMES',d['status'],d['completed_point_boxes'],[(r['id'],r['rank_lower_bound']) for r in d['rows'] if r['rank_lower_bound']>=22],flush=True)

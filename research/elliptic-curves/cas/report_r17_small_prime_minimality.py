#!/usr/bin/env python3
"""Universal non13 minimality at primes5through997 for primitive R17 parameters."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_r17_scaling_prime_support as support
import classify_r17_other_small_prime_scalings as local
import verify_r17_other_small_prime_scalings as independent
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint

ROOT=support.ROOT
ART=support.ART
OUT=ART/'r17_small_prime_minimality_v1.json'


def expected():
    a=support.expected();b=local.expected()
    if cert.read(support.OUT)!=json.loads(json.dumps(a)) or cert.read(local.OUT)!=json.loads(json.dumps(b)):raise ArithmeticError('resultant or residue proof differs')
    independent.main();primes=[p for p in _primes_up_to(997) if p>=5 and p!=13]
    pairs={(r['family'],r['prime']):r for r in b['rows']};rows=[]
    if len(pairs)!=26 or any(r['status']!='NO_REMOVABLE_SCALE' for r in pairs.values()):raise ArithmeticError('all26 exact local exclusions required')
    for r in a['rows']:
        if not int(r['homogeneous_resultant']):raise ArithmeticError('nonzero homogeneous resultant required')
        valuations=dict(r['trial_prime_factors']);direct=[];residue=[]
        for p in primes:
            if valuations.get(p,0)<4:direct.append(p)
            else:
                if (r['family'],p) not in pairs:raise ArithmeticError('unaccounted scaling-prime candidate')
                residue.append(p)
        rows.append({'family':r['family'],'excluded_by_resultant_valuation':direct,'excluded_by_complete_residue_tree':residue,'all_non13_primes_5_through997_minimal':True})
    paths=[Path(__file__).resolve(),Path(support.__file__).resolve(),Path(local.__file__).resolve(),Path(independent.__file__).resolve(),support.OUT,local.OUT]
    for path in [support.D/'independent.supervisor.json',local.D/'independent.supervisor.json']:
        s=cert.read(path);paths.append(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('independent exact verification required')
    return {'schema':'elliptic-curves.r17-small-prime-minimality.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'primes_checked':primes,'prime_family_pairs':6*len(primes),'pairs_excluded_by_residue_trees':26,'pairs_excluded_by_resultant_valuation':6*len(primes)-26,'claim_boundary':'For every primitive integer parameter pair in each of the six compact R17 families, no prime5through997 other than13 admits p^4|A_h and p^6|B_h. For nonsingular fibres this proves that the displayed integral short model is minimal at each such prime. The resultant adjugate makes valuation at least4 necessary, and all26 surviving candidate family/prime pairs have complete independent finite-residue exclusions. The separately certified13 classification is unchanged. Primes2,3, all primes above997, exact conductor, whole-curve rank, improved population density and new curves are not settled by this theorem.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();d=expected()
    if args.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('universal small-prime minimality report differs')
    else:
        if OUT.exists():raise FileExistsError('preserve universal small-prime proof')
        checkpoint(OUT,d)
    print('UNIVERSAL SMALL-PRIME MINIMALITY',d['prime_family_pairs'],'FAMILY/PRIME PAIRS PASS',flush=True)

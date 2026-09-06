#!/usr/bin/env python3
"""Concrete first26-point candidate, with novelty explicitly awaiting batch closure."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as arithmetic
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/compact192-r17-pari-v1/07ca9-001';INPUT=ART/'compact192_r17_07ca9_001_mod2_v1.json';OUT=ART/'compact192_first26_candidate_v1.json';SAGE=ART/'compact192_first26_candidate.sage'

def expected():
    cloud=cert.read(INPUT);v=cert.read(D/'verification.json');raw=cert.read(D/'result.json')
    if v['status']!='PASS' or v['input_sha256']!=cert.hashed(D/'result.json') or v['cloud_sha256']!=cert.hashed(INPUT) or cloud['rank_lower_bound']!=26 or raw['family']!='07ca9':raise ArithmeticError('fixed verified26 candidate required')
    short=tuple(map(cert.F,cloud['curve']));model=arithmetic.integral(short);inv=cert.weierstrass_invariants(model);points=[]
    for P in cloud['independent_points']:
        x,y=map(cert.F,P);X=x-inv['b2']/12;points.append((X,y-(model[0]*X+model[2])/2))
    if short!=(0,0,0,-inv['c4']/48,-inv['c6']/864) or any(not cert.is_on_weierstrass_curve(model,P) for P in points):raise ArithmeticError('exact integral model transport differs')
    back=[(x+inv['b2']/12,y+(model[0]*x+model[2])/2) for x,y in points];old=cloud['rank_certificate'];proof=checked_rank(short,back,[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
    if [list(map(str,P)) for P in back]!=cloud['independent_points'] or json.loads(json.dumps(proof))!=old:raise ArithmeticError('transported independent subgroup differs')
    minimal=arithmetic.minimality(model)
    sources={**arithmetic.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),INPUT,D/'verification.json',D/'result.json']}}
    return {'schema':'elliptic-curves.compact192-first26-candidate.v1','status':'PASS','sources':sources,'id':'07ca9-001','family':raw['family'],'parameter':raw['parameter'],'minimal_curve':list(map(str,model)),'points':[list(map(str,P)) for P in points],'discovery_curve':cloud['curve'],'discovery_points':cloud['independent_points'],'rank_lower_bound':26,'rank_certificate_on_discovery_curve':proof,'minimality':minimal,'novelty_status':'NOT_YET_CHECKED_AGAINST_CATALOGUE_FOR_THIS_BATCH','claim_boundary':'A concrete globally minimal integral model and26 exactly independent rational points from the completed first26-point candidate in the ongoing compact192 batch. Point and model proofs are complete; catalogue comparison and inventory promotion await terminal batch replay. This is not a new-curve, exact-rank, record, conductor or universal-novelty assertion.'}

def export(d):
    return '# Certified rank lower bound26. Catalogue novelty remains unverified for this batch.\n# Exact point-independence and minimality proof: compact192_first26_candidate_v1.json.\nfrom sage.all import QQ, EllipticCurve\nE = EllipticCurve(QQ, '+repr(d['minimal_curve'])+')\npoints = [E([QQ(x), QQ(y)]) for x, y in '+repr(d['points'])+']\nassert len(points) == 26\nprint(E)\nprint("26 exactly certified independent rational points; novelty pending")\n'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected();text=export(d)
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)) or SAGE.read_text()!=text:raise ArithmeticError('candidate model/point export differs')
    else:
        if OUT.exists() or SAGE.exists():raise FileExistsError('preserve first26 candidate evidence')
        checkpoint(OUT,d);SAGE.write_text(text)
    print('EXACT MINIMAL26 CANDIDATE',d['parameter'],d['minimality']['invariant_gcd'],d['novelty_status'],flush=True)

#!/usr/bin/env sage-python
"""A bounded local conductor upper bound for the new minimal rank27 curve."""
import sys
from pathlib import Path
from sage.all import QQ,ZZ,EllipticCurve,prime_range
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/paired-rank27-conductor-bound-v1';INPUT=ROOT/'artifacts/generated-results/elliptic-curves/new_paired_rank27_proof_v1.json';DB=ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json'

def main():
    p=cert.read(D/'protocol.json')
    if p['source_sha256']!=cert.hashed(Path(__file__).resolve()) or p['input_sha256']!=cert.hashed(INPUT) or p['database_sha256']!=cert.hashed(DB):raise ArithmeticError('bounded conductor input differs')
    if (D/'result.json').exists():raise FileExistsError('preserve bounded conductor output')
    d=cert.read(INPUT);E=EllipticCurve([QQ(a) for a in d['minimal_curve']]);delta=abs(ZZ(E.discriminant()));remaining=delta;known=ZZ(1);rows=[]
    if E.c4().gcd(E.c6())!=16 or delta.valuation(2)!=8:raise ArithmeticError('global minimality criterion changed')
    for q in prime_range(2,10001):
        v=int(remaining.valuation(q))
        if not v:continue
        remaining//=q**v
        if q==2:ld=E.local_data(q,proof=True);f=int(ld.conductor_valuation());kind=str(ld.kodaira_symbol())
        else:
            if E.c4()%q==0:raise ArithmeticError('odd bad prime not multiplicative')
            f=1;kind='multiplicative'
        known*=q**f;rows.append({'prime':str(q),'discriminant_valuation':v,'conductor_valuation':f,'reduction':kind})
    upper=known*remaining;db=cert.read(DB)['curves'];listed=sorted((ZZ(r['conductor']),r['id']) for r in db if r['rank_lower_bound']>=27 and str(r['conductor']).isdigit())
    checkpoint(D/'result.json',{'status':'COMPLETE_BOUNDED_LOCAL_AUDIT','sources':{str(path.relative_to(ROOT)):cert.hashed(path) for path in (Path(__file__).resolve(),INPUT,DB,D/'protocol.json')},'curve':d['minimal_curve'],'discriminant':str(E.discriminant()),'local_data':rows,'unprocessed_cofactor':str(remaining),'conductor_upper_bound':str(upper),'exact_conductor':remaining==1,'listed_rank27_minimum':{'id':listed[0][1],'conductor':str(listed[0][0])},'upper_bound_beats_listed_minimum':bool(upper<listed[0][0]),'missing_rank27_conductors':[r['id'] for r in db if r['rank_lower_bound']>=27 and not r['conductor']],'claim_boundary':'The conductor divides this upper bound. All residual primes have multiplicative reduction; the cofactor is not factored or asserted prime. A larger upper bound does not exclude a smaller actual conductor. No exact conductor or conductor record follows.'})
    print('NEW27 CONDUCTOR UPPER BOUND',len(str(upper)),'digits; residual',len(str(remaining)),'digits; beats listed',upper<listed[0][0],flush=True)
if __name__=='__main__':main()

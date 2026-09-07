#!/usr/bin/env python3
"""Fixed literature-motivated split-node bonus on saved retained models."""
import argparse
from pathlib import Path
from math import log
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'higher-split-node-score-v1'
OLD=LOCAL/'higher32768-r17-extended-v1/result.json';PRODUCT=LOCAL/'higher-r17-product-score-v1/result.json';NEW=LOCAL/'higher32768-product-first-extended-v1/result.json';EXCLUDED=[LOCAL/n/'protocol.json' for n in ('higher24-r17-pari-v1','product24-r17-pari-v1','productfirst24-r17-pari-v1','scaled13_24-r17-pari-v1')]
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',OLD,PRODUCT,NEW,*EXCLUDED)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve split-node scoring protocol')
    gate=ART/'product_first_portable_replay_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('completed unmodified product-first experiment required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-split-node-score.v1','sources':sources(),'gate_sha256':cert.hashed(gate),'prime_bound':65521,'selection_prime_bound':32749,'bonus_constant':[7,5],'per_family':4,'source_reference':'https://arxiv.org/pdf/2003.00077 section7 question2','gate':'The two unmodified product-score experiments complete without a high-rank addition. The cited primary source reports a heuristic bonus log(c*(p-1)/p) at split multiplicative primes, with1.2<=c<=1.68. That feature is absent from both completed product experiments. Test the fixed convenient value c=7/5 on the saved10482-address union, without fitting c to point results. The concurrent13-scaling cohort is excluded uniformly by its already frozen addresses, without waiting for or reading its outcomes.','formula':'Sum the existing quantized product score at good displayed primes and round(log(7*(p-1)/(5*p))*10^12) at each split-nodal scaled display. Primes through32749 select; higher primes through65521 are validation only, excluded from ties. No Tamagawa-number weighting or inferred rank is used.','local_algebra':'For each prime5..65521 dividing the integer short-model discriminant, divide A by p^4 and B by p^6 repeatedly where integral. If the scaled display is nonsingular, fail rather than silently omit a required good-prime trace. If A is a unit, the double root r=-3B/(2A) satisfies A=-3r^2,B=2r^3; the node splits precisely when3r is a nonzero square. Record all residues and the Euler-criterion value. A=0 gives a cuspidal display and no bonus. No complete factorization or conductor computation.','selection':'Merge the original S1-retained and product-first-retained pools by family and rational parameter, verify overlap score/model agreement, and exclude the union of addresses in four frozen higher-height point cohorts. Choose four per family by the augmented product score, existing good-prime count, denominator and signed numerator. No point, measured rank, catalogue label, public parameter or result-dependent refill enters the ordering.','limits':{'maximum_models':10482,'wall_seconds':180,'rss_bytes':1073741824,'workers':1,'checkpoint_rows':128},'future_scope':'This computes local algebra, a fixed heuristic score and proposed addresses only. A point cohort requires a separate frozen protocol after complete replay. No new trace or parameter scan is performed.','boundaries':'The literature heuristic is not a rank-prediction theorem and was developed on different families. Exact nodal reduction and numerical score replay do not establish point existence, rank upper bounds, Tamagawa factors, saturation, conductor or universal novelty. The union is still a retained subset of the original finite population.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen split-node inputs changed')
    return p
def pool():
    old=cert.read(OLD)['rows'];product={(r['family'],r['parameter']):r for r in cert.read(PRODUCT)['rows']};rows={}
    for r in old:
        key=r['family'],r['parameter'];q=product[key];rows[key]={k:r[k] for k in ('family','parameter','numerator','denominator','model','combined_good')};rows[key].update(product_selection_units=q['product_selection_units'],product_validation_units=q['product_validation_units'])
    for r in cert.read(NEW)['rows']:
        key=r['family'],r['parameter'];q={k:r[k] for k in ('family','parameter','numerator','denominator','model','combined_good')};q.update(product_selection_units=r['combined_selection_units'],product_validation_units=r['validation_units'])
        if key in rows and rows[key]!=q:raise ArithmeticError('overlap model/product score differs')
        rows[key]=q
    if len(rows)!=10482:raise ArithmeticError('fixed retained union differs')
    return [rows[k] for k in sorted(rows)]
def local(A,B,q):
    a,b=A,B;k=0
    while a%q**4==0 and b%q**6==0:a//=q**4;b//=q**6;k+=1
    a%=q;b%=q
    if (4*a*a*a+27*b*b)%q:raise ArithmeticError('unscaling reveals good reduction: missing trace must be handled explicitly')
    if not a:return {'prime':q,'scaling_exponent':k,'A_mod_p':a,'B_mod_p':b,'kind':'CUSPIDAL_DISPLAY'}
    r=-3*b*pow(2*a,-1,q)%q
    if (a+3*r*r)%q or (b-2*r**3)%q:raise ArithmeticError('double-root identity differs')
    test=pow(3*r%q,(q-1)//2,q)
    if test not in (1,q-1):raise ArithmeticError('nonzero tangent square test required')
    return {'prime':q,'scaling_exponent':k,'A_mod_p':a,'B_mod_p':b,'double_root':r,'tangent_square':3*r%q,'euler_criterion':test,'kind':'SPLIT_NODE' if test==1 else 'NONSPLIT_NODE'}
def run(check=False):
    p=protocol();out=D/'result.json';primes=[q for q in _primes_up_to(p['prime_bound']) if q>=5];bonus={q:round(log(7*(q-1)/(5*q))*10**12) for q in primes};excluded={(r['family'],r['parameter']) for path in EXCLUDED for r in cert.read(path)['rows']}
    if not check and out.exists():raise FileExistsError('preserve split-node scores')
    d={'schema':'elliptic-curves.higher-split-node-score-result.v1','status':'RUNNING','protocol_hash':digest(p),'rows':[]}
    for row in pool():
        A,B=map(int,row['model'][3:]);delta=-16*(4*A**3+27*B**2)
        if not delta:raise ArithmeticError('singular rational curve')
        places=[local(A,B,q) for q in primes if delta%q==0];low=sum(bonus[r['prime']] for r in places if r['kind']=='SPLIT_NODE' and r['prime']<=32749);high=sum(bonus[r['prime']] for r in places if r['kind']=='SPLIT_NODE' and r['prime']>32749)
        d['rows'].append({**row,'bad_places':places,'selection_bonus_units':low,'validation_bonus_units':high,'augmented_selection_units':row['product_selection_units']+low,'augmented_validation_units':row['product_validation_units']+high,'excluded_prior_address':(row['family'],row['parameter']) in excluded})
        if len(d['rows'])%128==0 and not check:checkpoint(out,d)
        if len(d['rows'])%1024==0:print('SPLIT-NODE SCORES',len(d['rows']),'of10482',flush=True)
    d['prospective_candidates']=[];d['selection']={}
    for family in sorted({r['family'] for r in d['rows']}):
        available=sorted((r for r in d['rows'] if r['family']==family and not r['excluded_prior_address']),key=lambda r:(-r['augmented_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))[:4]
        if len(available)!=4:raise ArithmeticError('fixed four-per-family selection unavailable')
        d['selection'][family]=[r['parameter'] for r in available];d['prospective_candidates'] += available
    d['status']='COMPLETE_FIXED_SPLIT_NODE_SCORING'
    if check:
        if cert.read(out)!=d:raise ArithmeticError('exact split-node score replay differs')
    else:checkpoint(out,d)
    print('SPLIT-NODE SCORING10482 AND FIXED24 REPLAY PASS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();prepare() if a.stage=='prepare' else run(a.stage=='replay')

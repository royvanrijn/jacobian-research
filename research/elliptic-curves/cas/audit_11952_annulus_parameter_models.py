#!/usr/bin/env python3
"""Independent homogeneous specialization and prospective equation-overlap audit."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
P=ROOT/'artifacts/local/elliptic-curves/11952-new-annulus-scores-v2/bank-protocol.json'
SCAN=ROOT/'artifacts/local/elliptic-curves/11952-new-annulus-v1/result.json'
OUT=ART/'annulus_11952_parameter_models_v1.json'
BANDS={1:(131072,524288)}

def stats(values):
    a=sorted(values);n=len(a)
    return {'minimum':a[0],'median':str(F(a[(n-1)//2]+a[n//2],2)),'maximum':a[-1]}

def expected():
    p=cert.read(P);scan=cert.read(SCAN)
    if len(p['rows'])!=32768 or scan['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION':raise ArithmeticError('fixed retained population required')
    families={r['family']:r for r in cert.read(spec.ATLAS)['families']}
    source={(s['id'],r['parameter']):r for s in scan['shards'] for r in s['rows']}
    if len(source)!=32768:raise ArithmeticError('duplicate retained address')
    seen={};pairs=[];groups={};keys=set();jcounts=Counter()
    for row in p['rows']:
        n=row['numerator'];d=row['denominator'];t=F(n,d);b=1;lo,hi=BANDS[b]
        if t.numerator!=n or t.denominator!=d or str(t)!=row['parameter'] or not lo<max(abs(n),d)<=hi:raise ArithmeticError('nonprimitive or out-of-band parameter')
        key=(row['slice_id'],row['parameter'])
        if key in keys or key not in source:raise ArithmeticError('score roster differs from retained addresses')
        keys.add(key);old=source[key]
        if any(row[k]!=old[k] for k in ('numerator','denominator','parameter','score_units','good_primes')):raise ArithmeticError('retained address/score changed')
        f=families[row['family']];model=[F(0)]*3
        for name,weight in [('A_coefficients_low_to_high',8),('B_coefficients_low_to_high',12)]:
            coeffs=list(map(F,f[name]))
            if len(coeffs)>weight+1:raise ArithmeticError('unexpected family degree')
            model.append(sum(c*n**k*d**(weight-k) for k,c in enumerate(coeffs)))
        if any(a.denominator!=1 for a in model) or list(map(str,model))!=row['model']:raise ArithmeticError('independent homogeneous specialization differs')
        a,bcoef=model[3:];delta=-16*(4*a**3+27*bcoef**2)
        if not delta:raise ArithmeticError('singular retained model')
        j=(-48*a)**3/delta;equation_matches=[]
        for other in seen.get(j,[]):
            if cert.isomorphic(model,other['model']):equation_matches.append(other['id'])
        if equation_matches:pairs.append({'id':row['id'],'matches':equation_matches})
        seen.setdefault(j,[]).append({'id':row['id'],'model':model});jcounts[b]+=int(not equation_matches)
        group=groups.setdefault((b,row['family']),{'height':[],'bits':[],'positive':0,'negative':0})
        group['height'].append(max(abs(n),d));group['bits'].append(max(abs(int(a)).bit_length() for a in model));group['positive' if n>0 else 'negative']+=1
    if keys!=set(source):raise ArithmeticError('missing retained equation')
    records=[{'band':b,'family':f,'curves':len(v['height']),'parameter_height':stats(v['height']),'largest_integral_coefficient_bits':stats(v['bits']),'positive':v['positive'],'negative':v['negative']} for (b,f),v in sorted(groups.items())]
    if len(records)!=1 or any(r['curves']!=32768 or r['positive']!=16384 or r['negative']!=16384 for r in records):raise ArithmeticError('balanced signed band/family counts differ')
    paths=[Path(__file__).resolve(),Path(cert.__file__),Path(spec.__file__),spec.ATLAS,P,SCAN,ROOT/'elliptic-curves/cas/research_runtime/store.py']
    return {'schema':'elliptic-curves.11952-annulus-parameter-models.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'curves_checked':32768,'groups':records,'distinct_j_invariants':len(seen),'distinct_Q_isomorphism_classes':32768-len(pairs),'within_roster_matches':pairs,'boundaries':'Independent exact homogeneous specialization and primitive annulus membership for the32768 frozen prospective models. Exact within-roster Q-isomorphism comparisons only; no catalogue or old point results read. Integral coefficient sizes are for these displayed models, not globally minimal heights. Different parameters are not automatically different curves; distinct classes here need not be absent from previous work. No new point, rank bound, selector modification or isolated portability claim.'}
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('outer model audit replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve outer model audit')
        checkpoint(OUT,d)
    print('NEW11952 ANNULUS MODEL AUDIT',d['curves_checked'],d['distinct_j_invariants'],d['distinct_Q_isomorphism_classes'],'PASS',flush=True)

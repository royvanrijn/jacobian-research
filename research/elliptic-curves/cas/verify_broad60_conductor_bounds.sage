#!/usr/bin/env sage-python
"""Independent PARI Tate replay of thethree-curve Sage generic local audit."""
import argparse,json,sys
from pathlib import Path
from sage.all import QQ,ZZ,pari,prime_range
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'broad60_conductor_bounds_v1.json'
OUT=ART/'broad60_conductor_bounds_pari_replay_v1.json'

def kodaira(k):
    if k in (1,2,3,4,-1,-2,-3,-4):return {1:'I0',2:'II',3:'III',4:'IV',-1:'I0*',-2:'II*',-3:'III*',-4:'IV*'}[k]
    return 'I'+str(abs(k)-4)+('*' if k<0 else '')

def expected():
    d=cert.read(INPUT)
    if d['status']!='PASS_COMPLETE_FIXED_LOCAL_AUDIT':raise ArithmeticError('terminal audit required')
    for name,h in d['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('source binding changed')
    inventory=[r for r in cert.read(ART/'broad60_mw16_results_v1.json')['curves'] if r['rank_lower_bound']>=22 and not r['icarm_matches'] and not r['previous_matches']]
    if len(inventory)!=3 or [r['id'] for r in inventory]!=[r['id'] for r in d['rows']]:raise ArithmeticError('whole inventory required')
    counts=[]
    for source,row in zip(inventory,d['rows']):
        if any(row[k]!=source[k] for k in ('id','curve','family','parameter','rank_lower_bound')):raise ArithmeticError('inventory row differs')
        E=pari.ellinit([QQ(a) for a in row['integral_curve']]);original=pari.ellinit([QQ(a) for a in source['curve']])
        # PARI's explicit c4,c6,Delta fields give an independent invariant check.
        if E[9]!=original[9] or E[10]!=original[10] or E[11]!=original[11]:raise ArithmeticError('integral invariant transport differs')
        delta=ZZ(E[11]);remaining=abs(delta);upper=ZZ(1);expected_primes=[]
        for p in prime_range(2,10001):
            if delta%p==0 or p in (2,3):expected_primes.append(int(p))
        if expected_primes!=[r['prime'] for r in row['local_data']]:raise ArithmeticError('local coverage differs')
        for local in row['local_data']:
            p=ZZ(local['prime']);raw=E.elllocalred(p);f=int(raw[0]);u=QQ(raw[2][0])
            e=int(delta.valuation(p));emin=e-12*int(u.valuation(p))
            if f!=local['conductor_exponent'] or kodaira(int(raw[1]))!=local['kodaira'] or e!=local['displayed_discriminant_valuation'] or emin!=local['minimal_discriminant_valuation']:
                raise ArithmeticError('independent PARI local Tate result differs')
            remaining//=p**e;upper*=p**f
        upper*=remaining
        if str(delta)!=row['discriminant'] or str(remaining)!=row['remaining_cofactor'] or str(upper)!=row['conductor_upper_bound']:
            raise ArithmeticError('independent conductor product differs')
        counts.append({'id':row['id'],'local_pairs':len(row['local_data']),'conductor_upper_bound':str(upper)})
    return {'schema':'elliptic-curves.broad60-conductor-pari-replay.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,ART/'broad60_mw16_results_v1.json')},
        'pari_version':str(pari.version()),'curves':3,'local_pairs':sum(r['local_pairs'] for r in counts),'rows':counts,
        'claim_boundary':'Independent PARI elllocalred agrees with every one of the Sage generic Tate local results and reconstructs all three new-curve conductor upper bounds. This does not factor the residual cofactors, verify external catalogue conductors, reprove ranks or assert conductor records.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('PARI replay record differs')
    else:
        if OUT.exists():raise FileExistsError('preserve independent conductor replay')
        checkpoint(OUT,d)
    print('PARI CONDUCTOR REPLAY',d['curves'],d['local_pairs'],'PASS',flush=True)

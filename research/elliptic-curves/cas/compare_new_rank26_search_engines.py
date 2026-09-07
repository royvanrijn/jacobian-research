#!/usr/bin/env python3
"""Compare exact affine square-hit sets on43 frozen PARI/GMP control boxes."""
import argparse,subprocess,time
from pathlib import Path
from math import isqrt
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from pointed_quartic_search import PointedQuarticSearch
ROOT=Path(__file__).resolve().parents[2];LOCAL=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=LOCAL/'new-rank26-engine-comparison-v1';INPUT=LOCAL/'new-rank26-pari43-v1/candidate-00/result.json';MAPS=LOCAL/'new-rank26-pari-maps-v1/maps.json'
OUTPUT=ART/'new_rank26_engine_comparison_v1.json'

def run():
    protocol=cert.read(D/'protocol.json');source=cert.read(INPUT);maps=cert.read(MAPS)
    if cert.hashed(INPUT)!=protocol['input_sha256'] or cert.hashed(MAPS)!=protocol['maps_sha256'] or cert.hashed(Path(__file__).resolve())!=protocol['source_sha256']:raise ArithmeticError('frozen engine comparison differs')
    if OUTPUT.exists():raise FileExistsError('preserve comparison attempt')
    result={'schema':'elliptic-curves.new-rank26-fixed-box-engine-comparison.v1','status':'RUNNING','protocol_hash':digest(protocol),'rows':[],'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,MAPS)}};checkpoint(OUTPUT,result)
    for i,(row,modelmap) in enumerate(zip(source['charts'],maps['rows'])):
        r=row['search']
        if r['completed_denominator']!=100000 or r['height_bound']!=100000:raise ArithmeticError('GMP box incomplete')
        P,Q=(modelmap[k] for k in ('reduced_P','reduced_Q'));poly=lambda a:'+'.join(f'({v})*x^{j}' for j,v in enumerate(a))
        program='C=['+poly(P)+','+poly(Q)+'];gettime();R=hyperellratpoints(C,100000);print("MS|",gettime());for(i=1,#R,print("X|",R[i][1]));quit\n'
        started=time.monotonic()
        try:
            completed=subprocess.run(['/usr/bin/gp','-q','-s','256000000'],input=program,text=True,capture_output=True,timeout=3)
            checkpoint(D/f'chart-{i:02}.json',{'program':program,'stdout':completed.stdout,'stderr':completed.stderr,'returncode':completed.returncode})
        except subprocess.TimeoutExpired as e:
            checkpoint(D/f'chart-{i:02}.json',{'program':program,'timeout_seconds':3,'stdout':str(e.stdout),'stderr':str(e.stderr)})
            result['rows'].append({'chart':i,'status':'CENSORED'});checkpoint(OUTPUT,result);continue
        if completed.returncode or '***' in completed.stderr:raise ArithmeticError('PARI fixed-box search failed')
        xs={cert.F(line[2:]) for line in completed.stdout.splitlines() if line.startswith('X|')};hits={cert.F(int(a),int(b)) for a,b,_ in r['primitive_square_hits']}
        # Validate every returned rational coordinate on the exact final quartic.
        coefficients=tuple(map(int,r['coefficients']));checked=[]
        for x in sorted(xs):
            n,d=x.numerator,x.denominator
            if max(abs(n),d)>100000:raise ArithmeticError('PARI coordinate outside fixed box')
            f=sum(coefficients[j]*n**j*d**(4-j) for j in range(5))
            if f<0 or isqrt(f)**2!=f:raise ArithmeticError('PARI coordinate not an exact final quartic square')
            checked.append([str(n),str(d),str(isqrt(f))])
        result['rows'].append({'chart':i,'status':'PASS' if xs==hits else 'EXACT_HIT_SET_MISMATCH','pari_wall_seconds':time.monotonic()-started,'pari_search_ms':int(next(line[3:] for line in completed.stdout.splitlines() if line.startswith('MS|'))),'affine_square_coordinates':checked,'gmp_hit_count':len(hits),'pari_hit_count':len(xs),'pari_only':list(map(str,sorted(xs-hits))),'gmp_only':list(map(str,sorted(hits-xs)))})
        checkpoint(OUTPUT,result);print('FIXED BOX ENGINE COMPARISON',i+1,result['rows'][-1]['status'],flush=True)
    result['status']='COMPLETE_DECLARED_COMPARISON';result['claim_boundary']='New-curve finite boxes only. Exact affine square-coordinate equality where PASS; infinity remains separately checked by the shared engine. A speed comparison does not establish universal completeness, future speedups, rank bounds or new curves.';checkpoint(OUTPUT,result)

def check():
    result=cert.read(OUTPUT);source=cert.read(INPUT)
    if any(cert.hashed(ROOT/p)!=h for p,h in result['sources'].items()):raise ArithmeticError('comparison source changed')
    for row in result['rows']:
        if row['status']!='PASS':raise ArithmeticError('comparison is not uniformly PASS')
        r=source['charts'][row['chart']]['search'];saved=cert.read(D/f"chart-{row['chart']:02}.json");xs={cert.F(line[2:]) for line in saved['stdout'].splitlines() if line.startswith('X|')};coefficients=tuple(map(int,r['coefficients']))
        expected=[]
        for x in sorted(xs):
            n,d=x.numerator,x.denominator;f=sum(coefficients[j]*n**j*d**(4-j) for j in range(5));root=isqrt(f)
            if root*root!=f or max(abs(n),d)>100000:raise ArithmeticError('square witness differs')
            expected.append([str(n),str(d),str(root)])
        if expected!=row['affine_square_coordinates'] or xs!={cert.F(int(n),int(d)) for n,d,_ in r['primitive_square_hits']}:raise ArithmeticError('hit set differs')
    if len(result['rows'])!=43 or result['status']!='COMPLETE_DECLARED_COMPARISON':raise ArithmeticError('incomplete comparison')
    print('REPLAYED IDENTICAL AFFINE HIT SETS',len(result['rows']),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else run()

#!/usr/bin/env python3
"""Compare saved generic17 and fixed own27 charts on every public extra witness certified mod2."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
from audit_recorded_point_mod2_rank_v3 import insert
from memory_rank_certificate import checked_rank
from search_observability import multiply
from research_runtime.store import checkpoint,digest

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PUBLIC=ART/'inventory188_public28_reproduction_v1.json'
GEOMETRY=ART/'inventory188_own27_geometry_visibility_v1.json'
SOURCE=ART/'full11952_late64_r17_results_v1.json'
MAPS=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1/maps.json'
OUT=ART/'inventory188_two_witness_chart_comparison_v1.json'


def expected():
    public=cert.read(PUBLIC);geometry=cert.read(GEOMETRY)
    if public['status']!='PASS' or geometry['status']!='PASS_EXACT_GEOMETRY_AND_COORDINATES':raise ArithmeticError('completed point and map proofs required')
    for data in (public,geometry):
        if any(cert.hashed(ROOT/n)!=h for n,h in data['sources'].items()):raise ArithmeticError('bound proof changed')
    signatures=[r for s in public['union_signatures'] for r in s['rows']];indices=[]
    for i in range(28):
        pivots={}
        for row in signatures:insert(pivots,row[:27]+[row[27+i]])
        if len(pivots)==28:indices.append(i)
    if indices!=[26,27]:raise ArithmeticError('fixed complete mod2-certified extra-witness roster differs')
    model=tuple(map(cert.F,public['curve']));old=[tuple(map(cert.F,P)) for P in public['points'][:27]]
    primes=[s['prime'] for s in public['union_signatures']];tp=public['rank_certificate']['no_rational_2_torsion_prime']
    proofs={}
    for i in indices:
        P=tuple(map(cert.F,public['transported_public_points'][i]));proofs[str(i)]=checked_rank(model,old+[P],primes,tp)
    r=next(r for r in cert.read(SOURCE)['curves'] if r['parameter']==public['parameter']);rawpath=ROOT/r['discovery_witness']['path']
    if cert.hashed(rawpath)!=r['discovery_witness']['sha256']:raise ArithmeticError('old charts changed')
    arms={'generic17':[],'own27':[]}
    for c in cert.read(rawpath)['charts']:
        s=c['search'];q=s['pointed_chart'];den=cert.F(q['point_denominator_root']);scale=cert.F(q['curve_coordinate_scale']);k=cert.F(q['shift_mod_denominator_squared'])
        M=multiply((den/scale,k/(den*scale),0,1),multiply(tuple(map(cert.F,q['unimodular_horizontal_matrix'])),tuple(map(cert.F,s['horizontal_matrix']))))
        arms['generic17'].append(((cert.F(s['base_point']['x']),cert.F(s['base_point']['y'])),M))
    for r in cert.read(MAPS)['rows']:
        raw=list(map(cert.F,r['raw_coefficients']));arms['own27'].append(((-raw[2]/6,-raw[1]/8),tuple(map(cert.F,r['matrix']))))
    if any(len(rows)!=49 for rows in arms.values()):raise ArithmeticError('equal49-chart arms required')
    observations=[];minima={}
    for arm,charts in arms.items():
        for i in indices:
            P=tuple(map(cert.F,public['transported_public_points'][i]));local=[]
            for chart,((x,y),(a,b,c,d)) in enumerate(charts):
                for sign in (-1,1):
                    if P[0]==x:raise ArithmeticError('extra point cannot equal old-subgroup endpoint')
                    t=(sign*P[1]+y)/(P[0]-x);num=d*t-b;den=a-c*t
                    if not den:coord=['1','0'];height=1;infinity=True
                    else:
                        s=num/den;coord=[str(s.numerator),str(s.denominator)];height=max(abs(s.numerator),s.denominator);infinity=False
                    row={'arm':arm,'public_point_index':i,'chart':chart,'sign':sign,'coordinate':coord,'height':height,'at_infinity':infinity,'inside_height125000_or_infinity':infinity or height<=125000}
                    observations.append(row);local.append(row)
            minima[arm+':'+str(i)]=min(local,key=lambda r:(r['height'],r['chart'],r['sign']))
    paths=[Path(__file__).resolve(),Path(cert.__file__),ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/search_observability.py',ROOT/'elliptic-curves/cas/audit_recorded_point_mod2_rank_v3.py',PUBLIC,GEOMETRY,SOURCE,MAPS,rawpath]
    return {'schema':'elliptic-curves.inventory188-two-witness-chart-comparison.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'witness_rank_certificates':proofs,
            'observations':observations,'minima':minima,
            'visible_observations_by_arm':{a:sum(r['inside_height125000_or_infinity'] for r in observations if r['arm']==a) for a in arms},
            'claim_boundary':'Exact coordinates for both public representatives independently certified to extend the old27 subgroup by the retained mod2 signatures, both signs, on equal49 generic17 and own27 maps. This is a retrospective two-witness, one-curve diagnostic. The other26 public points are not classified as lying in the old rational span. No point search, prospective sensitivity, upper bound, untested-translate exclusion or general chart-policy theorem.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if digest(cert.read(OUT))!=digest(d):raise ArithmeticError('two-witness chart replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve matched chart diagnostic')
        checkpoint(OUT,d)
    print('TWO WITNESS CHART COMPARISON',d['visible_observations_by_arm'],{k:v['height'] for k,v in d['minima'].items()})

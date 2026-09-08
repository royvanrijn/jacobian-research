#!/usr/bin/env python3
"""Integral equations and exact transported points for the five new26 candidates."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as arithmetic
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'compact192_r17_results_v1.json'
OUT=ART/'compact192_rank26_models_v1.json'
SAGE=ART/'new_compact192_rank26_curves.sage'
IDS=('07ca9-001','07ca9-006','07ca9-011','103b2-029','11952-010')


def expected():
    source=cert.read(INPUT)
    selected=[r for r in source['curves'] if r['rank_lower_bound']>=26 and not r['icarm_matches'] and not r['previous_matches']]
    if tuple(r['id'] for r in selected)!=IDS or any(r['rank_lower_bound']!=26 for r in selected):raise ArithmeticError('fixed five newly unmatched26 roster differs')
    stages=LOCAL=ROOT/'artifacts/local/elliptic-curves/compact192-r17-pari-v1/post-batch'
    paths=[Path(__file__).resolve(),INPUT]
    for name in ('certify','proof-replay'):
        path=stages/(name+'.supervisor.json');s=cert.read(path);paths.append(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('complete post-terminal point proofs and catalogue replay required')
    rows=[]
    for r in selected:
        short=tuple(map(cert.F,r['curve']));model=arithmetic.integral(short);inv=cert.weierstrass_invariants(model);points=[]
        for P in r['points']:
            x,y=map(cert.F,P);X=x-inv['b2']/12;points.append((X,y-(model[0]*X+model[2])/2))
        back=[(x+inv['b2']/12,y+(model[0]*x+model[2])/2) for x,y in points]
        if short!=(0,0,0,-inv['c4']/48,-inv['c6']/864) or [list(map(str,P)) for P in back]!=r['points'] or any(not cert.is_on_weierstrass_curve(model,P) for P in points):raise ArithmeticError('exact integral point transport differs')
        proof=r['rank_certificate'];actual=checked_rank(short,back,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.loads(json.dumps(actual))!=proof or len(points)!=26:raise ArithmeticError('independent26-point proof differs')
        minimality=arithmetic.minimality(model)
        if any(cert.isomorphic(model,q['ainvs']) for q in source['catalogue']['equations']) or any(cert.isomorphic(model,q['curve']) for q in source['previous_equations']):raise ArithmeticError('catalogue or previous model match after transport')
        rows.append({'id':r['id'],'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':26,'minimal_curve':list(map(str,model)),'points':[list(map(str,P)) for P in points],'discovery_curve':r['curve'],'discovery_points':r['points'],'rank_certificate_on_discovery_curve':proof,'minimality':minimality})
    if any(cert.isomorphic(r['minimal_curve'],s['minimal_curve']) for i,r in enumerate(rows) for s in rows[:i]):raise ArithmeticError('five rational-isomorphism-distinct models required')
    return {'schema':'elliptic-curves.compact192-rank26-models.v1','status':'PASS','sources':{**arithmetic.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}},'curves':rows,'catalogue':source['catalogue'],'claim_boundary':'Five globally minimal integral equations and26 exactly independent rational points on each, with exact transports to their discovery point certificates. All five are mutually nonisomorphic overQ, absent from the pinned593-equation catalogue and549 previous measured equations. This neither proves exact rank, conductor, full saturation, a record nor universal novelty. Full batch inventory and isolated replay are separate.'}


def export(data):
    lines=['# Five distinct curves with certified rank at least26.',
           '# Exact independence, model minimality and pinned catalogue comparisons:',
           '# compact192_rank26_models_v1.json. No exact-rank or universal-novelty claim.',
           'from sage.all import QQ, EllipticCurve','curves = {}']
    for row in data['curves']:
        lines += [f"E = EllipticCurve(QQ, {row['minimal_curve']!r})",
                  f"points = [E([QQ(x), QQ(y)]) for x, y in {row['points']!r}]",
                  'assert len(points) == 26',
                  f"curves[{row['id']!r}] = (E, points)",
                  f"print({row['id']!r}, E, 'rank >=26')"]
    return '\n'.join(lines)+'\n'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();data=expected();script=export(data)
    if args.check:
        if cert.read(OUT)!=json.loads(json.dumps(data)) or SAGE.read_text()!=script:raise ArithmeticError('five model certificate or Sage export differs')
    else:
        if OUT.exists() or SAGE.exists():raise FileExistsError('preserve five26 model export')
        checkpoint(OUT,data);SAGE.write_text(script)
    print('FIVE EXACT MINIMAL26 MODELS PASS',[(r['id'],r['minimality']['invariant_gcd']) for r in data['curves']],flush=True)

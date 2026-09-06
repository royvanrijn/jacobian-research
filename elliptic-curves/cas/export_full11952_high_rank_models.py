#!/usr/bin/env python3
"""Integral equations and independent points for unmatched full11952 bounds>=27."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as arithmetic
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'full11952_64_r17_results_v1.json';OUT=ART/'full11952_high_rank_models_v1.json';SAGE=ART/'new_full11952_high_rank_curves.sage'

def expected():
    source=cert.read(INPUT);portable=cert.read(ART/'full11952_64_point_portable_replay_v1.json')
    if portable['status']!='PASS' or portable['logical_stages']!=194:raise ArithmeticError('all194 isolated point-proof stages required')
    selected=[r for r in source['curves'] if r['rank_lower_bound']>=27 and not r['icarm_matches'] and not r['previous_matches']]
    paths=[Path(__file__).resolve(),INPUT,ART/'full11952_64_point_portable_replay_v1.json'];rows=[]
    for r in selected:
        short=tuple(map(cert.F,r['curve']));model=arithmetic.integral(short);inv=cert.weierstrass_invariants(model);points=[]
        for P in r['points']:
            x,y=map(cert.F,P);X=x-inv['b2']/12;points.append((X,y-(model[0]*X+model[2])/2))
        back=[(x+inv['b2']/12,y+(model[0]*x+model[2])/2) for x,y in points]
        if short!=(0,0,0,-inv['c4']/48,-inv['c6']/864) or [list(map(str,P)) for P in back]!=r['points'] or any(not cert.is_on_weierstrass_curve(model,P) for P in points):raise ArithmeticError('exact integral transport differs')
        proof=r['rank_certificate'];actual=checked_rank(short,back,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.loads(json.dumps(actual))!=proof or len(points)!=r['rank_lower_bound']:raise ArithmeticError('independent point proof differs')
        # Minimality is an additional cheap certificate, never a prerequisite for the rank proof.
        try:minimality={'status':'PROVED_GLOBAL_MINIMAL','certificate':arithmetic.minimality(model)}
        except ArithmeticError as exc:minimality={'status':'UNKNOWN','reason':str(exc)}
        if any(cert.isomorphic(model,q['ainvs']) for q in source['catalogue']['equations']) or any(cert.isomorphic(model,q['curve']) for q in source['previous_equations']):raise ArithmeticError('catalogue/prior match after exact transport')
        rows.append({'id':r['id'],'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':r['rank_lower_bound'],'integral_curve':list(map(str,model)),'points':[list(map(str,P)) for P in points],'discovery_curve':r['curve'],'discovery_points':r['points'],'rank_certificate_on_discovery_curve':proof,'minimality':minimality})
    if any(cert.isomorphic(r['integral_curve'],s['integral_curve']) for i,r in enumerate(rows) for s in rows[:i]):raise ArithmeticError('mutual rational nonisomorphism required')
    return {'schema':'elliptic-curves.full11952-high-rank-models.v1','status':'PASS','sources':{**arithmetic.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}},'curves':rows,'catalogue':source['catalogue'],'claim_boundary':'All cohort equations with certified lower bound at least27 and no match among the pinned593 catalogue and789 prior equations, with integral equations and exactly transported independent points. A minimality claim is made only when its separate certificate says PROVED_GLOBAL_MINIMAL. Exact ranks, conductors, full saturation, records and universal novelty remain unproved.'}

def export(data):
    lines=['# Full11952 curves with certified lower bounds at least27.',
           '# Exact proofs and pinned comparisons: full11952_high_rank_models_v1.json.',
           '# No exact-rank, conductor, record or universal-novelty claim.',
           'from sage.all import QQ, EllipticCurve','curves = {}']
    for row in data['curves']:
        lines += [f"E = EllipticCurve(QQ, {row['integral_curve']!r})",f"points = [E([QQ(x), QQ(y)]) for x, y in {row['points']!r}]",f"assert len(points) == {row['rank_lower_bound']}",f"curves[{row['id']!r}] = (E, points)",f"print({row['id']!r}, E, 'rank >=', {row['rank_lower_bound']})"]
    return '\n'.join(lines)+'\n'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected();script=export(d)
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)) or SAGE.read_text()!=script:raise ArithmeticError('model certificate or Sage export differs')
    else:
        if OUT.exists() or SAGE.exists():raise FileExistsError('preserve high-rank model export')
        checkpoint(OUT,d);SAGE.write_text(script)
    print('FULL11952 HIGH-RANK MODELS',[(r['id'],r['rank_lower_bound'],r['minimality']['status']) for r in d['curves']],flush=True)

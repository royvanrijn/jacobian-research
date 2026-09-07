#!/usr/bin/env python3
"""Independent exact replay of the marked auxiliary rank-three point subgroup."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict
import certify_compact_r17_candidates as cert
import audit_retained_cloud_modl as odd
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'native_rank3_carrier_marked_point_v1.json';CARRIER=ART/'rank_jump_minimal_native_block_carrier_v1.json';OUT=ART/'native_rank3_carrier_marked_point_replay_v1.json'
def short(model,P=None):
    a1,a2,a3,a4,a6=map(F,model);b=a2+a1*a1/4;c=a4+a1*a3/2;d=a6+a3*a3/4;curve=(F(0),F(0),F(0),c-b*b/3,d-b*c/3+2*b*b*b/27)
    if P is None:return curve
    x,y=map(F,P);return (x+b/3,y+(a1*x+a3)/2)
def expected():
    d=cert.read(INPUT);src=cert.read(CARRIER);g=src['geometry']
    if d['status']!='PASS' or any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('marked proof source changed')
    e,b1,c,b3,a=map(F,reversed(g['quartic_coefficients']));q=cert.square_root(e)
    if q is None or not q:raise ArithmeticError('rational marked quartic point required')
    base=(b1/q,c-b1*b1/(4*q*q),2*q*b3,-4*q*q*a,a*(b1*b1-4*q*q*c));x=b1*b1/(4*q*q)-c;opps=[(x,F(0)),(x,-base[0]*x-base[2])]
    if list(map(str,base))!=d['pointed_curve'] or str(x)!=d['opposite_x'] or set(opps)!={tuple(map(F,P)) for P in d['opposite_points']} or any(not cert.is_on_weierstrass_curve(base,P) for P in opps):raise ArithmeticError('rational opposite-point construction differs')
    J=tuple(map(F,g['minimal_Jacobian_model']));points=[tuple(map(F,P)) for P in d['Jacobian_points']]
    if list(map(str,J))!=d['Jacobian_model'] or d['Jacobian_points'][:2]!=src['descent']['points'] or any(not cert.is_on_weierstrass_curve(J,P) for P in points):raise ArithmeticError('auxiliary point identity differs')
    source=short(base);target=short(J);v,w=cert.weierstrass_invariants(source),cert.weierstrass_invariants(target);u=cert.square_root(v['c6']*w['c4']/(v['c4']*w['c6']))
    if u is None or not u or v['c4']!=u**4*w['c4'] or v['c6']!=u**6*w['c6']:raise ArithmeticError('rational Jacobian scale differs')
    transported={short(J,P) for P in points[2:]};old=[short(base,P) for P in opps]
    if not any(transported=={(x/u**2,y/(sign*u)**3) for x,y in old} for sign in (1,-1)):raise ArithmeticError('opposite point transport differs')
    shortpoints=[short(J,P) for P in points]
    if list(map(str,target))!=d['short_curve'] or [[str(v) for v in P] for P in shortpoints]!=d['short_points']:raise ArithmeticError('short auxiliary model differs')
    rows=[]
    for audit in d['audits']:
        ell=audit['modulus'];signatures=[]
        if ell not in (3,5) or not odd.ml.no_rational_l_torsion_reduction_certificate(target,audit['no_rational_ell_torsion_prime'],ell):raise ArithmeticError('odd torsion certificate differs')
        for old in audit['signatures']:
            actual=odd.finite.signature(target,shortpoints,old['prime'],ell)
            if json.dumps(asdict(actual),sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('exact finite signature differs')
            signatures.extend(actual.rows)
        pivots=odd.finite.pivots(signatures,ell)
        if pivots!=audit['independent_column_indices'] or len(pivots)!=audit['finite_column_rank'] or len(pivots)!=3:raise ArithmeticError('three independent auxiliary columns required')
        rows.append({'modulus':ell,'rank_lower_bound':len(pivots),'independent_indices':pivots})
    paths=[Path(__file__).resolve(),INPUT,CARRIER,Path(cert.__file__)]
    return {'schema':'elliptic-curves.native-carrier-marked-point-replay.v1','status':'PASS','sources':{**odd.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}},'audits':rows,'rank_lower_bound':3,'three_independent_Jacobian_points':[d['Jacobian_points'][i] for i in rows[0]['independent_indices']],'claim_boundary':'Independent rational formulas reconstruct both marked quartic points and their Jacobian transports; exact finite reductions modulo3 and5 prove three independent auxiliary points despite rational2-torsion. This closes the explicit rank-three subspace gap, not saturation or a full integral generator claim. No new original fibre or original-rank increase is asserted by this auxiliary proof.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('marked-point independent replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve independent marked-point proof')
        checkpoint(OUT,d)
    print('INDEPENDENT EXPLICIT AUXILIARY RANK3 POINT PROOF PASS',flush=True)

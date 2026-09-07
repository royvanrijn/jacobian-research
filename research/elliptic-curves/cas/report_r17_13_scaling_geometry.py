#!/usr/bin/env python3
"""Bind exact integral parameter charts and the universal13-scale classification."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_r17_integral13_charts as charts
import classify_r17_13_scaling as classification
import verify_r17_13_scaling_classification as independent
from research_runtime.store import checkpoint
ROOT=charts.ROOT;ART=charts.ART;OUT=ART/'r17_13_scaling_geometry_v1.json'

def expected():
    c=charts.expected();d=classification.expected()
    if cert.read(charts.OUT)!=json.loads(json.dumps(c)) or cert.read(classification.OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('exact chart/classification input differs')
    independent.main();rows=[]
    for r in d['rows']:
        q=next(q for q in c['rows'] if q['family']==r['family'] and q['integral_after_curve_scale13'])
        rows.append({'family':r['family'],'unique_nonminimal_residue_mod13':r['universal_residue_mod13'],'parameter_matrix':q['matrix'],'maximum_removable_13_exponent':1,'scaled_models_still_bad_at13':r['all_eligible_scaled_models_remain_bad_at13'],'square_box_weighted_bound_ratio':q['weighted_coefficient_bound_ratio'],'signed_permutation_self_presentations':q['signed_permutation_self_presentations']})
    paths=[Path(__file__).resolve(),Path(charts.__file__).resolve(),Path(classification.__file__).resolve(),Path(independent.__file__).resolve(),ROOT/'elliptic-curves/cas/verify_r17_integral13_charts.sage',charts.INPUT,charts.OUT,classification.OUT,classification.D/'protocol.json']
    for name in ('symbolic','independent-residues'):
        path=classification.D/(name+'.supervisor.json');s=cert.read(path);paths.append(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('independent exact check required')
    return {'schema':'elliptic-curves.r17-13-scaling-geometry.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'projective_cells_checked':84,'universal_integral_charts':6,'maximum_modular_witness_depth':3,'claim_boundary':'For each displayed compact R17 binary pair and every primitive integer parameter pair,13^4|A_h and13^6|B_h hold exactly in the single stated projective residue cell. One and only one13 scaling is removable there; no primitive pair admits13^8|A_h and13^12|B_h. For nonsingular fibres these are exact13-adic minimality statements, and the reduced model remains bad at13. The determinant13 parameter changes and explicit x/13^2,y/13^3 homogeneous scaling preserve the Q(t) family and its generic17-point subgroup. This proves neither a new family, improved search population, global model minimality, conductor, whole-curve rank nor a new near-record curve.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('integral13 geometry aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve integral13 geometry aggregate')
        checkpoint(OUT,d)
    print('UNIVERSAL13-SCALING CLASSIFICATION AND SIX EXACT PARAMETER CHARTS PASS',flush=True)

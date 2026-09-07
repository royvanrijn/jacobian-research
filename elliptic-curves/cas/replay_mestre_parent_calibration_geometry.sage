#!/usr/bin/env sage-python
"""Exact parity/norm transports and rational maps of the frozen specialized sample."""
import sys
from decimal import Decimal,localcontext
from pathlib import Path
from sage.all import matrix,ZZ,QQ,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import mestre_parent_calibration as control
import certify_compact_r17_candidates as cert
from replay_retention24_geometry import geometry,cloud_check,tuples
from research_runtime.store import digest
from mestre_parent_adapter import specialize

def main():
    p=control.protocol();rank=control.ROW['initial_rank'];D=control.D;maps=cert.read(D/'maps.json');data=cert.read(D/'result.json');seed=cert.read(control.SEED);initial=tuples(seed['points'])
    model,images=specialize(control.ROW['outer_u'],control.ROW['fibre_T'])
    if tuple(map(cert.F,seed['curve']))!=model or tuples(seed['covariant_images'])!=images:raise ArithmeticError('coherent parent/fibre images differ')
    E=EllipticCurve(QQ,[QQ(str(c)) for c in model]);image_points=[E([QQ(str(c)) for c in P]) for P in images]
    cloud=[E([QQ(c) for c in P]) for P in seed['divisor_cloud']]
    if cloud[0]!=image_points[0] or any(2*cloud[i]!=image_points[i]-image_points[0] for i in range(1,14)):raise ArithmeticError('exact divisor halving identities differ')
    if [seed['divisor_cloud'][i] for i in seed['independent_column_indices']]!=seed['points']:raise ArithmeticError('seed subset differs')
    if maps['protocol_hash']!=digest(p) or data['protocol_hash']!=digest(p) or len(initial)!=rank or tuples(data['initial_state']['state']['reductions']['points'])!=initial:raise ArithmeticError('known-only subgroup/protocol differs')
    g=matrix(ZZ,maps['rounded_gram']);u=matrix(ZZ,maps['change_of_basis']);h=matrix(ZZ,maps['reduced_gram'])
    with localcontext() as ctx:
        ctx.prec=110
        expected=[[int((Decimal(v)*1000000).to_integral_value()) for v in row] for row in maps['metric_gram']]
    if g!=matrix(ZZ,expected) or not g.is_symmetric() or not g.is_positive_definite() or abs(u.det())!=1 or h!=u*g*u.transpose():raise ArithmeticError('rounded positive metric or unimodular transport differs')
    if [r['parity'] for r in maps['sample']]!=control.masks(p):raise ArithmeticError('fixed2047 masks differ')
    for r in maps['sample']:
        w=matrix(ZZ,1,rank,r['representative']);v=matrix(ZZ,1,rank,r['reduced_representative'])
        if w!=v*u or int((w*g*w.transpose())[0,0])!=r['metric_norm'] or any((int(w[0,j])-(r['parity']>>j))%2 for j in range(rank)):raise ArithmeticError('exact representative parity/norm differs')
    selected=sorted(maps['sample'],key=lambda r:(-r['metric_norm'],r['parity']))[:49]
    if maps['centres']!=selected or len(maps['rows'])!=49:raise ArithmeticError('fixed selected roster differs')
    geometry(data,maps,initial,{**p,'maps_path':D/'maps.json'})
    cloud_check(data,data['charts'],control.ART/('mestre_parent_calibration_'+control.ROW['id'].replace('-','_')+'_mod2_v1.json'),D/'result.json')
    if data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=49:raise ArithmeticError('complete status on short roster')
    print('EXACT2047 PARITIES;49 SELECTED;RATIONAL MAP/POINT PROVENANCE',len(data['charts']),'SEARCHED',flush=True)
if __name__=='__main__':
    for index in range(len(control.protocol()['rows'])):control.configure(index);main()

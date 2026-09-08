#!/usr/bin/env python3
"""Preserve the completed generic-height and Picard replay inputs together."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
H=LOCAL/'mestre-generic-height-v1';C=LOCAL/'mestre-picard-count-v1'
OUT=ART/'mestre_468_replay_bundle_v1.json'

def compute():
    hp=cert.read(H/'protocol.json');cp=cert.read(C/'protocol.json');paths=[Path(__file__).resolve(),H/'protocol.json',H/'ledger.json',C/'protocol.json',C/'ledger.json']
    if cert.read(H/'ledger.json')['status']!='PASS' or cert.read(C/'ledger.json')['status']!='PASS':raise ArithmeticError('completed source ledgers required')
    rows=[]
    for u in hp['outer_parameters']:
        h=H/('u'+u)/'height.json';s=LOCAL/'mestre-parent-calibration-v1'/('u'+u+'-unit')/'seed.json';counts=[];paths += [h,s]
        for r in cp['rows']:
            if r['outer_u']==u:
                p=C/r['id']/'counts.json';counts.append(cert.read(p));paths.append(p)
        rows.append({'outer_u':u,'generic_heights':cert.read(h),'specialized_seed':cert.read(s),'finite_surface_counts':counts})
    certificates={}
    for name in ('mestre_generic_height_independent_v1.json','mestre_base_involution_v1.json','mestre_parent_picard_and_saturation_v1.json','mestre_geometric_discriminant_v2.json','mestre_parent_fibre_geometry_v1.json'):
        p=ART/name;certificates[name]=cert.read(p);paths.append(p)
    return {'schema':'elliptic-curves.mestre-468-replay-bundle.v1','status':'PASS','rows':rows,'certificates':certificates,
      'height_protocol':hp,'point_count_protocol':cp,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
      'scope':'Complete retained generic equations, heights, component profiles, exact group relation words, finite-field counts, specialized saturation bases and derived Picard/discriminant certificates. No new computation or claim beyond those certificates; raw ledgers remain preserved.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve replay bundle')
        checkpoint(OUT,result)
    print('PASS6 DETERMINANT468 PARENT REPLAY INPUTS')

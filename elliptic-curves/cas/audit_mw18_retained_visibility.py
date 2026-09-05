#!/usr/bin/env python3
"""Retrospective exact visibility audit; never enumerates or changes a search.

Audits both signs of the eleven published MW17 complement representatives.
One of their directions has become generic on each cover. Pointwise visibility
is deliberately not identified with recovery of the quotient rational span.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt, log2
from pathlib import Path
import zipfile

from search_observability import point_visibility
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results'
BUNDLE = ART / 'elliptic-curves/mw18_deep_centre_comparison_v1.zip'
PUBLIC = ART / 'elliptic-curves/icarm_curve_refresh_475_573_overview_v1.json'
COVERS = ART / 'elkies-k3-r17-extreme-anchored-mw18-covers-v1.json'


def transport(model, points, raw):
    a1, a2, a3, a4, a6 = map(Q, model)
    c2, c4, c6 = a2+a1*a1/4, a4+a1*a3/2, a6+a3*a3/4
    A, B = c4-c2*c2/3, c6-c2*c4/3+2*c2**3/27
    rawA, rawB = map(Q, raw[-2:])
    square = B*rawA/(rawB*A)
    scale = Q(isqrt(square.numerator), isqrt(square.denominator))
    if scale**2 != square or A != scale**4*rawA or B != scale**6*rawB:
        raise ArithmeticError('public/native rational isomorphism failed')
    output = []
    for x, y in points:
        x, y = Q(x), Q(y)
        xr, yr = (x+c2/3)/scale**2, (y+(a1*x+a3)/2)/scale**3
        if yr**2 != xr**3+rawA*xr+rawB:
            raise ArithmeticError('transported public point is off curve')
        output.append(dict(x=str(xr), y=str(yr)))
    return output


def run(output):
    public = {r['id']: r for r in json.loads(PUBLIC.read_text())['snapshot']['records']}
    z = zipfile.ZipFile(BUNDLE)
    protocol = json.loads(z.read('protocol.json'))
    covers = json.loads(COVERS.read_text())
    # Historical eleven-point complement is the retained control's explicit list.
    hp = ART/'elliptic-curves/elkies_2026_bisection_specialization_controls_v1.json'
    hist = next(r for r in json.loads(hp.read_text())['fibres'] if r['parameter']=='-9529/5471')
    result = dict(schema='elliptic-curves.mw18-retained-visibility.v1', retrospective_only=True,
        inputs={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
                for p in (BUNDLE,PUBLIC,COVERS,hp,Path(__file__),Path(__file__).with_name('search_observability.py'),
                          Path(__file__).with_name('elkies_rank28.py'))},
        protocol_hash=protocol['protocol_hash'], cases={},
        classification={'inside':100000,'modest_extension_at_most':1000000},
        claim_boundary='Signed published representatives only, not all generic translates or rational spans. '
        'Large required height diagnoses the retained coordinates, not a proof that no alternative small representative exists. '
        'The historical gates and all ranks are unchanged.')
    for case_id, case in sorted(protocol['cases'].items()):
        if case_id.startswith('historical'):
            model, pts = GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS
            labels = hist['public_complement']['source_point_indices_one_based']
        else:
            cid = int(case_id.split('-')[1]); row = public[cid]
            model, pts = row['ainvs'], row['points']; labels = list(range(18,29))
        all_points = transport(model, pts, case['curve'])
        oracle = {f'P{i}':all_points[i-1] for i in labels}
        cr = dict(oracle=oracle, policies={})
        for policy in protocol['policies']:
            cell = json.loads(z.read(f'anchor-trial/cells/{case_id}--{policy}.json'))
            counts = Counter(); minima = {}; in_box = []; discrepancies = []
            for label, p in oracle.items():
                rows=[]
                for sign in (1,-1):
                    signed=dict(x=p['x'],y=str(sign*Q(p['y'])))
                    for i, ch in enumerate(cell['charts']):
                        v=point_visibility(ch['search'],signed)
                        counts[v['status']]+=1
                        entry=dict(point=label,sign=sign,chart=i,**v)
                        if v.get('in_completed_box'):in_box.append(entry)
                        if v['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(entry)
                        if v.get('minimum_affine_height') is not None:rows.append(entry)
                best=min(rows,key=lambda r:(r['minimum_affine_height'],r['sign'],r['chart']))
                h=best['minimum_affine_height']
                best.update(height_ratio_to_claimed_box=str(Q(h,100000)), log2_height=log2(h),
                    classification='INSIDE' if h<=100000 else 'MODEST_10X' if h<=1000000 else 'LARGE_COORDINATE_COST')
                minima[label]=best
            cr['policies'][policy]=dict(chart_count=len(cell['charts']),certified_recovered_gain=cell['certified_gain'],
                status_counts=dict(counts),minimum_by_signed_representative=minima,
                in_completed_coverage=in_box,coverage_discrepancies=discrepancies,
                classification_counts=dict(Counter(v['classification'] for v in minima.values())))
            print('VISIBILITY',case_id,policy,dict(counts),cr['policies'][policy]['classification_counts'],flush=True)
        result['cases'][case_id]=cr
    result['status']='PASS_NO_COVERAGE_DISCREPANCY' if not any(p['coverage_discrepancies'] for c in result['cases'].values() for p in c['policies'].values()) else 'COVERAGE_DISCREPANCY'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ART/'elliptic-curves/mw18_retained_visibility_v1.json')
    run(parser.parse_args().output)

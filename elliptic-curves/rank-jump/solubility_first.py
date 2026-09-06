#!/usr/bin/env python3
"""Parameter-only exact splitting first, quotient diagnostics second. No search."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import json
from math import isqrt, lcm
from pathlib import Path
import retrospective as r
from cover_experiment import evaluate, sqrtq

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/'SOLUBILITY_FIRST_PROTOCOL.json'
GEOMETRY = r.OUT/'rank_jump_solubility_first_geometry_v1.json.gz'
INPUT = r.OUT/'rank_jump_solubility_first_inputs_v1.json'
OUTPUT = r.OUT/'rank_jump_solubility_first_v1.json'
COHORT = r.OUT/'rank_jump_completed_cohort_inputs_v1.json'
ATLAS = r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
MODEL = r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
ATLAS_SHA = '78e037dc4170955b8f79ddce4d1d3e0c0d3e9bb8f9614644c59ccc7d605226c4'


def integer_form(coeffs):
    q = list(map(F, coeffs)); assert len(q) == 3
    d = lcm(*(x.denominator for x in q))
    # f=d^2*q preserves rational square values, including their signs.
    return [int(x*d*d) for x in q], d


def split_rows(geometry, parameters):
    rows = []
    for par in parameters:
        t = F(par['published_parameter']); a, b = t.numerator, t.denominator
        hits, zeros = [], []
        for cover in geometry:
            c0,c1,c2 = cover['integer_quadratic']
            n = c0*b*b+c1*a*b+c2*a*a
            if n == 0: zeros.append(cover['label'])
            elif n > 0:
                z = isqrt(n)
                if z*z == n:
                    hits.append({'label': cover['label'], 'square_root': str(F(z,b*cover['denominator_scale']))})
        rows.append(par | {'nonzero_square_hits': hits, 'branch_degeneracies': zeros})
    return rows


def capture():
    raw = ATLAS.read_bytes(); assert r.digest(raw) == ATLAS_SHA
    atlas = json.loads(raw)['bisections']
    covers = sorted((c for c in atlas if c['residual_chord']['construction_chart']=='finite'), key=lambda c:c['label'])
    assert len(covers) == 39119
    geometry = []
    for c in covers:
        q,d = integer_form(c['residual_chord']['q_coefficients'])
        geometry.append({'label':c['label'], 'integer_quadratic':q, 'denominator_scale':d})
    parameters = [{'source_id':x['source_id'], 'compact_parameter':x['parameter'],
                   'published_parameter':str(-(F(x['parameter'])+50)/26)}
                  for x in r.read(COHORT)['rows'] if x['phase']=='initial' and '08234' in x['family']]
    assert len(parameters)==32
    # No rank or point is passed to the square-value operation.
    splits = split_rows(geometry, parameters)
    labels = {h['label'] for row in splits for h in row['nonzero_square_hits']}
    model = r.read(MODEL)
    r.write_new(INPUT, {'schema':'rank-jump.solubility-first-inputs.v1',
                       'atlas_sha256':ATLAS_SHA, 'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
                       'source_cohort_sha256':r.digest(COHORT.read_bytes()),
                       'parameters':parameters,
                       'A':model['A_coefficients_low_to_high'], 'B':model['B_coefficients_low_to_high'],
                       'split_lift_maps':{c['label']:c['lifted_section'] for c in covers if c['label'] in labels},
                       'excluded_inverted_chart':['orbit-0c54f'],
                       'boundary':'Map selection uses only parameter square values, never exceptional points or ranks.'})
    payload = json.dumps(geometry, separators=(',',':'), sort_keys=True).encode()
    with GEOMETRY.open('xb') as f: f.write(gzip.compress(payload,mtime=0))
    print('Frozen',len(geometry),'covers and',len(parameters),'parameters; soluble maps',len(labels),flush=True)


def compute():
    inp = r.read(INPUT)
    assert inp['protocol_sha256'] == r.digest(PROTOCOL.read_bytes())
    assert inp['source_cohort_sha256'] == r.digest(COHORT.read_bytes())
    geometry = json.loads(gzip.decompress(GEOMETRY.read_bytes()))
    assert len(geometry)==39119
    splits = split_rows(geometry,inp['parameters'])
    # Only after exact solubility has been decided read observed rank/point data.
    cohort = {x['source_id']:x for x in r.read(COHORT)['rows'] if x['phase']=='initial' and '08234' in x['family']}
    rows = []
    for split in splits:
        old = cohort[split['source_id']]; t = F(split['published_parameter'])
        A,B = (evaluate(list(map(F,inp[k])),t) for k in ('A','B'))
        model, points = r.short(old['model'],old['points'])
        # Recover an exact short-model scaling, then verify both coefficients.
        u2 = sqrtq(F(model[3])/A); assert u2 is not None
        u = sqrtq(u2); assert u is not None and F(model[4]) == u**6*B
        primes = [s['prime'] for s in old['rank_certificate']['signatures']]
        blocks = [(p,r.roots_at(model[3],model[4],p)) for p in primes]
        generic = [r.point_signature(model,P,blocks) for P in points[:17]]
        retained = [r.point_signature(model,P,blocks) for P in points]
        assert r.rank(generic)==17 and r.rank(retained)==len(points)
        lifted, signatures = [], []
        for hit in split['nonzero_square_hits']:
            c = inp['split_lift_maps'][hit['label']]; root = F(hit['square_root'])
            v = {k:evaluate(list(map(F,c[k+'_coefficients'])),t) for k in ('x0','x1','y0','y1')}
            x,y = v['x0']+root*v['x1'],v['y0']+root*v['y1']
            assert y*y == x**3+A*x+B
            P = [str(u**2*x),str(u**3*y)]
            r.short(model,[P])
            sig = r.point_signature(model,P,blocks); signatures.append(sig)
            lifted.append(hit | {'point':P,'Kummer_fingerprint':sig})
        q = r.rank(generic+signatures)-17
        rows.append(split | {'observed_quotient_rank':len(points)-17, 'full_quotient_rank':'UNKNOWN',
                            'split_cover_count':len(lifted), 'lifted_points':lifted,
                            'soluble_quotient_rank_lower_bound':q,
                            'soluble_quotient_rank_upper_bound':len(lifted),
                            'certified_multi_direction_block':q>=2,
                            'rank_lower_bound_with_retained_points':r.rank(retained+signatures),
                            'short_model_scaling_u':str(u),
                            'generic_rank':17, 'finite_characters_may_miss_independence':True})
    return {'schema':'rank-jump.solubility-first.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,GEOMETRY,INPUT,COHORT,Path(__file__),HERE/'retrospective.py',HERE/'cover_experiment.py')},
            'layer':'solubility','square_tests':len(geometry)*len(rows), 'rows':rows,
            'split_count_distribution':dict(sorted(Counter(str(x['split_cover_count']) for x in rows).items())),
            'certified_blocks':[x['source_id'] for x in rows if x['certified_multi_direction_block']],
            'boundary':'A fixed generic construction on completed retrospective fibres. No full Selmer/rank assertion and no negative solubility conclusion outside this dictionary.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','build','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    else:
        data=compute()
        if a.mode=='check':assert r.read(OUTPUT)==data
        else:r.write_new(OUTPUT,data)
        print('PASS',data['split_count_distribution'],'blocks',data['certified_blocks'])
        for x in data['rows']:
            print(x['source_id'],'observed',x['observed_quotient_rank'],'splits',x['split_cover_count'],'certified quotient',x['soluble_quotient_rank_lower_bound'])

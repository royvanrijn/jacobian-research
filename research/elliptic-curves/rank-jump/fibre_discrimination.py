#!/usr/bin/env python3
"""Frozen multi-fibre solubility comparison; no curve point construction."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import json
from itertools import combinations
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import retrospective as r

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/'FIBRE_DISCRIMINATION_PROTOCOL.json'
INPUT = r.OUT/'rank_jump_fibre_discrimination_inputs_v1.json'
GEOMETRY = r.OUT/'rank_jump_fibre_discrimination_geometry_v1.json.gz'
OUTPUT = r.OUT/'rank_jump_fibre_discrimination_v1.json'
CSV = r.OUT/'rank_jump_fibre_discrimination_v1.csv'
SOURCES = ['rank_jump_retrospective_report_v1.json', 'compact_r17_initial_measurements_v1.json',
           'rank_jump_completed_cohort_panel_v1.json', 'full11952_64_r17_results_v1.json',
           'full11952_late64_r17_results_v1.json']


def hash_file(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()


def capture():
    oldgeo = r.OUT/'rank_jump_solubility_first_geometry_v1.json.gz'
    geometry = {'published-R17': json.loads(gzip.decompress(oldgeo.read_bytes()))}
    atlas = r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
    assert hash_file(atlas) == '78e037dc4170955b8f79ddce4d1d3e0c0d3e9bb8f9614644c59ccc7d605226c4'
    # Stream one record; avoid loading the 242MB atlas while another agent works.
    buf = None; extra = None
    with atlas.open() as f:
        for line in f:
            if line == '    {\n': buf = []
            if buf is not None: buf.append(line)
            if buf is not None and line in ('    },\n', '    }\n'):
                if any('"label": "orbit-0c54f"' in s for s in buf):
                    extra = json.loads(''.join(buf).rstrip().rstrip(',')); break
                buf = None
    assert extra and extra['residual_chord']['construction_chart'] == 'inverted_at_infinity'
    q = list(map(F, extra['residual_chord']['q_coefficients']))
    den = lcm(*(v.denominator for v in q))
    geometry['published-R17'].append({'label': extra['label'], 'integer_quadratic': [int(v*den*den) for v in q], 'denominator_scale': den})
    geometry['published-R17'].sort(key=lambda c: c['label'])
    assert len(geometry['published-R17']) == 39120
    alt = r.ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json'
    geometry['11952'] = []
    for c in r.read(alt)['bisections']:
        assert c['branch']['denominator_coefficients'] == ['1']
        q = list(map(F, c['branch']['numerator_coefficients'])); assert len(q) == 3
        den = lcm(*(v.denominator for v in q))
        geometry['11952'].append({'label': c['label'], 'integer_quadratic': [int(v*den*den) for v in q], 'denominator_scale': den})
    assert len(geometry['11952']) == 1024
    compact = r.OUT/'compact_six_r17_atlas_v1.json'
    family = next(f for f in r.read(compact)['families'] if f['family'] == '11952')
    matrix = list(map(F, family['base_matrix_a_b_c_d']))
    assert matrix[0]*matrix[3] != matrix[1]*matrix[2]
    rows = []
    def add(source, id, fam, param, generic, rank, cohort, phase, boxes=None, score=None, origin=None, coverage=None):
        assert rank >= generic
        native = None; dictionary = None
        if param is not None:
            t = F(param)
            if fam == 'published-R17': dictionary, native = 'published-R17', t
            elif fam == '08234': dictionary, native = 'published-R17', -(t+50)/26
            elif fam == '11952':
                a,b,c,d = matrix
                if c*t+d: dictionary, native = '11952', (a*t+b)/(c*t+d)
        rows.append({'observation_id': source+':'+id, 'id': id, 'family': fam, 'parameter': param,
                     'dictionary': dictionary, 'native_parameter': str(native) if native is not None else None,
                     'generic_subgroup_rank': generic, 'retained_rank_lower_bound': rank,
                     'retained_quotient_rank': rank-generic, 'full_curve_rank': 'UNKNOWN',
                     'cohort': cohort, 'phase': phase, 'boxes': boxes, 'coverage_note': coverage,
                     'score_units': score, 'origin': origin, 'source': source})
    for old in r.read(r.OUT/SOURCES[0])['panel']:
        add(SOURCES[0], old['id'], old['family'], old.get('parameter'), old['generic_subgroup_rank_exact'],
            old['certified_independent_subgroup_rank_exact'], old.get('cohort'), 'historic_or_prior',
            old.get('chart_count'), old.get('score_units'), old['origin'], 'inherited exposure; not matched to later trials')
    for old in r.read(r.OUT/SOURCES[1])['rows']:
        add(SOURCES[1], old['cohort']+':'+str(old['index']), 'published-R17', old['parameter'], 17,
            old['rank_certificate']['rank_lower_bound'], old['cohort'], 'initial_measurement', old['chart_count'],
            old['score_units'], 'scored_observation', 'full_box_count='+str(old['full_box_count']))
    for old in r.read(r.OUT/SOURCES[2])['rows']:
        if old['phase'] != 'initial': continue
        add(SOURCES[2], old['source_id'], old['family'], old['parameter'], old['generic_subgroup_rank_exact'],
            old['witness_subgroup_rank_exact'], old['cohort'], 'completed_initial', old['completed_initial_boxes'],
            origin='completed_initial', coverage='all declared initial boxes complete')
    for source in SOURCES[3:]:
        for old in r.read(r.OUT/source)['curves']:
            assert old['completed_boxes'] == old['declared_charts']
            assert len(old['generic_points']) == 17
            assert old['rank_lower_bound'] == len(old['points']) == old['rank_certificate']['rank_lower_bound']
            add(source, old['id'], old['family'], old['parameter'], 17, old['rank_lower_bound'], source,
                'completed_initial', old['completed_boxes'], origin='completed_scored_cohort', coverage='all declared initial boxes complete')
    assert len(rows) <= 600
    inputs = {'schema': 'rank-jump.fibre-discrimination-inputs.v1', 'rows': rows,
              'native11952_parameter_matrix': list(map(str, matrix)), 'native11952_parameter_convention': family['parameter_convention'],
              'sources': {str(p.relative_to(r.ROOT)): hash_file(p) for p in [*(r.OUT/n for n in SOURCES), oldgeo, atlas, alt, compact]},
              'protocol_sha256': hash_file(PROTOCOL),
              'scope': 'Retrospective labels and censored search outcomes; no new point or parameter search.'}
    r.write_new(INPUT, inputs)
    with GEOMETRY.open('xb') as f: f.write(gzip.compress(json.dumps(geometry, sort_keys=True, separators=(',',':')).encode(), mtime=0))
    print('CAPTURED', len(rows), 'observations', Counter(x['dictionary'] for x in rows), flush=True)


def resultant(q, w):
    c,b,a = q; f,e,d = w
    return (a*f-c*d)**2-(a*e-b*d)*(b*f-c*e)


def compute():
    inputs = r.read(INPUT); assert inputs['protocol_sha256'] == hash_file(PROTOCOL)
    geometry = json.loads(gzip.decompress(GEOMETRY.read_bytes()))
    # This pass uses only parameters and generic equations, not rank/point labels.
    addresses = sorted({(row['dictionary'], row['native_parameter']) for row in inputs['rows'] if row['dictionary']})
    assert sum(len(geometry[fam]) for fam,t in addresses) <= 10000000
    blocks = {}; tests = 0
    for fam, parameter in addresses:
        t = F(parameter); a,b = t.numerator,t.denominator; hits=[]; zeros=[]
        for cover in geometry[fam]:
            c0,c1,c2 = cover['integer_quadratic']; v = c0*b*b+c1*a*b+c2*a*a
            tests += 1
            if v == 0: zeros.append(cover['label'])
            elif v > 0 and isqrt(v)**2 == v: hits.append(cover)
        forms=[]
        for c in hits:
            content = gcd(*c['integer_quadratic']); q=[v//content for v in c['integer_quadratic']]
            forms.append({'label':c['label'], 'primitive_form':q, 'removed_content':str(content),
                          'content_is_square': isqrt(content)**2 == content})
        pairs=[{'indices':[i,j], 'resultant':str(resultant(q['primitive_form'],w['primitive_form']))}
               for i,q in enumerate(forms) for j,w in enumerate(forms[:i])]
        n=len(hits)
        disjoint=all(int(x['resultant']) for x in pairs) and all(q['primitive_form'][1]**2 != 4*q['primitive_form'][0]*q['primitive_form'][2] for q in forms)
        degree=2**n if disjoint else 'UNKNOWN'
        genus=(0 if n == 0 else 1+2**(n-1)*(n-2)) if disjoint else 'UNKNOWN'
        blocks[fam+':'+parameter]={'dictionary':fam,'native_parameter':parameter,'compatible_cover_count':n,
            'branch_degeneracies':zeros,'compatible_forms':forms,'pair_resultants':pairs,
            'disjoint_branch_divisors':bool(disjoint),'simultaneous_carrier_degree':degree,'simultaneous_carrier_genus':genus,
            'split_triple_count':comb(n,3),'split_quartet_count':comb(n,4),
            'split_triples_with_some_fourth':comb(n,3) if n>=4 else 0,
            'triple_fourth_incidences':comb(n,3)*(n-3) if n>=4 else 0,
            'specialized_defect_on_compatible_subset':0,
            'collision_support_and_realizable_span':'UNKNOWN_PENDING_ARITHMETIC_AUDIT' if n>=2 else 'NOT_APPLICABLE',
            'rational_relation_component_count':'UNKNOWN_NO_COMMON_RELATION_DICTIONARY' if n>=2 else 'NOT_APPLICABLE'}
    rows=[]
    for old in inputs['rows']:
        key=old['dictionary']+':'+old['native_parameter'] if old['dictionary'] else None
        block=blocks.get(key)
        rows.append(old | {'block_key':key, 'dictionary_size':len(geometry[old['dictionary']]) if key else None,
            'dictionary_coverage':'COMPLETE_NATIVE_ATLAS' if old['dictionary']=='published-R17' else 'PARTIAL_1024' if key else 'UNTESTED_NO_VALID_ATLAS_OR_PARAMETER',
            'compatible_cover_count':block['compatible_cover_count'] if block else None,
            'simultaneous_carrier_genus':block['simultaneous_carrier_genus'] if block else None,
            'simultaneous_carrier_degree':block['simultaneous_carrier_degree'] if block else None,
            'split_quartet_count':block['split_quartet_count'] if block else None})
    summaries=[]
    groups=sorted({(x['dictionary'],str(x['cohort'])) for x in rows if x['dictionary']})
    for fam,cohort in groups:
        subset=[x for x in rows if x['dictionary']==fam and str(x['cohort'])==cohort]
        for gain in sorted({x['retained_quotient_rank'] for x in subset}):
            ss=[x for x in subset if x['retained_quotient_rank']==gain]
            summaries.append({'dictionary':fam,'cohort':cohort,'gain':gain,'observations':len(ss),
                'cover_count_distribution':dict(sorted(Counter(str(x['compatible_cover_count']) for x in ss).items()))})
    paired=[]
    def height(t):t=F(t);return max(abs(t.numerator),t.denominator)
    for high in rows:
        if not high['dictionary'] or high['retained_quotient_rank']<7:continue
        lows=[x for x in rows if x['family']==high['family'] and x['cohort']==high['cohort']
              and x['phase']==high['phase'] and x['boxes']==high['boxes'] and x['retained_quotient_rank']==0]
        if not lows:continue
        def key(low):
            h,l=height(high['parameter']),height(low['parameter'])
            return F(max(h,l),min(h,l)), abs(F(high['parameter'])-F(low['parameter'])),low['observation_id']
        low=min(lows,key=key)
        paired.append({'high':high['observation_id'],'low':low['observation_id'], 'height_ratio':str(key(low)[0]),
                       'high_gain':high['retained_quotient_rank'],'cover_counts':[high['compatible_cover_count'],low['compatible_cover_count']],
                       'same_family_cohort_phase_boxes':True,'score_matched':False})
    return {'schema':'rank-jump.fibre-discrimination.v1','status':'PASS','rows':rows,'blocks':blocks,
            'square_tests':tests,'unique_tested_addresses':len(blocks),'grouped_results':summaries,'height_matched_pairs':paired,
            'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (INPUT,GEOMETRY,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
            'boundary':r.read(PROTOCOL)['failure_semantics']}


def csv_text(data):
    import csv,io
    keys=['observation_id','family','parameter','native_parameter','generic_subgroup_rank','retained_rank_lower_bound',
          'retained_quotient_rank','full_curve_rank','cohort','phase','boxes','score_units','dictionary_coverage',
          'dictionary_size','compatible_cover_count','simultaneous_carrier_genus','simultaneous_carrier_degree','split_quartet_count']
    stream=io.StringIO(); w=csv.DictWriter(stream,fieldnames=keys,lineterminator='\n');w.writeheader()
    for row in data['rows']:w.writerow({k:row[k] for k in keys})
    return stream.getvalue()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['capture','build','check']);args=parser.parse_args()
    if args.mode=='capture':capture()
    else:
        data=compute()
        if args.mode=='build':
            r.write_new(OUTPUT,data)
            with CSV.open('x') as f:f.write(csv_text(data))
        else:assert r.read(OUTPUT)==data and CSV.read_text()==csv_text(data)
        print('PASS',len(data['rows']),'observations',data['unique_tested_addresses'],'addresses',data['square_tests'],'square tests',flush=True)
        for x in data['grouped_results']: print(x,flush=True)

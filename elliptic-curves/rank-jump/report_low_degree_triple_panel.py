#!/usr/bin/env python3
"""Join inherited retained gains only after the generic relation evaluation."""
import argparse
import csv
import io
from pathlib import Path
import retrospective as r
from fibre_discrimination import hash_file

HERE=Path(__file__).resolve().parent
PANEL=r.OUT/'rank_jump_low_degree_triple_panel_v1.json'
VERIFIED=r.OUT/'rank_jump_low_degree_triple_panel_verification_v1.json'
PAIR_REPORT=r.OUT/'rank_jump_degree_one_relation_discrimination_v1.json'
OUTPUT=r.OUT/'rank_jump_low_degree_triple_discrimination_v1.json'
CSV=r.OUT/'rank_jump_low_degree_triple_discrimination_v1.csv'


def compute():
    panel=r.read(PANEL);verified=r.read(VERIFIED);pairs=r.read(PAIR_REPORT)
    assert verified['status']==pairs['status']=='PASS'
    for doc in (panel,verified,pairs):
        for path,sha in doc['bindings'].items():assert hash_file(r.ROOT/path)==sha
    vv={x['block_key']:x for x in verified['rows']};pp={x['block_key']:x for x in pairs['rows']};rows=[]
    for x in panel['rows']:
        v=vv[x['block_key']];p=pp[x['block_key']]
        for k in v:assert v[k]==x[k]
        rows.append({'parameter':x['parameter'],'block_key':x['block_key'],'generic_subgroup_rank':17,
            'best_retained_quotient_rank':p['best_retained_quotient_rank'],'full_curve_rank':'UNKNOWN',
            'compatible_cover_count':x['compatible_cover_count'],'co_split_triples':x['co_split_triples'],
            'eligible_total_degree_6':x['eligible_by_total_degree']['6'],'eligible_total_degree_8':x['eligible_by_total_degree']['8'],
            'degree_6_incidence_count':x['incidences_by_total_degree']['6'],'degree_8_incidence_count':x['incidences_by_total_degree']['8'],
            'pair_incidence_count':p['rational_degree_one_components_at_parameter_mod_conjugation'],
            'pair_relation_rank':x['pair_relation_rank'],'triple_relation_rank':x['triple_relation_rank'],
            'combined_relation_rank':x['combined_relation_rank'],'extra_relation_rank_beyond_pairs':x['extra_relation_rank_beyond_pairs'],
            'native_quotient_upper_bound':x['native_quotient_upper_bound'],'observation_ids':p['observation_ids']})
    totals={k:sum(x[k] for x in rows) for k in ('co_split_triples','eligible_total_degree_6','eligible_total_degree_8',
                                             'degree_6_incidence_count','degree_8_incidence_count')}
    by_gain={str(g):{'addresses':sum(x['best_retained_quotient_rank']==g for x in rows),
                    **{k:sum(x[k] for x in rows if x['best_retained_quotient_rank']==g) for k in totals}}
             for g in sorted({x['best_retained_quotient_rank'] for x in rows})}
    return {'schema':'rank-jump.low-degree-triple-discrimination.v1','status':'PASS','rows':rows,'totals':totals,
        'by_best_retained_gain':by_gain,
        'bindings':{str(p.relative_to(r.ROOT)):hash_file(p) for p in (PANEL,VERIFIED,PAIR_REPORT,Path(__file__),HERE/'retrospective.py',HERE/'fibre_discrimination.py')},
        'boundary':'All labelled rational incidences at fixed parameters for this degree-six/eight integral-translate dictionary, modulo simultaneous conjugation. Counts are not distinct elliptic points, all global rational factors, or independent directions. Full ranks are UNKNOWN; a retained zero is censored. No new parameter or point search. The A1/MW16 and untransported historic-extreme coverage gaps persist.'}


def csv_text(data):
    keys=[k for k in data['rows'][0] if k not in ('block_key','observation_ids')]
    out=io.StringIO();w=csv.DictWriter(out,fieldnames=keys,lineterminator='\n');w.writeheader()
    for row in data['rows']:w.writerow({k:row[k] for k in keys})
    return out.getvalue()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':
        r.write_new(OUTPUT,d)
        with CSV.open('x') as f:f.write(csv_text(d))
    else:assert r.read(OUTPUT)==d and CSV.read_text()==csv_text(d)
    print('PASS',len(d['rows']),'addresses',d['totals'],flush=True)

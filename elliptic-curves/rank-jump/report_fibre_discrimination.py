#!/usr/bin/env python3
"""Join separately verified panel metrics without dropping unavailable families."""
import argparse
from collections import Counter
import csv
import io
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
CENSUS=r.OUT/'rank_jump_fibre_discrimination_v1.json'
VERIFY=r.OUT/'rank_jump_fibre_discrimination_verification_v1.json'
OUTPUT=r.OUT/'rank_jump_fibre_discrimination_summary_v1.json'
CSV=r.OUT/'rank_jump_fibre_discrimination_metrics_v1.csv'


def compute():
    c=r.read(CENSUS);v=r.read(VERIFY);arithmetic={x['block_key']:x for x in v['rows']};rows=[]
    for old in c['rows']:
        row=dict(old);block=c['blocks'].get(old['block_key']);detail=arithmetic.get(old['block_key'])
        if detail:
            row.update({k:detail[k] for k in ['collision_support_size','realizable_defect_span_exact']})
        elif block is not None:
            row['collision_support_size']=0
            row['realizable_defect_span_exact']=0 if block['compatible_cover_count']==1 else 'NOT_APPLICABLE'
        else:
            row['collision_support_size']=row['realizable_defect_span_exact']='UNTESTED'
        for key in ['split_triple_count','split_triples_with_some_fourth','triple_fourth_incidences','rational_relation_component_count']:
            row[key]=block[key] if block else 'UNTESTED'
        rows.append(row)
    unique={}
    for row in rows:
        if row['block_key'] and (row['block_key'] not in unique or row['retained_quotient_rank']>unique[row['block_key']]['retained_quotient_rank']):
            unique[row['block_key']]=row
    groups=[]
    for fam in ['published-R17','11952']:
        for lo,hi in [(0,0),(1,4),(5,8),(9,100)]:
            subset=[x for x in unique.values() if x['dictionary']==fam and lo<=x['retained_quotient_rank']<=hi]
            groups.append({'dictionary':fam,'best_retained_gain_bin':[lo,hi], 'unique_parameter_addresses':len(subset),
                           'cover_count_distribution':dict(sorted(Counter(str(x['compatible_cover_count']) for x in subset).items()))})
    summary={'schema':'rank-jump.fibre-discrimination-summary.v1','status':'PASS','rows':rows,
             'observations':len(rows),'tested_observations':sum(x['dictionary'] is not None for x in rows),
             'unique_tested_parameter_addresses':len(unique),'untested_observations':sum(x['dictionary'] is None for x in rows),
             'unique_address_groups':groups,'all_nontrivial_blocks_have_maximal_defect_span':all(x['realizable_defect_span_exact']==x['compatible_cover_count']-1 for x in arithmetic.values()),
             'nontrivial_blocks':len(arithmetic),'factored_pair_resultants':sum(x['pair_count'] for x in arithmetic.values()),
             'independent_local_witnesses':sum(len(x['independent_local_witnesses']) for x in arithmetic.values()),
             'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (CENSUS,VERIFY,Path(__file__),HERE/'retrospective.py')},
             'boundary':'Address groups use the best retained subgroup observation, not an exact rank. Different addresses might represent isomorphic curves; no cross-family statistical pooling or independence claim. Combinatorial extension counts do not add a feature beyond cover count. Relation-component census remains uncomputed.'}
    return summary


def csv_text(d):
    keys=['observation_id','family','parameter','generic_subgroup_rank','retained_rank_lower_bound','retained_quotient_rank','full_curve_rank',
          'cohort','phase','boxes','score_units','dictionary_coverage','dictionary_size','compatible_cover_count',
          'realizable_defect_span_exact','collision_support_size','rational_relation_component_count',
          'split_triple_count','split_triples_with_some_fourth','triple_fourth_incidences','split_quartet_count',
          'simultaneous_carrier_genus','simultaneous_carrier_degree']
    s=io.StringIO();w=csv.DictWriter(s,fieldnames=keys,lineterminator='\n');w.writeheader()
    for row in d['rows']:w.writerow({k:row[k] for k in keys})
    return s.getvalue()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':
        r.write_new(OUTPUT,d)
        with CSV.open('x') as f:f.write(csv_text(d))
    else:assert r.read(OUTPUT)==d and CSV.read_text()==csv_text(d)
    print({k:v for k,v in d.items() if k not in ['rows','bindings']})

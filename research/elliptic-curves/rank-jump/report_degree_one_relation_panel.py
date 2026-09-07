#!/usr/bin/env python3
"""Join the uniform relation census to existing gain labels after computation."""
import argparse
import csv
import io
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
CENSUS=r.OUT/'rank_jump_fibre_discrimination_v1.json'
PANEL=r.OUT/'rank_jump_degree_one_relation_panel_v1.json'
VERIFIED=r.OUT/'rank_jump_degree_one_relation_panel_verification_v1.json'
TRIPLES=r.OUT/'rank_jump_trace_zero_triple_panel_v1.json'
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_degree_one_relation_discrimination_v1.json'
CSV=r.OUT/'rank_jump_degree_one_relation_discrimination_v1.csv'


def compute():
    census=r.read(CENSUS);panel=r.read(PANEL);verified=r.read(VERIFIED);triples=r.read(TRIPLES);inp=r.read(INPUT)
    proof={x['block_key']:x for x in verified['rows']};tr={x['block_key']:x for x in triples['rows']}
    rows=[]
    for old in panel['rows']:
        key=old['block_key'];observations=[x for x in census['rows'] if x['block_key']==key]
        gain=max(x['retained_quotient_rank'] for x in observations)
        v=proof[key]
        canonical=[]
        for edge in old['signed_quotient_edges']:
            # The producer's generic_word is the translate for the chosen BRANCHES.
            # Subtract native traces when writing an equation in canonical +root points.
            word=edge['generic_word'].copy()
            for label,sign in zip(edge['labels'],edge['signs']):
                if sign==-1:word=[a-b for a,b in zip(word,inp['covers'][label]['trace'])]
            canonical.append(edge|{'branch_sum_translate_word':edge['generic_word'],
                                   'canonical_signed_relation_generic_word':word})
        rows.append({'parameter':old['parameter'],'block_key':key,'generic_subgroup_rank':17,
            'best_retained_quotient_rank':gain,'full_curve_rank':'UNKNOWN',
            'compatible_cover_count':v['covers'],'eligible_degree_one_pairs':v['eligible_degree_one_pairs'],
            'rational_degree_one_components_at_parameter_mod_conjugation':v['incidences'],
            'pair_relation_rank':v['signed_relation_rank'],'native_quotient_upper_bound':v['native_quotient_upper_bound'],
            'co_split_triples':tr[key]['co_split_triples'],'zero_residual_triples':len(tr[key]['trace_even_triples']),
            'triple_total_intersection_degree_lower_bound':6 if tr[key]['co_split_triples'] else 'NOT_APPLICABLE',
            'canonical_signed_relations':canonical,'observation_ids':[x['observation_id'] for x in observations]})
    return {'schema':'rank-jump.degree-one-relation-discrimination.v1','status':'PASS','rows':rows,
            'degree_one_incidence_count':sum(x['rational_degree_one_components_at_parameter_mod_conjugation'] for x in rows),
            'addresses_with_degree_one_incidence':sum(x['rational_degree_one_components_at_parameter_mod_conjugation']>0 for x in rows),
            'eligible_pair_count':sum(x['eligible_degree_one_pairs'] for x in rows),
            'relation_incidence_counts_by_best_retained_gain':{str(g):sum(x['rational_degree_one_components_at_parameter_mod_conjugation'] for x in rows if x['best_retained_quotient_rank']==g) for g in sorted({x['best_retained_quotient_rank'] for x in rows})},
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (CENSUS,PANEL,VERIFIED,TRIPLES,INPUT,Path(__file__),HERE/'retrospective.py')},
            'boundary':'A uniform degree-one dictionary on the complete published-R17 atlas, evaluated at frozen nonbranch fibres. Counts are not all low-degree rational components. Pair relations bound only the constructed native subgroup. The triple lower bound concerns total intersection degree; linear factors of higher-degree schemes remain possible.'}


def csv_text(data):
    keys=['parameter','generic_subgroup_rank','best_retained_quotient_rank','full_curve_rank','compatible_cover_count',
          'eligible_degree_one_pairs','rational_degree_one_components_at_parameter_mod_conjugation',
          'pair_relation_rank','native_quotient_upper_bound','co_split_triples','zero_residual_triples','triple_total_intersection_degree_lower_bound']
    s=io.StringIO();w=csv.DictWriter(s,fieldnames=keys,lineterminator='\n');w.writeheader()
    for row in data['rows']:w.writerow({k:row[k] for k in keys})
    return s.getvalue()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':
        r.write_new(OUTPUT,d)
        with CSV.open('x') as f:f.write(csv_text(d))
    else:assert r.read(OUTPUT)==d and CSV.read_text()==csv_text(d)
    print('PASS',len(d['rows']),'addresses;',d['degree_one_incidence_count'],'incidences on',d['addresses_with_degree_one_incidence'],'addresses')
    print(d['relation_incidence_counts_by_best_retained_gain'])

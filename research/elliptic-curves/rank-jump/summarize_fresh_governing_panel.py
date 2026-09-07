#!/usr/bin/env python3
"""Join masked arithmetic to retrospective labels only after all workers finish."""
import argparse
import csv
from pathlib import Path
import retrospective as r
import fresh_governing_panel as base
import fresh_governing_completion as completion
import fresh_governing_octics as octics
import verify_fresh_governing_panel as verifier

OUTPUT=r.OUT/'rank_jump_fresh_governing_comparison_v1.json'
CSV=OUTPUT.with_suffix('.csv')


def build(check=False):
    raw=r.read(base.OUTPUT);patch={x['token']:x for x in r.read(completion.OUTPUT)['rows']}
    labels=r.read(base.MANIFEST);metadata={x['token']:x for x in labels['rows']}
    octic={x['token']:x for x in r.read(octics.OUTPUT)['rows']}
    verified=r.read(verifier.OUTPUT);assert verified['status']=='PASS'
    rows=[]
    for old in raw['rows']:
        token=old['token'];label=metadata[token];local=old['local'];repair=patch.get(token,{})
        chosen=repair if repair else local;k=chosen.get('strict_generic_dimension')
        upper=chosen.get('strict_generic_dimension_upper_bound',k)
        M=chosen.get('minus_twist_CT_matrix',local.get('minus_twist_CT_matrix'))
        rank=chosen.get('minus_twist_CT_rank')
        o=octic[token];assert o['status']=='PASS'
        rows.append({**label,'inherited_strict_dimension':k,'inherited_strict_dimension_upper_bound':upper,
          'inherited_strict_nonzero_class_count':2**k-1 if k is not None else None,
          'inherited_unramified_extension_degree_over_cubic':2**k if k is not None else None,
          'inherited_basis_pair_count':k*(k-1)//2 if k is not None else None,
          'inherited_CT_switch_rank':rank,'inherited_CT_hyperbolic_planes':rank//2 if rank is not None else None,
          'inherited_restricted_CT_radical_dimension':k-rank if k is not None and rank is not None else None,
          'inherited_original_CT_rank':0 if k is not None else None,
          'inherited_CT_switch_matrix':M,
          'fixed_generic_pair_class_field_degree':o['joint_class_field_degree'],
          'fixed_generic_pair_governing_field_degree':o['governing_field_degree'],
          'fixed_generic_pair_inert_primes_checked':len(o['inert_prime_table']),
          'additional_Selmer_classes':'UNKNOWN','additional_quotient_CT':'UNKNOWN','full_CT':'UNKNOWN',
          'independent_class_group_status':old['independent']['status'],
          'independent_class_group_failure':'PARI_256MiB_STACK_LIMIT' if old['independent'].get('reason')=='worker failure' else old['independent'].get('reason'),
          'native_cover_count':'NOT_COMPUTED_IN_THIS_CT_PANEL',
          'simultaneous_solubility_carrier':'UNKNOWN'})
    lookup={x['token']:x for x in rows};pairs=[]
    for pair in labels['pairs']:
        h,l=lookup[pair['high']],lookup[pair['low']]
        pairs.append({**pair,'high_inherited_strict_dimension':h['inherited_strict_dimension'],
          'low_inherited_strict_dimension':l['inherited_strict_dimension'],
          'high_CT_switch_rank':h['inherited_CT_switch_rank'],'low_CT_switch_rank':l['inherited_CT_switch_rank'],
          'field_degree_discriminates':False,'full_quotient_CT_comparison':'UNKNOWN'})
    paths=(Path(__file__),base.INPUT,base.MANIFEST,base.OUTPUT,completion.OUTPUT,octics.OUTPUT,verifier.OUTPUT)
    result={'schema':'rank-jump.fresh-governing-comparison.v1','bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
      'rows':rows,'pairs':pairs,'feature_layers':{
        'inherited_strict_dimension':'incidence, inherited subspace only',
        'inherited_CT_switch_rank':'solubility obstruction on the -1 twist, not original rank gain',
        'fixed_generic_pair_governing_field_degree':'universal cochain structure; non-discriminating under these hypotheses',
        'coefficient_bits_and_completed_boxes':'matching/exposure metadata, not incidence'},
      'protocol_deviation':'The export used parameter-height bit length >13 for the large11952 control pool, whereas the prose protocol said height >4096. Thus4286/1881 used the compact pool. No case was rematched after arithmetic; comparisons are descriptive, not a preregistered aggregate test.',
      'selection_scope':'Retrospective, labels used to choose cases. Reused controls are not independent repetitions. Low gains are censored. Arithmetic workers receive no labels or exceptional coordinates.',
      'claim_boundary':'The requested additional-quotient governing/CT discriminator remains uncomputed. This is a verified inherited baseline and a precise independent-class bottleneck, not a rank predictor.'}
    if check:assert r.read(OUTPUT)==result
    else:
        r.write_new(OUTPUT,result)
        fields=['token','id','family','parameter','retained_rank_lower_bound','generic_rank','observed_quotient_rank','inherited_strict_dimension','inherited_strict_dimension_upper_bound','inherited_CT_switch_rank','inherited_restricted_CT_radical_dimension','fixed_generic_pair_governing_field_degree','additional_quotient_CT']
        with CSV.open('x',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print('PASS comparison:',len(rows),'fibres;',sum(x['inherited_strict_dimension'] is not None for x in rows),'exact inherited kernels; additional quotient UNKNOWN on all rows')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();build(a.mode=='check')

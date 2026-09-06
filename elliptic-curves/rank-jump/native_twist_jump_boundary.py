#!/usr/bin/env python3
"""Translate a fixed generic twist bound into a retrospective jump boundary."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
from cover_experiment import evaluate, sqrtq

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_native_twist_frobenius_inputs_v1.json'
MOMENTS=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
PARITY=r.OUT/'rank_jump_native_twist_moment_parity_v1.json'
VERIFY=r.OUT/'rank_jump_native_twist_frobenius_verification_v2.json'
FIBRES=r.OUT/'rank_jump_solubility_first_v1.json'
OUTPUT=r.OUT/'rank_jump_native_twist_jump_boundary_v1.json'


def compute():
    inp=r.read(INPUT);mom=r.read(MOMENTS);par=r.read(PARITY);verification=r.read(VERIFY);fibres=r.read(FIBRES)
    for data in (inp,mom,par,verification,fibres):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    case=next(row for row in fibres['rows'] if row['source_id']=='08234-009')
    root=sqrtq(evaluate(inp['q'],F(case['published_parameter'])));assert root and root==F(4307419,1937)
    assert 'orbit-1795d' in [x['label'] for x in case['nonzero_square_hits']]
    original=next(row for row in mom['rows'] if row['id']=='original')
    lower=original['arithmetic_generic_rank_lower_bound'];upper=par['rows'][0]['parity_refined_upper_bound']
    assert lower==upper==case['generic_rank']==17
    extra_lower=mom['rows'][1]['arithmetic_generic_rank_lower_bound']
    extra_upper=par['rows'][1]['parity_refined_upper_bound'];assert extra_lower==1 and extra_upper==7
    witness=case['generic_rank']+case['observed_quotient_rank'];assert witness==25
    return {'schema':'rank-jump.native-twist-jump-boundary.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,MOMENTS,PARITY,VERIFY,FIBRES,Path(__file__),HERE/'retrospective.py',HERE/'cover_experiment.py')},
            'source_id':case['source_id'],'published_parameter':case['published_parameter'],
            'primitive_cover_square_root':str(root),'marked_generic_rank':17,
            'native_twist_arithmetic_generic_rank_bounds':[extra_lower,extra_upper],
            'full_quadratic_pullback_arithmetic_generic_rank_bounds':[lower+extra_lower,upper+extra_upper],
            'retained_witness_rank':witness,'retained_witness_quotient_rank':case['observed_quotient_rank'],
            'witness_directions_outside_full_specialized_pullback_span_at_least':witness-(upper+extra_upper),
            'original_fibre_full_rank':'UNKNOWN','native_twist_exact_rank':'UNKNOWN',
            'argument':'Character decomposition gives full pullback rank=17+twist rank<=24. At either rational lift of this smooth parameter, the specialized image has rank<=24. The retained independent witness subgroup has rank25, so its image modulo that entire specialized rational span has dimension at least1.',
            'boundary':'This excludes only the claim that every retained direction is supplied by generic sections on this one quadratic cover. Other covers, product characters and specialization-only directions remain possible.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=compute()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS: full generic quadratic pullback rank18..24; witnessed fibre rank25; at least1 direction outside its specialized span')

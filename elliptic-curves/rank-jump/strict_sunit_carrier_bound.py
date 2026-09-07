#!/usr/bin/env python3
"""Retrospective S-unit carrier bounds from already certified Artin matrices."""
import argparse
from pathlib import Path
import retrospective as r

SOURCE=r.OUT/'rank_jump_half_ideal_artin_completion_v1.json'
OUTPUT=r.OUT/'rank_jump_strict_sunit_carrier_bound_v1.json'


def compute():
    rows=[]
    for old in r.read(SOURCE)['rows']:
        M=old['matrix_rows'];n=len(M);assert all(len(row)==n for row in M)
        # Exhaustive kernel, independently of the original matrix-rank routine.
        kernel=[v for v in range(1<<n) if all(sum(row[j]*((v>>j)&1) for j in range(n))%2==0 for row in M)]
        dim=(len(kernel)).bit_length()-1;assert len(kernel)==1<<dim and n-dim==old['artin_rank']
        rows.append({'case_index':old['case_index'],'id':old['id'],'known_strict_dimension':n,
            'artin_rank':n-dim,'right_kernel_masks':kernel,
            'strict_S_unit_dimension_in_known_space_upper_bound':dim,
            'known_strict_quotient_by_S_units_dimension_lower_bound':n-dim})
    return {'schema':'rank-jump.strict-sunit-carrier-bound.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),SOURCE)},
        'boundary':'Oracle-derived retrospective diagnostic only. Existing Artin matrices come from known exceptional points. They are not inputs to the separate equation-only principalization experiment or any prospective selector. A zero Artin kernel vector is not certified S-unit; no claim about S-units plus non-strict generic corrections.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print([(row['id'],row['strict_S_unit_dimension_in_known_space_upper_bound']) for row in result['rows']])

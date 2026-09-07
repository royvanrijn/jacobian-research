#!/usr/bin/env python3
"""Generic-only exact consequences and an explicitly artificial counting baseline."""
import argparse
from fractions import Fraction
from math import prod
from pathlib import Path
import retrospective as r
import generic_local_contraction as local

OUTPUT=r.OUT/'rank_jump_generic_contraction_consequences_v1.json'


def gaussian(n,k):
    if k<0 or k>n:return 0
    return prod(2**(n-i)-1 for i in range(k))//prod(2**(k-i)-1 for i in range(k))


def surjective_count(ell,g,d):
    return 2**(d*(ell-g))*gaussian(ell-d,g-d)


def build(check=False):
    data=r.read(local.OUTPUT);assert data['bindings']==local.bindings();rows=[]
    for row in data['rows']:
        assert row['status']=='PASS'
        G,O,T=row['joint_generic_signatures'],row['original_product_basis'],row['twist_product_basis']
        replay=local.linear_certificate(G,O,T)
        assert all(row[k]==v for k,v in replay.items())
        ell=row['local_product_dimension'];g=row['joint_generic_dimension']
        d=row['local_change_dimension'];e=row['generic_quotient_dimension'];m=row['generic_point_count']
        assert d==e
        fraction=Fraction(surjective_count(ell,g,d),gaussian(ell,g))
        rows.append({'case_index':row['case_index'],'generic_rank_from_retained_family_certificate':m,
            'generic_local_image_dimension':g,'local_product_dimension':ell,
            'exact_Selmer_dimension_drop':d,
            'generic_common_Selmer_dimension':m-d,
            'full_Selmer_identity':'Sel2(original) = Sel2(minus twist) + generic Kummer subgroup',
            'quotient_identity':'Sel2(original)/G = Sel2(minus twist)/(G intersect Sel2(minus twist))',
            'toy_uniform_subspace_surjectivity_fraction':str(fraction),
            'toy_model_only':'Uniform g-dimensional subspaces of an ell-dimensional binary space; no claim that arithmetic generic images follow this distribution.',
            'boundary':'The original/twist Selmer dimension difference is exact. Neither individual dimension, rational rank, Sha dimension nor prospective performance is inferred.'})
    result={'schema':'rank-jump.generic-contraction-consequences.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),local.OUTPUT)},'rows':rows}
    if check:assert r.read(OUTPUT)==result;print('PASS generic contraction consequences')
    else:r.write_new(OUTPUT,result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();build(args.mode=='check')

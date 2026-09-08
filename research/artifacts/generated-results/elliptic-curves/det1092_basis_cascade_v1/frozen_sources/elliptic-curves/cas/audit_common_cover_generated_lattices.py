#!/usr/bin/env python3
"""Reduce the actual common-cover subgroups, keeping their closures separate."""
import argparse
import json
from pathlib import Path
import search_low_height_mw_sublattices as base


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    curves=[];inputs={}
    for curve in (245,302):
        path=base.OUT/f'common_cover_mw_sublattices_v1_{curve}_selection.json'
        inputs[str(path.relative_to(base.ROOT))]=base.digest(path)
        selected=base.read(path)
        cloudpath=base.OUT/f'low_height_mw_sublattices_v1_{curve}_cloud.json.gz'
        inputs[str(cloudpath.relative_to(base.ROOT))]=base.digest(cloudpath)
        h=base.matrix_literal(base.read(cloudpath)['height_gram']);records=[]
        for candidate in selected['finalists']:
            r=candidate['rank'];b=base.gp_matrix(candidate['basis_rows'])
            program=f"""
default(realprecision,80);H={h};B={b};G=B*H*B~;U=qflllgram(G);L=U~*B;R=L*H*L~;
if(abs(matdet(U))!=1 || mathnf(B~)!=mathnf(L~),error("changed generated lattice"));
if(abs(matdet(G)/({candidate['determinant']})-{candidate['generated_index_in_primitive_closure']}^2)>1e-45,error("index determinant identity"));
Q=qfminim(R,,,2);print("DET|",matdet(G));print("MIN|",Q[2]);
print("H|",vector({r},j,R[j,j]));print("B|",vector({r},j,Vec(L[j,])));
print("G|",vector({r},j,Vec(R[j,])));
"""
            lines=base.run_gp(program,timeout=120);data=dict(x.split('|',1) for x in lines)
            records.append({'candidate_index':candidate['candidate_index'],'rank':r,
                'actual_generated_determinant':data['DET'],'actual_generated_minimum':data['MIN'],
                'actual_generated_lll_heights':json.loads(data['H']),
                'actual_generated_reduced_basis_rows':json.loads(data['B']),
                'actual_generated_reduced_gram':json.loads(data['G']),
                'index_in_primitive_closure':candidate['generated_index_in_primitive_closure']})
        curves.append({'curve':curve,'candidates':records})
    inputs[str(Path(__file__).relative_to(base.ROOT))]=base.digest(Path(__file__))
    report={'status':'PASS_GENERATED_SUBGROUP_INDEX_AND_REDUCTION_REPLAY','curves':curves,'inputs':inputs,
            'boundary':'Actual subgroup determinants and heights, separate from primitive-closure selection statistics. Numerical 80-digit heights, exact integer basis changes.'}
    path=base.OUT/'common_cover_generated_lattices_v1.json'
    if args.check:
        assert base.read(path)==report
    else:
        base.save(path,report)
    print(report['status'])


if __name__=='__main__':
    main()

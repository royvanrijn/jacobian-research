#!/usr/bin/env python3
"""Independent exact relation and returned-coordinate check for the known28 recovery."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import inventory188_fixed49_point_control as control
from half_lattice_pointed_sieve import linear_combination_python
from search_observability import point_visibility, prepare_chart
from research_runtime.store import checkpoint

ART=control.ART;ROOT=control.ROOT
RELATION=ART/'inventory188_recovered_direction_relation_v1.json'
NEAREST=ART/'inventory188_nearest_translate_visibility_v2.json'
OUT=ART/'inventory188_recovery_explanation_v1.json'


def expected():
    relation=control.engine.cert.read(RELATION);near=control.engine.cert.read(NEAREST)
    data=control.engine.cert.read(control.D/'result.json');cloud=control.engine.cert.read(control.MOD2)
    assert relation['status']=='PASS_EXACT_GROUP_IDENTITY' and relation['denominator']==1
    model=tuple(map(Q,relation['curve']));basis=[tuple(map(Q,p)) for p in relation['basis']]
    found=tuple(map(Q,relation['recovered_point']));word=relation['integer_coefficients']
    assert word[-1]==1 and linear_combination_python(model,basis,word)==found
    assert relation['recovered_point']==cloud['independent_points'][27]
    assert relation['basis'][:27]==cloud['independent_points'][:27]
    first=next(r for r in data['charts'] if r['rank_lower_bound']>27)
    assert first['index']==4
    record={**first['search'],'completed_denominator':125000}
    visibility=point_visibility(record,found)
    assert visibility['status']=='VISIBLE_AND_RECORDED' and visibility['recorded_square_hit']
    maps=control.engine.cert.read(control.D/'maps.json')
    M=prepare_chart(record)[5];N=list(map(Q,maps['rows'][4]['matrix']))
    i=next(i for i,v in enumerate(N) if v)
    assert all(M[j]*N[i]==N[j]*M[i] for j in range(4)), 'actual chart must agree projectively with prior map'
    all_observations=[]
    for r in data['charts']:
        c={**r['search'],'completed_denominator':125000}
        for sign in (-1,1):
            all_observations.append({'chart':r['index'],'sign':sign,
                'visibility':point_visibility(c,(found[0],sign*found[1]))})
    assert all(v['visibility']['status']!='VISIBLE_NOT_RECORDED' for v in all_observations)
    H=near['rows'][0]['rounded_full_gram']
    norm=sum(word[i]*H[i][j]*word[j] for i in range(28) for j in range(28))
    nearest=near['rows'][0]['representatives'][1]['old_word']+[1]
    near_norm=sum(nearest[i]*H[i][j]*nearest[j] for i in range(28) for j in range(28))
    paths=[RELATION,NEAREST,control.D/'result.json',control.D/'maps.json',control.MOD2,Path(__file__).resolve()]
    return {'schema':'elliptic-curves.inventory188-recovery-explanation.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):control.engine.cert.hashed(p) for p in paths},
        'first_gain_chart_one_based':5,'recovered_point':relation['recovered_point'],
        'recovered_equals_public26_plus_old_word':word[:27],
        'first_gain_visibility':visibility,'all_signed_recovered_observations':all_observations,
        'rounded_metric_scale':1000000,'recovered_point_metric_norm':norm,'nearest_translate_metric_norm':near_norm,
        'same_added_subgroup_as_old27_plus_public26':True,
        'boundary':'The returned representative generates the same added subgroup as public26. This is a successful known28 detector control, not a new curve or29th direction. The lower-metric representative was outside the boxes; minimizing canonical height alone need not minimize the searched coordinate.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args();r=expected()
    if args.check:assert r==control.engine.cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve recovery explanation')
        checkpoint(OUT,r)
    print('PASS chart5:',r['first_gain_visibility'],'rounded heights',r['recovered_point_metric_norm'],r['nearest_translate_metric_norm'])

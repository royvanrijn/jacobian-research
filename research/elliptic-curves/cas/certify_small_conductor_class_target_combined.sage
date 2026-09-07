#!/usr/bin/env sage-python
"""Combine audited norm relations with the known points' principal parity relations."""
import argparse
from pathlib import Path
import runpy
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint

ROOT,ART,cert=batch.ROOT,batch.ART,batch.cert
SOURCE=Path(__file__).resolve()
LOWER_SOURCE=ROOT/'elliptic-curves/cas/certify_small_conductor_class_lower16.sage'
LOWER=ART/'small_conductor_class_lower16_v1.json'


def expected(family,wave):
    suffix={'box':'','strip':'_strips','protected':'_protected'}[family]
    source=ROOT/('elliptic-curves/cas/pursue_small_conductor_class_target'+suffix+'.sage')
    runner=runpy.run_path(str(source))
    runner['audit_all'](wave,write=False)
    lower=runpy.run_path(str(LOWER_SOURCE))['expected']()
    if lower!=cert.read(LOWER):raise ArithmeticError('point parity proof differs')
    matrix=runner['prior_state'](wave+1)
    before=matrix.report()
    lookup={(str(c['p']),c['hnf']):i for i,c in enumerate(matrix.base['columns'])}
    relations=[]
    for j in range(22):
        factors=[[lookup[(r['p'],r['hnf'])],1] for r in lower['bad_prime_valuations'] if r['valuations'][j]%2]
        # All omitted valuations are proved even in the lower-bound certificate.
        matrix.add(factors)
        relations.append({'point_index':j,'odd_valuation_columns':[c for c,e in factors]})
    after=matrix.report();dim=matrix.dimension()
    if dim<16:raise ArithmeticError('upper bound contradicts proved class lower bound16')
    return {'schema':'elliptic-curves.small-conductor-class-target-combined.v1','status':'PASS',
        'family':family,'wave':wave,
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,source,runner['audit_path'](wave),LOWER_SOURCE,LOWER]},
        'point_principal_parity_relations':relations,'before_point_relations':before,'matrix':after,
        'independent_gain_from_point_relations':before['quotient_dimension']-dim,
        'unconditional_class_two_rank_lower_bound':16,'conditional_on_grh_class_two_rank_upper_bound':dim,
        'unconditional_curve_rank_lower_bound':22,'conditional_on_grh_curve_rank_upper_bound':2*((dim+7)//2),
        'conditional_on_grh_exact_class_two_rank':16 if dim==16 else None,
        'conditional_on_grh_exact_curve_rank':22 if dim==16 else None,
        'argument':'Each beta_i=4*x(P_i)-theta generates a principal ideal. Its valuations are exactly checked at bad primes and proved even at all others. Hence its listed odd-valuation ideal product is zero in Cl(K)/2, even without factoring the omitted square ideal. Append these valid small-base parity relations to the audited norm matrix and independently check the supported-intersection rank.',
        'claim_boundary':'Known points supply additional principal-ideal relations, not new rational points. The class lower bound16 is unconditional. Generation and the class/curve upper bounds remain GRH-conditional.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=['box','strip','protected'],required=True);p.add_argument('--wave',type=int,required=True);p.add_argument('--check',action='store_true');a=p.parse_args()
    out=ART/('small_conductor_class_target_combined_%s_%03d_v1.json'%(a.family,a.wave))
    result=expected(a.family,a.wave)
    if a.check:
        if cert.read(out)!=result:raise ArithmeticError('combined certificate differs')
    else:
        if out.exists():raise FileExistsError('preserve certificate')
        checkpoint(out,result)
    print('COMBINED CLASS BOUND',result['conditional_on_grh_class_two_rank_upper_bound'],'POINT-RELATION GAIN',result['independent_gain_from_point_relations'],'PASS',flush=True)

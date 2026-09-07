#!/usr/bin/env sage-python
"""Test the <=9 ternary lifts of each fixed MW16 input's mod-3 kernel.

Only exact rational zero sums prove relations. The resulting subgroup rank
does not bound the full elliptic-curve rank.
"""
import argparse
from itertools import product
from pathlib import Path
import sys
from sage.all import EllipticCurve,QQ,GF,Matrix

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
INPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_input_audit_v1.json'
ODD=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_odd_independence_v1.json'
OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_short_relations_v1.json'


def build():
    original={r['fibration_id']:r for r in cert.read(INPUT)['rows']};rows=[]
    for odd in cert.read(ODD)['rows']:
        source=original[odd['fibration_id']];attempt=odd['attempts'][0]
        if attempt['modulus']!=3:raise ArithmeticError('expected fixed mod3 kernel')
        matrix=Matrix(GF(3),[row for s in attempt['signatures'] for row in s['rows']])
        kernel=matrix.right_kernel();basis=kernel.basis()
        if kernel.dimension()>2:raise ArithmeticError('kernel exceeds declared nine-vector gate')
        E=EllipticCurve(QQ,list(map(QQ,source['curve'])));points=[E(list(map(QQ,p))) for p in source['points']]
        tested=[];relations=[];seen=set()
        for coefficients in product(range(3),repeat=kernel.dimension()):
            if not any(coefficients):continue
            v=sum((c*b for c,b in zip(coefficients,basis)),kernel.zero())
            word=tuple(int(x) if int(x)<2 else -1 for x in v)
            if next(x for x in word if x)<0:word=tuple(-x for x in word)
            if word in seen:continue
            seen.add(word)
            P=sum((c*p for c,p in zip(word,points) if c),E(0))
            relation=bool(P.is_zero())
            tested.append({'word':list(word),'exact_zero':relation})
            if relation:relations.append(list(word))
        relation_rank=int(Matrix(QQ,relations).rank()) if relations else 0
        subset=source['independent_subset_indices'];proof=source['independent_subset_certificate']
        model=tuple(map(cert.F,source['curve']));pts=tuple(tuple(map(cert.F,source['points'][i])) for i in subset)
        replay=cert.checked_rank(model,pts,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if cert.json.loads(cert.json.dumps(replay))!=proof:raise ArithmeticError('independent subset replay differs')
        lower=len(subset);upper=16-relation_rank
        if lower>upper:raise ArithmeticError('relations contradict independent subset')
        row={'fibration_id':source['fibration_id'],'parameter':'1','mod3_kernel_dimension':int(kernel.dimension()),
            'tested_ternary_lifts_up_to_sign':tested,'exact_relation_matrix_rank':relation_rank,
            'displayed_subgroup_rank_interval':[lower,upper],
            'status':'PASS_EXACT_DISPLAYED_SUBGROUP_RANK' if lower==upper else 'BOUNDED_RELATION_AUDIT_INCONCLUSIVE'}
        rows.append(row);print('MW16 EXACT RELATIONS',source['fibration_id'],relation_rank,'subgroup interval',lower,upper,flush=True)
    return {'schema':'elliptic-curves.compact-mw16-short-relations.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (INPUT,ODD,Path(__file__).resolve(),Path(cert.__file__).resolve())},
        'rows':rows,'scope':'Exact short relations among the sixteen displayed generic specializations at two fixed t=1 fibres. The independent subsets and relation matrix bound only this displayed subgroup. No upper bound on the full curve rank or generic family rank.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--check',type=Path)
    a=p.parse_args()
    if not a.check and a.output.exists():raise FileExistsError('preserve exact relation audit')
    data=build()
    if a.check:
        if cert.read(a.check)!=data:raise ArithmeticError('exact relation replay differs')
    else:checkpoint(a.output,data)

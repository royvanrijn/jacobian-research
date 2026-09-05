#!/usr/bin/env python3
"""Complete the fixed five-family t=1 audit, retaining inconclusive ranks.

This replaces the first audit's all-or-nothing export, not its frozen run.
The parameter and finite-prime budget are unchanged; no refill is allowed.
"""
import argparse
from dataclasses import asdict
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from mod2_reduction_independence import gf2_rank

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_input_audit_v1.json'


def sources():
    paths=(spec.ATLAS,Path(__file__).resolve(),Path(spec.__file__).resolve(),Path(cert.__file__).resolve(),
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def independent_columns(signatures):
    matrix=[row for signature in signatures for row in signature.rows]
    selected=[]
    for j in range(16):
        trial=selected+[j]
        if gf2_rank(([r[i] for i in trial] for r in matrix),len(trial))>len(selected):
            selected=trial
    return selected


def record(family,stored=None):
    model,points=spec.specialize(family,'1')
    if stored is None:
        signatures=cert.find_mod2_reduction_certificate(model,points,prime_bound=1000)
    else:
        signatures=tuple(cert.mod2_reduction_signature(model,points,s['prime']) for s in stored['signatures'])
    columns=independent_columns(signatures)
    primes=[s.prime for s in signatures]
    proof=cert.checked_rank(model,tuple(points[i] for i in columns),primes,
                            None if stored is None else stored['independent_subset_certificate']['no_rational_2_torsion_prime'])
    return {'fibration_id':family['fibration_id'],'presentation_id':family['presentation_id'],
        'parameter':'1','curve':list(map(str,model)),'points':[list(map(str,p)) for p in points],
        'status':'PASS_16_INDEPENDENT_POINTS' if len(columns)==16 else 'FULL_16_INDEPENDENCE_UNKNOWN_FROM_FIXED_TEST',
        'finite_quotient_rank':len(columns),'signatures':[asdict(s) for s in signatures],
        'independent_subset_indices':columns,'independent_subset_certificate':proof,
        'short_coefficient_bits':max(max(abs(q.numerator).bit_length(),q.denominator.bit_length()) for q in model)}


def main(output,check):
    families=cert.read(spec.ATLAS)['families']
    if check:
        data=cert.read(check)
        if data['sources']!=sources() or len(data['rows'])!=len(families):
            raise ArithmeticError('fixed input bindings changed')
        for family,row in zip(families,data['rows']):
            fresh=record(family,row)
            # JSON normalization converts signature tuples to lists.
            if cert.json.loads(cert.json.dumps(fresh))!=row:
                raise ArithmeticError('fixed input replay differs')
            print('REPLAYED MW16 INPUT',row['fibration_id'],row['status'],flush=True)
        return
    if output.exists(): raise FileExistsError('preserve audit')
    data={'schema':'elliptic-curves.compact-five-mw16-input-audit.v1','sources':sources(),'rows':[],
          'scope':'Exactly t=1 on all five families. Sixteen exported sections checked for membership, finite quotient primes through1000, torsion witness through200. An independent subset is certified; full16 independence stays unknown if the fixed test cannot prove it. No extra point search or new-curve claim.'}
    for family in families:
        row=record(family);data['rows'].append(row)
        cert.write(output,data)
        print('AUDITED MW16 INPUT',row['fibration_id'],row['finite_quotient_rank'],row['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--check',type=Path)
    a=p.parse_args();main(a.output,a.check)

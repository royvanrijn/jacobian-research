#!/usr/bin/env sage-python
"""Bounded numerical proposals, then exact relations for all supplied images.

This proves a subgroup rank, not the rank of the whole specialized curve.
No point search, fibre selection, or change to the frozen pilot.
"""
import argparse,json,sys
from fractions import Fraction
from pathlib import Path
from sage.all import QQ,ZZ,RealField,EllipticCurve,matrix,vector,pari,lcm
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import mestre_parent_calibration as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
OUT=batch.ART/'mestre_parent_calibration_input_span_v1.json'

def compute(propose):
    p=batch.protocol();rows=[];saved=None if propose else cert.read(OUT)
    RF=RealField(384);pari.default('realprecision',110)
    bindings={str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__))}
    for index,row in enumerate(p['rows']):
        source=batch.BATCH/row['id']/'seed.json';seed=cert.read(source)
        E=EllipticCurve(QQ,seed['curve']);cloud=[E([QQ(c) for c in P]) for P in seed['divisor_cloud']]
        selected=seed['independent_column_indices'];basis=[cloud[i] for i in selected]
        if len(cloud)!=14 or len(basis)!=11:raise ArithmeticError('fixed14 cloud and11 seed required')
        missing=[i for i in range(14) if i not in selected];relations=[]
        if propose:
            raw=pari(E).ellheightmatrix([list(P.xy()) for P in cloud],precision=384)
            H=matrix(RF,14,14,lambda i,j:RF(str(raw[i,j])))
            G=H.matrix_from_rows_and_columns(selected,selected)
        for k,target in enumerate(missing):
            if propose:
                values=G.solve_right(vector(RF,[H[i,target] for i in selected]))
                coefficients=[QQ(str(Fraction(str(v)).limit_denominator(64))) for v in values]
                den=lcm([c.denominator() for c in coefficients]);word=[ZZ(den*c) for c in coefficients]
                relation={'target_index':target,'target_multiplier':int(den),'basis_coefficients':list(map(int,word))}
            else:relation=saved['rows'][index]['relations'][k]
            if relation['target_index']!=target or not 1<=relation['target_multiplier']<=lcm(range(1,65)):raise ArithmeticError('relation target or bound differs')
            total=sum((ZZ(a)*P for a,P in zip(relation['basis_coefficients'],basis)),E(0))
            if len(relation['basis_coefficients'])!=11 or total!=ZZ(relation['target_multiplier'])*cloud[target]:raise ArithmeticError('bounded rational relation proposal failed exact group check')
            relations.append(relation)
        proof=seed['rank_certificate']
        from memory_rank_certificate import checked_rank
        actual=checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],
          [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if actual['rank_lower_bound']!=11:raise ArithmeticError('independent seed lower bound differs')
        rows.append({'id':row['id'],'basis_indices':selected,'relations':relations,'supplied_image_subgroup_rank':11})
        bindings[str(source.relative_to(ROOT))]=cert.hashed(source)
        print(row['id'],'EXACT SUPPLIED-IMAGE SUBGROUP RANK11',flush=True)
    return {'schema':'elliptic-curves.mestre-parent-calibration-input-span.v1','status':'PASS','rows':rows,
      'sources':bindings,'proposal_precision_bits':384,'coefficient_denominator_bound':64,
      'scope':'Three exact rational group relations per14-image cloud, plus certified11 independent seed points, prove supplied-image subgroup rank exactly11 at each of twelve frozen fibres. Numeric heights only propose relations; replay uses exact group arithmetic. No full-curve or generic-rank upper bound, saturation claim, or elliptic point search.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args()
    if not a.check and OUT.exists():raise FileExistsError('preserve relation certificate')
    result=compute(not a.check)
    if a.check:assert result==cert.read(OUT)
    else:checkpoint(OUT,result)

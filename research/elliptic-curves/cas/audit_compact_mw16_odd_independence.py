#!/usr/bin/env python3
"""Bounded odd-prime audit of the two retained t=1 mod-2 deficiencies.

Try ell=3, then ell=5 only if necessary, at reduction primes through 1000.
No new parameters, points, saturation search, or rank upper bounds.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
import mod_l_reduction_independence as odd
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_input_audit_v1.json'
OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_odd_independence_v1.json'


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in
        (INPUT,Path(__file__).resolve(),Path(odd.__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py')}


def duplicate_relations(points):
    return [{'indices':[j,i],'sign':sign} for i,(x,y) in enumerate(points)
            for j,(xx,yy) in enumerate(points[:i]) for sign in (1,-1) if x==xx and y==sign*yy]


def torsion_witness(model,ell):
    for p in odd._primes_up_to(200):
        if p in (2,ell):continue
        try:signature=odd.mod_l_reduction_signature(model,[],p,ell)
        except ValueError:continue
        if signature.group_order%ell:
            return {'prime':p,'group_order':signature.group_order}
    return None


def build(output):
    if output.exists():raise FileExistsError('preserve bounded audit')
    inputs=[r for r in cert.read(INPUT)['rows'] if r['finite_quotient_rank']<16]
    if [r['fibration_id'] for r in inputs]!=['a1-fibration-02','a1-fibration-05']:
        raise ArithmeticError('fixed two-input roster changed')
    data={'schema':'elliptic-curves.compact-mw16-odd-independence.v1','sources':sources(),
        'protocol':{'quotient_moduli':[3,5],'reduction_prime_bound':1000,'torsion_prime_bound':200,
                    'stop':'Stop a fibre at first full16 certificate with torsion exclusion.',
                    'scope':'Two retained t=1 inputs only; exact duplicate/sign check then bounded finite odd quotients. No new points, parameters, rational-halving search or rank upper bound.'},'rows':[]}
    checkpoint(output,data)
    for source in inputs:
        model=tuple(map(F,source['curve']));points=tuple(tuple(map(F,p)) for p in source['points'])
        row={'fibration_id':source['fibration_id'],'parameter':'1','duplicate_or_sign_relations':duplicate_relations(points),'attempts':[],
             'status':'INDEPENDENCE_UNKNOWN','rank_lower_bound':None}
        for ell in (3,5):
            signatures=odd.find_mod_l_reduction_certificate(model,points,modulus=ell,prime_bound=1000)
            rank=odd.combined_mod_l_rank(signatures,16,ell)
            witness=torsion_witness(model,ell)
            row['attempts'].append({'modulus':ell,'finite_quotient_rank':rank,'signatures':[asdict(s) for s in signatures],
                                    'no_rational_ell_torsion_witness':witness})
            print('MW16 ODD',source['fibration_id'],'ell',ell,'finite rank',rank,flush=True)
            if rank==16 and witness:
                row.update(status='PASS_16_INDEPENDENT_POINTS',rank_lower_bound=16);break
        data['rows'].append(row);checkpoint(output,data)


def check(path):
    data=cert.read(path)
    if data['sources']!=sources():raise ArithmeticError('sources changed')
    inputs={r['fibration_id']:r for r in cert.read(INPUT)['rows']}
    if [r['fibration_id'] for r in data['rows']]!=['a1-fibration-02','a1-fibration-05']:
        raise ArithmeticError('incomplete or changed two-input roster')
    for row in data['rows']:
        source=inputs[row['fibration_id']];model=tuple(map(F,source['curve']));points=tuple(tuple(map(F,p)) for p in source['points'])
        if row['parameter']!='1' or duplicate_relations(points)!=row['duplicate_or_sign_relations']:
            raise ArithmeticError('point relation audit differs')
        success=False
        for index,attempt in enumerate(row['attempts']):
            ell=attempt['modulus']
            if success or index>=2 or ell!=(3,5)[index]:raise ArithmeticError('attempt order differs')
            signatures=[]
            for stored in attempt['signatures']:
                if not 2<stored['prime']<=1000 or stored['prime']==ell:raise ArithmeticError('prime outside fixed bound')
                sig=odd.mod_l_reduction_signature(model,points,stored['prime'],ell)
                if cert.json.loads(cert.json.dumps(asdict(sig)))!=stored:raise ArithmeticError('finite quotient differs')
                signatures.append(sig)
            rank=odd.combined_mod_l_rank(signatures,16,ell)
            if rank!=attempt['finite_quotient_rank']:raise ArithmeticError('finite rank differs')
            witness=attempt['no_rational_ell_torsion_witness']
            if witness:
                p=witness['prime']
                if not 2<p<=200 or p==ell:raise ArithmeticError('invalid torsion witness prime')
                sig=odd.mod_l_reduction_signature(model,[],p,ell)
                if sig.group_order!=witness['group_order'] or sig.group_order%ell==0:raise ArithmeticError('torsion exclusion failed')
            success=rank==16 and witness is not None
        if (row['status']=='PASS_16_INDEPENDENT_POINTS')!=success or row['rank_lower_bound']!=(16 if success else None):
            raise ArithmeticError('independence claim differs')
        print('REPLAYED ODD MW16',row['fibration_id'],row['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--check',type=Path)
    a=p.parse_args();check(a.check) if a.check else build(a.output)

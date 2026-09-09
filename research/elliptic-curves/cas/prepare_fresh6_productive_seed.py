#!/usr/bin/env python3
"""Prepare first-M18 packets and exact winning generic parents; no point search."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from research_runtime.store import checkpoint
from memory_rank_certificate import checked_rank
from future_point_admission import FinitePointAdmission
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
import compact_atlas_specialization as spec

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(folder, certificate, curve_id):
    from sage.all import ZZ, matrix, pari
    folder.mkdir(exist_ok=False)
    bundle=json.loads(certificate.read_text())
    row=next(r for r in bundle['records'] if r['id']==curve_id)
    proof=row['packet'];family_name=row['family'];parameter=row['parameter'];expected_rank=18
    if proof['rank_lower_bound']!=18:raise ArithmeticError('M18 required')
    terminal=ROOT/row['evidence']['terminal.json']['path']
    verified=ROOT/row['evidence']['verified.json']['path']
    chartpath=terminal.parent/f"chart-{row['charts']-1:04d}.json"
    paths=[certificate,spec.ATLAS,terminal,verified,chartpath]
    for name in ('terminal.json','verified.json'):
        if sha(ROOT/row['evidence'][name]['path'])!=row['evidence'][name]['sha256']:raise ArithmeticError('evidence changed')
    t=json.loads(terminal.read_text());v=json.loads(verified.read_text())
    if v['status']!='PASS_INDEPENDENT_FIRST_SEED_REPLAY' or v['terminal_sha256']!=sha(terminal) or t['last_chart_sha256']!=sha(chartpath):raise ArithmeticError('replay seal differs')
    family=next(r for r in json.loads(spec.ATLAS.read_text())['families'] if r['family']==family_name)
    model=tuple(map(F,proof['curve']));basis=tuple(tuple(map(F,p)) for p in proof['points'])
    original,generic=spec.specialize(family,parameter)
    if original!=model or generic!=basis[:17] or len(basis)!=18:raise ArithmeticError('native seed differs')
    old=proof['proof'];fresh=checked_rank(model,basis,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(fresh))!=old:raise ArithmeticError('independence differs')
    protocol={'family':family_name,'parameter':parameter,'initial_rank':18,'generic_rank':17,'prime_bound':1000,'cvp_node_limit':2000000,'point_searches':0,
      'rule':'Exact winning parent of the first certified M18 point; V3 will diversify extensions from this observed productive parent. No higher points enter the seed.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'sources':{str((CAS/n).relative_to(ROOT)):sha(CAS/n) for n in ('prepare_fresh6_productive_seed.py','compact_atlas_specialization.py','memory_rank_certificate.py','visibility_lattice_fast.py','visibility_lattice_v2.py')}}
    checkpoint(folder/'protocol.json',protocol)
    checkpoint(folder/'seed-M18.json',dict(proof,generic_rank=17))
    chart=json.loads(chartpath.read_text());word=chart['mapping']['centre']['representative']
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    initial=json.loads((terminal.parent/'seed.json').read_text())
    state=certified_state(model,basis[:17],initial['proof'])
    search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=chart['mapping']['coordinate_policy'])
    returned=backend.replay(search,chart['mapping'],chart['search'])
    if basis[17] not in returned:raise ArithmeticError('new seed point absent from winning witness')
    mask=sum((int(w)%2)<<j for j,w in enumerate(word))
    gains=[{'arm':'first_M18','chart_index':chart['index'],'mask':mask,'before':17,'after':18,'proof':fresh,'points':proof['points'],'generic_word':word}]
    checkpoint(folder/'gains.json',gains)
    gram = [[F(x) for x in row] for row in family['generic_height_gram']]
    if any((2*x).denominator != 1 for row in gram for x in row):
        raise ArithmeticError('unexpected generic height denominator')
    g = matrix(ZZ, [[int(2*x) for x in row] for row in gram])
    u = matrix(ZZ, pari(g).qflllgram()).transpose()
    if abs(u.det()) != 1:
        raise ArithmeticError('nonunimodular generic LLL')
    inv = u.inverse()
    fast, reference = IntegerExactParity((u*g*u.transpose()).rows()), ExactParity((u*g*u.transpose()).rows())
    masks = sorted({row['mask'] for row in gains})
    seeds = {row['mask']: row for row in gains}
    rows, checks = [], []
    for mask in masks:
        w = matrix(ZZ, 1, 17, seeds[mask]['generic_word'])
        reduced_seed = tuple(map(int, (w*inv).row(0)))
        residue = tuple(x % 2 for x in reduced_seed)
        import numpy as np
        words0, _ = fast.babai(np.asarray([residue], dtype=np.int64))
        reduced_seed = tuple(map(int, words0[0]))
        cert = fast.solve(residue, reduced_seed, protocol['cvp_node_limit'])
        if cert != reference.solve(residue, reduced_seed, protocol['cvp_node_limit']):
            raise ArithmeticError('independent generic CVP replay differs')
        words = []
        for v in cert['minima']:
            z = tuple(map(int, (matrix(ZZ, 1, 17, v)*u).row(0)))
            if next(x for x in z if x) < 0:
                z = tuple(-x for x in z)
            if sum((x % 2) << j for j, x in enumerate(z)) != mask:
                raise ArithmeticError('generic parity transport differs')
            words.append(z)
        rows.append({'mask': mask, 'word': list(min(words)), 'norm': cert['norm']})
        checks.append({'mask': mask, 'reduced_seed': list(reduced_seed), 'proof': cert})
    checkpoint(folder/'generic-cvp-proofs.json', {'gram': [list(map(int, r)) for r in g.rows()],
        'LLL': [list(map(int, r)) for r in u.rows()], 'checks': checks})
    bank = {'status': 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'dimension': 17,
        'generic_gram_scale': 2, 'rows': rows, 'shells': sorted({r['norm'] for r in rows}),
        'protocol_sha256': sha(folder/'protocol.json'), 'gains_sha256': sha(folder/'gains.json'),
        'generic_cvp_sha256': sha(folder/'generic-cvp-proofs.json'),
        'claim_boundary': 'Exact historically productive subset, not the complete generic shell bank.'}
    checkpoint(folder/'anchor-bank.json', bank)
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('preparation binding changed')
    checkpoint(folder/'prepared.json', {'status': 'PASS_NATIVE_R17_SEED_AND_PRODUCTIVE_ANCHORS',
        'point_searches': 0, 'anchor_count': len(rows),
        'files': {p.name: sha(p) for p in sorted(folder.glob('*.json'))}})
    print('R17_SEED_PREPARED', curve_id, expected_rank, len(rows), masks, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder', type=Path, required=True)
    p.add_argument('--certificate', type=Path, default=ART/'fresh6_first_m18_v1/result.json')
    p.add_argument('--curve-id', required=True)
    a = p.parse_args()
    run(a.folder.resolve(), a.certificate.resolve(), a.curve_id)

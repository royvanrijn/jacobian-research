#!/usr/bin/env python3
"""Reference replay of expanded CVPs, reusing only previously sealed exact checks."""
import argparse
import json
from pathlib import Path
import numpy as np
from sage.all import ZZ, matrix
from visibility_lattice_v2 import ExactParity
from research_runtime.store import checkpoint
from run_future_visibility_upper_audit import guard, sha, read, LOCAL


def run(folder):
    output = folder/'reference-cvps.json'
    if output.exists():
        raise FileExistsError('preserve reference replay')
    seal, protocol = read(folder/'sealed.json'), read(folder/'protocol.json')
    guard(protocol)
    path = folder/'landscape/selection.json'
    if seal['selection_sha256'] != sha(path) or seal['protocol_sha256'] != sha(folder/'protocol.json'):
        raise ArithmeticError('expanded selection is unsealed')
    data = read(path)
    prior = LOCAL/'future-v3-curve90-masked-v1'
    checked = read(prior/'evaluation.json')
    guard(read(prior/'protocol.json'))
    if (checked['status'] != 'PASS_REFERENCE_CVP_AND_MASKED_VISIBILITY_AUDIT'
            or checked['selection_sha256'] != sha(prior/'landscape/selection.json')):
        raise ArithmeticError('previous independent CVP checks unavailable')
    old = read(prior/'landscape/selection.json')
    cache = {}
    if old['rounded_gram'] == data['rounded_gram'] and old['LLL'] == data['LLL']:
        for ai, anchor in enumerate(old['anchors']):
            npz = prior/'landscape'/f'anchor-{ai:02d}-full.npz'
            if sha(npz) != anchor['full_scores_sha256']:
                raise ArithmeticError('old checked CVP inputs changed')
            with np.load(npz, allow_pickle=False) as arrays:
                for row in anchor['refined']:
                    j = row['extension']
                    cache[tuple(map(int, arrays['reduced_residues'][j])),
                          tuple(map(int, arrays['babai_words'][j]))] = row['cvp']
    g, u = (matrix(ZZ, data[k]) for k in ('rounded_gram', 'LLL'))
    solver = ExactParity((u*g*u.transpose()).rows())
    result = {'status': 'RUNNING', 'selection_sha256': sha(path), 'reused': 0, 'fresh': 0,
        'previous_evaluation_sha256': sha(prior/'evaluation.json'),
        'source_sha256': sha(Path(__file__)), 'point_searches': 0}
    checkpoint(output, result)
    for ai, anchor in enumerate(data['anchors']):
        npz = folder/'landscape'/f'anchor-{ai:02d}-full.npz'
        if sha(npz) != anchor['full_scores_sha256']:
            raise ArithmeticError('expanded score arrays changed')
        with np.load(npz, allow_pickle=False) as arrays:
            for row in anchor['refined']:
                j = row['extension']
                p, w = (tuple(map(int, arrays[key][j])) for key in ('reduced_residues', 'babai_words'))
                if (p, w) in cache:
                    proof = cache[p, w]
                    result['reused'] += 1
                else:
                    proof = json.loads(json.dumps(solver.solve(p, w, protocol['policy']['exact_cvp_node_limit'])))
                    result['fresh'] += 1
                if proof != row['cvp']:
                    raise ArithmeticError('exact expanded CVP replay differs')
        result['anchors'] = ai+1
        checkpoint(output, result)
        print('UPPER_REFERENCE', ai+1, '/', len(data['anchors']), flush=True)
    guard(protocol)
    if sha(path) != result['selection_sha256'] or sha(prior/'evaluation.json') != result['previous_evaluation_sha256']:
        raise ArithmeticError('replay input changed')
    result['status'] = 'PASS_REFERENCE_EXPANDED_CVPS'
    result['claim_boundary'] = 'Original rational CVP solver or identical previously replayed input; no new search.'
    checkpoint(output, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    run(parser.parse_args().folder)

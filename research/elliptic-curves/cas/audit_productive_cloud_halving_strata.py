#!/usr/bin/env python3
"""Fixed size-stratified eight-point diagnostic for finite-span UNKNOWNs in a completed V3 cloud."""
import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
from sage.all import EllipticCurve, QQ
from split_seed_descent import build_frame, Classifier, decode_point
from future_point_admission import FinitePointAdmission
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]


def run(folder, output):
    read = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    output.mkdir(exist_ok=False)
    terminal, verified = read(folder/'terminal.json'), read(folder/'verified.json')
    assert verified['status'] == 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY'
    assert verified['terminal_sha256'] == sha(folder/'terminal.json')
    paths = [folder/'terminal.json', folder/'verified.json', Path(__file__),
             Path(__file__).with_name('split_seed_descent.py')]
    candidates = set()
    for path in sorted(folder.glob('epoch-*/cloud.json')):
        paths.append(path)
        for chart in read(path)['charts']:
            for p in chart['search']['finite_curve_points']:
                candidates.add((F(p['x']), abs(F(p['y']))))
    model = tuple(map(F, terminal['curve']))
    basis = tuple(tuple(map(F, p)) for p in terminal['points'])
    candidates -= {(x, abs(y)) for x,y in basis}
    height = lambda p: max(abs(v.numerator) for v in p) * max(v.denominator for v in p)
    ordered = sorted(candidates, key=lambda p:(height(p), p))
    chosen = [ordered[i] for i in sorted({j*(len(ordered)-1)//7 for j in range(8)})] if ordered else []
    protocol = {'point_searches': 0, 'sample_size': 8, 'max_steps': 4,
        'selection': 'Eight evenly spaced order statistics (including extremes) by max-absolute-numerator times max-denominator score, '
                     'with coordinate tie-break, among sign-normalized nonbasis cloud points.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'points': [list(map(str,p)) for p in chosen],
        'claim_boundary': 'Fixed diagnostic only. A finite-span UNKNOWN is not dependence. '
                          'Unexpected independence needs a separate replay before promotion.'}
    checkpoint(output/'protocol.json', protocol)
    admission = FinitePointAdmission(model, basis)
    E = EllipticCurve(QQ, list(map(str,model)))
    frame = build_frame(E, [E(list(map(str,p))) for p in basis])
    checkpoint(output/'frame.json',frame)
    classifier = Classifier(frame)
    rows = []
    for p in chosen:
        observation = admission.consider(p)
        result = classifier.classify(E(list(map(str,p))), max_steps=4)
        if result['status'] == 'INHERITED_RATIONAL_SPAN':
            assert int(result['relation_multiplier'])*E(list(map(str,p))) == sum(
                (int(n)*Q for n,Q in zip(result['relation_word'], classifier.basis)), E(0))
        rows.append({'point':list(map(str,p)), 'finite_observation':observation, 'result':result})
        checkpoint(output/'progress.json', {'rows':rows})
        print('HALVING_DIAGNOSTIC',len(rows),result['status'],result['steps'],flush=True)
    assert all(sha(ROOT/n)==h for n,h in protocol['inputs'].items())
    checkpoint(output/'result.json', {'status':'COMPLETE_FIXED_HALVING_DIAGNOSTIC',
        'protocol_sha256':sha(output/'protocol.json'), 'frame_sha256':sha(output/'frame.json'),
        'rows':rows, 'point_searches':0})


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.folder.resolve(),a.output.resolve())

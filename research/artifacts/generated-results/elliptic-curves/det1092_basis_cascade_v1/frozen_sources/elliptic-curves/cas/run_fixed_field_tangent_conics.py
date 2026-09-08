#!/usr/bin/env sage
"""Replay the cubic tangent-conic gate for genuine radical-cover lifts.

Default is an exact offline audit. --collect packages retained experiments;
--prepare extracts their inputs; --lattice replays the eight bounded local
conic reductions. No missing point or resource failure is an obstruction.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

from sage.all import GF, QQ, ZZ, NumberField, PolynomialRing, diagonal_matrix, matrix, vector
from sage.version import version as sage_version
import run_fixed_field_radical_covers as models

ROOT = models.ROOT
WORK = ROOT / 'artifacts/local/fixed-field-genuine-lift'
SUMMARY = ROOT / 'artifacts/generated-results/elliptic-curves/fixed_field_tangent_conics_v1.json'
EVIDENCE = SUMMARY.with_name('fixed_field_tangent_conics_evidence_v1.json.gz')
REMOTE = ('eight-on-target-1047173', 'cassels-conic-scalar',
          'conic-1047173-seed2', 'conic-596921-seed1', 'conic-450876-seed1',
          'conic-596921-quiet', 'conic-450876-quiet', 'best-reduced-conic')


def remote_status(raw):
    root = ET.fromstring(raw)
    text = '\n'.join(''.join(z.itertext()) for z in root.findall('.//results/line'))
    if 'memory limit' in text:
        return 'MEMORY_LIMIT_NO_WITNESS'
    if 'time limit' in raw:
        return 'TIME_LIMIT_NO_WITNESS'
    if root.findall('.//warning') or 'error:' in text.lower():
        return 'INCOMPLETE_NO_WITNESS'
    return 'REQUIRES_MANUAL_WITNESS_REPLAY'


def audit(evidence):
    old_summary = json.loads(models.SUMMARY.read_text())
    assert models.sha(models.EVIDENCE) == old_summary['evidence_sha256']
    assert evidence['models_evidence_sha256'] == models.sha(models.EVIDENCE)
    for path, digest in old_summary['source_hashes'].items():
        assert models.sha(ROOT / path) == digest
    records = json.loads(gzip.decompress(models.EVIDENCE.read_bytes()))['covers']
    _, _, A, B, E = models.base.context(models.base.SOURCE, -1)
    R = PolynomialRing(QQ, 't')
    q = R(evidence['field_polynomial'])
    assert q.degree() == 3 and q.is_monic()
    prime = ZZ(evidence['irreducibility_prime'])
    assert prime.is_prime() and q.change_ring(GF(prime)).is_irreducible()
    K = NumberField(q, 'theta')
    theta = K.gen()
    anchor_theta = 9*theta + 6
    assert anchor_theta**3 + A*anchor_theta + B == 0
    decode = lambda coefficients: K(list(map(QQ, coefficients)))
    assert len(evidence['conics']) == 3
    assert {r['mask'] for r in evidence['conics']} == set(models.approved_masks())
    conics = {}
    for row in evidence['conics']:
        rec = next(r for r in records if r['cover']['mask'] == row['mask']
                   and not r['cover']['translated_by_universal_point'])
        H, _ = models.verify_quadric_model(rec['cover'], rec['quadric_model'], E)
        lam = decode(row['lambda'])
        singular = lam*H[0] + H[1]
        assert singular.det() == 0 and singular.rank() == 3
        indices = row['indices']
        assert len(indices) == len(set(indices)) == 3 and set(indices) <= set(range(4))
        C = matrix(K, 3, [decode(z) for z in row['conic_matrix']])
        assert C == singular.matrix_from_rows_and_columns(indices, indices)
        assert C.det()
        # A nonzero principal minor gives projection from the vertex of
        # this quadric cone onto the displayed nonsingular plane conic.
        missing = next(iter(set(range(4))-set(indices)))
        kernel = singular.right_kernel().basis()[0]
        assert kernel[missing]
        conics[row['mask']] = C
    C = conics[1047173]
    start = evidence['norm_conic_start']
    a0, b0, d = decode(start['a']), decode(start['b']), QQ(start['ordinate_scale'])
    assert a0 == -C[1, 1]/2 and not C[0, 1] and not C[0, 2] and C[0, 0] == 2
    assert b0 == (-C[2, 2]/2-C[1, 2]**2/(4*a0))*d*d
    beta = -b0/a0
    base_map = matrix(K, [[0, 0, a0], [1, C[1, 2]*d/(2*a0), 0], [0, d, 0]])
    assert base_map.det()
    assert base_map.transpose()*C*base_map == 2*a0*diagonal_matrix(K, [-1, beta, a0])
    assert len(evidence['local_runs']) == 8
    assert {r['seed'] for r in evidence['local_runs']} == set(range(1, 9))
    norms = []
    points = []
    for trial in evidence['local_runs']:
        checkpoint = trial['checkpoint']
        result = trial['result']
        assert result['seed'] == trial['seed'] and result['target_mask'] == 1047173
        assert result['point_or_sha'] == 'UNKNOWN' and not result['failures']
        a, b = decode(checkpoint['a']), decode(checkpoint['b'])
        T = matrix(K, 3, [decode(z) for z in checkpoint['map_to_initial_norm_conic']])
        assert T.det()
        composed = base_map*T
        transformed = composed.transpose()*C*composed
        multiplier = transformed[0, 0]
        assert multiplier and transformed == multiplier*diagonal_matrix(K, [1, -a, -b])
        assert result['iterations'] == len(checkpoint['steps']) == len(result['steps'])
        if result['point_on_auxiliary_conic'] is not None:
            p = vector(K, [decode(z) for z in result['point_on_auxiliary_conic']])
            assert any(p) and p*C*p == 0
            points.append(trial['seed'])
        norms.append((max(abs(a.norm()), abs(b.norm())), trial['seed'], abs(a.norm()), abs(b.norm())))
    assert len(evidence['remote_attempts']) == len(REMOTE)
    assert {r['name'] for r in evidence['remote_attempts']} == set(REMOTE)
    for attempt in evidence['remote_attempts']:
        assert attempt['status'] == remote_status(attempt['xml'])
        assert attempt['status'] != 'REQUIRES_MANUAL_WITNESS_REPLAY', \
            'a completed CAS job requires examination before publication'
    probe = evidence['indefinite_probe']['result']
    assert probe['target_mask'] == 1047173 and probe['seed'] == 1
    assert probe['point_or_sha'] == 'UNKNOWN' and not probe['failures']
    assert probe['point_on_auxiliary_conic'] is None, \
        'a probe witness requires exact map replay before publication'
    best = min(norms)
    return {'exact_cubic_tangent_conics': 3, 'exact_reduced_conic_maps': 8,
            'verified_auxiliary_conic_points': len(points), 'genuine_higher_covers': 0,
            'best_trial_seed': best[1], 'best_coefficient_absolute_norms': list(map(str, best[2:])),
            'curve_rank_lower_bound': 1,
            'point_or_sha': {str(m): 'UNKNOWN' for m in models.MASKS}}


def collect(work):
    conics = [json.loads((work/f'conic-{m}-seed{s}.json').read_text())
              for m, s in [(1047173, 2), (596921, 1), (450876, 1)]]
    evidence = {'models_evidence_sha256': models.sha(models.EVIDENCE),
                'field_polynomial': conics[0]['field_polynomial'],
                'irreducibility_prime': 0, 'conics': conics,
                'norm_conic_start': json.loads((work/'lagrange_start.json').read_text()),
                'local_runs': [], 'remote_attempts': [], 'replay_inputs': {}}
    from sage.all import GF, prime_range
    q = PolynomialRing(QQ, 't')(evidence['field_polynomial'])
    evidence['irreducibility_prime'] = int(next(p for p in prime_range(100) if q.change_ring(GF(p)).is_irreducible()))
    for seed in range(1, 9):
        evidence['local_runs'].append({'seed': seed,
            'checkpoint': json.loads((work/f'conic-multistart-{seed}-checkpoint.json').read_text()),
            'result': json.loads((work/f'conic-multistart-{seed}-result.json').read_text()),
            'log': (work/f'conic-multistart-{seed}.log').read_text()})
    for name in REMOTE:
        raw = (work/f'{name}.xml').read_text()
        evidence['remote_attempts'].append({'name': name, 'input': (work/f'{name}.m').read_text(),
                                            'xml': raw, 'status': remote_status(raw)})
    evidence['indefinite_probe'] = {
        'result': json.loads((work/'conic-indefinite-1-result.json').read_text()),
        'log': (work/'conic-indefinite-1.log').read_text()}
    for name in ('conic_multistart.py', 'conic_indefinite.py',
                 'cassels-conic.json', 'lagrange_start.json'):
        evidence['replay_inputs'][name] = (work/name).read_text()
    arithmetic = audit(evidence)
    EVIDENCE.write_bytes(gzip.compress((json.dumps(evidence, sort_keys=True)+'\n').encode(), mtime=0))
    models.base.save(SUMMARY, {'schema': 'elliptic-curves.fixed-field-tangent-conics.v1',
        'status': 'GENUINE_LIFT_NOT_YET_CONSTRUCTED', 'arithmetic': arithmetic,
        'checker_sha256': models.sha(__file__), 'evidence_sha256': models.sha(EVIDENCE),
        'software': {'sage': sage_version, 'magma': '2.29-9'},
        'limits': {'remote_seconds_per_job': 60, 'local_starts': 8,
                   'local_seconds_per_start': 45, 'local_inner_seconds': 40,
                   'local_max_iterations': 18, 'local_lattice_variants_per_iteration': 64,
                   'indefinite_probe_inner_seconds': 110},
        'claim_boundary': ['The tangent conics are known soluble by the Hasse principle; no point was found.',
                           'Conic reductions and resource failures do not decide a target torsor.',
                           'No calibrated point search on a genuine higher cover was possible.']})
    print(json.dumps(arithmetic))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, default=WORK)
    parser.add_argument('--collect', action='store_true')
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--lattice', action='store_true')
    args = parser.parse_args()
    if args.collect:
        collect(args.workdir)
        return
    summary = json.loads(SUMMARY.read_text())
    assert models.sha(__file__) == summary['checker_sha256']
    assert models.sha(EVIDENCE) == summary['evidence_sha256']
    evidence = json.loads(gzip.decompress(EVIDENCE.read_bytes()))
    assert audit(evidence) == summary['arithmetic']
    if args.prepare:
        args.workdir.mkdir(parents=True, exist_ok=True)
        for name, content in evidence['replay_inputs'].items():
            if name.endswith('.py'):
                content = content.replace("D=Path('artifacts/local/fixed-field-genuine-lift')",
                                          'D=Path(__file__).resolve().parent')
            (args.workdir/name).write_text(content)
        for attempt in evidence['remote_attempts']:
            (args.workdir/(attempt['name']+'.m')).write_text(attempt['input'])
    if args.lattice:
        # In Sage 10.9 ideal.element_1_mod unnecessarily requests pari_bnf.
        # The retained solver calls nf.idealaddtoone directly. Enforce that
        # no Sage class-group request slips into this positive-witness path.
        wrapper = ('import runpy,sys\n'
                   'from sage.rings.number_field.number_field import NumberField_generic\n'
                   'def forbidden(*args,**kwargs):\n'
                   ' raise RuntimeError("Class-group call forbidden in conic lattice replay")\n'
                   'NumberField_generic.pari_bnf=forbidden\n'
                   'sys.argv=sys.argv[1:]\n'
                   'runpy.run_path(sys.argv[0],run_name="__main__")\n')
        for seed in range(1, 9):
            with (args.workdir/f'conic-multistart-{seed}.log').open('w') as log:
                subprocess.run([sys.executable, '-u', '-c', wrapper,
                                str(args.workdir/'conic_multistart.py'), str(seed)],
                               cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=45, check=True)
            trial = next(r for r in evidence['local_runs'] if r['seed'] == seed)
            actual = json.loads((args.workdir/f'conic-multistart-{seed}-checkpoint.json').read_text())
            assert actual == trial['checkpoint'], 'replayed exact checkpoint changed'
            actual = json.loads((args.workdir/f'conic-multistart-{seed}-result.json').read_text())
            # Wall-clock time varies; every arithmetic result must match.
            actual.pop('seconds')
            expected = {k: v for k, v in trial['result'].items() if k != 'seconds'}
            assert actual == expected, 'replayed arithmetic result changed'
    print('PASS_THREE_TANGENT_CONICS_EIGHT_EXACT_REDUCTIONS; NO_GENUINE_LIFT; ALL_TARGETS_UNKNOWN')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Freeze/replay only the four corrected, certified-M17 split admissions.

freeze uses Python; case uses Sage. No campaign or V3 worker is launched here.
Receipts are immutable; execution failures must retain their supervisor logs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2'
PILOT = ROOT/'artifacts/local/elliptic-curves/euclidean-seed-foundry-v1'
ROSTER = [(4, '2', 82931), (17, '-4/3', 65035), (43, '-5/7', 82931), (46, '2/7', 30223)]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(ok, message):
    if not ok:
        raise ArithmeticError(message)


def save(path, value):
    value = json.loads(json.dumps(value))
    if path.exists():
        check(read(path) == value, 'immutable receipt differs: '+str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True)+'\n')


def freeze():
    root = PILOT/'runtime/research'
    plan, manifest = read(PILOT/'plan.json'), read(PILOT/'manifest.json')
    check(sha(PILOT/'plan.json') == manifest['plan_sha256'], 'original plan seal')
    inputs = {}

    def retain(source, name):
        dest = OUT/'inputs'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copyfile(source, dest)
        check(sha(source) == sha(dest), 'frozen input changed')
        inputs[str(dest.relative_to(OUT))] = dict(sha256=sha(dest), origin=str(source.relative_to(ROOT)))

    retain(PILOT/'plan.json', 'plan.json')
    retain(PILOT/'manifest.json', 'manifest.json')
    retain(root/'seed-inputs/parent.json', 'parent.json')
    check(sha(root/'seed-inputs/parent.json') == manifest['files']['seed-inputs/parent.json'], 'parent manifest seal')
    for index, parameter, mask in ROSTER:
        trial = root/'seed-trials'/str(index)
        job = trial/'admit'
        request, old, seed = map(read, [job/'request.json', job/'admission.json', job/'seed.json'])
        check(request['parameter'] == parameter and request['candidates'][0]['mask'] == mask,
              'corrected four-case roster')
        check(len(request['candidates']) == len(old['candidates']) == 1, 'one frozen split')
        check(old['candidates'][0]['status'] == 'UNKNOWN_FINITE_COLUMN_IN_SPAN', 'original UNKNOWN')
        check(seed['rank_lower_bound'] == seed['generic_rank'] == len(seed['points']) == 17, 'sealed M17 seed')
        check(seed == read(job/'generic.json'), 'no subsequent points in seed')
        check(old['seed_sha256'] == sha(job/'seed.json') and old['request_sha256'] == sha(job/'request.json'),
              'original admission seals')
        check(request['incidence_sha256'] == sha(trial/'incidence.json'), 'incidence seal')
        supervision = read(job/'supervision-receipt.json')
        check(supervision['outcome'] == 'completed' and supervision['returncode'] == 0, 'original worker completed')
        for name, digest in supervision['output_sha256'].items():
            check(sha(job/name) == digest, 'original supervisor output seal')
        for name in ('request.json', 'admission.json', 'seed.json', 'seed-gate.json', 'supervision-receipt.json'):
            retain(job/name, f'{index}/{name}')
        retain(trial/'incidence.json', f'{index}/incidence.json')
        for name in ('request.json', 'result.json', 'supervision-receipt.json'):
            retain(root/'seed-traces'/str(mask)/name, f'traces/{mask}/{name}')
    save(OUT/'protocol.json', dict(schema='euclidean-four-split-admission.v2',
        roster=[dict(index=i, parameter=t, mask=m) for i, t, m in ROSTER], inputs=inputs,
        max_successful_halves=8, classification_wall_seconds=120,
        independent_replay_wall_seconds=180, rss_bytes=1610612736, maximum_workers=2,
        footprint='Exactly all places in each archived M17 seed proof, in saved order; no extra places.',
        nonhalving_proof='Primitive duplication quartic with a root-free projective reduction at a sealed footprint prime, or exact rational factorization replay.',
        boundary='Corrected user roster: four certified-M17 split-point UNKNOWNs only. Five generic-gate failures excluded. No new parameters/orbits, later V3 points, point search, heights or rank upper bounds.',
        statuses=['INHERITED_RATIONAL_SPAN', 'NEW_INDEPENDENT_DIRECTION', 'UNKNOWN']))
    print('FROZEN_FOUR_CERTIFIED_M17_UNKNOWNS', flush=True)


def case(index):
    from sage.all import QQ
    from sage.env import SAGE_VERSION
    from prospective_split_admission import ProspectiveAdmission, frame_from_packet
    protocol = read(OUT/'protocol.json')
    for name, info in protocol['inputs'].items():
        check(sha(OUT/name) == info['sha256'], 'sealed input differs')
    row = next(r for r in protocol['roster'] if r['index'] == index)
    folder = OUT/'inputs'/str(index)
    seed, request, old = map(read, [folder/'seed.json', folder/'request.json', folder/'admission.json'])
    frame = frame_from_packet(seed)
    save(OUT/'cases'/str(index)/'frame.json', frame)
    engine = ProspectiveAdmission(frame)
    t = QQ(row['parameter']); d = t.denominator()
    unscaled = list(map(QQ, request['candidates'][0]['point']))
    point = list(map(str, [unscaled[0]*d**4, unscaled[1]*d**6]))
    check(point == old['candidates'][0]['point'], 'exact short-model transport')
    result = engine.consider(point, max_steps=protocol['max_successful_halves'])
    sources = ['prospective_split_admission.py', 'split_seed_descent.py',
               'verify_det1092_funnel_small_conic_seed.sage', Path(__file__).name]
    result.update(schema='euclidean-exact-case-receipt.v2', **row,
        protocol_sha256=sha(OUT/'protocol.json'), frame_sha256=sha(OUT/'cases'/str(index)/'frame.json'),
        seed_sha256=sha(folder/'seed.json'), sage_version=SAGE_VERSION,
        source_hashes={str((CAS/name).relative_to(ROOT)): sha(CAS/name) for name in sources})
    save(OUT/'cases'/str(index)/'receipt.json', result)
    print(index, row['parameter'], row['mask'], result['status'], result.get('reason'), result['steps'], flush=True)


def run():
    from concurrent.futures import ThreadPoolExecutor
    from research_runtime.supervisor import Limits, run as supervise
    protocol = read(OUT/'protocol.json')
    sage = shutil.which('sage')
    check(sage is not None, 'Sage required')

    def attempt(row):
        index = row['index']
        folder = OUT/'execution'/str(index)
        folder.mkdir(parents=True, exist_ok=True)
        command = [sage, '-python', str(Path(__file__).resolve()), 'case', '--index', str(index)]
        record = supervise(command, limits=Limits(protocol['classification_wall_seconds'], protocol['rss_bytes']),
                           log_path=folder/'classification.log', cwd=ROOT)
        save(folder/'classification.json', record)
        print(index, record['outcome'], record['returncode'], flush=True)
        check(record['outcome'] == 'completed' and record['returncode'] == 0, 'classification worker failed; retain receipt')

    with ThreadPoolExecutor(max_workers=protocol['maximum_workers']) as pool:
        list(pool.map(attempt, protocol['roster']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('freeze', 'run', 'case'))
    parser.add_argument('--index', type=int, choices=[r[0] for r in ROSTER])
    args = parser.parse_args()
    if args.mode == 'freeze':
        freeze()
    elif args.mode == 'run':
        run()
    else:
        check(args.index is not None, 'case index required')
        case(args.index)

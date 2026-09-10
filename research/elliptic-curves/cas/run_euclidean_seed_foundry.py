#!/usr/bin/env python3
"""Prospective norm-ten seed sieve -> certified admission -> unchanged V3.

Opt-in, bounded, resumable campaign. Never edits or stops a running foundry.
Python prepares/plans; Sage is required only for the subprocess workers.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from euclidean_seed_sieve import address, conic_from_trace, evaluate, incidence, ordered_orbits

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DEFAULT = ROOT/'artifacts/local/elliptic-curves/euclidean-seed-foundry-v1'
SUFFIXES = {'.py', '.sage', '.gp', '.c', '.cpp', '.h', '.sh'}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    """Immutable receipts. All writers hold the campaign lease or unique jobs."""
    path = Path(path)
    value = json.loads(json.dumps(value))
    if path.exists():
        if read(path) != value:
            raise ArithmeticError('immutable receipt differs: '+str(path))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    tmp.replace(path)


def normalized_parent(reduced, gram):
    def polynomial(record):
        den = list(map(F, record['denominator']))
        if not den or not den[0] or any(den[1:]):
            raise ValueError('parent invariant is not polynomial')
        return [str(F(c)/den[0]) for c in record['numerator']]
    invariants = [polynomial(a) for a in reduced['a_invariants']]
    if any(any(map(F, a)) for a in invariants[:3]):
        raise ValueError('reduced det1092 parent must be short')
    def rational(record):
        return {'numerator_coefficients_low_to_high': list(map(str, record['numerator'])),
                'denominator_coefficients_low_to_high': list(map(str, record['denominator']))}
    sections = [{'basis_index': i, 'X': rational(p[0]), 'Y': rational(p[1])}
                for i, p in enumerate(reduced['basis_weierstrass_coordinates'])]
    if len(sections) != 17 or len(gram) != 17 or any(len(r) != 17 for r in gram):
        raise ValueError('expected the pinned seventeen-section parent')
    return {'family': 'det1092', 'generic_rank_lower_bound': 17,
            'A_coefficients_low_to_high': invariants[3], 'B_coefficients_low_to_high': invariants[4],
            'sections': sections, 'generic_height_gram': gram}


def prepare(args):
    if not (1 <= args.workers <= 4 and args.parameters > 0 and args.parameter_offset >= 0
            and args.control_every > 0 and args.trace_seconds > 0):
        raise ValueError('invalid bounded campaign settings')
    sage = shutil.which(args.sage)
    if sage is None or not Path('/usr/bin/gp').exists():
        raise RuntimeError('Sage and /usr/bin/gp are required; no campaign was created')
    reduced_path = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    original_path = ART/'curve302_recovered_mw17_parent_v1.json'
    table_path = ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
    parent = normalized_parent(read(reduced_path), read(original_path)['generic_height_gram'])
    with table_path.open() as stream:
        orbits = ordered_orbits(list(csv.DictReader(stream, delimiter='\t')), args.orbit_offset, args.orbits)
    gram = [[F(x) for x in r] for r in parent['generic_height_gram']]
    for r in orbits:
        w = r['word']
        if len(w) != 17 or sum(w[i]*gram[i][j]*w[j] for i in range(17) for j in range(17)) != 10:
            raise ArithmeticError('orbit word fails exact generic norm-ten gate')
    folder = args.folder.resolve()
    folder.mkdir(parents=True, exist_ok=False)
    root = folder/'runtime/research'
    # The workers see only frozen code and a generic-only projection of inputs.
    for directory in (CAS, ROOT/'elliptic-curves/ecsearch', ROOT/'elkies-k3/scripts'):
        for p in directory.rglob('*'):
            if p.is_file() and p.suffix in SUFFIXES and '__pycache__' not in p.parts:
                dest = root/p.relative_to(ROOT)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(p, dest)
    for p in (ROOT/'elliptic-curves').glob('*.py'):
        dest = root/p.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, dest)
    save(root/'seed-inputs/parent.json', parent)
    plan = {'schema': 'euclidean-seed-foundry.v1', 'root': str(root), 'sage': sage,
            'workers': args.workers, 'trace_seconds': args.trace_seconds,
            'admission_seconds': 900, 'cascade_seconds': 7200, 'rss_bytes': 3*1024**3,
            'point_calls': 64, 'height': 125000, 'control_every': args.control_every,
            'orbits': orbits, 'orbit_offset': args.orbit_offset,
            'parameters': [{'index': i, 'parameter': address(i), 'control': i % args.control_every == 0}
                           for i in range(args.parameter_offset, args.parameter_offset+args.parameters)],
            'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'input_hashes': {str(p.relative_to(ROOT)): sha(p) for p in (reduced_path, original_path, table_path)},
            'claim_boundary': 'A bounded prospective sieve, not an atlas closure, exact-rank claim, or demonstrated speedup. '
                              'Finite-column failures stay UNKNOWN. Seed and control costs are separate.'}
    save(folder/'plan.json', plan)
    save(folder/'manifest.json', {'plan_sha256': sha(folder/'plan.json'),
        'files': {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*')) if p.is_file()},
        'executables': {sage: sha(sage), '/usr/bin/gp': sha('/usr/bin/gp')}})
    print(json.dumps({'status': 'PREPARED_NOT_LAUNCHED', 'folder': str(folder),
                      'orbits': len(orbits), 'parameters': len(plan['parameters'])}))


def guard(folder):
    manifest, plan = read(folder/'manifest.json'), read(folder/'plan.json')
    if sha(folder/'plan.json') != manifest['plan_sha256']:
        raise ArithmeticError('frozen plan changed')
    root = Path(plan['root'])
    for name, digest in manifest['files'].items():
        if sha(root/name) != digest:
            raise ArithmeticError('frozen source/input changed: '+name)
    for name, digest in manifest['executables'].items():
        if sha(name) != digest:
            raise ArithmeticError('CAS executable changed: '+name)
    return plan, root


def build_trace(job):
    """Cold generic trace, checked again by reversed manual group arithmetic."""
    from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve
    request = read(job/'request.json')
    parent = read(ROOT/request['parent'])
    if sha(ROOT/request['parent']) != request['parent_sha256']:
        raise ArithmeticError('parent input changed')
    R = PolynomialRing(QQ, 't'); K = R.fraction_field()
    A, B = (R(parent[k]) for k in ('A_coefficients_low_to_high', 'B_coefficients_low_to_high'))
    def decode(r):
        return K(R(r['numerator_coefficients_low_to_high']))/R(r['denominator_coefficients_low_to_high'])
    E = EllipticCurve(K, [A, B])
    basis = [E([decode(s['X']), decode(s['Y'])]) for s in parent['sections']]
    w = request['word']
    P = -sum((ZZ(n)*Q for n, Q in zip(w, basis)), E(0))
    if P.is_zero():
        raise ArithmeticError('zero trace')
    def plus(p, q):
        if p is None: return q
        if q is None: return p
        x, y = p; u, v = q
        if x == u and y == -v: return None
        slope = (v-y)/(u-x) if x != u else (3*x*x+A)/(2*y)
        z = slope*slope-x-u
        return z, slope*(x-z)-y
    check = None
    for n, Q in reversed(list(zip(w, basis))):
        term = (Q[0], -Q[1] if n > 0 else Q[1])
        multiplier = abs(n)
        while multiplier:
            if multiplier & 1: check = plus(check, term)
            term = plus(term, term); multiplier >>= 1
    if check != (P[0], P[1]):
        raise ArithmeticError('independent rational-function trace replay failed')
    den = P[0].denominator().monic()
    if den.degree() != 6 or not den.is_square():
        save(job/'result.json', {'status': 'UNRESOLVED_TRACE_CHART', 'mask': request['mask'],
                                'reason': 'degree-three finite pole chart unavailable'})
        return
    h = R(den.sqrt()).monic()
    nx, ny = R(P[0]*h**2), R(P[1]*h**3)
    conic = conic_from_trace(*(list(map(str, f)) for f in (A, B, h, nx, ny)))
    save(job/'result.json', {'status': 'PASS_EXACT_CONIC_MAP', 'mask': request['mask'], 'word': w,
        'parent_sha256': request['parent_sha256'], 'request_sha256': sha(job/'request.json'),
        'conic': conic, 'checks': ['Sage trace', 'reversed manual trace', 'Fraction polynomial identities']})


def admit(job):
    """Only independently replayed finite rank certificates enter the cascade."""
    import parent_foundry_worker as worker
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    request = read(job/'request.json'); parent_path = ROOT/request['parent']
    if sha(parent_path) != request['parent_sha256']:
        raise ArithmeticError('parent seal differs')
    if sha(ROOT/request['incidence']) != request['incidence_sha256']:
        raise ArithmeticError('incidence provenance changed')
    cells = read(ROOT/request['incidence'])['cells']
    if request['candidates'] != [r for r in cells if r['status'] == 'SPLIT_REQUIRES_ADMISSION']:
        raise ArithmeticError('candidate roster differs from sealed incidence')
    parent = read(parent_path); packet = worker.generic_seed(parent, request['parameter'], job)
    if packet is None:
        save(job/'admission.json', {'status': 'UNRESOLVED_GENERIC_SPECIALIZATION', 'candidates': []})
        return
    model = tuple(map(F, packet['curve'])); points = [tuple(map(F, p)) for p in packet['points']]
    frame = FinitePointAdmission(model, points, prime_bound=1000)
    d = F(request['parameter']).denominator
    decisions = []
    for row in request['candidates']:
        x, y = map(F, row['point']); point = (x*d**4, y*d**6)
        if len(frame.points) >= 32:
            decisions.append({'mask': row['mask'], 'status': 'NOT_TESTED_TARGET_REACHED'})
        else:
            decisions.append({'mask': row['mask'], 'point': list(map(str, point)), **frame.consider(point)})
    rank = len(frame.points)
    if rank > packet['rank_lower_bound']:
        packet.update(points=[list(map(str, p)) for p in frame.points], rank_lower_bound=rank,
            proof=checked_rank(model, tuple(frame.points), frame.primes, packet['proof']['no_rational_2_torsion_prime']))
        worker.verify(packet, parent, request['parameter'])
    save(job/'seed.json', packet)
    save(job/'admission.json', {'status': 'PASS_CERTIFIED_SEED' if rank > 17 else 'NO_CERTIFIED_NEW_DIRECTION',
        'seed_rank_lower_bound': rank, 'generic_rank': 17, 'candidates': decisions,
        'seed_sha256': sha(job/'seed.json'), 'request_sha256': sha(job/'request.json'),
        'boundary': 'UNKNOWN_FINITE_COLUMN_IN_SPAN is not dependence. No inferred gain from a split.'})


def bounded(plan, job, command, seconds):
    from research_runtime.supervisor import run as supervise, Limits
    from run_parent_foundry import env
    receipt = job/'supervision-receipt.json'
    if receipt.exists():
        return read(receipt)
    child_env = env()
    child_env['PATH'] = str(Path(plan['sage']).parent)+os.pathsep+child_env.get('PATH', '')
    result = supervise(command, limits=Limits(seconds, plan['rss_bytes']),
                       log_path=job/'worker.log', cwd=Path(plan['root']), env=child_env)
    result['output_sha256'] = {p.name: sha(p) for name in
        ('result.json', 'admission.json', 'seed.json', 'packet.json', 'verified.json')
        if (p := job/name).exists()}
    save(receipt, result)
    return result


def successful(job, name):
    receipt = read(job/'supervision-receipt.json')
    good = receipt.get('outcome') == 'completed' and receipt.get('returncode') == 0 and (job/name).exists()
    if good and any(not (job/n).exists() or sha(job/n) != h for n, h in receipt.get('output_sha256', {}).items()):
        raise ArithmeticError('completed worker output changed: '+str(job))
    return good and name in receipt.get('output_sha256', {})


def run(folder):
    import fcntl
    with (folder/'campaign.lock').open('a') as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan, root = guard(folder)
        script = root/'elliptic-curves/cas/run_euclidean_seed_foundry.py'
        parent_path = root/'seed-inputs/parent.json'; parent_hash = sha(parent_path)
        parent = read(parent_path)
        base = {'parent': str(parent_path.relative_to(root)), 'parent_sha256': parent_hash}
        def build(row):
            if (folder/'STOP').exists(): return
            job = root/'seed-traces'/str(row['mask'])
            save(job/'request.json', {**base, **row})
            bounded(plan, job, [plan['sage'], '-python', str(script), '_build', '--folder', str(job)], plan['trace_seconds'])
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            list(pool.map(build, plan['orbits']))
        bank = []
        for row in plan['orbits']:
            job = root/'seed-traces'/str(row['mask'])
            if (job/'supervision-receipt.json').exists() and successful(job, 'result.json'):
                r = read(job/'result.json')
                if r['status'] == 'PASS_EXACT_CONIC_MAP':
                    if r['parent_sha256'] != parent_hash or r['request_sha256'] != sha(job/'request.json'):
                        raise ArithmeticError('conic receipt seal differs')
                    replay = conic_from_trace(parent['A_coefficients_low_to_high'],
                        parent['B_coefficients_low_to_high'], *(r['conic'][k] for k in ('h', 'nx', 'ny')))
                    if replay != r['conic']:
                        raise ArithmeticError('cached polynomial map replay differs')
                    bank.append(r)
        def process(slot):
            if (folder/'STOP').exists(): return
            trial = root/'seed-trials'/str(slot['index'])
            t = F(slot['parameter'])
            A, B = (evaluate(parent[k], t) for k in ('A_coefficients_low_to_high', 'B_coefficients_low_to_high'))
            if 4*A**3+27*B**2 == 0:
                save(trial/'incidence.json', {**slot, 'status': 'SINGULAR_FIBRE', 'cells': []})
                return
            cells = [{'mask': r['mask'], **incidence(r['conic'], t)} for r in bank]
            save(trial/'incidence.json', {**slot, 'status': 'SCANNED', 'cells': cells,
                'orbits_planned': len(plan['orbits']), 'orbits_constructed': len(bank)})
            candidates = [r for r in cells if r['status'] == 'SPLIT_REQUIRES_ADMISSION']
            seed = None
            if candidates:
                job = trial/'admit'
                save(job/'request.json', {**base, 'parameter': str(t), 'candidates': candidates,
                    'incidence': str((trial/'incidence.json').relative_to(root)),
                    'incidence_sha256': sha(trial/'incidence.json')})
                bounded(plan, job, [plan['sage'], '-python', str(script), '_admit', '--folder', str(job)], plan['admission_seconds'])
                if successful(job, 'admission.json') and read(job/'admission.json')['status'] == 'PASS_CERTIFIED_SEED':
                    seed = job/'seed.json'
            arms = [('control', None)] if slot['control'] else []
            if seed is not None and read(seed)['rank_lower_bound'] < 32:
                arms.append(('seeded', seed))
            for arm, packet in arms:
                job = trial/arm
                req = {**base, 'kind': 'panel', 'family': 'det1092', 'parameter': str(t),
                       'allowance': plan['point_calls'], 'height': plan['height'], 'bank_index': 0}
                if packet is not None:
                    req.update(packet=str(packet.relative_to(root)), packet_sha256=sha(packet))
                save(job/'request.json', req)
                bounded(plan, job, [plan['sage'], '-python', str(root/'elliptic-curves/cas/parent_foundry_worker.py'),
                                    '--job', str(job)], plan['cascade_seconds'])
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            list(pool.map(process, plan['parameters']))
        status(folder)


def status(folder):
    plan = read(folder/'plan.json'); root = Path(plan['root'])
    results = {'schema': 'euclidean-seed-foundry.status.v1', 'parameters_planned': len(plan['parameters']),
               'orbits_planned': len(plan['orbits']), 'trace_outcomes': {}, 'slots': [],
               'boundary': plan['claim_boundary']}
    for r in plan['orbits']:
        job = root/'seed-traces'/str(r['mask'])
        if not (job/'supervision-receipt.json').exists(): outcome = 'NOT_COMPLETED'
        elif successful(job, 'result.json'): outcome = read(job/'result.json')['status']
        else: outcome = 'UNRESOLVED_'+str(read(job/'supervision-receipt.json').get('outcome'))
        results['trace_outcomes'][outcome] = results['trace_outcomes'].get(outcome, 0)+1
    for slot in plan['parameters']:
        trial = root/'seed-trials'/str(slot['index']); row = {**slot, 'status': 'NOT_COMPLETED'}
        if (trial/'incidence.json').exists():
            data = read(trial/'incidence.json'); row.update(status=data['status'],
                square_splits=sum(c['status'] == 'SPLIT_REQUIRES_ADMISSION' for c in data['cells']))
        for kind, filename in (('admit', 'admission.json'), ('control', 'result.json'), ('seeded', 'result.json')):
            job = trial/kind
            if (job/'supervision-receipt.json').exists():
                receipt = read(job/'supervision-receipt.json')
                row[kind] = {'supervision': receipt, 'result': read(job/filename) if successful(job, filename) else None}
        results['slots'].append(row)
    # Live summary is replaceable, unlike evidence receipts. Never aggregates
    # seeded and control arms as independent fibres, or timeouts as nulls.
    tmp = folder/'STATUS.json.tmp'; tmp.write_text(json.dumps(results, indent=2, sort_keys=True)+'\n')
    tmp.replace(folder/'STATUS.json')
    print(json.dumps({'trace_outcomes': results['trace_outcomes'], 'slots': len(results['slots']),
                      'status_file': str(folder/'STATUS.json')}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'run', 'status', '_build', '_admit'))
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    parser.add_argument('--sage', default='sage'); parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--orbits', type=int, default=128); parser.add_argument('--orbit-offset', type=int, default=0)
    parser.add_argument('--parameters', type=int, default=48); parser.add_argument('--parameter-offset', type=int, default=0)
    parser.add_argument('--control-every', type=int, default=4); parser.add_argument('--trace-seconds', type=int, default=25)
    args = parser.parse_args(); args.folder = args.folder.resolve()
    if args.mode == 'prepare': prepare(args)
    elif args.mode == 'run':
        _, root = guard(args.folder)
        frozen_script = root/'elliptic-curves/cas/run_euclidean_seed_foundry.py'
        if frozen_script.resolve() != Path(__file__).resolve():
            subprocess.run([sys.executable, str(frozen_script), 'run', '--folder', str(args.folder)], check=True)
        else: run(args.folder)
    elif args.mode == 'status': status(args.folder)
    elif args.mode == '_build': build_trace(args.folder)
    else: admit(args.folder)


if __name__ == '__main__':
    main()

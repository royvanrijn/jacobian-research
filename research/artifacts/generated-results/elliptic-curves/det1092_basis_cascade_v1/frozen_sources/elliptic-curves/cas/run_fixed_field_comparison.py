#!/usr/bin/env python3
"""Frozen, checkpointed fixed-field CT comparison and radical-only point pilot.

The local-pilot certificate and its replayer certify the independent anchor
classes and W_u. This adapter reuses those intersections, the shared labelled
arithmetic context, norm-cover maps and Fisher pairing. No BNF is requested.
Discovery runs in the shared supervisor; replay does no point/conic search.
"""
import argparse
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import shutil
import sys
import time

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
sys.path.insert(0, str(CAS))
from research_runtime.store import FactStore, checkpoint, digest
from research_runtime.supervisor import Limits, run as supervise

SOURCE = ROOT/'artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json'
PROTOCOL = ROOT/'elliptic-curves/data/fixed_field_comparison_v1.json'
SUMMARY = ROOT/'artifacts/generated-results/elliptic-curves/fixed_field_comparison_v1.json'


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def extension_gate(protocol, store):
    """The second panel cannot run before the complete first comparison replays."""
    if 'extension_freeze' not in protocol:
        return
    path = ROOT/protocol['extension_freeze']
    if file_hash(path) != protocol['extension_freeze_sha256']:
        raise ArithmeticError('extension selection was changed')
    freeze = read(path)
    if protocol['parameters'] != freeze['parameters']:
        raise ArithmeticError('extension parameter scope changed')
    initial = read(ROOT/freeze['same_pairing_and_point_limits_as'])
    for key in ('entry_seconds', 'ct_curve_seconds', 'point_entry_seconds', 'point_curve_seconds',
                'point_cover_cap', 'point_height', 'local_witness_node_cap', 'rss_bytes',
                'pari_stack_bytes', 'seed', 'cover_policy'):
        if protocol[key] != initial[key]:
            raise ArithmeticError('extension budget changed')
    seen = set()
    for path, parent in freeze['first_stage_verified_evidence'].items():
        if file_hash(ROOT/path) != parent['sha256']:
            raise ArithmeticError('initial comparison evidence changed')
        evidence = read(ROOT/path)
        result = verify(evidence, store)
        if result['profile'] != parent['profile']:
            raise ArithmeticError('initial comparison profile changed')
        seen.add(evidence['parameter_u'])
    if seen != {'0', '-2', '1', '2'}:
        raise ArithmeticError('initial comparison incomplete')


def setup(protocol, parameter, store, retained=None):
    from sage.all import QQ, pari
    from research_runtime.arithmetic import TwoTorsionContext
    from research_runtime.sage_arithmetic import SageArithmetic
    from research_runtime.sage_subspace import SageSubspaceBackend
    from research_runtime.subspace import GlobalSquareclasses, local_intersection
    source_path = ROOT/protocol['local_source']
    if file_hash(source_path) != protocol['local_source_sha256']:
        raise ArithmeticError('local source binding changed')
    source = read(source_path)
    if source['status'] != 'PASS_EXACT_FULL_SPAN_LOCAL_INTERSECTIONS_NO_CLASS_GROUP':
        raise ArithmeticError('uncertified local input')
    row = next(r for r in source['runs'] if r['parameter_u'] == parameter)
    if not row['all_local_kummer_images_complete']:
        raise ArithmeticError('incomplete local input')
    masks = list(local_intersection(20, [r['known_span_quotient_rows'] for r in
        row['finite_local_conditions']+[row['real_local_condition']]]))
    if masks != [r['mask'] for r in row['W_u_basis']] or len(masks) != row['W_u_dimension']:
        raise ArithmeticError('local intersection mismatch')
    if retained is not None:
        if retained['protocol_hash'] != digest(protocol) or retained['parameter_u'] != parameter:
            raise ArithmeticError('misbound comparison evidence')
        # Cover, root and support witnesses are checked directly below; only
        # model transport and the shared order are needed as cached setup.
        snapshot = retained['arithmetic_facts']
        store.import_snapshot({**snapshot, 'facts': [item for item in snapshot['facts']
            if item['record']['key']['namespace'] in ('arithmetic-context', 'two-torsion/order')]})
    pari.setrand(protocol['seed'])
    pari.allocatemem(64_000_000, protocol['pari_stack_bytes'], silent=True)
    arithmetic = SageArithmetic(store)
    algebra = TwoTorsionContext(source['anchor']['base_polynomial_ascending'])
    context = arithmetic.prepare_congruent(row['raw_curve_ainvariants'], algebra, (0, 1, QQ(parameter)),
        factor_primes=row['complete_finite_place_support'], discover=retained is None)
    arithmetic.field(algebra, factor_primes=[r['prime'] for r in source['anchor']['base_discriminant_factorization']],
                     discover=retained is None)
    classes = GlobalSquareclasses(algebra.key,
        source['anchor']['known_kummer_basis_beta_power_coordinates'], protocol['local_source_sha256'])
    backend = SageSubspaceBackend(arithmetic, context, None,
        local_candidate_cap=protocol['local_witness_node_cap'], cover_policy=protocol['cover_policy'])
    if retained is not None and retained['context'] != json.loads(json.dumps(context.record())):
        raise ArithmeticError('context mismatch')
    # The source theorem certifies independence; these direct checks also bind
    # each representative to the displayed rational anchor point and norm.
    from research_runtime.arithmetic import CurveModel
    anchor = CurveModel(source['anchor']['short_model_ainvariants'])
    for point, beta in zip(source['anchor']['known_points_on_short_model'], classes.representatives):
        if not anchor.contains(point) or tuple(map(QQ, beta)) != (QQ(point[0]), -1, 0):
            raise ArithmeticError('anchor point/class mismatch')
        if QQ(pari.nfeltnorm(backend.nf, backend.element(beta))) != QQ(point[1])**2:
            raise ArithmeticError('anchor norm mismatch')
    return source, row, masks, context, classes, backend


def bounded(seconds, function):
    import signal
    def expired(signum, frame):
        raise TimeoutError('declared entry budget exhausted')
    old = signal.signal(signal.SIGALRM, expired)
    signal.alarm(seconds)
    try:
        return function()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def summarize(masks, ct, searches, parameter):
    from research_runtime.binary import BinaryBasis
    from research_runtime.subspace import restricted_radical
    radical = restricted_radical(masks, ct['matrix'], global_dimension=20)
    realized = BinaryBasis(20)
    for search in searches:
        if search.get('points'):
            realized, _ = realized.append(search['mask'])
    dimension = 20 if parameter == '0' else realized.rank
    return {**radical, 'certified_realized_dimension': dimension,
            'realized_dimension_is_lower_bound': parameter != '0',
            'profile': [len(masks), radical['pairing_rank'], radical['radical_dimension'], dimension]}


def search_queue(radical, cap):
    from research_runtime.binary import combine
    generators = radical['radical_global_masks']
    # Enumerate a bounded prefix by coordinate weight; never materialize 2^r.
    for count, indices in enumerate(indices for weight in range(1, len(generators)+1)
                                   for indices in combinations(range(len(generators)), weight)):
        if count == cap:
            break
        yield combine(sum(1 << i for i in indices), generators)


def worker(protocol, parameter, work, store, phase):
    from sage.all import QQ, pari
    from research_runtime.binary import combine
    source, row, masks, context, classes, backend = setup(protocol, parameter, store)
    identity = {'protocol_hash': digest(protocol), 'parameter_u': parameter, 'context': context.record()}
    binding = work/'binding.json'
    if binding.exists() and read(binding) != json.loads(json.dumps(identity)):
        raise ArithmeticError('checkpoint directory belongs to another experiment')
    checkpoint(binding, identity)
    covers = {}
    def cover(mask):
        if mask not in covers:
            print('COVER_START', parameter, mask, flush=True)
            c = backend.cover(context, classes, mask)
            backend.verify_cover(context, classes, mask, c)
            covers[mask] = c
            print('COVER_DONE', parameter, mask, flush=True)
        return covers[mask]
    def pair(trio, path):
        cs = [cover(m) for m in trio]
        if path.exists():
            result = backend._pair(classes, trio, cs, retained=read(path))
        else:
            result = backend._pair(classes, trio, cs)
            # Independent PARI Hilbert check in addition to shared elementary formula.
            for term in result['local_terms']:
                if term['place'] != 'infinity' and int(pari.hilbert(QQ(cs[1]['quartic'][4]),
                        QQ(term['gamma_value']), term['place'])) != term['hilbert_symbol']:
                    raise ArithmeticError('Hilbert implementation disagreement')
            checkpoint(path, result)
        return result
    if phase == 'ct':
        pairs, failures = [], []
        for i, j in combinations(range(len(masks)), 2):
            start = time.monotonic()
            try:
                result = bounded(protocol['entry_seconds'], lambda: pair(
                    (masks[i], masks[j], masks[i]^masks[j]), work/f'pair-{i}-{j}.json'))
                pairs.append(result)
                print('PAIR_DONE', parameter, i, j, result['value'], round(time.monotonic()-start, 3), flush=True)
            except Exception as exc:
                failures.append({'i': i, 'j': j, 'error': f'{type(exc).__name__}: {exc}'})
                checkpoint(work/'unknown.json', failures)
                print('PAIR_UNKNOWN', parameter, i, j, failures[-1]['error'], flush=True)
        if failures:
            raise ArithmeticError('incomplete matrix; no radical or point search permitted')
        rows = [[0]*len(masks) for _ in masks]
        for (i, j), result in zip(combinations(range(len(masks)), 2), pairs):
            rows[i][j] = rows[j][i] = result['value']
        ct = {'matrix': rows, 'pairs': pairs,
              'sum_covers': [covers[m] for m in sorted(set(covers)-set(masks))]}
        backend.verify_ct(context, classes, masks, [covers[m] for m in masks], ct)
        summary = summarize(masks, ct, [], parameter)
        if parameter == '0' and summary['profile'] != [20, 0, 20, 20]:
            raise ArithmeticError('anchor correctness control failed')
        # Reverse order, a sum, and a radical generator test symmetry/bilinearity.
        checks = []
        rad = summary['radical_coordinates']
        for a, b in [(2, 1), (3, 4)]+([(rad[0], 1)] if rad and rad[0] != 1 else []):
            left, right = combine(a, masks), combine(b, masks)
            result = bounded(protocol['entry_seconds'], lambda: pair((left, right, left^right), work/f'cross-{a}-{b}.json'))
            expected = sum(((a>>i)&1)*rows[i][j]*((b>>j)&1) for i in range(len(masks)) for j in range(len(masks))) % 2
            if result['value'] != expected:
                raise ArithmeticError('symmetry/bilinearity cross-check failed')
            checks.append({'coordinates': [a, b], 'pair': result})
        evidence = {**identity, 'protocol': protocol, 'status': 'VERIFIED_RESTRICTED_CT_RADICAL',
                    'admissible_masks': masks, 'covers': [covers[m] for m in sorted(covers)],
                    'ct': ct, 'cross_checks': checks, 'searches': [], 'summary': summary,
                    'arithmetic_facts': store.snapshot()}
        checkpoint(work/'ct.json', evidence)
        print('PROFILE', parameter, summary['profile'], flush=True)
        return
    evidence = read(work/'ct.json')
    verify(evidence, store)
    covers.update({c['mask']: c for c in evidence['covers']})
    searches = []
    if parameter != '0':
        for mask in search_queue(evidence['summary'], protocol['point_cover_cap']):
            path = work/f'search-{mask}.json'
            def solve():
                c = cover(mask)
                q = backend.R(list(map(QQ, c['quartic'])))
                # Remove the known rational square scale, to avoid wasting the
                # height budget on the ordinate of an invariant-normalized model.
                den = q.denominator()
                content = abs(QQ(pari.content(q*den**2)))
                factors = backend.arithmetic.factor_integer(content, discover=True)
                square = QQ(1)
                for p, e in factors:
                    square *= QQ(p)**(e//2)
                scale = square/den
                reduced = q/scale**2
                raw = pari.hyperellratpoints(reduced, protocol['point_height'])
                points = [(QQ(p[0]), QQ(1), QQ(p[1])*scale) for p in raw]
                if q[4].is_square():
                    points.append((QQ(1), QQ(0), q[4].sqrt()))
                hits = []
                for point in points:
                    recovered = backend.point_from_cover(context, classes, mask, c, point)
                    if recovered is None:
                        raise ArithmeticError('unexpected identity on nonzero independent class')
                    hits.append({'quartic_point': list(map(str, point)), 'raw_point': list(recovered)})
                return {'mask': mask, 'status': 'BOUNDED_SEARCH_COMPLETE',
                        'search_quartic': list(map(str, [reduced[i] for i in range(5)])),
                        'quartic_y_scale': str(scale), 'height': protocol['point_height'], 'points': hits}
            if path.exists():
                result = read(path)
            else:
                try:
                    result = bounded(protocol['point_entry_seconds'], solve)
                except Exception as exc:
                    result = {'mask': mask, 'status': 'UNKNOWN', 'error': f'{type(exc).__name__}: {exc}', 'points': []}
                checkpoint(path, result)
            searches.append(result)
            print('SEARCH_DONE', parameter, mask, result['status'], len(result['points']), flush=True)
    evidence.update({'covers': [covers[m] for m in sorted(covers)], 'searches': searches,
        'summary': summarize(masks, evidence['ct'], searches, parameter), 'arithmetic_facts': store.snapshot(),
        'point_phase_complete': True})
    verify(evidence, store)
    checkpoint(work/'result.json', evidence)
    print('PROFILE', parameter, evidence['summary']['profile'], flush=True)


def verify(evidence, store):
    from sage.all import QQ
    from research_runtime.binary import combine
    protocol, parameter = evidence['protocol'], evidence['parameter_u']
    source, row, masks, context, classes, backend = setup(protocol, parameter, store, evidence)
    if evidence['status'] != 'VERIFIED_RESTRICTED_CT_RADICAL' or evidence['admissible_masks'] != masks:
        raise ArithmeticError('incomplete or wrong subspace evidence')
    covers = {}
    for c in evidence['covers']:
        if c['mask'] in covers:
            raise ArithmeticError('duplicate cover')
        backend.verify_cover(context, classes, c['mask'], c)
        covers[c['mask']] = c
    backend.verify_ct(context, classes, masks, [covers[m] for m in masks], evidence['ct'])
    rows = evidence['ct']['matrix']
    for check in evidence['cross_checks']:
        a, b = check['coordinates']
        left, right = combine(a, masks), combine(b, masks)
        result = backend._pair(classes, (left, right, left^right), [covers[m] for m in (left, right, left^right)], retained=check['pair'])
        expected = sum(((a>>i)&1)*rows[i][j]*((b>>j)&1) for i in range(len(masks)) for j in range(len(masks))) % 2
        if result['value'] != expected:
            raise ArithmeticError('invalid cross-check')
    summary = summarize(masks, evidence['ct'], evidence['searches'], parameter)
    if summary != evidence['summary'] or (parameter == '0' and summary['profile'] != [20, 0, 20, 20]):
        raise ArithmeticError('profile or positive control mismatch')
    queue = [] if parameter == '0' else list(search_queue(summary, protocol['point_cover_cap']))
    if evidence.get('point_phase_complete') and [r['mask'] for r in evidence['searches']] != queue:
        raise ArithmeticError('incomplete point-search schedule')
    for search in evidence['searches']:
        mask = search['mask']
        if mask not in queue or search['status'] not in ('UNKNOWN', 'BOUNDED_SEARCH_COMPLETE'):
            raise ArithmeticError('point search escaped the restricted radical')
        if search['status'] == 'UNKNOWN':
            if search['points']:
                raise ArithmeticError('unknown search cannot certify points')
            continue
        c = covers[mask]
        q = backend.R(list(map(QQ, c['quartic'])))
        if backend.R(list(map(QQ, search['search_quartic'])))*QQ(search['quartic_y_scale'])**2 != q or search['height'] != protocol['point_height']:
            raise ArithmeticError('wrong search model or bound')
        for point in search['points']:
            recovered = backend.point_from_cover(context, classes, mask, c, point['quartic_point'])
            if recovered is None or list(recovered) != point['raw_point']:
                raise ArithmeticError('invalid rational realization')
    return summary


def verify_panel(document, store):
    """Replay the complete frozen comparison, preserving the baseline dependency."""
    from research_runtime.witnesses import retained_source
    for path, expected in document['source_hashes'].items():
        if Path(path).suffix in ('.py', '.sage'):
            retained_source(ROOT, path, expected)
        elif file_hash(ROOT/path) != expected:
            raise ArithmeticError(f'changed panel source: {path}')
    expected_parameters = ['-3', '-2', '-1', '0', '1', '2', '3']
    if [row['parameter_u'] for row in document['curves']] != expected_parameters:
        raise ArithmeticError('incomplete frozen panel')
    entries = covers = crosses = searches = 0
    for row in document['curves']:
        path = ROOT/row['evidence']
        if file_hash(path) != row['evidence_sha256']:
            raise ArithmeticError('panel evidence hash mismatch')
        if row['parameter_u'] == '-1':
            # This row is a retained theorem dependency, not a new campaign.
            baseline = read(path)['arithmetic']
            profile = [18, baseline['pairing_rank'], baseline['restricted_radical_dimension'], 0]
            if row['profile'] != profile or row['role'] != 'retained baseline':
                raise ArithmeticError('incorrect u=-1 comparator')
            continue
        evidence = read(path)
        result = verify(evidence, store)
        if not evidence.get('point_phase_complete') or evidence['parameter_u'] != row['parameter_u']:
            raise ArithmeticError('unfinished or misbound curve experiment')
        if result['profile'] != row['profile'] or result['radical_global_masks'] != row['radical_anchor_masks']:
            raise ArithmeticError('panel result mismatch')
        entries += len(evidence['ct']['pairs'])
        covers += len(evidence['covers'])
        crosses += len(evidence['cross_checks'])
        searches += len(evidence['searches'])
        if 'extension_freeze' in evidence['protocol']:
            extension_gate(evidence['protocol'], store)
    counts = dict(matrix_entries=entries, cover_maps=covers, cross_checks=crosses, point_searches=searches)
    if document['new_computation_counts'] != counts:
        raise ArithmeticError('incorrect computation counts')
    if document['success_criterion_met'] is not False or any(
            row['profile'][3] for row in document['curves'] if row['parameter_u'] != '0'):
        raise ArithmeticError('this negative-comparison certificate was changed')
    return {'status': 'PASS', 'profiles': {row['parameter_u']: row['profile'] for row in document['curves']}, **counts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=PROTOCOL)
    parser.add_argument('--u', type=str)
    parser.add_argument('--phase', choices=['ct', 'points'], default='ct')
    parser.add_argument('--workdir', type=Path, default=ROOT/'artifacts/local/fixed-field-comparison-v1')
    parser.add_argument('--verify', type=Path)
    parser.add_argument('--check', action='store_true', help='replay the retained frozen-panel summary')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    protocol = read(args.protocol)
    if args.check:
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            result = verify_panel(read(SUMMARY), FactStore(directory))
        print(json.dumps(result, indent=2))
        return
    if args.verify:
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            result = verify(read(args.verify), FactStore(directory))
        print(json.dumps(result, indent=2))
        return
    if args.u not in protocol['parameters']:
        parser.error('parameter is outside the frozen panel')
    work = args.workdir.resolve()/('u'+args.u.replace('-', 'm').replace('/', 'd'))
    work.mkdir(parents=True, exist_ok=True)
    if args.worker:
        store = FactStore(args.workdir.resolve()/'cache')
        extension_gate(protocol, store)
        worker(protocol, args.u, work, store, args.phase)
        return
    limits = Limits(protocol['ct_curve_seconds'] if args.phase == 'ct' else protocol['point_curve_seconds'],
                    protocol['rss_bytes'], pari_stack_bytes=protocol['pari_stack_bytes'])
    command = [shutil.which('sage'), '-python', str(Path(__file__).resolve()), '--worker',
               '--protocol', str(args.protocol.resolve()), '--u', args.u, '--phase', args.phase,
               '--workdir', str(args.workdir.resolve())]
    result = supervise(command, limits=limits, log_path=work/f'{args.phase}.log',
                       result_path=work/('ct.json' if args.phase == 'ct' else 'result.json'),
                       checkpoint_path=work/f'{args.phase}.supervisor.json')
    print(json.dumps(result, indent=2))
    if result['outcome'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()

"""Successor epochs: expose the latest admitted generator block first."""
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from math import ceil
from pathlib import Path
import resource
import signal
import time

from finite_cancellation_corpus import ROOT, canonical, digest, write
from cancellation_scheduler_cpu import cpu, atomic

CAS = Path(__file__).resolve().parent


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


from cancellation_basis_epoch import native_rank, rebuild


def run(case, arm, plan, dest, protocol_hash):
    from finite_cancellation_features import alarm
    from cysignals.signals import AlarmInterrupt
    from cancellation_basis_early_order import visit_order
    from cancellation_basis_exploration import Exposure, prepare
    from future_point_admission import FinitePointAdmission
    from half_lattice_pointed_sieve import linear_combination
    from memory_rank_certificate import checked_rank
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import execute, replay
    from v3_warm_engine import certified_state

    need(arm in plan['arms'], 'unsealed arm')
    need(not (dest/'events.json').exists() and not (dest/'result.json').exists(), 'preserve started arm')
    dest.mkdir(parents=True, exist_ok=True)
    resource.setrlimit(resource.RLIMIT_CPU, (plan['hard_process_cpu_seconds']-5, plan['hard_process_cpu_seconds']))
    signal.signal(signal.SIGALRM, alarm)
    events = []; calls = []; epochs = []; exposure = Exposure()
    components = {'rank_certification': 0., 'bank': 0., 'preparation': 0., 'point_calls': 0., 'cloud_reconciliation': 0.}
    packet = case['seed']; initial = len(packet['points']); curve = tuple(map(F, packet['curve']))
    admission = FinitePointAdmission(curve, [tuple(map(F, p)) for p in packet['points']], prime_bound=plan['prime_bound'])
    tick = cpu(); rank = native_rank(packet)
    need(rank['rank'] == initial, 'starting rank differs')
    write(dest/'initial-rank.json', {'status': 'PASS', 'rank': rank, 'seed_sha256': digest(canonical(packet))})
    components['rank_certification'] += cpu()-tick
    mapper = SourceFileLoader('basis_epoch_factor_free', str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000, silent=True)
    status = 'DECLARED_BOXES_EXHAUSTED'; first_cloud_gain = 0; first_gain_cpu = None
    epoch = 0; bank = None; unknowns = 0; fresh_start = initial

    def record(event):
        events.append(event); atomic(dest/'events.json', events)

    def remaining(): return plan['working_cpu_seconds']-cpu()

    def open_bank():
        nonlocal bank, unknowns, status
        folder = dest/f'epoch-{epoch:02d}'
        tick = cpu()
        record({'kind': 'bank_start', 'epoch': epoch, 'rank': len(packet['points']),
                'seed_sha256': digest(canonical(packet))})
        signal.alarm(max(1, min(plan['bank_wall_seconds'], ceil(max(0, remaining())))))
        try:
            bank, receipt = rebuild(packet, case, plan, folder)
        except (TimeoutError, AlarmInterrupt, ArithmeticError, ValueError, AssertionError) as error:
            bank = None; unknowns += 1; status = 'BANK_UNKNOWN'
            record({'kind': 'bank_unknown', 'epoch': epoch, 'error': type(error).__name__+': '+str(error), 'cpu_seconds': cpu()-tick})
            return False
        finally:
            signal.alarm(0); components['bank'] += cpu()-tick
        epochs.append({'epoch': epoch, 'rank': len(packet['points']), 'verification_sha256': sha(folder/'verification.json'),
                       'bank_sha256': sha(folder/'bank.json'), 'cpu_seconds': receipt['cpu_seconds']})
        record({'kind': 'bank_ready', **epochs[-1], 'fresh_generator_start': fresh_start,
                'visit_order': visit_order(bank['centres'], initial, fresh_start, epoch > 0)})
        return True

    if remaining() <= 0:
        status = 'CPU_CAP_BEFORE_BANK'
    elif open_bank():
        while bank is not None:
            search_packet = bank['seed']
            basis = [tuple(map(F, p)) for p in search_packet['points']]
            state = certified_state(curve, basis, search_packet['proof'])
            refresh = False
            for ci in visit_order(bank['centres'], initial, fresh_start, epoch > 0):
                centre = bank['centres'][ci]
                if remaining() <= 0:
                    status = 'CPU_CAP'; break
                tick = cpu(); word = centre['representative']
                signal.alarm(max(1, min(10, ceil(remaining()))))
                try:
                    prepared = prepare(curve, basis, word, mapper)
                except (TimeoutError, AlarmInterrupt, ArithmeticError, ValueError, AssertionError) as error:
                    unknowns += 1
                    record({'kind': 'preparation_unknown', 'epoch': epoch, 'centre': ci,
                            'error': type(error).__name__+': '+str(error)})
                    continue
                finally:
                    signal.alarm(0); components['preparation'] += cpu()-tick
                path = dest/f'epoch-{epoch:02d}/prepared-{ci:03d}.json'
                write(path, prepared)
                record({'kind': 'prepared', 'epoch': epoch, 'centre': ci, 'file': str(path.relative_to(dest)), 'sha256': sha(path)})
                anchor = list(map(str, linear_combination(curve, basis, word)))
                for mi, model in enumerate(prepared['models']):
                    mapping = model['mapping']; height = plan['height']
                    if remaining() <= 0:
                        status = 'CPU_CAP'; break
                    if len(calls) >= plan['maximum_calls']:
                        status = 'CALL_CAP'; break
                    if exposure.seen(anchor, mapping, height):
                        record({'kind': 'exact_completed_duplicate', 'epoch': epoch, 'centre': ci, 'model': mi})
                        continue
                    witness = exposure.inspect(anchor, mapping, prepared['models'][0]['mapping'], height)
                    search = PointedQuarticSearch(state=state, centre={'coefficients': word}, coordinate_policy=mapping['coordinate_policy'])
                    tick = cpu()
                    transcript, points = execute(search, mapping, height, min(plan['point_wall_seconds'], max(.001, remaining())), plan['gp_sha256'])
                    elapsed = cpu()-tick; components['point_calls'] += elapsed
                    tick = cpu()
                    returned = replay(search, mapping, transcript)
                    need(returned == points, 'immediate exact transcript replay differs')
                    before = len(admission.points); decisions = []
                    for point in returned:
                        decisions.append({'point': list(map(str, point)), **admission.consider(point)})
                    gain = len(admission.points)-before
                    components['cloud_reconciliation'] += cpu()-tick
                    complete = transcript['status'] == 'bounded_search_complete'
                    # Backend infinity is independently checked even on a timeout.
                    # Such points may be certified; no finite completed prefix is claimed.
                    rank_path = None
                    if gain:
                        tick = cpu()
                        proof = checked_rank(curve, admission.points, admission.primes, packet['proof']['no_rational_2_torsion_prime'])
                        packet = {'curve': case['seed']['curve'], 'points': [list(map(str, p)) for p in admission.points], 'proof': proof}
                        rank = native_rank(packet)
                        need(rank['rank'] == len(admission.points), 'complete cloud independent rank differs')
                        rank_path = dest/f'rank-{len(calls):03d}.json'
                        write(rank_path, {'packet': packet, 'independent_rank': rank})
                        components['rank_certification'] += cpu()-tick
                        if first_gain_cpu is None:
                            first_gain_cpu = cpu(); first_cloud_gain = gain
                    exposure.observe(anchor, mapping, height, complete)
                    path = dest/f'call-{len(calls):03d}.json'
                    write(path, {'epoch': epoch, 'centre': ci, 'model': mi, 'word': word,
                                 'mapping': mapping, 'search': transcript, 'admission': decisions,
                                 'gain': gain, 'rank_after': len(admission.points), 'tail_witness': witness,
                                 'call_cpu_seconds': elapsed, 'cpu_after_certification': cpu(),
                                 'rank_file': rank_path.name if rank_path else None,
                                 'rank_sha256': sha(rank_path) if rank_path else None})
                    calls.append({'file': path.name, 'sha256': sha(path), 'epoch': epoch, 'gain': gain,
                                  'status': transcript['status'], 'tail_witness': witness is not None})
                    record({'kind': 'call', **calls[-1]})
                    if len(admission.points) >= initial+plan['target_directions']:
                        status = 'CERTIFIED_TARGET'; break
                    if gain and arm == 'basis_refresh':
                        refresh = True; break
                if refresh or status in ('CPU_CAP', 'CALL_CAP', 'CERTIFIED_TARGET'):
                    break
            if status in ('CPU_CAP', 'CALL_CAP', 'CERTIFIED_TARGET') or not refresh:
                break
            if remaining() <= 0:
                status = 'CPU_CAP_BEFORE_REFRESH'; break
            fresh_start = len(search_packet['points'])
            epoch += 1
            # Entire reconciled basis, new bank and new local job order. Only
            # exact completed boxes survive; no old radial exposure is reused.
            if not open_bank():
                break

    write(dest/'final-rank.json', {'packet': packet, 'independent_rank': rank})
    atomic(dest/'events.json', events)
    search_end = cpu()
    result = {'case': case['id'], 'family': case['family'], 'stratum': case['stratum'], 'arm': arm,
              'status': status, 'initial_rank': initial, 'rank_lower_bound': len(packet['points']),
              'new_directions': len(packet['points'])-initial,
              'success': len(packet['points']) >= initial+plan['target_directions'],
              'first_cloud_gain': first_cloud_gain, 'first_gain_cpu_seconds': first_gain_cpu,
              'later_cloud_directions': len(packet['points'])-initial-first_cloud_gain,
              'calls': calls, 'epochs': epochs, 'unknowns': unknowns,
              'components': components, 'working_phase_cpu_seconds': search_end,
              'events_sha256': sha(dest/'events.json'), 'final_rank_sha256': sha(dest/'final-rank.json'),
              'protocol_sha256': protocol_hash}
    from verify_cancellation_basis_early import verify_arm
    verification = verify_arm(case, arm, dest, plan, result)
    write(dest/'independent-verification.json', verification)
    result.update(independent_verification_sha256=sha(dest/'independent-verification.json'),
                  final_replay_cpu_seconds=cpu()-search_end, total_internal_cpu_seconds=cpu())
    write(dest/'result.json', result)
    print(json.dumps({k: result[k] for k in ('case', 'arm', 'status', 'new_directions', 'later_cloud_directions', 'unknowns', 'total_internal_cpu_seconds')}), flush=True)

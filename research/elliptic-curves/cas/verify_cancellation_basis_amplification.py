#!/usr/bin/env sage -python
"""Replay basis/box transitions, full clouds, and independent rank certificates."""
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import time

from finite_cancellation_corpus import ROOT, canonical, digest, write

CAS = Path(__file__).resolve().parent


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify_arm(case, arm, dest, plan, result):
    from cancellation_basis_epoch import native_rank
    from cancellation_basis_exploration import Exposure, prepare
    from future_point_admission import FinitePointAdmission
    from half_lattice_pointed_sieve import linear_combination
    from pari_pointed_backend import replay
    from pointed_quartic_search import PointedQuarticSearch
    from run_complement_seed_v3 import normalized_selection
    from v3_warm_engine import certified_state
    tick = time.process_time(); seed = case['seed']; current = seed
    curve = tuple(map(F, seed['curve'])); initial = len(seed['points'])
    admission = FinitePointAdmission(curve, [tuple(map(F, p)) for p in seed['points']], prime_bound=plan['prime_bound'])
    initial_rank = native_rank(seed)
    need(initial_rank == read(dest/'initial-rank.json')['rank'], 'initial independent rank differs')
    mapper = SourceFileLoader('basis_replay_mapper', str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000, silent=True)
    exposure = Exposure(); epoch = -1; bank = None; pending_bank = False
    next_anchor = 0; next_model = 0; prepared = None; anchor = None; word = None
    calls = 0; gains = 0; first_gain = None; unknowns = 0; banks = 0; duplicates = 0
    required_refresh = False
    events = read(dest/'events.json')
    need(sha(dest/'events.json') == result['events_sha256'], 'event bytes differ')
    for event in events:
        need(gains < plan['target_directions'], 'events after completed target')
        kind = event['kind']
        if required_refresh:
            need(kind == 'bank_start', 'search continued on obsolete basis after gain')
        if kind == 'bank_start':
            need(not pending_bank and event['epoch'] == epoch+1, 'epoch order differs')
            need(epoch < 0 or arm == 'basis_refresh' and required_refresh, 'unjustified rebuild')
            epoch += 1; pending_bank = True; required_refresh = False
            need(event['rank'] == len(current['points']) and event['seed_sha256'] == digest(canonical(current)), 'bank does not use complete cloud')
        elif kind == 'bank_ready':
            need(pending_bank and event['epoch'] == epoch, 'bank receipt order differs')
            folder = dest/f'epoch-{epoch:02d}'; receipt = read(folder/'verification.json')
            need(sha(folder/'verification.json') == event['verification_sha256'], 'bank receipt changed')
            need(receipt['status'] == 'PASS_INDEPENDENT_BANK', 'unverified bank')
            for filename, field in (('bank.json', 'bank_sha256'), ('seed.json', 'seed_sha256'),
                                    ('producer/selection.json', 'producer_sha256'), ('reference/selection.json', 'reference_sha256')):
                need(sha(folder/filename) == receipt[field], 'landscape bytes differ')
            need(receipt['parent_bank_sha256'] == digest(canonical(case['parent_bank'])), 'parent subset changed')
            bank = read(folder/'bank.json'); producer = read(folder/'producer/selection.json'); reference = read(folder/'reference/selection.json')
            need(normalized_selection(producer) == normalized_selection(reference), 'independent rational CVP differs')
            need(bank['seed'] == current and producer['basis'] == current['points'] == reference['basis'], 'incompatible enlarged bank')
            need(bank['centres'] == producer['centres'][:plan['maximum_centres']], 'anchor selection differs')
            basis = [tuple(map(F, p)) for p in current['points']]
            state = certified_state(curve, basis, current['proof'])
            pending_bank = False; next_anchor = 0; next_model = 0; prepared = None; banks += 1
        elif kind == 'bank_unknown':
            need(pending_bank and event['epoch'] == epoch, 'unexpected bank failure')
            pending_bank = False; bank = None; unknowns += 1
        elif kind in ('prepared', 'preparation_unknown'):
            need(bank is not None and not pending_bank and event['epoch'] == epoch, 'preparation before compatible bank')
            need(event['centre'] == next_anchor, 'anchor order differs')
            need(prepared is None or next_model == len(prepared['models']), 'neighbour exploration silently skipped')
            word = bank['centres'][next_anchor]['representative']; next_anchor += 1
            if kind == 'preparation_unknown':
                unknowns += 1; prepared = None; continue
            path = dest/event['file']; need(sha(path) == event['sha256'], 'prepared map changed')
            prepared = read(path)
            need(prepared == prepare(curve, basis, word, mapper), 'model-independent preparation differs')
            anchor = list(map(str, linear_combination(curve, basis, word))); next_model = 0
        elif kind in ('call', 'exact_completed_duplicate'):
            need(prepared is not None and bank is not None and not pending_bank, 'call before preparation')
            model = prepared['models'][next_model]; mapping = model['mapping']
            if kind == 'exact_completed_duplicate':
                need(event['epoch'] == epoch and event['centre'] == next_anchor-1 and event['model'] == next_model, 'duplicate schedule differs')
                need(exposure.seen(anchor, mapping, plan['height']), 'non-exact box pruned')
                duplicates += 1; next_model += 1; continue
            row = read(dest/event['file']); record = result['calls'][calls]
            need(record == {k: event[k] for k in record} and sha(dest/event['file']) == event['sha256'], 'call receipt differs')
            need((row['epoch'], row['centre'], row['model']) == (epoch, next_anchor-1, next_model), 'call order differs')
            need(row['mapping'] == mapping and row['word'] == word and row['search']['height_bound'] == plan['height'], 'box binding differs')
            need(not exposure.seen(anchor, mapping, plan['height']), 'completed duplicate searched')
            witness = exposure.inspect(anchor, mapping, prepared['models'][0]['mapping'], plan['height'])
            need(row['tail_witness'] == witness and event['tail_witness'] == (witness is not None), 'tail exposure evidence differs')
            search = PointedQuarticSearch(state=state, centre={'coefficients': word}, coordinate_policy=mapping['coordinate_policy'])
            returned = replay(search, mapping, row['search']); before = len(admission.points); decisions = []
            for point in returned:
                decisions.append({'point': list(map(str, point)), **admission.consider(point)})
            gain = len(admission.points)-before
            need(decisions == row['admission'] and gain == row['gain'] == event['gain'], 'whole-cloud reconciliation differs')
            need(row['rank_after'] == len(admission.points), 'rank trace differs')
            if gain:
                path = dest/row['rank_file']; need(sha(path) == row['rank_sha256'], 'cloud certificate changed')
                packet = read(path); current = packet['packet']
                need(current['curve'] == seed['curve'] and current['points'] == [list(map(str, p)) for p in admission.points], 'cloud certificate drops or changes points')
                need(native_rank(current) == packet['independent_rank'], 'independent cloud rank differs')
                if first_gain is None: first_gain = gain
            else:
                need(row['rank_file'] is None and row['rank_sha256'] is None, 'unexpected rank packet')
            gains += gain
            required_refresh = bool(gain and arm == 'basis_refresh' and gains < plan['target_directions'])
            exposure.observe(anchor, mapping, plan['height'], row['search']['status'] == 'bounded_search_complete')
            calls += 1; next_model += 1
        else:
            raise ArithmeticError('unrecognized event: '+kind)
    need(not pending_bank, 'unreceipted bank')
    final = read(dest/'final-rank.json')
    need(sha(dest/'final-rank.json') == result['final_rank_sha256'] and final['packet'] == current, 'final subgroup differs')
    need(native_rank(current) == final['independent_rank'], 'final native rank differs')
    need(calls == len(result['calls']) and banks == len(result['epochs']) and unknowns == result['unknowns'], 'accounting differs')
    need(gains == result['new_directions'] and initial+gains == result['rank_lower_bound'], 'uncapped reward differs')
    need(result['first_cloud_gain'] == (first_gain or 0) and result['later_cloud_directions'] == gains-(first_gain or 0), 'later gain accounting differs')
    need(result['success'] == (gains >= plan['target_directions']), 'target status differs')
    if required_refresh:
        need(result['status'] in ('CPU_CAP', 'CPU_CAP_BEFORE_REFRESH'), 'unperformed refresh lacks budget boundary')
    return {'status': 'PASS', 'calls': calls, 'epochs': banks, 'exact_duplicates': duplicates,
            'rank': final['independent_rank'], 'new_directions': gains, 'unknowns': unknowns,
            'tail_witness_calls': sum(r['tail_witness'] for r in result['calls']),
            'cpu_seconds': time.process_time()-tick,
            'boundary': 'All full clouds, exact maps, box choices and enlarged subgroup transitions replay. Rational CVP comparison precedes every bank use; independent Sage quotient arithmetic certifies each rank increase.'}


def check():
    from cancellation_basis_amplification import OUT, RAW, guard
    plan, inputs = guard(); supervision = read(OUT/'supervision.json'); rows = []
    need(supervision['status'] == 'COMPLETE' and len(supervision['records']) == 2*plan['cases'], 'incomplete comparison')
    cases = {c['id']: c for c in inputs}
    for receipt in supervision['records']:
        dest = RAW/'arms'/receipt['case']/receipt['arm']; result = read(dest/'result.json')
        need(receipt['status'] == 'COMPLETE' and sha(dest/'result.json') == receipt['result_sha256'], 'process receipt differs')
        need(result['protocol_sha256'] == receipt['protocol_sha256'] == sha(OUT/'protocol.json'), 'protocol differs')
        need(sha(dest/'independent-verification.json') == result['independent_verification_sha256'], 'independent receipt differs')
        verification = read(dest/'independent-verification.json')
        need(verification['status'] == 'PASS' and verification['new_directions'] == result['new_directions'], 'uncertified result')
        need(receipt['charged_cpu_seconds'] >= result['total_internal_cpu_seconds']-.05, 'CPU accounting understates worker')
        need(sha(dest/'events.json') == result['events_sha256'], 'events changed')
        for call in result['calls']:
            need(sha(dest/call['file']) == call['sha256'], 'call changed')
        rows.append({'case': receipt['case'], 'arm': receipt['arm'], 'new_directions': result['new_directions'],
                     'charged_cpu_seconds': receipt['charged_cpu_seconds']})
    audit = {'status': 'PASS_RETAINED_BASIS_AMPLIFICATION_RECEIPTS', 'rows': rows,
             'protocol_sha256': sha(OUT/'protocol.json'), 'point_search_calls_added': 0}
    write(OUT/'audit.json', audit)
    print(json.dumps({'status': audit['status'], 'arms': len(rows)}), flush=True)


if __name__ == '__main__':
    check()

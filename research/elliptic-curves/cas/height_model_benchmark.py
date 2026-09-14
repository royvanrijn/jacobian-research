"""Frozen three-arm Curve302 experiment; no new selector or centre search.

The measured endpoint includes fresh-process map/rank replay. Target-location
diagnostics run separately, only after all arm endpoints have been sealed.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
import math
import os
from pathlib import Path
import sys
import time

from research_runtime.store import checkpoint, digest

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ARMS = ('v3_factor_free', 'bounded_single', 'bounded_pair')


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    from hashlib import sha256
    return sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def guard(folder):
    """Application read guard, not operating-system isolation."""
    roots = (folder.resolve(), ROOT.resolve())
    def audit(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes)):
            path = Path(args[0]).resolve()
            if 'artifacts' in path.parts and not any(path.is_relative_to(r) for r in roots):
                raise PermissionError('historical/sibling artifact read denied: '+str(path))
    sys.addaudithook(audit)


def scaled_height(height, ratio):
    """Least integer H with H^4 >= height^4 * ratio; exact ceiling."""
    ratio = F(ratio)
    require(0 < ratio <= 1, 'pair must match or improve the single coverage bound')
    lo, hi = 1, height
    while lo < hi:
        mid = (lo+hi)//2
        if mid**4*ratio.denominator >= height**4*ratio.numerator:
            hi = mid
        else:
            lo = mid+1
    return lo


def selection_from_portfolio(portfolio, packets):
    """Select only within the frozen original three-model vocabulary."""
    single = min(range(3), key=lambda i: (F(packets[i]['finite_gcd_divisor'])/
                                         F(packets[i]['real']['lower']), i))
    Dsingle = F(packets[single]['finite_gcd_divisor'])/F(packets[single]['real']['lower'])
    pairs = []
    for indices in combinations(range(3), 2):
        D = max(min(F(c['D_by_model'][i]) for i in indices) for c in portfolio['joint_cells'])
        pairs.append((D, indices))
    Dpair, pair = min(pairs)
    require(Dpair <= Dsingle, 'joint pair bound lost the best single model')
    order = sorted(pair, key=lambda i: (F(packets[i]['finite_gcd_divisor'])/
                                       F(packets[i]['real']['lower']), i))
    return {'single': single, 'pair': order, 'D_single': str(Dsingle),
            'D_pair': str(Dpair), 'single_height': 125000,
            'pair_height': scaled_height(125000, Dpair/Dsingle),
            'pair_squared_address_ratio': str(4*Dpair/Dsingle),
            'rule': 'Minimum uniform C/L single; minimum joint max-cell/min-model pair. '
                    'Ties by model index. Pair order by individual C/L. Exact ceiling height.'}


def prepare(folder, read_guard=True):
    if read_guard:
        guard(folder)
    from pointed_minimal_neighbours import build as neighbours
    from pointed_height_bounds import build
    from verify_pointed_height_bounds import verify
    from verify_pointed_height_portfolio import run as portfolio
    data = read(folder/'preconditioned_full-input.json')
    models, construction = neighbours(data)
    checkpoint(folder/'neighbours/construction.json', construction)
    require(len(models) == 2, 'frozen neighbour vocabulary incomplete')
    require(len({m['neighbour_witness']['prime'] for m in models}) == 2,
            'frozen independent-prime portfolio formula does not support this chart')
    paths = [folder/'preconditioned_full']
    for i, model in enumerate(models):
        path = folder/'neighbours'/f'neighbour-{i}'
        checkpoint(Path(str(path)+'-input.json'), model)
        paths.append(path)
    packets = []
    for data, path in zip([data]+models, paths):
        packet = build(data)
        checkpoint(Path(str(path)+'-bounds.json'), packet)
        checkpoint(Path(str(path)+'-verified.json'), verify(packet, ROOT))
        packets.append(packet)
    result = portfolio(folder, ROOT)
    checkpoint(folder/'portfolio.json', result)
    choice = selection_from_portfolio(result, packets)
    choice['models'] = [str(Path(str(p)+'-input.json').relative_to(folder)) for p in paths]
    choice['bounds'] = [str(Path(str(p)+'-bounds.json').relative_to(folder)) for p in paths]
    checkpoint(folder/'selection.json', choice)
    print('PREPARED', choice['single'], choice['pair'], choice['pair_height'], flush=True)


def bound_input(seed, mapping, policy, protocol):
    return {'case': seed['id'], 'policy': policy, 'curve': seed['curve'],
            'mapping': mapping, 'search_height': 125000,
            'certified_primes': protocol['certified_primes'],
            'prime_certificate': protocol['prime_certificate'], 'source_bindings': {}}


def support(folder):
    """Cold support discovery from E alone, bounded externally and fully charged."""
    guard(folder)
    from cypari2 import Pari
    from sage.all import ZZ, prime_range
    curve=read(folder/'input.json')['payload']['curve']
    A,B=map(F,curve[3:])
    discriminant=-16*(4*A**3+27*B**2)/6**12
    require(discriminant.denominator==1 and discriminant>0, 'unsupported minimal discriminant normalization')
    pari=Pari();pari.allocatemem(256000000,silent=True)
    factors=pari(str(discriminant.numerator)).factor()
    values=[(int(factors[i,0]),int(factors[i,1])) for i in range(factors.nrows())]
    require(math.prod(p**e for p,e in values)==discriminant, 'cold factor product differs')
    for p,e in values:
        require(e>0 and ZZ(p).is_prime(proof=True), 'unproved cold prime factor')
    primes=sorted(set(p for p,e in values)|set(map(int,prime_range(2,168))))
    data={'schema':'height-model-cold-discriminant-support.v1', 'proved_primes':list(map(str,primes)),
          'records':{'302':{'factorizations':{'DISCRIMINANT_FACTORIZATION':[[str(p),e] for p,e in values]}}},
          'curve':curve,'method':'Fresh PARI integer factorization followed by Sage unconditional prime proofs; no factor hints.'}
    path=folder/'support/primes.json';checkpoint(path,data)
    checkpoint(folder/'support/receipt.json', {'certified_primes':primes,
        'prime_certificate':{'path':os.path.relpath(path,ROOT),'sha256':sha(path)}})
    print('COLD_SUPPORT_CERTIFIED',len(values),flush=True)


def worker(folder):
    guard(folder)
    from memory_rank_certificate import checked_rank
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from lean_preconditioned_map_receipts import obtain
    from research_runtime.supervisor import run, Limits
    import pari_pointed_backend as backend
    started = time.monotonic()
    protocol, raw = read(folder/'protocol.json'), read(folder/'input.json')
    require(raw['sha256'] == digest(raw['payload']), 'changed input projection')
    seed = raw['payload']; arm = protocol['arm']
    model = tuple(map(F, seed['curve']))
    points = tuple(tuple(map(F, p)) for p in seed['points'])
    proof = checked_rank(model, points, seed['primes'], seed['torsion_prime'])
    checkpoint(folder/'initial-rank.json', proof)
    state = certified_state(model, points, proof)
    admission = FinitePointAdmission(model, points, prime_bound=1000)
    rows, gains = [], False
    status = 'NO_GAIN_IN_COMPLETED_BANK'
    if arm != ARMS[0]:
        checkpoint(folder/'progress.json', {'status':'COLD_DISCRIMINANT_PREPARATION','point_calls':0})
        outcome=run([protocol['map_python'],'-python',str(Path(__file__)),'support','--folder',str(folder)],
                    limits=Limits(protocol['support_seconds'],protocol['map_rss_bytes']),
                    log_path=folder/'support.log',checkpoint_path=folder/'support-supervisor.json')
        if outcome['outcome']!='completed' or outcome['returncode']!=0:
            checkpoint(folder/'result.json',{'status':'COLD_SUPPORT_UNRESOLVED','initial_rank':len(points),
                'rows':[],'arm':arm,'support_outcome':outcome['outcome'],'returncode':outcome['returncode'],
                'input_sha256':sha(folder/'input.json'),'worker_wall_seconds':time.monotonic()-started})
            return
        protocol.update(read(folder/'support/receipt.json'))
        checkpoint(folder/'bound-data-protocol.json',protocol)
    for ci, centre in enumerate(seed['centres']):
        if time.monotonic()-started >= protocol['search_phase_seconds']:
            status = 'CENSORED_ARM_ALLOWANCE'; break
        key = 'factor_free' if arm == ARMS[0] else 'preconditioned_full'
        mapping, receipt, limited = obtain(folder, ci, key, model, points, centre, state, protocol)
        if limited:
            rows.append({'centre': ci, 'status': limited, 'map_receipt': receipt}); continue
        if arm == ARMS[0]:
            work = [(key, mapping, 125000, None)]
        else:
            bank = folder/f'bank-{ci:04d}'
            checkpoint(bank/'preconditioned_full-input.json', bound_input(seed, mapping, key, protocol))
            remaining = protocol['search_phase_seconds']-(time.monotonic()-started)
            if remaining < 2:
                status = 'CENSORED_ARM_ALLOWANCE'; break
            record = run([protocol['map_python'], '-python', str(Path(__file__)),
                          'prepare', '--folder', str(bank)],
                         limits=Limits(min(20, remaining), protocol['map_rss_bytes']),
                         log_path=bank/'prepare.log', checkpoint_path=bank/'supervisor.json')
            if record['outcome'] != 'completed' or record['returncode'] != 0:
                rows.append({'centre': ci, 'status': 'INCOMPLETE_MODEL_SELECTION',
                             'preparation_outcome': record['outcome'], 'returncode': record['returncode'],
                             'bank': bank.name})
                # A failed frozen formula is a result, not an invitation to repair it mid-arm.
                continue
            choice = read(bank/'selection.json')
            indices = [choice['single']] if arm == ARMS[1] else choice['pair']
            height = choice['single_height'] if arm == ARMS[1] else choice['pair_height']
            work = [(str(i), read(bank/choice['models'][i])['mapping'], height,
                     str((bank/choice['bounds'][i]).relative_to(folder))) for i in indices]
        for key, mapping, height, bound in work:
            if time.monotonic()-started >= protocol['search_phase_seconds']:
                status = 'CENSORED_ARM_ALLOWANCE'; break
            search = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                         coordinate_policy=mapping['coordinate_policy'])
            transcript, returned = backend.execute(search, mapping, height,
                                                   protocol['seconds_per_search'], protocol['gp_sha256'])
            require(backend.replay(search, mapping, transcript) == returned, 'point-map replay differs')
            name = f'chart-{len(rows):04d}.json'
            checkpoint(folder/name, {'centre_index': ci, 'mapping': mapping, 'search': transcript,
                                     'bound': bound})
            row = {'centre': ci, 'model': key, 'height': height, 'chart': name,
                   'sha256': sha(folder/name), 'status': transcript['status'], 'admissions': {}}
            rows.append(row)
            for point in returned:
                outcome = admission.consider(point)['status']
                row['admissions'][outcome] = row['admissions'].get(outcome, 0)+1
                if outcome == 'INDEPENDENT_FINITE_COLUMN':
                    proof = checked_rank(model, admission.points, admission.primes, seed['torsion_prime'])
                    checkpoint(folder/'gain.json', {'curve': seed['curve'],
                        'points': [list(map(str, p)) for p in admission.points], 'proof': proof,
                        'centre': ci, 'chart': name, 'chart_sha256': sha(folder/name)})
                    status = 'GAIN_PENDING_INDEPENDENT_REPLAY'; gains = True; break
            checkpoint(folder/'progress.json', {'status': status, 'centre': ci,
                'point_calls': sum('chart' in r for r in rows), 'rows': rows,
                'worker_wall_seconds': time.monotonic()-started})
            if gains: break
        if gains or status == 'CENSORED_ARM_ALLOWANCE': break
    if status == 'NO_GAIN_IN_COMPLETED_BANK' and any(r['status'] != 'bounded_search_complete' for r in rows):
        status = 'NO_GAIN_WITH_INCOMPLETE_CHARTS'
    checkpoint(folder/'result.json', {'status': status, 'initial_rank': len(points),
        'rows': rows, 'arm': arm, 'input_sha256': sha(folder/'input.json'),
        'worker_wall_seconds': time.monotonic()-started})


def replay(folder):
    guard(folder)
    from memory_rank_certificate import checked_rank
    from v3_warm_engine import certified_state, load
    from pointed_quartic_search import PointedQuarticSearch
    from high_rank_foundry_certificate import compact_packet
    import pari_pointed_backend as backend
    seed, result = read(folder/'input.json')['payload'], read(folder/'result.json')
    model = tuple(map(F, seed['curve'])); points = tuple(tuple(map(F, p)) for p in seed['points'])
    proof = checked_rank(model, points, seed['primes'], seed['torsion_prime'])
    require(digest(proof) == digest(read(folder/'initial-rank.json')), 'starting rank replay differs')
    state = certified_state(model, points, proof)
    outputs = {}
    for row in result['rows']:
        if 'chart' not in row: continue
        require(sha(folder/row['chart']) == row['sha256'], 'changed chart')
        data = read(folder/row['chart']); centre = seed['centres'][row['centre']]
        search = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                     coordinate_policy=data['mapping']['coordinate_policy'])
        outputs[row['chart']] = backend.replay(search, data['mapping'], data['search'])
    gained = (folder/'gain.json').exists()
    target = read(folder/'gain.json') if gained else {'curve': seed['curve'], 'points': seed['points'], 'proof': proof}
    if gained:
        new = tuple(tuple(map(F, p)) for p in target['points'])
        require(new[:-1] == points and new[-1] in outputs[target['chart']], 'gain not from declared search')
        require(target['chart_sha256'] == sha(folder/target['chart']), 'gain chart changed')
        actual = checked_rank(model, new, [s['prime'] for s in target['proof']['signatures']], seed['torsion_prime'])
        require(digest(actual) == digest(target['proof']), 'gain certificate differs')
    compact = compact_packet({'curve': target['curve'], 'points': target['points'],
                              'proof': target['proof'], 'rank_lower_bound': len(target['points'])})
    packet = {**compact, 'signatures': compact['proof']['signatures'],
              'independent_column_indices': list(range(len(target['points']))),
              'rank_certificate': compact['proof']}
    checkpoint(folder/'independent-rank-input.json', packet)
    checker = load('independent_sage_finite_groups', CAS/'verify_factor_free_rank.sage')
    independent = checker.check_one(folder/'independent-rank-input.json')
    checkpoint(folder/'verified.json', {'status': 'PASS_POINT_MAP_AND_TWO_FINITE_IMPLEMENTATIONS',
        'independent_gain_verified': gained, 'rank_lower_bound': len(target['points']),
        'result_sha256': sha(folder/'result.json'), 'sage_enumeration': independent})


def diagnose(folder):
    """All-arm sealed output evaluation; never imported by search selection."""
    guard(folder)
    from memory_rank_certificate import checked_rank
    from v3_warm_engine import certified_state
    from lean_preconditioned_map_receipts import obtain
    from pointed_quartic_search import PointedQuarticSearch
    from pointed_height_bounds import build
    from verify_pointed_height_bounds import verify
    from search_observability import point_visibility
    from half_lattice_pointed_sieve import linear_combination
    from sage.all import EllipticCurve, QQ
    sealed = read(folder/'all-arms-sealed.json')
    rows = []
    for entry in sealed['rows']:
        if not entry['success']: continue
        arm = folder/entry['folder']; seed = read(arm/'input.json')['payload']
        gain = read(arm/'gain.json'); protocol = read(arm/'protocol.json')
        ci = gain['centre']; target = tuple(map(F, gain['points'][-1]))
        out = folder/'diagnostics'/seed['id']/protocol['arm']; out.mkdir(parents=True, exist_ok=False)
        model = tuple(map(F, seed['curve'])); points = tuple(tuple(map(F, p)) for p in seed['points'])
        proof = checked_rank(model, points, seed['primes'], seed['torsion_prime'])
        state = certified_state(model, points, proof); centre = seed['centres'][ci]
        maps = []
        for policy in ('factor_free', 'preconditioned_full'):
            mapping, _, limited = obtain(out, ci, policy, model, points, centre, state, protocol)
            if limited:
                rows.append({'source_arm': entry['folder'], 'status': 'DIAGNOSTIC_MAP_TIMEOUT', 'model': policy}); continue
            data = bound_input(seed, mapping, policy, protocol)
            if policy == 'preconditioned_full':
                checkpoint(out/'preconditioned_full-input.json', data)
                prepare(out, read_guard=False)
                choice = read(out/'selection.json')
                maps.extend((read(out/p), read(out/b)) for p,b in zip(choice['models'], choice['bounds']))
            else:
                packet = build(data); verify(packet, ROOT)
                checkpoint(out/'factor_free-bounds.json', packet); maps.append((data, packet))
        require(len(maps) == 4, 'incomplete diagnostic model set')
        anchor = tuple(map(F, maps[0][1]['anchor']))
        R = linear_combination(model, (target, anchor), (2, -1))
        E = EllipticCurve(QQ, list(map(str, model)))
        sageR = 2*E(list(map(str,target)))-E(list(map(str,anchor)))
        require((R is None and sageR.is_zero()) or
                (R is not None and tuple(map(str,R)) == tuple(str(c) for c in sageR.xy())), 'independent 2P-Q transport differs')
        Hx = 1 if R is None else max(abs(R[0].numerator), R[0].denominator)
        point_rows = []
        for i, (data, packet) in enumerate(maps):
            search = PointedQuarticSearch(state=state, centre={'coefficients':centre['representative']},
                                         coordinate_policy=data['mapping']['coordinate_policy'])
            location = point_visibility(search.chart_record(), target)
            n,d = map(int, location['coordinate']); H = max(abs(n), abs(d))
            D = F(packet['finite_gcd_divisor'])/F(packet['real']['lower'])
            ratio = F(H**4, Hx)/D
            require(ratio <= 1, 'winning point violates certified sufficient bound')
            point_rows.append({'model':data['policy'], 'coordinate':[str(n),str(d)],
                'exact_parameter_height': str(H), 'D':str(D), 'actual_over_sufficient_fourth_power':str(ratio),
                'h_parameter_minus_quarter_hx': math.log(H)-math.log(Hx)/4,
                'B2': math.log(D.numerator)/4-math.log(D.denominator)/4,
                'inside_height_125000': H<=125000,
                'inside_pair_box': H<=choice['pair_height'] if i else None,
                'parameter_infinity': d==0})
        pair = [point_rows[i+1] for i in choice['pair']]
        best_actual = min(F(r['exact_parameter_height'])**4 for r in pair)
        require(best_actual <= F(choice['D_pair'])*Hx, 'joint sufficient bound fails on actual target')
        result = {'source_arm':entry['folder'], 'centre':ci, 'target':gain['points'][-1],
                  'R_2P_minus_Q':None if R is None else list(map(str,R)), 'exact_x_height':str(Hx),
                  'models':point_rows, 'selection':choice,
                  'pair_actual_over_sufficient_fourth_power':str(best_actual/(F(choice['D_pair'])*Hx)),
                  'single_guarantees_target': F(Hx) <= F(125000**4)/F(choice['D_single']),
                  'pair_guarantees_target': F(Hx) <= F(choice['pair_height']**4)/F(choice['D_pair']),
                  'best_single_accessible': int(point_rows[choice['single']+1]['exact_parameter_height'])<=125000,
                  'pair_accessible': any(int(r['exact_parameter_height'])<=choice['pair_height'] for r in pair),
                  'boundary':'Literal sealed winning point and its winning anchor, transported into every model. '
                             'Different arm winners need not be the same rational point or anchor. Post-seal cost separate.'}
        checkpoint(out/'target-accessibility.json', result); rows.append(result)
    checkpoint(folder/'target-accessibility.json', {'status':'PASS_SEALED_OUTPUT_ACCESSIBILITY', 'rows':rows})


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('prepare','support','worker','replay','diagnose'))
    p.add_argument('--folder', type=Path, required=True)
    a=p.parse_args();globals()[a.mode](a.folder.resolve())

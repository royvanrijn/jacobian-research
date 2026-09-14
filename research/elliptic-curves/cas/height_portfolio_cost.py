"""Exact box allocation over an existing certified joint-state table.

This finite postprocessor constructs no models or bounds and performs no point
search. Costs are supplied estimates (the retained example uses H^2 units).
"""
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(value):
    require(type(value) in (int, str, Fraction), 'costs and bounds must be exact rationals')
    return Fraction(value)


def ceil_fourth_root(value):
    value = exact(value)
    require(value > 0, 'a height threshold must be positive')
    lo, hi = 0, 1
    while hi**4 * value.denominator < value.numerator:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**4 * value.denominator >= value.numerator:
            hi = mid
        else:
            lo = mid
    return hi


def instance(table, x_height, activation_costs, address_costs, common_cost=0):
    rows = [[exact(d) for d in row] for row in table]
    require(rows and rows[0], 'empty state/model bank')
    count = len(rows[0])
    require(all(len(row) == count and all(d > 0 for d in row) for row in rows),
            'invalid state table')
    X = exact(x_height)
    fixed = list(map(exact, activation_costs))
    rates = list(map(exact, address_costs))
    common = exact(common_cost)
    require(X >= 1, 'elliptic multiplicative height limit must be at least one')
    require(len(fixed) == len(rates) == count, 'model/cost count differs')
    require(common >= 0 and all(a >= 0 for a in fixed) and all(c > 0 for c in rates),
            'cost model must be monotone, with positive address rates')
    thresholds = [[ceil_fourth_root(X * d) for d in row] for row in rows]
    return thresholds, fixed, rates, common


def allocation(heights, thresholds, fixed, rates, common):
    covered = [[i for i, H in enumerate(heights) if H and H >= row[i]]
               for row in thresholds]
    if not all(covered):
        return None
    active = [i for i, H in enumerate(heights) if H]
    selected_cost = sum((fixed[i] + rates[i] * heights[i]**2 for i in active), Fraction(0))
    return {'heights': list(heights), 'active_models': active,
            'state_covering_models': covered, 'state_witness_models': [r[0] for r in covered],
            'selected_cost': str(selected_cost), 'objective': str(common + selected_cost)}


def choice_key(row):
    return (Fraction(row['objective']), len(row['active_models']), row['heights'])


def solve(table, x_height, activation_costs, address_costs, common_cost=0, node_cap=125):
    thresholds, fixed, rates, common = instance(
        table, x_height, activation_costs, address_costs, common_cost)
    count = len(fixed)
    options = [sorted({0, *(row[i] for row in thresholds)}) for i in range(count)]
    nodes = prod(map(len, options))
    require(nodes <= node_cap, 'declared allocation cap exceeded; no automatic enlargement')
    feasible = []
    for heights in product(*options):
        result = allocation(heights, thresholds, fixed, rates, common)
        if result is not None:
            feasible.append(result)
    require(feasible, 'no covered allocation')
    best_by_size = {}
    uniform_by_size = {}
    for size in range(1, count + 1):
        eligible = [r for r in feasible if len(r['active_models']) <= size]
        if eligible:
            best_by_size[str(size)] = min(eligible, key=choice_key)
        uniform = []
        for indices in combinations(range(count), size):
            H = max(min(row[i] for i in indices) for row in thresholds)
            result = allocation([H if i in indices else 0 for i in range(count)],
                                thresholds, fixed, rates, common)
            require(result is not None, 'uniform construction failed coverage')
            uniform.append(result)
        uniform_by_size[str(size)] = min(uniform, key=choice_key)
    return {'thresholds': thresholds, 'height_options': options, 'enumeration_nodes': nodes,
            'feasible_allocations': len(feasible), 'optimum': min(feasible, key=choice_key),
            'best_by_maximum_model_count': best_by_size,
            'best_uniform_by_exact_model_count': uniform_by_size}


def load_protocol(folder):
    protocol = json.loads((folder / 'protocol.json').read_text())
    source = ROOT / protocol['state_table']['path']
    require(hashlib.sha256(source.read_bytes()).hexdigest() == protocol['state_table']['sha256'],
            'state table bytes changed')
    for path, sha in protocol['source_lock'].items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha,
                'locked implementation changed: ' + path)
    packet = json.loads(source.read_text())
    require(packet['status'] == 'PASS_INDEPENDENT_LOCAL_PORTFOLIO_REPLAY',
            'table has no retained portfolio replay')
    require(packet['model_order'] == protocol['model_order'], 'model ordering changed')
    table = [[Fraction(d) for d in row['D_by_model']] for row in packet['joint_cells']]
    best_single = min(max(row[i] for row in table) for i in range(len(table[0])))
    require(best_single == Fraction(packet['best_single']['D']), 'single-model bound changed')
    reference = protocol['height_target']
    require(reference['rule'] == 'same certified range as best single at reference box',
            'unsupported height selection rule')
    H = reference['reference_height']
    require(type(H) is int and H > 0, 'invalid reference height')
    return protocol, packet, table, Fraction(H**4) / best_single


def run(folder):
    protocol, packet, table, X = load_protocol(folder)
    cost = protocol['cost_model']
    result = solve(table, X, cost['activation'], cost['per_H_squared'],
                   cost['common'], protocol['allocation_node_cap'])
    result.update(status='EXACT_COST_ALLOCATION_PRODUCED', x_height_limit=str(X),
                  model_order=packet['model_order'],
                  protocol_sha256=hashlib.sha256((folder / 'protocol.json').read_bytes()).hexdigest(),
                  boundary=protocol['boundary'])
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'preserve previous output')
    start, elapsed = time.process_time(), time.monotonic()
    result = run(args.folder.resolve())
    result['meter'] = {'cpu_seconds': time.process_time() - start,
                       'elapsed_seconds': time.monotonic() - elapsed,
                       'scope': 'postprocessor only; excludes input mathematics, preparation and search'}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'nodes': result['enumeration_nodes'],
                      'optimum': result['optimum'], 'meter': result['meter']}, indent=2))

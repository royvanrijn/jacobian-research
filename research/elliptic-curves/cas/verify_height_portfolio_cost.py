"""Check allocation optimality by independent assignments of states to models.

Unlike the producer's box-product enumeration, this enumerates m^s assignments.
The original mathematical bound/transport checker replays the supplied table.
No point search, new reduction, factorization or benchmark is performed.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def assignment_minima(thresholds, fixed, rates, common):
    """Return minima by *at most* k models, using a different finite reduction."""
    models, states = len(fixed), len(thresholds)
    require(models**states <= 125, 'independent assignment cap exceeded')
    minima = {}
    for assignment in product(range(models), repeat=states):
        H = [0] * models
        for s, i in enumerate(assignment):
            H[i] = max(H[i], thresholds[s][i])
        cost = common + sum((fixed[i] + rates[i]*h*h for i, h in enumerate(H) if h), F(0))
        used = sum(h > 0 for h in H)
        for k in range(used, models + 1):
            if k not in minima or cost < minima[k]:
                minima[k] = cost
    return minima


def verify(folder):
    p = json.loads((folder / 'protocol.json').read_text())
    r = json.loads((folder / 'allocation.json').read_text())
    require(r['status'] == 'EXACT_COST_ALLOCATION_PRODUCED', 'unexpected producer status')
    require(r['protocol_sha256'] == hashlib.sha256((folder / 'protocol.json').read_bytes()).hexdigest(),
            'protocol hash differs')
    for path, sha in p['source_lock'].items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == sha, 'source binding differs')
    source = ROOT / p['state_table']['path']
    require(hashlib.sha256(source.read_bytes()).hexdigest() == p['state_table']['sha256'],
            'state table binding differs')
    original = json.loads(source.read_text())
    # The legacy four-cell construction assumes two distinct neighbour primes.
    primes = [row['prime'] for row in original['partitions']]
    require(len(primes) == len(set(primes)) == 2, 'distinct-prime table prerequisite absent')
    from verify_pointed_height_portfolio import run as replay_table
    upstream = replay_table(source.parent, ROOT)
    for key in ('status', 'packet_bindings', 'model_order', 'partitions', 'joint_cells', 'subsets', 'best_single'):
        require(upstream[key] == original[key], 'upstream table replay differs: ' + key)
    models = len(original['model_order'])
    require(r['model_order'] == p['model_order'] == original['model_order'], 'model order differs')
    D = [[F(d) for d in row['D_by_model']] for row in original['joint_cells']]
    require(p['height_target']['rule'] == 'same certified range as best single at reference box',
            'unsupported height rule')
    reference_height = p['height_target']['reference_height']
    require(type(reference_height) is int and reference_height > 0, 'invalid reference height')
    X = F(reference_height**4) / min(
        max(row[i] for row in D) for i in range(models))
    require(r['x_height_limit'] == str(X), 'height target differs')
    thresholds = r['thresholds']
    require(len(thresholds) == len(D), 'state count differs')
    for ds, hs in zip(D, thresholds):
        require(len(hs) == models, 'model count differs')
        for d, h in zip(ds, hs):
            require(type(h) is int and h > 0 and (h-1)**4 < X*d <= h**4,
                    'incorrect exact threshold')
    costs = p['cost_model']
    fixed, rates, common = list(map(F, costs['activation'])), list(map(F, costs['per_H_squared'])), F(costs['common'])
    require(len(fixed) == len(rates) == models and common >= 0 and all(a >= 0 for a in fixed)
            and all(c > 0 for c in rates), 'invalid monotone costs')
    expected_options = [sorted({0, *(row[i] for row in thresholds)}) for i in range(models)]
    require(r['height_options'] == expected_options, 'finite breakpoints differ')
    nodes = 1
    for options in expected_options:
        nodes *= len(options)
    require(r['enumeration_nodes'] == nodes <= p['allocation_node_cap'], 'enumeration cap/count differs')

    def check(row):
        H = row['heights']
        require(len(H) == models and all(type(h) is int and h >= 0 for h in H), 'invalid box')
        active = [i for i, h in enumerate(H) if h]
        covers = [[i for i in active if H[i]**4 >= X*ds[i]] for ds in D]
        require(all(covers) and covers == row['state_covering_models'], 'incomplete or wrong coverage')
        require(row['active_models'] == active and row['state_witness_models'] == [v[0] for v in covers],
                'coverage witness differs')
        selected = sum((fixed[i] + rates[i]*H[i]**2 for i in active), F(0))
        require(F(row['selected_cost']) == selected and F(row['objective']) == common+selected,
                'allocation cost differs')
        return common + selected

    minima = assignment_minima(thresholds, fixed, rates, common)
    best = check(r['optimum'])
    require(best == minima[models], 'not an optimal allocation')
    entries = r['best_by_maximum_model_count']
    require(set(entries) == {str(k) for k in range(1, models+1)}, 'incomplete constrained comparisons')
    for k in range(1, models+1):
        require(len(entries[str(k)]['active_models']) <= k, 'model limit exceeded')
        require(check(entries[str(k)]) == minima[k], 'constrained optimum differs')
    # Recheck every uniform comparison independently, by all nonempty bit masks.
    uniform = {}
    for mask in range(1, 1 << models):
        active = [i for i in range(models) if mask & (1 << i)]
        H = max(min(row[i] for i in active) for row in thresholds)
        cost = common + sum((fixed[i] + rates[i]*H*H for i in active), F(0))
        k = len(active)
        uniform[k] = min(uniform.get(k, cost), cost)
    require(set(r['best_uniform_by_exact_model_count']) == set(entries), 'incomplete uniform comparisons')
    for key, row in r['best_uniform_by_exact_model_count'].items():
        require(len(row['active_models']) == int(key), 'uniform model count differs')
        require(len(set(h for h in row['heights'] if h)) == 1, 'claimed uniform boxes differ')
        require(check(row) == uniform[int(key)], 'uniform optimum differs')
    return {'status': 'PASS_INDEPENDENT_ASSIGNMENT_AND_BOUND_REPLAY',
            'protocol_sha256': r['protocol_sha256'],
            'allocation_sha256': hashlib.sha256((folder / 'allocation.json').read_bytes()).hexdigest(),
            'state_count': len(D), 'model_count': models, 'box_allocations': nodes,
            'assignments_checked': models**len(D),
            'optimal_heights': r['optimum']['heights'], 'objective': str(best),
            'boundaries': ['Optimal only for the supplied sufficient-state table and supplied cost model.',
                           'The retained cost units are an address proxy; setup/search CPU is unmeasured.',
                           'No rational point, new direction, target accessibility or speedup is inferred.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    require(args.output is None or not args.output.exists(), 'preserve previous receipt')
    start, elapsed = time.process_time(), time.monotonic()
    result = verify(args.folder.resolve())
    result['meter'] = {'cpu_seconds': time.process_time()-start,
                       'elapsed_seconds': time.monotonic()-elapsed,
                       'scope': 'independent allocation and retained bound replay; no search or cold preparation'}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))

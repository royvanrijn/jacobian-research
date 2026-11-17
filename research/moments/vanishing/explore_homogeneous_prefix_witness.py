"""Bounded Newton exploration, never a proof of exact moment vanishing.

P=Z^3+b W^3+c Z^2 T+d W^2 T+e(Z^2W-2ZT^2)
  +f(ZW^2-2WT^2)+h(T^3-3ZWT)+ZWT has Delta P=2T.
Thus no member is pure-null, using the retained homogeneous classification.
The six equations are M2,...,M12. All moments are normalized by (3n+1)!!.
"""
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np

P = Path(__file__).resolve().parent
BASE = [((3, 0, 0), 1), ((1, 1, 1), 1)]
BASIS = [[((0, 3, 0), 1)], [((2, 0, 1), 1)], [((0, 2, 1), 1)],
         [((2, 1, 0), 1), ((1, 0, 2), -2)],
         [((1, 2, 0), 1), ((0, 1, 2), -2)],
         [((0, 0, 3), 1), ((1, 1, 1), -3)]]


def dfact(n):
    return math.prod(range(n, 0, -2))


def maps(n, shift=(0, 0, 0)):
    total = 3*n+sum(shift)
    assert total % 2 == 0
    denominator = dfact(total+1)
    entries = []
    for i in range(3*n+1):
        j = i+shift[0]-shift[1]
        k = 3*n-i-j
        if 0 <= j <= 3*n and k >= 0 and (k+shift[2]) % 2 == 0:
            weight = math.factorial(i+shift[0])*dfact(k+shift[2]-1)/denominator
            entries.append((i, j, weight))
    return (np.array([v[0] for v in entries], dtype=int),
            np.array([v[1] for v in entries], dtype=int),
            np.array([v[2] for v in entries]))


MOMENT_MAPS = {n: maps(n) for n in range(2, 15, 2)}
DERIVATIVE_MAPS = {(n, j): [(factor, maps(n-1, shift)) for shift, factor in family]
                   for n in range(2, 13, 2) for j, family in enumerate(BASIS)}


def evaluate(parameters, with_m14=False):
    terms = list(BASE)
    for value, family in zip(parameters, BASIS):
        terms.extend((shift, value*factor) for shift, factor in family)
    values = []
    jacobian = []
    poly = np.ones((1, 1), dtype=complex)
    for n in range(1, 15 if with_m14 else 13):
        if n % 2 == 0 and n <= 12:
            jacobian.append([n*sum(factor*np.dot(poly[mp[0], mp[1]], mp[2])
                                    for factor, mp in DERIVATIVE_MAPS[n, j])
                             for j in range(6)])
        size = 3*n+1
        next_poly = np.zeros((size, size), dtype=complex)
        for (i, j, _), coefficient in terms:
            old_size = poly.shape[0]
            next_poly[i:i+old_size, j:j+old_size] += coefficient*poly
        poly = next_poly
        if n % 2 == 0:
            mp = MOMENT_MAPS[n]
            values.append(np.dot(poly[mp[0], mp[1]], mp[2]))
    return np.array(values), np.array(jacobian)


def main():
    start = time.monotonic()
    rng = np.random.default_rng(20261115)
    attempts = []
    candidates = []
    output = P/'homogeneous_prefix_witness_exploration.json'
    report = {'completed': False, 'kind': 'numerical exploration only',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seed': 20261115, 'maximum_starts': 96, 'maximum_iterations': 90,
              'wall_clock_limit_seconds': 110,
              'family': 'Z^3+ZWT+b W^3+c Z^2 T+d W^2 T+e(Z^2W-2ZT^2)+f(ZW^2-2WT^2)+h(T^3-3ZWT)',
              'target_moments': [2, 4, 6, 8, 10, 12],
              'normalization': 'M_n/(3n+1)!!', 'attempts': attempts, 'candidates': candidates,
              'exact_witness_certified': False}

    def checkpoint():
        report['runtime_seconds'] = round(time.monotonic()-start, 3)
        output.write_text(json.dumps(report, indent=2)+'\n')

    for index in range(96):
        if time.monotonic()-start > 110:
            break
        scale = [0.15, 0.4, 1.0][index % 3]
        parameters = scale*(rng.normal(size=6)+1j*rng.normal(size=6))
        best = float('inf')
        reason = 'iterations'
        for iteration in range(90):
            values, jac = evaluate(parameters)
            if not np.all(np.isfinite(values)) or not np.all(np.isfinite(jac)):
                reason = 'nonfinite'; break
            norm = float(np.linalg.norm(values, ord=np.inf))
            best = min(best, norm)
            if norm < 2e-15:
                condition = float(np.linalg.cond(jac))
                all_values, _ = evaluate(parameters, with_m14=True)
                candidates.append({'start': index, 'parameters': [[float(v.real), float(v.imag)] for v in parameters],
                                   'normalized_M2_to_M14': [[float(v.real), float(v.imag)] for v in all_values],
                                   'Jacobian_condition': condition, 'exact': False})
                reason = 'candidate'; break
            row_scale = np.maximum(np.max(abs(jac), axis=1), 1e-30)
            try:
                step = np.linalg.solve(jac/row_scale[:, None], values/row_scale)
            except np.linalg.LinAlgError:
                reason = 'singular'; break
            if np.max(abs(step)) > 15:
                step *= 15/np.max(abs(step))
            accepted = False
            for damping in [1, .5, .25, .125, .0625, .03125, .015625]:
                trial = parameters-damping*step
                if np.max(abs(trial)) > 30:
                    continue
                residual, _ = evaluate(trial)
                if np.linalg.norm(residual/row_scale, ord=2) < np.linalg.norm(values/row_scale, ord=2):
                    parameters = trial; accepted = True; break
            if not accepted:
                reason = 'line_search'; break
        attempts.append({'start': index, 'iterations': iteration+1, 'best_absolute_residual': best, 'reason': reason})
        checkpoint()
        if len(candidates) >= 3:
            break
        if index % 8 == 0:
            print(json.dumps({'starts': index+1, 'candidates': len(candidates), 'elapsed': report['runtime_seconds']}), flush=True)
    report['completed'] = True
    report['terminal_reason'] = 'three_candidates' if len(candidates) >= 3 else ('wall_clock' if len(attempts) < 96 else 'start_budget')
    checkpoint()
    print(json.dumps({'completed': True, 'attempts': len(attempts), 'candidates': len(candidates),
                      'runtime_seconds': report['runtime_seconds'], 'exact_witness_certified': False}), flush=True)


if __name__ == '__main__':
    main()

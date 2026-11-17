"""Bounded numerical search for M2=...=M18=0 on an odd cubic slice.

This is only a candidate generator. A later interval certificate is required.
The cubic part has Laplacian2T, and three free linear terms are added to
the retained homogeneous six-parameter family.
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
         [((0, 0, 3), 1), ((1, 1, 1), -3)],
         [((1, 0, 0), 1)], [((0, 1, 0), 1)], [((0, 0, 1), 1)]]


def dfact(n):
    return math.prod(range(n, 0, -2))


def maps(power, moment, shift=(0, 0, 0)):
    denominator = dfact(3*moment+1)
    entries = []
    for i in range(3*power+1):
        j = i+shift[0]-shift[1]
        if j < 0:
            continue
        for k in range(max(0, power-i-j), 3*power-i-j+1):
            if (i+j+k-power) % 2 or (k+shift[2]) % 2:
                continue
            weight = math.factorial(i+shift[0])*dfact(k+shift[2]-1)/denominator
            entries.append((i, j, k, weight))
    return tuple(np.array([v[i] for v in entries], dtype=int if i < 3 else float) for i in range(4))


MOMENT_MAPS = {n: maps(n, n) for n in range(2, 21, 2)}
DERIVATIVE_MAPS = {(n, j): [(factor, maps(n-1, n, shift)) for shift, factor in family]
                   for n in range(2, 19, 2) for j, family in enumerate(BASIS)}


def evaluate(parameters, maximum=18):
    terms = dict(BASE)
    for value, family in zip(parameters, BASIS):
        for shift, factor in family:
            terms[shift] = terms.get(shift, 0)+value*factor
    values = []
    jacobian = []
    poly = np.ones((1, 1, 1), dtype=complex)
    for n in range(1, maximum+1):
        if n % 2 == 0 and n <= 18:
            jacobian.append([n*sum(factor*np.dot(poly[mp[0], mp[1], mp[2]], mp[3])
                                    for factor, mp in DERIVATIVE_MAPS[n, j])
                             for j in range(9)])
        size = 3*n+1
        next_poly = np.zeros((size, size, size), dtype=complex)
        old_size = poly.shape[0]
        for (i, j, k), coefficient in terms.items():
            next_poly[i:i+old_size, j:j+old_size, k:k+old_size] += coefficient*poly
        poly = next_poly
        if n % 2 == 0:
            mp = MOMENT_MAPS[n]
            values.append(np.dot(poly[mp[0], mp[1], mp[2]], mp[3]))
    return np.array(values), np.array(jacobian)


def main():
    start = time.monotonic()
    rng = np.random.default_rng(20261116)
    hom = json.loads((P/'homogeneous_prefix_witness_certificate.json').read_text())
    hom_center = np.array([complex(float(v[0]), float(v[1])) for v in hom['center']])
    attempts = []
    candidates = []
    output = P/'odd_prefix_witness_exploration.json'
    assert not output.exists(), 'Preserve completed attempts; use a new version for another campaign.'
    report = {'completed': False, 'kind': 'numerical exploration only',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seed': 20261116, 'maximum_starts': 24, 'maximum_iterations': 65,
              'wall_clock_limit_seconds': 150,
              'family': 'Z^3+ZWT+b W^3+c Z^2 T+d W^2 T+e(Z^2W-2ZT^2)+f(ZW^2-2WT^2)+h(T^3-3ZWT)+uZ+vW+wT',
              'coefficient_order': ['b', 'c', 'd', 'e', 'f', 'h', 'u', 'v', 'w'],
              'target_moments': list(range(2, 19, 2)), 'normalization': 'M_n/(3n+1)!!',
              'homogeneous_seed_sha256': hashlib.sha256((P/'homogeneous_prefix_witness_certificate.json').read_bytes()).hexdigest(),
              'attempts': attempts, 'candidates': candidates, 'exact_witness_certified': False}

    def checkpoint():
        report['runtime_seconds'] = round(time.monotonic()-start, 3)
        output.write_text(json.dumps(report, indent=2)+'\n')

    for index in range(24):
        if time.monotonic()-start > 150:
            break
        if index < 6:
            parameters = np.concatenate([hom_center, np.zeros(3, dtype=complex)])
            if index:
                parameters += .1*(rng.normal(size=9)+1j*rng.normal(size=9))
        else:
            parameters = .3*(rng.normal(size=9)+1j*rng.normal(size=9))
        best = float('inf')
        reason = 'iterations'
        for iteration in range(65):
            if time.monotonic()-start > 150:
                reason = 'wall_clock'; break
            values, jac = evaluate(parameters)
            row_scale = np.maximum(np.max(abs(jac), axis=1), 1e-100)
            if not np.all(np.isfinite(values)) or not np.all(np.isfinite(jac)):
                reason = 'nonfinite'; break
            norm = float(np.max(abs(values/row_scale)))
            best = min(best, norm)
            try:
                step = np.linalg.solve(jac/row_scale[:, None], values/row_scale)
            except np.linalg.LinAlgError:
                reason = 'singular'; break
            stepnorm = float(np.max(abs(step)))
            if norm < 1e-10 and stepnorm < 2e-6:
                all_values, _ = evaluate(parameters, maximum=20)
                candidates.append({'start': index, 'parameters': [[float(v.real), float(v.imag)] for v in parameters],
                                   'normalized_M2_to_M20': [[float(v.real), float(v.imag)] for v in all_values],
                                   'row_scaled_Jacobian_condition': float(np.linalg.cond(jac/row_scale[:, None])),
                                   'Newton_step_sup': stepnorm, 'exact': False})
                reason = 'candidate'; break
            if stepnorm > 8:
                step *= 8/stepnorm
            accepted = False
            for damping in [1, .5, .25, .125, .0625, .03125, .015625]:
                trial = parameters-damping*step
                if np.max(abs(trial)) > 25:
                    continue
                residual, _ = evaluate(trial)
                if np.linalg.norm(residual/row_scale, ord=2) < np.linalg.norm(values/row_scale, ord=2):
                    parameters = trial; accepted = True; break
            if not accepted:
                reason = 'line_search'; break
        attempts.append({'start': index, 'iterations': iteration+1, 'best_row_scaled_residual': best,
                         'reason': reason, 'last_parameters': [[float(v.real), float(v.imag)] for v in parameters]})
        checkpoint()
        print(json.dumps({'starts': index+1, 'candidates': len(candidates), 'best': best, 'reason': reason,
                          'elapsed': report['runtime_seconds']}), flush=True)
        if candidates:
            break
    report['completed'] = True
    report['terminal_reason'] = 'candidate' if candidates else ('wall_clock' if report['runtime_seconds'] >= 150 else 'start_budget')
    checkpoint()
    print(json.dumps({'completed': True, 'attempts': len(attempts), 'candidates': len(candidates),
                      'runtime_seconds': report['runtime_seconds'], 'exact_witness_certified': False}), flush=True)


if __name__ == '__main__':
    main()

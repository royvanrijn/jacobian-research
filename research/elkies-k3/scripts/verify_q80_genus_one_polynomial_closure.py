#!/usr/bin/env python3
"""Check source-bound finite premises of the written Q80 genus-one closure."""
import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-genus-one-polynomial-closure-v1'
NOTE = 'elkies-k3/Q80_GENUS_ONE_POLYNOMIAL_CLOSURE_2026-09-14.md'
POLE = 'artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1'
MIXED = 'artifacts/generated-results/elkies-k3-q80-genus-one-k0-orthogonal-contact-v1'
HELPER = 'elkies-k3/scripts/verify_q80_k0_pole_boundary.py'
SOURCES = {
    'k4': 'elkies-k3/COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md',
    'k3': 'elkies-k3/Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md',
    'k2': 'elkies-k3/Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md',
    'k1': 'elkies-k3/Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md',
    'k0_integral': 'elkies-k3/Q80_GENUS_ONE_K0_QUARTIC_GATE_2026-09-14.md',
    'k0_poles': 'elkies-k3/Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md',
    'k0_mixed': 'elkies-k3/Q80_GENUS_ONE_K0_MIXED_VALUATION_2026-09-14.md',
    'k0_two_negative': NOTE,
}


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bindings():
    paths = {NOTE, HELPER, *SOURCES.values(), str(Path(__file__).relative_to(ROOT))}
    for base, records in [(POLE, ['input.json', 'receipt.json']),
                          (MIXED, ['input.json', 'result.json', 'independent-replay.json'])]:
        paths.update(base + '/' + name for name in records)
        packet = read(ROOT / base / 'input.json')
        paths.update(packet['bindings'])
        paths.update(packet.get('preserved_preflight', {}))
    return {path: digest(ROOT / path) for path in sorted(paths)}


def freeze(out):
    out.mkdir(parents=True, exist_ok=True)
    packet = {
        'schema': 1, 'prime': 131, 'cpu_seconds': 20, 'memory_bytes': 2 * 1024**3,
        'scope': 'Finite nodal premises, old exact receipts, divisor-degree arithmetic and fifteen source-linked allocations. The local analytic closure requires the written proof.',
        'bindings': bindings(), 'coverage_sources': SOURCES,
    }
    with (out / 'input.json').open('x') as stream:
        json.dump(packet, stream, indent=2, sort_keys=True)
        stream.write('\n')


def verify(out):
    started = time.process_time()
    packet = read(out / 'input.json')
    require(packet['schema'] == 1 and packet['prime'] == 131, 'literal prime/schema')
    require(packet['cpu_seconds'] == 20 and packet['memory_bytes'] == 2 * 1024**3, 'declared finite limits')
    require(packet['bindings'] == bindings() and packet['coverage_sources'] == SOURCES, 'immutable written proof and premise bytes')
    spec = importlib.util.spec_from_file_location('q80_closure_pole_premises', ROOT / HELPER)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    pole = helper.verify(ROOT / POLE)
    old_pole = read(ROOT / POLE / 'receipt.json')
    require(old_pole['status'] == 'PASS' and old_pole['input_sha256'] == digest(ROOT / POLE / 'input.json'), 'retained pole input/receipt')
    for key in ['height_pole_table', 'excluded_k0_j_by_written_proof', 'remaining_k0_j',
                'remaining_pole_quadratic_reduction', 'new_analytic_proof_independently_replayed']:
        require(pole[key] == old_pole[key], 'unchanged pole component ' + key)
    require(old_pole['remaining_boundary'] == 'UNKNOWN', 'historical pole boundary retained')
    mixed = helper.check_prior_packet(ROOT / MIXED)
    replay = read(ROOT / MIXED / 'independent-replay.json')
    require(replay['status'] == 'PASS' and replay['input_sha256'] == digest(ROOT / MIXED / 'input.json')
            and replay['result_sha256'] == digest(ROOT / MIXED / 'result.json'), 'completed mixed finite replay binding')
    require(mixed['mixed_valuation_excluded_by_written_proof'] and mixed['remaining_valuation_pattern'] == 'both negative'
            and mixed['remaining_genus1_k0_j'] == [2] and mixed['remaining_boundary'] == 'UNKNOWN', 'unchanged mixed component boundary')
    require(replay['independent_finite_replay'] and not replay['analytic_proof_independently_verified'], 'finite versus analytic assurance')

    # Finite divisor arithmetic in the new proof. Its analytic premises are
    # not proved by this enumeration: they are explicitly sourced in NOTE.
    divisor_cases = []
    for pole_factor in ['F', 'G']:
        for alpha_zero_factor in ['F', 'G']:
            degree_F = (-2 if pole_factor == 'F' else 2) + (alpha_zero_factor == 'F')
            divisor_cases.append({'pole_factor': pole_factor, 'alpha_zero_factor': alpha_zero_factor,
                                  'degree_F': degree_F})
    selected = [row for row in divisor_cases if row['degree_F'] == 4 - 2]
    require(selected == [{'pole_factor': 'G', 'alpha_zero_factor': 'G', 'degree_F': 2}], 'unique divisor allocation with degree two')
    require(2 not in [-2, 3], 'coincident pole and nodal zero cannot give degree two')
    require(Fraction(1, 2)**2 == Fraction(1, 4) and pow(4, -1, 131) == 33, 'unit formal-doubling leading multiplier')
    require(pole['remaining_k0_j'] == [2] and pole['remaining_pole_quadratic_reduction'] == [88, 62, 1], 'quadratic pole premise')

    coverage = []
    for k in range(5):
        for j in range(5-k):
            if k:
                sources = [SOURCES['k' + str(k)]]
            else:
                sources = [SOURCES['k0_integral'], SOURCES['k0_poles']]
                if j == 2:
                    sources += [SOURCES['k0_mixed'], SOURCES['k0_two_negative']]
            coverage.append({'k': k, 'j': j, 'written_exclusion_sources': sources})
    require(len(coverage) == 15 and len({(row['k'], row['j']) for row in coverage}) == 15, 'all fifteen factor allocations')
    return {
        'status': 'PASS', 'input_sha256': digest(out / 'input.json'),
        'checker_sha256': digest(Path(__file__)), 'written_proof_sha256': digest(ROOT / NOTE),
        'scope': packet['scope'], 'nodal_premise_replay': pole['nodal_premise_replay'],
        'retained_pole_receipt_sha256': digest(ROOT / POLE / 'receipt.json'),
        'retained_mixed_replay_sha256': digest(ROOT / MIXED / 'independent-replay.json'),
        'divisor_cases': divisor_cases, 'forced_divisor_case': selected[0],
        'written_allocation_coverage': coverage,
        'genus1_polynomial_chart_excluded_over_Q_by_written_proof': True,
        'new_two_negative_case_excluded_over_Q131_by_written_proof': True,
        'analytic_proof_independently_verified': False, 'formal_verification': False,
        'external_review': False, 'large_censuses_rerun': False, 'denominator_sweep_run': False,
        'old_norm8_singular_replay_gap_upgraded': False, 'positive_mw17_target_complete': False,
        'cpu_seconds': time.process_time() - started,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['freeze', 'check'])
    parser.add_argument('--directory', type=Path, default=DEFAULT)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (20, 25))
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    if args.mode == 'freeze':
        freeze(args.directory)
    else:
        result = verify(args.directory)
        if args.receipt:
            with args.receipt.open('x') as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write('\n')
        print(json.dumps(result), flush=True)

#!/usr/bin/env python3
"""Independently replay the new49-pencil rational singular-member exclusion.

Uses finite integer arithmetic, the existing independent trace checker, and
the expanded quartic discriminant, rather than the producer's I,J formula.
The inherited norm8 exclusion and the written descent remain dependencies.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-genus-zero-singular-gate-v1'
OLD = ROOT / 'artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1'
SPEC = importlib.util.spec_from_file_location('moving_replay', ROOT / 'elkies-k3/scripts/verify_q80_norm12_moving_pencils.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
read, require, ev = V.read, V.require, V.ev


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def quartic_discriminant(coeff, p):
    """Discriminant of a binary quartic, including its roots at infinity."""
    require(len(coeff) == 5, 'binary quartic needs five coefficients')
    e,d,c,b,a = coeff
    return (256*a**3*e**3 - 192*a*a*b*d*e*e - 128*a*a*c*c*e*e
            + 144*a*a*c*d*d*e - 27*a*a*d**4 + 144*a*b*b*c*e*e
            - 6*a*b*b*d*d*e - 80*a*b*c*c*d*e + 18*a*b*c*d**3
            + 16*a*c**4*e - 4*a*c**3*d*d - 27*b**4*e*e
            + 18*b**3*c*d*e - 4*b**3*d**3 - 4*b*b*c**3*e
            + b*b*c*c*d*d) % p


def root_multiplicity(coeff, root, p):
    current, count = coeff[:], 0
    while current and not current[-1]:
        current.pop()
    require(current, 'zero discriminant is not a root certificate')
    while len(current) > 1 and ev(current, root, p) == 0:
        quotient = [0]*(len(current)-1)
        carry = current[-1]
        for j in range(len(current)-2,-1,-1):
            quotient[j] = carry
            carry = (current[j]+root*carry) % p
        require(carry == 0, 'root division remainder')
        current, count = quotient, count+1
    return count


def check_discriminant(row, frame, p):
    coeff = row['discriminant_coefficients']
    require(len(coeff) == 25 and all(type(c) is int and 0 <= c < p for c in coeff), 'degree24 binary discriminant coefficients')
    require(any(coeff), 'identically zero discriminant')
    require(p > 24, 'distinct interpolation sites required')
    # Both expressions have degree <=24. 25 agreements are an exact identity.
    for u in range(25):
        quartic = [ev(r,u,p) for r in frame['branch_matrix']]
        require(ev(coeff,u,p) == quartic_discriminant(quartic,p), 'discriminant interpolation identity')
    infinity_value = quartic_discriminant([r[4] for r in frame['branch_matrix']],p)
    require(coeff[24] == infinity_value, 'parameter infinity discriminant')
    finite = [u for u in range(p) if ev(coeff,u,p) == 0]
    roots = [[u,root_multiplicity(coeff,u,p)] for u in finite]
    infinity_order = 24-max(i for i,c in enumerate(coeff) if c)
    require(roots == row['finite_roots_with_multiplicity'], 'full finite root list and multiplicities')
    require(infinity_order == row['infinity_multiplicity'], 'lost projective infinity root')
    return not roots and not infinity_order


def verify(path):
    start = time.monotonic()
    packet, result = read(path/'input.json'), read(path/'result.json')
    require(packet['checker_sha256'] == digest(Path(__file__)), 'checker source binding')
    producer = ROOT/'elkies-k3/scripts/certify_q80_genus_zero_singular_gate.sage'
    require(packet['producer_sha256'] == digest(producer), 'producer source binding')
    for name, expected in packet['bindings'].items():
        require(digest(ROOT/name) == expected, 'input binding '+name)
    required = [str((OLD/name).relative_to(ROOT)) for name in (
        'input.json.gz','frames-521.json.gz','frames-523.json.gz',
        'prime-521.json.gz','prime-523.json.gz','independent-replay.json')]
    required += ['artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json',
                 'artifacts/generated-results/elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json']
    require(set(required) <= set(packet['bindings']), 'missing retained proof binding')
    require(packet['primes'] == [521,523] and packet['pencil_indices'] == list(range(49)), 'complete49-pencil input')
    require(result['input_sha256'] == digest(path/'input.json'), 'result input binding')
    require(result['status'] == 'PASS' and result['norm12_pencils_excluded'] == 49, 'completion status')
    require(result['positive_target_complete'] is False and result['whole_family_genus_bound'] == 'UNKNOWN', 'scope widened beyond witnesses')
    source = read(OLD/'input.json.gz')
    parity = read(ROOT/required[-2])['invariant_trace_parity']
    require(parity['isotropic_minimum_norm_partition'] == {'4':1313,'8':63917,'12':49}, 'retained parity partition')
    require(source['words'] == [r['section_basis_w'] for r in parity['deep_norm12_classes']], 'all49 trace words')
    prior = read(ROOT/required[-1])
    require(prior['minimum_translation_class_count'] == 63917 and prior['covered_half_open_range'] == [0,63917] and prior['candidate_count'] == 0, 'inherited norm8 scope')
    require(read(OLD/'independent-replay.json')['new_trace_frames_verified'] == 98, 'inherited frame replay receipt')
    remaining, witnesses = set(range(49)), {}
    count, parameters = 0, 0
    summaries = []
    require([s['prime'] for s in result['stages']] == packet['primes'], 'missing or extra prime stage')
    for meta in result['stages']:
        p = meta['prime']
        name = 'prime-%d.json' % p
        require(meta['path'] == name and digest(path/name) == meta['sha256'], 'prime stage binding')
        stage = read(path/name)
        require(stage['prime'] == p and stage['input_sha256'] == digest(path/'input.json'), 'prime input attachment')
        frames_path = OLD/('frames-%d.json.gz' % p)
        old_frames = read(frames_path)
        require(old_frames['prime'] == p and old_frames['input_sha256'] == digest(OLD/'input.json.gz'), 'frame attachment')
        require(stage['frames_sha256'] == digest(frames_path) == read(OLD/('prime-%d.json.gz' % p))['frames_sha256'], 'retained frame hash')
        require([r['index'] for r in old_frames['records']] == list(range(49)), 'frame completeness')
        selected = sorted(remaining)
        require(stage['selected_indices'] == selected and [r['index'] for r in stage['records']] == selected, 'adaptive coverage')
        frames = {r['index']:r for r in old_frames['records']}
        ctx = V.trace_context(source,p)
        V.validate_norm12([frames[i] for i in selected],ctx)
        excluded = []
        for row in stage['records']:
            i = row['index']
            if check_discriminant(row,frames[i],p):
                excluded.append(i)
                witnesses[str(i)] = p
        remaining.difference_update(excluded)
        require(stage['excluded_indices'] == excluded and stage['remaining_indices'] == sorted(remaining), 'exact stage conclusion')
        require((meta['selected'],meta['excluded'],meta['remaining']) == (len(selected),len(excluded),len(remaining)), 'stage counts')
        count += len(selected)
        parameters += len(selected)*(p+1)
        summaries.append({'prime':p,'selected':len(selected),'excluded':len(excluded),'remaining':len(remaining)})
    require(not remaining and witnesses == result['witnesses'] and len(witnesses) == 49, 'complete witness assignment')
    rows = result['genus_zero_strata']
    require([(r['k'],r['j']) for r in rows] == [(k,j) for k in range(3) for j in range(5-k)], 'all12 degree strata')
    for r in rows:
        k,j = r['k'],r['j']
        require(r['degree_bounds'] == [k,2-k,j,4-k-j,5-j,1+k+j,4,0,0], 'genus-zero degree allocation')
        I = k+2*j
        require((r['intersection'],r['cross_height'],r['sum_height'],r['difference_height']) == (I,4-I,24-2*I,8+2*I), 'genus-zero height allocation')
    return {'status':'PASS', 'input_sha256':digest(path/'input.json'),
        'result_sha256':digest(path/'result.json'), 'new_discriminant_replay_independent':True,
        'trace_frames_reverified':count, 'degree24_discriminants_verified':count,
        'projective_parameter_tests':parameters, 'norm12_pencils_excluded':49,
        'stages':summaries, 'positive_target_complete':False,
        'inherited_norm8_independent_replay':False,
        'written_descent_formally_verified':False,
        'elapsed_seconds':round(time.monotonic()-start,6)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT)
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    result = verify(args.input)
    if args.record:
        with args.record.open('x') as stream:
            json.dump(result,stream,indent=2,sort_keys=True)
            stream.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)

#!/usr/bin/env sage-python
"""Replay the bounded descent diagnostic, theorem target and exact norm relations."""
import argparse
import json
from pathlib import Path
import runpy

import certify_small_conductor_selmer_rank_target as parity
import prepare_small_conductor_norm_form as forms
import pilot_small_conductor_norm_smoothness as pilot

ROOT, ART = forms.ROOT, forms.ART
OUT = ART / 'small_conductor_descent_shortcut_v1.json'


def expected():
    cert = forms.target.original.cert
    for module, path in [(forms.target, forms.target.OUTPUT), (parity, parity.OUTPUT), (forms, forms.OUT)]:
        if json.loads(json.dumps(module.expected())) != cert.read(path):
            raise ArithmeticError('exact component certificate differs: ' + str(path))
    pilot.calculate(check=True)
    relation_checker = ROOT / 'elliptic-curves/cas/audit_small_conductor_norm_relations.sage'
    relations = runpy.run_path(str(relation_checker))
    actual = relations['expected']()
    if actual != cert.read(relations['OUTPUT']):
        raise ArithmeticError('principal-ideal relations differ')
    stages = []
    inputs = {Path(__file__).resolve(), forms.target.OUTPUT, parity.OUTPUT, forms.OUT,
              pilot.OUT, relations['OUTPUT'], relation_checker}
    for stage in ['field', 'cold_rankinit', 'hinted_rankinit', 'prepared_bnf']:
        path = forms.D / (stage + '.json')
        row = cert.read(path)
        log, program = forms.D / (stage + '.log'), forms.D / (stage + '.gp')
        if cert.hashed(log) != row['log_sha256'] or cert.hashed(program) != row['program_sha256']:
            raise ArithmeticError('profile evidence hash differs')
        required = 'completed' if stage == 'field' else 'strict_wall_timeout'
        if row['supervision']['outcome'] != required or row['unconditional_rank_upper_bound'] is not None:
            raise ArithmeticError('profile outcome differs')
        stages.append({'stage': stage, 'outcome': required,
                       'wall_seconds': row['supervision']['wall_seconds'],
                       'limit_seconds': row['supervision']['timeout_seconds']})
        inputs.update([path, log, program])
    return {'schema': 'elliptic-curves.small-conductor-descent-shortcut.v1', 'status': 'PASS',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in sorted(inputs)},
            'profile_stages': stages,
            'profile_replay_scope': 'Retained program/log hashes and bounded outcomes are checked; slow profiling is not rerun. The old gp_error flag is ignored because normal PARI debug lines contain ***.',
            'class_two_rank_lower_bound': 15,
            'sufficient_class_two_rank_upper_bound_for_exact_rank22': 16,
            'maximum_coefficient_bits': cert.read(forms.OUT)['maximum_coefficient_bits'],
            'smoothness_arms': cert.read(pilot.OUT)['arms'],
            'exact_principal_relations': len(actual['relations']),
            'additional_mod2_relation_rank': actual['additional_relation_rank'],
            'factor_base_generation_certified': False,
            'class_rank_upper_bound': None, 'curve_rank_upper_bound': None,
            'claim_boundary': 'Known-theorem specialization, exact integral norm coordinates, and a fixed finite relation-yield experiment on MW16 family05 at3/17. No full descent, class upper bound, new point, exact rank, or general runtime guarantee.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = expected()
    cert = forms.target.original.cert
    if args.check:
        if result != cert.read(OUT):
            raise ArithmeticError('aggregate certificate differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve aggregate certificate')
        cert.write(OUT, result)
    print('DESCENT SHORTCUT REPLAY PASS: criterion and18 relations; rank upper bound UNKNOWN')

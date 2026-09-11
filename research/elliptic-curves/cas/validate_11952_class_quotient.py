#!/usr/bin/env python3
"""Validate the frozen large-stack lane, distinguishing diagnostics from errors.

PARI prints both its stack-size announcement and Bach-constant diagnostic with
***. The frozen launcher's conservative classifier rejects those too. This
separate, hash-bound validator permits ONLY those two exact diagnostic forms;
it preserves the original log and result and never relaxes arithmetic gates.
"""
import argparse
import json
from pathlib import Path
import re

import run_11952_class_quotient as lane
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run


def normalize_diagnostics(log):
    allowed = [r'\s*\*\*\*\s+Warning: new maximum stack size = 8589934592 \(8192\.000 Mbytes\)\.',
               r'\*\*\* Bach constant: [0-9]+\.[0-9]+']
    return '\n'.join(line for line in log.splitlines()
                     if not any(re.fullmatch(pattern, line) for pattern in allowed))


def classify(stage, log, outcome, binary_present, provisional=None):
    result = lane.classify(stage, normalize_diagnostics(log),
                           'rss_limit' if outcome == 'strict_rss_limit' else outcome,
                           binary_present, provisional)
    return result


def check_protocol():
    p = lane.read(lane.OUT / 'protocol.json')
    for name, expected in p['input_sha256'].items():
        assert lane.digest(Path(name)) == expected, name
    for stage, expected in p['program_sha256'].items():
        assert lane.digest(lane.OUT / (stage+'.gp')) == expected
    assert p['wall_seconds'] == {'provisional': 600, 'quotient': 300}
    assert p['pari_stack_bytes'] == lane.STACK and p['rss_bytes'] == lane.RSS
    return p


def validate(stage, check=False):
    protocol = check_protocol()
    raw = lane.read(lane.OUT / (stage+'.json'))
    assert raw['supervision'] == lane.read(lane.OUT / (stage+'.supervisor.json'))
    assert raw['supervision']['command'] == [str(lane.GP), '-q', '-f', '-s', str(lane.STACK), str(lane.OUT / (stage+'.gp'))]
    for key, value in [('wall_seconds', protocol['wall_seconds'][stage]),
                       ('rss_bytes', lane.RSS), ('pari_stack_bytes', lane.STACK)]:
        assert raw['supervision']['limits'][key] == value
    log_path = lane.OUT / (stage+'.log')
    assert raw['log_sha256'] == lane.digest(log_path)
    assert raw['protocol_sha256'] == lane.digest(lane.OUT / 'protocol.json')
    binary = lane.LOCAL / 'provisional_bnf.bin'
    if raw['bnf_sha256'] is not None:
        assert binary.exists() and raw['bnf_sha256'] == lane.digest(binary)
    else:
        assert not binary.exists()
    provisional = validate('provisional', check) if stage == 'quotient' else None
    result = classify(stage, log_path.read_text(), raw['supervision']['outcome'], binary.exists(), provisional)
    result.update(stage=stage, original_result_sha256=lane.digest(lane.OUT / (stage+'.json')),
                  validator_sha256=lane.digest(Path(__file__)),
                  protocol_sha256=lane.digest(lane.OUT / 'protocol.json'),
                  log_sha256=lane.digest(log_path), bnf_sha256=raw['bnf_sha256'],
                  original_supervision_outcome=raw['supervision']['outcome'],
                  independence='Separate output validation, not independent BNF arithmetic; exact rank depends on the pinned local criterion and independent point replay.')
    path = lane.OUT / (stage+'.validation.json')
    if path.exists():
        assert lane.read(path) == result
    elif check:
        raise FileNotFoundError(path)
    else:
        checkpoint(path, result)
    return result


def run_quotient():
    p = check_protocol()
    provisional = validate('provisional')
    assert provisional['status'] == 'PROVISIONAL_GRH_SUFFICIENT'
    assert lane.memory_available() >= lane.RSS + 2*2**30
    for suffix in ['json', 'supervisor.json', 'log']:
        if (lane.OUT / ('quotient.'+suffix)).exists():
            raise FileExistsError('no repeat of quotient stage')
    checkpoint(lane.OUT / 'quotient-admission.json', {
        'validator_sha256': lane.digest(Path(__file__)),
        'protocol_sha256': lane.digest(lane.OUT / 'protocol.json'),
        'provisional_validation_sha256': lane.digest(lane.OUT / 'provisional.validation.json'),
        'g_comp': provisional['provisional_class_2rank'],
        'reason': 'Completed sufficient provisional BNF; only two known nonerror diagnostics excluded. Frozen GP code, memory and time limits unchanged.'})
    supervision = run([str(lane.GP), '-q', '-f', '-s', str(lane.STACK), str(lane.OUT / 'quotient.gp')],
                      limits=Limits(p['wall_seconds']['quotient'], lane.RSS, pari_stack_bytes=lane.STACK),
                      log_path=lane.OUT / 'quotient.log', checkpoint_path=lane.OUT / 'quotient.supervisor.json', cwd=lane.ROOT)
    checkpoint(lane.OUT / 'quotient.json', {
        'stage': 'quotient', 'supervision': supervision,
        'protocol_sha256': lane.digest(lane.OUT / 'protocol.json'),
        'log_sha256': lane.digest(lane.OUT / 'quotient.log'),
        'bnf_sha256': lane.digest(lane.LOCAL / 'provisional_bnf.bin')})
    return validate('quotient')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--run-quotient', action='store_true')
    args = parser.parse_args()
    if args.run_quotient:
        print(json.dumps(run_quotient(), indent=2))
    else:
        found = False
        for stage in ['provisional', 'quotient']:
            if (lane.OUT / (stage+'.json')).exists():
                found = True
                print(json.dumps(validate(stage, args.check), indent=2))
        if not found:
            raise SystemExit('No completed stage record yet; no validation certificate issued.')

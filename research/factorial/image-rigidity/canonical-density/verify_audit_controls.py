#!/usr/bin/env python3
"""Replay optimized checkers and reject deliberately corrupted private copies."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
import json
import re
import shutil
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parent
CASES = [
    ('verify_density_exclusion.py', 'canonical_density_exclusion.md', 'density_exclusion_verification.json'),
    ('verify_inverse_reachability.py', 'inverse_canonical_reachability.md', 'inverse_reachability_verification.json'),
    ('verify_moment_pairs.py', 'moment_pair_obstructions.md', 'moment_pair_verification.json'),
]
INPUTS = [name for row in CASES for name in row[:2]] + ['moment_pair_certificates.json']

def require(value, message):
    if not value:
        raise ArithmeticError(message)

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def run(checker, output):
    return subprocess.run([sys.executable, '-O', str(checker), '--output', str(output)],
                          capture_output=True, text=True)

def fixture_copy(directory):
    for name in INPUTS:
        shutil.copyfile(BASE / name, directory / name)

def negative(label, checker, mutate, expected):
    with tempfile.TemporaryDirectory(prefix='sfc-audit-') as name:
        directory = Path(name)
        fixture_copy(directory)
        mutate(directory)
        output = directory / 'unexpected_success.json'
        result = run(directory / checker, output)
        require(result.returncode != 0 and expected in result.stderr,
                'Wrong rejection for ' + label + ': ' + result.stdout + result.stderr)
        require(not output.exists(), 'Success receipt written for corrupted case')
        return dict(test=label, status='REJECTED', optimized=True,
                    error=result.stderr.strip().splitlines()[-1], success_receipt_written=False)

def changed_proof(name):
    def change(directory):
        p = directory / name
        p.write_text(p.read_text() + '\nDeliberately changed proof copy.\n')
    return change

def replace_code(old, new, checker):
    def change(directory):
        p = directory / checker
        content = p.read_text()
        require(content.count(old) == 1, 'Mutation locator must be unique')
        p.write_text(content.replace(old, new))
    return change

def changed_certificate(kind):
    def change(directory):
        path = directory / 'moment_pair_certificates.json'
        data = json.loads(path.read_text())
        row = data['generic_counterpairs'][0]
        if kind == 'moment':
            row['upper_coefficients_descending'][1] = str(
                Fraction(row['upper_coefficients_descending'][1]) + 1)
        else:
            row['bezout_upper_multiplier_ascending'][0] += 1
        path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
        checker = directory / 'verify_moment_pairs.py'
        content, count = re.subn(r"(?m)^CERTIFICATE_SHA256='[a-f0-9]+'$",
                                 "CERTIFICATE_SHA256='" + digest(path) + "'",
                                 checker.read_text())
        require(count == 1, 'Isolated fixture pin locator')
        checker.write_text(content)
    return change

before = {name: digest(BASE / name) for name in INPUTS}
replays = []
for checker, proof, receipt in CASES:
    with tempfile.TemporaryDirectory(prefix='sfc-replay-') as name:
        output = Path(name) / 'optimized.json'
        result = run(BASE / checker, output)
        require(result.returncode == 0, result.stdout + result.stderr)
        require(output.read_bytes() == (BASE / receipt).read_bytes(),
                'Normal/optimized receipt mismatch: ' + checker)
        replays.append(dict(checker=checker, normal_optimized_byte_identical=True))

mutations = [negative('altered proof: ' + proof, checker, changed_proof(proof),
                      'proof authentication') for checker, proof, receipt in CASES]
mutations.append(negative('moment coefficient with copied fixture pin refreshed',
                         'verify_moment_pairs.py', changed_certificate('moment'),
                         'independent factorial-moment reconstruction'))
mutations.append(negative('Bezout coefficient with copied fixture pin refreshed',
                         'verify_moment_pairs.py', changed_certificate('bezout'),
                         'modular Bezout identity'))
mutations.append(negative('wrong backward leading coefficient',
                         'verify_inverse_reachability.py',
                         replace_code('m*(m+D)', 'm*(m+D+1)', 'verify_inverse_reachability.py'),
                         'backward top-coefficient cancellation'))
mutations.append(negative('wrong generic collision factor', 'verify_moment_pairs.py',
                         replace_code('B=F(9,748)', 'B=F(9,749)', 'verify_moment_pairs.py'),
                         'genuine generic collision'))
after = {name: digest(BASE / name) for name in INPUTS}
require(before == after, 'Production input was changed')
result = dict(status='PASS', source_sha256=before, source_bytes_unchanged=True,
              replays=replays, deliberate_mutations_rejected=len(mutations),
              mutations=mutations, external_review=False, formal_verification=False,
              scope='Internal exact replay and negative controls; no unrestricted SFC closure.')
(BASE / 'portable_verification.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
print(json.dumps(dict(status='PASS', optimized_replays=len(replays),
                     deliberate_mutations_rejected=len(mutations))))

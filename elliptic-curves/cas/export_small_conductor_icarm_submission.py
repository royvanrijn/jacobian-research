#!/usr/bin/env python3
"""Export the four ICARM form fields from the verified rank22 certificate."""
import argparse
from pathlib import Path

import certify_small_conductor_curve as proof_checker

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
PROOF = ART / 'small_conductor_rank22_proof_v1.json'
OUTPUT = ART / 'small_conductor_rank22_icarm_submission.md'
POINTS = ART / 'small_conductor_rank22_icarm_points.txt'


def export(check=False):
    proof = proof_checker.cert.read(PROOF)
    proof_checker.verify(proof)
    coefficients = ', '.join(proof['integral_model'])
    points = '\n'.join(', '.join(p) for p in proof['integral_points']) + '\n'
    primes = ', '.join(p for p, _ in proof['discriminant_factorization'])
    commentary = (
        'Specialization of the compact A1/MW16 family a1-fibration-05 at parameter 3/17. '
        'The supplied 22 rational points have an exact finite-reduction independence certificate, '
        'proving rank at least 22; exact rank is unknown. '
        'The displayed integral equation is globally minimal. '
        'The complete minimal-discriminant factorization is proved with recursive Lucas primality certificates. '
        'All bad primes except 17 have multiplicative reduction and conductor exponent 1; '
        'at 17 the minimal discriminant valuation is 4 and the conductor exponent is 2. '
        'The exact conductor is ' + proof['conductor'] + '. '
        'Local research identifier: new-20260905-36. '
        'Exact point, conductor and primality certificates are available.'
    )
    text = (
        '# ICARM submission: MW16 at 3/17\n\n'
        'Open https://elliptic-rank.icarm.cloud/ and log in with GitHub. '
        'Paste each block into its corresponding form field.\n\n'
        '## a-invariants\n\n```text\n' + coefficients + '\n```\n\n'
        '## points\n\nAll 22 points below are on the displayed integral model. '
        'Paste the complete block, one `x, y` pair per line.\n\n```text\n' + points + '```\n\n'
        '## primes of bad reduction\n\nPaste the ten distinct primes, without their discriminant exponents.\n\n```text\n'
        + primes + '\n```\n\n'
        '## commentary\n\n```text\n' + commentary + '\n```\n\n'
        '## Expected result\n\n'
        '- Certified rank lower bound: at least 22.\n'
        '- Exact conductor: `' + proof['conductor'] + '`.\n'
        '- Exact rank: unknown.\n\n'
        'The form checks the points and their independence and records the conductor from the supplied bad primes. '
        'The certificate is supporting evidence; the form itself asks for the four fields above. '
        'A public proof URL can be added to the commentary when one is available.\n\n'
        'Source: `small_conductor_rank22_proof_v1.json`, SHA256 `' + proof_checker.cert.hashed(PROOF) + '`.\n\n'
        'Regenerate from the repository root:\n\n```sh\n'
        'python3 elliptic-curves/cas/export_small_conductor_icarm_submission.py\n```\n'
    )
    # Verify that the literal form payload round-trips to the certified model,
    # transported points and complete set of bad primes.
    assert [v.strip() for v in coefficients.split(',')] == proof['integral_model']
    assert [[v.strip() for v in line.split(',')] for line in points.splitlines()] == proof['integral_points']
    assert [v.strip() for v in primes.split(',')] == [p for p, _ in proof['discriminant_factorization']]
    for path, content in [(OUTPUT, text), (POINTS, points)]:
        if check:
            if path.read_text() != content:
                raise ArithmeticError('submission field export differs: ' + str(path))
        else:
            path.write_text(content)
    print('VERIFIED ICARM FIELDS: five coefficients,22 points,10 bad primes and commentary')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    export(parser.parse_args().check)
